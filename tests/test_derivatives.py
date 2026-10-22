#!/usr/bin/env python3

import os
import unittest
import numpy as np
from scipy.io import netcdf_file
from booz_xform import Booz_xform, InputArrays, OutputArrays

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

INPUTS = ['iota', 'rmnc', 'rmns', 'zmnc', 'zmns', 'lmnc', 'lmns',
          'bmnc', 'bmns', 'bsubumnc', 'bsubumns', 'bsubvmnc', 'bsubvmns']
OUTPUTS = ['bmnc_b', 'bmns_b', 'rmnc_b', 'rmns_b', 'zmnc_b', 'zmns_b',
           'numnc_b', 'numns_b', 'gmnc_b', 'gmns_b', 'Boozer_G', 'Boozer_I']

# (wout file, compute_surfs, mboz, nboz): one stellarator-symmetric and one
# asymmetric case.
CASES = [('wout_li383_1.4m.nc', [3, 20, 40], 10, 6),
         ('wout_up_down_asymmetric_tokamak.nc', [1, 5], 8, 0)]


def setup(filename, surfs, mboz, nboz):
    b = Booz_xform()
    b.verbose = 0
    b.read_wout(os.path.join(TEST_DIR, filename))
    b.compute_surfs = surfs
    b.mboz = mboz
    b.nboz = nboz
    b.run()
    return b


def present(obj, names):
    return {k: np.array(getattr(obj, k)) for k in names if np.size(getattr(obj, k)) > 0}


def to_struct(cls, arrays):
    s = cls()
    for k, v in arrays.items():
        setattr(s, k, v)
    return s


def random_like(arrays, rng):
    """Random perturbations scaled with the magnitude of each entry."""
    return {k: rng.standard_normal(v.shape) * (np.abs(v) + 1e-3 * np.abs(v).max())
            for k, v in arrays.items()}


def dot(a, b):
    return sum(np.vdot(a[k], b[k]) for k in a)


class DerivativeTest(unittest.TestCase):
    def test_jvp_matches_finite_differences(self):
        """jvp() agrees with 4th-order central differences of run()."""
        for filename, surfs, mboz, nboz in CASES:
            b = setup(filename, surfs, mboz, nboz)
            x0 = present(b, INPUTS)
            t = random_like(x0, np.random.default_rng(0))
            dy = present(b.jvp(to_struct(InputArrays, t)), OUTPUTS)

            def outputs_at(eps):
                for k in x0:
                    setattr(b, k, x0[k] + eps * t[k])
                b.run()
                return present(b, OUTPUTS)

            h = 1e-4
            yp, ym = outputs_at(h), outputs_at(-h)
            yp2, ym2 = outputs_at(2 * h), outputs_at(-2 * h)
            self.assertEqual(set(dy), set(yp))
            for k in dy:
                fd = (8 * (yp[k] - ym[k]) - (yp2[k] - ym2[k])) / (12 * h)
                scale = np.abs(dy[k]).max()
                np.testing.assert_allclose(dy[k], fd, rtol=0, atol=1e-8 * scale,
                                           err_msg=f'{filename} {k}')

    def test_vjp_is_transpose_of_jvp(self):
        for filename, surfs, mboz, nboz in CASES:
            b = setup(filename, surfs, mboz, nboz)
            rng = np.random.default_rng(1)
            t = random_like(present(b, INPUTS), rng)
            dy = present(b.jvp(to_struct(InputArrays, t)), OUTPUTS)
            ybar = {k: rng.standard_normal(v.shape) for k, v in dy.items()}
            xbar = present(b.vjp(to_struct(OutputArrays, ybar)), INPUTS)
            self.assertEqual(set(xbar), set(t))
            lhs, rhs = dot(ybar, dy), dot(xbar, t)
            self.assertAlmostEqual(lhs / rhs, 1.0, delta=1e-12, msg=filename)

    def test_init_from_vmec_derivatives(self):
        """init_from_vmec_jvp() reproduces init_from_vmec() on perturbed
        wout arrays, and init_from_vmec_vjp() is its transpose."""
        for filename, _, _, _ in CASES:
            b = Booz_xform()
            b.verbose = 0
            b.read_wout(os.path.join(TEST_DIR, filename))
            with netcdf_file(os.path.join(TEST_DIR, filename), mmap=False) as f:
                ns = int(f.variables['ns'][()])
                wout = {k: f.variables[k][()].copy().T for k in INPUTS
                        if k != 'iota' and k in f.variables}
                wout['iota'] = f.variables['iotas'][()].copy()
            if not b.asym:
                for k in ['rmns', 'zmnc', 'lmnc', 'bmns', 'bsubumns', 'bsubvmns']:
                    wout.pop(k, None)
            rng = np.random.default_rng(2)
            t = random_like(wout, rng)
            d = present(b.init_from_vmec_jvp(to_struct(InputArrays, t)), INPUTS)

            def internal(arrays):
                full = {k: arrays.get(k, np.zeros((1, 1))) for k in INPUTS}
                b.init_from_vmec(ns, full['iota'], full['rmnc'], full['rmns'],
                                 full['zmnc'], full['zmns'], full['lmnc'], full['lmns'],
                                 full['bmnc'], full['bmns'], full['bsubumnc'],
                                 full['bsubumns'], full['bsubvmnc'], full['bsubvmns'])
                return present(b, INPUTS)

            # init_from_vmec is linear, so applying it to the tangent gives the jvp.
            ref = internal(t)
            self.assertEqual(set(d), set(ref))
            for k in d:
                np.testing.assert_allclose(d[k], ref[k], rtol=1e-14, atol=0,
                                           err_msg=f'{filename} {k}')

            xbar = {k: rng.standard_normal(v.shape) for k, v in d.items()}
            wbar = present(b.init_from_vmec_vjp(to_struct(InputArrays, xbar)), INPUTS)
            self.assertEqual(set(wbar), set(t))
            self.assertAlmostEqual(dot(xbar, d) / dot(wbar, t), 1.0, delta=1e-13,
                                   msg=filename)

    def test_empty_members_are_zero_and_shapes_are_checked(self):
        b = setup(*CASES[0])
        rng = np.random.default_rng(3)
        t = InputArrays()
        t.bmnc = rng.standard_normal(np.array(b.bmnc).shape)
        dy = b.jvp(t)
        # Only |B| was perturbed, so the geometry and nu are unchanged.
        np.testing.assert_array_equal(dy.rmnc_b, 0)
        np.testing.assert_array_equal(dy.numns_b, 0)
        np.testing.assert_array_equal(dy.Boozer_G, 0)
        self.assertGreater(np.abs(dy.bmnc_b).max(), 0)
        self.assertEqual(np.size(dy.bmns_b), 0)

        t.bmnc = np.zeros((2, 2))
        with self.assertRaises(RuntimeError):
            b.jvp(t)
        with self.assertRaises(RuntimeError):
            Booz_xform().vjp(OutputArrays())


if __name__ == '__main__':
    unittest.main()
