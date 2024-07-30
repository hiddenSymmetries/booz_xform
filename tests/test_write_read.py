#!/usr/bin/env python3

import unittest
import os
import numpy as np
from booz_xform import Booz_xform

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class WriteReadTest(unittest.TestCase):
    def test_write_read(self):
        """
        Write a boozmn file, then read in the contents of that file to a
        second Booz_xform object. The data in the two objects should match.
        """
        configurations = ['circular_tokamak',
                          'up_down_asymmetric_tokamak',
                          'li383_1.4m',
                          'LandremanSenguptaPlunk_section5p3']

        for configuration in configurations:
            # Try a few different array sizes:
            for compute_surfs in [[0], [0, 5], [5, 10, 15]]:
                wout_filename = 'wout_' + configuration + '.nc'
                boozmn_filename = 'boozmn_new_' + configuration + '.nc'
                b1 = Booz_xform()
                for flux in [True,False]:
                    b1.read_wout(os.path.join(TEST_DIR, wout_filename),flux)
                    b1.compute_surfs = compute_surfs
                    b1.run()
                    b1.write_boozmn(boozmn_filename)

                    # Read the results into a new object
                    b2 = Booz_xform()
                    b2.read_boozmn(boozmn_filename)

                    self.assertEqual(b1.asym, b2.asym)
                    self.assertEqual(b1.nfp, b2.nfp)
                    self.assertEqual(b1.mboz, b2.mboz)
                    self.assertEqual(b1.nboz, b2.nboz)
                    self.assertEqual(b1.mnboz, b2.mnboz)
                    self.assertEqual(b1.ns_b, b2.ns_b)
                    np.testing.assert_equal(b1.xm_b, b2.xm_b)
                    np.testing.assert_equal(b1.xn_b, b2.xn_b)
                    np.testing.assert_equal(b1.compute_surfs, b2.compute_surfs)

                    rtol = 1e-15
                    atol = 1e-15
                    np.testing.assert_allclose(b1.bmnc_b, b2.bmnc_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.bmns_b, b2.bmns_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.gmnc_b, b2.gmnc_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.gmns_b, b2.gmns_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.rmnc_b, b2.rmnc_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.rmns_b, b2.rmns_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.zmnc_b, b2.zmnc_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.zmns_b, b2.zmns_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.numnc_b, b2.numnc_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.numns_b, b2.numns_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.s_b, b2.s_b, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.iota, b2.iota, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.Boozer_G_all, b2.Boozer_G_all, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.Boozer_I_all, b2.Boozer_I_all, rtol=rtol, atol=atol)
                    np.testing.assert_allclose(b1.phi,b2.phi, rtol=rtol, atol=atol)
                    if flux:
                        np.testing.assert_allclose(b1.phip,b2.phip, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.chi,b2.chi, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.pres,b2.pres, rtol=rtol, atol=atol)
                    else:
                        np.testing.assert_equal(b2.phip,0)
                        assert len(b1.phip) == 0
                        np.testing.assert_equal(b2.chi,0)
                        assert len(b1.chi) == 0
                    os.remove(boozmn_filename)

if __name__ == '__main__':
    unittest.main()
