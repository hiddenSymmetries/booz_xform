#!/usr/bin/env python3

import unittest
import os
import tempfile
import numpy as np
from booz_xform import Booz_xform
from scipy.io import netcdf_file

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class WriteReadTest(unittest.TestCase):
    def test_flux_arrays_preserve_vmec_values_and_units(self):
        """
        VMEC phi and chi are signed fluxes in Wb, not values divided by
        2*pi. Check their native import and boozmn serialization against an
        independent NetCDF reader so a hidden normalization cannot survive a
        Booz_xform-only write/read comparison.
        """
        wout_filename = os.path.join(TEST_DIR, 'wout_circular_tokamak.nc')
        with netcdf_file(wout_filename, 'r', mmap=False) as wout:
            phi_vmec = np.array(wout.variables['phi'].data, copy=True)
            chi_vmec = np.array(wout.variables['chi'].data, copy=True)
            phipf_vmec = np.array(wout.variables['phipf'].data, copy=True)
        phip_expected = -phipf_vmec / (2 * np.pi)

        b1 = Booz_xform()
        b1.read_wout(wout_filename, True)
        np.testing.assert_array_equal(b1.phi, phi_vmec)
        np.testing.assert_array_equal(b1.chi, chi_vmec)
        np.testing.assert_allclose(b1.phip, phip_expected,
                                   rtol=2e-16, atol=0)

        b1.compute_surfs = [0, 5]
        b1.run()
        with tempfile.TemporaryDirectory() as tmpdir:
            boozmn_filename = os.path.join(tmpdir, 'boozmn_flux_units.nc')
            b1.write_boozmn(boozmn_filename)

            with netcdf_file(boozmn_filename, 'r', mmap=False) as boozmn:
                phi_b = boozmn.variables['phi_b']
                chi_b = boozmn.variables['chi_b']
                phip_b = boozmn.variables['phip_b']
                np.testing.assert_array_equal(phi_b.data, phi_vmec)
                np.testing.assert_array_equal(chi_b.data, chi_vmec)
                phip_expected[0] = 0
                np.testing.assert_allclose(phip_b.data, phip_expected,
                                           rtol=2e-16, atol=0)
                self.assertEqual(phi_b.units, b'Tesla * meter^2')
                self.assertEqual(chi_b.units, b'Tesla * meter^2')
                self.assertEqual(phip_b.units, b'Tesla * meter^2')

            b2 = Booz_xform()
            b2.read_boozmn(boozmn_filename)
            np.testing.assert_array_equal(b2.phi, phi_vmec)
            np.testing.assert_array_equal(b2.chi, chi_vmec)
            np.testing.assert_allclose(b2.phip, phip_expected,
                                       rtol=2e-16, atol=0)

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
                    if flux:
                        np.testing.assert_allclose(b1.phip,b2.phip, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.chi,b2.chi, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.pres,b2.pres, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.phi,b2.phi, rtol=rtol, atol=atol)
                        np.testing.assert_allclose(b1.toroidal_flux, b2.toroidal_flux, rtol=rtol, atol=atol)
                    else:
                        np.testing.assert_equal(b2.phip,0)
                        assert len(b1.phip) == 0
                        np.testing.assert_equal(b2.chi,0)
                        assert len(b1.chi) == 0
                        np.testing.assert_equal(b2.pres,0)
                        assert len(b1.pres) == 0
                        np.testing.assert_equal(b2.phi,0)
                        assert len(b1.phi) == 0
                        np.testing.assert_equal(b1.toroidal_flux,0)
                        np.testing.assert_equal(b2.toroidal_flux,0)
                    os.remove(boozmn_filename)

if __name__ == '__main__':
    unittest.main()
