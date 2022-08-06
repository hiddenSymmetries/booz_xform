#!/usr/bin/env python3

import unittest
import os
import numpy as np
from scipy.io import netcdf_file
from booz_xform import Booz_xform

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class RegressionTest(unittest.TestCase):
    def test_regression(self):

        configurations = ['circular_tokamak',
                          'up_down_asymmetric_tokamak',
                          'li383_1.4m',
                          'LandremanSenguptaPlunk_section5p3']
        
        for configuration in configurations:
            wout_filename = 'wout_' + configuration + '.nc'
            boozmn_filename = 'boozmn_' + configuration + '.nc'
            boozmn_new_filename = 'boozmn_new_' + configuration + '.nc'
            f = netcdf_file(os.path.join(TEST_DIR, boozmn_filename),
                                   'r', mmap=False)
            b = Booz_xform()
            b.read_wout(os.path.join(TEST_DIR, wout_filename))
            # Transfer parameters from the reference file to the new
            # calculation
            b.mboz = f.variables['mboz_b'][()]
            b.nboz = f.variables['nboz_b'][()]
            b.compute_surfs = f.variables['jlist'][()] - 2

            b.run()

            # Compare 2D arrays
            vars = ['bmnc_b', 'rmnc_b', 'zmns_b', 'numns_b', 'gmnc_b']
            asym = bool(f.variables['lasym__logical__'][()])
            if asym:
                vars += ['bmns_b', 'rmns_b', 'zmnc_b', 'numnc_b', 'gmns_b']

            rtol = 1e-12
            atol = 1e-12
            for var in vars:
                # gmnc_b is misspelled in the fortran version
                var_ref = var
                if var == 'gmnc_b':
                    var_ref = 'gmn_b'

                # Handle the issue that we now use the variable nu,
                # whereas the boozmn format uses the variable
                # p = -nu.
                sign = 1
                if var[:2] == 'nu':
                    sign = -1
                    var_ref = 'p' + var[2:]

                # Reference values:
                arr1 = f.variables[var_ref][()]
                # Newly computed values:
                arr2 = getattr(b, var).transpose()

                print('abs diff in ' + var + ':', np.max(np.abs(arr1 - sign * arr2)))
                np.testing.assert_allclose(arr1, sign * arr2,
                                           rtol=rtol, atol=atol)

            # Now compare some values written to the boozmn files.
            b.write_boozmn(boozmn_new_filename)
            f2 = netcdf_file(boozmn_new_filename)

            vars = f.variables.keys()
            # These variables will not match:
            exclude = ['rmax_b', 'rmin_b', 'betaxis_b', 'version', 'pres_b', 'beta_b', 'phip_b']
            for var in vars:
                if var in exclude:
                    continue
                # Reference values:
                arr1 = f.variables[var][()]
                # Newly computed values:
                arr2 = f2.variables[var][()]
                print('abs diff in ' + var + ':', np.max(np.abs(arr1 - arr2)))
                np.testing.assert_allclose(arr1, arr2,
                                           rtol=rtol, atol=atol)
            
            f.close()

if __name__ == '__main__':
    unittest.main()
