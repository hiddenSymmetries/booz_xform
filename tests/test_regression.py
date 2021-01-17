#!/usr/bin/env python3

import unittest
import os
import numpy as np
from scipy.io import netcdf
from booz_xform import Booz_xform

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class RegressionTest(unittest.TestCase):
    def test_li383(self):
        f = netcdf.netcdf_file(os.path.join(TEST_DIR, 'boozmn_li383_1.4m.nc'),
                               'r', mmap=False)
        b = Booz_xform()
        b.read_wout(os.path.join(TEST_DIR, 'wout_li383_1.4m.nc'))
        # Transfer parameters from the reference file to the new
        # calculation
        b.mboz = f.variables['mboz_b'][()]
        b.nboz = f.variables['nboz_b'][()]
        b.jlist = f.variables['jlist'][()]

        b.run()

        vars = ['bmnc_b', 'rmnc_b', 'zmns_b', 'pmns_b', 'gmnc_b']
        
        rtol = 1e-12
        atol = 1e-12
        for var in vars:
            # gmnc_b is misspelled in the fortran version
            var_ref = var
            if var == 'gmnc_b':
                var_ref = 'gmn_b'
                
            # Reference values:
            arr1 = f.variables[var_ref][()]
            # Newly computed values:
            arr2 = np.array(getattr(b, var)).transpose()
            
            print('abs diff in ' + var + ':', np.max(np.abs(arr1 - arr2)))
            np.testing.assert_allclose(arr1, arr2,
                                       rtol=rtol, atol=atol)
        f.close()

if __name__ == '__main__':
    unittest.main()
