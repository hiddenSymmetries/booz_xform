#!/usr/bin/env python3

import os
import unittest
import booz_xform as bx

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class PlotTest(unittest.TestCase):
    def test_surfplot(self):
        """
        Try calling surfplot(). Because plt.show() is not called, no plots
        will actually appear on the screen.
        """
        configurations = ['circular_tokamak',
                          'up_down_asymmetric_tokamak',
                          'li383_1.4m',
                          'LandremanSenguptaPlunk_section5p3']

        for configuration in configurations:
            wout_filename = 'wout_' + configuration + '.nc'
            b = bx.Booz_xform()
            b.read_wout(os.path.join(TEST_DIR, wout_filename))
            b.compute_surfs = [2, 4]
            b.run()
            bx.surfplot(b)

if __name__ == '__main__':
    unittest.main()
