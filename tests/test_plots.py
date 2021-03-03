#!/usr/bin/env python3

import os
import unittest
import numpy as np
import matplotlib.pyplot as plt
import booz_xform as bx

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_files')

class PlotTest(unittest.TestCase):
    def test_calling_plots(self):
        """
        Try calling the plotting functions. Because plt.show() is not
        called, no plots will actually appear on the screen.
        """
        configurations = ['circular_tokamak',
                          'up_down_asymmetric_tokamak',
                          'li383_1.4m',
                          'LandremanSenguptaPlunk_section5p3']

        for configuration in configurations:
            for which_input in [0, 1]:
                if which_input == 0:
                    # Supply a filename string
                    b = os.path.join(TEST_DIR,
                                     "boozmn_" + configuration + ".nc")
                else:
                    # Supply a Booz_xform instance that was used to
                    # drive the transformation:
                    wout_filename = 'wout_' + configuration + '.nc'
                    b = bx.Booz_xform()
                    b.read_wout(os.path.join(TEST_DIR, wout_filename))
                    b.compute_surfs = [2, 15]
                    b.run()

                # Try a few different options for each plot type:
                bx.surfplot(b)
                bx.surfplot(b, js=1, fill=False, cmap=plt.cm.jet)

                bx.symplot(b)
                bx.symplot(b, marker='.', sqrts=True, log=False)

                bx.modeplot(b)
                bx.modeplot(b, marker='.', sqrts=True, log=False,
                            legend_args={'fontsize':6})

if __name__ == '__main__':
    unittest.main()
