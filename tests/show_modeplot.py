#!/usr/bin/env python3

import matplotlib.pyplot as plt
import booz_xform as bx

wout_filename = 'test_files/wout_li383_1.4m.nc'
b = bx.Booz_xform()
b.read_wout(wout_filename)
#b.compute_surfs = [2, 4]
b.run()

plt.figure()
bx.modeplot(b, nmodes=5, legend_args={"fontsize":7})

plt.figure()
bx.modeplot(b, marker='.', sqrts=True, log=False, B0=False)

plt.show()
