#!/usr/bin/env python3

import matplotlib.pyplot as plt
import booz_xform as bx

wout_filename = 'test_files/wout_li383_1.4m.nc'
b = bx.Booz_xform()
b.read_wout(wout_filename)
#b.compute_surfs = [2, 4]
b.run()

plt.figure()
bx.symplot(b)

plt.figure()
bx.symplot(b, marker='.', sqrts=True, log=False)

plt.show()
