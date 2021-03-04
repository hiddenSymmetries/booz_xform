#!/usr/bin/env python3

import matplotlib.pyplot as plt
import booz_xform as bx

if False:
    # Use a Booz_xform object that was just used to drive the transformation:
    wout_filename = 'test_files/wout_li383_1.4m.nc'
    b = bx.Booz_xform()
    b.read_wout(wout_filename)
    #b.compute_surfs = [2, 4]
    b.run()
else:
    # Use a boozmn_*.nc file from a previous transformation:
    b = 'test_files/boozmn_li383_1.4m_manySurfs.nc'

plt.figure()
bx.modeplot(b, nmodes=6)

plt.figure()
bx.modeplot(b, sqrts=True, log=False, B0=False,
            legend_args={"fontsize":10})

plt.show()
