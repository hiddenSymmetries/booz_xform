#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import booz_xform as bx

wout_filename = 'test_files/wout_li383_1.4m.nc'
b = bx.Booz_xform()
b.read_wout(wout_filename)
b.compute_surfs = [47]
b.run()
bx.surfplot(b)
plt.tight_layout()

plt.figure()
bx.surfplot(b, fill=False, cmap=plt.cm.jet, levels=np.arange(1.3, 2.0, 0.05))
plt.tight_layout()

plt.show()
