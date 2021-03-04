#!/usr/bin/env python3

import booz_xform as bx

wout_filename = 'test_files/wout_li383_1.4m.nc'
b = bx.Booz_xform()
b.read_wout(wout_filename)
#b.mboz = 64
#b.nboz = 64
b.compute_surfs = [47]
b.run()
#fig = bx.wireplot(b, surf=False)
#fig.show()
bx.wireplot(b, refine=6, orig=False).show()
#bx.wireplot(b, refine=6, surf=False).show()
