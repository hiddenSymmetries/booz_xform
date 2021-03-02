#!/usr/bin/env python

import numpy as np
# No error is raised if matplotlib cannot be imported here, because
# matplotlib is not needed for core operation of booz_xform, so we
# don't want to make it a mandatory requirement.
try:
    import matplotlib.pyplot as plt
except:
    pass

def hw():
    """ For testing. """
    print("Hello world from surfplot")

def surfplot(bx,
             js = 0,
             ntheta = 50,
             nphi = 90,
             ncontours = 25):
    """
    Plot |B| on a surface vs the Boozer poloidal and toroidal angles.

    Args:
      bx (Booz_xform): A Booz_xform instance.
    """

    theta1d = np.linspace(0, 2 * np.pi, ntheta)
    phi1d = np.linspace(0, 2 * np.pi / bx.nfp, nphi)
    phi, theta = np.meshgrid(phi1d, theta1d)

    B = np.zeros_like(phi)
    for jmn in range(len(bx.xm_b)):
        m = bx.xm_b[jmn]
        n = bx.xn_b[jmn]
        angle = m * theta - n * phi
        B += bx.bmnc_b[jmn, js] * np.cos(angle)
        if bx.asym:
            B += bx.bmns_b[jmn, js] * np.sin(angle)

    plt.rcParams.update({'font.size': 16})
    plt.contourf(phi, theta, B, ncontours)
    cbar = plt.colorbar()
    plt.xlabel(r'Boozer toroidal angle $\varphi$')
    plt.ylabel(r'Boozer poloidal angle $\theta$')
    plt.title('|B| on surface {}'.format(bx.compute_surfs[js]))
