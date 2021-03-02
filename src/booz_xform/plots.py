#!/usr/bin/env python

import numpy as np
# No error is raised if matplotlib cannot be imported here, because
# matplotlib is not needed for core operation of booz_xform, so we
# don't want to make it a mandatory requirement.
try:
    import matplotlib.pyplot as plt
except:
    pass

def surfplot(bx,
             js = 0,
             ntheta = 50,
             nphi = 90,
             ncontours = 25,
             **kwargs):
    """
    Plot :math:`|B|` on a surface vs the Boozer poloidal and toroidal angles.

    Args:
      bx (Booz_xform): The Booz_xform instance to plot.
      js (int): The index among the output surfaces to plot.
      ntheta (int): Number of grid points in the poloidal angle.
      nphi (int): Number of grid points in the toroidal angle.
      ncontours (int): Number of contours to show.
      kwargs: Any additional key-value pairs to pass to matplotlib's ``contourf`` command.

    This function can generate figures like this:

    .. image:: surfplot.png
       :width: 400

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
    plt.contourf(phi, theta, B, ncontours, **kwargs)
    cbar = plt.colorbar()
    plt.xlabel(r'Boozer toroidal angle $\varphi$')
    plt.ylabel(r'Boozer poloidal angle $\theta$')
    plt.title('|B| on surface $s$={:.4} [Tesla]'.format(bx.s_b[js]))


def symplot(bx,
            max_m = 20,
            max_n = 20,
            ymin = None,
            sqrts = False,
            log = True,
            B0 = True,
            legend_loc = "best",
            **kwargs):
    """
    Plot the radial variation of all the Fourier modes of :math:`|B|`
    in Boozer coordinates. Color is used to group modes with
    :math:`m=0` and/or :math:`n=0`.

    Args:
      bx (Booz_xform): The Booz_xform instance to plot.
      max_m (int): Maximum poloidal mode number to include in the plot.
      max_n (int): Maximum toroidal mode number (divided by nfp) to include in the plot.
      ymin (float): Lower limit for the y-axis. Only used if ``log==True``.
      sqrts (bool): If true, the x axis will be sqrt(toroidal flux) instead of toroidal flux.
      log (bool): Whether to use a logarithmic y axis.
      B0 (bool): Whether to include the m=n=0 mode in the figure.
      legend_loc (str): Location of the legend.
      kwargs: Any additional key-value pairs to pass to matplotlib's ``plot`` command.
    """

    background_color = 'b'
    QA_color = [0, 0.7, 0]
    mirror_color = [0.7, 0.5, 0]
    helical_color = [1, 0, 1]

    # If ymin is not specified, pick a default value such that the
    # plot mostly shows the largest modes, not all the modes down to
    # machine precision.
    if ymin is None:
        ymin = np.max(bx.bmnc_b) * 1e-4
    
    nmodes = len(bx.xm_b)

    # This next line needs to be corrected!
    s = np.linspace(0, 1, len(bx.compute_surfs))

    if sqrts:
        rad = np.sqrt(bx.s_b)
    else:
        rad = bx.s_b

    def my_abs(x):
        if log:
            return np.abs(x)
        else:
            return x

    # Draw a reference line at 0.
    if not log:
        plt.plot([0, 1], [0, 0], ':k')
        
    # First, plot just the 1st mode of each type, so the legend looks nice.
    if B0:
        for imode in range(nmodes):
            if bx.xn_b[imode] == 0 and bx.xm_b[imode] == 0:
                plt.plot(rad, my_abs(bx.bmnc_b[imode, :]), color=background_color,
                             label='m = 0, n = 0 (Background)', **kwargs)
                break
    for imode in range(nmodes):
        if bx.xn_b[imode] == 0 and bx.xm_b[imode] != 0:
            plt.plot(rad, my_abs(bx.bmnc_b[imode, :]), color=QA_color,
                         label=r'm $\ne$ 0, n = 0 (Quasiaxisymmetric)', **kwargs)
            break
    for imode in range(nmodes):
        if bx.xn_b[imode] != 0 and bx.xm_b[imode] == 0:
            plt.plot(rad, my_abs(bx.bmnc_b[imode, :]), color=mirror_color,
                         label=r'm = 0, n $\ne$ 0 (Mirror)', **kwargs)
            break
    for imode in range(nmodes):
        if bx.xn_b[imode] != 0 and bx.xm_b[imode] != 0:
            plt.plot(rad, my_abs(bx.bmnc_b[imode, :]), color=helical_color,
                         label=r'm $\ne$ 0, n $\ne$ 0 (Helical)', **kwargs)
            break

    plt.legend(fontsize=9, loc=legend_loc)
    # Now that the legend is made, plot all modes

    for imode in range(nmodes):
        if np.abs(bx.xm_b[imode]) > max_m:
            continue
        if np.abs(bx.xn_b[imode]) > max_n * bx.nfp:
            continue
        if bx.xn_b[imode] == 0:
            if bx.xm_b[imode] == 0:
                mycolor = background_color
                if not B0:
                    continue
            else:
                mycolor = QA_color
        else:
            if bx.xm_b[imode] == 0:
                mycolor = mirror_color
            else:
                mycolor = helical_color
        plt.plot(rad, my_abs(bx.bmnc_b[imode, :]), color=mycolor, **kwargs)

    if sqrts:
        plt.xlabel('$r/a$ = sqrt(Normalized toroidal flux)')
    else:
        plt.xlabel('$s$ = Normalized toroidal flux')
    plt.title('Fourier harmonics of |B| in Boozer coordinates [Tesla]')
    plt.xlim([0, 1])
    if log:
        plt.yscale("log")
        plt.gca().set_ylim(bottom=ymin)
