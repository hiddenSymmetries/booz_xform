#!/usr/bin/env python

import numpy as np
# No error is raised if matplotlib cannot be imported here, because
# matplotlib is not needed for core operation of booz_xform, so we
# don't want to make it a mandatory requirement.
try:
    import matplotlib.pyplot as plt
except:
    pass

from ._booz_xform import Booz_xform

def handle_b_input(b):
    if isinstance(b, str):
        filename = b
        b2 = Booz_xform()
        b2.read_boozmn(filename)
        return b2
    
    elif isinstance(b, Booz_xform):
        return b
    
    else:
        raise ValueError("b argument must be a booz_xform.Booz_xform instance or string")

def surfplot(b,
             js = 0,
             fill = True,
             ntheta = 50,
             nphi = 90,
             ncontours = 25,
             **kwargs):
    """
    Plot :math:`|B|` on a surface vs the Boozer poloidal and toroidal angles.

    Args:
      b (Booz_xform, str): The Booz_xform instance to plot,
        or a filename of a boozmn_*.nc file.
      js (int): The index among the output surfaces to plot.
      fill (bool): Whether the contours are filled, i.e.
        whether to use plt.contourf vs plt.contour.
      ntheta (int): Number of grid points in the poloidal angle.
      nphi (int): Number of grid points in the toroidal angle.
      ncontours (int): Number of contours to show.
      kwargs: Any additional key-value pairs to pass to matplotlib's 
        ``contourf`` or ``contour`` command.

    This function can generate figures like this:

    .. image:: surfplot.png
       :width: 400

    """
    b = handle_b_input(b)
    
    theta1d = np.linspace(0, 2 * np.pi, ntheta)
    phi1d = np.linspace(0, 2 * np.pi / b.nfp, nphi)
    phi, theta = np.meshgrid(phi1d, theta1d)

    modB = np.zeros_like(phi)
    for jmn in range(len(b.xm_b)):
        m = b.xm_b[jmn]
        n = b.xn_b[jmn]
        angle = m * theta - n * phi
        modB += b.bmnc_b[jmn, js] * np.cos(angle)
        if b.asym:
            modB += b.bmns_b[jmn, js] * np.sin(angle)

    plt.rcParams.update({'font.size': 16})
    if fill:
        plt.contourf(phi, theta, modB, ncontours, **kwargs)
    else:
        plt.contour(phi, theta, modB, ncontours, **kwargs)
        
    cbar = plt.colorbar()
    plt.xlabel(r'Boozer toroidal angle $\varphi$')
    plt.ylabel(r'Boozer poloidal angle $\theta$')
    plt.title('|B| on surface $s$={:.4} [Tesla]'.format(b.s_b[js]))


def symplot(b,
            max_m = 20,
            max_n = 20,
            ymin = None,
            sqrts = False,
            log = True,
            B0 = True,
            legend_args = {"loc":"best"},
            **kwargs):
    """
    Plot the radial variation of all the Fourier modes of :math:`|B|`
    in Boozer coordinates. Color is used to group modes with
    :math:`m=0` and/or :math:`n=0`.

    Args:
      b (Booz_xform, str): The Booz_xform instance to plot,
        or a filename of a boozmn_*.nc file.
      max_m (int): Maximum poloidal mode number to include in the plot.
      max_n (int): Maximum toroidal mode number (divided by nfp) to include in the plot.
      ymin (float): Lower limit for the y-axis. Only used if ``log==True``.
      sqrts (bool): If true, the x axis will be sqrt(toroidal flux) instead of toroidal flux.
      log (bool): Whether to use a logarithmic y axis.
      B0 (bool): Whether to include the m=n=0 mode in the figure.
      legend_args (dict): Any arguments to pass to ``plt.legend()``.
         Useful for setting the legend font size and location.
      kwargs: Any additional key-value pairs to pass to matplotlib's ``plot`` command.
    """

    b = handle_b_input(b)

    background_color = 'b'
    QA_color = [0, 0.7, 0]
    mirror_color = [0.7, 0.5, 0]
    helical_color = [1, 0, 1]

    # If ymin is not specified, pick a default value such that the
    # plot mostly shows the largest modes, not all the modes down to
    # machine precision.
    if ymin is None:
        ymin = np.max(b.bmnc_b) * 1e-4
    
    mnmax = len(b.xm_b)

    if sqrts:
        rad = np.sqrt(b.s_b)
    else:
        rad = b.s_b

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
        for imode in range(mnmax):
            if b.xn_b[imode] == 0 and b.xm_b[imode] == 0:
                plt.plot(rad, my_abs(b.bmnc_b[imode, :]), color=background_color,
                             label='m = 0, n = 0 (Background)', **kwargs)
                break
    for imode in range(mnmax):
        if b.xn_b[imode] == 0 and b.xm_b[imode] != 0:
            plt.plot(rad, my_abs(b.bmnc_b[imode, :]), color=QA_color,
                         label=r'm $\ne$ 0, n = 0 (Quasiaxisymmetric)', **kwargs)
            break
    for imode in range(mnmax):
        if b.xn_b[imode] != 0 and b.xm_b[imode] == 0:
            plt.plot(rad, my_abs(b.bmnc_b[imode, :]), color=mirror_color,
                         label=r'm = 0, n $\ne$ 0 (Mirror)', **kwargs)
            break
    for imode in range(mnmax):
        if b.xn_b[imode] != 0 and b.xm_b[imode] != 0:
            plt.plot(rad, my_abs(b.bmnc_b[imode, :]), color=helical_color,
                         label=r'm $\ne$ 0, n $\ne$ 0 (Helical)', **kwargs)
            break

    plt.legend(**legend_args)
    # Now that the legend is made, plot all modes

    for imode in range(mnmax):
        if np.abs(b.xm_b[imode]) > max_m:
            continue
        if np.abs(b.xn_b[imode]) > max_n * b.nfp:
            continue
        if b.xn_b[imode] == 0:
            if b.xm_b[imode] == 0:
                mycolor = background_color
                if not B0:
                    continue
            else:
                mycolor = QA_color
        else:
            if b.xm_b[imode] == 0:
                mycolor = mirror_color
            else:
                mycolor = helical_color
        plt.plot(rad, my_abs(b.bmnc_b[imode, :]), color=mycolor, **kwargs)

    if sqrts:
        plt.xlabel('$r/a$ = sqrt(Normalized toroidal flux)')
    else:
        plt.xlabel('$s$ = Normalized toroidal flux')
    plt.title('Fourier harmonics of |B| in Boozer coordinates [Tesla]')
    plt.xlim([0, 1])
    if log:
        plt.yscale("log")
        plt.gca().set_ylim(bottom=ymin)

        
def modeplot(b,
             nmodes = 10,
             ymin = None,
             sqrts = False,
             log = True,
             B0 = True,
             legend_args = {"loc":"best"},
             **kwargs):
    """
    Plot the radial variation of the Fourier modes of :math:`|B|` in
    Boozer coordinates. The plot includes only the largest few modes,
    based on their magnitude at the outermost surface for which data
    are available.

    Args:
      b (Booz_xform, str): The Booz_xform instance to plot,
        or a filename of a boozmn_*.nc file.
      nmodes (int): How many modes to include
      ymin (float): Lower limit for the y-axis. Only used if ``log==True``.
      sqrts (bool): If true, the x axis will be sqrt(toroidal flux) instead of toroidal flux.
      log (bool): Whether to use a logarithmic y axis.
      B0 (bool): Whether to include the m=n=0 mode in the figure.
      legend_args (dict): Any arguments to pass to ``plt.legend()``.
         Useful for setting the legend font size and location.
      kwargs: Any additional key-value pairs to pass to matplotlib's ``plot`` command.
    """
    b = handle_b_input(b)

    # If ymin is not specified, pick a default value such that the
    # plot mostly shows the largest modes, not all the modes down to
    # machine precision.
    if ymin is None:
        ymin = np.max(b.bmnc_b) * 1e-4

    # Either keep or discard the m=n=0 mode:
    assert b.xm_b[0] == 0
    assert b.xn_b[0] == 0
    if B0:
        data = b.bmnc_b
        xm = b.xm_b
        xn = b.xn_b
    else:
        data = b.bmnc_b[1:, :]
        xm = b.xm_b[1:]
        xn = b.xn_b[1:]

    # Sort by the absolute value of the mode amplitudes at the
    # outermost radius. Minus sign gives a decreasing instead of
    # increasing sort.
    sorted_indices = np.argsort(-np.abs(data[:, -1]))
    indices = sorted_indices[:nmodes]

    if sqrts:
        rad = np.sqrt(b.s_b)
    else:
        rad = b.s_b

    def my_abs(x):
        if log:
            return np.abs(x)
        else:
            return x

    # Draw a reference line at 0.
    if not log:
        plt.plot([0, 1], [0, 0], ':k')
        
    for jmode in range(nmodes):
        index = indices[jmode]
        plt.plot(rad, my_abs(data[index, :]),
                 label='$m=${}, $n=${}'.format(xm[index], xn[index]),
                 **kwargs)

    plt.legend(**legend_args)

    if sqrts:
        plt.xlabel('$r/a$ = sqrt(Normalized toroidal flux)')
    else:
        plt.xlabel('$s$ = Normalized toroidal flux')
    plt.title('Fourier harmonics of |B| in Boozer coordinates [Tesla]')
    plt.xlim([0, 1])
    if log:
        plt.yscale("log")
        plt.gca().set_ylim(bottom=ymin)
