import numpy as np
from ._booz_xform import Booz_xform as Booz_xform_cpp


def _calculate_fourier_series(m, n, cos_mn_ampl, sin_mn_ampl, theta, phi):

    out = np.zeros_like(phi)

    for jmn in range(len(m)):

        m_i = m[jmn]
        n_i = n[jmn]
        angle = m_i * theta - n_i * phi

        out += cos_mn_ampl[jmn] * np.cos(angle)
        if not sin_mn_ampl.size == 0:
            out += sin_mn_ampl[jmn] * np.sin(angle)

    return out


class Booz_xform(Booz_xform_cpp):

    def calculate_modB_boozer_on_surface(self, js, theta, phi):
        """
        Calculates :math:`|B|` on a surface in Boozer poloidal and toroidal
        angles.
        Args:
          js (int): The index among the output surfaces.
          phi (array-like): The toroidal angle values.
          theta (array-like): The poloidal angle values.
        """

        phi = np.asanyarray(phi, dtype=float)
        theta = np.asanyarray(theta, dtype=float)

        cos_ampl = self.bmnc_b[:, js]

        if self.asym:
            sin_ampl = self.bmns_b[:, js]
        else:
            sin_ampl = np.array([])

        return _calculate_fourier_series(self.xm_b,
                                         self.xn_b,
                                         cos_ampl,
                                         sin_ampl,
                                         theta,
                                         phi,
                                         )
