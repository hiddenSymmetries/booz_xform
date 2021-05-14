import numpy as np
from ._booz_xform import Booz_xform as Booz_xform_cpp
from .plots import surfplot, symplot, modeplot, wireplot



class Booz_xform(Booz_xform_cpp):

    def calculate_modB_boozer(self, js, phi, theta):
        """
        Calculates :math:`|B|` on a surface in Boozer poloidal and toroidal angles.
        Args:
          js (int): The index among the output surfaces to plot.
          phi (array-like): The toroidal angle values.
          theta (array-like): The poloidal angle values.
        """

        phi = np.asanyarray(phi, dtype=np.float)
        theta = np.asanyarray(theta, dtype=np.float)

        modB = np.zeros_like(phi)

        for jmn in range(len(self.xm_b)):
            m = self.xm_b[jmn]
            n = self.xn_b[jmn]
            angle = m * theta - n * phi
            modB += self.bmnc_b[jmn, js] * np.cos(angle)
            if self.asym:
                modB += self.bmns_b[jmn, js] * np.sin(angle)
        return modB


