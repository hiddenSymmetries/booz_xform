Theory and numerical implementation
===================================

Theory
******

The code is based on an unpublished 1995 note by Steve Hirshman,
"Transformation from VMEC to Boozer coordinates".
(The method is not specific to VMEC, and can be applied to other equilibrium representations.)
The calculation has two steps. In the first step, the difference is
computed between the original toroidal and poloidal angles and the Boozer angles.
In the second step, this information is used to transform other quantities such as
:math:`R`, :math:`Z`, and :math:`B` that are known as functions of the original
angles to functions of the new angles.

The Boozer poloidal angle :math:`\theta_B` and toroidal angle :math:`\zeta_B` are defined by the relations

.. math::

   \vec{B} = \nabla\psi\times\nabla\theta_B + \iota(\psi) \nabla\zeta_B\times\nabla\psi

   \vec{B} = \beta(\psi,\theta_B,\zeta_B)\nabla\psi + I(\psi)\nabla\theta_B + G(\psi)\nabla\zeta_B

where :math:`2 \pi \psi` is the toroidal flux, and :math:`\iota` is the rotational transform.
