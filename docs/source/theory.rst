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

We assume that at the start, each flux surface is parameterized by a poloidal angle :math:`\theta_0`
and a toroidal angle :math:`\zeta_0`, which need not be straight-field-line angles. The toroidal
angle need not be the standard toroidal angle (although it often is). We assume that
we know the quantity :math:`\lambda` by which the poloidal angle can be shifted to obtain
the straight-field-line angle :math:`\theta^*=\theta_0 + \lambda` compatible with :math:`\zeta_0`.
Thus, the magnetic field :math:`\vec{B}` can be written as

.. math::

   \vec{B} = \nabla\psi\times\nabla\theta^* + \iota(\psi) \nabla\zeta_0\times\nabla\psi

where :math:`2 \pi \psi` is the toroidal flux, and :math:`\iota` is the rotational transform.
The Boozer poloidal angle :math:`\theta_B` and toroidal angle :math:`\zeta_B` are also straight field line
angles, so

.. math::

   \vec{B} = \nabla\psi\times\nabla\theta_B + \iota(\psi) \nabla\zeta_B\times\nabla\psi.

Moreover, Boozer coordinates satisfy

.. math::
   
   \vec{B} = \beta(\psi,\theta_B,\zeta_B)\nabla\psi + I(\psi)\nabla\theta_B + G(\psi)\nabla\zeta_B



.. note:: In Hirshman's note "Transformation from VMEC to Boozer coordinates", the angle difference
	  is defined as :math:`p = \zeta_B - \zeta_0`. However in the fortran ``booz_xform`` code,
	  a minus sign appears on line 83 of ``boozer.f`` that reverses the sign, so the ``p`` quantity
	  saved in ``boozmn_*.nc`` files is in fact :math:`\zeta_0 - \zeta_B`.
