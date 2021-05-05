Theory and numerical implementation
===================================

Theory
******

The code is based on an unpublished 1995 note by Steve Hirshman,
:download:`"Transformation from VMEC to Boozer coordinates".<Transformation_from_VMEC_to_Boozer_coordinates.pdf>`
(The method is not specific to VMEC, and can be applied to other equilibrium representations.)
The calculation has two steps. In the first step, the difference is
computed between the original toroidal and poloidal angles and the Boozer angles.
In the second step, this information is used to transform other quantities such as
:math:`R`, :math:`Z`, and :math:`B` that are known as functions of the original
angles to functions of the new angles.

1. Determining the toroidal angle difference
--------------------------------------------

We assume that at the start, each flux surface is parameterized by a poloidal angle :math:`\theta_0`
and a toroidal angle :math:`\zeta_0`, which need not be straight-field-line angles. The toroidal
angle need not be the standard toroidal angle (although it typically is). We assume that
we know the quantity :math:`\lambda` by which the poloidal angle can be shifted to obtain
the straight-field-line angle :math:`\theta^*=\theta_0 + \lambda` compatible with :math:`\zeta_0`.
Thus, the magnetic field :math:`\vec{B}` can be written as

.. math:: \vec{B} = \nabla\psi\times\nabla\theta^* + \iota(\psi) \nabla\zeta_0\times\nabla\psi
   :label: straight_field_line

where :math:`2 \pi \psi` is the toroidal flux, and :math:`\iota` is the rotational transform.
The Boozer poloidal angle :math:`\theta_B` and toroidal angle :math:`\zeta_B` are also straight field line
angles, so

.. math::
   :label: straight_field_line_B

   \vec{B} = \nabla\psi\times\nabla\theta_B + \iota(\psi) \nabla\zeta_B\times\nabla\psi.

Moreover, Boozer coordinates satisfy

.. math::
   :label: Boozer
   
   \vec{B} = \beta(\psi,\theta_B,\zeta_B)\nabla\psi + I(\psi)\nabla\theta_B + G(\psi)\nabla\zeta_B.

We define :math:`\nu` to be the difference between the two toroidal angles:

.. math::
   :label: angle_diff

   \zeta_B = \zeta_0 + \nu.

Plugging this expression into :eq:`straight_field_line_B` and equating the result with :eq:`straight_field_line`
gives :math:`\nabla\psi\times\nabla\left( \theta_B - \theta_0^* - \iota \nu \right) = 0`.
Therefore :math:`\theta_B - \theta_0^* - \iota \nu = y(\psi)` for some flux function :math:`y(\psi)`.
We are free to choose :math:`y(\psi)=0`, which amounts to specifying the origin of :math:`\theta_B` on each surface.
Thus,

.. math::
   :label: theta_diff

   \theta_B = \theta_0 + \lambda + \iota \nu.

.. note:: In Hirshman's note "Transformation from VMEC to Boozer coordinates", the angle difference
	  is defined as :math:`p = \zeta_B - \zeta_0`. However in the fortran ``booz_xform`` code,
	  a minus sign appears on line 83 of ``boozer.f`` that reverses the sign, so the ``p`` quantity
	  saved in ``boozmn_*.nc`` files is in fact :math:`\zeta_0 - \zeta_B`. In the C++/python ``booz_xform``
	  code, the sign convention :math:`\nu = \zeta_B - \zeta_0` is used consistently throughout the
	  code and output files.

Plugging :eq:`angle_diff`-:eq:`theta_diff` into the covariant Boozer representation :eq:`Boozer` gives

.. math::
   :label: B_intermediate

   \vec{B} = G \nabla \zeta_0 + (\nabla \theta_0 + \nabla \lambda) I
   + (G + \iota I) \nabla \nu
     + \left(\beta + I \nu \frac{d \iota}{d\psi}\right) \nabla \psi.

Note that the covariant components of :math:`\vec{B}` in the original coordinates are

.. math::
   :label: covar_components

   B_{\theta_0} = \vec{B} \cdot \frac{\partial\vec{r}}{\partial\theta_0}
   = \vec{B}\cdot\frac{\nabla\zeta_0\times\nabla\psi}{\nabla\psi\cdot\nabla\theta_0\cdot\nabla\zeta_0}, \\
   
   B_{\zeta_0} = \vec{B} \cdot \frac{\partial\vec{r}}{\partial\zeta_0}
   = \vec{B}\cdot\frac{\nabla\psi\times\nabla\theta_0}{\nabla\psi\cdot\nabla\theta_0\cdot\nabla\zeta_0},

where :math:`\vec{r}` is the position vector, and the dual relations have been used to get the
right expressions from the central ones. Plugging :eq:`B_intermediate` into these expressions gives

.. math ::
   :label: find_nu1

   B_{\theta_0} = \left( 1 + \frac{\partial\lambda}{\partial\theta_0}\right) I
   + (G + \iota I) \frac{\partial\nu}{\partial\theta_0}

and

.. math ::
   :label: find_nu2

   B_{\zeta_0} = G + I \frac{\partial\lambda}{\partial\zeta_0}
   + (G + \iota I) \frac{\partial\nu}{\partial\zeta}.

Equations :eq:`find_nu1` and :eq:`find_nu2` determine :math:`\nu` up to a flux function (i.e. a constant on each surface).
From :eq:`angle_diff`, it can be seen that this constant effectively determines the origin of the :math:`\zeta_B` coordinate.
We are free to fix this constant by requiring that the average of :math:`\nu` 
over :math:`\theta_0` and :math:`\zeta_0` be zero.
We can write a double Fourier series for :math:`\nu`,

.. math ::
   :label: nu_Fourier

   \nu(\psi,\theta_0,\zeta_0) = \sum_{m,n} \left[
   \hat{\nu}_{m,n}^s(\psi) \sin(m \theta_0 - n \zeta_0)
   + \hat{\nu}_{m,n}^c(\psi) \cos(m \theta_0 - n \zeta_0) \right],

with analogous series for :math:`\lambda`, :math:`B_{\theta_0}`, and :math:`B_{\zeta_0}`.
The :math:`m=n=0` modes of :eq:`find_nu1`-:eq:`find_nu2` then give

.. math::
   :label: GI

   I = \hat{B}_{\theta_0, 0, 0}^c, \;\;
   G = \hat{B}_{\zeta_0, 0, 0}^c,

allowing the flux functions :math:`I` and :math:`G` to be computed from the input data.
When at least one of :math:`m` or :math:`n` is nonzero, the :math:`\cos(m \theta_0 - n \zeta_0)`
modes of :eq:`find_nu1`-:eq:`find_nu2` are

.. math::
   :label: nus1

   \hat{\nu}_{m,n}^s = \frac{1}{G+\iota I} \left(
   \frac{\hat{B}_{\theta_0, m, n}^c}{m} - I \hat{\lambda}_{m,n}^s \right)

and

.. math::
   :label: nus2

   \hat{\nu}_{m,n}^s = \frac{1}{G+\iota I} \left(
   \frac{-\hat{B}_{\zeta_0, m, n}^c}{n} - I \hat{\lambda}_{m,n}^s \right),

and the :math:`\sin(m \theta_0 - n \zeta_0)` modes of :eq:`find_nu1`-:eq:`find_nu2` are

.. math::
   :label: nuc1

   \hat{\nu}_{m,n}^c = \frac{1}{G+\iota I} \left(
   \frac{-\hat{B}_{\theta_0, m, n}^s}{m} - I \hat{\lambda}_{m,n}^c \right)

and

.. math::
   :label: nuc2

   \hat{\nu}_{m,n}^c = \frac{1}{G+\iota I} \left(
   \frac{\hat{B}_{\zeta_0, m, n}^s}{n} - I \hat{\lambda}_{m,n}^c \right).

Equations :eq:`nus1`-:eq:`nuc2` enable :math:`\nu` to be calculated from the input data.
To see that :eq:`nus1`-:eq:`nus2` are consistent with each other, the curl of
:math:`\vec{B} = B_\psi\nabla\psi + B_{\theta_0}\nabla\theta_0 + B_{\zeta_0}\nabla\zeta_0`
can be plugged into the MHD equilibrium property :math:`\vec{J}\cdot\nabla\psi=0`,
yielding :math:`\partial B_{\zeta_0}/\partial\theta_0 - \partial B_{\theta_0}/\partial\zeta_0=0`.
In Fourier space, then,

.. math::
   :label: radial_current

   m \hat{B}_{\zeta_0,m,n}^s = -n \hat{B}_{\zeta_0,m,n}^s, \;\;
   m \hat{B}_{\zeta_0,m,n}^c = -n \hat{B}_{\zeta_0,m,n}^c.

The modes of :math:`\nu` with :math:`m \ne 0` can be computed from :eq:`nus1` and :eq:`nuc1`.
The modes of :math:`\nu` with :math:`n \ne 0` can be computed from :eq:`nus2` and :eq:`nuc2`.
For modes of :math:`\nu` with both :math:`n \ne 0` and :math:`m \ne 0`, we can use either
:eq:`nus1` or :eq:`nus2`, since the expressions are equivalent by :eq:`radial_current`;
for the same reason we can use either :eq:`nuc1` or :eq:`nuc2`.
The mode of :math:`\nu` with :math:`m=n=0` can be set to zero, for as described above, this
choice amounts to specifying the origin of :math:`\zeta_B`.

The results :eq:`nus1`-:eq:`nuc2` can be summarized as

.. math::
   :label: w

   \nu = \frac{w - I \lambda}{G + \iota I}

where :math:`w` has a Fourier series analogous to :eq:`nu_Fourier`,
with :math:`\hat{w}_{0,0}^c=0` (assuming :math:`\lambda_{0,0}^c=0`) and

.. math::
   :label: w_modes
	   
   \hat{w}_{m,n}^s = \hat{B}_{\theta_0,m,n}^c / m = -\hat{B}_{\zeta_0,m,n}^c / n, \\
   \hat{w}_{m,n}^c = -\hat{B}_{\theta_0,m,n}^s / m = \hat{B}_{\zeta_0,m,n}^s / n.

   
2. Transforming other quantities
--------------------------------

Now that :math:`\nu` is determined, the remaining task is to express
other scalar quantities like :math:`B` as functions of the new angles instead of the old angles.
We let :math:`\Omega` denote any scalar quantity that we wish to transform in this way,
including the cylindrical coordinates :math:`R` and :math:`Z`, and :math:`\nu` itself. We write
a double Fourier expansion, using a bar instead of hat to denote that the independent variables
are now the Boozer angles:

.. math::
   :label: Fourier_Boozer

   \Omega(\psi,\theta_B,\zeta_B)
   = \sum_{m,n} \left[
   \bar{\Omega}_{m,n}^s(\psi) \sin(m \theta_B - n \zeta_B)
   + \bar{\Omega}_{m,n}^c(\psi) \cos(m \theta_B - n \zeta_B) \right].

The amplitudes satisfy

.. math::
   :label: Fourier_Boozer2

   \bar{\Omega}_{m,n}^s = \frac{1}{4\pi^2} \int_0^{2\pi}d\theta_B
   \int_0^{2\pi}d\zeta_B \; \Omega \sin(m \theta_B - n \zeta_B), \\
   \bar{\Omega}_{m,n}^c = \frac{1}{4\pi^2} \int_0^{2\pi}d\theta_B
   \int_0^{2\pi}d\zeta_B \; \Omega \cos(m \theta_B - n \zeta_B).

Changing the integration variables to the old angles,

.. math::
   :label: change_variables

   \bar{\Omega}_{m,n}^s = \frac{1}{4\pi^2} \int_0^{2\pi}d\theta_0
   \int_0^{2\pi}d\zeta_0 \frac{\partial(\theta_B,\zeta_B)}{\partial(\theta_0,\zeta_0)}
   \Omega \sin( m[\theta_0 + \lambda + \iota \nu] - n[\zeta_0 + \nu]), \\
   \bar{\Omega}_{m,n}^c = \frac{1}{4\pi^2} \int_0^{2\pi}d\theta_0
   \int_0^{2\pi}d\zeta_0 \frac{\partial(\theta_B,\zeta_B)}{\partial(\theta_0,\zeta_0)}
   \Omega \cos( m[\theta_0 + \lambda + \iota \nu] - n[\zeta_0 + \nu]),

where :eq:`angle_diff`-:eq:`theta_diff` have been applied. The Jacobian appearing in
:eq:`change_variables` is

.. math::
   :label: Jacobian

   \frac{\partial(\theta_B,\zeta_B)}{\partial(\theta_0,\zeta_0)}
   &=\frac{\partial\theta_B}{\partial\theta_0} \frac{\partial\zeta_B}{\partial\zeta_0}
   -\frac{\partial\theta_B}{\partial\zeta_0} \frac{\partial\zeta_B}{\partial\theta_0} \\
   &=\left(1 + \frac{\partial\lambda}{\partial\theta_0}\right)
   \left(1 + \frac{\partial\nu}{\partial\zeta_0}\right)
   +\left(\iota - \frac{\partial\lambda}{\partial\zeta_0}\right)
   \frac{\partial\nu}{\partial\theta_0}.

It can be seen that :eq:`change_variables` and the last expression in :eq:`Jacobian`
express the new Fourier amplitudes we seek in terms of the old angles
and known quantities on these angles.




Implementation details
**********************

Consistent with the discussion following :eq:`radial_current`, the modes
of :math:`\nu` are computed via :math:`w` as follows. If :math:`m \ne 0`,
:eq:`nus1` and :eq:`nuc1` are used, equivalent to the first equality
on each line of :eq:`w_modes`.
If :math:`m = 0` but :math:`n \ne 0`,
:eq:`nus2` and :eq:`nuc2` are used, equivalent to the last expression
on each line of :eq:`w_modes`. This logic appears in ``surface_solve.cpp``.

The integrals in :eq:`change_variables` are computed using a uniform tensor-product
grid in :math:`\theta_0` and :math:`\zeta_0`.
