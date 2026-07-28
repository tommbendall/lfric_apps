.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_vertical_treatment:

Vertical Treatment
===================

The external reference data used for nudging is supplied on its own set of
hybrid-pressure levels, which do not in general coincide with the model's
levels. Before the reference fields can be used in the nudging increment,
they must therefore be interpolated onto the model's levels, and the
temperature field converted to potential temperature. A height-dependent
weight is also computed, to allow nudging to be restricted to a chosen
range of vertical levels.

Reference pressure profile
---------------------------

The reference dataset's hybrid-pressure levels are defined by coefficients
:math:`a_k`, :math:`b_k`, supplied for either 137 or 88 levels
(``nudge_data_levels``). Given the reference surface pressure,
:math:`p_{s,ref}`, the pressure on reference half-level :math:`k` is

.. math:: :label: eq:nudging_ref_pressure

   p_{ref,k} = a_k + b_k \, p_{s,ref} ,

and full-level values are formed by averaging :math:`a_k` and :math:`b_k`
between adjacent half levels.

Model pressure profile
------------------------

The corresponding pressure on the model's levels is obtained from the
model's Exner pressure, :math:`\Pi`, on both the ``W3`` and ``Wtheta``
function spaces:

.. math:: :label: eq:nudging_model_pressure

   p = p_0 \, \Pi^{1/\kappa}

where :math:`p_0` is the reference surface pressure (``p_zero``) and
:math:`\kappa = R_d / c_p` (``kappa``).

Temperature to potential temperature
--------------------------------------

The reference dataset supplies temperature, :math:`T_{ref}`, rather than
potential temperature. This is converted using the reference pressure
profile computed above:

.. math:: :label: eq:nudging_T_to_theta

   \theta_{ref} = T_{ref} \left( \frac{p_0}{p_{ref}} \right)^{\kappa} .

Vertical interpolation
------------------------

The reference fields (:math:`\theta_{ref}`, :math:`u_{ref}`,
:math:`v_{ref}`) are interpolated from the reference pressure levels onto
the model's pressure levels by linear interpolation in log-pressure:

.. math:: :label: eq:nudging_log_p_interp

   X_{ref}(p) = X_{ref,l} + \left( X_{ref,l+1} - X_{ref,l} \right)
   \frac{\ln p - \ln p_{ref,l}}{\ln p_{ref,l+1} - \ln p_{ref,l}}

where :math:`l` is the reference level immediately below the model level
being interpolated to (in pressure). To avoid unbounded extrapolation
above the top or below the bottom of the reference data, the interpolation
factor is limited to the range :math:`[-1, 2]`. The lowest model
:math:`\theta` level is set equal to the value interpolated for the level
immediately above it.

.. _nudging_science_vertical_taper:

Vertical tapering
-------------------

To allow nudging to be targeted at a chosen range of vertical levels, a
height-dependent weight, :math:`w_{taper}` (introduced in
:eq:`eq:nudging_weight`), is computed separately for the ``W3``
and ``Wtheta`` function spaces. This weight ramps linearly from 0 to 1
across a configurable band of levels at the bottom of the nudging region,
remains at 1 throughout the region, and ramps back down to 0 across a
configurable band at the top:

.. math:: :label: eq:nudging_vertical_taper

   w_{taper}(k) =
   \begin{cases}
     0 &
       k < k_{bot} - w_{bot} \\
     \dfrac{k - (k_{bot} - w_{bot})}{2 \, w_{bot}} &
       k_{bot} - w_{bot} \le k \le k_{bot} + w_{bot} \\
     1 &
       k_{bot} + w_{bot} < k < k_{top} - w_{top} \\
     1 - \dfrac{k - (k_{top} - w_{top})}{2 \, w_{top}} &
       k_{top} - w_{top} \le k \le k_{top} + w_{top} \\
     0 &
       k > k_{top} + w_{top}
   \end{cases}

where :math:`k` is the model-level index, :math:`k_{bot}` and
:math:`k_{top}` are the ``nudging_level_bottom`` and ``nudging_level_top``
configuration options, and :math:`w_{bot}` and :math:`w_{top}` are the
``nudging_width_bottom`` and ``nudging_width_top`` options, all specified
in terms of model levels.
