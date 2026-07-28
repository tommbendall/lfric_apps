.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_formulation:

Formulation
===========

Let :math:`X` denote a field being nudged (potential temperature
:math:`\theta`, or one of the horizontal wind components :math:`u`, :math:`v`),
and let :math:`X_{ref}` denote the corresponding reference field, having
already been vertically interpolated onto the model's levels (see
:ref:`nudging_science_vertical_treatment`). At each call to the nudging
scheme, an increment :math:`\Delta X` is calculated and added to the
model's tendency for :math:`X`.

Two methods are available for calculating :math:`\Delta X`, selected via the
``nudging_method`` configuration option.

.. _nudging_science_formulation_newtonian:

Newtonian relaxation
---------------------

When ``nudging_method='newtonian'``, the increment is simply proportional to
the local difference between the field and its reference value:

.. math:: :label: eq:nudging_newtonian

   \Delta X = w \left( X_{ref} - X \right)

where :math:`w` is a weight described in :ref:`nudging_science_weights`
below. This relaxes :math:`X` directly towards :math:`X_{ref}` at every
horizontal scale present in the reference data, with no spectral
selectivity.

.. _nudging_science_formulation_convolution:

Convolution
-----------

When ``nudging_method='convolution'``, the local difference between the
field and its reference value is first filtered using a convolution in
physical space, before being scaled by the same weight :math:`w`:

.. math:: :label: eq:nudging_convolution

   \Delta X = w \, \left( K * \left( X_{ref} - X \right) \right)

where :math:`K * (\cdot)` denotes convolution with a kernel :math:`K`,
computed over a stencil of cells of extent ``spectral_stencil_extent``
around each grid point (set wide enough to resolve the lowest retained
wavenumber). The kernel is constructed to approximate an ideal spherical
low-pass, or band-pass, filter retaining only total wavenumbers between
``spectral_kmin`` and ``spectral_kmax``:

.. math:: :label: eq:nudging_kernel

   K(\gamma) \propto \sum_{l=k_{min}}^{k_{max}} \frac{2l+1}{4\pi} P_l(\cos \gamma)
   \, \exp\left( -\frac{1}{2}\left(\frac{\gamma}{\sigma}\right)^2 \right)

where :math:`\gamma` is the great-circle (central) angle between a grid
point and a neighbouring point within the stencil, :math:`P_l` is the
Legendre polynomial of degree :math:`l`, and
:math:`\sigma = 2\pi / (1 + k_{max}/3)` is the width of a Gaussian envelope
applied to the sum. The envelope ensures the kernel decays smoothly to zero
at the edge of the stencil, avoiding sharp truncation which would otherwise
spuriously amplify some retained wavenumbers. The weights are normalised so
that they sum to one over the stencil, ensuring the convolution preserves
the mean of the filtered field.

This convolution acts only in the horizontal (it is a 2D operation applied
independently on each model level), and is only implemented for
lowest-order finite elements.

.. _nudging_science_weights:

Nudging weights
----------------

In both of the above methods, the increment is scaled by a single weight
:math:`w`, which combines a vertical tapering factor, a spin-up ramp, and a
relaxation-timescale factor:

.. math:: :label: eq:nudging_weight

   w = w_{taper} \; r_{spinup} \; \min\left(1, \frac{\Delta t}{\tau}\right)

where:

* :math:`w_{taper}` is a height-dependent factor, computed as described in
  :ref:`nudging_science_vertical_treatment`, that ramps the nudging on and
  off over a chosen range of vertical levels;
* :math:`r_{spinup}` is a factor that ramps linearly from 0 to 1 over the
  configured spin-up period, described below;
* :math:`\Delta t` is the model timestep and :math:`\tau` is the
  configured relaxation timescale, ``nudging_relax_time``. If
  ``nudging_relax_time`` is set to zero, this factor is instead set to 1
  (the reference state is relaxed towards fully within a single timestep).

Spin-up ramp
~~~~~~~~~~~~

To avoid introducing a shock into the model when nudging is first applied,
the nudging weight is ramped up gradually over a configurable spin-up
period, :math:`[t_{start}, t_{end}]`, set via ``nudging_spinup_start`` and
``nudging_spinup_end``:

.. math:: :label: eq:nudging_spinup

   r_{spinup} =
   \begin{cases}
     0 & t \le t_{start} \\
     \dfrac{t - t_{start}}{t_{end} - t_{start}} & t_{start} < t < t_{end} \\
     1 & t \ge t_{end}
   \end{cases}

where :math:`t` is the elapsed model time since the start of the run. If
``nudging_spinup_start`` and ``nudging_spinup_end`` are equal, the ramp is
skipped and :math:`r_{spinup} = 1` once :math:`t` exceeds this value.
