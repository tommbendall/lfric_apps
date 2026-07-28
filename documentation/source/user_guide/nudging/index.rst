.. -----------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

.. _nudging_user_index:

Nudging
=======

Nudging is a form of external forcing that relaxes model fields towards a
specified reference state, typically taken from analysis or reanalysis data.
It can be used, for example, to keep a forecast close to observed large-scale
conditions, or to spin up a model state consistent with an external data
source.

A full mathematical description of the nudging scheme is given in the
:ref:`Science Guide <nudging_science_index>`. This page describes how to
configure and use nudging.

Enabling nudging
----------------

Nudging is enabled per-field via the ``namelist:external_forcing`` options
``theta_forcing`` and ``wind_forcing``. Setting either (or both) of these to
``'nudging'`` activates the nudging scheme for potential temperature or the
horizontal wind components respectively:

.. code-block:: fortran

    &external_forcing
    theta_forcing = 'nudging',
    wind_forcing  = 'nudging',
    /

When nudging is enabled for at least one field, the ``namelist:nudging``
namelist becomes active and must be configured as described below.

.. attention::

    In the current implementation, the nudging algorithm always reads
    ``namelist:multires_coupling=nudging_mesh_name``, regardless of whether
    coarse-mesh nudging (``coarse_nudging``) is actually being used. Since the
    ``multires_coupling`` namelist is only present when
    ``namelist:formulation=use_multires_coupling`` is set to ``.true.``,
    this means that ``use_multires_coupling`` must currently be set to
    ``.true.`` and ``nudging_mesh_name`` must be given a value whenever
    nudging is enabled, even if ``coarse_nudging`` is ``.false.`` and
    nudging is intended to run on the same mesh as the dynamical core.
    Without this, the model will fail at run time when it attempts to use
    nudging.

Reference data files
---------------------

The external reference data used for nudging is read from a file specified
in the ``namelist:files`` namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``nudging_directory``
      - Path to the directory containing the nudging reference data file.
    * - ``nudging_filename``
      - Name of the file (within ``nudging_directory``) containing the
        reference temperature, wind and surface pressure fields.

Configuring the nudging scheme
-------------------------------

The behaviour of the nudging scheme is controlled by the
``namelist:nudging`` namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``nudging_method``
      - Selects the nudging algorithm: ``'newtonian'`` for pointwise
        relaxation towards the reference state, or ``'convolution'`` to
        relax only the large scales of the field, using a convolution in
        physical space that emulates a low-pass top-hat filter in spectral
        space (see the :ref:`Science Guide <nudging_science_index>`).
    * - ``nudging_relax_time``
      - Timescale, in seconds, over which the field (or its large-scale
        component) is relaxed towards the reference state. A value of zero
        sets the field to the reference state in a single time step.
    * - ``nudging_spinup_start`` / ``nudging_spinup_end``
      - Times, in seconds since the start of the run, over which the
        strength of nudging is ramped up linearly from zero to full
        strength. This allows nudging to be introduced gradually rather
        than as a step change.
    * - ``nudge_data_levels``
      - Number of vertical levels on which the reference temperature, wind
        and surface pressure data are provided.
    * - ``nudging_level_bottom`` / ``nudging_level_top``
      - Model levels defining the vertical extent over which nudging is
        applied at full strength.
    * - ``nudging_width_bottom`` / ``nudging_width_top``
      - Number of levels over which nudging is tapered linearly to zero,
        below ``nudging_level_bottom`` and above ``nudging_level_top``
        respectively.
    * - ``spectral_kmin`` / ``spectral_kmax``
      - (``nudging_method = 'convolution'`` only) Minimum and maximum
        wavenumbers retained by the top-hat spectral filter.
    * - ``spectral_stencil_extent``
      - (``nudging_method = 'convolution'`` only) Extent of the stencil
        used to construct the physical-space convolution.

Coarse-mesh nudging
--------------------

For efficiency, nudging may optionally be computed on a mesh coarser than
the one used by the dynamical core. This is controlled from the
``namelist:multires_coupling`` namelist:

.. list-table::
    :header-rows: 1
    :widths: 25 75

    * - Option
      - Description
    * - ``coarse_nudging``
      - If ``.true.``, nudging increments are computed on a coarser mesh and
        then mapped back onto the dynamics mesh, rather than being computed
        directly on the dynamics mesh.
    * - ``nudging_mesh_name``
      - Tag-name of the coarser mesh to use when ``coarse_nudging`` is
        ``.true.``. This mesh must also be listed in
        ``multires_coupling_mesh_tags``.
