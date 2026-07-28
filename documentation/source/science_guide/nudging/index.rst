.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_index:

Spectral Nudging
================

Spectral nudging is a technique for relaxing a free-running simulation
towards an external reference state (typically a reanalysis, or a driving
model in a regional/limited-area configuration), targeted at a chosen range
of horizontal scales and a chosen range of vertical levels. This allows the
large-scale state of the simulation to be kept close to the reference state,
while still permitting the model to generate its own small-scale
variability.

This section describes the formulation of spectral nudging implemented in
LFRic-Apps. For details of how to configure and enable nudging, see the
:ref:`nudging_user_index` section of the User Guide.

.. toctree::
    :maxdepth: 1

    overview
    formulation
    vertical_treatment
    multiresolution
