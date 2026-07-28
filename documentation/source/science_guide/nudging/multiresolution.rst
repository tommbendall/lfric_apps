.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _nudging_science_multiresolution:

Multi-resolution Nudging
==========================

Nudging can optionally be performed on a coarser mesh than the mesh used
for the rest of the model, controlled by the ``coarse_nudging``
configuration option. This can be useful when the reference data used for
nudging has a coarser effective resolution than the model, since it avoids
computing the nudging increment (and, for the convolution method, the
associated stencil operations) at the model's full resolution.

When ``coarse_nudging=.true.``, the nudging scheme proceeds as follows:

1. The field being nudged, :math:`X`, is mapped from the model's mesh onto
   the coarser mesh named by ``nudging_mesh_name``.
2. The nudging increment, :math:`\Delta X`, is calculated on the coarse
   mesh exactly as described in :ref:`nudging_science_formulation`, using
   a reference field that has been prepared directly on the coarse mesh.
3. The resulting coarse-mesh increment is prolongated back onto the
   model's mesh, using linear-order recovery, to give the increment that
   is actually applied to the model's tendency for :math:`X`.

When ``coarse_nudging=.false.``, no coarsening or prolongation takes
place, and the nudging increment is calculated directly on the model's own
mesh.

.. attention::
   Coarse-mesh nudging requires ``use_multires_coupling=.true.`` (in
   ``namelist:formulation``) in addition to ``coarse_nudging=.true.``, so
   that ``nudging_mesh_name`` is defined. See the note in the
   :ref:`nudging_user_index` section of the User Guide for a current
   limitation regarding this requirement.

The external reference fields used to compute :math:`\Delta X` (see
:ref:`nudging_science_vertical_treatment`) are themselves prepared on the
coarse mesh whenever ``coarse_nudging=.true.``: in particular, the Exner
pressure fields used to construct the model pressure profile are first
mapped onto the coarse mesh before being used in the vertical
interpolation.
