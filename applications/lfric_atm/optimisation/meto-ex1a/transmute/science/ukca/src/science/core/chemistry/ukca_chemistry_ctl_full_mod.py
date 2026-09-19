##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Summary
# =======
#
# This transformation introduces a chunking loop around the call to the ASAD
# solver in `ukca_chemistry_ctl_full_mod.F90`.  A single call to the ASAD
# solver is replaced with multiple calls, each of which operates on a "chunk"
# of the full domain. The chunk size is taken at compile time from the
# environment variable UKCA_FULL_CHUNK_SIZE.  If this variable is not set then
# the source code is passed through unmodified.
#
# Chunking is mainly achieved by slicing the arguments to the ASAD solver.
# However, the ASAD solver is dependent not only on its arguments but
# also on global ASAD arrays, several of which are read and written directly by
# UKCA full-domain mode.  When chunking is enabled, any ASAD arrays that were
# originally full-domain sized become chunk sized. To resolve this size change,
# the transformer needs to be told about all of these ASAD arrays via
# the following parameter.
#
#   * asad_vars:              a dict mapping the ASAD arrays (originally
#                             full-domain sized), which are accessed by UKCA
#                             full-domain mode, to their associated ranks
#
# Unfortunately, this parameter cannot be inferred automatically because the
# ASAD arrays are dynamically allocated (and hence we don't know which ones
# were originally full-domain sized at compile time). Any changes to UKCA
# full-domain mode must therefore ensure that asad_vars is updated, if
# necessary.
#
# For each variable in asad_vars, the transformer introduces a
# full-domain-sized counterpart.  Outside the chunking loop, it renames each
# access of a narrow (chunk sized) ASAD array into an access of its wide
# (full-domain sized) counterpart.  Inside the chunking loop, slices of the
# newly introduced wide arrays are copied into narrow ASAD arrays, then the
# ASAD solver is called, and then narrow ASAD arrays are copied back into
# slices of their wide counterparts.
#
# In addition to asad_vars, the transformer uses the following parameters.
#
#   * fulldom_size_name:      name of variable holding full-domain size
#   * asad_call_name:         name of the top-level ASAD solver routine
#
# OpenMP can also be added to the chunked loop by setting the environment
# variable UKCA_FULL_CHUNK_OMP to True. By default it will be turned on
# provided the chunk size is not equal to domain size, i.e. loop of length 1
#
# OpenMP parallelism is then added to the chunking loop using an omp parallell do
# directive to allow for top level parallelism on the ASAD solver. Importantly
# a call to ukca_reallocate_asad_arrays which reallocates the THREADPRIVATE
# arrays. This is done within the parallel region to account for the potential
# for chunk_size to be different between iterations, i.e. smaller last
# iteration. Dynamic scheduling has been selected based upon the number of
# solver iterations varying between chunks depending on the complexity of the
# chemistry.
#
# Example
# =======
#
# Given the program
#
#   subroutine main()
#     integer, parameter :: n = 256
#     integer :: arr1(n)
#     integer :: arr2(n, 2)
#     integer :: asad_arr1(32)
#     arr1(:) = foo
#     asad_arr1(:) = arr1(:) + 1
#     call asad_cdrive(arr1, arr2, n)
#     arr1(:) = arr1(:) + asad_arr1(:)
#     arr1(:) = arr1(:) + 1
#   end subroutine
#
# and the parameters
#
#   fulldom_size_name     = "n"
#   asad_vars             = {"asad_arr1" : 1}
#   asad_call_name        = "asad_cdrive"
#   UKCA_FULL_CHUNK_SIZE  = 32
#
# the following program is produced.
#
#   subroutine main()
#     integer, parameter :: n = 256
#     integer :: arr1(n)
#     integer :: arr2(n, 2)
#     integer :: asad_arr1(32)
#     integer :: full_asad_arr1(n)
#     integer :: chunk_begin, chunk_end, chunk_size
#
#     arr1(:) = foo
#     full_asad_arr1(:) = arr1(:) + 1
#     do chunk_begin = 1, n, 32
#       chunk_end = min(n, chunk_begin+chunk_size-1)
#       chunk_size = 1 + chunk_end - chunk_begin
#       asad_arr1(1:chunk_size) = full_asad_arr1(chunk_begin:chunk_end)
#       call asad_cdrive(arr1(chunk_begin:chunk_end),                       &
#                        arr2(chunk_begin:chunk_end, :),                    &
#                        chunk_size)
#       full_asad_arr1(chunk_begin:chunk_end) = asad_arr1(1:chunk_size)
#     end do
#     arr1(:) = arr1(:) + full_asad_arr1(:)
#     arr1(:) = arr1(:) + 1
#   end subroutine

# Imports
# =======

import logging
import os

from psyclone.psyir.nodes import (
    ArrayReference,
    Assignment,
    BinaryOperation,
    Call,
    IfBlock,
    IntrinsicCall,
    Literal,
    Loop,
    Reference,
    Routine,
    Schedule,
    UnaryOperation,
)
from psyclone.psyir.symbols import (
    ArrayType,
    ContainerSymbol,
    DataSymbol,
    ImportInterface,
    RoutineSymbol,
    ScalarType,
    Symbol,
)
from psyclone.psyir.transformations.reference2arrayrange_trans import (
    Reference2ArrayRangeTrans,
)
from psyclone.transformations import OMPParallelLoopTrans, TransformationError
from psyclone.version import __MAJOR__, __MICRO__, __MINOR__

# Conditonal imports
# ==================

psy_version = (__MAJOR__, __MINOR__, __MICRO__)

# Transformation Parameters
# =========================

# Name of variable holding the full-domain size
fulldom_size_name = "tot_n_pnts"

# ASAD arrays in use (and their ranks)
asad_vars = {"rk":  2, "sph2o": 1, "sphno3": 1, "tnd": 1,
             "y":   2, "za":    1, "dpd":    2, "dpw": 2,
             "prk": 2, "fpsc1": 1, "fpsc2":  1}

# Name of the top-level ASAD call
asad_call_name = "asad_cdrive"

# Name of the routine in which to apply the transformation
routine_name = "ukca_chemistry_ctl_full"

# Source and name of the reallocation routine
asad_realloc_routine_loc = ("ukca_chemistry_ctl_col_mod",
                        "ukca_reallocate_asad_arrays")


# Utility
# ==============
def get_bool_env(var_name: str, default: bool = False) -> bool:
    val = os.getenv(var_name)
    if val is None:
        return default
    return val.strip().lower() in ('1', 'true', 't', 'yes', 'y', 'on')


# Transformation
# ==============

def trans(psyir):
    desired_chunk_size = os.getenv("UKCA_FULL_CHUNK_SIZE")
    if desired_chunk_size is None:
        return
    elif desired_chunk_size == "FULL_DOMAIN":
        # Message to print (via umPrint) when chunking enabled
        message_text = ("UKCA full-domain chunking enabled with " +
                        "a chunk size equal to the size of the full " +
                        "domain")
        # We use None to represent the full-domain chunk size
        desired_chunk_size = None
    else:
        # Message to print (via umPrint) when chunking enabled
        message_text = ("UKCA full-domain chunking enabled with " +
                        "a chunk size of " + desired_chunk_size)
    use_omp = get_bool_env("UKCA_FULL_CHUNK_OMP", True)
    if desired_chunk_size is None and use_omp:
        logging.WARNING(
            "Turning off omp as chunk size is set to full domain size")
        use_omp = False

    # Locate correct routine within which to apply the transformation
    for routine in psyir.walk(Routine):
        if routine.name != "ukca_chemistry_ctl_full":
            continue
        for loop in routine.walk(Loop):
            try:
                # Parallelise the "DO l = 1, dim_ntp" loops
                if match_loop(loop, "l", "dim_ntp"):
                    omp_trans.apply(loop)

                # Parallelise the "DO jspf = 1, jpcspf" loop
                if match_loop(loop, "jspf", "jpcspf"):
                    omp_trans.apply(loop, force=True)

                # Parallelise the 3D chunking loop
                if match_loop(loop, "zi", "model_levels"):
                    # Find all "chunk_" arrays (to be marked as private)
                    privates = set()
                    for sym in loop.get_all_accessed_symbols():
                        if (sym.name.startswith("chunk_") and
                                isinstance(sym, DataSymbol) and
                                isinstance(sym.datatype, ArrayType)):
                            privates.add(sym)

                    # Apply the transformation
                    parent, position = loop.parent, loop.position
                    omp_trans.apply(loop, force=True, collapse=3)

                    # Mark explicitly private variables
                    if psy_version < (3, 3, 0):
                        loop.explicitly_private_symbols.update(privates)
                    else:
                        directive = parent.children[position]
                        directive.explicitly_private_symbols.update(
                            privates)

            except TransformationError as err:
                err_msg = ("ukca_chemistry_ctl_full_mod.py: Error: "
                           "could not apply OMP transformation "
                           f"to loop '{loop.variable.name}': "
                           f"{err.message_text}")
                raise TransformationError(err_msg) from err
