##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
This transformation introduces OpenMP directives around loops inside
UKCA chemistry full-domain mode.
"""

from psyclone.psyir.symbols import (
    ArrayType,
    ContainerSymbol,
    DataSymbol,
    ImportInterface,
    RoutineSymbol,
    ScalarType,
)
from psyclone.psyir.nodes import (
    BinaryOperation,
    Call,
    IfBlock,
    IntrinsicCall,
    Literal,
    Loop,
    Reference,
    Routine,
    UnaryOperation,
)
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError
)

# Location of the routine that (re)allocates the THREADPRIVATE ASAD arrays
asad_realloc_routine_loc = ("ukca_chemistry_ctl_col_mod",
                             "ukca_reallocate_asad_arrays")


def match_loop(loop: Loop, var_name: str, stop_name: str) -> bool:
    """
    Return true only if loop's variable is named var_name
    and loop's stop expression is a reference named stop_name.
    """
    return (loop.variable.name == var_name and
            isinstance(loop.stop_expr, Reference) and
            loop.stop_expr.name == stop_name)


def trans(psyir):
    """
    Add OpenMP directives to selected loops.
    """
    # All loops are dynamically scheduled
    omp_trans = OMPLoopTrans(omp_directive="paralleldo",
                             omp_schedule="dynamic")

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
                    for ref in loop.walk(Reference):
                        sym = ref.symbol
                        if (sym.name.startswith("chunk_") and
                                isinstance(sym, DataSymbol) and
                                isinstance(sym.datatype, ArrayType)):
                            privates.add(sym)

                    # Ensure the THREADPRIVATE ASAD scratch arrays (sph2o,
                    # za, rk, y, dpd, etc.) are (re)allocated to the chunk
                    # size on every thread before they are used. Relying
                    # solely on the one-off allocation in ukca_iniasad is
                    # not sufficient: that call only guarantees correct
                    # sizing for the threads active in *that* parallel
                    # region, not necessarily every thread that later
                    # executes this one. Without this, a thread whose
                    # THREADPRIVATE copies are unallocated (or wrongly
                    # sized) crashes the first time it reaches this loop.
                    sym_tab = routine.symbol_table
                    asad_realloc_mod_sym = sym_tab.find_or_create(
                        asad_realloc_routine_loc[0],
                        symbol_type=ContainerSymbol)
                    asad_realloc_routine = sym_tab.find_or_create(
                        asad_realloc_routine_loc[1],
                        symbol_type=RoutineSymbol,
                        interface=ImportInterface(asad_realloc_mod_sym))

                    chunk_n_pnts_sym = sym_tab.lookup("chunk_n_pnts")
                    sph2o_sym = sym_tab.lookup("sph2o")

                    realloc_call = Call.create(
                        asad_realloc_routine,
                        [Reference(chunk_n_pnts_sym)])

                    realloc_block = IfBlock.create(
                        BinaryOperation.create(
                            BinaryOperation.Operator.OR,
                            UnaryOperation.create(
                                UnaryOperation.Operator.NOT,
                                IntrinsicCall.create(
                                    IntrinsicCall.Intrinsic.ALLOCATED,
                                    [Reference(sph2o_sym)])),
                            BinaryOperation.create(
                                BinaryOperation.Operator.NE,
                                Reference(chunk_n_pnts_sym),
                                IntrinsicCall.create(
                                    IntrinsicCall.Intrinsic.SIZE,
                                    [Reference(sph2o_sym),
                                     ("dim",
                                      Literal("1",
                                              ScalarType.integer_type()))]))),
                        [realloc_call])

                    innermost_loop = loop.walk(Loop)[-1]
                    innermost_loop.loop_body.addchild(realloc_block, index=0)

                    # Apply the transformation
                    parent, position = loop.parent, loop.position
                    omp_trans.apply(loop, force=True, collapse=3)

                    # Mark explicitly private variables. PSyclone 3.3.0
                    # moved this attribute from the Loop node to the
                    # enclosing Directive node.
                    directive = parent.children[position]
                    directive.explicitly_private_symbols.update(privates)

            except TransformationError as err:
                err_msg = ("ukca_chemistry_ctl_full_mod.py: Error: "
                           "could not apply OMP transformation "
                           f"to loop '{loop.variable.name}': "
                           f"{err}")
                raise TransformationError(err_msg) from err
