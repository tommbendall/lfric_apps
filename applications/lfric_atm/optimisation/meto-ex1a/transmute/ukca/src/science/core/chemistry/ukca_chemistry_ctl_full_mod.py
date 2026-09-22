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
    DataSymbol,
)
from psyclone.psyir.nodes import (
    Loop,
    Reference,
    Routine
)
from psyclone.transformations import (
    OMPLoopTrans,
    TransformationError
)


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
