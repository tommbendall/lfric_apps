import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""

class vn32_t500(MacroUpgrade):
    """Upgrade macro for ticket #500 by Thomas Bendall."""

    BEFORE_TAG = "vn3.2"
    AFTER_TAG = "vn3.2_t500"

    def upgrade(self, config, meta_config=None):
        # Add the per-field swift_peregrin transport option. Default it to
        # .false. for every profile so that existing behaviour is preserved.
        if config.get(["namelist:transport"]) is not None:
            profile_size = self.get_setting_value(
                config, ["namelist:transport", "profile_size"]
            )
            self.add_setting(
                config,
                ["namelist:transport", "swift_peregrin"],
                profile_size + "*.false.",
            )

        return config, self.reports
