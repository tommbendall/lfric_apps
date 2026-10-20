import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


class vn31_tNNN(MacroUpgrade):
    """Upgrade macro for ticket #NNN by Thomas Bendall."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_tNNN"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        """Replace horizontal_transport_predictor with wind_transport_predictor."""
        old = ["namelist:formulation", "horizontal_transport_predictor"]
        new = ["namelist:formulation", "wind_transport_predictor"]
        value = self.get_setting_value(config, old)
        self.remove_setting(config, old)
        if value == ".true.":
            self.add_setting(config, new, "'no_etadot'")
        else:
            self.add_setting(config, new, "'3d'")
        return config, self.reports


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


class vn31_t238(MacroUpgrade):
    """Upgrade macro for ticket #238 by Thomas Bendall."""

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t238"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-driver
        self.add_setting(
            config, ["namelist:finite_element", "coord_space"], "'Wchi'"
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )
        self.add_setting(
            config,
            ["namelist:finite_element", "coord_order_multigrid"],
            coord_order,
        )

        return config, self.reports


class vn31_t180(MacroUpgrade):
    """Upgrade macro for ticket #180 by Thomas Bendall."""

    BEFORE_TAG = "vn3.1_t238"
    AFTER_TAG = "vn3.1_t180"

    def upgrade(self, config, meta_config=None):
        # Get values
        n_orog_smooth = self.get_setting_value(
            config, ["namelist:initialization", "n_orog_smooth"]
        )
        w0_mapping = self.get_setting_value(
            config, ["namelist:initialization", "w0_orography_mapping"]
        )
        coord_order = self.get_setting_value(
            config, ["namelist:finite_element", "coord_order"]
        )

        # Add new settings
        self.add_setting(
            config, ["namelist:orography", "n_orog_smooth"], n_orog_smooth
        )
        self.add_setting(
            config, ["namelist:orography", "w0_multigrid_mapping"], w0_mapping
        )
        self.add_setting(
            config, ["namelist:orography", "orography_order"], coord_order
        )

        # Remove old settings
        self.remove_setting(
            config, ["namelist:initialization", "n_orog_smooth"]
        )
        self.remove_setting(
            config, ["namelist:initialization", "w0_orography_mapping"]
        )

        return config, self.reports
class vn31_t111(MacroUpgrade):
    # Upgrade macro for #111 by Ian Boutle

    BEFORE_TAG = "vn3.1_t180"
    AFTER_TAG = "vn3.1_t111"

    def upgrade(self, config, meta_config=None):
        # Add settings
        nml="namelist:boundaries"
        self.add_setting(config,[nml,"lbc_bal_meth"], "'keep_rho'")
        self.add_setting(config,[nml,"lbc_sort_theta"], ".true.")
        nml="namelist:initialization"
        eos_height=self.get_setting_value(config,[nml, "model_eos_height"])
        self.remove_setting(config,[nml,"model_eos_height"])
        self.add_setting(config,[nml,"init_eos_height"], eos_height)
        self.add_setting(config,[nml,"init_exner_method"], "'hydrostatic'")
        self.add_setting(config,[nml,"init_sort_theta"], ".true.")
        return config, self.reports
