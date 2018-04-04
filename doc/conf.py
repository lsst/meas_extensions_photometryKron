"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.meas.extensions.photometryKron


_g = globals()
_g.update(build_package_configs(
    project_name='meas_extensions_photometryKron',
    version=lsst.meas.extensions.photometryKron.version.__version__))
