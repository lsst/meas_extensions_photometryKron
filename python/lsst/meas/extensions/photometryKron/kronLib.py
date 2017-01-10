from __future__ import absolute_import

from lsst.meas.base import BasePlugin, wrapSimpleAlgorithm
from ._photometryKron import KronFluxAlgorithm, KronFluxControl, KronAperture

__all__ = ["KronFluxAlgorithm", "KronFluxControl", "KronAperture", "KronFluxPlugin", "KronFluxForcedPlugin"]

KronFluxPlugin, KronFluxForcedPlugin = wrapSimpleAlgorithm(
    KronFluxAlgorithm,
    name = "ext_photometryKron_KronFlux",
    Control = KronFluxControl,
    executionOrder = BasePlugin.FLUX_ORDER,
    shouldApCorr = True,
    needsMetadata = True,
)
