#
# LSST Data Management System
# Copyright 2008-2016 AURA/LSST.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

from lsst.meas.base import BasePlugin, wrapSimpleAlgorithm
from .photometryKron import KronFluxAlgorithm, KronFluxControl, KronAperture

__all__ = ["KronFluxAlgorithm", "KronFluxControl", "KronAperture", "KronFluxPlugin", "KronFluxForcedPlugin"]

KronFluxPlugin, KronFluxForcedPlugin = wrapSimpleAlgorithm(
    KronFluxAlgorithm,
    name = "ext_photometryKron_KronFlux",
    Control = KronFluxControl,
    executionOrder = BasePlugin.FLUX_ORDER,
    shouldApCorr = True,
    needsMetadata = True,
)

from .version import *
