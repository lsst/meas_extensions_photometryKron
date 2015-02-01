// -*- lsst-c++ -*-

/*
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

%define kronLib_DOCSTRING
"
Interface to Kron magnitudes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.photometryKron", docstring=kronLib_DOCSTRING) kronLib

%{
#include "lsst/base.h"
#include "lsst/pex/logging.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/table.h"
#include "lsst/meas/base.h"
#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_BASE_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
#include "ndarray/swig/eigen.h"
%}

%include "ndarray.i"
%declareNumPyConverters(lsst::meas::base::CentroidCov);
%lsst_exceptions();

%include "lsst/p_lsstSwig.i"

%import "lsst/afw/math/mathLib.i"
%import "lsst/afw/image/imageLib.i"
%import "lsst/afw/detection/detectionLib.i"
%import "lsst/afw/table/tableLib.i"



%include "lsst/meas/base/constants.h"
%include "lsst/meas/base/exceptions.i"
%include "lsst/meas/base/utilities.i"
%include "lsst/meas/base/Algorithm.h"
%import "lsst/pex/config.h"

%{
#include "lsst/meas/extensions/photometryKron.h"
%}

%feature("notabstract") lsst::meas::extensions::photometryKron::KronFluxAlgorithm;
%include "lsst/meas/extensions/photometryKron.h"
