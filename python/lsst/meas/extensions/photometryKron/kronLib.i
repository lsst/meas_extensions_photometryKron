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
%module(package="lsst.meas.extensions.photometryKron.kronLib", docstring=kronLib_DOCSTRING) kronLib

%pythonnondynamic;
%naturalvar;  // use const reference typemaps

%include "lsst/p_lsstSwig.i"

%lsst_exceptions()

%{
#include "lsst/pex/policy.h"
#include "lsst/afw/geom.h"
#include "lsst/meas/extensions/detail/KronPhotometry.h"
%}

%import "lsst/afw/detection/detectionLib.i"

SWIG_SHARED_PTR_DERIVED(KronShapePtr,
                        lsst::afw::detection::Photometry, lsst::meas::algorithms::detail::KronPhotometry);

%include "lsst/meas/extensions/detail/KronPhotometry.h"

%inline %{
    PTR(lsst::meas::algorithms::detail::KronPhotometry)
    cast_Kron(PTR(lsst::afw::detection::Photometry) photom) {
        return boost::shared_dynamic_cast<lsst::meas::algorithms::detail::KronPhotometry>(photom);
    }
%}

%extend lsst::meas::algorithms::detail::KronPhotometry {
    %pythoncode {
    def getKronRadius(self):
        return self.getParameter()
    }
}
