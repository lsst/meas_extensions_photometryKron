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
/*
 * This module is only needed to allow us to call the constructor that registers KRON measurements
 */

%define kronLib_DOCSTRING
"
Interface to Kron magnitudes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.meas.extensions.photometryKron.kronLib", docstring=kronLib_DOCSTRING) kronLib

#define HAVE_ellipticalFootprint 1

#if HAVE_ellipticalFootprint
%{
#include "lsst/base.h"
#include "lsst/afw/detection/Footprint.h"
#include "lsst/afw/image/Utils.h"
#include "lsst/afw/geom/Point.h"
%}

%include "lsst/p_lsstSwig.i"

%include "lsst/base.h"
%import "lsst/afw/detection/Footprint.h"
%import "lsst/afw/image/Utils.h"
%import "lsst/afw/geom/Point.h"

%inline %{
namespace lsst {
namespace meas {
namespace algorithms {
PTR(lsst::afw::detection::Footprint)
ellipticalFootprint(lsst::afw::geom::Point2I const& center, //!< The center of the circle
                    double a,                               //!< Major axis (pixels)
                    double b,                               //!< Minor axis (pixels)
                    double theta,                           //!< angle of major axis from x-axis; (radians)
                    lsst::afw::geom::Box2I const& region    //!< Bounding box of MaskedImage footprint
                   );
}}}
%}
#endif
