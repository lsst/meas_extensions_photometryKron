/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include <cstdint>

#include "lsst/pex/config/python.h"  // defines LSST_DECLARE_CONTROL_FIELD
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/afw/table/Source.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/extensions/photometryKron.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace extensions {
namespace photometryKron {

namespace {

void declareKronFluxControl(py::module &mod) {
    py::class_<KronFluxControl> cls(mod, "KronFluxControl");

    cls.def(py::init<>());

    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, fixed);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, nSigmaForRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, nIterForRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, nRadiusForFlux);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, maxSincRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, minimumRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, enforceMinimumRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, useFootprintRadius);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, smoothingSigma);
    LSST_DECLARE_CONTROL_FIELD(cls, KronFluxControl, refRadiusName);
}

void declareKronFluxAlgorithm(py::module &mod) {
    py::class_<KronFluxAlgorithm, std::shared_ptr<KronFluxAlgorithm>, base::SimpleAlgorithm> cls(
            mod, "KronFluxAlgorithm");

    cls.def_static("getFlagDefinitions", &KronFluxAlgorithm::getFlagDefinitions,
                   py::return_value_policy::copy);
    cls.attr("FAILURE") = py::cast(KronFluxAlgorithm::FAILURE);
    cls.attr("EDGE") = py::cast(KronFluxAlgorithm::EDGE);
    cls.attr("BAD_SHAPE_NO_PSF") = py::cast(KronFluxAlgorithm::BAD_SHAPE_NO_PSF);
    cls.attr("NO_MINIMUM_RADIUS") = py::cast(KronFluxAlgorithm::NO_MINIMUM_RADIUS);
    cls.attr("NO_FALLBACK_RADIUS") = py::cast(KronFluxAlgorithm::NO_FALLBACK_RADIUS);
    cls.attr("BAD_RADIUS") = py::cast(KronFluxAlgorithm::BAD_RADIUS);
    cls.attr("USED_MINIMUM_RADIUS") = py::cast(KronFluxAlgorithm::USED_MINIMUM_RADIUS);
    cls.attr("USED_PSF_RADIUS") = py::cast(KronFluxAlgorithm::USED_PSF_RADIUS);
    cls.attr("SMALL_RADIUS") = py::cast(KronFluxAlgorithm::SMALL_RADIUS);
    cls.attr("BAD_SHAPE") = py::cast(KronFluxAlgorithm::BAD_SHAPE);

    cls.def(py::init<KronFluxAlgorithm::Control const &, std::string const &, afw::table::Schema &,
                     daf::base::PropertySet &>(),
            "ctrl"_a, "name"_a, "schema"_a, "metadata"_a);

    cls.def("measure", &KronFluxAlgorithm::measure, "measRecord"_a, "exposure"_a);
    cls.def("measureForced", &KronFluxAlgorithm::measureForced, "measRecord"_a, "exposure"_a, "refRecord"_a,
            "refWcs"_a);
    cls.def("fail", &KronFluxAlgorithm::fail, "measRecord"_a, "error"_a = NULL);
}

using PyKronAperture = py::class_<KronAperture>;

/**
 * Wrap templated methods of KronAperture
 *
 * @tparam ImageT  Image-like class, such as lsst::afw::image::MaskedImage<float>
 * @param cls  pybind11 class wrapping KronAperture
 */
template <typename ImageT>
void declareKronApertureTemplatedMethods(PyKronAperture &cls) {
    cls.def_static("determineRadius", &KronAperture::determineRadius<ImageT>, "image"_a, "axes"_a, "center"_a,
                   "ctrl"_a);
    cls.def("measureFlux", &KronAperture::measureFlux<ImageT>, "image"_a, "nRadiusForFlux"_a,
            "maxSincRadius"_a);
}

void declareKronAperture(py::module &mod) {
    PyKronAperture cls(mod, "KronAperture");

    cls.def(py::init<afw::geom::Point2D const &, afw::geom::ellipses::BaseCore const &, float>(), "center"_a,
            "core"_a, "radiusForRadius"_a = std::nanf(""));
    cls.def(py::init<afw::table::SourceRecord const &, float>(), "source"_a,
            "radiusForRadius"_a = std::nanf(""));
    cls.def(py::init<afw::table::SourceRecord const &, afw::geom::AffineTransform const &, double, float>(),
            "reference"_a, "refToMeas"_a, "radius"_a, "radiusForRadius"_a = std::nanf(""));

    cls.def_static("getKronAxes", &KronAperture::getKronAxes, "shape"_a, "transformation"_a, "radius"_a);

    cls.def("getX", &KronAperture::getX);
    cls.def("getY", &KronAperture::getY);
    cls.def("getRadiusForRadius", &KronAperture::getRadiusForRadius);
    cls.def("getCenter", &KronAperture::getCenter);
    cls.def("getAxes", (afw::geom::ellipses::Axes & (KronAperture::*)()) & KronAperture::getAxes,
            py::return_value_policy::reference_internal);
    cls.def("transform", &KronAperture::transform, "trans"_a);

    declareKronApertureTemplatedMethods<afw::image::MaskedImage<float>>(cls);
}

}  // namespace lsst::meas::extensions::<anonymous>

PYBIND11_PLUGIN(_photometryKron) {
    py::module mod("_photometryKron", "Python wrapper for PhotometryKroh.h");

    declareKronFluxControl(mod);
    declareKronFluxAlgorithm(mod);
    declareKronAperture(mod);

    return mod.ptr();
}
}
}
}
}  // namespace lsst::meas::extensions::photometryKron
