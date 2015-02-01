// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_BASE_KronFlux_h_INCLUDED
#define LSST_MEAS_BASE_KronFlux_h_INCLUDED

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace extensions { namespace photometryKron {
namespace {
    struct KronAperture;
}

/**
 *  @brief A C++ control class to handle KronFluxAlgorithm's configuration
 */
class KronFluxControl {
public:
    LSST_CONTROL_FIELD(fixed, bool,
                       "if true, use existing shape and centroid measurements instead of fitting");
    LSST_CONTROL_FIELD(nSigmaForRadius, double,
                       "Multiplier of rms size for aperture used to initially estimate the Kron radius");
    LSST_CONTROL_FIELD(nIterForRadius, int, "Number of times to iterate when setting the Kron radius");
    LSST_CONTROL_FIELD(nRadiusForFlux, double, "Number of Kron radii for Kron flux");
    LSST_CONTROL_FIELD(maxSincRadius, double,
                       "Largest aperture for which to use the slow, accurate, sinc aperture code");
    LSST_CONTROL_FIELD(minimumRadius, double,
                       "Minimum Kron radius (if == 0.0 use PSF's Kron radius). "
                       "Ignored if enforceMinimumRadius is false");
    LSST_CONTROL_FIELD(enforceMinimumRadius, bool, "If true check that the Kron radius exceeds some minimum");
    LSST_CONTROL_FIELD(useFootprintRadius, bool,
                       "Use the Footprint size as part of initial estimate of Kron radius");
    LSST_CONTROL_FIELD(smoothingSigma, double,
                       "Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K");
    KronFluxControl() :
        fixed(false),
        nSigmaForRadius(6.0),
        nIterForRadius(1),
        nRadiusForFlux(2.5),
        maxSincRadius(10.0),
        minimumRadius(0.0),
        enforceMinimumRadius(true),
        useFootprintRadius(false),
        smoothingSigma(-1.0)
    {}
};
/**
 *  @brief A measurement algorithm that estimates flux using Kron photometry 
 */
class KronFluxAlgorithm : public base::SimpleAlgorithm {
public:

    enum {
        FAILURE=base::FlagHandler::FAILURE,
        RADIUS,
        SMALL_RADIUS, 
        N_FLAGS
    };

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef KronFluxControl Control;

    KronFluxAlgorithm(Control const & ctrl, std::string const & name, afw::table::Schema & schema);

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void fail(
        afw::table::SourceRecord & measRecord,
        meas::base::MeasurementError * error=NULL
    ) const;

private:

    void _applyAperture(
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const& exposure,
        KronAperture const& aperture
        ) const;

    void _applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const;
    std::string const & _name;
    Control _ctrl;
    //algorithms::ScaledFlux::KeyTuple _fluxCorrectionKeys;
    meas::base::FluxResultKey _fluxResultKey;
    afw::table::Key<float> _radiusKey;
    afw::table::Key<float> _radiusForRadiusKey;
    meas::base::FlagHandler _flagHandler;
    meas::base::SafeCentroidExtractor _centroidExtractor;
};

}}}} // namespace lsst::meas::extensions::shapeHSM

#endif // !LSST_MEAS_BASE_KronFlux_h_INCLUDED
