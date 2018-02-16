// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2015 LSST Corporation.
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
#ifndef LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
#define LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H

#include <memory>
#include <cmath>

#include "lsst/pex/config.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/geom/SkyWcs.h"
#include "lsst/meas/base/Algorithm.h"
#include "lsst/meas/base/FluxUtilities.h"
#include "lsst/meas/base/CentroidUtilities.h"
#include "lsst/meas/base/FlagHandler.h"
#include "lsst/meas/base/InputUtilities.h"

namespace lsst { namespace meas { namespace extensions { namespace photometryKron {

class KronAperture;

/**
 *  @brief C++ control object for Kron flux.
 *
 *  @sa KronFluxConfig.
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
                       "Minimum Kron radius (if == 0.0 use PSF's Kron radius) if enforceMinimumRadius. "
                       "Also functions as fallback aperture radius if set.");
    LSST_CONTROL_FIELD(enforceMinimumRadius, bool, "If true check that the Kron radius exceeds some minimum");
    LSST_CONTROL_FIELD(useFootprintRadius, bool,
                       "Use the Footprint size as part of initial estimate of Kron radius");
    LSST_CONTROL_FIELD(smoothingSigma, double,
                       "Smooth image with N(0, smoothingSigma^2) Gaussian while estimating R_K");
    LSST_CONTROL_FIELD(refRadiusName, std::string,
                       "Name of field specifying reference Kron radius for forced measurement");

    KronFluxControl() :
        fixed(false),
        nSigmaForRadius(6.0),
        nIterForRadius(1),
        nRadiusForFlux(2.5),
        maxSincRadius(10.0),
        minimumRadius(0.0),
        enforceMinimumRadius(true),
        useFootprintRadius(false),
        smoothingSigma(-1.0),
        refRadiusName("ext_photometryKron_KronFlux_radius")
    {}
};

/**
 *  @brief A measurement algorithm that estimates flux using Kron photometry
 */
class KronFluxAlgorithm : public base::SimpleAlgorithm {
public:

    // Structures and routines to manage flaghandler
    static meas::base::FlagDefinitionList const & getFlagDefinitions();
    static meas::base::FlagDefinition const FAILURE;
    static meas::base::FlagDefinition const EDGE;
    static meas::base::FlagDefinition const BAD_SHAPE_NO_PSF;
    static meas::base::FlagDefinition const NO_MINIMUM_RADIUS;
    static meas::base::FlagDefinition const NO_FALLBACK_RADIUS;
    static meas::base::FlagDefinition const BAD_RADIUS;
    static meas::base::FlagDefinition const USED_MINIMUM_RADIUS;
    static meas::base::FlagDefinition const USED_PSF_RADIUS;
    static meas::base::FlagDefinition const SMALL_RADIUS;
    static meas::base::FlagDefinition const BAD_SHAPE;

    /// A typedef to the Control object for this algorithm, defined above.
    /// The control object contains the configuration parameters for this algorithm.
    typedef KronFluxControl Control;

    KronFluxAlgorithm(
        Control const & ctrl, 
        std::string const & name, 
        afw::table::Schema & schema, 
        daf::base::PropertySet & metadata
    );

    virtual void measure(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure
    ) const;

    virtual void measureForced(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::geom::SkyWcs const & refWcs
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

    PTR(KronAperture) _fallbackRadius(afw::table::SourceRecord& source, double const R_K_psf,
                                      pex::exceptions::Exception& exc) const;

    std::string _name;
    Control _ctrl;
    meas::base::FluxResultKey _fluxResultKey;
    afw::table::Key<float> _radiusKey;
    afw::table::Key<float> _radiusForRadiusKey;
    afw::table::Key<float> _psfRadiusKey;
    meas::base::FlagHandler _flagHandler;
    meas::base::SafeCentroidExtractor _centroidExtractor;
};

class KronAperture {
public:
    KronAperture(afw::geom::Point2D const& center, afw::geom::ellipses::BaseCore const& core,
                 float radiusForRadius=std::nanf("")) :
        _center(center),
        _axes(core),
        _radiusForRadius(radiusForRadius)
        {}

    explicit KronAperture(afw::table::SourceRecord const& source, float radiusForRadius=std::nanf("")) :
        _center(afw::geom::Point2D(source.getX(), source.getY())),
        _axes(source.getShape()),
        _radiusForRadius(radiusForRadius)
        {}

    KronAperture(afw::table::SourceRecord const& reference, afw::geom::AffineTransform const& refToMeas,
                 double radius, float radiusForRadius=std::nanf("")) :
        _center(refToMeas(reference.getCentroid())),
        _axes(getKronAxes(reference.getShape(), refToMeas.getLinear(), radius)),
        _radiusForRadius(radiusForRadius)
        {}

    /// Accessors
    double getX() const { return _center.getX(); }
    double getY() const { return _center.getY(); }
    float getRadiusForRadius() const { return _radiusForRadius; }

    afw::geom::Point2D const& getCenter() const { return _center; }

    afw::geom::ellipses::Axes & getAxes() { return _axes; }

    afw::geom::ellipses::Axes const& getAxes() const { return _axes; }


    /// Determine the Kron Aperture from an image
    ///
    /// Determines the object Kron aperture, using the shape from source.getShape()
    /// (e.g. SDSS's adaptive moments)
    template<typename ImageT>
    static PTR(KronAperture) determineRadius(
        ImageT const& image,  ///< Image to measure
        afw::geom::ellipses::Axes axes,  ///< Shape of aperture
        afw::geom::Point2D const& center,   ///< Centre of source
        KronFluxControl const& ctrl  ///< control the algorithm
        );

    /// Photometer within the Kron Aperture on an image
    template<typename ImageT>
    std::pair<double, double> measureFlux(
        ImageT const& image,  ///< Image to measure
        double const nRadiusForFlux,  ///< Kron radius multiplier
        double const maxSincRadius  ///< largest radius that we use sinc apertyres
        ) const;

    /// Transform a Kron Aperture to a different frame
    PTR(KronAperture) transform(afw::geom::AffineTransform const& trans) const {
        afw::geom::Point2D const center = trans(getCenter());
        afw::geom::ellipses::Axes const axes(*getAxes().transform(trans.getLinear()).copy());
        return std::make_shared<KronAperture>(center, axes);
    }

    /// Determine Kron axes from a reference image
    static
    afw::geom::ellipses::Axes getKronAxes(
        afw::geom::ellipses::Axes const& shape,
        afw::geom::LinearTransform const& transformation,
        double const radius
        );

private:
    afw::geom::Point2D const _center;     // Center of aperture
    afw::geom::ellipses::Axes _axes;      // Ellipse defining aperture shape
    float _radiusForRadius;               // Radius used to estimate the Kron radius
};

}}}} // namespace lsst::meas::extensions::photometryKron

#endif // !LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
