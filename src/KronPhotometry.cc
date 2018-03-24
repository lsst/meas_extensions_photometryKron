// -*- LSST-C++ -*-
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

#include <numeric>
#include <cmath>
#include <functional>
#include "boost/math/constants/constants.hpp"
#include "lsst/utils/Cache.h"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/SpanSet.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/math/KernelFunctions.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/base.h"
#include "lsst/meas/base/ApertureFlux.h"

#include "lsst/meas/extensions/photometryKron.h"

namespace lsst {
namespace meas {
namespace extensions {
namespace photometryKron {

namespace {
base::FlagDefinitionList flagDefinitions;
} // end anonymous

base::FlagDefinition const KronFluxAlgorithm::FAILURE = flagDefinitions.addFailureFlag( "general failure flag, set if anything went wrong");
base::FlagDefinition const KronFluxAlgorithm::EDGE = flagDefinitions.add("flag_edge", "bad measurement due to image edge");
base::FlagDefinition const KronFluxAlgorithm::BAD_SHAPE_NO_PSF = flagDefinitions.add("flag_bad_shape_no_psf", "bad shape and no PSF");
base::FlagDefinition const KronFluxAlgorithm::NO_MINIMUM_RADIUS = flagDefinitions.add("flag_no_minimum_radius", "minimum radius could not enforced: no minimum value or PSF");
base::FlagDefinition const KronFluxAlgorithm::NO_FALLBACK_RADIUS = flagDefinitions.add("flag_no_fallback_radius", "no minimum radius and no PSF provided");
base::FlagDefinition const KronFluxAlgorithm::BAD_RADIUS = flagDefinitions.add("flag_bad_radius", "bad Kron radius");
base::FlagDefinition const KronFluxAlgorithm::USED_MINIMUM_RADIUS = flagDefinitions.add("flag_used_minimum_radius", "used the minimum radius for the Kron aperture");
base::FlagDefinition const KronFluxAlgorithm::USED_PSF_RADIUS = flagDefinitions.add("flag_used_psf_radius", "used the PSF Kron radius for the Kron aperture");
base::FlagDefinition const KronFluxAlgorithm::SMALL_RADIUS = flagDefinitions.add("flag_small_radius", "measured Kron radius was smaller than that of the PSF");
base::FlagDefinition const KronFluxAlgorithm::BAD_SHAPE = flagDefinitions.add("flag_bad_shape", "shape for measuring Kron radius is bad; used PSF shape");

base::FlagDefinitionList const & KronFluxAlgorithm::getFlagDefinitions() {
    return flagDefinitions;
}

LSST_EXCEPTION_TYPE(BadKronException, pex::exceptions::RuntimeError,
                    lsst::meas::extensions::photometryKron::BadKronException);

namespace {

template <typename MaskedImageT>
        class FootprintFlux {
public:
    explicit FootprintFlux() : _sum(0.0), _sumVar(0.0) {}

    /// @brief Reset everything for a new Footprint
    void reset() {
        _sum = _sumVar = 0.0;
    }
    void reset(afw::detection::Footprint const&) {}

    /// @brief method called for each pixel by applyFunctor
    void operator()(afw::geom::Point2I const & pos,
                    typename MaskedImageT::Image::Pixel const & ival,
                    typename MaskedImageT::Variance::Pixel const & vval) {
        _sum += ival;
        _sumVar += vval;
    }

    /// Return the Footprint's flux
    double getSum() const { return _sum; }

    /// Return the variance of the Footprint's flux
    double getSumVar() const { return _sumVar; }

private:
    double _sum;
    double _sumVar;
};

/************************************************************************************************************/
///
/// Find the first elliptical moment of an object
///
/// We know the shape and orientation of the ellipse to use, specified by ab and theta
/// If we take theta to be 0 the major axis is aligned along the x axis.  In this case
/// we may define the elliptical radius as sqrt(x^2 + (y*a/b)^2) of a point (x, y).
///
/// In other words, it's the length of the major axis of the ellipse of specified shape that passes through
/// the point
///
template <typename MaskedImageT, typename WeightImageT>
class FootprintFindMoment {
public:
    FootprintFindMoment(MaskedImageT const& mimage, ///< The image the source lives in
                        afw::geom::Point2D const& center, // center of the object
                        double const ab,                // axis ratio
                        double const theta // rotation of ellipse +ve from x axis
        ) : _xcen(center.getX()), _ycen(center.getY()),
                           _ab(ab),
                           _cosTheta(::cos(theta)),
                           _sinTheta(::sin(theta)),
                           _sum(0.0), _sumR(0.0),
#if 0
                           _sumVar(0.0), _sumRVar(0.0),
#endif
                           _imageX0(mimage.getX0()), _imageY0(mimage.getY0())
        {}

    /// @brief Reset everything for a new Footprint
    void reset() {}
    void reset(afw::detection::Footprint const& foot) {
        _sum = _sumR = 0.0;
#if 0
        _sumVar = _sumRVar = 0.0;
#endif

        MaskedImageT const& mimage = this->getImage();
        afw::geom::Box2I const& bbox(foot.getBBox());
        int const x0 = bbox.getMinX(), y0 = bbox.getMinY(), x1 = bbox.getMaxX(), y1 = bbox.getMaxY();

        if (x0 < _imageX0 || y0 < _imageY0 ||
            x1 >= _imageX0 + mimage.getWidth() || y1 >= _imageY0 + mimage.getHeight()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::OutOfRangeError,
                              (boost::format("Footprint %d,%d--%d,%d doesn't fit in image %d,%d--%d,%d")
                               % x0 % y0 % x1 % y1
                               % _imageX0 % _imageY0
                               % (_imageX0 + mimage.getWidth() - 1) % (_imageY0 + mimage.getHeight() - 1)
                              ).str());
        }
    }

    /// @brief method called for each pixel by applyFunctor
    void operator()(afw::geom::Point2I const & pos, typename MaskedImageT::Image::Pixel const & ival) {
        double x = static_cast<double>(pos.getX());
        double y = static_cast<double>(pos.getY());
        double const dx = x - _xcen;
        double const dy = y - _ycen;
        double const du =  dx*_cosTheta + dy*_sinTheta;
        double const dv = -dx*_sinTheta + dy*_cosTheta;

        double r = ::hypot(du, dv*_ab); // ellipsoidal radius
#if 1
        if (::hypot(dx, dy) < 0.5) {    // within a pixel of the centre
            /*
             * We gain significant precision for flattened Gaussians by treating the central pixel specially
             *
             * If the object's centered in the pixel (and has constant surface brightness) we have <r> == eR;
             * if it's at the corner <r> = 2*eR; we interpolate between these exact results linearily in the
             * displacement.  And then add in quadrature which is also a bit dubious
             *
             * We could avoid all these issues by estimating <r> using the same trick as we use for
             * the sinc fluxes; it's not clear that it's worth it.
             */

            double const eR = 0.38259771140356325; // <r> for a single square pixel, about the centre
            r = ::hypot(r, eR*(1 + ::hypot(dx, dy)/afw::geom::ROOT2));
        }
#endif

        _sum += ival;
        _sumR += r*ival;
#if 0
        typename MaskedImageT::Variance::Pixel vval = iloc.variance(0, 0);
        _sumVar += vval;
        _sumRVar += r*r*vval;
#endif
    }

    /// Return the Footprint's <r_elliptical>
    double getIr() const { return _sumR/_sum; }

#if 0
    /// Return the variance of the Footprint's <r>
//    double getIrVar() const { return _sumRVar/_sum - getIr()*getIr(); } // Wrong?
    double getIrVar() const { return _sumRVar/(_sum*_sum) + _sumVar*_sumR*_sumR/::pow(_sum, 4); }
#endif

    /// Return whether the measurement might be trusted
    bool getGood() const { return _sum > 0 && _sumR > 0; }

private:
    double const _xcen;                 // center of object
    double const _ycen;                 // center of object
    double const _ab;                   // axis ratio
    double const _cosTheta, _sinTheta;  // {cos,sin}(angle from x-axis)
    double _sum;                        // sum of I
    double _sumR;                       // sum of R*I
#if 0
    double _sumVar;                     // sum of Var(I)
    double _sumRVar;                    // sum of R*R*Var(I)
#endif
    int const _imageX0, _imageY0;       // origin of image we're measuring

};
} // end anonymous namespace

afw::geom::ellipses::Axes KronAperture::getKronAxes(
    afw::geom::ellipses::Axes const& shape,
    afw::geom::LinearTransform const& transformation,
    double const radius
    )
{
    afw::geom::ellipses::Axes axes(shape);
    axes.scale(radius/axes.getDeterminantRadius());
    return axes.transform(transformation);
}

template<typename ImageT>
PTR(KronAperture) KronAperture::determineRadius(
    ImageT const& image,
    afw::geom::ellipses::Axes axes,
    afw::geom::Point2D const& center,
    KronFluxControl const& ctrl
    )
{
    //
    // We might smooth the image because this is what SExtractor and Pan-STARRS do.  But I don't see much gain
    //
    double const sigma = ctrl.smoothingSigma; // Gaussian width of smoothing sigma to apply
    bool const smoothImage = sigma > 0;
    int kSize = smoothImage ? 2*int(2*sigma) + 1 : 1;
    afw::math::GaussianFunction1<afw::math::Kernel::Pixel> gaussFunc(smoothImage ? sigma : 100);
    afw::math::SeparableKernel kernel(kSize, kSize, gaussFunc, gaussFunc);
    bool const doNormalize = true, doCopyEdge = false;
    afw::math::ConvolutionControl convCtrl(doNormalize, doCopyEdge);
    double radius0 = axes.getDeterminantRadius();
    double radius = std::numeric_limits<double>::quiet_NaN();
    float radiusForRadius = std::nanf("");
    for (int i = 0; i < ctrl.nIterForRadius; ++i) {
        axes.scale(ctrl.nSigmaForRadius);
        radiusForRadius = axes.getDeterminantRadius(); // radius we used to estimate R_K
        //
        // Build an elliptical Footprint of the proper size
        //
        afw::detection::Footprint foot(afw::geom::SpanSet::fromShape(
            afw::geom::ellipses::Ellipse(axes, center)));
        afw::geom::Box2I bbox = !smoothImage ?
            foot.getBBox() :
            kernel.growBBox(foot.getBBox()); // the smallest bbox needed to convolve with Kernel
        bbox.clip(image.getBBox());
        ImageT subImage(image, bbox, afw::image::PARENT, smoothImage);
        if (smoothImage) {
            afw::math::convolve(subImage, ImageT(image, bbox, afw::image::PARENT, false), kernel, convCtrl);
        }
        //
        // Find the desired first moment of the elliptical radius, which corresponds to the major axis.
        //
        FootprintFindMoment<ImageT, afw::detection::Psf::Image> iRFunctor(
            subImage, center, axes.getA()/axes.getB(), axes.getTheta()
        );

        try {
            foot.getSpans()->applyFunctor(
                iRFunctor, *(subImage.getImage()));
        } catch(lsst::pex::exceptions::OutOfRangeError &e) {
            if (i == 0) {
                LSST_EXCEPT_ADD(e, "Determining Kron aperture");
            }
            break;                      // use the radius we have
        }

        if (!iRFunctor.getGood()) {
            throw LSST_EXCEPT(BadKronException, "Bad integral defining Kron radius");
        }

        radius = iRFunctor.getIr()*sqrt(axes.getB()/axes.getA());
        if (radius <= radius0) {
            break;
        }
        radius0 = radius;

        axes.scale(radius/axes.getDeterminantRadius()); // set axes to our current estimate of R_K
        iRFunctor.reset();
    }

    return std::make_shared<KronAperture>(center, axes, radiusForRadius);
}

// Photometer an image with a particular aperture
template<typename ImageT>
std::pair<double, double> photometer(
    ImageT const& image, // Image to measure
    afw::geom::ellipses::Ellipse const& aperture, // Aperture in which to measure
    double const maxSincRadius // largest radius that we use sinc apertures to measure
    )
{
    afw::geom::ellipses::Axes const& axes = aperture.getCore();
    if (axes.getB() > maxSincRadius) {
        FootprintFlux<ImageT> fluxFunctor;
        auto spans = afw::geom::SpanSet::fromShape(aperture);
        spans->applyFunctor(
                fluxFunctor, *(image.getImage()), *(image.getVariance()));
        return std::make_pair(fluxFunctor.getSum(), ::sqrt(fluxFunctor.getSumVar()));
    }
    try {
        base::ApertureFluxResult fluxResult = base::ApertureFluxAlgorithm::computeSincFlux<float>(image, aperture);
        return std::make_pair(fluxResult.flux, fluxResult.fluxSigma);
    } catch(pex::exceptions::LengthError &e) {
        LSST_EXCEPT_ADD(e, (boost::format("Measuring Kron flux for object at (%.3f, %.3f);"
                                          " aperture radius %g,%g theta %g")
                            % aperture.getCenter().getX() % aperture.getCenter().getY()
                            % axes.getA() % axes.getB() % afw::geom::radToDeg(axes.getTheta())).str());
        throw e;
    }
}


double calculatePsfKronRadius(
    CONST_PTR(afw::detection::Psf) const& psf, // PSF to measure
    afw::geom::Point2D const& center, // Centroid of source on parent image
    double smoothingSigma=0.0         // Gaussian sigma of smoothing applied
    )
{
    assert(psf);
    double const radius = psf->computeShape(center).getDeterminantRadius();
    // For a Gaussian N(0, sigma^2), the Kron radius is sqrt(pi/2)*sigma
    return ::sqrt(afw::geom::PI/2)*::hypot(radius, std::max(0.0, smoothingSigma));
}

template<typename ImageT>
std::pair<double, double> KronAperture::measureFlux(
    ImageT const& image,
    double const nRadiusForFlux,
    double const maxSincRadius
    ) const
{
    afw::geom::ellipses::Axes axes(getAxes()); // Copy of ellipse core, so we can scale
    axes.scale(nRadiusForFlux);
    afw::geom::ellipses::Ellipse const ellip(axes, getCenter());

    return photometer(image, ellip, maxSincRadius);
}

/************************************************************************************************************/

/**
 * @brief A class that knows how to calculate fluxes using the KRON photometry algorithm
 *
 * @ingroup meas/algorithms
 */
KronFluxAlgorithm::KronFluxAlgorithm(
    KronFluxControl const & ctrl,
    std::string const & name,
    afw::table::Schema & schema,
    daf::base::PropertySet & metadata
) : _name(name),
    _ctrl(ctrl),
    _fluxResultKey(
        meas::base::FluxResultKey::addFields(schema, name, "flux from Kron Flux algorithm")
    ),
    _radiusKey(schema.addField<float>(name + "_radius", "Kron radius (sqrt(a*b))")),
    _radiusForRadiusKey(schema.addField<float>(name + "_radius_for_radius",
                            "radius used to estimate <radius> (sqrt(a*b))")),
    _psfRadiusKey(schema.addField<float>(name + "_psf_radius", "Radius of PSF")),
    _centroidExtractor(schema, name, true)
{
    _flagHandler = meas::base::FlagHandler::addFields(schema, name, getFlagDefinitions());
    metadata.add(name + "_nRadiusForFlux", ctrl.nRadiusForFlux);
}

void KronFluxAlgorithm::fail(
    afw::table::SourceRecord & measRecord,
    meas::base::MeasurementError * error
) const {
    _flagHandler.handleFailure(measRecord, error);
}

void KronFluxAlgorithm::_applyAperture(
    afw::table::SourceRecord & source,
    afw::image::Exposure<float> const& exposure,
    KronAperture const& aperture
    ) const
{
    double const rad = aperture.getAxes().getDeterminantRadius();
    if (rad < std::numeric_limits<double>::epsilon()) {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            BAD_RADIUS.doc,
            BAD_RADIUS.number
        );
    }

    std::pair<double, double> result;
    try {
        result = aperture.measureFlux(exposure.getMaskedImage(), _ctrl.nRadiusForFlux, _ctrl.maxSincRadius);
    } catch (pex::exceptions::LengthError const& e) {
        // We hit the edge of the image; there's no reasonable fallback or recovery
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            EDGE.doc,
            EDGE.number
        );
    } catch(lsst::pex::exceptions::OutOfRangeError &e) {
	    throw LSST_EXCEPT(
                meas::base::MeasurementError,
                EDGE.doc,
                EDGE.number
            );
    }

    // set the results in the source object
    meas::base::FluxResult fluxResult;
    fluxResult.flux = result.first;
    fluxResult.fluxSigma = result.second;
    source.set(_fluxResultKey, fluxResult);
    source.set(_radiusKey, aperture.getAxes().getDeterminantRadius());
    //
    //  REMINDER:  In the old code, the psfFactor is calculated using getPsfFactor,
    //  and the values set for _fluxCorrectionKeys.  See old meas_algorithms version.
}

void KronFluxAlgorithm::_applyForced(
        afw::table::SourceRecord & source,
        afw::image::Exposure<float> const & exposure,
        afw::geom::Point2D const & center,
        afw::table::SourceRecord const & reference,
        afw::geom::AffineTransform const & refToMeas
    ) const
{
    float const radius = reference.get(reference.getSchema().find<float>(_ctrl.refRadiusName).key);
    KronAperture const aperture(reference, refToMeas, radius);
    _applyAperture(source, exposure, aperture);
    if (exposure.getPsf()) {
        source.set(_psfRadiusKey, calculatePsfKronRadius(exposure.getPsf(), center, _ctrl.smoothingSigma));
    }
}

void KronFluxAlgorithm::measure(
                      afw::table::SourceRecord & source,
                      afw::image::Exposure<float> const& exposure
                     ) const {
    afw::geom::Point2D center = _centroidExtractor(source, _flagHandler);

    // Did we hit a condition that fundamentally prevented measuring the Kron flux?
    // Such conditions include hitting the edge of the image and bad input shape, but not low signal-to-noise.
    bool bad = false;

    afw::image::MaskedImage<float> const& mimage = exposure.getMaskedImage();

    double R_K_psf = -1;
    if (exposure.getPsf()) {
        R_K_psf = calculatePsfKronRadius(exposure.getPsf(), center, _ctrl.smoothingSigma);
    }

    //
    // Get the shape of the desired aperture
    //
    afw::geom::ellipses::Axes axes;
    if (!source.getShapeFlag()) {
        axes = source.getShape();
    } else {
        bad = true;
        if (!exposure.getPsf()) {
            throw LSST_EXCEPT(
                meas::base::MeasurementError,
                BAD_SHAPE_NO_PSF.doc,
                BAD_SHAPE_NO_PSF.number
            );
        }
        axes = exposure.getPsf()->computeShape();
        _flagHandler.setValue(source, BAD_SHAPE.number, true);
    }
    if (_ctrl.useFootprintRadius) {
        afw::geom::ellipses::Axes footprintAxes(source.getFootprint()->getShape());
        // if the Footprint's a disk of radius R we want footRadius == R.
        // As <r^2> = R^2/2 for a disk, we need to scale up by sqrt(2)
        footprintAxes.scale(::sqrt(2));

        double radius0 = axes.getDeterminantRadius();
        double const footRadius = footprintAxes.getDeterminantRadius();

        if (footRadius > radius0*_ctrl.nSigmaForRadius) {
            radius0 = footRadius/_ctrl.nSigmaForRadius; // we'll scale it up by nSigmaForRadius
            axes.scale(radius0/axes.getDeterminantRadius());
        }
    }

    PTR(KronAperture) aperture;
    if (_ctrl.fixed) {
        aperture.reset(new KronAperture(source));
    } else {
        try {
            aperture = KronAperture::determineRadius(mimage, axes, center, _ctrl);
        } catch (pex::exceptions::OutOfRangeError& e) {
            // We hit the edge of the image: no reasonable fallback or recovery possible
            throw LSST_EXCEPT(
                meas::base::MeasurementError,
                EDGE.doc,
                EDGE.number
            );
        } catch (BadKronException& e) {
            // Not setting bad=true because we only failed due to low S/N
            aperture = _fallbackRadius(source, R_K_psf, e);
        } catch(pex::exceptions::Exception& e) {
            bad = true; // There's something fundamental keeping us from measuring the Kron aperture
            aperture = _fallbackRadius(source, R_K_psf, e);
        }
    }

    /*
     * Estimate the minimum acceptable Kron radius as the Kron radius of the PSF or the
     * provided minimum radius
     */

    // Enforce constraints on minimum radius
    double rad = aperture->getAxes().getDeterminantRadius();
    if (_ctrl.enforceMinimumRadius) {
        double newRadius = rad;
        if (_ctrl.minimumRadius > 0.0) {
            if (rad < _ctrl.minimumRadius) {
                newRadius = _ctrl.minimumRadius;
                _flagHandler.setValue(source, USED_MINIMUM_RADIUS.number, true);
            }
        } else if (!exposure.getPsf()) {
            throw LSST_EXCEPT(
                meas::base::MeasurementError,
                NO_MINIMUM_RADIUS.doc,
                NO_MINIMUM_RADIUS.number
            );
        } else if (rad < R_K_psf) {
            newRadius = R_K_psf;
            _flagHandler.setValue(source, USED_PSF_RADIUS.number, true);
        }
        if (newRadius != rad) {
            aperture->getAxes().scale(newRadius/rad);
            _flagHandler.setValue(source, SMALL_RADIUS.number, true); // guilty after all
        }
    }

    _applyAperture(source, exposure, *aperture);
    source.set(_radiusForRadiusKey, aperture->getRadiusForRadius());
    source.set(_psfRadiusKey, R_K_psf);
    if (bad) _flagHandler.setValue(source, FAILURE.number, true);
}

void KronFluxAlgorithm::measureForced(
        afw::table::SourceRecord & measRecord,
        afw::image::Exposure<float> const & exposure,
        afw::table::SourceRecord const & refRecord,
        afw::geom::SkyWcs const & refWcs
    ) const {
    afw::geom::Point2D center = _centroidExtractor(measRecord, _flagHandler);
    auto xytransform = afw::geom::makeWcsPairTransform(refWcs, *exposure.getWcs());
    _applyForced(measRecord, exposure, center, refRecord,
                    linearizeTransform(*xytransform, refRecord.getCentroid())
                );

}


PTR(KronAperture) KronFluxAlgorithm::_fallbackRadius(afw::table::SourceRecord& source, double const R_K_psf,
                                            pex::exceptions::Exception& exc) const
{
    _flagHandler.setValue(source, BAD_RADIUS.number, true);
    double newRadius;
    if (_ctrl.minimumRadius > 0) {
        newRadius = _ctrl.minimumRadius;
        _flagHandler.setValue(source, USED_MINIMUM_RADIUS.number, true);
    } else if (R_K_psf > 0) {
        newRadius = R_K_psf;
        _flagHandler.setValue(source, USED_PSF_RADIUS.number, true);
    } else {
        throw LSST_EXCEPT(
            meas::base::MeasurementError,
            NO_FALLBACK_RADIUS.doc,
            NO_FALLBACK_RADIUS.number
        );
    }
    PTR(KronAperture) aperture(new KronAperture(source));
    aperture->getAxes().scale(newRadius/aperture->getAxes().getDeterminantRadius());
    return aperture;
}


#define INSTANTIATE(TYPE) \
template PTR(KronAperture) KronAperture::determineRadius<afw::image::MaskedImage<TYPE> >( \
    afw::image::MaskedImage<TYPE> const&, \
    afw::geom::ellipses::Axes, \
    afw::geom::Point2D const&, \
    KronFluxControl const& \
    ); \
template std::pair<double, double> KronAperture::measureFlux<afw::image::MaskedImage<TYPE> >( \
    afw::image::MaskedImage<TYPE> const&, \
    double const, \
    double const \
    ) const;

INSTANTIATE(float);

}}}} // namespace lsst::meas::extensions::photometryKron
