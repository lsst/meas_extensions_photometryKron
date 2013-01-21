// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "boost/math/constants/constants.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/image/Exposure.h"
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/ellipses.h"
#include "lsst/meas/algorithms/Photometry.h"

#include "lsst/meas/extensions/photometryKron.h"

namespace lsst {
namespace afwDet = afw::detection;
namespace afwEllipses = afw::geom::ellipses;

namespace meas {
namespace extensions {
namespace photometryKron {
namespace {

/**
 * @brief A class that knows how to calculate fluxes using the KRON photometry algorithm
 *
 * @ingroup meas/algorithms
 */
class KronFlux : public algorithms::FluxAlgorithm {
public:

    KronFlux(KronFluxControl const & ctrl, afw::table::Schema & schema) :
        algorithms::FluxAlgorithm(
            ctrl, schema,
            "Kron photometry: photometry with aperture set to some multiple of <radius>"
            "determined within some multiple of the source size"
        ),
        _radiusKey(schema.addField<float>(ctrl.name + ".radius", "Kron radius (sqrt(a*b))")),
        _radiusForRadiusKey(schema.addField<float>(ctrl.name + ".radiusForRadius",
                                          "Radius used to estimate <radius> (sqrt(a*b))")),
        _badRadiusKey(schema.addField<afw::table::Flag>(ctrl.name + ".flags.radius", "Bad Kron radius"))
    {}

private:
    
    template <typename PixelT>
    void _apply(
        afw::table::SourceRecord & source,
        afw::image::Exposure<PixelT> const & exposure,
        afw::geom::Point2D const & center
    ) const;

    LSST_MEAS_ALGORITHM_PRIVATE_INTERFACE(KronFlux);

    afw::table::Key<float> _radiusKey;
    afw::table::Key<float> _radiusForRadiusKey;
    afw::table::Key<afw::table::Flag> _badRadiusKey;
};

/************************************************************************************************************/

template <typename MaskedImageT, typename WeightImageT>
class FootprintFindMoment : public afwDet::FootprintFunctor<MaskedImageT> {
public:
    FootprintFindMoment(MaskedImageT const& mimage, ///< The image the source lives in
                        afw::geom::Point2D const& center, // center of the object
                        double const ab,                // axis ratio
                        double const theta // rotation of ellipse +ve from x axis
                       ) : afwDet::FootprintFunctor<MaskedImageT>(mimage),
                           _xcen(center.getX()), _ycen(center.getY()),
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
    void reset(afwDet::Footprint const& foot) {
        _sum = _sumR = 0.0;
#if 0
        _sumVar = _sumRVar = 0.0;
#endif

        MaskedImageT const& mimage = this->getImage();
        afw::geom::Box2I const& bbox(foot.getBBox());
        int const x0 = bbox.getMinX(), y0 = bbox.getMinY(), x1 = bbox.getMaxX(), y1 = bbox.getMaxY();

        if (x0 < _imageX0 || y0 < _imageY0 ||
            x1 >= _imageX0 + mimage.getWidth() || y1 >= _imageY0 + mimage.getHeight()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::OutOfRangeException,
                              (boost::format("Footprint %d,%d--%d,%d doesn't fit in image %d,%d--%d,%d")
                               % x0 % y0 % x1 % y1
                               % _imageX0 % _imageY0
                               % (_imageX0 + mimage.getWidth() - 1) % (_imageY0 + mimage.getHeight() - 1)
                              ).str());
        }
    }
    
    /// @brief method called for each pixel by apply()
    void operator()(typename MaskedImageT::xy_locator iloc, ///< locator pointing at the image pixel
                    int x,                                  ///< column-position of pixel
                    int y                                   ///< row-position of pixel
                   ) {
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

        typename MaskedImageT::Image::Pixel ival = iloc.image(0, 0);
        _sum += ival;
        _sumR += r*ival;
#if 0
        typename MaskedImageT::Variance::Pixel vval = iloc.variance(0, 0);
        _sumVar += vval;
        _sumRVar += r*r*vval;
#endif
    }

    /// Return the Footprint's <r>
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


struct KronAperture {
    KronAperture(afw::geom::Point2D const& center, afw::geom::ellipses::Axes const& ellipse) :
        _x(center.getX()), _y(center.getY()), _ellipse(ellipse) {}
    explicit KronAperture(afw::table::SourceRecord const& source) :
        _x(source.getX()), _y(source.getY()),
        _ellipse(source.getShape()) {}

    /// Accessors
    double getX() const { return _x; }
    double getY() const { return _y; }
    afw::geom::ellipses::Axes getEllipse() const { return _ellipse; }

    /// Determine the Kron Aperture from an image
    template<typename ImageT>
    static PTR(KronAperture) determine(ImageT const& image, afw::table::SourceRecord const& source,
                                       afw::geom::Point2D const& center,
                                       KronFluxControl const& ctrl, float *radiusForRadius
                                      );

    /// Photometer within the Kron Aperture on an image
    template<typename ImageT>
    std::pair<double, double> measure(ImageT const& image, // Image to measure
                                      double nRadiusForFlux // Kron radius multiplier
                                     ) const;

    /// Transform a Kron Aperture to a different frame
    PTR(KronAperture) transform(afw::geom::AffineTransform const& trans) const {
        afw::geom::Point2D const center = trans(afw::geom::Point2D(_x, _y));
        afw::geom::ellipses::Axes const ellipse(_ellipse.transform(trans.getLinear()));
        return boost::make_shared<KronAperture>(center, ellipse);
    }

private:
    double _x, _y;                      // Centre
    afw::geom::ellipses::Axes _ellipse;          // Ellipse defining aperture shape
};

/*
 * Estimate the object Kron aperture, using the shape from source.getShape() (e.g. SDSS's adaptive moments)
 */
template<typename ImageT>
PTR(KronAperture) KronAperture::determine(ImageT const& image, // Image to measure
                                          afw::table::SourceRecord const& source, // Source with measurements
                                          afw::geom::Point2D const& center, // Centre of source
                                          KronFluxControl const& ctrl,      // control the algorithm
                                          float *radiusForRadius            // radius used to estimate radius
                                         )
{
    afw::geom::ellipses::Axes axes(source.getShape());
    afw::geom::ellipses::Axes footprintAxes(source.getFootprint()->getShape());
    footprintAxes.scale(2);             // <r^2> = 1/2 for a disk

    double radius0 = axes.getDeterminantRadius();
    double const footRadius = footprintAxes.getDeterminantRadius();
    if (ctrl.useFootprintRadius && footRadius > radius0*ctrl.nSigmaForRadius) {
        radius0 = footRadius/ctrl.nSigmaForRadius; // we'll scale it up by nSigmaForRadius
        axes.scale(radius0/axes.getDeterminantRadius());
    }
    double radius = std::numeric_limits<double>::quiet_NaN();
    for (int i = 0; i < ctrl.nIterForRadius; ++i) {
        axes.scale(ctrl.nSigmaForRadius);
        *radiusForRadius = axes.getDeterminantRadius(); // radius we used to estimate R_K
        
        FootprintFindMoment<ImageT, afwDet::Psf::Image> iRFunctor(
                                                   image, center, axes.getA()/axes.getB(), axes.getTheta()
                                                                 );
        
        // Build an elliptical Footprint of the proper size
        afwDet::Footprint foot(afw::geom::ellipses::Ellipse(axes, center));

        try {
            iRFunctor.apply(foot);
        } catch(lsst::pex::exceptions::OutOfRangeException &e) {
            if (i == 0) {
                LSST_EXCEPT_ADD(e, "Determining Kron aperture");
            }
            break;                      // use the radius we have
        }
        
        if (!iRFunctor.getGood()) {
            throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                              "Bad footprint when determining Kron aperture");
        }
        
        radius = iRFunctor.getIr();
        if (radius <= radius0) {
            break;
        }
        radius0 = radius;

        axes.scale(radius/axes.getDeterminantRadius()); // set axes to our current estimate of R_K
    }

    return boost::make_shared<KronAperture>(
                 center, afw::geom::ellipses::Axes(radius, radius*axes.getB()/axes.getA(), axes.getTheta()));
}

template<typename ImageT>
std::pair<double, double> KronAperture::measure(ImageT const& image, // Image of interest
                                                double nRadiusForFlux // Kron radius multiplier
    ) const
{
    try {
        double const r2 = nRadiusForFlux*_ellipse.getA(); // outer radius
        double const ellip = 1.0 - _ellipse.getB()/_ellipse.getA();
        return algorithms::photometry::calculateSincApertureFlux(
            image, _x, _y, 0.0, r2, _ellipse.getTheta(), ellip
        );
    } catch(pex::exceptions::LengthErrorException &e) {
        LSST_EXCEPT_ADD(e, (boost::format("Measuring Kron flux for object at (%.3f, %.3f);"
                                          " aperture radius %g,%g theta %g")
                            % _x % _y % _ellipse.getA() % _ellipse.getB() %
                            afw::geom::radToDeg(_ellipse.getTheta())).str());
        throw e;
    }
}

/************************************************************************************************************/

template <typename PixelT>
void KronFlux::_apply(
    afw::table::SourceRecord & source,
    afw::image::Exposure<PixelT> const& exposure,
    afw::geom::Point2D const & center
) const {
    source.set(getKeys().flag, true); // bad unless we get all the way to success at the end
    afw::image::MaskedImage<PixelT> const& mimage = exposure.getMaskedImage();

    KronFluxControl const & ctrl = static_cast<KronFluxControl const &>(this->getControl());

    if (source.getShapeFlag()) {        // the shape's bad; give up now
        return;
    }
    
    CONST_PTR(KronAperture) aperture;
    float radiusForRadius = std::numeric_limits<double>::quiet_NaN();
    if (ctrl.fixed) {
        aperture.reset(new KronAperture(source));
    } else {
        try {
            aperture = KronAperture::determine(mimage, source, center,
                                               ctrl.nSigmaForRadius, ctrl.nIterForRadius, &radiusForRadius);
        } catch(pex::exceptions::Exception& e) {
            return;
        }
    }
    source.set(_badRadiusKey, false);

    std::pair<double, double> result = aperture->measure(mimage, ctrl.nRadiusForFlux);
    source.set(getKeys().meas, result.first);
    source.set(getKeys().err, result.second);
    source.set(_radiusKey, aperture->getEllipse().getDeterminantRadius());
    source.set(_radiusForRadiusKey, radiusForRadius);
    source.set(getKeys().flag, false);
}

LSST_MEAS_ALGORITHM_PRIVATE_IMPLEMENTATION(KronFlux);

} // anonymous namespace

PTR(algorithms::AlgorithmControl) KronFluxControl::_clone() const {
    return boost::make_shared<KronFluxControl>(*this);
}

PTR(algorithms::Algorithm) KronFluxControl::_makeAlgorithm(
    afw::table::Schema & schema,
    PTR(daf::base::PropertyList) const &
) const {
    return boost::make_shared<KronFlux>(*this, boost::ref(schema));
}


}}}} // namespace lsst::meas::extensions::photometryKron
