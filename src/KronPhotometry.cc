// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "boost/math/constants/constants.hpp"
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Angle.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/afw/detection/AperturePhotometry.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/ellipses.h"


#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDet = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;
namespace afwEllipse = lsst::afw::geom::ellipses;
namespace afwCoord = lsst::afw::coord;

/*
 * Do we have the ellipticalFootprint function?  In sufficiently recent versions of afw it's available
 * as a Footprint constructor;  when this is in a cut version this function should be deleted
 */
#define HAVE_ellipticalFootprint 1      // n.b. if == 0 you can clean up kronLib.i too

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief A class that knows how to calculate fluxes using the KRON photometry algorithm
 *
 * @ingroup meas/algorithms
 */
template<typename ExposureT>
class KronPhotometer : public Algorithm<afwDet::Photometry, ExposureT>
{
public:
    typedef Algorithm<afwDet::Photometry, ExposureT> AlgorithmT;
    typedef boost::shared_ptr<KronPhotometer> Ptr;
    typedef boost::shared_ptr<KronPhotometer const> ConstPtr;

    /// Ctor
    KronPhotometer(bool fixed=false, double nSigmaForRadius=6.0, double nRadiusForFlux=2.0, 
                   double background=0.0, double shiftmax=10.0) :
        AlgorithmT(), _fixed(fixed), _nSigmaForRadius(nSigmaForRadius),
        _nRadiusForFlux(nRadiusForFlux), _background(background), _shiftmax(shiftmax) {}

    virtual std::string getName() const { return "KRON"; }

    virtual PTR(AlgorithmT) clone() const {
        return boost::make_shared<KronPhotometer<ExposureT> >(_nSigmaForRadius, _nRadiusForFlux,
                                                              _background, _shiftmax);
    }

    virtual void configure(pexPolicy::Policy const& policy) {
        if (policy.isBool("fixed")) {
            _fixed = policy.getBool("fixed");
        }
        if (policy.isDouble("nSigmaForRadius")) {
            _nSigmaForRadius = policy.getDouble("nSigmaForRadius");
        }
        if (policy.isDouble("nRadiusForFlux")) {
            _nRadiusForFlux = policy.getDouble("nRadiusForFlux");
        }
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        }
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        }
    }

    virtual PTR(afwDet::Photometry) measureNull(void) const {
        double const NaN = std::numeric_limits<double>::quiet_NaN();
        return boost::make_shared<afwDet::AperturePhotometry>(NaN, NaN, NaN);
    }

    virtual PTR(afwDet::Photometry) measureSingle(afwDet::Source const&, afwDet::Source const&,
                                                  ExposurePatch<ExposureT> const&) const;

private:
    bool _fixed;
    double _nSigmaForRadius;
    double _nRadiusForFlux;
    double _background;
    double _shiftmax;
};


#if HAVE_ellipticalFootprint
/************************************************************************************************************/
/**
 * Create an elliptical Footprint.
 */
PTR(afwDet::Footprint)
ellipticalFootprint(afwGeom::Point2I const& center, //!< The center of the circle
                    double a,                       //!< Major axis (pixels)
                    double b,                       //!< Minor axis (pixels)
                    double theta,                   //!< angle of major axis from x-axis; (radians)
                    afwGeom::Box2I const& region=afwGeom::Box2I() //!< Bounding box of MaskedImage footprint
                   )
{
    PTR(afwDet::Footprint) foot(new afwDet::Footprint);
    foot->setRegion(region);
    
    int const xc = center[0];           // x-centre
    int const yc = center[1];           // y-centre

    double const c = ::cos(theta);
    double const s = ::sin(theta);

    double const c0 = a*a*s*s + b*b*c*c;
    double const c1 = c*s*(a*a - b*b)/c0;
    double const c2 = a*b/c0;

    double const ymax = ::sqrt(c0) + 1; // max extent of ellipse in y-direction
    //
    // We go to quite a lot of annoying trouble to ensure that all pixels that are within or intercept
    // the ellipse are included in the Footprint
    //
    double x1, x2, y;
    for (int i = -ymax; i <= ymax; ++i) {
        double const dy = (i > 0) ? -0.5 : 0.5;
        y = i + dy;              // chord at top of pixel (above centre)
        if (c0 > y*y) {
            x1 = y*c1 - c2*std::sqrt(c0 - y*y);
            x2 = y*c1 + c2*std::sqrt(c0 - y*y);
        } else {
            x1 = x2 = y*c1;
        }

        y = i - dy;                     // chord at bottom of pixel (above centre)
        if (c0 > y*y) {
            double tmp = y*c1 - c2*std::sqrt(c0 - y*y);
            if (tmp < x1) {
                x1 = tmp;
            }

            tmp = y*c1 + c2*std::sqrt(c0 - y*y);
            if (tmp > x2) {
                x2 = tmp;
            }
        }

        foot->addSpan(yc + i, xc + x1 + 0.5, xc + x2 + 0.5);
    }

    return foot;
}
#endif


/************************************************************************************************************/

namespace {

template <typename MaskedImageT, typename WeightImageT>
class FootprintFindMoment : public afwDet::FootprintFunctor<MaskedImageT> {
public:
    FootprintFindMoment(MaskedImageT const& mimage, ///< The image the source lives in
                        double const xcen, double const ycen, // center of the object
                        double const ab,                      // axis ratio
                        double const theta                    // rotation of ellipse +ve from x axis
                       ) : afwDet::FootprintFunctor<MaskedImageT>(mimage),
                           _xcen(xcen), _ycen(ycen),
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
        afwGeom::Box2I const& bbox(foot.getBBox());
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
        double const dx = (x - _imageX0) - _xcen;
        double const dy = (y - _imageY0) - _ycen;
        double const du =  dx*_cosTheta + dy*_sinTheta;
        double const dv = -dx*_sinTheta + dy*_cosTheta;

        double r = ::hypot(du, dv*_ab); // ellipsoidal radius
#if 1
        if (dx*dx + dy*dy < 0.25) {     // within a pixel of the centre
            /*
             * We gain significant precision for flattened Gaussians by treating the central pixel specially
             *
             * If the object's centered in the pixel (and has constant surface brightness) we have <r> == eR;
             * if it's at the corner <r> = 2*eR; we interpolate between these exact results linearily in the
             * displacement.
             */
            
            double const eR = 0.38259771140356325; // <r> for a single square pixel, about the centre
            r = (eR/_ab)*(1 + afwGeom::ROOT2*::hypot(::fmod(du, 1), ::fmod(dv, 1)));
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
    KronAperture(double x, double y, afwEllipse::Axes const& ellipse) : _x(x), _y(y), _ellipse(ellipse) {}

    /// Accessors
    double getX() const { return _x; }
    double getY() const { return _y; }
    afwEllipse::Axes getEllipse() const { return _ellipse; }

    /// Determine the Kron Aperture from an image
    template<typename ImageT>
    static PTR(KronAperture) determine(ImageT const& image, // Image to measure
                                       afwDet::Source const& source, // Source with measurements
                                       double xcen, double ycen, // Centre of source
                                       double nSigmaForRadius,   // Multiplier for Kron radius
                                       double background, // Background to remove
                                       double shiftmax // Maximum shift permitted
        );

    /// Photometer within the Kron Aperture on an image
    template<typename ImageT>
    std::pair<double, double> measure(ImageT const& image, // Image to measure
                                      double nRadiusForFlux // Kron radius multiplier
        ) const;

    /// Transform a Kron Aperture to a different frame
    PTR(KronAperture) transform(afwGeom::AffineTransform const& trans) const {
        afwGeom::Point2D const center = trans(afwGeom::Point2D(_x, _y));
        afwEllipse::Axes const ellipse(_ellipse.transform(trans.getLinear()));
        return boost::make_shared<KronAperture>(center.getX(), center.getY(), ellipse);
    }

private:
    double _x, _y;                      // Centre
    afwEllipse::Axes _ellipse;          // Ellipse defining aperture shape
};

/*
 * Estimate the object's moments using the SDSS adaptive moments algorithm
 */
template<typename ImageT>
PTR(KronAperture) KronAperture::determine(ImageT const& image, // Image to measure
                                          afwDet::Source const& source, // Source with measurements
                                          double xcen, double ycen, // Centre of source
                                          double nSigmaForRadius,   // Multiplier for Kron radius
                                          double background, // Background to remove
                                          double shiftmax // Maximum shift permitted
    )
{
    /*
     * Estimate the object's moments using the SDSS adaptive moments algorithm
     */
    double const NaN = std::numeric_limits<double>::quiet_NaN();
    double Ixx = NaN; // <xx>
    double Ixy = NaN; // <xy>
    double Iyy = NaN; // <yy>
    short flags = 0;                    // Status flags
    try {
        CONST_PTR(afwDet::Measurement<afwDet::Shape>) shape = source.getShape();

        Ixx = shape->find("SDSS")->getIxx();
        Ixy = shape->find("SDSS")->getIxy();
        Iyy = shape->find("SDSS")->getIyy();
        flags = shape->find("SDSS")->getShapeStatus();
    } catch (pexExceptions::Exception& e) {
        detail::SdssShapeImpl shapeImpl;
        
        if (!detail::getAdaptiveMoments(image, background, xcen, ycen, shiftmax, &shapeImpl)) {
            std::string const& msg = "Failed to estimate adaptive moments while measuring KRON radius";
            LSST_EXCEPT_ADD(e, msg);
            throw e;
        }
        Ixx = shapeImpl.getIxx();
        Ixy = shapeImpl.getIxy();
        Iyy = shapeImpl.getIyy();
        flags = shapeImpl.getFlags();
    }
    if (flags & (Flags::SHAPE_MAXITER | Flags::SHAPE_UNWEIGHTED | Flags::SHAPE_UNWEIGHTED_BAD)) {
        // Don't trust the adaptive moment: they could make us take forever measuring a very large aperture
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Unable to measure adaptive moments");
    }

    /*
     * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
     * <x^2> + <y^2> = <u^2> + <v^2>
     * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
     * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
     */
    double const Iuu_p_Ivv = Ixx + Iyy;                             // <u^2> + <v^2>
    double const Iuu_m_Ivv = ::sqrt(::pow(Ixx - Iyy, 2) + 4*::pow(Ixy, 2)); // <u^2> - <v^2>
    double const Iuu = 0.5*(Iuu_p_Ivv + Iuu_m_Ivv);                         // (major axis)^2; a
    double const Ivv = 0.5*(Iuu_p_Ivv - Iuu_m_Ivv);                         // (minor axis)^2; b
    double const theta = 0.5*::atan2(2*Ixy, Ixx - Iyy);                     // angle of a +ve from x axis

    double const a = nSigmaForRadius*::sqrt(Iuu);
    double const b = nSigmaForRadius*::sqrt(Ivv);

    FootprintFindMoment<ImageT, afwDet::Psf::Image> iRFunctor(image, xcen, ycen, a/b, theta);

    // Build an elliptical Footprint of the proper size
    afwGeom::Point2I center(xcen + 0.5, ycen + 0.5);
#if HAVE_ellipticalFootprint
    PTR(afwDet::Footprint) foot(ellipticalFootprint(center, a, b, theta));
    iRFunctor.apply(*foot);
#else
    afwDet::Footprint foot(center, a, b, theta);
    iRFunctor.apply(foot);
#endif

    if (!iRFunctor.getGood()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Bad footprint when determining Kron aperture");
    }

    double const radius = iRFunctor.getIr();

    return boost::make_shared<KronAperture>(xcen, ycen, afwEllipse::Axes(radius, radius * b / a, theta));
}

template<typename ImageT>
std::pair<double, double> KronAperture::measure(ImageT const& image, // Image of interest
                                                double nRadiusForFlux // Kron radius multiplier
    ) const
{
    try {
        double const r2 = nRadiusForFlux * _ellipse.getA() * _ellipse.getA();
        double const ellip = 1.0 - _ellipse.getB()/_ellipse.getA();
        return photometry::calculateSincApertureFlux(image, _x, _y, 0.0, r2, _ellipse.getTheta(), ellip);
    } catch(pexExceptions::LengthErrorException &e) {
        LSST_EXCEPT_ADD(e, (boost::format("Measuring Kron flux for object at (%.3f, %.3f);"
                                          " aperture radius %g,%g theta %g")
                            % _x % _y % _ellipse.getA() % _ellipse.getB() %
                            afwGeom::radToDeg(_ellipse.getTheta())).str());
        throw e;
    }
}


} // anonymous namespace

/************************************************************************************************************/

/**
 * Calculate the desired Kron radius and flux
 */
template<typename ExposureT>
PTR(afwDet::Photometry) KronPhotometer<ExposureT>::measureSingle(
    afwDet::Source const& target,
    afwDet::Source const& source,
    ExposurePatch<ExposureT> const& patch
    ) const
{
    typedef typename ExposureT::MaskedImageT MaskedImageT;

    CONST_PTR(ExposureT) exposure = patch.getExposure();
    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = patch.getCenter().getX() - mimage.getX0(); // column position in image pixel coords
    double const ycen = patch.getCenter().getY() - mimage.getY0();  // row position

    CONST_PTR(KronAperture) aperture = KronAperture::determine(mimage, source, xcen, ycen, _nSigmaForRadius, 
                                                               _background, _shiftmax);
    std::pair<double, double> const& result = aperture->measure(mimage, _nRadiusForFlux);
    double const flux = result.first;
    double const fluxErr = result.second;
    double const radius = aperture->getEllipse().getA();
    return boost::make_shared<afwDet::AperturePhotometry>(flux, fluxErr, radius);
}


/// Declare the existence of a "KRON" algorithm to MeasurePhotometry
LSST_DECLARE_ALGORITHM(KronPhotometer, afwDet::Photometry);

}}}
