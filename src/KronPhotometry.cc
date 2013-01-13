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
#include "lsst/afw/table/Source.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/coord/Coord.h"
#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/ellipses.h"

#include "lsst/meas/algorithms/Photometry.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"

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
        _radiusKey(schema.addField<double>(ctrl.name + ".radius", "Kron radius (sqrt(a*b))")),
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

    afw::table::Key<double> _radiusKey;
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
            r = (eR/_ab)*(1 + afw::geom::ROOT2*::hypot(::fmod(du, 1), ::fmod(dv, 1)));
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
    static PTR(KronAperture) determine(ImageT const& image, // Image to measure
                                       afw::table::SourceRecord const& source, // Source with measurements
                                       afw::geom::Point2D const& center, // Center of source
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
 * Estimate the object's moments using the SDSS adaptive moments algorithm
 */
template<typename ImageT>
PTR(KronAperture) KronAperture::determine(ImageT const& image, // Image to measure
                                          afw::table::SourceRecord const& source, // Source with measurements
                                          afw::geom::Point2D const& center, // Centre of source
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

    typedef algorithms::detail::SdssShapeImpl ShapeImpl;

    ShapeImpl shapeImpl;
    
    if (!algorithms::detail::getAdaptiveMoments(image, background, center.getX(), center.getY(),
                                                shiftmax, &shapeImpl)) {
        std::string msg = "Failed to estimate adaptive moments while measuring KRON radius";
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, msg);        
    }
    Ixx = shapeImpl.getIxx();
    Ixy = shapeImpl.getIxy();
    Iyy = shapeImpl.getIyy();

    if (shapeImpl.getFlag(ShapeImpl::MAXITER)
        || shapeImpl.getFlag(ShapeImpl::UNWEIGHTED)
        || shapeImpl.getFlag(ShapeImpl::UNWEIGHTED_BAD)
    ) {
        // Don't trust the adaptive moment: they could make us take forever measuring a very large aperture
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Unable to measure adaptive moments");
    }

    afw::geom::ellipses::Axes axes(afw::geom::ellipses::Quadrupole(Ixx, Iyy, Ixy));
    axes.scale(nSigmaForRadius);

    FootprintFindMoment<ImageT, afwDet::Psf::Image> iRFunctor(
        image, center, axes.getA() / axes.getB(), axes.getTheta()
    );

    // Build an elliptical Footprint of the proper size
    afwDet::Footprint foot(afw::geom::ellipses::Ellipse(axes, center));

    iRFunctor.apply(foot);

    if (!iRFunctor.getGood()) {
        throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException,
                          "Bad footprint when determining Kron aperture");
    }

    double const radius = iRFunctor.getIr();

    return boost::make_shared<KronAperture>(
        center, afw::geom::ellipses::Axes(radius, radius * axes.getB() / axes.getA(), axes.getTheta())
    );
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

#if 0
    // XXX Do we have to worry about this?
    double const xcen = patch.getCenter().getX() - mimage.getX0(); // column position in image pixel coords
    double const ycen = patch.getCenter().getY() - mimage.getY0();  // row position
#endif

    KronFluxControl const & ctrl = static_cast<KronFluxControl const &>(this->getControl());

    if (source.getShapeFlag()) {        // the shape's bad; give up now
        return;
    }
    
    CONST_PTR(KronAperture) aperture;
    if (ctrl.fixed) {
        aperture.reset(new KronAperture(source));
    } else {
        try {
            aperture = KronAperture::determine(mimage, source, center, ctrl.nSigmaForRadius,
                                               ctrl.background, ctrl.shiftmax);
        } catch(pex::exceptions::Exception& e) {
            return;
        }
    }
    source.set(_badRadiusKey, false);

    std::pair<double, double> result = aperture->measure(mimage, ctrl.nRadiusForFlux);
    source.set(getKeys().meas, result.first);
    source.set(getKeys().err, result.second);
    source.set(_radiusKey, aperture->getEllipse().getDeterminantRadius());
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
