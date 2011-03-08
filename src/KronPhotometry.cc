// -*- LSST-C++ -*-
#include <numeric>
#include <cmath>
#include <functional>
#include "lsst/pex/exceptions.h"
#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math/Integrate.h"
#include "lsst/meas/algorithms/Measure.h"

#include "lsst/afw/detection/Psf.h"
#include "lsst/afw/detection/Photometry.h"
#include "lsst/meas/algorithms/detail/SdssShape.h"
#include "lsst/meas/algorithms/Photometry.h"

namespace pexExceptions = lsst::pex::exceptions;
namespace pexLogging = lsst::pex::logging;
namespace afwDetection = lsst::afw::detection;
namespace afwGeom = lsst::afw::geom;
namespace afwImage = lsst::afw::image;
namespace afwMath = lsst::afw::math;

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief A class that knows how to calculate fluxes using the KRON photometry algorithm
 * @ingroup meas/algorithms
 */
class KronPhotometry : public afwDetection::Photometry
{
    enum { RADIUS = Photometry::NVALUE,
           NVALUE = RADIUS + 1 };

public:
    typedef boost::shared_ptr<KronPhotometry> Ptr;
    typedef boost::shared_ptr<KronPhotometry const> ConstPtr;

    /// Ctor
    KronPhotometry(double radius, double flux, double fluxErr=std::numeric_limits<double>::quiet_NaN()) :
        afwDetection::Photometry() {
        init();                         // This allocates space for everything in the schema

        set<FLUX>(flux);
        set<FLUX_ERR>(fluxErr);
        set<RADIUS>(radius);
    }

    /// Add desired fields to the schema
    virtual void defineSchema(afwDetection::Schema::Ptr schema ///< our schema; == _mySchema
                     ) {
        Photometry::defineSchema(schema);
        schema->add(afwDetection::SchemaEntry("radius", RADIUS, afwDetection::Schema::DOUBLE, 1, "pixels"));
    }

    double getRadius(int) const {
        return get<RADIUS, double>();
    }

    static bool doConfigure(lsst::pex::policy::Policy const& policy)
    {
        if (policy.isDouble("nSigmaForRad")) {
            _nSigmaForRad = policy.getDouble("nSigmaForRad");
        }
        if (policy.isDouble("background")) {
            _background = policy.getDouble("background");
        } 
        if (policy.isDouble("shiftmax")) {
            _shiftmax = policy.getDouble("shiftmax");
        } 
        
        return true;
    }
    template<typename ImageT>
    static Photometry::Ptr doMeasure(CONST_PTR(ImageT),
                                     CONST_PTR(afwDetection::Peak),
                                     CONST_PTR(afwDetection::Source)
                                    );

private:
    static double _nSigmaForRad;
    static double _background;
    static double _shiftmax;

    KronPhotometry(void) : afwDetection::Photometry() { }
    LSST_SERIALIZE_PARENT(afwDetection::Photometry)
};

LSST_REGISTER_SERIALIZER(KronPhotometry)

double KronPhotometry::_nSigmaForRad = 6.0;   // Size of aperture (in sigma) to estimate Kron radius
double KronPhotometry::_background = 0.0;     // the frame's background level
double KronPhotometry::_shiftmax = 10;        // Max allowed centroid shift

/************************************************************************************************************/

namespace {
    
template<typename ImageT>
std::pair<double, double>
getKronFlux(
        ImageT const& mimage,           // the data to process
        double background,               // background level
        double xcen, double ycen,         // centre of object
        double shiftmax                  // max allowed centroid shift
               )
{
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    detail::SdssShapeImpl shapeImpl;

    if (!detail::getAdaptiveMoments(mimage, background, xcen, ycen, shiftmax, &shapeImpl)) {
        ;                               // Should set a flag here
    } else {
        /*
         * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
         * <x^2> + <y^2> = <u^2> + <v^2>
         * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
         * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
         */
        double const Mxx = shapeImpl.getIxx(); // <x^2>
        double const Mxy = shapeImpl.getIxy(); // <xy>
        double const Myy = shapeImpl.getIyy(); // <y^2>
        
        double const Muu_p_Mvv = Mxx + Myy;                             // <u^2> + <v^2>
        double const Muu_m_Mvv = ::sqrt(::pow(Mxx - Myy, 2) + 4*::pow(Mxy, 2)); // <u^2> - <v^2>
        double const Muu = 0.5*(Muu_p_Mvv + Muu_m_Mvv);
        double const Mvv = 0.5*(Muu_p_Mvv - Muu_m_Mvv);
        
        double const scale = 2*M_PI*::sqrt(Muu*Mvv);
        flux = scale*shapeImpl.getI0();
        fluxErr = scale*shapeImpl.getI0Err();
    }

    return std::make_pair(flux, fluxErr);
}
}

/************************************************************************************************************/
/**
 * Calculate the desired kron flux
 */
template<typename ExposureT>
afwDetection::Photometry::Ptr KronPhotometry::doMeasure(CONST_PTR(ExposureT) exposure,
                                                            CONST_PTR(afwDetection::Peak) peak,
                                                            CONST_PTR(afwDetection::Source) source
                                                           )
{
    double radius = std::numeric_limits<double>::quiet_NaN();
    double flux = std::numeric_limits<double>::quiet_NaN();
    double fluxErr = std::numeric_limits<double>::quiet_NaN();

    if (!peak) {
        return boost::make_shared<KronPhotometry>(radius, flux, fluxErr);
    }

    typedef typename ExposureT::MaskedImageT MaskedImageT;
    typedef typename MaskedImageT::Image Image;
    typedef typename Image::Pixel Pixel;
    typedef typename Image::Ptr ImagePtr;

    MaskedImageT const& mimage = exposure->getMaskedImage();
    
    double const xcen = peak->getFx() - mimage.getX0(); ///< object's column position in image pixel coords
    double const ycen = peak->getFy() - mimage.getY0();  ///< object's row position
    /*
     * Estimate the object's moments using the SDSS adaptive moments algorithm
     */
    double Ixx = std::numeric_limits<double>::quiet_NaN(); // <xx>
    double Ixy = std::numeric_limits<double>::quiet_NaN(); // <xy>
    double Iyy = std::numeric_limits<double>::quiet_NaN(); // <yy>
    try {
        if (source) {
            CONST_PTR(lsst::afw::detection::Measurement<lsst::afw::detection::Shape>)
                shape = source->getShape();
            
            Ixx = shape->find("SDSS")->getIxx();
            Ixy = shape->find("SDSS")->getIxy();
            Iyy = shape->find("SDSS")->getIyy();
        } else {
            throw LSST_EXCEPT(lsst::pex::exceptions::NotFoundException, "No Source");
        }
    } catch (lsst::pex::exceptions::Exception& e) {
        detail::SdssShapeImpl shapeImpl;
        
        if (!detail::getAdaptiveMoments(mimage, _background, xcen, ycen, _shiftmax, &shapeImpl)) {
            LSST_EXCEPT_ADD(e, "Measuring KRON flux");
            throw e;
        }
        Ixx = shapeImpl.getIxx();
        Ixy = shapeImpl.getIxy();
        Iyy = shapeImpl.getIyy();
    }
#if 0
    std::cout << "SDSS moments " << " " << Ixx << " " << Ixy << " " << Iyy << std::endl;
#endif
    /*
     * The shape is an ellipse that's axis-aligned in (u, v) [<uv> = 0] after rotation by theta:
     * <x^2> + <y^2> = <u^2> + <v^2>
     * <x^2> - <y^2> = cos(2 theta)*(<u^2> - <v^2>)
     * 2*<xy>        = sin(2 theta)*(<u^2> - <v^2>)
     */
    double const Iuu_p_Ivv = Ixx + Iyy;                             // <u^2> + <v^2>
    double const Iuu_m_Ivv = ::sqrt(::pow(Ixx - Iyy, 2) + 4*::pow(Ixy, 2)); // <u^2> - <v^2>
    double const Iuu = 0.5*(Iuu_p_Ivv + Iuu_m_Ivv);
    double const Ivv = 0.5*(Iuu_p_Ivv - Iuu_m_Ivv);
    double const theta = 0.5*::atan2(2*Ixy, Ixx - Iyy);
    
    radius = ::sqrt(Iuu);               // major axis
    flux = ::sqrt(Ivv);                 // minor axis
    fluxErr = theta*180/M_PI;

#if 0
    {
        std::pair<double, double> flux_fluxErr = getKronFlux(mimage, _background, xcen, ycen, _shiftmax);
        flux = flux_fluxErr.first;
        fluxErr = flux_fluxErr.second;
    }
#endif

    return boost::make_shared<KronPhotometry>(radius, flux, fluxErr);
}

/*
 * Declare the existence of a "KRON" algorithm to MeasurePhotometry
 *
 * \cond
 */
#define MAKE_PHOTOMETRYS(TYPE)                                          \
    MeasurePhotometry<afwImage::Exposure<TYPE> >::declare("KRON", \
        &KronPhotometry::doMeasure<afwImage::Exposure<TYPE> >, \
        &KronPhotometry::doConfigure \
    )

namespace {
    volatile bool isInstance[] = {
        MAKE_PHOTOMETRYS(float)
#if 0
        ,MAKE_PHOTOMETRYS(double)
#endif
    };
}
    
// \endcond

}}}
