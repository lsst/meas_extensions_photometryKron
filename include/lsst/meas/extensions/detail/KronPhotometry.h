#if !defined(LSST_MEAS_EXTENSIONS_DETAIL_KRONPHOTOMETRY_H)
#define LSST_MEAS_EXTENSIONS_DETAIL_KRONPHOTOMETRY_H 1

#include "lsst/afw/detection/Photometry.h"

namespace lsst {
    namespace afw {
        namespace detection {
            class Peak;
            class Source;
        }
    }
    namespace meas {
        namespace algorithms {
            namespace detail {
/**
 * @brief A class that knows how to calculate fluxes using the KRON photometry algorithm
 * @ingroup meas/algorithms
 */
class KronPhotometry : public lsst::afw::detection::Photometry
{
    enum { RADIUS=Photometry::NVALUE,
           NVALUE   = RADIUS   + 1 };

public:
    typedef boost::shared_ptr<KronPhotometry> Ptr;
    typedef boost::shared_ptr<KronPhotometry const> ConstPtr;

    /// Ctor
    KronPhotometry(double radius, double flux, double fluxErr=std::numeric_limits<double>::quiet_NaN());

    /// Add desired fields to the schema
    virtual void defineSchema(lsst::afw::detection::Schema::Ptr schema);

    virtual double getParameter(int=0) const;

    static bool doConfigure(lsst::pex::policy::Policy const& policy);

    template<typename ImageT>
    static Photometry::Ptr doMeasure(CONST_PTR(ImageT),
                                     CONST_PTR(lsst::afw::detection::Peak),
                                     CONST_PTR(lsst::afw::detection::Source)
                                    );

private:
    static double _nSigmaForRad;
    static double _background;
    static double _shiftmax;

    KronPhotometry(void) : lsst::afw::detection::Photometry() { }
    LSST_SERIALIZE_PARENT(lsst::afw::detection::Photometry)
};
            }
        }
    }
}       
#endif
