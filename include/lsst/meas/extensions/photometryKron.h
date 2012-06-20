// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
#define LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H

#include "lsst/meas/algorithms/FluxControl.h"

namespace lsst { namespace meas { namespace extensions { namespace photometryKron {

/**
 *  @brief C++ control object for Gaussian flux.
 *
 *  @sa GaussianFluxConfig.
 */
class KronFluxControl : public algorithms::FluxControl {
public:

    LSST_CONTROL_FIELD(fixed, bool,
                       "if true, use existing shape and centroid measurements instead of fitting");
    LSST_CONTROL_FIELD(background, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(shiftmax, double, "FIXME! NEVER DOCUMENTED!");
    LSST_CONTROL_FIELD(nSigmaForRadius, double, "Number of sigma to set Kron radius");
    LSST_CONTROL_FIELD(nRadiusForFlux, double, "Number of Kron radii for Kron flux");

    KronFluxControl() : 
        algorithms::FluxControl("flux.kron"), fixed(false), background(0.0), shiftmax(10.0),
        nSigmaForRadius(6.0), nRadiusForFlux(2.0)
    {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::photometryKron

#endif // !LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
