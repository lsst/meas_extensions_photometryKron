// -*- lsst-c++ -*-
#ifndef LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
#define LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H

#include "lsst/meas/algorithms/FluxControl.h"

namespace lsst { namespace meas { namespace extensions { namespace photometryKron {

/**
 *  @brief C++ control object for Kron flux.
 *
 *  @sa KronFluxConfig.
 */
class KronFluxControl : public algorithms::FluxControl {
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

    KronFluxControl() : 
        algorithms::FluxControl("flux.kron"), fixed(false),
        nSigmaForRadius(6.0),
        nIterForRadius(1),
        nRadiusForFlux(2.5),
        maxSincRadius(10.0),
        minimumRadius(0.0),
        enforceMinimumRadius(true),
        useFootprintRadius(false),
        smoothingSigma(-1.0)
    {}

private:
    virtual PTR(algorithms::AlgorithmControl) _clone() const;
    virtual PTR(algorithms::Algorithm) _makeAlgorithm(
        afw::table::Schema & schema, PTR(daf::base::PropertyList) const & metadata
    ) const;
};

}}}} // namespace lsst::meas::extensions::photometryKron

#endif // !LSST_MEAS_EXTENSIONS_PHOTOMETRY_KRON_H
