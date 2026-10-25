#ifndef GOLLUMFIT_FASTWEIGHTER_H_
#define GOLLUMFIT_FASTWEIGHTER_H_

#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include <PhysTools/autodiff.h>

#include "Event.h"
#include "GollumParameters.h"
#include "ReducedSplines.h"
#include "SystematicsCache.h"

/**
* @file FastWeighter.h
* @brief Closed-form event weighter (value and analytic gradient).
*
* Computes the same weight as the composed weighter chain built by
* sterile::WeighterMaker (see analysisWeighting.h):
*
*   w = convNorm * (conv + prompt + astro) * icegrad            [enableTotalNorm]
*   w = (convNorm * conv + prompt + astro) * icegrad             [otherwise]
*
*   conv   = H_c * D_c * T_c * (1+adu*ad)(1+klu*kl) * F,  F = convFlux + sum_i h_i*had_i + sum_j cr_j*cr_j (FLT_MAX if negative)
*   prompt = promptNorm * H_p * D_p * T_p * promptFlux
*   astro  = astroNorm  * H_a * D_a * T_a * astroFlux * nn * tilt
*   icegrad = prod_k (1 + g_k * icegrad_k)
*
* with H (hole ice) and D (DOM efficiency) = 10^(rate(p) - rate(nominal)),
* T (attenuation) a spline in nuxs/nubarxs, nn the nu/nubar balance and tilt the
* broken power law in the primary energy. The spline factors come from the
* per-event reduced splines in SystematicsCache, so a weight costs a few dozen
* flops instead of 18 three-dimensional spline evaluations, and the gradient
* with respect to the 38 parameters is written in closed form instead of being
* propagated through 38-wide dual numbers.
*/

namespace gollumfit {

struct FastWeighterParams {
    double convNorm, promptNorm, adu, klu;
    double had[10];        ///< hadronic yield parameters (indices 4..13)
    double cr[6];          ///< cosmic-ray spectrum parameters (14..19)
    double icegrad[9];     ///< ice gradient parameters (20..28)
    double domEff;         ///< 29
    double holeIce;        ///< 30
    double astroNorm;      ///< 31
    double astroDeltaGamma, astroDeltaGammaSec; ///< 32, 33
    double astroPivot;     ///< 34 (log10 of the break energy)
    double nuNubarRatio;   ///< 35
    double nuxs, nubarxs;  ///< 36, 37

    template <typename T>
    static double val(const T& x) { return static_cast<double>(x); }
    template <unsigned int N>
    static double val(const phys_tools::autodiff::FD<N>& x) { return x.value(); }

    template <typename DataType>
    explicit FastWeighterParams(const std::vector<DataType>& p) {
        if (p.size() < 38)
            throw std::runtime_error("FastWeighter: expected 38 fit parameters");
        convNorm = val(p[0]); promptNorm = val(p[1]); adu = val(p[2]); klu = val(p[3]);
        for (int i = 0; i < 10; ++i) had[i] = val(p[4 + i]);
        for (int i = 0; i < 6; ++i) cr[i] = val(p[14 + i]);
        for (int i = 0; i < 9; ++i) icegrad[i] = val(p[20 + i]);
        domEff = val(p[29]); holeIce = val(p[30]); astroNorm = val(p[31]);
        astroDeltaGamma = val(p[32]); astroDeltaGammaSec = val(p[33]); astroPivot = val(p[34]);
        nuNubarRatio = val(p[35]); nuxs = val(p[36]); nubarxs = val(p[37]);
    }
};

/// Shared implementation: evaluates the weight and, if grad != nullptr, its 38 partial derivatives.
class FastWeighterCore {
public:
    static constexpr unsigned int kNParams = 38;
    static constexpr double kLn10 = 2.302585092994046;

    FastWeighterCore(const SystematicsCache& cache, const FastWeighterParams& p, bool enableTotalNorm)
        : cache_(&cache), p_(p), enableTotalNorm_(enableTotalNorm) {
        // basis functions along the parameter axis, once per evaluation
        for (const ReducedAxisSpline& t : cache.domTables()) domBasis_.push_back(t.basis(p.domEff));
        for (const ReducedAxisSpline& t : cache.holeTables()) holeBasis_.push_back(t.basis(p.holeIce));
        for (const ReducedAxisSpline& t : cache.attTables()) {
            attBasisNu_.push_back(t.basis(p.nuxs));
            attBasisNubar_.push_back(t.basis(p.nubarxs));
        }
    }

    template <bool WithGrad>
    inline double evaluate(const Event& e, double* grad) const {
        const SystematicsCache& cache = *cache_;
        const size_t idx = cache.index(e);
        const SystematicsCache::EventInfo& info = cache.info(idx);
        const SplineReal* row = cache.row(idx);
        const SystematicsCache::Layout& L = cache.layout();
        const FastWeighterParams& p = p_;

        // --- ice gradients: A = prod_k (1 + g_k c_k)
        const double ig[9] = {e.cachedIceGrad0, e.cachedIceGrad1, e.cachedIceGrad2, e.cachedIceGrad3, e.cachedIceGrad4,
                              e.cachedIceGrad5, e.cachedIceGrad6, e.cachedIceGrad7, e.cachedIceGrad8};
        double igf[9], A = 1.0;
        for (int k = 0; k < 9; ++k) { igf[k] = 1.0 + p.icegrad[k] * ig[k]; A *= igf[k]; }

        // --- atmospheric flux shape factors
        const double adf = 1.0 + p.adu * e.cachedAtmDensity;
        const double klf = 1.0 + p.klu * e.cachedKaonLosses;
        const double B = adf * klf;

        // --- conventional flux with hadronic / cosmic-ray gradients
        const double had[10] = {e.cachedHadronicHEkp, e.cachedHadronicHEkm, e.cachedHadronicVHE1pip, e.cachedHadronicVHE1pim,
                                e.cachedHadronicVHE3kp, e.cachedHadronicVHE3km, e.cachedHadronicVHE3pip, e.cachedHadronicVHE3pim,
                                e.cachedHadronicVHE3p, e.cachedHadronicVHE3n};
        const double cr[6] = {e.cachedCosmicRay1, e.cachedCosmicRay2, e.cachedCosmicRay3,
                              e.cachedCosmicRay4, e.cachedCosmicRay5, e.cachedCosmicRay6};
        double hadronic = 0.0, cosmic = 0.0;
        for (int i = 0; i < 10; ++i) hadronic += p.had[i] * had[i];
        for (int j = 0; j < 6; ++j) cosmic += p.cr[j] * cr[j];
        double F = e.cachedConvWeight + cosmic + hadronic;
        const bool Fneg = (F < 0.0);
        if (Fneg)
            F = static_cast<double>(std::numeric_limits<float>::max()); // unphysical region sentinel, no gradient

        // --- detector systematics and attenuation per flux component
        double D[3], dD[3], H[3], dH[3], T[3], dT[3];
        for (unsigned int c = 0; c < 3; ++c) {
            detectorFactor(info.domTable[c], info.domNominal[c], row + L.domOffset[c], domBasis_, D[c], dD[c]);
            detectorFactor(info.holeTable[c], info.holeNominal[c], row + L.holeOffset[c], holeBasis_, H[c], dH[c]);
            const int8_t at = info.attTable[c];
            if (at == SystematicsCache::kTableNone) { T[c] = 1.0; dT[c] = 0.0; }
            else if (at == SystematicsCache::kTableZero) { T[c] = 0.0; dT[c] = 0.0; }
            else {
                const ReducedAxisSpline::Basis& bs = info.isNeutrino ? attBasisNu_[at] : attBasisNubar_[at];
                if (!bs.valid) { T[c] = 0.0; dT[c] = 0.0; }
                else {
                    T[c] = ReducedAxisSpline::eval(row + L.attOffset[c], bs);
                    if (T[c] < 0.0)
                        throw std::runtime_error("GollumFit::attenuation weighter making a neg weight.");
                    dT[c] = WithGrad ? ReducedAxisSpline::evalDeriv(row + L.attOffset[c], bs) : 0.0;
                }
            }
        }

        // --- astrophysical shape: nu/nubar balance and broken power-law tilt
        const double nn = info.isAntiparticle ? p.nuNubarRatio : 2.0 - p.nuNubarRatio;
        const double lnE = info.logPrimaryEnergy;
        const double lnEp = p.astroPivot * kLn10;
        const bool highBranch = (lnE > lnEp);
        const double dg = highBranch ? p.astroDeltaGammaSec : p.astroDeltaGamma;
        const double tilt = std::exp(-dg * (lnE - lnEp));

        // --- components (without their normalisations)
        const double convNoB = H[0] * D[0] * T[0] * F;
        const double conv = convNoB * B;
        const double promptBase = H[1] * D[1] * T[1] * e.cachedPromptWeight;
        const double astroNoTilt = H[2] * D[2] * T[2] * e.cachedAstroWeight * nn;
        const double astroBase = astroNoTilt * tilt;
        const double prompt = p.promptNorm * promptBase;
        const double astro = p.astroNorm * astroBase;

        const double nC = p.convNorm;
        const double nP = enableTotalNorm_ ? p.convNorm : 1.0;
        const double nA = nP;
        const double S = nC * conv + nP * prompt + nA * astro;
        const double w = S * A;

        if (WithGrad) {
            double* g = grad;
            g[0] = (enableTotalNorm_ ? (conv + prompt + astro) : conv) * A;
            g[1] = nP * promptBase * A;
            const double convA = nC * A;
            const double convNoBF = H[0] * D[0] * T[0];
            g[2] = convA * convNoBF * F * e.cachedAtmDensity * klf;
            g[3] = convA * convNoBF * F * e.cachedKaonLosses * adf;
            const double convNoF = Fneg ? 0.0 : convA * convNoBF * B;
            for (int i = 0; i < 10; ++i) g[4 + i] = convNoF * had[i];
            for (int j = 0; j < 6; ++j) g[14 + j] = convNoF * cr[j];
            // d A / d g_k = c_k * prod_{j != k} (1 + g_j c_j), via prefix/suffix products
            double prefix[10], suffix[10];
            prefix[0] = 1.0;
            for (int k = 0; k < 9; ++k) prefix[k + 1] = prefix[k] * igf[k];
            suffix[9] = 1.0;
            for (int k = 8; k >= 0; --k) suffix[k] = suffix[k + 1] * igf[k];
            for (int k = 0; k < 9; ++k) g[20 + k] = S * ig[k] * prefix[k] * suffix[k + 1];
            // DOM efficiency and hole ice: each factor is 10^(r(p)-r0), d/dp = factor * ln10 * r'(p)
            g[29] = A * (nC * H[0] * T[0] * B * F * dD[0]
                         + nP * p.promptNorm * H[1] * T[1] * e.cachedPromptWeight * dD[1]
                         + nA * p.astroNorm * H[2] * T[2] * e.cachedAstroWeight * nn * tilt * dD[2]);
            g[30] = A * (nC * D[0] * T[0] * B * F * dH[0]
                         + nP * p.promptNorm * D[1] * T[1] * e.cachedPromptWeight * dH[1]
                         + nA * p.astroNorm * D[2] * T[2] * e.cachedAstroWeight * nn * tilt * dH[2]);
            g[31] = nA * astroBase * A;
            const double dTiltDg = -(lnE - lnEp) * tilt;
            const double astroTiltCoef = nA * p.astroNorm * astroNoTilt * A;
            g[32] = highBranch ? 0.0 : astroTiltCoef * dTiltDg;
            g[33] = highBranch ? astroTiltCoef * dTiltDg : 0.0;
            g[34] = astroTiltCoef * (dg * kLn10 * tilt);
            g[35] = nA * p.astroNorm * H[2] * D[2] * T[2] * e.cachedAstroWeight * tilt * A * (info.isAntiparticle ? 1.0 : -1.0);
            const double dTsum = A * (nC * H[0] * D[0] * B * F * dT[0]
                                      + nP * p.promptNorm * H[1] * D[1] * e.cachedPromptWeight * dT[1]
                                      + nA * p.astroNorm * H[2] * D[2] * e.cachedAstroWeight * nn * tilt * dT[2]);
            g[36] = info.isNeutrino ? dTsum : 0.0;
            g[37] = info.isNeutrino ? 0.0 : dTsum;
        }
        return w;
    }

private:
    const SystematicsCache* cache_;
    FastWeighterParams p_;
    bool enableTotalNorm_;
    std::vector<ReducedAxisSpline::Basis> domBasis_, holeBasis_, attBasisNu_, attBasisNubar_;

    /// DOM-efficiency / hole-ice factor 10^(rate(p) - rate(nominal)) and its derivative wrt p.
    static inline void detectorFactor(int8_t table, SplineReal nominal, const SplineReal* coef,
                                      const std::vector<ReducedAxisSpline::Basis>& bases, double& f, double& df) {
        if (table == SystematicsCache::kTableNone) { f = 1.0; df = 0.0; return; }
        if (table == SystematicsCache::kTableZero || nominal == SystematicsCache::kNominalSentinel) { f = 0.0; df = 0.0; return; }
        const ReducedAxisSpline::Basis& bs = bases[table];
        if (!bs.valid) { f = 0.0; df = 0.0; return; } // parameter outside the spline range: event ignored
        const SplineReal r = ReducedAxisSpline::eval(coef, bs);
        if (r == 0.f) { f = 0.0; df = 0.0; return; }
        f = std::exp(kLn10 * (static_cast<double>(r) - static_cast<double>(nominal)));
        df = f * kLn10 * static_cast<double>(ReducedAxisSpline::evalDeriv(coef, bs));
    }
};

/// Weight only.
class FastWeighter {
public:
    using result_type = double;
    FastWeighter(std::shared_ptr<const SystematicsCache> cache, const FastWeighterParams& p, bool enableTotalNorm)
        : cache_(std::move(cache)), core_(*cache_, p, enableTotalNorm) {}
    inline double operator()(const Event& e) const { return core_.evaluate<false>(e, nullptr); }
private:
    std::shared_ptr<const SystematicsCache> cache_; ///< keeps the table alive
    FastWeighterCore core_;
};

/// Weight with its gradient as a 38-dimensional dual number.
class FastGradWeighter {
public:
    using result_type = phys_tools::autodiff::FD<38>;
    FastGradWeighter(std::shared_ptr<const SystematicsCache> cache, const FastWeighterParams& p, bool enableTotalNorm)
        : cache_(std::move(cache)), core_(*cache_, p, enableTotalNorm) {}
    inline result_type operator()(const Event& e) const {
        double g[38];
        const double w = core_.evaluate<true>(e, g);
        result_type r(w);
        for (unsigned int i = 0; i < 38; ++i)
            r.setDerivative(i, g[i]);
        return r;
    }
private:
    std::shared_ptr<const SystematicsCache> cache_;
    FastWeighterCore core_;
};

} // namespace gollumfit

#endif /* GOLLUMFIT_FASTWEIGHTER_H_ */
