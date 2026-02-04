/**
 * @file GollumFitDiverWrapper.h
 * @brief Wrapper for the Diver differential evolution optimizer.
 *
 * Provides integration of the Diver global optimizer with GollumFit
 * for finding the global minimum of the likelihood surface before
 * local refinement with MinLLH.
 */

#ifndef GOLLUMFITDIVERWRAPPER_H_
#define GOLLUMFITDIVERWRAPPER_H_

#include "GollumFit.h"
#include "diver.hpp"

#include <vector>
#include <memory>
#include <random>
#include <limits>
#include <cmath>
#include <string>
#include <atomic>
#include <iostream>
#include <unordered_map>
#include <algorithm>

namespace gollumfit {
namespace diver {

/// Penalty value returned for invalid parameter combinations
static constexpr double BAD_OBJ = 1e12;

/// @name Diagnostic Counters
/// @brief Thread-safe counters for tracking optimizer evaluations
/// @{
static std::atomic<long long> n_eval{0};  ///< Total objective function evaluations
static std::atomic<long long> n_good{0};  ///< Successful (finite) evaluations
static std::atomic<long long> n_bad{0};   ///< Failed/invalid evaluations
/// @}

/// Current free indices for objective function (set by Dive before calling cdiver)
static std::vector<int> current_free_indices;

/**
 * @brief Returns default indices of free parameters for Diver optimization.
 *
 * These 9 parameters are the most important for global basin-finding.
 * Other nuisance parameters are held fixed at their seeded values during
 * Diver exploration and then freed for local refinement with MinLLH.
 *
 * @return Reference to static vector of parameter indices
 */
inline const std::vector<int>& default_free_indices() {
    static const std::vector<int> idx = {
        0,   // convNorm
        8,   // hadronicVHE3kp
        9,   // hadronicVHE3km
        29,  // domEfficiency
        30,  // holeiceForward
        31,  // astroNorm
        32,  // astroDeltaGamma
        33,  // astroDeltaGammaSec
        34   // astroPivot
    };
    return idx;
}

/**
 * @brief Get free indices to use - custom if provided, otherwise default.
 *
 * @param custom_indices Custom indices from SteeringParams (empty = use default)
 * @return Vector of parameter indices to optimize
 */
inline std::vector<int> get_free_indices(const std::vector<int>& custom_indices) {
    if (custom_indices.empty()) {
        return default_free_indices();
    }
    return custom_indices;
}

/**
 * @struct PriorPair
 * @brief Holds Gaussian prior parameters for a single nuisance parameter.
 */
struct PriorPair {
    double mu;   ///< Prior mean
    double sig;  ///< Prior standard deviation
};

/**
 * @brief Builds a lookup map from parameter index to prior (mu, sigma).
 *
 * @param p Priors object containing all prior center/width values
 * @return Map from parameter index to PriorPair
 */
inline std::unordered_map<int, PriorPair> make_prior_map(const Priors& p)
{
    std::unordered_map<int, PriorPair> m;

    m[0]  = {p.convNormCenter,           p.convNormWidth};
    m[2]  = {p.zenithCorrectionCenter,   p.zenithCorrectionWidth};
    m[3]  = {p.kaonLossesCenter,         p.kaonLossesWidth};
    m[4]  = {p.hadronicHEkpCenter,       p.hadronicHEkpWidth};
    m[5]  = {p.hadronicHEkmCenter,       p.hadronicHEkmWidth};
    m[6]  = {p.hadronicVHE1pipCenter,    p.hadronicVHE1pipWidth};
    m[7]  = {p.hadronicVHE1pimCenter,    p.hadronicVHE1pimWidth};
    m[8]  = {p.hadronicVHE3kpCenter,     p.hadronicVHE3kpWidth};
    m[9]  = {p.hadronicVHE3kmCenter,     p.hadronicVHE3kmWidth};
    m[10] = {p.hadronicVHE3pipCenter,    p.hadronicVHE3pipWidth};
    m[11] = {p.hadronicVHE3pimCenter,    p.hadronicVHE3pimWidth};
    m[12] = {p.hadronicVHE3pCenter,      p.hadronicVHE3pWidth};
    m[13] = {p.hadronicVHE3nCenter,      p.hadronicVHE3nWidth};
    m[14] = {p.cosmicRay1Center,         p.cosmicRay1Width};
    m[15] = {p.cosmicRay2Center,         p.cosmicRay2Width};
    m[16] = {p.cosmicRay3Center,         p.cosmicRay3Width};
    m[17] = {p.cosmicRay4Center,         p.cosmicRay4Width};
    m[18] = {p.cosmicRay5Center,         p.cosmicRay5Width};
    m[19] = {p.cosmicRay6Center,         p.cosmicRay6Width};
    m[20] = {p.icegrad0Center,           p.icegrad0Width};
    m[21] = {p.icegrad1Center,           p.icegrad1Width};
    m[22] = {p.icegrad2Center,           p.icegrad2Width};
    m[23] = {p.icegrad3Center,           p.icegrad3Width};
    m[24] = {p.icegrad4Center,           p.icegrad4Width};
    m[25] = {p.icegrad5Center,           p.icegrad5Width};
    m[26] = {p.icegrad6Center,           p.icegrad6Width};
    m[27] = {p.icegrad7Center,           p.icegrad7Width};
    m[28] = {p.icegrad8Center,           p.icegrad8Width};

    m[29] = {p.domEfficiencyCenter,      p.domEfficiencyWidth};
    m[30] = {p.holeiceForwardCenter,     p.holeiceForwardWidth};

    m[31] = {p.astroNormCenter,          p.astroNormWidth};
    m[32] = {p.astroDeltaGammaCenter,    p.astroDeltaGammaWidth};
    m[33] = {p.astroDeltaGammaSecCenter, p.astroDeltaGammaSecWidth};
    m[36] = {p.nuxsCenter,               p.nuxsWidth};
    m[37] = {p.nubarxsCenter,            p.nubarxsWidth};

    return m;
}

/**
 * @brief Objective function for Diver optimizer.
 *
 * Evaluates the likelihood at the given parameter point. Only the free
 * parameters (from free_indices()) are varied; others remain at seed values.
 *
 * @param params Array of free parameter values from Diver
 * @param D Dimension (number of free parameters)
 * @param fcall Function call counter (incremented by this function)
 * @param quit Set to true to signal early termination (not used here)
 * @param validvector Whether params are within bounds
 * @param context Pointer to GollumFit shared_ptr
 * @return Negative log-likelihood value, or BAD_OBJ if invalid
 */
inline double objective(
    double params[],
    const int D,
    int &fcall,
    bool &quit,
    const bool validvector,
    void*& context
){
    ++fcall;
    ++n_eval;
    quit = false;

    if (!validvector) {
        ++n_bad;
        return BAD_OBJ;
    }

    auto* gf_ptr = static_cast<std::shared_ptr<GollumFit>*>(context);
    std::shared_ptr<GollumFit> gf = *gf_ptr;

    const auto& idx = current_free_indices;

    FitParameters fp = gf->GetFitParametersSeed()[0];
    std::vector<double> full = gf->ConvertFitParameters(fp);

    for (int j = 0; j < D; ++j)
        full[idx[j]] = params[j];

    fp = gf->ConvertVecToFitParameters(full);

    double llh = gf->EvalLLH(fp, true);

    if (!std::isfinite(llh) || llh >= BAD_OBJ) {
        ++n_bad;
        return BAD_OBJ;
    }

    ++n_good;
    return llh;
}

/**
 * @brief Run the Diver differential evolution global optimizer.
 *
 * Performs global optimization over the free parameters to find the
 * best-fit region before local refinement with MinLLH. Diver parameters
 * (maxgen, NP, F, Cr, seed) are read from SteeringParams.
 *
 * Initial population is sampled from truncated Gaussian priors.
 *
 * @param gf GollumFit instance (must have likelihood problem constructed)
 * @param seed_params Initial parameter values (used for fixed parameters)
 * @param bounds Parameter bounds for optimization
 * @param priors Prior distributions for initial population sampling
 * @param path Output directory for Diver diagnostic files
 * @return FitResult containing best-fit parameters and likelihood
 */
inline FitResult Dive(
    std::shared_ptr<GollumFit> gf,
    FitParameters& seed_params,
    const FitParametersBound& bounds,
    const Priors& priors,
    const std::string& path
){
    FitResult result;

    // Get diver parameters from steering
    const auto& steer = gf->GetSteeringParams();
    const int maxgen = steer.diver_maxgen;
    const int NP = steer.diver_NP;
    const double F_in = steer.diver_F;
    const double Cr_in = steer.diver_Cr;
    const int diver_seed = steer.diver_seed;

    // Get free indices (custom if provided, otherwise default 9)
    const auto idx = get_free_indices(steer.diver_free_indices);
    current_free_indices = idx;  // Set for objective function
    const int D = static_cast<int>(idx.size());

    std::cout << "[DIVER] Using " << D << " free parameters: ";
    for (int i : idx) std::cout << i << " ";
    std::cout << std::endl;

    /* ---------- full bounds (must match FitParameters ordering) ---------- */
    std::vector<std::pair<double,double>> full_bounds = {
        {bounds.convNormMin, bounds.convNormMax},                             // 0
        {bounds.promptNormMin, bounds.promptNormMax},                         // 1
        {bounds.zenithCorrectionMin, bounds.zenithCorrectionMax},             // 2
        {bounds.kaonLossesMin, bounds.kaonLossesMax},                         // 3
        {bounds.hadronicHEkpMin, bounds.hadronicHEkpMax},                     // 4
        {bounds.hadronicHEkmMin, bounds.hadronicHEkmMax},                     // 5
        {bounds.hadronicVHE1pipMin, bounds.hadronicVHE1pipMax},               // 6
        {bounds.hadronicVHE1pimMin, bounds.hadronicVHE1pimMax},               // 7
        {bounds.hadronicVHE3kpMin, bounds.hadronicVHE3kpMax},                 // 8
        {bounds.hadronicVHE3kmMin, bounds.hadronicVHE3kmMax},                 // 9
        {bounds.hadronicVHE3pipMin, bounds.hadronicVHE3pipMax},               // 10
        {bounds.hadronicVHE3pimMin, bounds.hadronicVHE3pimMax},               // 11
        {bounds.hadronicVHE3pMin, bounds.hadronicVHE3pMax},                   // 12
        {bounds.hadronicVHE3nMin, bounds.hadronicVHE3nMax},                   // 13
        {bounds.cosmicRay1Min, bounds.cosmicRay1Max},                         // 14
        {bounds.cosmicRay2Min, bounds.cosmicRay2Max},                         // 15
        {bounds.cosmicRay3Min, bounds.cosmicRay3Max},                         // 16
        {bounds.cosmicRay4Min, bounds.cosmicRay4Max},                         // 17
        {bounds.cosmicRay5Min, bounds.cosmicRay5Max},                         // 18
        {bounds.cosmicRay6Min, bounds.cosmicRay6Max},                         // 19
        {bounds.icegrad0Min, bounds.icegrad0Max},                             // 20
        {bounds.icegrad1Min, bounds.icegrad1Max},                             // 21
        {bounds.icegrad2Min, bounds.icegrad2Max},                             // 22
        {bounds.icegrad3Min, bounds.icegrad3Max},                             // 23
        {bounds.icegrad4Min, bounds.icegrad4Max},                             // 24
        {bounds.icegrad5Min, bounds.icegrad5Max},                             // 25
        {bounds.icegrad6Min, bounds.icegrad6Max},                             // 26
        {bounds.icegrad7Min, bounds.icegrad7Max},                             // 27
        {bounds.icegrad8Min, bounds.icegrad8Max},                             // 28
        {bounds.domEfficiencyMin, bounds.domEfficiencyMax},                   // 29
        {bounds.holeiceForwardMin, bounds.holeiceForwardMax},                 // 30
        {bounds.astroNormMin, bounds.astroNormMax},                           // 31
        {bounds.astroDeltaGammaMin, bounds.astroDeltaGammaMax},               // 32
        {bounds.astroDeltaGammaSecMin, bounds.astroDeltaGammaSecMax},         // 33
        {bounds.astroPivotMin, bounds.astroPivotMax},                         // 34
        {bounds.NeutrinoAntineutrinoRatioMin, bounds.NeutrinoAntineutrinoRatioMax}, // 35
        {bounds.nuxsMin, bounds.nuxsMax},                                     // 36
        {bounds.nubarxsMin, bounds.nubarxsMax}                                // 37
    };

    std::vector<double> lower(D), upper(D);
    for (int j = 0; j < D; ++j) {
        lower[j] = full_bounds[idx[j]].first;
        upper[j] = full_bounds[idx[j]].second;
    }

    /* ---------- initial guesses (layout required by Diver example) ---------- */
    const auto prior_map = make_prior_map(priors);

    std::mt19937 rng(12345);

    auto trunc_gauss = [&](double mu, double sig, double lo, double hi){
        std::normal_distribution<double> n(mu, sig);
        double x;
        do { x = n(rng); } while (x < lo || x > hi);
        return x;
    };

    // We include nGuesses = NP guesses (one per slot). You can reduce this if desired.
    const int nGuesses = NP;

    // Diver expects 2D array flattened with FIRST index = parameter, SECOND = individual
    // i.e. initial_guesses[param + D*individual] is NOT what Diver expects.
    // It expects: initial_guesses[ param_index * nGuesses + individual_index ].
    std::vector<double> initial_guesses(D * nGuesses, 0.0);

    std::cout << "\n[DIVER INITIAL GUESSES] (printing first individual only)\n";
    for (int j = 0; j < D; ++j) {
        int i = idx[j];
        double lo = lower[j], hi = upper[j];

        double v0;
        if (i == 34) {
            v0 = std::uniform_real_distribution<double>(lo, hi)(rng);
        } else {
            auto it = prior_map.find(i);
            if (it == prior_map.end()) v0 = std::uniform_real_distribution<double>(lo, hi)(rng);
            else v0 = trunc_gauss(it->second.mu, it->second.sig, lo, hi);
        }
        std::cout << "param idx " << i << " = " << v0 << "\n";
    }

    // Fill all individuals
    for (int g = 0; g < nGuesses; ++g) {
        for (int j = 0; j < D; ++j) {
            int i = idx[j];
            double lo = lower[j], hi = upper[j];

            double v;
            if (i == 34) {
                v = std::uniform_real_distribution<double>(lo, hi)(rng);
            } else {
                auto it = prior_map.find(i);
                if (it == prior_map.end()) v = std::uniform_real_distribution<double>(lo, hi)(rng);
                else v = trunc_gauss(it->second.mu, it->second.sig, lo, hi);
            }

            // param-major layout:
            initial_guesses[j * nGuesses + g] = v;
        }
    }

    /* ---------- context pointer ---------- */
    std::shared_ptr<GollumFit>* ctx_ptr = &gf;
    void* ctx = static_cast<void*>(ctx_ptr);

    /* ---------- Diver best-fit outputs ---------- */
    const int nDerived = 0;
    std::vector<double> bestFitParams(D, 0.0);
    std::vector<double> bestFitDerived(1, 0.0); // dummy; unused if nDerived==0

    const int nDiscrete = 0;
    const int discrete_dummy[1] = {0}; // dummy; unused if nDiscrete==0
    const bool partitionDiscrete = false;

    const int nF = 1;
    const double F_arr[nF] = {F_in};
    const double Cr = Cr_in;

    const double lambda = 0.6;    
    const bool current = false;
    const bool expon = false;
    const int bndry = 3;
    const bool jDE = true;
    const bool lambdajDE = true;
    const double convthresh = 1e-3;
    const int convsteps = 10;
    const bool removeDuplicates = true;
    const int savecount = 1;
    const bool resume = false;
    const bool disableIO = false;
    const bool outputRaw = true;
    const bool outputSam = false;
    const int init_pop_strategy = 2;
    const bool discard_unfit_points = false;
    const int max_init_attempts = 100000;
    const double max_acceptable_val = BAD_OBJ;
    const int seed = diver_seed;
    const int verbose = 1;

    // IMPORTANT: cdiver path argument expects const char[]
    // so we pass path.c_str()
    double llh = cdiver(
        objective,
        D,
        lower.data(),
        upper.data(),
        path.c_str(),
        nDerived,
        bestFitParams.data(),
        bestFitDerived.data(),
        nDiscrete,
        (nDiscrete > 0 ? discrete_dummy : nullptr),
        partitionDiscrete,
        maxgen,
        NP,
        nF,
        F_arr,
        Cr,
        lambda,
        current,
        expon,
        bndry,
        jDE,
        lambdajDE,
        convthresh,
        convsteps,
        removeDuplicates,
        savecount,
        resume,
        disableIO,
        outputRaw,
        outputSam,
        init_pop_strategy,
        nGuesses,
        (nGuesses > 0 ? initial_guesses.data() : nullptr),
        discard_unfit_points,
        max_init_attempts,
        max_acceptable_val,
        seed,
        ctx,
        verbose
    );

    std::cout << "\n[DIVER DIAGNOSTICS]\n";
    std::cout << "Total evals : " << n_eval << "\n";
    std::cout << "Good evals  : " << n_good << "\n";
    std::cout << "Bad evals   : " << n_bad << "\n";
    std::cout << "Good frac   : " << (n_eval ? double(n_good)/double(n_eval) : 0.0) << "\n";

    std::cout << "\n[DIVER BEST-FIT FREE PARAMS]\n";
    for (int j = 0; j < D; ++j) {
        std::cout << "  idx " << idx[j] << " = " << bestFitParams[j] << "\n";
    }

    // Insert best-fit free params into full vector (seed for fixed params)
    FitParameters fp0 = gf->GetFitParametersSeed()[0];
    std::vector<double> full_best = gf->ConvertFitParameters(fp0);
    for (int j = 0; j < D; ++j) {
        full_best[idx[j]] = bestFitParams[j];
    }

    result.params = gf->ConvertVecToFitParameters(full_best);
    result.likelihood = llh;
    return result;
}

} // namespace diver
} // namespace gollumfit

#endif
