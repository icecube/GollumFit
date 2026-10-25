#ifndef GOLLUMFIT_SYSTEMATICSCACHE_H_
#define GOLLUMFIT_SYSTEMATICSCACHE_H_

#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <LeptonWeighter/ParticleType.h>
#include <photospline/splinetable.h>

#include "Event.h"
#include "GollumEnumDefinitions.h"
#include "GollumTools.h"
#include "ReducedSplines.h"

/**
* @file SystematicsCache.h
* @brief Per-event side table of reduced detector-systematics splines.
*
* For every simulated event and every flux component this stores the reduced
* (1-D) coefficients of the DOM-efficiency, hole-ice and attenuation splines,
* plus the log-rate at the nominal DOM-efficiency / hole-ice values, so that a
* likelihood evaluation only has to evaluate a 1-D spline per event and factor.
* The table is built once (lazily) from the loaded spline maps and the event
* containers registered by GollumFit.
*/

namespace gollumfit {

class SystematicsCache {
public:
    using DOMMapType = std::map<std::pair<FluxComponent, Topology>, std::shared_ptr<photospline::splinetable<>>>;
    using HoleIceMapType = DOMMapType;
    using AttenuationMapType = std::map<std::pair<FluxComponent, LW::ParticleType>, std::shared_ptr<photospline::splinetable<>>>;

    static constexpr unsigned int kComponents = 3;
    /// Nominal detector parameters of the MC (the committed weighters hard-code these).
    static constexpr double kDOMEfficiencyNominal = 1.27;
    static constexpr double kHoleIceNominal = -1.0;
    /// Marks a factor that is identically 1 (component not corrected).
    static constexpr int8_t kTableNone = -1;
    /// Marks a factor that is identically 0 (event outside the spline support).
    static constexpr int8_t kTableZero = -2;
    /// Sentinel for "nominal rate is zero: ignore the event" (as in the setCache helpers).
    static constexpr float kNominalSentinel = -std::numeric_limits<float>::max();

    static constexpr FluxComponent kComponentOf[kComponents] = {FluxComponent::atmConv, FluxComponent::atmPrompt, FluxComponent::diffuseAstro};

    /// Per-event bookkeeping (the reduced coefficients live in a separate flat array).
    struct EventInfo {
        SplineReal domNominal[kComponents];  ///< log10-rate at the nominal DOM efficiency, or kNominalSentinel
        SplineReal holeNominal[kComponents]; ///< log10-rate at the nominal hole-ice parameter, or kNominalSentinel
        float logPrimaryEnergy;         ///< ln(primaryEnergy), for the astrophysical tilt
        int8_t domTable[kComponents];   ///< index into domTables, kTableNone or kTableZero
        int8_t holeTable[kComponents];
        int8_t attTable[kComponents];
        bool isNeutrino;                ///< uses the nuxs (true) or nubarxs (false) attenuation scale
        bool isAntiparticle;            ///< (int)primaryType < 0
    };

    /// Offsets (in floats) of the reduced coefficients inside one event row.
    struct Layout {
        unsigned int nkDom = 0, nkHole = 0, nkAtt = 0;
        unsigned int domOffset[kComponents], holeOffset[kComponents], attOffset[kComponents];
        unsigned int stride = 0;
    };

    SystematicsCache(const DOMMapType& dom, const HoleIceMapType& hole, const AttenuationMapType& att)
        : domMap_(dom), holeMap_(hole), attMap_(att) {
        for (const auto& kv : domMap_) {
            domKeys_.push_back(kv.first);
            domTables_.emplace_back(*kv.second);
            layout_.nkDom = std::max(layout_.nkDom, domTables_.back().ncoeffs());
        }
        for (const auto& kv : holeMap_) {
            holeKeys_.push_back(kv.first);
            holeTables_.emplace_back(*kv.second);
            layout_.nkHole = std::max(layout_.nkHole, holeTables_.back().ncoeffs());
        }
        for (const auto& kv : attMap_) {
            attKeys_.push_back(kv.first);
            attTables_.emplace_back(*kv.second);
            layout_.nkAtt = std::max(layout_.nkAtt, attTables_.back().ncoeffs());
        }
        if (domTables_.size() > 100 || holeTables_.size() > 100 || attTables_.size() > 100)
            throw std::runtime_error("SystematicsCache: too many spline tables");
        unsigned int off = 0;
        for (unsigned int c = 0; c < kComponents; ++c) { layout_.domOffset[c] = off; off += layout_.nkDom; }
        for (unsigned int c = 0; c < kComponents; ++c) { layout_.holeOffset[c] = off; off += layout_.nkHole; }
        for (unsigned int c = 0; c < kComponents; ++c) { layout_.attOffset[c] = off; off += layout_.nkAtt; }
        layout_.stride = off;
        for (unsigned int c = 0; c < kComponents; ++c) {
            domComponentPresent_[c] = componentInMap(domMap_, kComponentOf[c]);
            holeComponentPresent_[c] = componentInMap(holeMap_, kComponentOf[c]);
        }
    }

    /// Register a contiguous range of events and build their rows.
    void addEvents(const Event* begin, size_t n) {
        Range r;
        r.begin = begin;
        r.n = n;
        r.offset = infos_.size();
        ranges_.push_back(r);
        infos_.resize(infos_.size() + n);
        rows_.resize(rows_.size() + static_cast<size_t>(n) * layout_.stride);
        for (size_t i = 0; i < n; ++i)
            build(begin[i], infos_[r.offset + i], &rows_[(r.offset + i) * layout_.stride]);
    }

    size_t size() const { return infos_.size(); }
    const Layout& layout() const { return layout_; }
    const std::vector<ReducedAxisSpline>& domTables() const { return domTables_; }
    const std::vector<ReducedAxisSpline>& holeTables() const { return holeTables_; }
    const std::vector<ReducedAxisSpline>& attTables() const { return attTables_; }

    /// Index of a registered event (by address).
    inline size_t index(const Event& e) const {
        const Event* p = &e;
        for (const Range& r : ranges_)
            if (p >= r.begin && p < r.begin + r.n)
                return r.offset + static_cast<size_t>(p - r.begin);
        throw std::runtime_error("SystematicsCache: event is not part of a registered simulation set");
    }
    inline const EventInfo& info(size_t idx) const { return infos_[idx]; }
    inline const SplineReal* row(size_t idx) const { return &rows_[idx * layout_.stride]; }

private:
    struct Range { const Event* begin; size_t n; size_t offset; };

    DOMMapType domMap_;
    HoleIceMapType holeMap_;
    AttenuationMapType attMap_;
    std::vector<std::pair<FluxComponent, Topology>> domKeys_, holeKeys_;
    std::vector<std::pair<FluxComponent, LW::ParticleType>> attKeys_;
    std::vector<ReducedAxisSpline> domTables_, holeTables_, attTables_;
    bool domComponentPresent_[kComponents], holeComponentPresent_[kComponents];
    Layout layout_;
    std::vector<Range> ranges_;
    std::vector<EventInfo> infos_;
    std::vector<SplineReal> rows_;

    template <typename Map>
    static bool componentInMap(const Map& m, FluxComponent c) {
        for (const auto& kv : m)
            if (kv.first.first == c)
                return true;
        return false;
    }
    template <typename Key>
    static int8_t findKey(const std::vector<Key>& keys, const Key& k) {
        for (size_t i = 0; i < keys.size(); ++i)
            if (keys[i] == k)
                return static_cast<int8_t>(i);
        return kTableNone;
    }

    /// Reduce one detector spline (DOM efficiency or hole ice) for one event.
    /// Mirrors DOMEffWeighter / holeiceWeighter / DOMEfficiencySetter semantics:
    /// factor 1 if the component is not corrected, an error if the component
    /// is corrected but the (component,topology) table is missing, and the
    /// event is ignored (weight 0) if the nominal rate is 0 or the event lies
    /// outside the table support.
    static void reduceDetector(const Event& e, FluxComponent comp, bool componentPresent,
                               const std::vector<std::pair<FluxComponent, Topology>>& keys,
                               const std::vector<ReducedAxisSpline>& tables, double nominal,
                               const char* what, int8_t& tableIdx, SplineReal& nominalRate, SplineReal* out) {
        if (!componentPresent) {
            tableIdx = kTableNone;
            nominalRate = 0.f;
            return;
        }
        tableIdx = findKey(keys, std::make_pair(comp, static_cast<Topology>(e.topology)));
        if (tableIdx == kTableNone)
            throw std::runtime_error(std::string(what) + " correction for " + GetFluxComponentName(comp) +
                                     " flux component and for topology " + GetTopologyName(static_cast<Topology>(e.topology)) + " not found.");
        const ReducedAxisSpline& t = tables[tableIdx];
        if (!t.reduce(std::log10(e.energy), std::cos(e.zenith), out)) {
            tableIdx = kTableZero;
            nominalRate = kNominalSentinel;
            return;
        }
        ReducedAxisSpline::Basis bs = t.basis(nominal);
        SplineReal r = bs.valid ? ReducedAxisSpline::eval(out, bs) : 0.;
        if (!std::isfinite(r))
            throw std::runtime_error(std::string(what) + " spline evaluates to a non-finite value at the nominal point.");
        nominalRate = (r == 0.f) ? kNominalSentinel : r;
    }

    void build(const Event& e, EventInfo& info, SplineReal* row) const {
        info.logPrimaryEnergy = static_cast<float>(std::log(static_cast<double>(e.primaryEnergy)));
        info.isNeutrino = tools::isNeutrino(e.primaryType);
        info.isAntiparticle = static_cast<int>(e.primaryType) < 0;
        const bool attenuationSkipped = std::cos(e.primaryZenith) > 0.1;
        for (unsigned int c = 0; c < kComponents; ++c) {
            reduceDetector(e, kComponentOf[c], domComponentPresent_[c], domKeys_, domTables_, kDOMEfficiencyNominal,
                           "DOM efficiency", info.domTable[c], info.domNominal[c], row + layout_.domOffset[c]);
            reduceDetector(e, kComponentOf[c], holeComponentPresent_[c], holeKeys_, holeTables_, kHoleIceNominal,
                           "HoleIce", info.holeTable[c], info.holeNominal[c], row + layout_.holeOffset[c]);
            // attenuation: factor 1 when no table matches or the event is down-going
            int8_t idx = findKey(attKeys_, std::make_pair(kComponentOf[c], e.primaryType));
            if (idx == kTableNone || attenuationSkipped) {
                info.attTable[c] = kTableNone;
                continue;
            }
            const ReducedAxisSpline& t = attTables_[idx];
            if (!t.reduce(std::log10(e.primaryEnergy), std::cos(e.primaryZenith), row + layout_.attOffset[c]))
                info.attTable[c] = kTableZero; // photospline evaluates to 0 outside the support
            else
                info.attTable[c] = idx;
        }
    }
};

} // namespace gollumfit

#endif /* GOLLUMFIT_SYSTEMATICSCACHE_H_ */
