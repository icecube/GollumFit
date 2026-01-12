#ifndef ANALYSISWEIGHTING_H_
#define ANALYSISWEIGHTING_H_

#include <cassert>
#include <type_traits>
#include <limits>

#include <PhysTools/likelihood/weighting.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/autodiff.h>
#include <nuSQuIDS/marray.h>
#include <nuSQuIDS/tools.h>
#include <photospline/splinetable.h>
#include <math.h>

#include "GollumParameters.h"
#include "GollumTools.h"
#include "utils.h"

/**
* @file analysisWeighting.h
* @brief Definitions of the individual weighters used in the likelihood evaluation
* and the weighter maker to perform the actual weighting of events.
*/

using namespace photospline;
using namespace phys_tools::autodiff;

//================================================================================
// FUNDAMENTAL WEIGHTERS
//================================================================================

/**
* @brief A template struct for weighting values cached in events.
*
* This weighter extracts values of a specified member from Event objects,
* casting them to a specified result type.
*
* @tparam T The type to which the cached value will be cast when retrieved.
* @tparam Event The type of the event object from which the cached value will be extracted.
* @tparam U The type of the cached value within the Event object.
*/
template<typename T, typename Event, typename U>
struct cachedValueWeighter : public phys_tools::GenericWeighter<cachedValueWeighter<T,Event,U>>{
    private:
      /**
      * @brief Pointer to the member of Event that holds the cached value.
      */
      U Event::* cachedPtr;
    public:
      /**
      * @brief The type to which the cached value will be cast.
      */
      using result_type=T;
      /**
      * @brief Constructor
      * @param ptr Pointer to the member of Event that holds the cached value.
      */
      cachedValueWeighter(U Event::* ptr):cachedPtr(ptr){}
      /**
      * @brief Retrieves and casts the cached value from an Event object.
      *
      * @param e The Event object from which to extract the cached value.
      * @return The cached value cast to the specified result type.
      */
      result_type operator()(const Event& e) const{
          return(result_type(e.*cachedPtr));
      }
};

//================================================================================
// ATMOSPHERIC FLUX WEIGHTERS
//================================================================================

//Weight particles and antiparticles differently
/**
* @brief A template struct to assign weights to particles and antiparticles differently.
*
* The weighter uses a balance factor to determine the weights for particles and antiparticles.
*
* @tparam Event The type of the event object.
* @tparam T The numeric type of the weight and balance factor.
*/
template<typename Event, typename T>
struct antiparticleWeighter : public phys_tools::GenericWeighter<antiparticleWeighter<Event,T>>{
    private:
      /**
      * @brief The balance factor for weighting particles and antiparticles.
      *
      * A balance factor of 1 assigns equal weight to particles and antiparticles.
      * A balance factor of 0 assigns zero weight to particles and double weight to antiparticles.
      * A balance factor of 2 assigns double weight to particles and zero weight to antiparticles.
      */
      T balance;
    public:
      /**
      * @brief The result type of the weight computation.
      */
      using result_type=T;

      /**
      * @brief Constructor
      *
      * @param b The balance factor to be used for weighting.
      */
      antiparticleWeighter(T b):
          balance(b){}

      /**
      * @brief Calculates the weight for an event based on its primary type.
      *
      * If the event's primary type indicates an antiparticle (negative), the balance factor is used as the weight.
      * Otherwise, the weight is calculated as 2 - balance, effectively inverting the balance for particles.
      *
      * @param e The event object for which to calculate the weight.
      * @return The calculated weight as the result type T.
      */
      result_type operator()(const Event& e) const{
          return((int)e.primaryType<0 ? balance : 2-balance);
      }
};


// Tilt a spectrum by an incremental powerlaw index about a precomputed median energy
/**
* @brief A weighter struct template to tilt an event spectrum by incremental power law indices.
* 
* This weighter applies a tilt to a spectrum using two different power-law indices, which
* are applied below and above a specified median energy threshold. The tilt changes the spectral index
* of the event energies, effectively accentuating or attenuating the contribution of events
* with energies higher or lower than the median energy.
* 
* @tparam Event The type of the event object which must contain a primaryEnergy member.
* @tparam T The numeric type for energy and weight calculations.
*/
template<typename Event, typename T>
struct brokenpowerlawTiltWeighter : public phys_tools::GenericWeighter<brokenpowerlawTiltWeighter<Event,T>>{
    private:
      T medianlog10Energy; ///>The (base 10) log of the median energy about which the spectrum is tilted.
      T deltaIndex1; ///>The incremental change in spectral index for energies below the median energy.
      T deltaIndex2; ///>The incremental change in spectral index for energies above the median energy.
    public:
        using result_type=T;
        /**
        * @brief Constructor
        * 
        * @param me The logarithm (base 10) of the median energy.
        * @param dg1 The change in spectral index for energies below the median energy.
        * @param dg2 The change in spectral index for energies above the median energy.
        */
        brokenpowerlawTiltWeighter(T me, T dg1, T dg2):
            medianlog10Energy(me),deltaIndex1(dg1),deltaIndex2(dg2){}

        /**
        * @brief Calculates the weight for an event according to the broken power law tilt.
        * 
        * The weight is determined by the event's primary energy relative to the median energy.
        * If the primary energy is greater than the median, the weight is calculated using
        * deltaIndex2; otherwise, deltaIndex1 is used.
        * 
        * The weight is calculated by exponentiating; for example:
        * @code{.cpp}
        * weight = pow(e.primaryEnergy/medianEnergy,-deltaIndex2)
        * @endcode
        * 
        * 
        * @param e The event object for which to calculate the weight.
        * @return The calculated weight as the result type T.
        */
        result_type operator()(const Event& e) const{
            result_type medianEnergy = pow(10.,medianlog10Energy);
            result_type weight = (e.primaryEnergy>medianEnergy) ? pow(e.primaryEnergy/medianEnergy,-deltaIndex2) : pow(e.primaryEnergy/medianEnergy,-deltaIndex1);
            return(weight);
        }
};

//================================================================================
// CONV FLUX WEIGHTER (PENALTY NEGATIVE)
//================================================================================

/**
* @brief A weighter to calculate the weighted flux for an event.
*
* The ConvFluxWeigther combines several contributions to the flux of an event,
* including hadronic contributions and cosmic ray components.
* Each type of contribution is weighted by a set of coefficients that are passed
* to the constructor.
*
* @tparam Event The type of the event object which must contain cached values for the
*             different flux components and hadronic and cosmic ray contributions.
* @tparam T The numeric type for the coefficients and the result type for the flux.
*/
template<typename Event, typename T> 
    struct ConvFluxWeigther : public phys_tools::GenericWeighter<ConvFluxWeigther<Event,T>>{
        private:
            T hekp; ///< High energy K+ (158 GeV) (DAEMONFLUX parameter - hardonic yield)
            T hekm; ///< High energy K- (158 GeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe1pip; ///< very high energy pi+ (20 TeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe1pim; ///< very high energy pi- (20 TeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3kp; ///< very high energy K+ (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3km; ///< very high energy K- (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3pip; ///< very high energy pi+ (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3pim; ///< very high energy pi- (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3p; ///< very high energy p (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T vhe3n; ///< very high energy n (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            T cr1; ///< Global Spline Fit 1 (DAEMONFLUX parameter - cosmic ray spectrum)
            T cr2; ///< Global Spline Fit 2 (DAEMONFLUX parameter - cosmic ray spectrum)
            T cr3; ///< Global Spline Fit 3 (DAEMONFLUX parameter - cosmic ray spectrum)
            T cr4; ///< Global Spline Fit 4 (DAEMONFLUX parameter - cosmic ray spectrum)
            T cr5; ///< Global Spline Fit 5 (DAEMONFLUX parameter - cosmic ray spectrum)
            T cr6; ///< Global Spline Fit 6 (DAEMONFLUX parameter - cosmic ray spectrum)
        public:
            using result_type=T;///< The result type of the flux computation.

            /**
            * @brief Constructor.
            * 
            * Assigns inputs to private struct variables.
            *
            * @param hekp_ High energy K+ (158 GeV) (DAEMONFLUX parameter - hardonic yield)
            * @param hekm_ High energy K- (158 GeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe1pip_ very high energy pi+ (20 TeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe1pim_ very high energy pi- (20 TeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3kp_ very high energy K+ (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3km_ very high energy K- (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3pip_ very high energy pi+ (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3pim_ very high energy pi- (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3p_ very high energy p (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param vhe3n_ very high energy n (2 PeV) (DAEMONFLUX parameter - hardonic yield)
            * @param cr1_ Global Spline Fit 1 (DAEMONFLUX parameter - cosmic ray spectrum)
            * @param cr2_ Global Spline Fit 2 (DAEMONFLUX parameter - cosmic ray spectrum)
            * @param cr3_ Global Spline Fit 3 (DAEMONFLUX parameter - cosmic ray spectrum)
            * @param cr4_ Global Spline Fit 4 (DAEMONFLUX parameter - cosmic ray spectrum)
            * @param cr5_ Global Spline Fit 5 (DAEMONFLUX parameter - cosmic ray spectrum)
            * @param cr6_ Global Spline Fit 6 (DAEMONFLUX parameter - cosmic ray spectrum)
            */
            ConvFluxWeigther(T hekp_, T hekm_, 
                             T vhe1pip_, T vhe1pim_, 
                             T vhe3kp_, T vhe3km_, 
                             T vhe3pip_, T vhe3pim_, 
                             T vhe3p_, T vhe3n_, 
                             T cr1_, T cr2_, 
                             T cr3_, T cr4_, 
                             T cr5_, T cr6_):
                                 hekp(hekp_),  hekm(hekm_), 
                                 vhe1pip(vhe1pip_),  vhe1pim(vhe1pim_), 
                                 vhe3kp(vhe3kp_),  vhe3km(vhe3km_), 
                                 vhe3pip(vhe3pip_),  vhe3pim(vhe3pim_), 
                                 vhe3p(vhe3p_),  vhe3n(vhe3n_), 
                                 cr1(cr1_),  cr2(cr2_), 
                                 cr3(cr3_),  cr4(cr4_), 
                                 cr5(cr5_),  cr6(cr6_) {}

            /**
            * @brief Calculates the weighted flux for an event.
            *
            * This operator sums the conventional flux with the hadronic and cosmic ray
            * contributions, each multiplied by their respective coefficients. If the resulting
            * flux is negative (which is unphysical), it returns the maximum representable float
            * value as a fallback.
            *
            * @param e The event object for which to calculate the weighted flux.
            * @return The weighted flux or max float value if the flux is unphysical.
            */
            result_type operator()(const Event& e) const{
              result_type convFlux = e.cachedConvWeight;
              result_type hadronic = hekp    * e.cachedHadronicHEkp +
                                      hekm    * e.cachedHadronicHEkm +
                                      vhe1pip * e.cachedHadronicVHE1pip +
                                      vhe1pim * e.cachedHadronicVHE1pim +
                                      vhe3kp  * e.cachedHadronicVHE3kp +
                                      vhe3km  * e.cachedHadronicVHE3km +
                                      vhe3pip * e.cachedHadronicVHE3pip +
                                      vhe3pim * e.cachedHadronicVHE3pim +
                                      vhe3p   * e.cachedHadronicVHE3p +
                                      vhe3n   * e.cachedHadronicVHE3n ;
              result_type cr       = cr1     * e.cachedCosmicRay1 +
                                      cr2     * e.cachedCosmicRay2 +
                                      cr3     * e.cachedCosmicRay3 +
                                      cr4     * e.cachedCosmicRay4 +
                                      cr5     * e.cachedCosmicRay5 +
                                      cr6     * e.cachedCosmicRay6 ;
              result_type flux = convFlux + cr + hadronic;
              if ( flux<0 ) {
                  // std::cout << "Unphysical region " << flux << std::endl;
                  // std::cout << e << std::endl;
                  // std::cout << hekp << " " << hekm << " " 
                  //           << vhe1pip << " " << vhe1pim << " " 
                  //           << vhe3kp << " " << vhe3km << " " 
                  //           << vhe3pip << " " << vhe3pim << " " 
                  //           << vhe3p << " " << vhe3n << " " 
                  //           << cr1 << " " << cr2 << " " 
                  //           << cr3 << " " << cr4 << " " 
                  //           << cr5 << " " << cr6 << std::endl;
                  return result_type(std::numeric_limits<float>::max());
              }
              return flux;
            }
    };

//================================================================================
// DOM EFFICIENCY WEIGHTER
//================================================================================

using DOMMapType=std::map<std::pair<gollumfit::FluxComponent,gollumfit::Topology>,std::shared_ptr<splinetable<>>>;


/**
* @brief Weighter struct to calculate the digital optical module (DOM) efficiency weights for events.
*
* This struct uses a map of DOM efficiencies to calculate a correction factor based on the event's properties
* and predefined efficiency of the DOM. The correction factor is derived from a spline table and is applied to
* the cached efficiency value for the event. The weight is the power of 10 of the difference between the
* correction factor and the cached efficiency.
*
* @tparam Event The event type.
* @tparam DataType The numeric type for the DOM efficiency values (e.g., float, double). This type should
*                support operations with doubles for the weight calculation.
*/
template<typename Event, typename DataType>
    struct DOMEffWeighter : public phys_tools::GenericWeighter<DOMEffWeighter<Event,DataType>>{
        private:
            const DOMMapType& dom_efficiency_map_;
            const DataType dom_efficiency_;
            const gollumfit::FluxComponent flux_component_;
            const bool enforce_needed;
        private:
            bool IsComponentInMap(const DOMMapType& mapa, gollumfit::FluxComponent componente){
              for(auto elo : mapa){
                if(elo.first.first == componente) return true;
              }
              return false;
            }
        public:

            /**
            * @brief Constructor.
            *
            * @param dom_efficiency_map_ Reference to the map associating flux component and topology pairs with
            *                            their corresponding DOM efficiency spline tables.
            * @param dom_efficiency_ The efficiency value used in the correction calculation.
            * @param flux_component_ The flux component type for which the weight will be calculated.
            */
            DOMEffWeighter(const DOMMapType& dom_efficiency_map_, DataType dom_efficiency_, gollumfit::FluxComponent flux_component_):
                dom_efficiency_map_(dom_efficiency_map_),
                dom_efficiency_(dom_efficiency_),
                flux_component_(flux_component_),
                enforce_needed((not dom_efficiency_map_.empty()) and IsComponentInMap(dom_efficiency_map_,flux_component_))
            {}

            using result_type=double;

            /**
            * @brief Calculates the DOM efficiency weight for a given event.
            *
            * If the DOM efficiency map is not empty and the flux component is found within the map, the function
            * calculates a correction factor using a spline table and returns the efficiency weight. Otherwise,
            * a default weight of 1.0 is returned.
            * 
            * The calculated DOM efficiency weight is given as `pow(10.,rate-cache)`, where `rate` is the 
            * interpolated correction from the spline table and `cache` is the cached weight of the event.
            *
            * @param e The event for which to calculate the weight.
            * @return The calculated DOM efficiency weight or 1.0 if no correction is needed or possible.
            * 
            */
            result_type operator()(const Event& e) const{
                if(dom_efficiency_map_.empty() or not enforce_needed){
                  return result_type(1.);
                }
                // double cache;
                // if(flux_component_ == gollumfit::FluxComponent::atmConv)
                //   cache = e.cachedDOMEffConv;
                // else if (flux_component_ == gollumfit::FluxComponent::atmPrompt)
                //   cache = e.cachedDOMEffPrompt;
                // else if (flux_component_ == gollumfit::FluxComponent::diffuseAstro)
                //   cache = e.cachedDOMEffAstro;
                // else
                //   throw std::runtime_error("DOM efficiency correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");

                double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_};

                auto domefficiencycorrection = dom_efficiency_map_.find(std::make_pair(flux_component_,static_cast<gollumfit::Topology>(e.topology)));
                if(domefficiencycorrection == dom_efficiency_map_.end())
                    throw std::runtime_error("DOM efficiency correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");
                double rate=(*(*domefficiencycorrection).second)(coordinates);

                coordinates[2]=1.27;
                double cache = (*(*domefficiencycorrection).second)(coordinates);

                if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                  // this is out of parameter space, ignore the event
                  return 0.0;
                }
                return(pow(10.,rate-cache));
            }
    };


/**
* @brief Weighter struct to calculate DOM efficiency weights suitable for automatic differentiation.
*
* This struct specialization is similar to the DOMEffWeighter for simple numeric types but is capable of
* handling @c phys_tools::autodiff::FD types, which carry both value and derivative information. This allows
* for the computation of gradients alongside the efficiency weight, which is useful for optimization algorithms
* that require derivative information.
*
* @tparam Event The event type.
* @tparam Dim The dimensionality of the automatic differentiation. This should correspond to the number of
*             variables with respect to which derivatives will be taken.
*/
template<typename Event,unsigned int Dim>
    struct DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<DOMEffWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
        private:
            const DOMMapType& dom_efficiency_map_;
            const phys_tools::autodiff::FD<Dim> dom_efficiency_;
            const gollumfit::FluxComponent flux_component_;
            unsigned int didx;
            const bool enforce_needed;
        private:
            bool IsComponentInMap(const DOMMapType& mapa, gollumfit::FluxComponent componente){
              for(auto elo : mapa){
                if(elo.first.first == componente) return true;
              }
              return false;
            }
        public:

            /**
            * @brief Constructor. Requires map of DOM efficiencies, a specific efficiency value,
            *       a flux component type, and an index indicating the variable with respect to which to differentiate.
            *
            * @param dom_efficiency_map_ Reference to the map associating flux component and topology pairs with
            *                            their corresponding DOM efficiency spline tables.
            * @param dom_efficiency_ The @c phys_tools::autodiff::FD efficiency value used in the correction calculation.
            * @param flux_component_ The flux component type for which the weight will be calculated.
            */
            DOMEffWeighter(const DOMMapType& dom_efficiency_map_, phys_tools::autodiff::FD<Dim> dom_efficiency_, gollumfit::FluxComponent flux_component_):
                dom_efficiency_map_(dom_efficiency_map_),
                dom_efficiency_(dom_efficiency_),
                flux_component_(flux_component_),
                enforce_needed((not dom_efficiency_map_.empty()) and IsComponentInMap(dom_efficiency_map_,flux_component_))
            {
                //const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(dom_efficiency_);
                for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(dom_efficiency_.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                }
            }

            /**
            * @brief Calculates the DOM efficiency weight and its derivative for a given event.
            *
            * If the DOM efficiency map is not empty and the flux component is within the map, the function
            * calculates a correction factor using a spline table, its derivative, and returns them wrapped
            * in a @c phys_tools::autodiff::FD object. Otherwise, a default weight of 1.0 is returned.
            *
            * @param e The event for which to calculate the weight.
            * @return The calculated DOM efficiency weight wrapped in an autodiff type, or 1.0 if no correction
            *         is needed or possible.
            */
            using result_type=phys_tools::autodiff::FD<Dim>;
            result_type operator()(const Event& e) const{
                if(dom_efficiency_map_.empty() or not enforce_needed){
                  return result_type(1.);
                }

                // double cache;
                // if(flux_component_ == gollumfit::FluxComponent::atmConv)
                //   cache = e.cachedDOMEffConv;
                // else if (flux_component_ == gollumfit::FluxComponent::atmPrompt)
                //   cache = e.cachedDOMEffPrompt;
                // else if (flux_component_ == gollumfit::FluxComponent::diffuseAstro)
                //   cache = e.cachedDOMEffAstro;
                // else
                //   throw std::runtime_error("DOM efficiency correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");

                double rate, derivative;
                double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_.value()};
                int centers[3];

                auto domefficiencycorrection = dom_efficiency_map_.find(std::make_pair(flux_component_,static_cast<gollumfit::Topology>(e.topology)));
                if(domefficiencycorrection == dom_efficiency_map_.end())
                    throw std::runtime_error("DOM efficiency correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");
                rate=(*(*domefficiencycorrection).second)(coordinates);

                if(not (*(*domefficiencycorrection).second).searchcenters(coordinates,centers))
                    throw std::runtime_error("Out of phase space: domeff spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(dom_efficiency_.value()) + ").");
                derivative=(*(*domefficiencycorrection).second).ndsplineeval(coordinates,centers,1<<2);
                derivative*=dom_efficiency_.derivative(didx);

                coordinates[2]=1.27;
                double cache = (*(*domefficiencycorrection).second)(coordinates);

                if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                  // this is out of parameter space, ignore the event
                  result_type r(0.0);
                  r.setDerivative(didx,0.0);
                  return r;
                }

                result_type r(rate);
                r.setDerivative(didx,derivative);
                return(pow(10.,r-cache));
            }
    };

//used for initializing the per-event dom-eff related caches
/**
* @brief Struct to initialize the per-event DOM efficiency related caches.
*
* This utility is used to precompute and store event-specific DOM efficiency values
* by querying a map of spline tables based on the event properties. Each event's cache is
* updated with a log-transformed efficiency value obtained from the corresponding spline
* for the event's topology and energy.
*
* @tparam Event The event type that must contain energy, zenith, and topology properties,
*               as well as fields for caching DOM efficiency values.
*/
template<typename Event>
    struct DOMEfficiencySetter{
        const DOMMapType& dom_efficiency_map_;///< Reference to the map of DOM efficiency spline tables.
        const double dom_efficiency_;///< The global DOM efficiency to be used for all events.

        /**
        * @brief Constructor.
        *
        * @param dom_efficiency_map_ The map of spline tables corresponding to different combinations
        *                            of flux components and topologies.
        * @param dom_efficiency_ The DOM efficiency value to be used in the spline evaluation.
        */
        DOMEfficiencySetter(const DOMMapType& dom_efficiency_map_, double dom_efficiency_):
            dom_efficiency_map_(dom_efficiency_map_),dom_efficiency_(dom_efficiency_){}

        /**
        * @brief Sets the cache of an event with the precomputed DOM efficiency value.
        *
        * The function calculates the efficiency value from the DOM efficiency map using the event's
        * energy, zenith, and the provided DOM efficiency. It then caches this value in the appropriate
        * field of the event. If the calculated cache value is zero, it is replaced with the lowest
        * possible float value to effectively ignore the event in subsequent processing.
        *
        * @param e Reference to the event whose cache will be updated.
        *
        * @throw std::runtime_error If the spline evaluation yields an efficiency value that is not finite,
        *                         or if the flux component for the DOM efficiency is not implemented.
        */
        void setCache(Event& e) const{
            assert(not dom_efficiency_map_.empty());
            double coordinates[3]={log10(e.energy),cos(e.zenith),dom_efficiency_};
            for(auto element : dom_efficiency_map_){
              if(element.first.second == static_cast<gollumfit::Topology>(e.topology)){
                double cache = (*element.second)(coordinates);
                assert(std::isfinite(cache));
                if(cache < -1.e5){
                  std::cout << e << std::endl;
                  throw std::runtime_error("Even out of parameter space.");
                }
                if(cache == 0){
                  // the cache if the log of the rate; this cases the event to be ignored.
                  cache = -std::numeric_limits<float>::max();
                }
                if(element.first.first == gollumfit::FluxComponent::atmConv)
                  e.cachedDOMEffConv = cache;
                else if(element.first.first == gollumfit::FluxComponent::atmPrompt)
                  e.cachedDOMEffPrompt = cache;
                else if(element.first.first == gollumfit::FluxComponent::diffuseAstro)
                  e.cachedDOMEffAstro= cache;
                else
                  throw std::runtime_error("Unimplemented DOMEff component");
              }
            }
        }
    };


//================================================================================
// HOLEICE WEIGHTER
//================================================================================

using HoleIceMapType=std::map<std::pair<gollumfit::FluxComponent,gollumfit::Topology>,std::shared_ptr<splinetable<>>>;

/**
* @brief Struct to weight events based on hole ice corrections using a map of spline tables.
* 
* @tparam Event The type of the event data.
* @tparam DataType The type of the data used in the hole ice correction, typically a float or double.
*/
template<typename Event, typename DataType>
    struct holeiceWeighter : public phys_tools::GenericWeighter<holeiceWeighter<Event,DataType>>{
        private:
            const HoleIceMapType& hole_ice_map_;///< Map of flux components and topologies to spline tables.
            DataType hole_ice_forward_;///< Hole ice model parameter.
            const gollumfit::FluxComponent flux_component_;///< The flux component for the correction.
            const bool enforce_needed;///< Flag to determine if correction is needed.
        private:
            /**
            * @brief Checks if the specified flux component is in the hole ice map.
            * 
            * @param mapa The hole ice map to check.
            * @param componente The flux component to look for.
            * @return true If the component is found.
            * @return false If the component is not found.
            */
            bool IsComponentInMap(const HoleIceMapType& mapa, gollumfit::FluxComponent componente){
              for(auto elo : mapa){
                if(elo.first.first == componente) return true;
              }
              return false;
            }
        public:
            /**
            * @brief Constructor. 
            * 
            * @param hole_ice_map A reference to the map of spline tables for corrections.
            * @param hole_ice_forward A hole ice model parameter (p2). 
            * @param flux_component The flux component (e.g., atmConv, atmPrompt, diffuseAstro).
            */
            holeiceWeighter(const HoleIceMapType& hole_ice_map_, DataType hole_ice_forward_, gollumfit::FluxComponent flux_component_):
                hole_ice_map_(hole_ice_map_),
                hole_ice_forward_(hole_ice_forward_),
                flux_component_(flux_component_),
                enforce_needed((not hole_ice_map_.empty()) and IsComponentInMap(hole_ice_map_,flux_component_))
            {}

            using result_type=double;///< The result type of the weight calculation.

            /**
            * @brief Apply the hole ice correction to an event.
            * 
            * If hole ice corrections are not necessary or if the map is empty, returns a weight of 1.0.
            * Otherwise, it calculates a correction based on the event's properties and the corresponding spline table.
            * 
            * @param e The event to which the correction will be applied.
            * @return The calculated weight for the event, after applying the hole ice correction.
            * @exception std::runtime_error Thrown if the correction cannot be found for the given event properties.
            */
            result_type operator()(const Event& e) const{
                if(hole_ice_map_.empty() or not enforce_needed)
                  return 1.;

                // double cache;
                // if(flux_component_ == gollumfit::FluxComponent::atmConv)
                //   cache = e.cachedHoleIceConv;
                // else if (flux_component_ == gollumfit::FluxComponent::atmPrompt)
                //   cache = e.cachedHoleIcePrompt;
                // else if (flux_component_ == gollumfit::FluxComponent::diffuseAstro)
                //   cache = e.cachedHoleIceAstro;
                // else
                //   throw std::runtime_error("HoleIce correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");

                double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_};

                auto holeicecorrection = hole_ice_map_.find(std::make_pair(flux_component_,static_cast<gollumfit::Topology>(e.topology)));
                if(holeicecorrection == hole_ice_map_.end())
                    throw std::runtime_error("HoleIce correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");
                double rate=(*(*holeicecorrection).second)(coordinates);

                coordinates[2]=-1.;
                double cache = (*(*holeicecorrection).second)(coordinates);
                
                if(cache == 0 and rate > -1.e5)
                  throw std::runtime_error("HoleIce weight out of parameter space.");

                if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                  // this is out of parameter space, ignore the event
                  return 0.0;
                }
                return(pow(10.,rate-cache));
            }
    };


/**
* @brief Weighter struct to calculate hole ice weights suitable for automatic differentiation.
*
* This struct specialization is similar to the holeiceWeighter for simple numeric types but is capable of
* handling @c phys_tools::autodiff::FD types, which carry both value and derivative information. This allows
* for the computation of gradients alongside the weight, which is useful for optimization algorithms
* that require derivative information.
*
* @tparam Event The event type.
* @tparam Dim The dimensionality of the automatic differentiation. This should correspond to the number of
*             variables with respect to which derivatives will be taken.
*/
template<typename Event,unsigned int Dim>
    struct holeiceWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<holeiceWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
        private:
            const HoleIceMapType& hole_ice_map_;///< Reference to the map containing spline tables for hole ice corrections.
            phys_tools::autodiff::FD<Dim> hole_ice_forward_;///< The forward parameter for hole ice, capable of automatic differentiation.
            const gollumfit::FluxComponent flux_component_;///< The flux component enum indicating the type of neutrino flux.
            unsigned int didx; ///< Index of the derivative in the automatic differentiation.
            const bool enforce_needed;///< Flag indicating whether the correction should be applied.
        private:
            /**
            * @brief Checks if the specified flux component is present in the hole ice map.
            * 
            * @param mapa The map containing the hole ice correction data.
            * @param componente The flux component to search for.
            * @return true if the component is found, false otherwise.
            */
            bool IsComponentInMap(const HoleIceMapType& mapa, gollumfit::FluxComponent componente){
              for(auto elo : mapa){
                if(elo.first.first == componente) return true;
              }
              return false;
            }
        public:
            /**
            * @brief Constructs a new holeiceWeighter specialized for automatic differentiation.
            * 
            * Initializes the weighter with the provided hole ice map, forward parameter, and flux component.
            * It also determines the index of the active derivative, which is required for automatic differentiation.
            * 
            * @param hole_ice_map Reference to the hole ice correction map.
            * @param hole_ice_forward Hole ice forward parameter capable of automatic differentiation.
            * @param flux_component The specific flux component for which corrections are to be applied.
            */
            holeiceWeighter(const HoleIceMapType& hole_ice_map_, phys_tools::autodiff::FD<Dim> hole_ice_forward_, gollumfit::FluxComponent flux_component_):
                hole_ice_map_(hole_ice_map_),
                hole_ice_forward_(hole_ice_forward_),
                flux_component_(flux_component_),
                enforce_needed((not hole_ice_map_.empty()) and IsComponentInMap(hole_ice_map_,flux_component_)) {
                //const unsigned int n=phys_tools::autodiff::detail::dimensionExtractor<phys_tools::autodiff::FD,Dim,double>::nVars(hole_ice_forward_);
                for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                    if(hole_ice_forward_.derivative(i)!=0){
                        didx=i;
                        break;
                    }
                }
            }

            using result_type=phys_tools::autodiff::FD<Dim>;
            /**
            * @brief Applies the hole ice correction to an event, capable of returning both value and derivatives.
            * 
            * If hole ice corrections are not necessary or if the map is empty, returns a weight of 1.0.
            * Otherwise, it calculates a correction based on the event's properties and the corresponding spline table.
            * The correction also includes the derivative information needed for automatic differentiation.
            * 
            * @param e The event to which the correction will be applied.
            * @return The calculated weight for the event, which includes derivative information.
            * @exception std::runtime_error Thrown if the correction cannot be found for the given event properties or if the event is out of the valid phase space.
            */
            result_type operator()(const Event& e) const{
                if(hole_ice_map_.empty() or not enforce_needed)
                  return result_type(1.);

                // double cache;
                // if(flux_component_ == gollumfit::FluxComponent::atmConv)
                //   cache = e.cachedHoleIceConv;
                // else if (flux_component_ == gollumfit::FluxComponent::atmPrompt)
                //   cache = e.cachedHoleIcePrompt;
                // else if (flux_component_ == gollumfit::FluxComponent::diffuseAstro)
                //   cache = e.cachedHoleIceAstro;
                // else
                //   throw std::runtime_error("HoleIce correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");

                double rate, derivative;
                double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_.value()};
                int centers[3];

                auto holeicecorrection = hole_ice_map_.find(std::make_pair(flux_component_,static_cast<gollumfit::Topology>(e.topology)));
                if(holeicecorrection == hole_ice_map_.end())
                    throw std::runtime_error("HoleIce correction for " + gollumfit::GetFluxComponentName(flux_component_) + " flux component and for topology " + gollumfit::GetTopologyName(static_cast<gollumfit::Topology>(e.topology)) + " not found.");
                rate=(*(*holeicecorrection).second)(coordinates);

                if(not (*(*holeicecorrection).second).searchcenters(coordinates,centers))
                    throw std::runtime_error("Out of phase space: holeice spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(hole_ice_forward_.value()) + ").");
                derivative=(*(*holeicecorrection).second).ndsplineeval(coordinates,centers,1<<2);
                derivative*=hole_ice_forward_.derivative(didx);

                coordinates[2]=-1.;
                double cache = (*(*holeicecorrection).second)(coordinates);

                if(rate == 0 or cache == -std::numeric_limits<float>::max()){
                  // this is out of parameter space, ignore the event
                  result_type r(0.0);
                  r.setDerivative(didx,0.0);
                  return r;
                }

                result_type r(rate);
                r.setDerivative(didx,derivative);
                return(pow(10.,r-cache));
            }
    };

//used for initializing the per-event dom-eff related caches
/**
* @brief Struct to initialize per-event caches for DOM efficiency related to hole ice corrections.
* 
* The HoleIceSetter is used to precompute and store correction factors for individual events
* based on their energy and zenith angle, avoiding redundant calculations during the
* application of hole ice corrections.
* 
* @tparam Event The event type which contains properties such as energy and zenith angle.
*/
template<typename Event>
    struct HoleIceSetter{
        const HoleIceMapType& hole_ice_map_;///< Reference to the hole ice correction map.
        const double hole_ice_forward_;///< Hole ice model parameter.
        /**
        * @brief Constructor. 
        * 
        * @param hole_ice_map Reference to the map of spline tables for hole ice corrections.
        * @param hole_ice_forward A parameter related to the hole ice model.
        */
        HoleIceSetter(const HoleIceMapType& hole_ice_map_, double hole_ice_forward_):
            hole_ice_map_(hole_ice_map_),hole_ice_forward_(hole_ice_forward_){}

        /**
        * @brief Sets the cache for hole ice corrections on an event.
        * 
        * Computes and caches the log of the correction rate for an event based on its properties.
        * If the cache evaluates to zero, the event is marked to be ignored by setting the cache
        * to the lowest possible float value.
        * 
        * @param e Reference to the event whose cache will be set.
        * @exception std::runtime_error Thrown if an unimplemented hole ice component is encountered.
        */
        void setCache(Event& e) const{
            double coordinates[3]={log10(e.energy),cos(e.zenith),hole_ice_forward_};
            assert(not hole_ice_map_.empty());
            for(auto element : hole_ice_map_){
              if(element.first.second == static_cast<gollumfit::Topology>(e.topology)){
                double cache = (*element.second)(coordinates);
                assert(std::isfinite(cache));
                if(cache == 0){
                  // the cache if the log of the rate; this cases the event to be ignored.
                  cache = -std::numeric_limits<float>::max();
                }
                if(element.first.first == gollumfit::FluxComponent::atmConv)
                  e.cachedHoleIceConv = cache;
                else if(element.first.first == gollumfit::FluxComponent::atmPrompt)
                  e.cachedHoleIcePrompt = cache;
                else if(element.first.first == gollumfit::FluxComponent::diffuseAstro)
                  e.cachedHoleIceAstro= cache;
                else
                  throw std::runtime_error("Unimplemented HoleIce component");
              }
            }
        }
    };


//================================================================================
// Attenuation WEIGHTER
//================================================================================

using AttenuationMapType=std::map<std::pair<gollumfit::FluxComponent,LW::ParticleType>,std::shared_ptr<splinetable<>>>;

template<typename Event, typename DataType>
    struct attenuationWeighter : public phys_tools::GenericWeighter<attenuationWeighter<Event,DataType>>{
        private:
            const AttenuationMapType & attenuation_uncertainty_spline_;
            gollumfit::FluxComponent flux_component_;
            DataType scale_nu_;
            DataType scale_nubar_;
        public:
            attenuationWeighter(const AttenuationMapType &  attenuation_uncertainty_spline_, gollumfit::FluxComponent flux_component_, DataType scale_nu_, DataType scale_nubar_):
                attenuation_uncertainty_spline_(attenuation_uncertainty_spline_),
                flux_component_(flux_component_),
                scale_nu_(scale_nu_),
                scale_nubar_(scale_nubar_)
            {}

            using result_type=double;
            result_type operator()(const Event& e) const{
              auto attenuation_ratio_spline = attenuation_uncertainty_spline_.find(std::make_pair(flux_component_,e.primaryType));
              if(attenuation_ratio_spline == attenuation_uncertainty_spline_.end())
                 return result_type(1.0);
              if(cos(e.primaryZenith)>0.1)
                 return result_type(1.0);
              DataType scale;
              if(isNeutrino(e.primaryType))
                scale = scale_nu_;
              else
                scale = scale_nubar_;
              double coordinates[3]={log10(e.primaryEnergy),cos(e.primaryZenith), scale};
              double correction = (*(*attenuation_ratio_spline).second)(coordinates);
              if(correction<0.0)
                throw std::runtime_error("GollumFit::attenuation weighter making a neg weight.");
              return correction;
            }
    };

template<typename Event,unsigned int Dim>
    struct attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>> : public phys_tools::GenericWeighter<attenuationWeighter<Event,phys_tools::autodiff::FD<Dim>>>{
        private:
            const AttenuationMapType &  attenuation_uncertainty_spline_;
            gollumfit::FluxComponent flux_component_;
            phys_tools::autodiff::FD<Dim> scale_nu_;
            phys_tools::autodiff::FD<Dim> scale_nubar_;
            unsigned int didx_nu_,didx_nubar_;
        public:
            attenuationWeighter(const AttenuationMapType &  attenuation_uncertainty_spline_, gollumfit::FluxComponent flux_component_,
                phys_tools::autodiff::FD<Dim> scale_nu_,
                phys_tools::autodiff::FD<Dim> scale_nubar_):
                attenuation_uncertainty_spline_(attenuation_uncertainty_spline_),
                flux_component_(flux_component_),
                scale_nu_(scale_nu_),
                scale_nubar_(scale_nubar_)
            {
              for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                if(scale_nu_.derivative(i)!=0){
                    didx_nu_=i;
                    break;
                }
              }
              for(unsigned int i=0; i<Dim; i++){ //this will break if Dim is Dynamic
                if(scale_nubar_.derivative(i)!=0){
                    didx_nubar_=i;
                    break;
                }
              }
            }

            using result_type=phys_tools::autodiff::FD<Dim>;
            result_type operator()(const Event& e) const{
              auto attenuation_ratio_spline = attenuation_uncertainty_spline_.find(std::make_pair(flux_component_,e.primaryType));
              if(attenuation_ratio_spline == attenuation_uncertainty_spline_.end())
                 return result_type(1.0);
              if(cos(e.primaryZenith)>0.1)
                 return result_type(1.0);
              double scale;
              unsigned int didx;
              if(isNeutrino(e.primaryType)){
                scale = scale_nu_.value();
                didx = didx_nu_;
              } else {
                scale = scale_nubar_.value();
                didx = didx_nubar_;
              }
              double coordinates[3]={log10(e.primaryEnergy),cos(e.primaryZenith), scale};
              double correction = (*(*attenuation_ratio_spline).second)(coordinates);

              int centers[3];
              if(not (*(*attenuation_ratio_spline).second).searchcenters(coordinates,centers))
                  throw std::runtime_error("Out of phase space: attenuation spline (" + std::to_string(log10(e.energy)) + ", " + std::to_string(cos(e.zenith)) + ", " + std::to_string(scale) + ").");
              double derivative=(*(*attenuation_ratio_spline).second).ndsplineeval(coordinates,centers,1<<2);

              result_type c(correction);
              c.setDerivative(didx,derivative);
              return(c);
            }
    };




//================================================================================
// THE WEIGHT MAKER
//================================================================================

// This function construct a Weighter object to weight each event independently.
//
// "Time can't be measured in days the way money is measured in pesos and centavos,
// because all pesos are equal, while every day, perhaps every hour, is different."
// JLB

namespace sterile {

/**
 * @brief Struct to construct the weighters for the fitting.
 * 
 * Handles the creation of weighters for all the nuisance parameters. 
 */
struct WeighterMaker{
    private:
        static constexpr double astroPivotEnergy=1.0e5; ///< astro pivot energy in GeV
        // DOM efficienfy splines;
        // Type aliases for maps associating spline tables.
        using DOMMapType=std::map<std::pair<gollumfit::FluxComponent,gollumfit::Topology>,std::shared_ptr<splinetable<>>>;
        using HoleIceMapType=std::map<std::pair<gollumfit::FluxComponent,gollumfit::Topology>,std::shared_ptr<splinetable<>>>;
        using AttenuationMapType=std::map<std::pair<gollumfit::FluxComponent, LW::ParticleType>,std::shared_ptr<splinetable<>>>;
        // Shared pointers to the maps holding the spline tables.
        std::shared_ptr<const DOMMapType> domefficiencySplines_;
        std::shared_ptr<const HoleIceMapType> holeiceSplines_;
        std::shared_ptr<const AttenuationMapType> attenuationSplines_;
        // bools to indicate whether splines are loaded correctly
        bool domeff_splines_loaded_ = false;
        bool holeice_splines_loaded_ = false;
        bool attenuation_splines_loaded_ = false;
        gollumfit::SteeringParams const * steering;
    public:
        // default constructor // bad bad
        /**
        * @brief Default constructor that initializes all members to their default states.
        * 
        * Typically not used.
        */
        WeighterMaker():
          domefficiencySplines_(nullptr),
          holeiceSplines_(nullptr),
          attenuationSplines_(nullptr),
          domeff_splines_loaded_(false),
          holeice_splines_loaded_(false),
          attenuation_splines_loaded_(false),
          steering(nullptr)
        {}

        /**
        * @brief Constructor that initializes members with provided spline maps and steering parameters.
        * 
        * @param dom_splines Map of DOM efficiency splines.
        * @param holeice_splines Map of hole ice splines.
        * @param attenuation_splines Map of attenuation splines.
        * @param steeringParams Steering parameters for the analysis.
        * @param datapaths Path to get the data, MC, and systematics splines.
        */
        WeighterMaker(const DOMMapType& dom_splines,
                      const HoleIceMapType& holeice_splines, 
                      const AttenuationMapType& attenuation_splines,
                      const gollumfit::SteeringParams& steeringParams,
                      const gollumfit::DataPaths& datapaths):
            domefficiencySplines_(std::make_shared<const DOMMapType>(dom_splines)),
            holeiceSplines_(std::make_shared<const HoleIceMapType>(holeice_splines)),
            attenuationSplines_(std::make_shared<const AttenuationMapType>(attenuation_splines)),
            domeff_splines_loaded_(true),
            holeice_splines_loaded_(true),
            attenuation_splines_loaded_(true),
            steering(&steeringParams)
        {
        }
        // copy constructor // CA
        /**
        * @brief Copy constructor.
        * 
        * @param guy The `WeighterMaker` object to copy from.
        */
        WeighterMaker(const WeighterMaker& guy):
          domefficiencySplines_(guy.domefficiencySplines_),
          holeiceSplines_(guy.holeiceSplines_),
          attenuationSplines_(guy.attenuationSplines_),
          domeff_splines_loaded_(guy.domeff_splines_loaded_),
          holeice_splines_loaded_(guy.holeice_splines_loaded_),
          attenuation_splines_loaded_(guy.attenuation_splines_loaded_),
          steering(guy.steering)
        {}
        // assign guy // CA
        /**
        * @brief Copy assignment operator.
        * 
        * @param other The `WeighterMaker` object to assign from.
        * @return A reference to the modified `WeighterMaker` object.
        */
        WeighterMaker& operator=(WeighterMaker& other){
          if(&other==this)
            return(*this);

          domefficiencySplines_ = other.domefficiencySplines_;
          holeiceSplines_ = other.holeiceSplines_;
          attenuationSplines_ = other.attenuationSplines_;
          domeff_splines_loaded_ = other.domeff_splines_loaded_;
          holeice_splines_loaded_ = other.holeice_splines_loaded_;
          attenuation_splines_loaded_ = other.attenuation_splines_loaded_;
          steering = other.steering;

          return(*this);
        }
        // move guy // CA
        /**
        * @brief Move assignment operator.
        * 
        * @param other The `WeighterMaker` object to move from.
        * @return A reference to the modified `WeighterMaker` object.
        */
        WeighterMaker& operator=(WeighterMaker&& other){
          if(&other==this)
            return(*this);

          domefficiencySplines_ = std::move(other.domefficiencySplines_);
          holeiceSplines_ = std::move(other.holeiceSplines_);
          attenuationSplines_ = std::move(other.attenuationSplines_);
          domeff_splines_loaded_ = std::move(other.domeff_splines_loaded_);
          holeice_splines_loaded_ = std::move(other.holeice_splines_loaded_);
          attenuation_splines_loaded_ = std::move(other.attenuation_splines_loaded_);
          steering = std::move(other.steering);

          return(*this);
        }
    public:
        // Methods to set the respective spline maps and mark them as loaded.
        void SetDOMEfficiencySplines(const DOMMapType& dom_splines){ domefficiencySplines_ = std::make_shared<const DOMMapType>(dom_splines); domeff_splines_loaded_ = true; }
        void SetHoleIceSplines(const HoleIceMapType& holeice_splines){ holeiceSplines_ = std::make_shared<const HoleIceMapType>(holeice_splines); holeice_splines_loaded_ = true; }
        void SetAttenuationSplines(const AttenuationMapType& attenuation_splines){ attenuationSplines_ = std::make_shared<const AttenuationMapType>(attenuation_splines); attenuation_splines_loaded_ = true; }
        void SetSteering(const gollumfit::SteeringParams& steeringParams){steering = &steeringParams;}

        /**
        * @brief Templated functor that calculates weights for events based on various components.
        * 
        * This function uses a vector of parameters to calculate the weights for different flux components and applies corrections for various physical effects.
        * 
        * @param params A vector of nuisance parameters used in the weight calculations.
        * @return A function that takes an `Event` and returns a weight of type `DataType`.
        */
        template<typename DataType>
            std::function<DataType(const Event&)> operator()(const std::vector<DataType>& params) const{

                //unpack things so we have legible names
                DataType convNorm                  = params[0];
                DataType promptNorm                = params[1];
                DataType adu                       = params[2];
                DataType klu                       = params[3];
                DataType hekp                      = params[4];
                DataType hekm                      = params[5];
                DataType vhe1pip                   = params[6];
                DataType vhe1pim                   = params[7];
                DataType vhe3kp                    = params[8];
                DataType vhe3km                    = params[9];
                DataType vhe3pip                   = params[10];
                DataType vhe3pim                   = params[11];
                DataType vhe3p                     = params[12];
                DataType vhe3n                     = params[13];
                DataType cr1                       = params[14];
                DataType cr2                       = params[15];
                DataType cr3                       = params[16];
                DataType cr4                       = params[17];
                DataType cr5                       = params[18];
                DataType cr6                       = params[19];
                DataType icegrad0                  = params[20];
                DataType icegrad1                  = params[21];
                DataType icegrad2                  = params[22];
                DataType icegrad3                  = params[23];
                DataType icegrad4                  = params[24];
                DataType icegrad5                  = params[25];
                DataType icegrad6                  = params[26];
                DataType icegrad7                  = params[27];
                DataType icegrad8                  = params[28];
                DataType deltaDomEff               = params[29];
                DataType holeiceForward            = params[30];
                DataType astroNorm                 = params[31];
                DataType astroDeltaGamma           = params[32];
                DataType astroDeltaGammaSec        = params[33];
                DataType astroPivot                = params[34];
                DataType NeutrinoAntineutrinoRatio = params[35];
                DataType nuxs                      = params[36];
                DataType nubarxs                   = params[37];


                // phys_tools::autodiff::FD<38> paux;
                // for (unsigned int i=0; i<params.size(); i++) { paux = params[i]; std::cout << paux.value() << "  "; }
                // std::cout << std::endl;

                using cachedWeighter=cachedValueWeighter<DataType,Event,double>;
                cachedWeighter astroFlux(&Event::cachedAstroWeight);
                cachedWeighter promptFlux(&Event::cachedPromptWeight);
                //flux
                cachedWeighter adu_wgt(&Event::cachedAtmDensity);
                cachedWeighter klu_wgt(&Event::cachedKaonLosses);
                //ice gradients
                cachedWeighter icegrad0_wgt(&Event::cachedIceGrad0);
                cachedWeighter icegrad1_wgt(&Event::cachedIceGrad1);
                cachedWeighter icegrad2_wgt(&Event::cachedIceGrad2);
                cachedWeighter icegrad3_wgt(&Event::cachedIceGrad3);
                cachedWeighter icegrad4_wgt(&Event::cachedIceGrad4);
                cachedWeighter icegrad5_wgt(&Event::cachedIceGrad5);
                cachedWeighter icegrad6_wgt(&Event::cachedIceGrad6);
                cachedWeighter icegrad7_wgt(&Event::cachedIceGrad7);
                cachedWeighter icegrad8_wgt(&Event::cachedIceGrad8);

                using neuaneu_t = antiparticleWeighter<Event,DataType>;
                neuaneu_t neuaneu_wgt(NeutrinoAntineutrinoRatio);

                using domEffW_t = DOMEffWeighter<Event,DataType>;
                domEffW_t convDOMEff(*domefficiencySplines_,deltaDomEff,gollumfit::FluxComponent::atmConv);
                domEffW_t promptDOMEff(*domefficiencySplines_,deltaDomEff,gollumfit::FluxComponent::atmPrompt);
                domEffW_t astroDOMEff(*domefficiencySplines_,deltaDomEff,gollumfit::FluxComponent::diffuseAstro);

                using holeIceW_t = holeiceWeighter<Event,DataType>;
                holeIceW_t convHoleIce(*holeiceSplines_,holeiceForward,gollumfit::FluxComponent::atmConv);
                holeIceW_t promptHoleIce(*holeiceSplines_,holeiceForward,gollumfit::FluxComponent::atmPrompt);
                holeIceW_t astroHoleIce(*holeiceSplines_,holeiceForward,gollumfit::FluxComponent::diffuseAstro);

                using attW_t = attenuationWeighter<Event,DataType>;
                attW_t  convAtt((*attenuationSplines_), (gollumfit::FluxComponent::atmConv), nuxs, nubarxs);
                attW_t  promptAtt((*attenuationSplines_), (gollumfit::FluxComponent::atmPrompt), nuxs, nubarxs);
                attW_t  astroAtt((*attenuationSplines_), (gollumfit::FluxComponent::diffuseAstro), nuxs, nubarxs);

                auto icegrad_wgt = (1.+icegrad0*icegrad0_wgt)*
                                   (1.+icegrad1*icegrad1_wgt)*
                                   (1.+icegrad2*icegrad2_wgt)*
                                   (1.+icegrad3*icegrad3_wgt)*
                                   (1.+icegrad4*icegrad4_wgt)*
                                   (1.+icegrad5*icegrad5_wgt)*
                                   (1.+icegrad6*icegrad6_wgt)*
                                   (1.+icegrad7*icegrad7_wgt)*
                                   (1.+icegrad8*icegrad8_wgt);


                auto atm_wgt = (1.+adu*adu_wgt)*
                               (1.+klu*klu_wgt);

                auto conv   = convHoleIce   * convDOMEff   * convAtt   * atm_wgt * ConvFluxWeigther<Event,DataType>(hekp,hekm,vhe1pip,vhe1pim,vhe3kp,vhe3km,vhe3pip,vhe3pim,vhe3p,vhe3n,cr1,cr2,cr3,cr4,cr5,cr6);
                auto prompt = promptNorm * promptHoleIce * promptDOMEff * promptAtt * promptFlux;
                auto astro  = astroNorm  * astroHoleIce  * astroDOMEff  * astroAtt  * astroFlux * neuaneu_wgt * brokenpowerlawTiltWeighter<Event,DataType>(astroPivot, astroDeltaGamma, astroDeltaGammaSec);
                
                if (steering->enableTotalNorm) { return convNorm*(conv+prompt+astro)*icegrad_wgt; }
                else { return (convNorm*conv+prompt+astro)*icegrad_wgt; }
            }
};

} // close sterile namespace

#endif /* ANALYSISWEIGHTING_H_ */
