#ifndef GOLLUMFIT_H_
#define GOLLUMFIT_H_

#include <sys/types.h>
#include <sys/stat.h>
#include <cstdio>
#include <stdlib.h>
#include <iterator>
#include <algorithm>
#include <functional>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>
#include <memory>

#include <LeptonWeighter/Weighter.h>
#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>
#include <PhysTools/likelihood/likelihood.h>
#include <PhysTools/histogram.h>
#include <PhysTools/bin_types.h>

#include "GollumParameters.h"
#include "Event.h"
#include "analysisWeighting.h"
#include "compactIO.h"
#include "GollumEnumDefinitions.h"
#include "GollumMCSet.h"
#include "GollumTools.h"

namespace gollumfit{

struct simpleLocalDataWeighter{
  double operator()(const Event& e) const{
    return(e.cachedWeight);
  }
};

struct simpleLocalDataWeighterConstructor{
    template<typename DataType>
    simpleLocalDataWeighter operator()(const std::vector<DataType> params) const{
        return simpleLocalDataWeighter();
    }
};


using namespace nusquids;
using namespace phys_tools::histograms;
using namespace phys_tools::likelihood;
using HistType = std::tuple< histogram<3,entryStoringBin<std::reference_wrapper<const Event>>> >;
using PriorIndices = std::tuple<phys_tools::likelihood::parameters<>,phys_tools::likelihood::parameters<4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19>,phys_tools::likelihood::parameters<20,21,22,23,24,25,26,27,28>>;
using BasicPrior = FixedSizePriorSet<GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     UniformPrior,
                                     GaussianPrior,
                                     GaussianPrior,
                                     GaussianPrior>;

using CPrior = ArbitraryPriorType<PriorIndices, BasicPrior, GaussianNDPrior<16>, GaussianNDPrior<9>>::type;
using LType=LikelihoodProblem<std::reference_wrapper<const Event>, HistType,simpleLocalDataWeighterConstructor,sterile::WeighterMaker,CPrior,SAYLikelihoodRelativeUncertaintyMod,38>;
using hist_marray=marray<double,3>;

template<typename ContainerType, typename HistType, typename BinnerType>
  void bin(const ContainerType& data, HistType& hist, const BinnerType& binner){
  for(const Event& event : data)
    binner(hist,event);
}

/**
* @class GollumFit
* @brief Main GollumFit object.
* 
* Key components include the loading of data and MC, appropriately setting event weights and re-weighting them,
* loading the required splines and resources needed for re-weighting, obtaining the expectation histogram given 
* nusiance parameters, calculating the likelihood, and performing likelihood minimization.
*/
class GollumFit {
  private:
    // All the local configuration variables
    SteeringParams  steeringParams_; ///
    DataPaths       dataPaths_;

    // to store events
    std::deque<Event> metaEvents_;
    std::deque<Event> mainSimulation_;
    std::deque<Event> sample_;

    // for fast mode only
    std::deque<Event> auxSimulation_;

    // histograms
    HistType dataHist_; // analysis data
    HistType simHist_; // analysis MC is

    // minimizing objects
    std::vector<FitParameters> fitSeed_;
    FitParameters asimovSeed_;
    FitParametersFlag fixedParams_;
    FitParametersBound boundParams_;
    Priors priors_;

    // weighter object
    sterile::WeighterMaker DFWM;
    std::shared_ptr<LW::Flux> fluxConv_,fluxPrompt_,fluxAstro_;
    std::shared_ptr<LW::CrossSectionFromSpline> xsw_;
    std::vector<std::shared_ptr<LW::Generator>> mcw_;
    std::shared_ptr<LW::Weighter> convFluxWeighter_;
    std::shared_ptr<LW::Weighter> promptFluxWeighter_;
    std::shared_ptr<LW::Weighter> astroFluxWeighter_;

    // Status flags
    bool xs_weighter_constructed_ = (false);
    bool flux_weighter_constructed_ = (false);
    bool lepton_weighter_constructed_ = (false);
    bool data_histogram_constructed_ = (false);
    bool simulation_histogram_constructed_ = (false);
    bool mc_generation_weighter_constructed_ = (false);
    bool simulation_loaded_ = (false);
    bool data_loaded_ = (false);
    bool likelihood_problem_constructed_ = (false);
    bool simulation_initialized_ = (false);
    bool fastmode_constructed_= (false);
    bool fitSeed_constructed_= (false);
    bool priors_constructed_= (false);
    bool fixedParams_constructed_= (false);
    bool boundParams_constructed_= (false);

    // splines flags
    bool domeff_spline_loaded_ = (false);
    bool holeice_spline_loaded_ = (false);
    bool attenuation_spline_loaded_ = (false);
    bool ice_gradient_spline_loaded_ = (false);
    bool atmospheric_density_spline_loaded_ = (false);
    bool atmospheric_kaonlosses_spline_loaded_ = (false);
    bool hadronic_spline_loaded_ = (false);
    bool cosmic_ray_spline_loaded_ = (false);



    // Ice Gradient splines
    std::map<std::pair<int,Topology>,std::shared_ptr<splinetable<>>> iceGradientSplines_;

    // DOM efficiency splines
    std::map<std::pair<FluxComponent,Topology>, std::shared_ptr<splinetable<>>> domefficiencySplines_;
    DOMEfficiencySetter<Event> * domEffSetter_;
    double domEffSpline_minValidValue_ = 0.5;
    double domEffSpline_maxValidValue_ = 1.5;

    // Hole ice splines
    std::map<std::pair<FluxComponent,Topology>, std::shared_ptr<splinetable<>>> holeIceSplines_;
    HoleIceSetter<Event> * holeIceSetter_;
    double holeIceSpline_minValidValue_ = -1.5;
    double holeIceSpline_maxValidValue_ =  0.5;

    // Attenuation splines
    std::map<std::pair<FluxComponent,LW::ParticleType>, std::shared_ptr<splinetable<>>> attenuationSplines_;
    double attenuationSpline_minValidValue_ = 0.5;
    double attenuationSpline_maxValidValue_ = 1.5;

    // atmopsheric density spline
    std::shared_ptr<splinetable<>> atmosphericDensityUncertaintySpline_ = nullptr;

    // kaon losses spline
    std::shared_ptr<splinetable<>> kaonLossesUncertaintySpline_ = nullptr;

    // hadronic splines
    std::map<HadronicParameter,std::shared_ptr<LW::nuSQUIDSAtmFlux<>>> hadronicFluxWeighter_;

    // cosmic ray splines
    std::map<CosmicRayParameter,std::shared_ptr<LW::nuSQUIDSAtmFlux<>>> cosmicRayFluxWeighter_;

    // likehood problem object
    std::shared_ptr<LType> prob_;

    // Nasty template part of fit function

    /**
    * @brief Performs a fit using the LBFGSB algorithm.
    * 
    * This function template takes a likelihood function and minimizer object and
    * uses the Limited-memory Broyden-Fletcher-Goldfarb-Shanno with Bounds (L-BFGS-B)
    * algorithm to minimize the likelihood function. The likelihood function is wrapped
    * in a @c BFGS_Function object to facilitate the minimization process.
    * 
    * @tparam @c LikelihoodType The type of the likelihood function object. This type must
    *         have an overloaded operator() that computes the likelihood given a set
    *       of parameters.
    * @param likelihood A reference to the likelihood function object.
    * @param minimizer A reference to the an instance of
    *                   @c phys_tools::lbfgsb::LBFGSB_Driver which will perform the minimization.
    *                   note that attributes of this object must be set, such as the list of variables to minimize over.
    * @return True if the minimization succeeded, false otherwise.
    */
    template<typename LikelihoodType>
      bool DoFitLBFGSB(LikelihoodType& likelihood, phys_tools::lbfgsb::LBFGSB_Driver& minimizer) const{
      using namespace phys_tools::likelihood;
      bool succeeded=minimizer.minimize(BFGS_Function<LikelihoodType>(likelihood));
      return succeeded;
    }

  protected:
    // to check events
    bool IsEventHealthy(Event& e, bool print_reason = true) const;

    /**
    * @brief Adjusts the fit seed values to ensure they are within the valid ranges of the spline functions.
    *
    * This function modifies the hole ice forward and DOM efficiency parameters of the fit seed
    * if they fall outside their respective valid ranges. It sets them to the midpoint of the valid range.
    *
    * @param fitSeed The fit seed parameters to be adjusted if necessary.
    */
    void ForceFitSeedSanity(FitParameters& fitSeed);
    
    /**
    * @brief Ensures that all fit seeds are within valid ranges.
    *
    * Iterates over all fit seeds and enforces them to fall within valid ranges by
    * calling ForceFitSeedSanity on each individual @c FitParameters object.
    */
    void ForceFitSeedSanity();
    
    /**
    * @brief Loads experimental data into the framework from a specified file.
    *
    * Events are loaded from a file and each event is processed by
    * assigning a cached weight (of 1) and topology before being added to the sample set.
    * If loading fails, an error message is displayed and an exception is thrown
    * if no data is loaded.
    * 
    * @throws std::runtime_error If the sample is not found after attempting to load data.
    */
    void LoadData();

    /**
    * @brief Loads MC simulation data into the framework.
    *
    * Loads full MC simulation data. It processes
    * each event, setting cached weights and properties that are required
    * for the fitting. To set some cached weights, it will read some necessary spline files.
    * If any issues occur during event processing, the
    * event is skipped and an error is displayed.
    * 
    * 
    * @throws std::runtime_error If the simulation data is not found after loading attempts.
    * @throws std::runtime_error If an impossible hadronic or cosmic ray parameter label is encountered.
    */
    void LoadMC();

    /**
    * @brief Loads FastMC from a compact file.
    * 
    * The function loads the FastMC from a compact format file with a '.meows'
    * extension into the class members.
    */
    void LoadCompact();

    /**
    * @brief Clears the data container.
    * 
    * This function is used to clear the container @c sample_ that holds the sample data within the class.
    */
    void ClearData();

    /**
    * @brief Clears the simulation data containers.
    * 
    * This function is used to clear the containers @c mainSimulation_ and @c metaEvents_
    * that hold the main simulation data and meta events within the class.
    */
    void ClearSimulation();

    // loading DOM efficiency splines
    /**
    * @brief Loads the DOM Efficiency splines for atmospheric and astrophysical neutrino flux components.
    * 
    * Initializes the DOM efficiency splines for atmospheric conventional, atmospheric prompt, and diffuse
    * astrophysical neutrino flux components based on the topology. The splines are loaded from FITS files
    * and inserted into the @c domefficiencySplines_ map. The function also determines and sets the
    * validity range of the loaded splines and updates the DFWM weighter object with these splines.
    * 
    * @pre The @c dataPaths_.domeff_spline_path must be set correctly to the directory containing the
    * spline FITS files.
    * 
    * @post The @c domefficiencySplines_ map is populated with the loaded splines. The validity range
    * for the efficiency splines is set, and the DFWM is updated.
    * 
    * @return true if the splines are loaded successfully, false otherwise.
    * 
    * @throws std::runtime_error If no splines are loaded (assertion fails).
    */
    bool LoadDOMEfficiencySplines();

    /**
    * @brief Loads the Hole Ice splines for atmospheric and astrophysical neutrino flux components.
    * 
    * Initializes the Hole Ice splines for atmospheric conventional, atmospheric prompt, and diffuse
    * astrophysical neutrino flux components based on the topology. The splines are loaded from FITS files
    * and inserted into the @c holeIceSplines_ map. The function also determines and sets the
    * validity range of the loaded splines and updates the DFWM weighter object with these splines.
    * 
    * @pre The @c dataPaths_.holeice_spline_path must be set correctly to the directory containing the
    * spline FITS files.
    * 
    * @post The @c holeIceSplines_ map is populated with the loaded splines. The validity range
    * for the Hole Ice splines is set, and the DFWM is updated.
    * 
    * @return true if the splines are loaded successfully, false otherwise.
    * 
    * @throws std::runtime_error If no splines are loaded (assertion fails).
    */
    bool LoadHoleIceResources();

    /**
    * @brief Loads the attenuation splines for different flux components and particle types.
    * 
    * This function initializes the attenuation splines for atmospheric conventional, atmospheric prompt,
    * and diffuse astrophysical flux components for neutrino and antineutrino types.
    * The splines are loaded from FITS files and stored in the @c attenuationSplines_ map. The function
    * also asserts that the correct number of splines have been loaded and sets the validity range
    * for the splines based on the first spline loaded. These splines are used in the framework to
    * account for the attenuation of neutrinos as they pass through the Earth due to cross section uncertainties.
    * 
    * @pre The @c dataPaths_.attenuation_spline_path must be set correctly to the directory containing
    * the attenuation spline FITS files.
    * 
    * @post The @c attenuationSplines_ map is populated with the loaded splines. An assertion checks
    * that 12 splines (3 flux components × 4 particle types) are loaded. The DFWM is updated with
    * the loaded splines.
    * 
    * @return true if the splines are loaded successfully, false otherwise.
    * 
    * @throws std::runtime_error If the assertion on the number of loaded splines fails.
    */
    bool LoadAttenuationResources();

    /**
    * @brief Loads the ice gradient splines used for weighting the effect of ice anisotropy.
    * 
    * This function reads from FITS files to load splines representing ice gradients, which are
    * used to simulate the effect of anisotropic ice on neutrino detection. The splines are indexed
    * by an integer and associated with a particular topology. The loaded splines
    * are inserted into the @c iceGradientSplines_ map. The selection of which splines to load is
    * dependent on the @c steeringParams_.selectionStart parameter.
    * 
    * In the latest edition of MEOWS, these splines are no longer used in favor of gradients. 
    * 
    * @pre The @c dataPaths_.ice_gradient_spline_path and @c steeringParams_.ice_gradient_filename
    * must be correctly set to point to the directory and names of the ice gradient spline FITS files.
    * 
    * @post The @c iceGradientSplines_ map is populated with the loaded splines.
    * 
    * @return true if the splines are loaded successfully, false otherwise.
    */
    bool LoadIceGradientResources();

    /**
    * @brief Loads the spline representing uncertainties in Kaon energy losses.
    *
    * This function loads a spline from a specified file path which models the uncertainties
    * in the energy losses of Kaons in the atmosphere.
    *
    * @pre @c dataPaths_.atmospheric_kaonlosses_spline_path must be set to the file path
    * of the Kaon losses uncertainty spline.
    *
    * @return true if the spline is loaded successfully, false otherwise.
    */
    bool LoadKaonEnergyLossesUncertaintiesResources();

    /**
    * @brief Loads the spline representing uncertainties in atmospheric density.
    *
    * This function loads a spline from a specified file path which models the uncertainties
    * in atmospheric density.
    *
    * @pre @c dataPaths_.atmospheric_density_spline_path must be set to the file path
    * of the atmospheric density uncertainty spline.
    *
    * @return true if the spline is loaded successfully, false otherwise.
    */
    bool LoadAtmosphericDensityUncertaintiesResources();

    /**
    * @brief Loads resources related to hadronic interaction uncertainties.
    *
    * This function initializes @c hadronicFluxWeighter_ with resources that model
    * the uncertainties in hadronic interactions. It iterates over a set of hadronic
    * interaction parameters, loading the corresponding resources.
    *
    * @pre @c dataPaths_.hadronic_spline_path must be set to the directory containing
    * the hadronic uncertainty resources, and @c steeringParams_.model_label must be set
    * to the appropriate model label.
    *
    * @post @c hadronicFluxWeighter_ is populated with the loaded resources.
    *
    * @return true if all resources are loaded successfully, false otherwise.
    */
    bool LoadHadronicUncertaintiesResources();

    /**
    * @brief Loads resources related to cosmic ray interaction uncertainties.
    *
    * This function initializes @c cosmicRayFluxWeighter_ with resources that model
    * the uncertainties in cosmic ray interactions. It iterates over a set of cosmic
    * ray parameters, loading the corresponding resources into the framework.
    *
    * @pre @c dataPaths_.cosmic_ray_spline_path must be set to the directory containing
    * the cosmic ray uncertainty resources, and @c steeringParams_.model_label must be set
    * to the appropriate model label.
    *
    * @post @c cosmicRayFluxWeighter_ is populated with the loaded resources.
    *
    * @return true if all resources are loaded successfully, false otherwise.
    */
    bool LoadCosmicRayUncertaintiesResources();
    
    // functions to construct the weighters

    /**
    * @brief Constructs the neutrino cross section weighter from spline files.
    *
    * Initializes the cross section weighter with neutrino and antineutrino cross section splines for
    * both charged current (CC) and neutral current (NC) interactions. The splines are loaded from the
    * specified file paths.
    *
    * @pre The file paths for cross section splines must be set correctly, such as @c dataPaths_.diff_neutrino_cc_xs_spline_path and so on.
    * @post xsw_ is populated with the cross section weighter and @c xs_weighter_constructed_ is set to true.
    */
    void ConstructCrossSectionWeighter();

    /**
    * @brief Constructs the atmospheric and astrophysical neutrino flux weighters.
    *
    * Initializes the flux weighters for conventional atmospheric, prompt atmospheric, and astrophysical
    * neutrinos. The corresponding nuSQuIDS flux objects are created from the specified file paths.
    *
    * @pre The file paths for nuSQuIDS atmospheric and astrophysical flux files must be set correctly 
    * as @c dataPaths_.conventional_nusquids_atmospheric_file and so on.
    * @post @c fluxConv_ , @c fluxPrompt_ , @c fluxAstro_ are populated with the respective flux weighters,
    *       @c flux_weighter_constructed_ is set to true, and @c simulation_initialized_ is set to false.
    */
    void ConstructFluxWeighter();

    /**
    * @brief Constructs a Monte Carlo generation weighter based on the simulation sets.
    *
    * This function creates a weighter for each generator in the Monte Carlo simulation sets,
    * which are determined based on the simulation parameters and paths provided.
    *
    * @pre @c steeringParams_.simToLoad and @c dataPaths_.mc_path must be set to identify the correct
    * simulation sets.
    * @post @c mcw_ is filled with the Monte Carlo generation weighters and @c mc_generation_weighter_constructed_
    * is set to true.
    * @throws @c std::runtime_error If the @c mcw_ is empty after attempting to fill it.
    */
    void ConstructMonteCarloGenerationWeighter();

    /**
    * @brief Constructs the lepton weighters for different neutrino flux types.
    *
    * Creates weighters for conventional, prompt, and astrophysical neutrinos based on the
    * previously constructed flux and cross section weighters, as well as the Monte Carlo
    * generation weighter.
    *
    * @pre @c mc_generation_weighter_constructed_ , @c flux_weighter_constructed_, and
    * @c xs_weighter_constructed_ must all be true. If any of these weighters is not constructed,
    * an exception will be thrown.
    * @post @c convFluxWeighter_ , @c promptFluxWeighter_ , @c astroFluxWeighter_ are populated with
    * their respective weighters and @c lepton_weighter_constructed_ is set to true.
    * @throws @c std::runtime_error If any of the required weighters are not constructed prior to
    * calling this function.
    */
    void ConstructLeptonWeighter();
    // functions to construct the histograms of data and simulation

    /**
    * @brief Constructs a histogram for the data.
    *
    * Initializes the data histogram using the steering parameters for energy and zenith dimensions.
    * The histogram is then filled with the loaded data. The energy range and zenith
    * range for the histogram axes are set according to the steering parameters.
    *
    * @pre @c data_loaded_ must be true, indicating that data has been loaded into the framework.
    * @post @c dataHist_ contains the histogram of the loaded data, with axes limits set, and
    *       @c data_histogram_constructed_ is set to true.
    * @throws @c std::runtime_error If no data has been loaded prior to calling this function.
    */
    void ConstructDataHistogram();

    /**
    * @brief Constructs a histogram for the MC simulation.
    *
    * Initializes the simulation histogram to match the structure of the data histogram and then
    * fills it with the loaded MC simulation data.
    *
    * @pre @c simulation_loaded_ must be true, indicating that simulation has been loaded into the framework.
    * @pre @c data_histogram_constructed_ must be true, indicating that the data histogram has been constructed and can be used as a template.
    * @post @c simHist_ contains the histogram of the loaded simulation data and
    *       @c simulation_histogram_constructed_ is set to true.
    * @throws @c std::runtime_error If no simulation has been loaded or if the data histogram has not been
    *       constructed prior to calling this function.
    */
    void ConstructSimulationHistogram();

  public:

    /**
    * @brief GollumFit object constructor.
    *
    * This constructor initializes a GollumFit instance by setting up the necessary
    * configuration and parameters for the fit. It encapsulates the configuration
    * within shared pointers for dynamic memory management and potential shared ownership.
    *
    * @param dataPaths Struct containing filepaths to be used in the fit.
    * @param steeringParams Struct containing steering parameters.
    */
    GollumFit(DataPaths dataPaths, SteeringParams steeringParams);

    /**
    * @brief Reconfigures the analysis framework with new data paths and steering parameters. The function template
    * called in the GollumFit constructor.
    * 
    * This function template is used to update the configuration of the analysis framework. It conditionally resets
    * and reloads various components such as splines, weighters, and histograms based on the new configuration provided.
    * If called as a constructor (template parameter constructor set to true), it will always reset; otherwise, 
    * it will reset only if the new configuration objects are different from the current ones.
    * 
    * @tparam constructor A boolean template parameter that determines whether the function should act as a constructor,
    *                     which forces a reset of configuration regardless of the new objects being different.
    * @param dataPaths Shared pointer to a dataPaths object containing new paths to data files.
    * @param steeringParams Shared pointer to a steeringParams object containing new steering parameters.
    * 
    * @pre The function should be called with valid shared pointers to dataPaths and steeringParams.
    *      If constructor is false, the function compares the addresses of the provided shared pointers to the
    *    current configuration to decide if a reset is necessary.
    * 
    * @post The configuration is updated with the new data paths and steering parameters. Components such as splines,
    *       weighters, and histograms are reloaded if necessary. The internal state flags like @c domeff_spline_loaded_ , etc.
    *       are updated to reflect the new configuration status.
    * 
    * @throws @c std::runtime_error If any of the conditions for reloading data or simulation are not met after the
    *         update, or if the new data paths are invalid.
    */
    template<bool constructor>
    void ReConfig(std::shared_ptr<DataPaths> dataPaths, std::shared_ptr<SteeringParams> steeringParams) {
      bool reset_steering = constructor || (bool(steeringParams) && (&steeringParams_ != steeringParams.get()));
      bool reset_data = constructor || (bool(dataPaths) && (&dataPaths_ != dataPaths.get()));

      // default minimizer seed
      fitSeed_ = {FitParameters()};

      std::cout << "reset_steering: " << reset_steering << std::endl;
      std::cout << "reset_data: " << reset_data << std::endl;

      std::shared_ptr<DataPaths> old_data = std::shared_ptr<DataPaths>(reset_data ? new DataPaths(dataPaths_) : nullptr);
      if(reset_data) {
        dataPaths_ = *dataPaths;
      }

      std::shared_ptr<SteeringParams> old_steering = std::shared_ptr<SteeringParams>(reset_steering ? new SteeringParams(steeringParams_) : nullptr);
      if(reset_steering) {
        steeringParams_ = *steeringParams;
      }

      if(reset_steering || reset_data) {

        if(dataPaths_.domeff_spline_path!="") domeff_spline_loaded_ = LoadDOMEfficiencySplines();
        if(dataPaths_.holeice_spline_path!="") holeice_spline_loaded_ = LoadHoleIceResources();
        if(dataPaths_.attenuation_spline_path!="") attenuation_spline_loaded_ = LoadAttenuationResources();
        if(dataPaths_.ice_gradient_spline_path!="") ice_gradient_spline_loaded_ = LoadIceGradientResources();
        if(dataPaths_.atmospheric_density_spline_path!="") atmospheric_density_spline_loaded_ = LoadAtmosphericDensityUncertaintiesResources();
        if(dataPaths_.atmospheric_kaonlosses_spline_path!="") atmospheric_kaonlosses_spline_loaded_ = LoadKaonEnergyLossesUncertaintiesResources();
        if(dataPaths_.hadronic_spline_path!="") hadronic_spline_loaded_ = LoadHadronicUncertaintiesResources();
        if(dataPaths_.cosmic_ray_spline_path!="") cosmic_ray_spline_loaded_ = LoadCosmicRayUncertaintiesResources();

        std::cout<<"Constructing DiffuseWeightMaker" <<std::endl;
        DFWM = sterile::WeighterMaker(domefficiencySplines_,
                                      holeIceSplines_,
                                      attenuationSplines_,
                                      steeringParams_,
                                      dataPaths_);

        if(dataPaths_.compact_file_path!=""){
            std::cout<<"Loading compact data" <<std::endl;
            LoadCompact();
        } else {
          std::cout<<"Loading Flux weighter" <<std::endl;
          ConstructFluxWeighter();
          std::cout<<"Loading XS" <<std::endl;
          ConstructCrossSectionWeighter();
          std::cout<<"Loading MC weighter" <<std::endl;
          ConstructMonteCarloGenerationWeighter();
          std::cout<<"Loading Lepton weighter" <<std::endl;
          ConstructLeptonWeighter();
          // std::cout<<"Loading data" <<std::endl;
          // assert(CheckDataPath(dataPaths_.data_path));
          // LoadData();
          std::cout<<"Loading MC" <<std::endl;
          assert(CheckDataPath(dataPaths_.mc_path));
          LoadMC();
        }

        // std::cout<<"Making data hist" <<std::endl;
        // ConstructDataHistogram();
        std::cout<<"Making sim hist" <<std::endl;
        ConstructSimulationHistogram();

      }
    }


    /**
    * @brief Checks if the given directory path exists and is accessible.
    *
    * This function verifies the existence and accessibility of a directory at the
    * specified path. If the path is not empty, it checks whether the directory
    * can be accessed and whether it is indeed a directory and not a file. If the
    * checks fail, an exception is thrown with a relevant error message.
    *
    * @param p The path to the directory to check.
    * @return True if the directory exists and is accessible, false otherwise.
    * @throws @c std::runtime_error If the path is not accessible or not a directory.
    */
    bool CheckDataPath(std::string p) const;

    /**
    * @brief Ensures that the specified file path points to a valid file.
    *
    * This function checks if the file at the given path exists and can be
    * opened. If the file is valid, the path is returned. If the file does not
    * exist or cannot be opened, an exception is thrown.
    *
    * @param FilePath The path to the file to check.
    * @return The path to the file if it exists and is accessible.
    * @throws @c std::runtime_error If the file does not exist or cannot be accessed.
    */
    std::string CheckedFilePath(std::string) const;

    /**
    * @brief Writes the FastMC result to a file.
    * 
    * The function writes the FastMC to a compact format file
    * with a '.meows' extension. 
    * 
    * @param file_path The base path where the compact file will be created.
    * @return Returns true on successful write operation.
    */
    bool WriteCompact(std::string) const;


    /**
    * @brief Constructs a "fast mode" or "FastMC" simulation using a provided scaling parameter.
    * 
    * This function is responsible for initializing a fast mode simulation, which is a
    * reduced representation of the full simulation dataset. The fast mode
    * allows for quicker computations at the cost of some accuracy.
    *
    * The function performs a series of checks to ensure that both the simulation and
    * data have been loaded before proceeding. It then constructs a meta histogram
    * with modified bin widths based on the meta scaling parameter. Subsequently, it
    * processes the simulation data to populate this meta histogram, effectively
    * reducing the dataset size.
    *
    * The function prints a warning message if the meta scaling parameter is set to
    * a value that could lead to inaccurate results. It throws an exception if the
    * meta scaling is greater than one, as this is likely to produce unreliable
    * outcomes.
    *
    * @param meta_scaling A double value that scales the bin width of the histograms.
    *                   Values greater than 1 are not allowed and will throw an exception.
    * 
    * @throws @c std::runtime_error If no simulation has been loaded.
    * @throws @c std::runtime_error If no data has been loaded.
    * @throws @c std::runtime_error If meta_scaling is greater than one.
    * 
    * @note The function outputs to the standard console the before and after sizes of
    *       the Monte Carlo (MC) event sets, indicating the reduction in dataset size.
    *       It also asserts that the resulting set of meta events is not empty.
    * 
    * Example python usage of creating a FastMC dataset: 
    * @code 
    * gollumfit = gf.GollumFit(datapaths,steering_params)
    * gollumfit.ConstructFastMode(float(meta_scaling))
    * gollumfit.WriteCompact(outpath)
    * @endcode
    */
    void ConstructFastMode(double meta_scaling);

    /**
    * @brief Sets up the likelihood problem.
    *
    * This function prepares the likelihood problem by checking that all necessary components
    * have been constructed. These include the data histogram, simulation histogram, fit seed,
    * priors, and fixed parameters. Once verified, it initializes the priors for the fit parameters
    * based on the provided prior settings. It then combines these priors with other settings
    * to construct the full likelihood problem, by calling @c phys_tools::likelihood::makeLikelihoodProblem 
    * which is encapsulated in a shared pointer to a likelihood problem object. 
    * This object is later used for likelihood evaluations.
    *
    * @pre The following conditions must be satisfied before constructing the likelihood problem:
    *    - The data histogram must be constructed ( @c data_histogram_constructed_ should be true).
    *    - The simulation histogram must be constructed ( @c simulation_histogram_constructed_ should be true).
    *    - The fitting seed must be constructed ( @c fitSeed_constructed_ should be true).
    *    - The priors must be constructed ( @c priors_constructed_ should be true).
    *    - The fixed parameters must be constructed ( @c fixedParams_constructed_ should be true).
    * 
    * @post The likelihood problem is constructed and encapsulated within @c prob_ , and
    *       the flag @c likelihood_problem_constructed_ is set to true, indicating that the
    *       likelihood problem is ready to be used for fitting.
    *
    * @throws @c std::runtime_error If any of the components required for constructing the likelihood
    *         problem are not yet constructed.
    */
    void ConstructLikelihoodProblem();

    // Converters between human and vector forms

    /**
    * @brief Converts a @c FitParameters object into a vector of doubles.
    *
    * This function maps the human-readable fit parameters from a structured
    * representation ( @c FitParameters ) into a flat vector of doubles. It is useful
    * for numerical operations that require array-like data structures.
    *
    * @param ns The FitParameters object to convert.
    * @return A vector of doubles representing the fit parameters.
    * @note Asserts that the size of the resulting vector is exactly 38.
    */
    std::vector<double> ConvertFitParameters(FitParameters p) const;

    /**
    * @brief Converts a @c FitParametersFlag object into a vector of bools.
    *
    * This function takes a set of fit parameter flags and converts them into
    * a vector of booleans. Each flag corresponds to a different fit parameter,
    * indicating whether it should be included in the fitting process.
    *
    * @param ns The FitParametersFlag object to convert.
    * @return A vector of booleans representing the fit parameter flags.
    * @note Asserts that the size of the resulting vector is exactly 38.
    */
    std::vector<bool> ConvertFitParametersFlag(FitParametersFlag pf) const;

    /**
    * @brief Converts a vector of doubles back into a human-readable @c FitParameters object.
    *
    * This function is the inverse of @c ConvertFitParameters, taking a vector of
    * doubles and converting it back into a structured @c FitParameters object.
    *
    * @param vecns The vector of doubles to convert.
    * @return A @c FitParameters object populated with the values from the vector.
    */
    FitParameters ConvertVecToFitParameters(std::vector<double> vecns) const;

    // functions to obtain distributions
    /**
    * @brief Generates a random realization of the expected event distribution based on input fit parameters and a seed.
    *
    * This function uses a random number generator to produce a statistical realization
    * of the event distribution. It uses a weighter constructed from the provided fit parameters
    * to weight the events from the main simulation. The generated sample is then binned to produce a histogram.
    *
    * @param fit_parameters The fit parameters to use for weighting the events.
    * @param seed The seed for the random number generator.
    * @return A multi-dimensional array representing the binned realization of events.
    * @throws @c std::runtime_error If no events are generated, which likely indicates an issue with the expected weights.
    */
    hist_marray GetRealization(std::vector<double> fit_parameters, int seed) const;

    /**
    * @brief Generates a random realization of the expected event distribution based on input fit parameters and a seed.
    *
    * This function is a convenience wrapper that converts the provided fit parameters (in a @c FitParameters object)
    * to a vector of doubles and then calls GetRealization with the converted parameters and a given seed.
    *
    * @param fp The fit parameters to be converted and used.
    * @param seed The seed for the random number generator.
    * @return A multi-dimensional array representing the binned realization of events.
    */
    hist_marray GetRealization(FitParameters fp, int seed) const;

    /**
    * @brief Computes the weighted expectation for fit parameters.
    *
    * This function calculates the expected distribution based on the provided fit parameters.
    * It constructs a weighter from the fit parameters and applies it to the simulation histogram
    * to get the expected distribution, or the expectation.
    *
    * @param fit_parameters The parameters to be used for constructing the weighter.
    * @return A multi-dimensional array, the histogram in energy and zenith for each topology.
    */
    hist_marray GetExpectation(std::vector<double> fit_parameters) const;

    /**
    * @brief Computes the weighted expectation for fit parameters.
    *
    * This function is a convenience wrapper that converts the provided fit parameters 
    * (in a @c FitParameters object)
    * to a vector of doubles and then calls GetExpectation.
    *
    * @param fp The fit parameters to be converted and used.
    * @return A multi-dimensional array, the histogram in energy and zenith for each topology.
    */
    hist_marray GetExpectation(FitParameters fp) const;

    /**
    * @brief Generates a simulated data realization based on input nuisance parameters
    *       using a specified random seed.
    *
    * This function simulates a data realization from a given set of nuisance parameters.
    * It calculates the expected number of events and uses a random number generator
    * to produce a realization with appropriate statistics.
    * Note that the output is an array of events, not a binned histogram.
    *
    * @param fit_parameters A vector of doubles containing the nuisance parameters.
    * @param seed An integer used to seed the random number generator for the realization.
    * @return A two-dimensional array holding all the events. The four columns
    * correspond to the energy, zenith, topology, and weight.
    */
    marray<double,2> GetRealizationEvents(std::vector<double> fit_parameters, int seed) const;

    /**
    * @brief Provides the expected event distribution for a given set of fit parameters.
    *
    * This function calculates the expected distribution of events based on the provided
    * fit parameters. Depending on whether the fast mode is constructed, it uses either
    * the meta simulation or the main simulation to compute the expected event distribution.
    * Note that the output is an array of events, not a binned histogram.
    *
    * @param fit_parameters A vector of doubles containing the fit parameters.
    * @return A two-dimensional array holding all the events. The four columns
    * correspond to the energy, zenith, topology, and weight.
    */
    marray<double,2> GetExpectationEvents(std::vector<double> fit_parameters) const;

    /**
    * @brief Checks if the expected event distribution is valid for the given fit parameters.
    *
    * This is an overloaded function that first converts the human-readable @c FitParameters
    * object into a vector of doubles before calling the main @c CheckExpectation function.
    *
    * @param fit_parameters A FitParameters object containing the fit parameters.
    * @return An integer, where 0 indicates a valid expectation and 1 indicates an invalid one.
    */
    int CheckExpectation(FitParameters fit_parameters) const;

    /**
    * @brief Checks if the expected event distribution is valid for the given fit parameters.
    *
    * This function verifies the validity of the expected event distribution by ensuring
    * that all weighted events have non-negative weights. The function returns an error code
    * if any event has a negative weight.
    *
    * @param fit_params A vector of doubles containing the fit parameters.
    * @return An integer, where 0 indicates a valid expectation and 1 indicates an invalid one.
    */
    int CheckExpectation(std::vector<double> fit_parameters) const;

    /**
    * @brief Computes the square of the weighted expectation for fit parameters.
    *
    * Similar to GetExpectation, but instead of the weight, it uses the square of the weight.
    * This can be useful for certain statistical analyses where the square of weights is required.
    *
    * @param fit_parameters The parameters to be used for constructing the square of the weighter.
    * @return A multi-dimensional array representing the squared expectation distribution.
    */
    hist_marray GetSquareExpectation(std::vector<double> fit_parameters) const;

    /**
    * @brief Computes the square of the weighted expectation for fit parameters.
    *
    * This function is a convenience wrapper that converts the provided fit parameters
    * (in a @c FitParameters object)
    * to a vector of doubles and then calls GetSquareExpectation.
    *
    * @param fit_params The fit parameters to be converted and used.
    * @return A multi-dimensional array representing the squared expected event distribution.
    */
    hist_marray GetSquareExpectation(FitParameters fit_params) const;

    /**
    * @brief Retrieves the bin edges for a specific dimension of a histogram.
    *
    * This function template extracts the bin edges for the specified dimension
    * of a histogram, which is part of a tuple of histograms. The histogram
    * from which to extract the edges is selected by the template parameter.
    *
    * @tparam hist_index The index of the histogram in the tuple from which to pull bin edges.
    * @param dim The dimension of the histogram from which to pull the bin edges.
    * @param hh The tuple containing the histogram.
    * @return A vector containing the edges of the bins in the specified dimension.
    */
    template<unsigned int hist_type>
    std::vector<double> PullBinEdges(int dim, const  HistType& h) const;

    /**
    * @brief General function to calculate a weighted expectation based on an arbitrary function.
    *
    * This function applies a user-provided function to the events in the simulation histogram
    * to calculate the weighted expectation. The function should compute the weight to be applied
    * to each event.
    *
    * @pre @c simulation_histogram_constructed_ must be true, indicating that the simulation histogram
    * has been successfully constructed.
    * @param f A function taking an Event object and returning a weight.
    * @return A multi-dimensional array representing the weighted expectation distribution.
    * @throws @c std::runtime_error If the simulation histogram has not been constructed, or if the
    * expectation is negative.
    */
    hist_marray GetWeightedExpectation(std::function<double(const Event &)> f) const;

    /**
    * @brief Constructs a 2D histogram of energy vs. zenith for the expected events based on fit parameters and topology.
    *
    * This function generates a histogram with energy and zenith dimensions from the main simulation.
    * If topology is non-negative, only events matching the topology are included.
    *
    * @param fit_params The fit parameters to be used in the weighter function.
    * @param topology The index of the topology to filter the events by. If negative, all topologies are included.
    * @return A 2D histogram binned energy and zenith.
    */
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyZenithExpectationHistogram(std::vector<double> fit_params, int topology = -1) const;

    /**
    * @brief Gets a weighter function from the provided fit parameters.
    *
    * This function converts fit parameters to a vector of doubles and uses them
    * to construct a weighter function that can be applied to Event objects.
    *
    * @param fp The fit parameters to be converted into a weighter function.
    * @return A function that takes an Event object and returns a weight based on the input fit parameters.
    */
    std::function<double(const Event&)> GetEventWeighter(FitParameters fp) const;
    
    /**
    * @brief Retrieves the distribution from the constructed data histogram.
    *
    * This function extracts the data distribution from the constructed data histogram and
    * returns it as a multi-dimensional array. Each cell in the array corresponds to the
    * sum of the weights of the events that fall into the corresponding histogram bin.
    *
    * @pre @c data_histogram_constructed_ must be true, indicating that the data histogram
    * has been successfully constructed.
    * @return A multi-dimensional array representing the summed weights in the data histogram.
    * @throws @c std::runtime_error If the data histogram has not been constructed.
    */
    hist_marray GetDataDistribution() const;

    /**
    * @brief Constructs a 2D histogram of energy vs. zenith for the expected events based on fit parameters and topology.
    *
    * This function is a convenience wrapper that converts the provided fit parameters (in a @c FitParameters object)
    * to a vector of doubles and then calls GetEnergyZenithExpectationHistogram with the converted parameters and topology.
    *
    * @param fit_params The fit parameters to be converted and used.
    * @param topology The index of the topology to filter the events by. If negative, all topologies are included.
    * @return A 2D histogram with bins for energy and zenith dimensions.
    */
    histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GetEnergyZenithExpectationHistogram(FitParameters fit_params, int topology = -1) const;

    // Methods to get out and set event samples

    /**
    * @brief Outputs the current data sample as a two-dimensional array.
    *
    * This function creates a two-dimensional array with the contents of the
    * current data sample. Each event's attributes are stored in the array's rows.
    * Note that the output is an array of events, not a binned histogram.
    *
    * @return A two-dimensional array representing the current data sample.
    */
    marray<double,2> GetDataEvents() const;

    /**
    * @brief Generates a simulated data realization based on input nuisance parameters
    *       using a specified random seed.
    *
    * This is an overloaded function that first converts the human-readable @c FitParameters
    * object into a vector of doubles before calling the main @c GetRealizationEvents function.
    * Note that the output is an array of events, not a binned histogram.
    *
    * @param nuisance A FitParameters object containing the nuisance parameters.
    * @param seed An integer used to seed the random number generator for the realization.
    * @return A two-dimensional array representing the simulated data realization.
    */
    marray<double,2> GetRealizationEvents(FitParameters nuisance, int seed) const;

    /**
    * @brief Provides the expected event distribution for a given set of fit parameters.
    *
    * This is an overloaded function that first converts the human-readable @c FitParameters
    * object into a vector of doubles before calling the main @c GetExpectationEvents function.
    * Note that the output is an array of events, not a binned histogram.
    *
    * @param fit_params A FitParameters object containing the fit parameters.
    * @return A two-dimensional array representing the expected event distribution.
    */
    marray<double,2> GetExpectationEvents(FitParameters nuisance) const;

    /**
    * @brief Given data, sums their weights, and updates the internal data histogram.
    *
    * This function processes an array of data events, each with several attributes,
    * and stores them into the internal @c sample_ structure. Additionally, it reconstructs the data
    * histogram based on the new sample. 
    *
    * @param Data A two-dimensional array of data where each row represents an event
    *             with attributes such as energy, zenith, topology, and cached weight.
    *             The input format should be the output of @c GetExpectationEvents .
    * @return The sum of weights of all data events.
    */
    double SetData(marray<double,2> Data);

    // Methods to get histogram binning
    /**
    * @brief Gets the energy bin edges from the Monte Carlo simulation histogram.
    *
    * This function uses the @c PullBinEdges template function to extract the
    * bin edges for the energy dimension (dimension 0) from the first histogram
    * in the simulation histogram tuple.
    *
    * @return A vector of doubles containing the energy bin edges of the MC histogram.
    */
    std::vector<double> GetEnergyBinsMC() const;

    /**
    * @brief Gets the zenith bin edges from the Monte Carlo simulation histogram.
    *
    * This function uses the @c PullBinEdges template function to extract the
    * bin edges for the zenith dimension (dimension 1) from the first histogram
    * in the simulation histogram tuple.
    *
    * @return A vector of doubles containing the zenith bin edges of the MC histogram.
    */
    std::vector<double> GetZenithBinsMC() const;

    /**
    * @brief Gets the topology bin edges from the Monte Carlo simulation histogram.
    *
    * This function uses the @c PullBinEdges template function to extract the
    * bin edges for the topology dimension (dimension 2) from the first histogram
    * in the simulation histogram tuple.
    *
    * @return A vector of doubles containing the topology bin edges of the MC histogram.
    */
    std::vector<double> GetTopologyBinsMC() const;

    // functions to check the status of the object

    /**
    * @brief Checks if the data is loaded.
    * @return The value of @c data_loaded_.
    */
    bool CheckDataLoaded() const                       {return data_loaded_;};
    /**
    * @brief Checks if the simulation is loaded.
    * @return The value of @c simulation_loaded_.
    */
    bool CheckSimulationLoaded() const                 {return simulation_loaded_;};
    /**
    * @brief Checks if cross section weighter is constructed.
    * @return The value of @c xs_weighter_constructed_.
    */
    bool CheckCrossSectionWeighterConstructed() const  {return xs_weighter_constructed_;};
    /**
    * @brief Checks if flux weighter is constructed.
    * @return The value of @c flux_weighter_constructed_.
    */
    bool CheckFluxWeighterConstructed() const          {return flux_weighter_constructed_;};
    /**
    * @brief Checks if lepton weighter is constructed.
    * @return The value of @c lepton_weighter_constructed_.
    */
    bool CheckLeptonWeighterConstructed() const        {return lepton_weighter_constructed_;};
    /**
    * @brief Checks if data histogram is constructed.
    * @return The value of @c data_histogram_constructed_.
    */
    bool CheckDataHistogramConstructed() const         {return data_histogram_constructed_;};
    /**
    * @brief Checks if simulation histogram is constructed.
    * @return The value of @c simulation_histogram_constructed_.
    */
    bool CheckSimulationHistogramConstructed() const   {return simulation_histogram_constructed_;};
    /**
    * @brief Checks if likelihood problem is constructed.
    * @return The value of @c likelihood_problem_constructed_.
    */
    bool CheckLikelihoodProblemConstruction() const    {return likelihood_problem_constructed_;};
    
    /**
    * @brief Reports the status of various components of the @c GollumFit object.
    *
    * This function prints to the standard output the status of data and simulation
    * loading, as well as the construction status of various weighters and histograms.
    * It is used for debugging and verification purposes.
    */
    void ReportStatus() const;

    // functions to evaluate the likelihood
    /**
    * @brief Evaluates the negative log-likelihood for a given set of nuisance parameters.
    *
    * This function computes the likelihood of the nuisance parameters optionally including the prior.
    * It is a wrapper that takes a vector of doubles as the nuisance parameters, and calls the 
    * @c evaluateLikelihood method of an object declared with @c phys_tools::likelihood::makeLikelihoodProblem .
    *
    * @param fp The vector of nuisance parameters to evaluate.
    * @param include_prior A boolean flag to include the prior in the likelihood evaluation.
    * @return The negative log-likelihood value.
    * @throws @c std::runtime_error If the likelihood problem has not been constructed.
    */
    double EvalLLH(std::vector<double> fp, bool include_prior) const;
    double operator()(std::vector<double> fp) const {return EvalLLH(fp,true);}

    /**
    * @brief Evaluates the negative log-likelihood for a given set of nuisance parameters.
    *
    * This function is an overloaded version that takes a FitParameters object, converts it to
    * a vector of doubles, and calls the primary EvalLLH function.
    *
    * @param fp The nuisance parameters to evaluate, encapsulated in a @c FitParameters object.
    * @param include_prior A boolean flag to include the prior in the likelihood evaluation.
    * @return The negative log-likelihood value.
    */
    double EvalLLH(FitParameters fp, bool include_prior) const;
    double operator()(FitParameters fp) const {return EvalLLH(fp,true);}

    /**
    * @brief Computes the gradient of the negative log-likelihood function.
    *
    * This function evaluates the gradient of the likelihood using automatic differentiation.
    * It is used for gradient-based optimization algorithms.
    *
    * @param v A vector containing the nuisance parameters wrapped in autodiff variables.
    * @return The gradient of the negative log-likelihood with respect to the nuisance parameters.
    */
    phys_tools::autodiff::FD<38> EvalLLHGradient(std::vector<phys_tools::autodiff::FD<38>> v) const;

    /**
    * @brief Minimize the negative log-likelihood for the set problem.
    *
    * This function performs the minimization of the negative log-likelihood 
    * over the parameter space defined by the user. It iterates over a set of seeds
    * and uses the LBFGSB algorithm to find the minimum. The function ensures that
    * the likelihood problem is constructed before proceeding with the optimization
    * and throws an exception if this is not the case.
    * The optimization process involves declaring the fit parameters with their initial
    * values, boundaries, and fixing certain parameters if required. 
    * 
    * @note A central @c phys_tools::lbfgsb::LBFGSB_Driver object is declared and minimization is performed. 
    * The results and relevant metadata are stored to a @c FitResult object.
    * 
    * @pre @c data_histogram_constructed_ must be true, indicating that the data histogram
    * has been successfully constructed.
    *
    * @return FitResult object containing the results of the fit. This includes the
    * success of the fit, the minimum NLL value, the best-fit parameters, and the
    * number of function and gradient evaluations performed during the fit.
    *
    * @throws @c std::runtime_error If the likelihood problem is not constructed before
    * the function is called.
    */
    FitResult MinLLH() const;

    // get/set functions
    /**
    * @brief Steering parameters getter function. 
    * @return The current @c SteeringParams object.
    */
    SteeringParams       GetSteeringParams()     { return steeringParams_;};
    /**
    * @brief Data paths getter function. 
    * @return The current @c DataPaths object.
    */
    DataPaths            GetDataPaths()          { return dataPaths_;};
    /**
    * @brief Fit parameter seed getter function. 
    * @return The current seed as a @c std::vector<FitParameters> object.
    */
    std::vector<FitParameters>        GetFitParametersSeed() { return fitSeed_;};
    /**
    * @brief Fit priors getter function. 
    * @return The current seed as a @c Priors object.
    */
    Priors               GetFitPriors()          { return priors_;};
    /**
    * @brief Fit parameters flag getter function.
    * @return The current seed as a @c FitParametersFlag object.
    */
    FitParametersFlag    GetFitParametersFlag()  { return fixedParams_;};
    /**
    * @brief Fit parameters bounds getter function.
    * @return The current seed as a @c FitParametersBound object.
    */
    FitParametersBound   GetFitParametersBound() { return boundParams_;};

    /**
    * @brief Sets the steering parameters.
    * 
    * This method updates the steering parameters and re-constructs the @c DFWM WeighterMaker object
    * with the new steering params and the existing splines and paths.
    * 
    * @param p The new @c SteeringParams object to set.
    */
    void SetSteeringParams(SteeringParams p) { 
      steeringParams_=p; 
      DFWM = sterile::WeighterMaker(domefficiencySplines_,
                                    holeIceSplines_,
                                    attenuationSplines_,
                                    p,
                                    dataPaths_);
    }
    /**
    * @brief Sets the seed for fit parameters.
    * 
    * This method updates the seeds for the fit parameters. It also ensures that
    * the seeds meet certain criteria by calling @c ForceFitSeedSanity and marks the
    * fit seed as constructed.
    * 
    * @param fp A vector of @c FitParameters objects to set as the new seeds.
    */
    void SetFitParametersSeed(std::vector<FitParameters> fp) { 
      fitSeed_= fp; 
      ForceFitSeedSanity(); 
      fitSeed_constructed_ = true;
    }
    /**
    * @brief Sets the flags for fit parameters.
    * 
    * This method updates the flags that determine which fit parameters should
    * be held fixed during the fitting process and marks the fixed parameters
    * as constructed.
    * 
    * @param fpg The new @c FitParametersFlag object to set.
    */
    void SetFitParametersFlag(FitParametersFlag fpg) { 
      fixedParams_ = fpg; 
      fixedParams_constructed_ = true;
    }
    /**
    * @brief Sets the bounds for fit parameters.
    * 
    * This method updates the bounds for the fit parameters and marks the bounds
    * as constructed.
    * 
    * @param fpb The new @c FitParametersBound object to set.
    */
    void SetFitParametersBound(FitParametersBound fpb) { 
      boundParams_ = fpb; 
      boundParams_constructed_ = true;
    }
    /**
    * @brief Sets the priors for the fit parameters.
    * 
    * This method updates the priors for the fit parameters and marks the priors
    * as constructed.
    * 
    * @param priors The new @c Priors object to set.
    */
    void SetFitParametersPriors(Priors priors) { 
      priors_ = priors; 
      priors_constructed_ = true;
    }

};

} // close namespace SterileSearch

#endif /* GOLLUMFIT_H_ */
