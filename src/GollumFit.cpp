#include <cmath>
#include <algorithm>
#include <iostream>

#include "GollumFit.h"
#include "FastMode.h"
#include "GollumMCSpecifications.h"

namespace gollumfit {

/*************************************************************************************************************
 * Constructor
 * **********************************************************************************************************/

GollumFit::GollumFit(DataPaths dataPaths, SteeringParams steeringParams){
  ReConfig<true>(std::shared_ptr<DataPaths>(new DataPaths(dataPaths)),
      std::shared_ptr<SteeringParams>(new SteeringParams(steeringParams)));
}

/*************************************************************************************************************
 * Implementation auxiliary functions
 * **********************************************************************************************************/

/**
* @brief Auxiliary function to bin events into a histogram.
*
* This function takes an event and bins it into the provided histogram structure.
* It performs several assertions on the event data to ensure the validity of the
* event before adding it to the histogram.
*
* @param h Reference to a HistType object. The histogram to bin into.
* @param e Constant reference to an Event object
* @throws std::logic_error If any of the assertions fail, indicating the event data is invalid.
*/
void binner (HistType& h, const Event& e) {
    
    assert(e.topology >= 0);
    assert(!std::isnan(e.energy));
    assert(!std::isnan(e.zenith));
    assert(!std::isnan(e.topology));
    std::get<0>(h).add(e.energy,cos(e.zenith),e.topology,amount(std::cref(e)));

};

/*************************************************************************************************************
 * Functions to read and write data
 * **********************************************************************************************************/

void GollumFit::LoadData(){

  using namespace phys_tools::tableio;
  using namespace H5Load::sterile;
  try{
    auto dataAction = [&](RecordID id, Event& e){
      e.cachedWeight=1.;
      e.topology = (unsigned int)Topology::track;// its a track
      sample_.push_back(e);
    };
    auto ic86Action=[&](RecordID id, Event& e){ dataAction(id,e); };
    readFile(CheckedFilePath(dataPaths_.data_path+"IC86.h5"),"MuExEnergy",-1,ic86Action);
  } catch(std::exception& ex){
    std::cerr << "Problem loading experimental data: " << ex.what() << std::endl;
  }
  std::cout << "Loaded " << sample_.size() << " experimental events" << std::endl;

  data_loaded_ = true;
  if(not data_loaded_)
    throw std::runtime_error("Sample not found.");

}


void GollumFit::LoadMC(){

  std::cout << "Begin Loading MC" << std::endl;

  using namespace phys_tools::tableio;

  mainSimulation_.clear();
  metaEvents_.clear();

  double livetime=steeringParams_.fullLivetime;

  std::vector<MCSet> simSetsToLoad = sterile::GetSimulationSets(steeringParams_.simToLoad, dataPaths_.mc_path+ "/STERILE");

  auto simAction=[&](RecordID id, Event& e, double number_of_files, const DOMEfficiencySetter<Event>& domEff, const HoleIceSetter<Event>& holeIce){

    try{

      LW::Event lw_e {
        e.primaryType,
        e.final_state_particle_0,
        e.final_state_particle_1,
        e.intX,
        e.intY,
        e.primaryEnergy,
        e.primaryAzimuth,
        e.primaryZenith,
        0.0,//e.x,
        0.0,//e.y,
        0.0,//e.z,
        0.0,//e.r,
        e.totalColumnDepth
      };

      double Scaling = livetime/number_of_files;

      e.cachedConvWeight=(*convFluxWeighter_)(lw_e)*Scaling;
      e.cachedPromptWeight=(*promptFluxWeighter_)(lw_e)*Scaling;
      e.cachedAstroWeight=(*astroFluxWeighter_)(lw_e)*Scaling;

      if (domeff_spline_loaded_) domEff.setCache(e);
      if (holeice_spline_loaded_) holeIce.setCache(e);
      if (ice_gradient_spline_loaded_) {
        double reco_coordinates[2]={cos(e.zenith),log10(e.energy)};
        Topology topo = static_cast<gollumfit::Topology>(e.topology);
        std::pair<int,Topology> ig0(0,topo);
        std::pair<int,Topology> ig1(1,topo);
        std::pair<int,Topology> ig2(2,topo);
        std::pair<int,Topology> ig3(3,topo);
        std::pair<int,Topology> ig4(4,topo);
        std::pair<int,Topology> ig5(5,topo);
        std::pair<int,Topology> ig6(6,topo);
        std::pair<int,Topology> ig7(7,topo);
        std::pair<int,Topology> ig8(8,topo);
        e.cachedIceGrad0 = (*(*iceGradientSplines_.find(ig0)).second)(reco_coordinates);
        e.cachedIceGrad1 = (*(*iceGradientSplines_.find(ig1)).second)(reco_coordinates);
        e.cachedIceGrad2 = (*(*iceGradientSplines_.find(ig2)).second)(reco_coordinates);
        e.cachedIceGrad3 = (*(*iceGradientSplines_.find(ig3)).second)(reco_coordinates);
        e.cachedIceGrad4 = (*(*iceGradientSplines_.find(ig4)).second)(reco_coordinates);
        e.cachedIceGrad5 = (*(*iceGradientSplines_.find(ig5)).second)(reco_coordinates);
        e.cachedIceGrad6 = (*(*iceGradientSplines_.find(ig6)).second)(reco_coordinates);
        e.cachedIceGrad7 = (*(*iceGradientSplines_.find(ig7)).second)(reco_coordinates);
        e.cachedIceGrad8 = (*(*iceGradientSplines_.find(ig8)).second)(reco_coordinates);
      }

      double true_coordinates[2]={log10(e.primaryEnergy),cos(e.primaryZenith)};
      if (atmospheric_density_spline_loaded_) e.cachedAtmDensity = (*atmosphericDensityUncertaintySpline_)(true_coordinates);
      if (atmospheric_kaonlosses_spline_loaded_) e.cachedKaonLosses = (*kaonLossesUncertaintySpline_)(true_coordinates);

      double oneWeightScaling = convFluxWeighter_->get_oneweight(lw_e)*Scaling;

      if (hadronic_spline_loaded_) {
        for(const auto& entry : hadronicFluxWeighter_){
          switch (entry.first) {
            case HadronicParameter::HEkp: 
              e.cachedHadronicHEkp = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::HEkm: 
              e.cachedHadronicHEkm = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE1pip: 
              e.cachedHadronicVHE1pip = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE1pim: 
              e.cachedHadronicVHE1pim = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3kp: 
              e.cachedHadronicVHE3kp = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3km: 
              e.cachedHadronicVHE3km = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3pip: 
              e.cachedHadronicVHE3pip = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3pim: 
              e.cachedHadronicVHE3pim = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3p: 
              e.cachedHadronicVHE3p = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case HadronicParameter::VHE3n: 
              e.cachedHadronicVHE3n = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            default:
              throw std::runtime_error("Impossible hadronic parameter label");
          }
        }
      }
      if (cosmic_ray_spline_loaded_) {
        for(const auto& entry : cosmicRayFluxWeighter_){
          switch (entry.first) {
            case CosmicRayParameter::GSF1 : 
              e.cachedCosmicRay1 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case CosmicRayParameter::GSF2 : 
              e.cachedCosmicRay2 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case CosmicRayParameter::GSF3 : 
              e.cachedCosmicRay3 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case CosmicRayParameter::GSF4 : 
              e.cachedCosmicRay4 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case CosmicRayParameter::GSF5 : 
              e.cachedCosmicRay5 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            case CosmicRayParameter::GSF6 : 
              e.cachedCosmicRay6 = (*(entry.second))(lw_e)*oneWeightScaling - e.cachedConvWeight;
              break;
            default:
              throw std::runtime_error("Impossible cosmicray parameter label");
          }
        }
      }

      mainSimulation_.push_back(e);

    } catch (std::exception & exception){
      std::cout << exception.what() << " Encountered impossible to weight event. Removing it from deque. CAD." << std::endl;
      std::cout << e << std::endl;
    }

  };

  for(auto simSet : simSetsToLoad){

    if(domeff_spline_loaded_) domEffSetter_=new DOMEfficiencySetter<Event>(domefficiencySplines_,simSet.unshadowedFraction);
    if(holeice_spline_loaded_) holeIceSetter_=new HoleIceSetter<Event>(holeIceSplines_,simSet.holeiceForward);

    auto callback=[&](RecordID id, Event& e){ simAction(id,e,simSet.number_of_files,*domEffSetter_,*holeIceSetter_); };
    
    std::cout << "Number of MC sets and splits: " << std::get<0>(simSet.split) << " " << std::get<1>(simSet.split) << std::endl;
    if(std::get<0>(simSet.split)){
      for(unsigned int split = 0; split < std::get<1>(simSet.split); split++){
        auto path=CheckedFilePath(simSet.path+"/"+simSet.filename+"/"+simSet.filename+"_"+std::to_string(split)+".h5");
        H5Load::sterile::readFile(path,steeringParams_.energyName,steeringParams_.selectionStart,callback);
      }
    }
    else {
      auto path=CheckedFilePath(simSet.path+"/"+simSet.filename);
      H5Load::sterile::readFile(path,steeringParams_.energyName,steeringParams_.selectionStart,callback);      
    }
  }

  std::cout << "Loaded " << mainSimulation_.size() << " events in main simulation set" << std::endl;
  assert(not mainSimulation_.empty());

  std::cout << "End Loading MC" << std::endl;

  simulation_loaded_ = true;
  if(not simulation_loaded_)
    throw std::runtime_error("Tag simulation not found.");
}

/*************************************************************************************************************
 * Functions to do load compact
 * **********************************************************************************************************/

bool GollumFit::WriteCompact(string file_path) const {
  // ALEJO: be aware that the sideband MC and data are not dump on the compact format.
  std::string file_path_compact = file_path + "/STERILE.meows";
  gollumfit::dump::splatData(file_path_compact,0,sample_,(fastmode_constructed_)?metaEvents_:mainSimulation_);
  return true;
}

void GollumFit::LoadCompact() {
  // ALEJO: be aware that the sideband MC and data are not dump on the compact format.
  std::string file_path_compact = dataPaths_.compact_file_path + "/STERILE.meows";
  gollumfit::dump::unsplatData(CheckedFilePath(file_path_compact),0,sample_,mainSimulation_);
  data_loaded_ = true;
  simulation_loaded_ = true;
}

void GollumFit::ClearData(){
  sample_.clear();
}

void GollumFit::ClearSimulation(){
  mainSimulation_.clear();
  metaEvents_.clear();
}

/*************************************************************************************************************
 * Functions to load to load DOM efficiency splines
 * **********************************************************************************************************/

bool GollumFit::LoadDOMEfficiencySplines(){
  domefficiencySplines_.clear();
  std::cout << "Loading DOM efficiency splines..." << std::endl;

  auto register_domefficiency_spline = [&](FluxComponent component, Topology topology) {
    domefficiencySplines_.insert({std::make_pair(component,topology),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.domeff_spline_path+"/domefficiency_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits")))});
  };

  for ( auto flux : {FluxComponent::atmConv,FluxComponent::atmPrompt,FluxComponent::diffuseAstro} ) {
    if(steeringParams_.selectionStart != -1){
      register_domefficiency_spline(flux,Topology::track);
      register_domefficiency_spline(flux,Topology::shower);
    }
    else {
      register_domefficiency_spline(flux,Topology::all);
    }
  }

  assert(domefficiencySplines_.size() != 0);

  // assuming all splines have the same validity ranges
  domEffSpline_minValidValue_ = domefficiencySplines_.begin()->second->lower_extent(2);
  domEffSpline_maxValidValue_ = domefficiencySplines_.begin()->second->upper_extent(2);

  std::cout << "DOM efficiency splines validity range set to [" + std::to_string(domEffSpline_minValidValue_) + ", " + std::to_string(domEffSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetDOMEfficiencySplines(domefficiencySplines_);

  return true;

}

/*************************************************************************************************************
 * Functions to load to load Hole Ice splines
 * **********************************************************************************************************/

bool GollumFit::LoadHoleIceResources(){
  holeIceSplines_.clear();
  std::cout << "Loading Hole Ice splines..." << std::endl;

  auto register_holeice_spline = [&](FluxComponent component, Topology topology) {
    holeIceSplines_.insert({std::make_pair(component,topology),
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.holeice_spline_path+"/holeice_spline_stacked_"+GetFluxComponentName(component)+"_"+GetTopologyName(topology)+".fits")))});
  };


  for ( auto flux : {FluxComponent::atmConv,FluxComponent::atmPrompt,FluxComponent::diffuseAstro} ) {
    if(steeringParams_.selectionStart != -1){
      register_holeice_spline(flux,Topology::track);
      register_holeice_spline(flux,Topology::shower);
    }
    else {
      register_holeice_spline(flux,Topology::all);
    }
  }
  
  assert(holeIceSplines_.size() != 0);

  // assuming all splines have the same validity ranges
  holeIceSpline_minValidValue_ = holeIceSplines_.begin()->second->lower_extent(2);
  holeIceSpline_maxValidValue_ = holeIceSplines_.begin()->second->upper_extent(2);

  std::cout << "Hole Ice splines forward parameter validity range set to [" + std::to_string(holeIceSpline_minValidValue_) + ", " + std::to_string(holeIceSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetHoleIceSplines(holeIceSplines_);

  return true;

}

/*************************************************************************************************************
 * Functions to load to load Attenuation splines
 * **********************************************************************************************************/

bool GollumFit::LoadAttenuationResources(){
  attenuationSplines_.clear();
  std::cout << "Loading attenuation splines..." << std::endl;

  auto register_attenuation_spline = [&](FluxComponent component, LW::ParticleType particle_type) {
    attenuationSplines_.insert({{component, particle_type},
        std::shared_ptr<splinetable<>>(new splinetable<>(CheckedFilePath(dataPaths_.attenuation_spline_path+"/attenuation_spline_"+GetFluxComponentName(component)+"_"+GetParticleName(particle_type)+".fits")))});
  };


  for ( auto flux : {FluxComponent::atmConv,FluxComponent::atmPrompt,FluxComponent::diffuseAstro} ) {
    for ( auto part : {LW::ParticleType::NuMu,LW::ParticleType::NuMuBar,LW::ParticleType::NuTau,LW::ParticleType::NuTauBar} ) {
      register_attenuation_spline(flux, part);
    }
  }
  assert(attenuationSplines_.size() == 12);

  // assuming all splines have the same validity ranges
  attenuationSpline_minValidValue_ = attenuationSplines_.begin()->second->lower_extent(2);
  attenuationSpline_maxValidValue_ = attenuationSplines_.begin()->second->upper_extent(2);

  std::cout << "Attenuation splines forward parameter validity range set to [" + std::to_string(attenuationSpline_minValidValue_) + ", " + std::to_string(attenuationSpline_maxValidValue_) + "]."<< std::endl;

  DFWM.SetAttenuationSplines(attenuationSplines_);

  return true;

}

/*************************************************************************************************************
 * Functions to load to load icegradients splines
 * **********************************************************************************************************/

bool GollumFit::LoadIceGradientResources(){
  iceGradientSplines_.clear();
  std::cout << "Loading Ice Gradients resources..."<< std::endl;

  auto register_icegrad_spline = [&](int idx, Topology topology) {
    iceGradientSplines_.insert({{idx,topology},std::make_shared<splinetable<>>(CheckedFilePath(dataPaths_.ice_gradient_spline_path+ "/" + steeringParams_.ice_gradient_filename[idx] + "_" + GetTopologyName(topology) + ".fits"))});
  };
  
  for ( unsigned int i=0; i<steeringParams_.ice_gradient_filename.size(); i++ ) {
    if(steeringParams_.selectionStart != -1){
      for ( auto topology : {Topology::track,Topology::shower} ) register_icegrad_spline(i, topology); 
    }
    else {
      register_icegrad_spline(i, Topology::all); 
    }
  }

  return true;

}

/*************************************************************************************************************
 * Functions to load to load Flux splines
 * **********************************************************************************************************/

bool GollumFit::LoadKaonEnergyLossesUncertaintiesResources(){
  std::cout << "Loading Kaon Energy Losses Uncertainty resources..." << std::endl;
  kaonLossesUncertaintySpline_ = std::make_shared<splinetable<>>(CheckedFilePath(dataPaths_.atmospheric_kaonlosses_spline_path));
  return true;
}

bool GollumFit::LoadAtmosphericDensityUncertaintiesResources(){
  std::cout << "Loading Atmospheric Density Uncertainty resources..."<< std::endl;
  atmosphericDensityUncertaintySpline_ = std::make_shared<splinetable<>>(CheckedFilePath(dataPaths_.atmospheric_density_spline_path));
  return true;
}

bool GollumFit::LoadHadronicUncertaintiesResources(){
  hadronicFluxWeighter_.clear();
  std::cout << "Loading Hadronic Uncertainty resources..."<< std::endl;

  std::string model_label = steeringParams_.model_label;

  auto register_flux = [&](string had_name) {
    hadronicFluxWeighter_.insert({GetHadronicParameter(had_name),
        std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.hadronic_spline_path+"/" + had_name + "_" + model_label + ".hdf5"))});
  };

  for ( auto had_name : steeringParams_.active_hadronic_parameters ) register_flux(had_name);            

  return true;

}


bool GollumFit::LoadCosmicRayUncertaintiesResources(){
  cosmicRayFluxWeighter_.clear();
  std::cout << "Loading Cosmic Ray Uncertainty resources..."<< std::endl;

  std::string model_label = steeringParams_.model_label;

  auto register_flux = [&](string cr_name) {
    cosmicRayFluxWeighter_.insert({GetCosmicRayParameter(cr_name),
        std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.cosmic_ray_spline_path+"/" + cr_name + "_" + model_label + ".hdf5"))});
  };

  for ( auto cr_name : steeringParams_.active_cosmicray_parameters ) register_flux(cr_name);            

  return true;

}


/*************************************************************************************************************
 * Functions to construct weighters
 * **********************************************************************************************************/

void GollumFit::ConstructCrossSectionWeighter(){
  xsw_ = std::make_shared<LW::CrossSectionFromSpline>(static_cast<std::string>(CheckedFilePath(dataPaths_.diff_neutrino_cc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_antineutrino_cc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_neutrino_nc_xs_spline_path)),
                                                      static_cast<std::string>(CheckedFilePath(dataPaths_.diff_antineutrino_nc_xs_spline_path)));
  xs_weighter_constructed_=true;
}

void GollumFit::ConstructFluxWeighter(){
  fluxConv_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.conventional_nusquids_atmospheric_file), false);
  fluxPrompt_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.prompt_nusquids_atmospheric_file), false);
  fluxAstro_ = std::make_shared<LW::nuSQUIDSAtmFlux<>>(CheckedFilePath(dataPaths_.astro_nusquids_file));

  // Enable averaged flavor evaluation for astro flux weighter
  if (steeringParams_.astroOscAvgScale.has_value()) {
    fluxAstro_->EnableAveragedEval(steeringParams_.astroOscAvgScale.value());
    std::cout << "Astro flux weighter will use averaged oscillation evaluation above " << steeringParams_.astroOscAvgScale.value() << " eV." << std::endl;
  } else {
    std::cout << "Warning: Astro flux weighter will NOT use averaged oscillation evaluation." << std::endl;
  }
  
  flux_weighter_constructed_= true;
  simulation_initialized_ = false;
}

void GollumFit::ConstructMonteCarloGenerationWeighter(){
  mcw_.clear();

  std::vector<MCSet> simSetsToLoad = sterile::GetSimulationSets(steeringParams_.simToLoad, dataPaths_.mc_path+ "/STERILE");
  for(auto simset: simSetsToLoad ){
    for(auto g : simset.generators){
      mcw_.emplace_back(g);
    }
  }

  assert(not mcw_.empty());
  mc_generation_weighter_constructed_=true;
}

void GollumFit::ConstructLeptonWeighter(){
  if(not mc_generation_weighter_constructed_)
    throw std::runtime_error("MonteCarlo generation weighter has to be constructed first.");
  if(not flux_weighter_constructed_)
    throw std::runtime_error("Flux weighter has to be constructed first.");
  if(not xs_weighter_constructed_)
    throw std::runtime_error("Cross section weighter has to be constructed first.");
  convFluxWeighter_ = std::make_shared<LW::Weighter>(fluxConv_,xsw_,mcw_);
  promptFluxWeighter_ = std::make_shared<LW::Weighter>(fluxPrompt_,xsw_,mcw_);
  astroFluxWeighter_ = std::make_shared<LW::Weighter>(fluxAstro_,xsw_,mcw_);
  lepton_weighter_constructed_=true;
}


/*************************************************************************************************************
 * Functions to construct histograms
 * **********************************************************************************************************/

void GollumFit::ConstructDataHistogram(){
  if(not data_loaded_)
    throw std::runtime_error("No data has been loaded. Cannot construct data histogram.");

  typedef std::remove_reference<decltype(std::get<0>(dataHist_))>::type Hist0;

  Hist0 h0(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                       LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth), // zenith dimension
                       LinearAxis(0,1)); // topology dimension

  dataHist_ = std::make_tuple(h0);

  auto& data0 = std::get<0>(dataHist_);
  data0.getAxis(0)->setLowerLimit(steeringParams_.minFitEnergy);
  data0.getAxis(0)->setUpperLimit(steeringParams_.maxFitEnergy);
  data0.getAxis(1)->setLowerLimit(steeringParams_.minCosth);
  data0.getAxis(1)->setUpperLimit(steeringParams_.maxCosth);

  // fill in the histogram with the data
  bin(sample_, dataHist_, binner);

  data_histogram_constructed_=true;
}

void GollumFit::ConstructSimulationHistogram(){
  if(not simulation_loaded_)
    throw std::runtime_error("No simulation has been loaded. Cannot construct simulation histogram.");

  typedef std::remove_reference<decltype(std::get<0>(simHist_))>::type Hist0;

  Hist0 h0(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                       LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth), // zenith dimension
                       LinearAxis(0,1)); // topology dimension

  simHist_ = std::make_tuple(h0);

  auto& sim0 = std::get<0>(simHist_);
  sim0.getAxis(0)->setLowerLimit(steeringParams_.minFitEnergy);
  sim0.getAxis(0)->setUpperLimit(steeringParams_.maxFitEnergy);
  sim0.getAxis(1)->setLowerLimit(steeringParams_.minCosth);
  sim0.getAxis(1)->setUpperLimit(steeringParams_.maxCosth);

  bin(mainSimulation_, simHist_, binner);

  simulation_histogram_constructed_=true;
}

/*************************************************************************************************************
 * Functions to obtain distributions
 * **********************************************************************************************************/

hist_marray GollumFit::GetDataDistribution() const {
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before asking for it.");

  const auto& dataHist = std::get<0>(dataHist_);

  marray<double,3> array {static_cast<size_t>(dataHist.getBinCount(2)),
                          static_cast<size_t>(dataHist.getBinCount(1)),
                          static_cast<size_t>(dataHist.getBinCount(0))};

  for(size_t it=0; it<dataHist.getBinCount(2); it++){ // topology
    for(size_t ic=0; ic<dataHist.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<dataHist.getBinCount(0); ie++){ // energy
        auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(dataHist(ie,ic,it));
        array[it][ic][ie] = 0;
        for(Event event : itc){
          array[it][ic][ie] += event.cachedWeight;
        }
      }
    }
  }

  return array;
}


hist_marray GollumFit::GetExpectation(std::vector<double> fit_parameters) const {
  auto weighter = DFWM(fit_parameters);
  std::function<double(const Event&)> f = [&weighter](const Event & e){return weighter(e);};
  return GetWeightedExpectation(f);
}


hist_marray GollumFit::GetSquareExpectation(std::vector<double> fit_parameters) const {
  auto weighter = DFWM(fit_parameters);
  std::function<double(const Event &)> f = [&weighter](const Event & e){return pow(weighter(e),2.);};
  return GetWeightedExpectation(f);
}

//AW: write a similar function as this one that uses a ML prediction for the expectation instead
hist_marray GollumFit::GetWeightedExpectation(std::function<double(const Event &)> f) const {
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before asking for distributions.");

  const auto& simHist = std::get<0>(simHist_);

  marray<double,3> array {static_cast<size_t>(simHist.getBinCount(2)),
                          static_cast<size_t>(simHist.getBinCount(1)),
                          static_cast<size_t>(simHist.getBinCount(0))};
  std::fill(array.begin(),array.end(),0);

  for(size_t it=0; it<simHist.getBinCount(2); it++){ // topology
    for(size_t ic=0; ic<simHist.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<simHist.getBinCount(0); ie++){ // energy
        auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(simHist(ie,ic,it));
        double expectation=0;
        for(auto event : itc.entries()){
          expectation+=f(event);
        }
        assert(expectation>=0.0 && "Expectation cannot be negative");

        array[it][ic][ie] = expectation;
      }
    }
  }
 
  return array;
}

hist_marray GollumFit::GetExpectation(FitParameters fit_params) const {
  return GetExpectation(ConvertFitParameters(fit_params));
}

hist_marray GollumFit::GetSquareExpectation(FitParameters fit_params) const {
  return GetSquareExpectation(ConvertFitParameters(fit_params));
}

hist_marray GollumFit::GetRealization(std::vector<double> fit_params, int seed) const {
  std::mt19937 rng;
  rng.seed(seed);

  std::cout << "construct weighter" << std::endl;
  auto weighter=DFWM(fit_params);

  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization=phys_tools::likelihood::generateSample(weights,mainSimulation_,expected,rng);
  auto fullRealizationHist = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(dataHist_)));
  bin(realization,fullRealizationHist,binner);
  auto& realizationHist = std::get<0>(fullRealizationHist);

  if(realization.size() == 0){
    throw std::runtime_error("No events generated. Expected events are "+std::to_string(expected));
  }

  marray<double,3> array {static_cast<size_t>(realizationHist.getBinCount(2)),
                          static_cast<size_t>(realizationHist.getBinCount(1)),
                          static_cast<size_t>(realizationHist.getBinCount(0))};
  std::fill(array.begin(),array.end(),0);

  for(size_t it=0; it<realizationHist.getBinCount(2); it++){ // topology
    for(size_t ic=0; ic<realizationHist.getBinCount(1); ic++){ // zenith
      for(size_t ie=0; ie<realizationHist.getBinCount(0); ie++){ // energy
        auto itc = static_cast<phys_tools::likelihood::entryStoringBin<std::reference_wrapper<const Event>>>(realizationHist(ie,ic,it));
        array[it][ic][ie] = itc.size();
      }
    }
  }

  return array;
}

hist_marray GollumFit::GetRealization(FitParameters nuisance, int seed) const {
  return GetRealization(ConvertFitParameters(nuisance),seed);
}

histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GollumFit::GetEnergyZenithExpectationHistogram(std::vector<double> fit_params, int topology) const {
  auto EnergyZenithHistogram = histogram<2,entryStoringBin<std::reference_wrapper<const Event>>>(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth), // energy dimension
                                                                                                 LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth)); // zenith dimension

  for(auto & e : mainSimulation_){
    if(topology<0)
      EnergyZenithHistogram.add(e.energy,cos(e.zenith),amount(std::cref(e)));
    else if (e.topology == (unsigned int)(topology))
      EnergyZenithHistogram.add(e.energy,cos(e.zenith),amount(std::cref(e)));
  }

  return EnergyZenithHistogram;
}

histogram<2,entryStoringBin<std::reference_wrapper<const Event>>> GollumFit::GetEnergyZenithExpectationHistogram(FitParameters fit_params, int topology) const {
  return GetEnergyZenithExpectationHistogram(ConvertFitParameters(fit_params),topology);
}

std::function<double(const Event&)> GollumFit::GetEventWeighter(FitParameters fp) const {
  return DFWM(ConvertFitParameters(fp));
}

/*************************************************************************************************************
 * Functions to construct likelihood problem and evaluate it
 * **********************************************************************************************************/

void GollumFit::ConstructLikelihoodProblem(){
  if(not data_histogram_constructed_)
    throw std::runtime_error("Data histogram needs to be constructed before likelihood problem can be formulated.");
  if(not simulation_histogram_constructed_)
    throw std::runtime_error("Simulation histogram needs to be constructed before likelihood problem can be formulated.");
  if(not fitSeed_constructed_)
    throw std::runtime_error("Fitting seed needs to be constructed before likelihood problem can be formulated.");
  if(not priors_constructed_)
    throw std::runtime_error("Priors needs to be constructed before likelihood problem can be formulated.");
  if(not fixedParams_constructed_)
    throw std::runtime_error("Fixed parameters needs to be constructed before likelihood problem can be formulated.");

  // standard flux parameters
  GaussianPrior convNormPrior                  (priors_.convNormCenter,priors_.convNormWidth);
  GaussianPrior promptNormPrior                (priors_.promptNormCenter,priors_.promptNormWidth);
  GaussianPrior zenithCorrectionPrior          (priors_.zenithCorrectionCenter,priors_.zenithCorrectionWidth);
  GaussianPrior kaonLossesPrior                (priors_.kaonLossesCenter,priors_.kaonLossesWidth);
  GaussianPrior hadronicHEkpPrior              (priors_.hadronicHEkpCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicHEkmPrior              (priors_.hadronicHEkmCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE1pipPrior           (priors_.hadronicVHE1pipCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE1pimPrior           (priors_.hadronicVHE1pimCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3kpPrior            (priors_.hadronicVHE3kpCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3kmPrior            (priors_.hadronicVHE3kmCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3pipPrior           (priors_.hadronicVHE3pipCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3pimPrior           (priors_.hadronicVHE3pimCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3pPrior             (priors_.hadronicVHE3pCenter,std::numeric_limits<double>::max());
  GaussianPrior hadronicVHE3nPrior             (priors_.hadronicVHE3nCenter,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay1Prior                (priors_.cosmicRay1Center,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay2Prior                (priors_.cosmicRay2Center,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay3Prior                (priors_.cosmicRay3Center,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay4Prior                (priors_.cosmicRay4Center,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay5Prior                (priors_.cosmicRay5Center,std::numeric_limits<double>::max());
  GaussianPrior cosmicRay6Prior                (priors_.cosmicRay6Center,std::numeric_limits<double>::max());
  // detector systematics
  GaussianPrior icegrad0Prior                  (priors_.icegrad0Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad1Prior                  (priors_.icegrad1Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad2Prior                  (priors_.icegrad2Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad3Prior                  (priors_.icegrad3Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad4Prior                  (priors_.icegrad4Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad5Prior                  (priors_.icegrad5Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad6Prior                  (priors_.icegrad6Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad7Prior                  (priors_.icegrad7Center,std::numeric_limits<double>::max());
  GaussianPrior icegrad8Prior                  (priors_.icegrad8Center,std::numeric_limits<double>::max());
  GaussianPrior domEfficiencyPrior             (priors_.domEfficiencyCenter,priors_.domEfficiencyWidth);
  GaussianPrior holeiceForwardPrior            (priors_.holeiceForwardCenter,priors_.holeiceForwardWidth);
  // astro systematics
  GaussianPrior astroNormPrior                 (priors_.astroNormCenter,priors_.astroNormWidth);
  GaussianPrior astroDeltaGammaPrior           (priors_.astroDeltaGammaCenter,priors_.astroDeltaGammaWidth);
  GaussianPrior astroDeltaGammaSecPrior        (priors_.astroDeltaGammaSecCenter,priors_.astroDeltaGammaSecWidth);
  UniformPrior  astroPivotPrior                (priors_.astroPivotMin,priors_.astroPivotMax);
  GaussianPrior NeutrinoAntineutrinoRatioPrior (priors_.NeutrinoAntineutrinoRatioCenter,priors_.NeutrinoAntineutrinoRatioWidth);
  // xsec systematics
  GaussianPrior nuxsPrior                      (priors_.nuxsCenter,priors_.nuxsWidth);
  GaussianPrior nubarxsPrior                   (priors_.nubarxsCenter,priors_.nubarxsWidth);


  auto basicpriors=makePriorSet(
    convNormPrior,
    promptNormPrior,
    zenithCorrectionPrior,
    kaonLossesPrior,
    hadronicHEkpPrior,
    hadronicHEkmPrior,
    hadronicVHE1pipPrior,
    hadronicVHE1pimPrior,
    hadronicVHE3kpPrior,
    hadronicVHE3kmPrior,
    hadronicVHE3pipPrior,
    hadronicVHE3pimPrior,
    hadronicVHE3pPrior,
    hadronicVHE3nPrior,
    cosmicRay1Prior,
    cosmicRay2Prior,
    cosmicRay3Prior,
    cosmicRay4Prior,
    cosmicRay5Prior,
    cosmicRay6Prior,
    icegrad0Prior,
    icegrad1Prior,
    icegrad2Prior,
    icegrad3Prior,
    icegrad4Prior,
    icegrad5Prior,
    icegrad6Prior,
    icegrad7Prior,
    icegrad8Prior,
    domEfficiencyPrior,
    holeiceForwardPrior,
    astroNormPrior,
    astroDeltaGammaPrior,
    astroDeltaGammaSecPrior,
    astroPivotPrior,
    NeutrinoAntineutrinoRatioPrior,
    nuxsPrior,
    nubarxsPrior
  );


  std::vector<double> FluxCorrCenter = { 
    priors_.hadronicHEkpCenter,
    priors_.hadronicHEkmCenter,
    priors_.hadronicVHE1pipCenter,
    priors_.hadronicVHE1pimCenter,
    priors_.hadronicVHE3kpCenter,
    priors_.hadronicVHE3kmCenter,
    priors_.hadronicVHE3pipCenter,
    priors_.hadronicVHE3pimCenter,
    priors_.hadronicVHE3pCenter,
    priors_.hadronicVHE3nCenter,
    priors_.cosmicRay1Center,
    priors_.cosmicRay2Center,
    priors_.cosmicRay3Center,
    priors_.cosmicRay4Center,
    priors_.cosmicRay5Center,
    priors_.cosmicRay6Center };

  std::vector<double> FluxCorrWidth = { 
    priors_.hadronicHEkpWidth,
    priors_.hadronicHEkmWidth,
    priors_.hadronicVHE1pipWidth,
    priors_.hadronicVHE1pimWidth,
    priors_.hadronicVHE3kpWidth,
    priors_.hadronicVHE3kmWidth,
    priors_.hadronicVHE3pipWidth,
    priors_.hadronicVHE3pimWidth,
    priors_.hadronicVHE3pWidth,
    priors_.hadronicVHE3nWidth,
    priors_.cosmicRay1Width,
    priors_.cosmicRay2Width,
    priors_.cosmicRay3Width,
    priors_.cosmicRay4Width,
    priors_.cosmicRay5Width,
    priors_.cosmicRay6Width };


  std::vector<double> IceCorrCenter = { 
    priors_.icegrad0Center,
    priors_.icegrad1Center,
    priors_.icegrad2Center,
    priors_.icegrad3Center,
    priors_.icegrad4Center,
    priors_.icegrad5Center,
    priors_.icegrad6Center,
    priors_.icegrad7Center,
    priors_.icegrad8Center };

  std::vector<double> IceCorrWidth = { 
    priors_.icegrad0Width,
    priors_.icegrad1Width,
    priors_.icegrad2Width,
    priors_.icegrad3Width,
    priors_.icegrad4Width,
    priors_.icegrad5Width,
    priors_.icegrad6Width,
    priors_.icegrad7Width,
    priors_.icegrad8Width };


  std::vector<std::vector<double>> FluxCorr;
  for (unsigned int i=0; i<FluxCorrWidth.size(); i++) {    
    std::vector<double> aux_corr;
    for (unsigned int j=0; j<FluxCorrWidth.size(); j++) aux_corr.push_back(priors_.flux_corr[i][j]);
    FluxCorr.push_back(aux_corr);   
  }

  std::vector<std::vector<double>> IceCorr;
  for (unsigned int i=0; i<IceCorrWidth.size(); i++) {    
    std::vector<double> aux_corr;
    for (unsigned int j=0; j<IceCorrWidth.size(); j++) aux_corr.push_back(priors_.ice_corr[i][j]);
    IceCorr.push_back(aux_corr);   
  }

  GaussianNDPrior<16> flux_prior(FluxCorrCenter, FluxCorrWidth, FluxCorr);
  GaussianNDPrior<9> ice_prior(IceCorrCenter, IceCorrWidth, IceCorr);

  auto llhpriors = makeArbitraryPriorSet<PriorIndices>(basicpriors,flux_prior,ice_prior);

  auto fitseedvec = ConvertFitParameters(fitSeed_.front());

  prob_ = std::make_shared<LType>(phys_tools::likelihood::makeLikelihoodProblem<std::reference_wrapper<const Event>,38>(
                                  dataHist_, {simHist_}, llhpriors, {0.0}, simpleLocalDataWeighterConstructor(), DFWM,
                                  phys_tools::likelihood::SAYLikelihoodRelativeUncertaintyMod(0), fitseedvec));
  prob_->setEvaluationThreadCount(steeringParams_.evalThreads);
  prob_->likelihoodFunction.SetSigmaOverMu(steeringParams_.uncertaintyModSigmaOverMu);

  likelihood_problem_constructed_=true;
}


double GollumFit::EvalLLH(std::vector<double> nuisance, bool include_prior) const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");
  return -prob_->evaluateLikelihood(nuisance,include_prior);
}

double GollumFit::EvalLLH(FitParameters nuisance, bool include_prior) const {
  return EvalLLH(ConvertFitParameters(nuisance),include_prior);
}

phys_tools::autodiff::FD<38> GollumFit::EvalLLHGradient(std::vector<phys_tools::autodiff::FD<38>> v) const {
  return -prob_->evaluateLikelihood(v);
}

void GollumFit::ForceFitSeedSanity(){
  for(auto & fitSeed: fitSeed_)
    ForceFitSeedSanity(fitSeed);
}

void GollumFit::ForceFitSeedSanity(FitParameters& fitSeed) {
  if( fitSeed.holeiceForward < holeIceSpline_minValidValue_ or fitSeed.holeiceForward > holeIceSpline_maxValidValue_ ){
    double proposed_holeice_seed = (holeIceSpline_maxValidValue_ + holeIceSpline_minValidValue_)/2.;
    std::cout << "Holeice forward seed value ("+std::to_string(fitSeed.holeiceForward)+ ") " + " outside of spline valid range [" + std::to_string(holeIceSpline_minValidValue_) + ", " + std::to_string(holeIceSpline_maxValidValue_) + "]. Changing it to " + std::to_string(proposed_holeice_seed) + "." << std::endl;
    fitSeed.holeiceForward = proposed_holeice_seed;
  }

  if( fitSeed.domEfficiency < domEffSpline_minValidValue_ or fitSeed.domEfficiency > domEffSpline_maxValidValue_ ){
    double proposed_dom_efficiency_seed = (domEffSpline_maxValidValue_ + domEffSpline_minValidValue_)/2.;
    std::cout << "DOM efficiency seed value ("+std::to_string(fitSeed.domEfficiency)+ ") " + " outside of spline valid range [" + std::to_string(domEffSpline_minValidValue_) + ", " + std::to_string(domEffSpline_maxValidValue_) + "]. Changing it to " + std::to_string(proposed_dom_efficiency_seed) + "." << std::endl;
    fitSeed.domEfficiency = proposed_dom_efficiency_seed;
  }

}

// make a copy of this function to ML-augment
FitResult GollumFit::MinLLH() const {
  if(not likelihood_problem_constructed_)
    throw std::runtime_error("Likelihood problem has not been constructed..");

  FitResult final_result;
  final_result.likelihood = std::numeric_limits<double>::max();

  for(auto fitSeed : fitSeed_){
    std::vector<double> seed=ConvertFitParameters(fitSeed);
    prob_->setSeed(seed);

    std::vector<unsigned int> fixedIndices;
    std::vector<bool> FixVec=ConvertFitParametersFlag(fixedParams_);
    for(size_t i=0; i!=FixVec.size(); i++)
        if(FixVec[i]) fixedIndices.push_back(i);

    phys_tools::lbfgsb::LBFGSB_Driver minimizer;
    minimizer.setGradientTolerance(steeringParams_.grad_tol);
    minimizer.setChangeTolerance(steeringParams_.change_tol);

    // parameter seed, gradient factor, min boundary, max boundary
    minimizer.addParameter(  seed[0], .001, boundParams_.convNormMin,                  boundParams_.convNormMax                  ); // conv norm
    minimizer.addParameter(  seed[1], .001, boundParams_.promptNormMin,                boundParams_.promptNormMax                ); // prompt norm
    minimizer.addParameter(  seed[2], .001, boundParams_.zenithCorrectionMin,          boundParams_.zenithCorrectionMax          ); // jordi delta
    minimizer.addParameter(  seed[3], .001, boundParams_.kaonLossesMin,                boundParams_.kaonLossesMax                ); // kaon losses
    minimizer.addParameter(  seed[4], .001, boundParams_.hadronicHEkpMin,              boundParams_.hadronicHEkpMax              ); // hadronic HE K plus
    minimizer.addParameter(  seed[5], .001, boundParams_.hadronicHEkmMin,              boundParams_.hadronicHEkmMax              ); // hadronic HE K minus
    minimizer.addParameter(  seed[6], .001, boundParams_.hadronicVHE1pipMin,           boundParams_.hadronicVHE1pipMax           ); // hadronic VHE1 pi plus
    minimizer.addParameter(  seed[7], .001, boundParams_.hadronicVHE1pimMin,           boundParams_.hadronicVHE1pimMax           ); // hadronic VHE1 pi minus
    minimizer.addParameter(  seed[8], .001, boundParams_.hadronicVHE3kpMin,            boundParams_.hadronicVHE3kpMax            ); // hadronic VHE3 K plus
    minimizer.addParameter(  seed[9], .001, boundParams_.hadronicVHE3kmMin,            boundParams_.hadronicVHE3kmMax            ); // hadronic VHE3 K minus
    minimizer.addParameter( seed[10], .001, boundParams_.hadronicVHE3pipMin,           boundParams_.hadronicVHE3pipMax           ); // hadronic VHE3 pi plus
    minimizer.addParameter( seed[11], .001, boundParams_.hadronicVHE3pimMin,           boundParams_.hadronicVHE3pimMax           ); // hadronic VHE3 pi minus
    minimizer.addParameter( seed[12], .001, boundParams_.hadronicVHE3pMin,             boundParams_.hadronicVHE3pMax             ); // hadronic VHE3 proton
    minimizer.addParameter( seed[13], .001, boundParams_.hadronicVHE3nMin,             boundParams_.hadronicVHE3nMax             ); // hadronic VHE3 neutron
    minimizer.addParameter( seed[14], .001, boundParams_.cosmicRay1Min,                boundParams_.cosmicRay1Max                ); // cosmic ray PCA1
    minimizer.addParameter( seed[15], .001, boundParams_.cosmicRay2Min,                boundParams_.cosmicRay2Max                ); // cosmic ray PCA2
    minimizer.addParameter( seed[16], .001, boundParams_.cosmicRay3Min,                boundParams_.cosmicRay3Max                ); // cosmic ray PCA3
    minimizer.addParameter( seed[17], .001, boundParams_.cosmicRay4Min,                boundParams_.cosmicRay4Max                ); // cosmic ray PCA4
    minimizer.addParameter( seed[18], .001, boundParams_.cosmicRay5Min,                boundParams_.cosmicRay5Max                ); // cosmic ray PCA5
    minimizer.addParameter( seed[19], .001, boundParams_.cosmicRay6Min,                boundParams_.cosmicRay6Max                ); // cosmic ray PCA6
    minimizer.addParameter( seed[20], .001, boundParams_.icegrad0Min,                  boundParams_.icegrad0Max                  ); // ice grad 0
    minimizer.addParameter( seed[21], .001, boundParams_.icegrad1Min,                  boundParams_.icegrad1Max                  ); // ice grad 1
    minimizer.addParameter( seed[22], .001, boundParams_.icegrad2Min,                  boundParams_.icegrad2Max                  ); // ice grad 2
    minimizer.addParameter( seed[23], .001, boundParams_.icegrad3Min,                  boundParams_.icegrad3Max                  ); // ice grad 3
    minimizer.addParameter( seed[24], .001, boundParams_.icegrad4Min,                  boundParams_.icegrad4Max                  ); // ice grad 4
    minimizer.addParameter( seed[25], .001, boundParams_.icegrad5Min,                  boundParams_.icegrad5Max                  ); // ice grad 5
    minimizer.addParameter( seed[26], .001, boundParams_.icegrad6Min,                  boundParams_.icegrad6Max                  ); // ice grad 6
    minimizer.addParameter( seed[27], .001, boundParams_.icegrad7Min,                  boundParams_.icegrad7Max                  ); // ice grad 7
    minimizer.addParameter( seed[28], .001, boundParams_.icegrad8Min,                  boundParams_.icegrad8Max                  ); // ice grad 8
    minimizer.addParameter( seed[29], .001, boundParams_.domEfficiencyMin,             boundParams_.domEfficiencyMax             ); // dom eff
    minimizer.addParameter( seed[30], .001, boundParams_.holeiceForwardMin,            boundParams_.holeiceForwardMax            ); // hole ice forward
    minimizer.addParameter( seed[31], .001, boundParams_.astroNormMin,                 boundParams_.astroNormMax                 ); // astro norm
    minimizer.addParameter( seed[32], .001, boundParams_.astroDeltaGammaMin,           boundParams_.astroDeltaGammaMax           ); // astro delta gamma
    minimizer.addParameter( seed[33], .001, boundParams_.astroDeltaGammaSecMin,        boundParams_.astroDeltaGammaSecMax        ); // second astro component parameters
    minimizer.addParameter( seed[34], .001, boundParams_.astroPivotMin,                boundParams_.astroPivotMax                ); // pivot point astro component 
    minimizer.addParameter( seed[35], .001, boundParams_.NeutrinoAntineutrinoRatioMin, boundParams_.NeutrinoAntineutrinoRatioMax ); // conv particle balance
    minimizer.addParameter( seed[36], .001, boundParams_.nuxsMin,                      boundParams_.nuxsMax                      ); // nuxs
    minimizer.addParameter( seed[37], .001, boundParams_.nubarxsMin,                   boundParams_.nubarxsMax                   ); // nubarxs

    minimizer.setHistorySize(20);

    for(auto idx : fixedIndices){
      minimizer.fixParameter(idx);
    }

    FitResult result;
    result.succeeded=DoFitLBFGSB(*prob_, minimizer);
    result.likelihood=minimizer.minimumValue();
    
    // printing out LH here! 
    std::cout << "LH: " << result.likelihood << std::endl; 
    
    result.params=ConvertVecToFitParameters(minimizer.minimumPosition());

    result.nEval+=minimizer.numberOfEvaluations();
    result.nGrad+=minimizer.numberOfEvaluations();

    if(result.likelihood < final_result.likelihood)
      final_result = result;

  }

  return final_result;
}

/*************************************************************************************************************
 * Implementation of the fast mode
 * **********************************************************************************************************/


void GollumFit::ConstructFastMode(double meta_scaling) {
  if((not simulation_loaded_ ) or mainSimulation_.empty())
    throw std::runtime_error("No simulation has been loaded.");
  // This piece of code is inspired in ideas implemented in LVTools by CA. A better implementation
  // was done, independently, by BJPJ on the Sterilizer.

  std::cout << "Using metascaling " << meta_scaling << ". Note: large metascaling can produce inaccurate results." << std::endl;
  if( meta_scaling > 1. )
    throw std::runtime_error("Using metascaling greater than one can give bad results. Reconsider.");

  using MetaHistType = histogram<3,entryStoringBin<std::reference_wrapper<const Event>>>;

  MetaHistType metaHist(LogarithmicAxis(steeringParams_.logEbinEdge, steeringParams_.logEbinWidth*meta_scaling), // energy dimension
                        LinearAxis(steeringParams_.cosThbinEdge, steeringParams_.cosThbinWidth*meta_scaling), // zenith dimension
                        LinearAxis(0,1)); // particle type dimension

  auto meta_binner = [](MetaHistType& h, const Event& e){
    h.add(e.primaryEnergy,cos(e.primaryZenith),int(e.primaryType),amount(std::cref(e)));
  };

  struct combiner {
      Event operator()(const std::vector<std::reference_wrapper<const Event>>& events) {
        return combine_events(events);
      }
  };

  metaEvents_ = gollumfit::fastmode::get_fastmode_events<Event>(metaHist, meta_binner, combiner(), simHist_);

  std::cout << "Went from this many MC events " << mainSimulation_.size() << " to " << metaEvents_.size() << std::endl;
  assert(!metaEvents_.empty());

  auto newZeroHist = std::make_tuple(makeEmptyHistogramCopy(std::get<0>(simHist_)));
  bin(metaEvents_, newZeroHist, binner);
  simHist_ = std::move(newZeroHist);

  fastmode_constructed_ = true;
}

/*************************************************************************************************************
 * Functions to set options in the class
 * **********************************************************************************************************/

// Check a directory exists and throw a relevant error otherwise.
bool GollumFit::CheckDataPath(std::string p) const {
 struct stat info;
 bool status=true;
 if(p!="")
   {
     if( stat(p.c_str(), &info) != 0 )
       {
         status=false;
         throw std::runtime_error("cannot access "+ p);
       }
     else if( !(info.st_mode & S_IFDIR) )
       {
         status=false;
         throw std::runtime_error("is not a directory: " +p);
       }
   }
 else{
   std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
   return false;
 }
 return status;
}

std::string GollumFit::CheckedFilePath(std::string FilePath) const {
  std::cout<<"Checking file from path "<<FilePath<<std::endl;
  try{
    std::ifstream thefile(FilePath);
    if(thefile.good())
      return FilePath;
    else
      throw std::runtime_error("File " + FilePath + " does not exist!");
  }
  catch(std::exception &re)
    {
      throw std::runtime_error("File " + FilePath + " does not exist!");
    }
}

// Given a human readable nuisance parameter set, make a nuisance vector
std::vector<double> GollumFit::ConvertFitParameters(FitParameters ns) const {
  std::vector<double> nuis;

  nuis.push_back(ns.convNorm);
  nuis.push_back(ns.promptNorm);
  nuis.push_back(ns.zenithCorrection);
  nuis.push_back(ns.kaonLosses);
  nuis.push_back(ns.hadronicHEkp);
  nuis.push_back(ns.hadronicHEkm);
  nuis.push_back(ns.hadronicVHE1pip);
  nuis.push_back(ns.hadronicVHE1pim);
  nuis.push_back(ns.hadronicVHE3kp);
  nuis.push_back(ns.hadronicVHE3km);
  nuis.push_back(ns.hadronicVHE3pip);
  nuis.push_back(ns.hadronicVHE3pim);
  nuis.push_back(ns.hadronicVHE3p);
  nuis.push_back(ns.hadronicVHE3n);
  nuis.push_back(ns.cosmicRay1);
  nuis.push_back(ns.cosmicRay2);
  nuis.push_back(ns.cosmicRay3);
  nuis.push_back(ns.cosmicRay4);
  nuis.push_back(ns.cosmicRay5);
  nuis.push_back(ns.cosmicRay6);
  nuis.push_back(ns.icegrad0);
  nuis.push_back(ns.icegrad1);
  nuis.push_back(ns.icegrad2);
  nuis.push_back(ns.icegrad3);
  nuis.push_back(ns.icegrad4);
  nuis.push_back(ns.icegrad5);
  nuis.push_back(ns.icegrad6);
  nuis.push_back(ns.icegrad7);
  nuis.push_back(ns.icegrad8);
  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.holeiceForward);
  nuis.push_back(ns.astroNorm);
  nuis.push_back(ns.astroDeltaGamma);
  nuis.push_back(ns.astroDeltaGammaSec);
  nuis.push_back(ns.astroPivot);
  nuis.push_back(ns.NeutrinoAntineutrinoRatio);
  nuis.push_back(ns.nuxs);
  nuis.push_back(ns.nubarxs);

  assert(nuis.size() == 38);

  return nuis;
}

// Given a human readable flag set, make a bool vector
std::vector<bool> GollumFit::ConvertFitParametersFlag(FitParametersFlag ns) const {
  std::vector<bool> nuis;

  nuis.push_back(ns.convNorm);
  nuis.push_back(ns.promptNorm);
  nuis.push_back(ns.zenithCorrection);
  nuis.push_back(ns.kaonLosses);
  nuis.push_back(ns.hadronicHEkp);
  nuis.push_back(ns.hadronicHEkm);
  nuis.push_back(ns.hadronicVHE1pip);
  nuis.push_back(ns.hadronicVHE1pim);
  nuis.push_back(ns.hadronicVHE3kp);
  nuis.push_back(ns.hadronicVHE3km);
  nuis.push_back(ns.hadronicVHE3pip);
  nuis.push_back(ns.hadronicVHE3pim);
  nuis.push_back(ns.hadronicVHE3p);
  nuis.push_back(ns.hadronicVHE3n);
  nuis.push_back(ns.cosmicRay1);
  nuis.push_back(ns.cosmicRay2);
  nuis.push_back(ns.cosmicRay3);
  nuis.push_back(ns.cosmicRay4);
  nuis.push_back(ns.cosmicRay5);
  nuis.push_back(ns.cosmicRay6);
  nuis.push_back(ns.icegrad0);
  nuis.push_back(ns.icegrad1);
  nuis.push_back(ns.icegrad2);
  nuis.push_back(ns.icegrad3);
  nuis.push_back(ns.icegrad4);
  nuis.push_back(ns.icegrad5);
  nuis.push_back(ns.icegrad6);
  nuis.push_back(ns.icegrad7);
  nuis.push_back(ns.icegrad8);
  nuis.push_back(ns.domEfficiency);
  nuis.push_back(ns.holeiceForward);
  nuis.push_back(ns.astroNorm);
  nuis.push_back(ns.astroDeltaGamma);
  nuis.push_back(ns.astroDeltaGammaSec);
  nuis.push_back(ns.astroPivot);
  nuis.push_back(ns.NeutrinoAntineutrinoRatio);
  nuis.push_back(ns.nuxs);
  nuis.push_back(ns.nubarxs);

  assert(nuis.size() == 38);

  return nuis;
}

// And go back to human readable
FitParameters GollumFit::ConvertVecToFitParameters(std::vector<double> vecns) const {

  FitParameters ns;
  ns.convNorm                  = vecns[0];
  ns.promptNorm                = vecns[1];
  ns.zenithCorrection          = vecns[2];
  ns.kaonLosses                = vecns[3];
  ns.hadronicHEkp              = vecns[4];
  ns.hadronicHEkm              = vecns[5];
  ns.hadronicVHE1pip           = vecns[6];
  ns.hadronicVHE1pim           = vecns[7];
  ns.hadronicVHE3kp            = vecns[8];
  ns.hadronicVHE3km            = vecns[9];
  ns.hadronicVHE3pip           = vecns[10];
  ns.hadronicVHE3pim           = vecns[11];
  ns.hadronicVHE3p             = vecns[12];
  ns.hadronicVHE3n             = vecns[13];
  ns.cosmicRay1                = vecns[14];
  ns.cosmicRay2                = vecns[15];
  ns.cosmicRay3                = vecns[16];
  ns.cosmicRay4                = vecns[17];
  ns.cosmicRay5                = vecns[18];
  ns.cosmicRay6                = vecns[19];
  ns.icegrad0                  = vecns[20];
  ns.icegrad1                  = vecns[21];
  ns.icegrad2                  = vecns[22];
  ns.icegrad3                  = vecns[23];
  ns.icegrad4                  = vecns[24];
  ns.icegrad5                  = vecns[25];
  ns.icegrad6                  = vecns[26];
  ns.icegrad7                  = vecns[27];
  ns.icegrad8                  = vecns[28];
  ns.domEfficiency             = vecns[29];
  ns.holeiceForward            = vecns[30];
  ns.astroNorm                 = vecns[31];
  ns.astroDeltaGamma           = vecns[32];
  ns.astroDeltaGammaSec        = vecns[33];
  ns.astroPivot                = vecns[34];
  ns.NeutrinoAntineutrinoRatio = vecns[35];
  ns.nuxs                      = vecns[36];
  ns.nubarxs                   = vecns[37];

  return ns;
}

// Report back the status of object construction
void GollumFit::ReportStatus() const {
  std::cout<< "Data loaded:                   " << CheckDataLoaded() <<std::endl;
  std::cout<< "Sim loaded:                    " << CheckSimulationLoaded()  <<std::endl;
  std::cout<< "XS weighter constructed:       " << CheckCrossSectionWeighterConstructed() <<std::endl;
  std::cout<< "Flux weighter constructed:     " << CheckFluxWeighterConstructed()<<std::endl;
  std::cout<< "Lepton weighter constructed:   " << CheckLeptonWeighterConstructed()<<std::endl;
  std::cout<< "Data histogram constructed:    " << CheckDataHistogramConstructed()<<std::endl;
  std::cout<< "Sim histogram constructed:     " << CheckSimulationHistogramConstructed()<<std::endl;
  std::cout<< "LLH problem constructed:       " << CheckLikelihoodProblemConstruction()<<std::endl;
}

/*************************************************************************************************************
 * Functions to get out event distributions and set data
 * **********************************************************************************************************/

double GollumFit::SetData(marray<double,2> Data) {
  double TotalWeight=0;
  sample_.clear();
  for(size_t i=0; i!=Data.extent(0); ++i)
    {
      Event e;
      e.energy       = Data[i][0];
      e.zenith       = Data[i][1];
      e.topology     = Data[i][2];
      e.cachedWeight = Data[i][3];
      TotalWeight+=Data[i][3];
      sample_.push_back(e);
    }
  // remaking data histogram
  std::cout<<"Remaking data hist" <<std::endl;
  ConstructDataHistogram();

  return TotalWeight;
}

marray<double,2> GollumFit::GetDataEvents() const {
  marray<double,2> ReturnVec { sample_.size(), 4} ;
  for(size_t i=0; i!=sample_.size(); ++i) {
      ReturnVec[i][0]=sample_[i].energy;
      ReturnVec[i][1]=sample_[i].zenith;
      ReturnVec[i][2]=sample_[i].topology;
      ReturnVec[i][3]=sample_[i].cachedWeight;
  }
  return ReturnVec;
}


marray<double,2> GollumFit::GetRealizationEvents( FitParameters nuisance, int seed) const {
  return GetRealizationEvents(ConvertFitParameters(nuisance), seed);
}

marray<double,2> GollumFit::GetRealizationEvents( std::vector<double> nuisance, int seed) const {
  std::mt19937 rng;
  rng.seed(seed);

  auto weighter=DFWM(nuisance);
  double expected=0;
  std::vector<double> weights;
  for(const Event& e : mainSimulation_){
    auto w=weighter(e);
    weights.push_back(w);
    expected+=w;
  }

  std::vector<Event> realization=phys_tools::likelihood::generateSample(weights,mainSimulation_,expected,rng);

  marray<double,2> ReturnVec { realization.size(), 4} ;

  for(size_t i=0; i!=realization.size(); ++i)
    {
      ReturnVec[i][0]=realization[i].energy;
      ReturnVec[i][1]=realization[i].zenith;
      ReturnVec[i][2]=realization[i].topology;
      ReturnVec[i][3]=1.;
    }
  return ReturnVec;
}

marray<double,2> GollumFit::GetExpectationEvents(FitParameters fit_params) const {
  return GetExpectationEvents(ConvertFitParameters(fit_params));
}

marray<double,2> GollumFit::GetExpectationEvents(std::vector<double> fit_params) const {
  const std::deque<Event> & sim = (fastmode_constructed_)?metaEvents_:mainSimulation_;
  marray<double,2> ReturnVec {sim.size(), 4} ;
  auto weighter = DFWM(fit_params);
  for(size_t i=0; i!=sim.size(); ++i) {
      ReturnVec[i][0]=sim[i].energy;
      ReturnVec[i][1]=sim[i].zenith;
      ReturnVec[i][2]=sim[i].topology;
      ReturnVec[i][3]=weighter(sim[i]);
  }
  return ReturnVec;
}

int GollumFit::CheckExpectation(FitParameters fit_params) const {
  return CheckExpectation(ConvertFitParameters(fit_params));
}

int GollumFit::CheckExpectation(std::vector<double> fit_params) const {
  const std::deque<Event> & sim = (fastmode_constructed_)?metaEvents_:mainSimulation_;
  auto weighter = DFWM(fit_params);
  for(size_t i=0; i!=sim.size(); ++i) {
      if ( weighter(sim[i])<0 ) return 1;
  }
  return 0;
}

/*************************************************************************************************************
 * Functions to get bin edges
 * **********************************************************************************************************/

// Given a histogram reference h, get bin edges in dimension dim
  template<unsigned int hist_index>
  std::vector<double> GollumFit::PullBinEdges(int dim, const HistType& hh) const{
    auto h = std::get<hist_index>(hh);
    std::vector<double> edges_i(h.getBinCount(dim));
    for(unsigned int j=0; j<h.getBinCount(dim); j++)
      edges_i[j]=h.getBinEdge(dim,j);
    edges_i.push_back(h.getBinEdge(dim,h.getBinCount(dim)-1)+h.getBinWidth(dim,h.getBinCount(dim)-1));
    return edges_i;
  }


  std::vector<double> GollumFit::GetEnergyBinsMC() const{
    return PullBinEdges<0>(0,simHist_);
  }

  std::vector<double> GollumFit::GetZenithBinsMC() const{
    return PullBinEdges<0>(1,simHist_);
  }

  std::vector<double> GollumFit::GetTopologyBinsMC() const{
    return PullBinEdges<0>(2,simHist_);
  }


} // close namespace gollumfit
