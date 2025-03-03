#include <vector>
#include <iostream>
#include <string>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>

#include <GollumFit/utils.h>
#include <GollumFit/json.hpp>
#include <GollumFit/GollumTools.h>

#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>

#include "tuple_hash.h"

using namespace nusquids;
using namespace nlohmann;

std::vector<double> ConvertToVector(const marray<double,1>& array){
  std::vector<double> vec(array.size());
  for(unsigned int i = 0; i < array.size(); i++){
    vec[i] = array[i];
  }
  return vec;
}

int main(int argc, char* argv[]){
  using splinetable=photospline::splinetable<>;

  // get mixing angles and square mass difference
  double dm41sq,th14,th24,th34,del14,del24;
  unsigned int index;
  std::string input_earth_path;
  std::string input_xsec_dir;
  std::string output_path;
  if(argc != 11){
      printf("ERROR:USAGE: program inearth inxsec outdir index dmsq41 [eVˆ2] th14 th24 th34 del14 del24 [rad] \n");
      exit(0);
  } else {
      input_earth_path = argv[1];
      input_xsec_dir   = argv[2];
      output_path      = argv[3];
      index            = atoi(argv[4]);
      dm41sq           = atof(argv[5]);
      th14             = atof(argv[6]);
      th24             = atof(argv[7]);
      th34             = atof(argv[8]);
      del14            = atof(argv[9]);
      del24            = atof(argv[10]);
  }
  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Indir Xsec: "<<input_xsec_dir<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  bool doStack = true;
  const unsigned int numneu = 3;
  const unsigned int number_of_scales = 25;
  const squids::Const units;
  nuSQUIDSAtm<> nus_atm(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350),numneu,both,true);

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path);
  nus_atm.Set_EarthModel(earth);

  nus_atm.Set_ProgressBar(false);
  nus_atm.Set_TauRegeneration(true);
  nus_atm.Set_EvalThreads(15);

  // Values from https://arxiv.org/pdf/1611.01514.pdf
  // set mixing angles and masses
  nus_atm.Set_MixingAngle(0,1,0.585732);
  nus_atm.Set_MixingAngle(0,2,0.147655);
  nus_atm.Set_MixingAngle(1,2,0.726057);
  nus_atm.Set_MixingAngle(0,3,th14);
  nus_atm.Set_MixingAngle(1,3,th24);
  nus_atm.Set_MixingAngle(2,3,th34);

  nus_atm.Set_SquareMassDifference(1,7.50e-05);
  nus_atm.Set_SquareMassDifference(2,0.00252);
  nus_atm.Set_SquareMassDifference(3,dm41sq);

  nus_atm.Set_CPPhase(0,2,0.0);
  nus_atm.Set_CPPhase(0,3,del14);
  nus_atm.Set_CPPhase(1,3,del24);

  double error = 1.0e-17;
  // setup integration settings
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm.Set_rel_error(error);
  nus_atm.Set_abs_error(error);

  double baseline_astro_normalization = 1.0e-18; // nu/GeV/s/cm^2/sr
  double baseline_astro_spectral_index = -2.5; // center of things
  auto fluxAstro_ = std::make_shared<LW::PowerLawFlux>(baseline_astro_normalization,
                                                       baseline_astro_spectral_index);

  // construct the kaon initial state
  marray<double,4> inistate {nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
  std::fill(inistate.begin(),inistate.end(),0);

  marray<double,1> cos_range = nus_atm.GetCosthRange();
  marray<double,1> e_range = nus_atm.GetERange();

  LW::Event scratch_lw_e;
   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       double cth = cos_range[ci];

       scratch_lw_e.energy=enu;
       scratch_lw_e.zenith=acos(cth);

       inistate[ci][ei][0][0] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][0][1] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][0][2] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][0][3] = 0.;

       inistate[ci][ei][1][0] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][1][1] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][1][2] = (*fluxAstro_)(scratch_lw_e);
       inistate[ci][ei][1][3] = 0.;
     }
   }

  nus_atm.Set_initial_state(inistate,flavor);

  nus_atm.EvolveState();

  // spline energy grid
  auto costh_grid = linspace(-1.1,0.2,300);
  auto e_grid = logspace(1.e2,1.e6,600);

  // photo spline settings
  std::ifstream settings_stream("./spline_settings.txt");
  json settings_json;
  settings_stream >> settings_json;

  // Define spline properties
  std::vector<double> smoothing = settings_json.at("smoothing");
  std::vector<unsigned int> splineOrder = settings_json.at("splineOrder");

  std::cout << json(smoothing) << std::endl;
  std::cout << json(splineOrder) << std::endl;

  std::vector<double> energyKnots = settings_json.at("energyKnots");
  std::vector<double> zenithKnots = settings_json.at("zenithKnots");

  std::transform(energyKnots.begin(),energyKnots.end(),energyKnots.begin(),(double(*)(double))&log10);
  std::sort(energyKnots.begin(), energyKnots.end());
  std::sort(zenithKnots.begin(), zenithKnots.end());

  std::cout << "Energy knots: ";
  for(auto e : energyKnots) {
      std::cout << e << " ";
  }
  std::cout << std::endl;

  std::cout << "zenith knots: ";
  for(auto e : zenithKnots) {
      std::cout << e << " ";
  }
  std::cout << std::endl;

  // nusquids magic
  std::map<std::string, unsigned int> neutrino_label_to_nusquids {{"neutrino",0}, {"antineutrino",1}};

  // scaling grid
  auto scalings = linspace(0.8,1.2,number_of_scales);

  // splines map
  std::map<std::tuple<double,std::string>,std::shared_ptr<photospline::splinetable<>>> splines;

  for(double scale : scalings){
    std::cout << "Calculating propagation for scaling " << scale << std::endl;

    std::vector<double> edges {1.e2*units.GeV,1.e8*units.GeV};
    std::vector<double> scale_segment {scale};
    auto scaled_xs = std::make_shared<CrossSectionLibrary>() ;
    scaled_xs->addTarget(proton, gollumfit::tools::ScaledNeutrinoCrossSections(input_xsec_dir+"/csms_proton.h5",edges,scale_segment));
    scaled_xs->addTarget(neutron,gollumfit::tools::ScaledNeutrinoCrossSections(input_xsec_dir+"/csms_neutron.h5",edges,scale_segment));
    scaled_xs->addTarget(electron,GlashowResonanceCrossSection());
    nuSQUIDSAtm<> nus_atm_scaled(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350),numneu,both,true,scaled_xs);
    nus_atm_scaled.Set_EarthModel(earth);

    nus_atm_scaled.Set_ProgressBar(false);
    nus_atm_scaled.Set_TauRegeneration(true);
    nus_atm_scaled.Set_EvalThreads(15);

    // Values from https://arxiv.org/pdf/1611.01514.pdf
    // set mixing angles and masses
    nus_atm_scaled.Set_MixingAngle(0,1,0.585732);
    nus_atm_scaled.Set_MixingAngle(0,2,0.147655);
    nus_atm_scaled.Set_MixingAngle(1,2,0.726057);
    nus_atm_scaled.Set_MixingAngle(0,3,th14);
    nus_atm_scaled.Set_MixingAngle(1,3,th24);
    nus_atm_scaled.Set_MixingAngle(2,3,th34);

    nus_atm_scaled.Set_SquareMassDifference(1,7.50e-05);
    nus_atm_scaled.Set_SquareMassDifference(2,0.00252);
    nus_atm_scaled.Set_SquareMassDifference(3,dm41sq);

    nus_atm_scaled.Set_CPPhase(0,2,0.0);
    nus_atm_scaled.Set_CPPhase(0,3,del14);
    nus_atm_scaled.Set_CPPhase(1,3,del24);

    // setup integration settings
    nus_atm_scaled.Set_GSL_step(gsl_odeiv2_step_rk4);
    nus_atm_scaled.Set_rel_error(error);
    nus_atm_scaled.Set_abs_error(error);

    nus_atm_scaled.Set_initial_state(inistate,flavor);

    nus_atm_scaled.EvolveState();

    for(std::string primary_type : {"neutrino","antineutrino"}){
      std::vector<std::tuple<double, double, size_t, size_t>> sparseHistEntries;
      for(size_t ic=0; ic<costh_grid.extent(0); ic++){ // zenith
          for(size_t ie=0; ie<e_grid.extent(0); ie++){ // energy
            double costh = costh_grid[ic];
            if ( costh < -1.0 )
              costh = -1.0;
            double exp = nus_atm_scaled.EvalFlavor(index,costh,e_grid[ie]*units.GeV,neutrino_label_to_nusquids[primary_type]);
            exp /= nus_atm.EvalFlavor(index,costh,e_grid[ie]*units.GeV,neutrino_label_to_nusquids[primary_type]);
            double weight = 1.;// relative error
            if(exp> 0) {
                sparseHistEntries.push_back(std::make_tuple((exp), weight, ie, ic));
            }
          }
      }
      auto sparseHist = std::make_shared<photospline::ndsparse>(sparseHistEntries.size(), 2);
      std::vector<double> weights;
      unsigned int coord[2];
      for(auto entry : sparseHistEntries) {
          coord[0] = std::get<2>(entry);
          coord[1] = std::get<3>(entry);
          sparseHist->insertEntry(std::get<0>(entry), coord);
          weights.push_back(std::get<1>(entry));
      }

      auto spline = std::make_shared<photospline::splinetable<>>();

      std::vector<std::vector<double>> binCoordinates(2);
      binCoordinates[0] = ConvertToVector(e_grid);
      std::transform(binCoordinates[0].begin(), binCoordinates[0].end(), binCoordinates[0].begin(), [](double e)->double{return log10(e);});
      binCoordinates[1] = ConvertToVector(costh_grid);

      spline->fit(*sparseHist, weights, binCoordinates, splineOrder, {energyKnots, zenithKnots}, smoothing, {2});

      splines[std::make_tuple(scale, primary_type)] = spline;
    }
  }

  // Save splines
  std::cout << "Saving splines..." << std::endl;
  std::vector<photospline::splinetable<> *> splinesToSave;

  for(std::string primary_type : {"neutrino","antineutrino"}){
    std::vector<photospline::splinetable<> *> splinesToStack;
    for(double scale : scalings){
        auto spline = splines[std::make_tuple(scale, primary_type)];
        splinesToStack.push_back(&(*spline));
    }

    if(doStack){
        std::cout << "Stacking splines..." << std::endl;
        photospline::splinetable<> spline(splinesToStack, ConvertToVector(scalings), 2);
        spline.write_fits(output_path+"/attenuation_spline_diffuseAstro_" + argv[4] + "_" + primary_type + ".fits");
    }
  }

  return 0;
}
