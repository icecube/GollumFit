#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>
#include <LeptonWeighter/Flux.h>
#include <LeptonWeighter/nuSQFluxInterface.h>

using namespace nusquids;

int main(int argc, char* argv[]){
  // get mixing angles and square mass difference
  double dm41sq,th14,th24,th34,del14,del24;
  unsigned int index;
  std::string input_earth_path;
  std::string output_path;
  if(argc != 10){
      printf("ERROR:USAGE: program inearth outdir index dmsq41 [eVˆ2] th14 th24 th34 del14 del24 [rad] \n");
      exit(0);
  } else {
      input_earth_path = argv[1];
      output_path      = argv[2];
      index            = atoi(argv[3]);
      dm41sq           = atof(argv[4]);
      th14             = atof(argv[5]);
      th24             = atof(argv[6]);
      th34             = atof(argv[7]);
      del14            = atof(argv[8]);
      del24            = atof(argv[9]);
  }
  if(output_path[output_path.length()-1]!='/'){
    output_path = output_path +'/';}

  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  double baseline_astro_normalization = 1.0e-18; // nu/GeV/s/cm^2/sr
  double baseline_astro_spectral_index = -2.5; // center of things
  auto fluxAstro_ = std::make_shared<LW::PowerLawFlux>(baseline_astro_normalization,
                                                       baseline_astro_spectral_index);

  const unsigned int numneu = 4;
  const squids::Const units;
  nuSQUIDSAtm<> nus_atm_astro(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350),numneu,both,true);

  nus_atm_astro.Set_TauRegeneration(true);

  nus_atm_astro.Set_IncludeOscillations(true);
  nus_atm_astro.Set_ProgressBar(false);

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path);
  nus_atm_astro.Set_EarthModel(earth);

  // Values from https://arxiv.org/pdf/1611.01514.pdf
  // set mixing angles and masses
  nus_atm_astro.Set_MixingAngle(0,1,0.585732);
  nus_atm_astro.Set_MixingAngle(0,2,0.147655);
  nus_atm_astro.Set_MixingAngle(1,2,0.726057);
  nus_atm_astro.Set_MixingAngle(0,3,th14);
  nus_atm_astro.Set_MixingAngle(1,3,th24);
  nus_atm_astro.Set_MixingAngle(2,3,th34);

  nus_atm_astro.Set_SquareMassDifference(1,7.50e-05);
  nus_atm_astro.Set_SquareMassDifference(2,0.00252);
  nus_atm_astro.Set_SquareMassDifference(3,dm41sq);

  nus_atm_astro.Set_CPPhase(0,2,0.0);
  nus_atm_astro.Set_CPPhase(0,3,del14);
  nus_atm_astro.Set_CPPhase(1,3,del24);

  double error = 1.0e-17;
  // setup integration settings
  nus_atm_astro.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm_astro.Set_rel_error(error);
  nus_atm_astro.Set_abs_error(error);


  // construct the astro initial state
   marray<double,4> inistate_astro {nus_atm_astro.GetNumCos(),nus_atm_astro.GetNumE(),2,numneu};
   std::fill(inistate_astro.begin(),inistate_astro.end(),0);

   marray<double,1> cos_range = nus_atm_astro.GetCosthRange();
   marray<double,1> e_range = nus_atm_astro.GetERange();

  LW::Event scratch_lw_e;
   for ( int ci = 0 ; ci < nus_atm_astro.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm_astro.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       double cth = cos_range[ci];

       scratch_lw_e.energy=enu;
       scratch_lw_e.zenith=acos(cth);

       inistate_astro[ci][ei][0][0] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][1] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][2] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][0][3] = 0.;

       inistate_astro[ci][ei][1][0] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][1] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][2] = (*fluxAstro_)(scratch_lw_e);
       inistate_astro[ci][ei][1][3] = 0.;
     }
   }

  nus_atm_astro.Set_initial_state(inistate_astro,flavor);

  nus_atm_astro.EvolveState();


   for ( int ci = 0 ; ci < nus_atm_astro.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm_astro.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < 3; flv++){
         if(nus_atm_astro.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_astro.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm_astro.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm_astro.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }


  nus_atm_astro.WriteStateHDF5(output_path+"/astro_"+std::to_string(index)+"_"+std::to_string(dm41sq)+"_"+std::to_string(th14)+"_"+std::to_string(th24)+"_"+std::to_string(th34)+"_"+std::to_string(del14)+"_"+std::to_string(del24)+".hdf5");

  std::cout << "finish" << std::endl;

  return 0;
}
