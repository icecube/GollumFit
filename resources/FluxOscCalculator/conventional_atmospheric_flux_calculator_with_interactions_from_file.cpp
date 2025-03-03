#include <vector>
#include <iostream>
#include <string>
#include <nuSQuIDS/nuSQuIDS.h>
#include <nuSQuIDS/marray.h>

using namespace nusquids;

int main(int argc, char* argv[]){
  // get mixing angles and square mass difference
  double dm41sq,th14,th24,th34,del14,del24;
  unsigned int index;
  std::string input_flux_path, input_earth_path;
  std::string output_path;
  if(argc != 11){
      printf("ERROR:USAGE: program influx inearth outdir index dmsq41 [eVˆ2] th14 th24 th34 del14 del24 [rad] \n");
      exit(0);
  } else {
      input_flux_path  = argv[1];
      input_earth_path = argv[2];
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

  std::cout<<"Inpath Flux: "<<input_flux_path<<std::endl;
  std::cout<<"Inpath Earth: "<<input_earth_path<<std::endl;
  std::cout<<"Outpath: "<<output_path<<std::endl;

  const unsigned int numneu = 4;
    
  const squids::Const units;
  nuSQUIDSAtm<> nus_atm(linspace(-1.,0.2,100),logspace(1.e2*units.GeV,1.e6*units.GeV,350),numneu,both,true);

  std::shared_ptr<EarthAtm> earth = std::make_shared<EarthAtm>(input_earth_path);
  nus_atm.Set_EarthModel(earth);

  nus_atm.Set_TauRegeneration(true);

  nus_atm.Set_ProgressBar(false);

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

  double error = 1.0e-15;
  // setup integration settings
  nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4);
  nus_atm.Set_rel_error(error);
  nus_atm.Set_abs_error(error);


  // loading kaon and pion flux files
  marray<double,2> input_flux = quickread(input_flux_path);

  // construct the kaon initial state
   marray<double,4> inistate {nus_atm.GetNumCos(),nus_atm.GetNumE(),2,numneu};
   std::fill(inistate.begin(),inistate.end(),0);

   marray<double,1> cos_range = nus_atm.GetCosthRange();
   marray<double,1> e_range = nus_atm.GetERange();

   assert( input_flux.extent(0) == nus_atm.GetNumCos()*nus_atm.GetNumE() );

   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei]/units.GeV;
       assert( std::fabs(enu - input_flux[ci*e_range.size() + ei][1] ) < 1.e-4 );
       double cth = cos_range[ci];
       assert( std::fabs(cth- input_flux[ci*e_range.size() + ei][0] ) < 1.e-4 );

       inistate[ci][ei][0][0] = input_flux[ci*e_range.size() + ei][2];
       inistate[ci][ei][0][1] = input_flux[ci*e_range.size() + ei][4];
       inistate[ci][ei][0][2] = 0.;
       inistate[ci][ei][0][3] = 0.;

       inistate[ci][ei][1][0] = input_flux[ci*e_range.size() + ei][3];
       inistate[ci][ei][1][1] = input_flux[ci*e_range.size() + ei][5];
       inistate[ci][ei][1][2] = 0.;
       inistate[ci][ei][1][3] = 0.;
     }
   }

  nus_atm.Set_initial_state(inistate,flavor);

  nus_atm.EvolveState();

   for ( int ci = 0 ; ci < nus_atm.GetNumCos(); ci++){
     for ( int ei = 0 ; ei < nus_atm.GetNumE(); ei++){
       double enu = e_range[ei];
       double cth = cos_range[ci];
       for(unsigned int flv = 0; flv < 3; flv++){
         if(nus_atm.EvalFlavor(flv,cth,enu,0) < 0)
           std::cout << "neg nu    propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,0) << std::endl;
         if(nus_atm.EvalFlavor(flv,cth,enu,1) < 0)
           std::cout << "neg nubar propagated fluxes: " << flv << " " << cth << " " << enu/units.GeV << " " << nus_atm.EvalFlavor(flv,cth,enu,1) << std::endl;
       }
     }
   }

  nus_atm.WriteStateHDF5(output_path+"/atmospheric_"+std::to_string(index)+"_"+std::to_string(dm41sq)+"_"+std::to_string(th14)+"_"+std::to_string(th24)+"_"+std::to_string(th34)+"_"+std::to_string(del14)+"_"+std::to_string(del24)+".hdf5");

  std::cout << "finish" << std::endl;

  return 0;
}
