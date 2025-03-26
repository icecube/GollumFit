#include "GollumTools.h"

/**
* @file GollumTools.cpp
* @brief Helper function definitions.
*/

namespace gollumfit {
  namespace tools {

    bool isNeutrino(LW::ParticleType p){
      return(p==LW::ParticleType::NuE 
             || p==LW::ParticleType::NuMu
             || p==LW::ParticleType::NuTau);
    }

    bool isAntineutrino(LW::ParticleType p){
      return(p==LW::ParticleType::NuEBar
             || p==LW::ParticleType::NuMuBar
             || p==LW::ParticleType::NuTauBar);
    }

    std::string CheckedFilePath(std::string FilePath, bool quiet) {
      if(!quiet) std::cout<<"Reading a file from path "<<FilePath<<std::endl;
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

  }// closing gollumfit::tools namespace
}// closing gollumfit namespace

