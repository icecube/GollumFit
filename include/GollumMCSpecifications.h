#ifndef GOLLUMMCSPECIFICATIONS_H_
#define GOLLUMMCSPECIFICATIONS_H_

#include <LeptonWeighter/Generator.h>
#include <LeptonWeighter/Event.h>
#include <LeptonWeighter/ParticleType.h>
#include "GollumTools.h"
#include "GollumEnumDefinitions.h"
#include "GollumMCSet.h"

/**
* @file GollumMCSpecifications.h
* @brief Functions to get MC simulation sets
*/

namespace gollumfit {
  namespace sterile {

    /**
    * @brief Constructs a map of simulation identifiers to MCSet objects.
    * 
    * This function reads simulation configuration and data from specified files 
    * and constructs MCSet objects with corresponding parameters. The function
    * expects the presence of .lic and .h5 files in the provided mc_dataPath
    * directory. These files are used to initialize the MCSet objects which are
    * then mapped to their respective string identifiers.
    * 
    * @param mc_dataPath The path to the directory containing the .lic and .h5 
    * files necessary for initializing MCSet objects.
    * @return std::map<std::string, MCSet> A map linking simulation set 
    * identifiers to their corresponding MCSet objects.
    * @throw std::runtime_error Throws if the file paths are invalid or required 
    * files are missing.
    */
    std::map<std::string,MCSet> GetSimInfo(std::string mc_dataPath)
    {
       std::map<std::string,MCSet> simInfo = {

        /*
        * OFFICIAL ARES SETS
        */
        {"BDT_Test_HE",MCSet(mc_dataPath,"BDT_19738",{true,1},992,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true)))
        },
        {"BDT_Split_HE",MCSet(mc_dataPath,"BDT_19738",{true,5},19738,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true)))
        },
        {"BDT_Tau",MCSet(mc_dataPath,"BDT_Tau.h5",{false,0},5000,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Tau.lic",true)))
        },
        {"Platinum_Split_HE",MCSet(mc_dataPath,"Platinum_98000",{true,98},98000/5.,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true)))
        },


      };
      return simInfo;
    }


    /**
    * @brief Retrieves a vector of MCSet objects based on a simulation tag.
    * 
    * This function allows users to retrieve various combinations of simulation 
    * sets by specifying a simulation tag. The function searches for the 
    * simulation sets within a map obtained from GetSimInfo(mc_dataPath). It can
    * handle specific combinations of simulation sets (e.g., combining low 
    * energy and high energy sets) or retrieve a single simulation set based on 
    * the provided tag.
    * 
    * @param simulation_tag A string identifier for the desired simulation set(s).
    * @param mc_dataPath The path to the directory containing the .lic and .h5 
    * files necessary for initializing MCSet objects.
    * @return std::vector<MCSet> A vector containing the requested MCSet objects.
    * @throw std::runtime_error Throws if the simulation tag does not match any 
    * predefined set or combination of sets, or if appending "_LE", "_HE", or 
    * "_EHE" does not result in valid sets.
    */
    std::vector<MCSet> GetSimulationSets(std::string simulation_tag, std::string mc_dataPath){
      auto simulations = GetSimInfo(mc_dataPath);

      if(simulation_tag=="Platinum_Split_Plus_Tau")
          return {simulations.at("Platinum_Split_HE"), simulations.at("Platinum_Tau")};

      if(simulation_tag=="BDT_Split_HE_Plus_Tau")
          return {simulations.at("BDT_Split_HE"), simulations.at("BDT_Tau")};

      if(simulation_tag=="BDT_Test_HE_Plus_Tau")
          return {simulations.at("BDT_Test_HE"), simulations.at("BDT_Tau")};

      if(simulation_tag=="BDT1_Split_HE_Plus_Tau")
          return {simulations.at("BDT1_Split_HE"), simulations.at("BDT1_Tau")};

      if(simulation_tag=="BDT1_Test_HE_Plus_Tau")
          return {simulations.at("BDT1_Test_HE"), simulations.at("BDT1_Tau")};

      // user specifies directly the simulation they want
      if(not (simulations.find(simulation_tag) == simulations.end()))
        return {simulations.at(simulation_tag)};

      // else we use LE and HE sets
      if((not (simulations.find(simulation_tag+"_LE") == simulations.end())) and (not (simulations.find(simulation_tag+"_HE") == simulations.end())) and (not (simulations.find(simulation_tag+"_EHE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE"),simulations.at(simulation_tag+"_HE"),simulations.at(simulation_tag+"_EHE")};
      else if ((not (simulations.find(simulation_tag+"_LE") == simulations.end())) and (not (simulations.find(simulation_tag+"_HE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE"),simulations.at(simulation_tag+"_HE")};
      else if ((not (simulations.find(simulation_tag+"_HE") == simulations.end())) and (not (simulations.find(simulation_tag+"_EHE") == simulations.end())))
        return {simulations.at(simulation_tag+"_HE"),simulations.at(simulation_tag+"_EHE")};
      else if ((not (simulations.find(simulation_tag+"_LE") == simulations.end())) and (not (simulations.find(simulation_tag+"_EHE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE"),simulations.at(simulation_tag+"_EHE")};
      else if ((not (simulations.find(simulation_tag+"_LE") == simulations.end())))
        return {simulations.at(simulation_tag+"_LE")};
      else if ((not (simulations.find(simulation_tag+"_HE") == simulations.end())))
        return {simulations.at(simulation_tag+"_HE")};
      else if ((not (simulations.find(simulation_tag+"_EHE") == simulations.end())))
        return {simulations.at(simulation_tag+"_EHE")};
      else
        throw std::runtime_error("GollumFit::Invalid simulation tag: "+ simulation_tag + ". Do not append EHE, HE or LE postfix to select them both.");
    }

  } // close sterile namespace

} // close gollumfit namespace

#endif /* GOLLUMMCSPECIFICATIONS_H_ */
