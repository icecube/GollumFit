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
* @brief Functions to get MC simulation sets, with runtime registration support.
*/

namespace gollumfit {
  namespace sterile {

    /**
    * @brief Singleton registry for MC simulation sets.
    *
    * Holds both individual MCSet entries and compound simulation tags.
    * Hardcoded defaults are loaded lazily via InitDefaults(). Custom sets
    * can be registered from C++ or Python before constructing GollumFit.
    */
    class MCSetRegistry {
    public:
      /// Get the singleton instance.
      inline static MCSetRegistry& Instance() {
        static MCSetRegistry instance;
        return instance;
      }

      /// Register a single MCSet with a unique name.
      void RegisterMCSet(const std::string& name, MCSet mcset) {
        simInfo_.insert_or_assign(name, std::move(mcset));
      }

      /// Register a compound simulation tag that resolves to multiple MCSet names.
      void RegisterSimulationTag(const std::string& tag, const std::vector<std::string>& mcset_names) {
        compoundTags_[tag] = mcset_names;
      }

      /// Get the full simulation info map.
      const std::map<std::string, MCSet>& GetSimInfo() const {
        return simInfo_;
      }

      /// Check if defaults have been initialized.
      bool DefaultsInitialized() const { return defaults_initialized_; }

      /// Clear all registered sets (both defaults and custom).
      void ClearAll() {
        simInfo_.clear();
        compoundTags_.clear();
        defaults_initialized_ = false;
      }

      /**
      * @brief Initialize the hardcoded default MC sets.
      *
      * @param mc_dataPath Path to the STERILE MC directory.
      */
      void InitDefaults(const std::string& mc_dataPath) {
        if (defaults_initialized_) return;

        // Individual MC sets
        RegisterMCSet("Cheetah_Test", MCSet(mc_dataPath,"Cheetah_Test_1.h5",{false,0},1,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Cheetah_Generation_data.lic",true))));

        RegisterMCSet("Cheetah_Mini", MCSet(mc_dataPath,"Cheetah_2997.h5",{false,0},2997,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Cheetah_Generation_data.lic",true))));

        RegisterMCSet("BDT_Inelasticity", MCSet(mc_dataPath,"BDT_Inelasticity_992.h5",{false,0},992,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true))));

        RegisterMCSet("BDT_Tau_Inelasticity", MCSet(mc_dataPath,"BDT_Tau_Inelasticity_1000.h5",{false,0},1000,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Tau.lic",true))));

        RegisterMCSet("BDT_Test_HE", MCSet(mc_dataPath,"BDT_19738",{true,1},992,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true))));

        RegisterMCSet("BDT_Split_HE", MCSet(mc_dataPath,"BDT_19738",{true,5},19738,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true))));

        RegisterMCSet("BDT_Tau", MCSet(mc_dataPath,"BDT_Tau.h5",{false,0},5000,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Tau.lic",true))));

        RegisterMCSet("Platinum_Split_HE", MCSet(mc_dataPath,"Platinum_98000",{true,98},98000/5.,1.27,-1,
                     LW::MakeGeneratorsFromLICFile(tools::CheckedFilePath(mc_dataPath + "/" + "Platinum_Generation_data.lic",true))));

        // Compound tags
        RegisterSimulationTag("Platinum_Split_Plus_Tau", {"Platinum_Split_HE", "Platinum_Tau"});
        RegisterSimulationTag("BDT_Inelasticity_Plus_Tau", {"BDT_Inelasticity", "BDT_Tau_Inelasticity"});
        RegisterSimulationTag("BDT_Split_HE_Plus_Tau", {"BDT_Split_HE", "BDT_Tau"});
        RegisterSimulationTag("BDT_Test_HE_Plus_Tau", {"BDT_Test_HE", "BDT_Tau"});
        RegisterSimulationTag("BDT1_Split_HE_Plus_Tau", {"BDT1_Split_HE", "BDT1_Tau"});
        RegisterSimulationTag("BDT1_Test_HE_Plus_Tau", {"BDT1_Test_HE", "BDT1_Tau"});

        defaults_initialized_ = true;
      }

      /**
      * @brief Retrieve simulation sets for a given tag.
      *
      * Checks compound tags first, then direct matches, then _LE/_HE/_EHE suffixes.
      *
      * @param simulation_tag Tag identifying the desired simulation set(s).
      * @return std::vector<MCSet> The requested MCSet objects.
      * @throw std::runtime_error If the tag cannot be resolved.
      */
      std::vector<MCSet> GetSimulationSets(const std::string& simulation_tag) const {
        // Check compound tags
        auto ct = compoundTags_.find(simulation_tag);
        if (ct != compoundTags_.end()) {
          std::vector<MCSet> result;
          for (const auto& name : ct->second) {
            auto it = simInfo_.find(name);
            if (it == simInfo_.end())
              throw std::runtime_error("MCSetRegistry: compound tag '" + simulation_tag +
                                       "' references unknown MCSet '" + name + "'");
            result.push_back(it->second);
          }
          return result;
        }

        // Direct match
        if (simInfo_.count(simulation_tag))
          return {simInfo_.at(simulation_tag)};

        // Try _LE/_HE/_EHE suffix combinations
        bool has_LE  = simInfo_.count(simulation_tag+"_LE")  > 0;
        bool has_HE  = simInfo_.count(simulation_tag+"_HE")  > 0;
        bool has_EHE = simInfo_.count(simulation_tag+"_EHE") > 0;

        if (has_LE && has_HE && has_EHE)
          return {simInfo_.at(simulation_tag+"_LE"), simInfo_.at(simulation_tag+"_HE"), simInfo_.at(simulation_tag+"_EHE")};
        else if (has_LE && has_HE)
          return {simInfo_.at(simulation_tag+"_LE"), simInfo_.at(simulation_tag+"_HE")};
        else if (has_HE && has_EHE)
          return {simInfo_.at(simulation_tag+"_HE"), simInfo_.at(simulation_tag+"_EHE")};
        else if (has_LE && has_EHE)
          return {simInfo_.at(simulation_tag+"_LE"), simInfo_.at(simulation_tag+"_EHE")};
        else if (has_LE)
          return {simInfo_.at(simulation_tag+"_LE")};
        else if (has_HE)
          return {simInfo_.at(simulation_tag+"_HE")};
        else if (has_EHE)
          return {simInfo_.at(simulation_tag+"_EHE")};
        else
          throw std::runtime_error("GollumFit::Invalid simulation tag: " + simulation_tag +
                                   ". Do not append EHE, HE or LE postfix to select them both.");
      }

    private:
      MCSetRegistry() = default;
      std::map<std::string, MCSet> simInfo_;
      std::map<std::string, std::vector<std::string>> compoundTags_;
      bool defaults_initialized_ = false;
    };


    // ---- Backward-compatible free functions ----

    /**
    * @brief Constructs a map of simulation identifiers to MCSet objects.
    * Delegates to MCSetRegistry singleton.
    */
    inline std::map<std::string,MCSet> GetSimInfo(std::string mc_dataPath)
    {
      auto& reg = MCSetRegistry::Instance();
      if (!reg.DefaultsInitialized()) reg.InitDefaults(mc_dataPath);
      return std::map<std::string,MCSet>(reg.GetSimInfo());
    }

    /**
    * @brief Retrieves a vector of MCSet objects based on a simulation tag.
    * Delegates to MCSetRegistry singleton.
    */
    inline std::vector<MCSet> GetSimulationSets(std::string simulation_tag, std::string mc_dataPath){
      auto& reg = MCSetRegistry::Instance();
      if (!reg.DefaultsInitialized()) reg.InitDefaults(mc_dataPath);
      return reg.GetSimulationSets(simulation_tag);
    }

  } // close sterile namespace

} // close gollumfit namespace

#endif /* GOLLUMMCSPECIFICATIONS_H_ */
