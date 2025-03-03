#include "GollumEnumDefinitions.h"

/**
* @file GollumEnumDefinitions.cpp
* @brief Definitions of enumerators for related sets of parameters, and associated functions to map them to string names.
*/

namespace gollumfit {

/**
* @brief Retrieves the string associated with a given FluxComponent.
*
* This function maps the enumerators of the FluxComponent enum to their corresponding
* string representations. This is useful for obtaining human-readable names for
* flux components.
*
* @param component The FluxComponent enumerator for which the name is requested.
* @return A string representing the name of the flux component.
*
* @throws std::runtime_error If the provided enumerator is not part of the registered flux components.
*/
std::string GetFluxComponentName(const FluxComponent& component){
  switch(component) {
    case FluxComponent::atmConv: return "atmConv";
    case FluxComponent::atmPrompt: return "atmPrompt";
    case FluxComponent::diffuseAstro: return "diffuseAstro";
    default: throw std::runtime_error("Flux component given by enumerator not registered. No label guy found.");
  }
}

/**
* @brief Overload of the output stream operator for FluxComponent.
*
* Allows for easy streaming of FluxComponent's name representation to an output stream such as std::cout.
*
* @param os Reference to the output stream.
* @param component The FluxComponent to be output.
* @return Reference to the output stream after the operation.
*/
std::ostream& operator<<(std::ostream& os, const FluxComponent& component){
  os << GetFluxComponentName(component);
  return os;
}

/**
* @brief Retrieves the string name associated with a given Topology enumerator.
*
* This function maps the enumerators of the Topology enum to their corresponding
* string representations. This is useful for obtaining human-readable names for
* event topologies.
*
* @param topology The Topology enumerator for which the name is requested.
* @return A string representing the name of the topology, or an empty string if the topology is not recognized.
*/
std::string GetTopologyName(const Topology& topology){
  switch(topology) {
    case Topology::shower: return "shower";
    case Topology::track: return "track";
    case Topology::doublebang: return "doublebang";
    case Topology::all: return "all";
    default: return std::string();
  }
}

/**
* @brief Overload of the output stream operator for Topology.
*
* Allows for easy streaming of Topology's name representation to an output stream such as std::cout.
*
* @param os Reference to the output stream.
* @param topology The Topology to be output.
* @return Reference to the output stream after the operation.
*/
std::ostream& operator<<(std::ostream& os, const Topology& topology){
  os << GetTopologyName(topology);
  return os;
}

/**
* @brief Retrieves the string name associated with a given LW::ParticleType enumerator.
*
* This function maps the enumerators of the LW::ParticleType enum to their corresponding
* string representations. This is useful for obtaining human-readable names for
* event topologies.
*
* @param particle The LW::ParticleType enumerator for which the name is requested.
* @return A string representing the name of the particle, or an empty string if the particle is not recognized.
*/
std::string GetParticleName(const LW::ParticleType& particle){
  switch(particle) {
    case LW::ParticleType::NuE: return "nue";
    case LW::ParticleType::NuEBar: return "nuebar";
    case LW::ParticleType::NuMu: return "numu";
    case LW::ParticleType::NuMuBar: return "numubar";
    case LW::ParticleType::NuTau: return "nutau";
    case LW::ParticleType::NuTauBar: return "nutaubar";
    default: return std::string();
  }
}

/**
* @brief Converts a string identifier to its corresponding HadronicParameter enumerator.
*
* Given a string representing the name of a hadronic model parameter, this function returns
* the corresponding enumerator from the HadronicParameter enum. If the name does not correspond
* to any known enumerator, an exception is thrown.
*
* @param had_name The string name of the hadronic parameter.
* @return The HadronicParameter enumerator corresponding to the input string.
* @throws std::runtime_error If the input string does not match any known parameter label.
*/
HadronicParameter GetHadronicParameter(std::string had_name){
  if      ( had_name=="he_K+" ) return HadronicParameter::HEkp;
  else if ( had_name=="he_K-" ) return HadronicParameter::HEkm;
  else if ( had_name=="vhe1_pi+" ) return HadronicParameter::VHE1pip;
  else if ( had_name=="vhe1_pi-" ) return HadronicParameter::VHE1pim;
  else if ( had_name=="vhe3_K+" ) return HadronicParameter::VHE3kp;
  else if ( had_name=="vhe3_K-" ) return HadronicParameter::VHE3km;
  else if ( had_name=="vhe3_pi+" ) return HadronicParameter::VHE3pip;
  else if ( had_name=="vhe3_pi-" ) return HadronicParameter::VHE3pim;
  else if ( had_name=="vhe3_p" ) return HadronicParameter::VHE3p;
  else if ( had_name=="vhe3_n" ) return HadronicParameter::VHE3n;
  else throw std::runtime_error("Impossible hadronic parameter label");
}


/**
* @brief Retrieves the string name of a HadronicParameter enumerator.
*
* This function maps the enumerators of the HadronicParameter enum to their corresponding
* string representations, which can be used for display or logging purposes.
*
* @param had_par The HadronicParameter enumerator for which the name is requested.
* @return A string representing the name of the hadronic parameter, or an empty string if the parameter is not recognized.
*/
std::string GetHadronicParameterName(const HadronicParameter&  had_par){
  switch(had_par) {
    case HadronicParameter::HEkp : return "he_K+";
    case HadronicParameter::HEkm : return "he_K-";
    case HadronicParameter::VHE1pip : return "vhe1_pi+";
    case HadronicParameter::VHE1pim : return "vhe1_pi-";
    case HadronicParameter::VHE3pip : return "vhe3_pi+";
    case HadronicParameter::VHE3pim : return "vhe3_pi-";
    case HadronicParameter::VHE3kp : return "vhe3_K+";
    case HadronicParameter::VHE3km : return "vhe3_K-";
    case HadronicParameter::VHE3p : return "vhe3_p";
    case HadronicParameter::VHE3n : return "vhe3_n";
    default: return std::string();
  }
}

/**
* @brief Overload of the output stream operator for HadronicParameter.
*
* Allows for easy streaming of a HadronicParameter's name representation to an output stream such as std::cout.
*
* @param os Reference to the output stream.
* @param had_par The HadronicParameter to be output.
* @return Reference to the output stream after the operation.
*/
std::ostream& operator<<(std::ostream& os, const HadronicParameter& had_par){
  os << GetHadronicParameterName(had_par);
  return os;
}

/**
* @brief Converts a string identifier to its corresponding CosmicRayParameter enumerator.
*
* Given a string representing the name of a cosmic ray model parameter, this function returns
* the corresponding enumerator from the CosmicRayParameter enum. If the name does not correspond
* to any known enumerator, an exception is thrown.
*
* @param cr_name The string name of the cosmic ray parameter.
* @return The CosmicRayParameter enumerator corresponding to the input string.
* @throws std::runtime_error If the input string does not match any known parameter label.
*/
CosmicRayParameter GetCosmicRayParameter(std::string cr_name){
  if      ( cr_name=="GSF_1") return CosmicRayParameter::GSF1;
  else if ( cr_name=="GSF_2") return CosmicRayParameter::GSF2;
  else if ( cr_name=="GSF_3") return CosmicRayParameter::GSF3;
  else if ( cr_name=="GSF_4") return CosmicRayParameter::GSF4;
  else if ( cr_name=="GSF_5") return CosmicRayParameter::GSF5;
  else if ( cr_name=="GSF_6") return CosmicRayParameter::GSF6;
  else throw std::runtime_error("Impossible cosmicray parameter label");
}


/**
* @brief Retrieves the string name of a CosmicRayParameter enumerator.
*
* This function maps the enumerators of the CosmicRayParameter enum to their corresponding
* string representations, which can be used for display or logging purposes.
*
* @param cr_par The CosmicRayParameter enumerator for which the name is requested.
* @return A string representing the name of the cosmic ray parameter, or an empty string if the parameter is not recognized.
*/
std::string GetCosmicRayParameterName(const CosmicRayParameter& cr_par){
  switch(cr_par) {
    case CosmicRayParameter::GSF1 : return "GSF_1";
    case CosmicRayParameter::GSF2 : return "GSF_2";
    case CosmicRayParameter::GSF3 : return "GSF_3";
    case CosmicRayParameter::GSF4 : return "GSF_4";
    case CosmicRayParameter::GSF5 : return "GSF_5";
    case CosmicRayParameter::GSF6 : return "GSF_6";
    default: return std::string();
  }
}


/**
* @brief Overload of the output stream operator for CosmicRayParameter.
*
* Allows the name associated with a CosmicRayParameter enumerator to be streamed directly to an output stream.
* This facilitates logging and printing of the cosmic ray parameter's name.
*
* @param os The output stream to which the name will be written.
* @param cr_par The cosmic ray parameter enumerator to be written to the output stream.
* @return Reference to the output stream.
*/
std::ostream& operator<<(std::ostream& os, const CosmicRayParameter& cr_par){
  os << GetCosmicRayParameterName(cr_par);
  return os;
}


} // close namespace
