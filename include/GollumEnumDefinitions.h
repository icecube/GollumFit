#ifndef GOLLUMENUMDEFINITIONS_H_
#define GOLLUMENUMDEFINITIONS_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <LeptonWeighter/ParticleType.h>

/**
* @file GollumEnumDefinitions.h
* @brief Declarations of enumerators for related sets of parameters, and associated functions to map them to string names.
*/

namespace gollumfit{

/**
* @enum FluxComponent
* @brief Enumerates the components of the flux model.
*
* This enumeration represents the different components of the neutrino flux such as
* atmospheric conventional (atmConv), atmospheric prompt (atmPrompt), and
* diffuse astrophysical (diffuseAstro) neutrinos.
*/
enum class FluxComponent {atmConv, atmPrompt, diffuseAstro};
std::string GetFluxComponentName(const FluxComponent& component);
std::ostream& operator<<(std::ostream& os, const FluxComponent& component);

/**
* @enum Topology
* @brief Enumerates the event topology types.
*
* This enumeration represents the different event topology types such as
* shower, track, doublebang, or all.
*/
enum class Topology {shower = 0, track = 1, doublebang = 2, all = 3};
std::string GetTopologyName(const Topology& topology);
std::ostream& operator<<(std::ostream& os, const Topology& model);

/**
* @enum DiffuseFitType
* @brief Enumerates the types of diffuse fit models.
*
* This enumeration represents the different diffuse fit models such as
* SinglePowerLaw and BrokenPowerLaw.
*/
enum class DiffuseFitType {SinglePowerLaw, BrokenPowerLaw};

/**
* @enum HadronicParameter
* @brief Enumerates hadronic interaction model parameters.
*
* This enumeration represents different parameters of hadronic interaction models.
*/
enum class HadronicParameter {
    HEkp, 
    HEkm, 
    VHE1pip, 
    VHE1pim, 
    VHE3kp, 
    VHE3km, 
    VHE3pip, 
    VHE3pim, 
    VHE3p, 
    VHE3n
};
HadronicParameter GetHadronicParameter(std::string had_name);
std::string GetHadronicParameterName(const HadronicParameter& had_par);
std::ostream& operator<<(std::ostream& os, const HadronicParameter& had_par);

/**
* @enum CosmicRayParameter
* @brief Enumerates cosmic ray model parameters.
*
* This enumeration represents different parameters used in cosmic ray models, specifically
* the Global Spline Fit (GSF) model parameters.
*/
enum class CosmicRayParameter {GSF1, GSF2, GSF3, GSF4, GSF5, GSF6};
CosmicRayParameter GetCosmicRayParameter(std::string cr_name);
std::string GetCosmicRayParameterName(const CosmicRayParameter& cr_par);
std::ostream& operator<<(std::ostream& os, const CosmicRayParameter& cr_par);

// this function should totally get moved to LeptonWeighter CA. Note this aschneider on the review thingy.
std::string GetParticleName(const LW::ParticleType& particle);

} // close namespace

#endif /* GOLLUMENUMDEFINITIONS_H_ */
