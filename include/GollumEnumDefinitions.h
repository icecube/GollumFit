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
* Each value corresponds to a secondary particle species and energy scale whose
* production cross section is varied as a nuisance parameter. The naming convention
* is `<energy regime>_<particle>`, where HE = high-energy (158 GeV), VHE1 = very-high-energy
* scale 1 (20 TeV), and VHE3 = very-high-energy scale 3 (2 PeV).
*/
enum class HadronicParameter {
    HEkp,    ///< High-energy K+ production at 158 GeV
    HEkm,    ///< High-energy K- production at 158 GeV
    VHE1pip, ///< Very-high-energy π+ production at 20 TeV
    VHE1pim, ///< Very-high-energy π- production at 20 TeV
    VHE3kp,  ///< Very-high-energy K+ production at 2 PeV
    VHE3km,  ///< Very-high-energy K- production at 2 PeV
    VHE3pip, ///< Very-high-energy π+ production at 2 PeV
    VHE3pim, ///< Very-high-energy π- production at 2 PeV
    VHE3p,   ///< Very-high-energy proton production at 2 PeV
    VHE3n    ///< Very-high-energy neutron production at 2 PeV
};
HadronicParameter GetHadronicParameter(std::string had_name);
std::string GetHadronicParameterName(const HadronicParameter& had_par);
std::ostream& operator<<(std::ostream& os, const HadronicParameter& had_par);

/**
* @enum CosmicRayParameter
* @brief Enumerates cosmic ray model parameters.
*
* Each value corresponds to one principal component of the Global Spline Fit (GSF)
* cosmic ray flux model. Varying these parameters shifts the cosmic ray energy spectrum
* within the uncertainties of the GSF model, propagating into the atmospheric neutrino flux.
*/
enum class CosmicRayParameter {
    GSF1, ///< GSF cosmic ray principal component 1
    GSF2, ///< GSF cosmic ray principal component 2
    GSF3, ///< GSF cosmic ray principal component 3
    GSF4, ///< GSF cosmic ray principal component 4
    GSF5, ///< GSF cosmic ray principal component 5
    GSF6  ///< GSF cosmic ray principal component 6
};
CosmicRayParameter GetCosmicRayParameter(std::string cr_name);
std::string GetCosmicRayParameterName(const CosmicRayParameter& cr_par);
std::ostream& operator<<(std::ostream& os, const CosmicRayParameter& cr_par);

std::string GetParticleName(const LW::ParticleType& particle);

} // namespace gollumfit

#endif /* GOLLUMENUMDEFINITIONS_H_ */
