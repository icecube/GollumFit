/*
 * GollumParameters.cpp
 *
 *  Created on: Nov 21, 2017
 *      Author: aschneider
 */

#include "GollumParameters.h"
#include "Event.h"

/**
* @file GollumParameters.cpp
*/

namespace gollumfit {

/*************************************************************************************************************
 * METAGOLLUM
 * **********************************************************************************************************/

std::ostream& operator<<(std::ostream& os, const FitParameters& fp){
  os << "convNorm                  : " << fp.convNorm << std::endl;
  os << "promptNorm                : " << fp.promptNorm << std::endl;
  os << "zenithCorrection          : " << fp.zenithCorrection << std::endl;
  os << "kaonLosses                : " << fp.kaonLosses << std::endl;
  os << "hadronicHEkp              : " << fp.hadronicHEkp << std::endl;
  os << "hadronicHEkm              : " << fp.hadronicHEkm << std::endl;
  os << "hadronicVHE1pip           : " << fp.hadronicVHE1pip << std::endl;
  os << "hadronicVHE1pim           : " << fp.hadronicVHE1pim << std::endl;
  os << "hadronicVHE3kp            : " << fp.hadronicVHE3kp << std::endl;
  os << "hadronicVHE3km            : " << fp.hadronicVHE3km << std::endl;
  os << "hadronicVHE3pip           : " << fp.hadronicVHE3pip << std::endl;
  os << "hadronicVHE3pim           : " << fp.hadronicVHE3pim << std::endl;
  os << "hadronicVHE3p             : " << fp.hadronicVHE3p << std::endl;
  os << "hadronicVHE3n             : " << fp.hadronicVHE3n << std::endl;
  os << "cosmicRay1                : " << fp.cosmicRay1 << std::endl;
  os << "cosmicRay2                : " << fp.cosmicRay2 << std::endl;
  os << "cosmicRay3                : " << fp.cosmicRay3 << std::endl;
  os << "cosmicRay4                : " << fp.cosmicRay4 << std::endl;
  os << "cosmicRay5                : " << fp.cosmicRay5 << std::endl;
  os << "cosmicRay6                : " << fp.cosmicRay6 << std::endl;
  os << "icegrad0                  : " << fp.icegrad0 << std::endl;
  os << "icegrad1                  : " << fp.icegrad1 << std::endl;
  os << "icegrad2                  : " << fp.icegrad2 << std::endl;
  os << "icegrad3                  : " << fp.icegrad3 << std::endl;
  os << "icegrad4                  : " << fp.icegrad4 << std::endl;
  os << "icegrad5                  : " << fp.icegrad5 << std::endl;
  os << "icegrad6                  : " << fp.icegrad6 << std::endl;
  os << "icegrad7                  : " << fp.icegrad7 << std::endl;
  os << "icegrad8                  : " << fp.icegrad8 << std::endl;
  os << "domEfficiency             : " << fp.domEfficiency << std::endl;
  os << "holeiceForward            : " << fp.holeiceForward << std::endl;
  os << "astroNorm                 : " << fp.astroNorm << std::endl;
  os << "astroDeltaGamma           : " << fp.astroDeltaGamma << std::endl;
  os << "astroDeltaGammaSec        : " << fp.astroDeltaGammaSec << std::endl;
  os << "astroPivot                : " << fp.astroPivot << std::endl;
  os << "NeutrinoAntineutrinoRatio : " << fp.NeutrinoAntineutrinoRatio << std::endl;
  os << "nuxs                      : " << fp.nuxs << std::endl;
  os << "nubarxs                   : " << fp.nubarxs << std::endl;
  return os;
}

} // namespace gollumfit
