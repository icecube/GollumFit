#include "Event.h"

/**
* @file Event.cpp
* @brief Event class
* 
* Defines operator overload to stream variables in the event class.
*/


std::ostream& operator<<(std::ostream& os, const Event& e){
    std::cout << " energy: " << e.energy << '\n'
    << " zenith: " << e.zenith << '\n'
    << " topology: " << (int)e.topology << '\n'
    << " primaryType: " << (int)e.primaryType << '\n'
    << " primaryZenith: " << e.primaryZenith << '\n'
    << " primaryAzimuth: " << e.primaryAzimuth << '\n'
    << " primaryEnergy: " << e.primaryEnergy << '\n'
    << " cachedWeight: " << e.cachedWeight << '\n'
    << " cachedConvWeight: " << e.cachedConvWeight << '\n'
    << " cachedPromptWeight: " << e.cachedPromptWeight << '\n'
    << " cachedAstroWeight: " << e.cachedAstroWeight << '\n'
    << " cachedDOMEffConv: " << e.cachedDOMEffConv<< '\n'
    << " cachedDOMEffPrompt: " << e.cachedDOMEffPrompt << '\n'
    << " cachedDOMEffAstro: " << e.cachedDOMEffAstro << '\n'
    << " cachedHoleIceConv: " << e.cachedHoleIceConv << '\n'
    << " cachedHoleIcePrompt: " << e.cachedHoleIcePrompt << '\n'
    << " cachedHoleIceAstro: " << e.cachedHoleIceAstro << '\n'
    << " cachedHadronicHEkp: " << e.cachedHadronicHEkp << '\n'
    << " cachedHadronicHEkm: " << e.cachedHadronicHEkm << '\n'
    << " cachedHadronicVHE1pip: " << e.cachedHadronicVHE1pip << '\n'
    << " cachedHadronicVHE1pim: " << e.cachedHadronicVHE1pim << '\n'
    << " cachedHadronicVHE3kp: " << e.cachedHadronicVHE3kp << '\n'
    << " cachedHadronicVHE3km: " << e.cachedHadronicVHE3km << '\n'
    << " cachedHadronicVHE3pip: " << e.cachedHadronicVHE3pip << '\n'
    << " cachedHadronicVHE3pim: " << e.cachedHadronicVHE3pim << '\n'
    << " cachedHadronicVHE3p: " << e.cachedHadronicVHE3p << '\n'
    << " cachedHadronicVHE3n: " << e.cachedHadronicVHE3n << '\n'
    << " cachedCosmicRay1: " << e.cachedCosmicRay1 << '\n'
    << " cachedCosmicRay2: " << e.cachedCosmicRay2 << '\n'
    << " cachedCosmicRay3: " << e.cachedCosmicRay3 << '\n'
    << " cachedCosmicRay4: " << e.cachedCosmicRay4 << '\n'
    << " cachedCosmicRay5: " << e.cachedCosmicRay5 << '\n'
    << " cachedCosmicRay6: " << e.cachedCosmicRay6 << '\n'
    << " cachedAtmDensity: " << e.cachedAtmDensity << '\n'
    << " cachedKaonLosses: " << e.cachedKaonLosses << '\n'
    << " cachedIceGrad0: " << e.cachedIceGrad0 << '\n'
    << " cachedIceGrad1: " << e.cachedIceGrad1 << '\n'
    << " cachedIceGrad2: " << e.cachedIceGrad2 << '\n'
    << " cachedIceGrad3: " << e.cachedIceGrad3 << '\n'
    << " cachedIceGrad4: " << e.cachedIceGrad4 << '\n'
    << " cachedIceGrad5: " << e.cachedIceGrad5 << '\n'
    << " cachedIceGrad6: " << e.cachedIceGrad6 << '\n'
    << " cachedIceGrad7: " << e.cachedIceGrad7 << '\n'
    << " cachedIceGrad8: " << e.cachedIceGrad8 << '\n';
    return(os);
}

herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data){
    std::set<std::string>* items=static_cast<std::set<std::string>*>(operator_data);
    items->insert(member_name);
    return(0);
}
