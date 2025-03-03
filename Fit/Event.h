#ifndef EVENT_H_
#define EVENT_H_

#include <iostream>
#include <set>

#include <boost/math/constants/constants.hpp>
#include <LeptonWeighter/ParticleType.h>
#include <PhysTools/tableio.h>
#include "GollumEnumDefinitions.h"
#include "FastMode.h"

/**
* @file Event.h
* @brief Event class
* 
* Defines the event class. Also defines a couple helper functions to combine events and read files.
*/


/**
* @brief Iterates over the members of a group, collecting each member's name.
* 
* This function is designed to be used as an operator with H5Giterate or similar
* functions that require an H5L_iterate_t callback. It will cast the operator_data
* to a std::set<std::string> and insert the member_name into the set. Typical
* usage is to collect all table names (or other member names) within a HDF5 group.
*
* @param group_id The identifier for the group being iterated over. This is usually
*                obtained from a call to one of the HDF5 library's group functions.
* @param member_name A pointer to the null-terminated name of the current member.
*                   The name is retrieved by the library calling this function.
* @param operator_data A pointer to user-defined data. In this case, it should be
*                      a pointer to a std::set<std::string> which collects the member names.
* 
* @return Returns zero to indicate success. If a non-zero value is returned, the 
*         iterating function will stop and return that value.
*
* @note The function always returns 0, indicating successful iteration. If you need
*       to handle errors or interruptions, additional logic should be added to handle
*       such cases.
* 
* Example usage:
* @code
* std::set<std::string> tableNames;
* herr_t status = H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
* if (status < 0) {
*    // Handle error
* }
* // tableNames now contains all the names of the members of the root group
* @endcode
* Presently collectTableNames is used within the function readFile, in the following way:
* @code
* H5File h5file(filePath);
* std::set<std::string> tables;
* H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
* @endcode
*
*/
herr_t collectTableNames(hid_t group_id, const char * member_name, void* operator_data);

/**
* @brief Basic event structure used throughout GollumFit.
* 
* 
*/
class Event{
public:

  // LW/LI compatibility
  LW::ParticleType final_state_particle_0;
  LW::ParticleType final_state_particle_1;
  LW::ParticleType primaryType;
  float primaryEnergy;
  float primaryAzimuth;
  float primaryZenith;
  float totalColumnDepth;
  float intX,intY;

  // meta-event quantities
  int num_events;

  // reconstructed quantities
  float energy;
  float zenith;

  // topology
  unsigned int topology;

  // cached quantities
  float cachedWeight;

  // cached fluxes
  double cachedConvWeight;
  double cachedPromptWeight;
  double cachedAstroWeight;

  // cached detector systematic
  // holeice cached
  double cachedHoleIceConv;
  double cachedHoleIcePrompt;
  double cachedHoleIceAstro;

  // domeff cached
  double cachedDOMEffConv;
  double cachedDOMEffPrompt;
  double cachedDOMEffAstro;

  //hadronic cached
  double cachedHadronicHEkp;            
  double cachedHadronicHEkm;            
  double cachedHadronicVHE1pip;            
  double cachedHadronicVHE1pim;            
  double cachedHadronicVHE3kp;            
  double cachedHadronicVHE3km;            
  double cachedHadronicVHE3pip;            
  double cachedHadronicVHE3pim;            
  double cachedHadronicVHE3p;            
  double cachedHadronicVHE3n;            

  //cosmic ray cached
  double cachedCosmicRay1;    
  double cachedCosmicRay2;    
  double cachedCosmicRay3;    
  double cachedCosmicRay4;    
  double cachedCosmicRay5;    
  double cachedCosmicRay6;    

  //flux cached
  double cachedAtmDensity;
  double cachedKaonLosses;

  //ice gradients cached
  double cachedIceGrad0;
  double cachedIceGrad1;
  double cachedIceGrad2;
  double cachedIceGrad3;
  double cachedIceGrad4;
  double cachedIceGrad5;
  double cachedIceGrad6;
  double cachedIceGrad7;
  double cachedIceGrad8;

  /**
  * @brief Default constructor.
  */
  Event():
    // LW/LI compatibility
    final_state_particle_0(LW::ParticleType::unknown),
    final_state_particle_1(LW::ParticleType::unknown),
    primaryType(LW::ParticleType::unknown),
    primaryEnergy(std::numeric_limits<float>::quiet_NaN()),
    primaryAzimuth(std::numeric_limits<float>::quiet_NaN()),
    primaryZenith(std::numeric_limits<float>::quiet_NaN()),
    totalColumnDepth(std::numeric_limits<float>::quiet_NaN()),
    intX(std::numeric_limits<float>::quiet_NaN()),
    intY(std::numeric_limits<float>::quiet_NaN()),
    // meta-event quantities
    num_events(1),
    // reconstruction quantities
    energy(std::numeric_limits<float>::quiet_NaN()),
    zenith(std::numeric_limits<float>::quiet_NaN()),
    // topology
    topology(0),
    // cached quantities
    cachedWeight(0),
    // flux caches
    cachedConvWeight(0),
    cachedPromptWeight(0),
    cachedAstroWeight(0),
    // holeice cache
    cachedHoleIceConv(0),
    cachedHoleIcePrompt(0),
    cachedHoleIceAstro(0),
    // dom efficiency cache
    cachedDOMEffConv(0),
    cachedDOMEffPrompt(0),
    cachedDOMEffAstro(0),
    //hadronic cached
    cachedHadronicHEkp(0),            
    cachedHadronicHEkm(0),            
    cachedHadronicVHE1pip(0),            
    cachedHadronicVHE1pim(0),            
    cachedHadronicVHE3kp(0),            
    cachedHadronicVHE3km(0),            
    cachedHadronicVHE3pip(0),            
    cachedHadronicVHE3pim(0),            
    cachedHadronicVHE3p(0),            
    cachedHadronicVHE3n(0),            
    //cosmic ray cached
    cachedCosmicRay1(0),    
    cachedCosmicRay2(0),    
    cachedCosmicRay3(0),    
    cachedCosmicRay4(0),    
    cachedCosmicRay5(0),    
    cachedCosmicRay6(0),    
    //flux cached
    cachedAtmDensity(0),
    cachedKaonLosses(0),
    //ice gradients cached
    cachedIceGrad0(0),
    cachedIceGrad1(0),
    cachedIceGrad2(0),
    cachedIceGrad3(0),
    cachedIceGrad4(0),
    cachedIceGrad5(0),
    cachedIceGrad6(0),
    cachedIceGrad7(0),
    cachedIceGrad8(0)
  {}

  /**
  * @brief Method to define == operator in a natural way to compare events.
  */
  bool operator==(const Event& e) {
    return  (((final_state_particle_0 == e.final_state_particle_0)) &&
            ((final_state_particle_1 == e.final_state_particle_1)) &&
            ((primaryType == e.primaryType)) &&
            ((primaryEnergy == e.primaryEnergy) || (std::isnan(primaryEnergy) && std::isnan(e.primaryEnergy))) &&
            ((primaryAzimuth == e.primaryAzimuth) || (std::isnan(primaryAzimuth) && std::isnan(e.primaryAzimuth))) &&
            ((primaryZenith == e.primaryZenith) || (std::isnan(primaryZenith) && std::isnan(e.primaryZenith))) &&
            ((totalColumnDepth == e.totalColumnDepth) || (std::isnan(totalColumnDepth) && std::isnan(e.totalColumnDepth))) &&
            ((intX == e.intX) || (std::isnan(intX) && std::isnan(e.intX))) &&
            ((intY == e.intY) || (std::isnan(intY) && std::isnan(e.intY))) &&
            ((num_events == e.num_events) || (std::isnan(num_events) && std::isnan(e.num_events))) &&
            ((energy == e.energy) || (std::isnan(energy) && std::isnan(e.energy))) &&
            ((zenith == e.zenith) || (std::isnan(zenith) && std::isnan(e.zenith))) &&
            ((topology == e.topology) || (std::isnan(topology) && std::isnan(e.topology))) &&
            ((cachedConvWeight == e.cachedConvWeight) || (std::isnan(cachedConvWeight) && std::isnan(e.cachedConvWeight))) &&
            ((cachedPromptWeight == e.cachedPromptWeight) || (std::isnan(cachedPromptWeight) && std::isnan(e.cachedPromptWeight))) &&
            ((cachedAstroWeight == e.cachedAstroWeight) || (std::isnan(cachedAstroWeight) && std::isnan(e.cachedAstroWeight))) &&
            ((cachedWeight == e.cachedWeight) || (std::isnan(cachedWeight) && std::isnan(e.cachedWeight))) &&
            ((cachedHoleIceConv == e.cachedHoleIceConv) || (std::isnan(cachedHoleIceConv) && std::isnan(e.cachedHoleIceConv))) &&
            ((cachedHoleIcePrompt == e.cachedHoleIcePrompt) || (std::isnan(cachedHoleIcePrompt) && std::isnan(e.cachedHoleIcePrompt))) &&
            ((cachedHoleIceAstro== e.cachedHoleIceAstro) || (std::isnan(cachedHoleIceAstro) && std::isnan(e.cachedHoleIceAstro))) &&
            ((cachedDOMEffConv == e.cachedDOMEffConv) || (std::isnan(cachedDOMEffConv) && std::isnan(e.cachedDOMEffConv))) &&
            ((cachedDOMEffPrompt == e.cachedDOMEffPrompt) || (std::isnan(cachedDOMEffPrompt) && std::isnan(e.cachedDOMEffPrompt))) &&
            ((cachedDOMEffAstro == e.cachedDOMEffAstro) || (std::isnan(cachedDOMEffAstro) && std::isnan(e.cachedDOMEffAstro))) &&
            ((cachedHadronicHEkp == e.cachedHadronicHEkp) || (std::isnan(cachedHadronicHEkp) && std::isnan(e.cachedHadronicHEkp))) &&            
            ((cachedHadronicHEkm == e.cachedHadronicHEkm) || (std::isnan(cachedHadronicHEkm) && std::isnan(e.cachedHadronicHEkm))) &&            
            ((cachedHadronicVHE1pip == e.cachedHadronicVHE1pip) || (std::isnan(cachedHadronicVHE1pip) && std::isnan(e.cachedHadronicVHE1pip))) &&            
            ((cachedHadronicVHE1pim == e.cachedHadronicVHE1pim) || (std::isnan(cachedHadronicVHE1pim) && std::isnan(e.cachedHadronicVHE1pim))) &&            
            ((cachedHadronicVHE3kp == e.cachedHadronicVHE3kp) || (std::isnan(cachedHadronicVHE3kp) && std::isnan(e.cachedHadronicVHE3kp))) &&            
            ((cachedHadronicVHE3km == e.cachedHadronicVHE3km) || (std::isnan(cachedHadronicVHE3km) && std::isnan(e.cachedHadronicVHE3km))) &&            
            ((cachedHadronicVHE3pip == e.cachedHadronicVHE3pip) || (std::isnan(cachedHadronicVHE3pip) && std::isnan(e.cachedHadronicVHE3pip))) &&            
            ((cachedHadronicVHE3pim == e.cachedHadronicVHE3pim) || (std::isnan(cachedHadronicVHE3pim) && std::isnan(e.cachedHadronicVHE3pim))) &&            
            ((cachedHadronicVHE3p == e.cachedHadronicVHE3p) || (std::isnan(cachedHadronicVHE3p) && std::isnan(e.cachedHadronicVHE3p))) &&            
            ((cachedHadronicVHE3n == e.cachedHadronicVHE3n) || (std::isnan(cachedHadronicVHE3n) && std::isnan(e.cachedHadronicVHE3n))) &&            
            ((cachedCosmicRay1 == e.cachedCosmicRay1) || (std::isnan(cachedCosmicRay1) && std::isnan(e.cachedCosmicRay1))) &&    
            ((cachedCosmicRay2 == e.cachedCosmicRay2) || (std::isnan(cachedCosmicRay2) && std::isnan(e.cachedCosmicRay2))) &&    
            ((cachedCosmicRay3 == e.cachedCosmicRay3) || (std::isnan(cachedCosmicRay3) && std::isnan(e.cachedCosmicRay3))) &&    
            ((cachedCosmicRay4 == e.cachedCosmicRay4) || (std::isnan(cachedCosmicRay4) && std::isnan(e.cachedCosmicRay4))) &&    
            ((cachedCosmicRay5 == e.cachedCosmicRay5) || (std::isnan(cachedCosmicRay5) && std::isnan(e.cachedCosmicRay5))) &&    
            ((cachedCosmicRay6 == e.cachedCosmicRay6) || (std::isnan(cachedCosmicRay6) && std::isnan(e.cachedCosmicRay6))) &&    
            ((cachedAtmDensity == e.cachedAtmDensity) || (std::isnan(cachedAtmDensity) && std::isnan(e.cachedAtmDensity))) &&
            ((cachedKaonLosses == e.cachedKaonLosses) || (std::isnan(cachedKaonLosses) && std::isnan(e.cachedKaonLosses))) &&
            ((cachedIceGrad0 == e.cachedIceGrad0) || (std::isnan(cachedIceGrad0) && std::isnan(e.cachedIceGrad0))) &&
            ((cachedIceGrad1 == e.cachedIceGrad1) || (std::isnan(cachedIceGrad1) && std::isnan(e.cachedIceGrad1))) &&
            ((cachedIceGrad2 == e.cachedIceGrad2) || (std::isnan(cachedIceGrad2) && std::isnan(e.cachedIceGrad2))) &&
            ((cachedIceGrad3 == e.cachedIceGrad3) || (std::isnan(cachedIceGrad3) && std::isnan(e.cachedIceGrad3))) &&
            ((cachedIceGrad4 == e.cachedIceGrad4) || (std::isnan(cachedIceGrad4) && std::isnan(e.cachedIceGrad4))) &&
            ((cachedIceGrad5 == e.cachedIceGrad5) || (std::isnan(cachedIceGrad5) && std::isnan(e.cachedIceGrad5))) &&
            ((cachedIceGrad6 == e.cachedIceGrad6) || (std::isnan(cachedIceGrad6) && std::isnan(e.cachedIceGrad6))) &&
            ((cachedIceGrad7 == e.cachedIceGrad7) || (std::isnan(cachedIceGrad7) && std::isnan(e.cachedIceGrad7))) &&
            ((cachedIceGrad8 == e.cachedIceGrad8) || (std::isnan(cachedIceGrad8) && std::isnan(e.cachedIceGrad8)))
            );
  }

  double Get_MetaWeight() const {
    double meta_weight = (cachedConvWeight + cachedPromptWeight + cachedAstroWeight);
    return meta_weight;
  }
 
  /**
  * @brief Method to define += operator for events.
  * 
  * Provides an interpretation for adding events. In particular, the weights are combined in a way such that
  * the new event is correctly weighted accounting for the contributions of the previous two events.
  * 
  */
  void operator+=(const Event& e) {
    // if the meta_weight is not 1, but its the oneWeight
    // the weighted fast mode is been turn on. This is a
    // super sharp knive. Appropiate magic needs to be done in the
    // division guy. CA.
    double meta_weight = e.Get_MetaWeight();
    assert((not std::isnan(meta_weight)) and (meta_weight >= 0.0));

    num_events += e.num_events;
    energy += e.energy*meta_weight;
    zenith += e.zenith*meta_weight;
    primaryEnergy += e.primaryEnergy*meta_weight;
    primaryZenith += e.primaryZenith*meta_weight;
    primaryAzimuth += e.primaryAzimuth*meta_weight;
    totalColumnDepth += e.totalColumnDepth*meta_weight;
    intX += e.intX*meta_weight;
    intY += e.intY*meta_weight;
    cachedAtmDensity += e.cachedAtmDensity*meta_weight;
    cachedKaonLosses += e.cachedKaonLosses*meta_weight;
    cachedIceGrad0 += e.cachedIceGrad0*meta_weight;
    cachedIceGrad1 += e.cachedIceGrad1*meta_weight;
    cachedIceGrad2 += e.cachedIceGrad2*meta_weight;
    cachedIceGrad3 += e.cachedIceGrad3*meta_weight;
    cachedIceGrad4 += e.cachedIceGrad4*meta_weight;
    cachedIceGrad5 += e.cachedIceGrad5*meta_weight;
    cachedIceGrad6 += e.cachedIceGrad6*meta_weight;
    cachedIceGrad7 += e.cachedIceGrad7*meta_weight;
    cachedIceGrad8 += e.cachedIceGrad8*meta_weight;
    // we accumulate the rates
    cachedWeight += e.cachedWeight;
    cachedConvWeight += e.cachedConvWeight;
    cachedPromptWeight += e.cachedPromptWeight;
    cachedAstroWeight += e.cachedAstroWeight;
    cachedHadronicHEkp += e.cachedHadronicHEkp;
    cachedHadronicHEkm += e.cachedHadronicHEkm;
    cachedHadronicVHE1pip += e.cachedHadronicVHE1pip;
    cachedHadronicVHE1pim += e.cachedHadronicVHE1pim;
    cachedHadronicVHE3kp += e.cachedHadronicVHE3kp;
    cachedHadronicVHE3km += e.cachedHadronicVHE3km;
    cachedHadronicVHE3pip += e.cachedHadronicVHE3pip;
    cachedHadronicVHE3pim += e.cachedHadronicVHE3pim;
    cachedHadronicVHE3p += e.cachedHadronicVHE3p;
    cachedHadronicVHE3n += e.cachedHadronicVHE3n;
    cachedCosmicRay1 += e.cachedCosmicRay1;
    cachedCosmicRay2 += e.cachedCosmicRay2;
    cachedCosmicRay3 += e.cachedCosmicRay3;
    cachedCosmicRay4 += e.cachedCosmicRay4;
    cachedCosmicRay5 += e.cachedCosmicRay5;
    cachedCosmicRay6 += e.cachedCosmicRay6;
    // we average the expectation from systematics variation in reco bin
    cachedHoleIceConv = log10(pow(10.,cachedHoleIceConv) + pow(10.,e.cachedHoleIceConv)*meta_weight);
    cachedHoleIcePrompt = log10(pow(10.,cachedHoleIcePrompt) + pow(10.,e.cachedHoleIcePrompt)*meta_weight);
    cachedHoleIceAstro = log10(pow(10.,cachedHoleIceAstro) + pow(10.,e.cachedHoleIceAstro)*meta_weight);
    cachedDOMEffConv = log10(pow(10.,cachedDOMEffConv) + pow(10.,e.cachedDOMEffConv)*meta_weight);
    cachedDOMEffPrompt = log10(pow(10.,cachedDOMEffPrompt) + pow(10.,e.cachedDOMEffPrompt)*meta_weight);
    cachedDOMEffAstro = log10(pow(10.,cachedDOMEffAstro) + pow(10.,e.cachedDOMEffAstro)*meta_weight);
  }

  /**
  * @brief Method to define /= operator for events.
  * 
  * Provides an interpretation for adding events. In particular, the weights are combined in a way such that
  * the new event is correctly weighted accounting for the contributions of the previous two events.
  * 
  */
  void operator/=(double d) {
    energy /= d;
    zenith /= d;
    primaryEnergy /= d;
    primaryZenith /= d;
    primaryAzimuth /= d;
    totalColumnDepth /= d;
    intX /= d;
    intY /= d;
    cachedAtmDensity /= d;
    cachedKaonLosses /= d;
    cachedIceGrad0 /= d;
    cachedIceGrad1 /= d;
    cachedIceGrad2 /= d;
    cachedIceGrad3 /= d;
    cachedIceGrad4 /= d;
    cachedIceGrad5 /= d;
    cachedIceGrad6 /= d;
    cachedIceGrad7 /= d;
    cachedIceGrad8 /= d;
    cachedHoleIceConv = log10(pow(10.,cachedHoleIceConv)/d);
    cachedHoleIcePrompt = log10(pow(10.,cachedHoleIcePrompt)/d);
    cachedHoleIceAstro = log10(pow(10.,cachedHoleIceAstro)/d);
    cachedDOMEffConv = log10(pow(10.,cachedDOMEffConv)/d);
    cachedDOMEffPrompt = log10(pow(10.,cachedDOMEffPrompt)/d);
    cachedDOMEffAstro= log10(pow(10.,cachedDOMEffAstro)/d);
  }

};

/**
* @brief Function to define the << operator to cout event details.
* 
*/
std::ostream& operator<<(std::ostream& os, const Event& e);

/**
* @brief Combines a vector of events into a single weighted meta event.
* 
* This function takes a vector of events and combines them by summing their attributes after
* normalizing each by its respective meta weight. The final combined meta event is then
* normalized by the total sum of weights. This function will throw an exception if the total
* weight is found to be negative, as this is an undefined condition for the combination process.
*
* @tparam dataType The data type of the events within the vector. This type must have an 
*       overloaded division operator to allow for normalization by weight, an overloaded
*         addition operator for summing events, and must support the Get_MetaWeight method 
*         to retrieve the specific weight of an event. Event objects will be made from this data type.
* @param events A constant reference to a std::vector containing the events to be combined.
* @return Event The combined meta event as an Event object.
* @exception std::runtime_error Thrown if the total weight after combination is negative.
*
* Usage example:
* @code
*   std::vector<dataType> eventVec = {...};
* Event combinedEvent = combine_events(eventVec);
* @endcode
*
* @note If the total weight is zero, the function will return the first event without any changes.
*/
template<typename dataType>
Event combine_events(const std::vector<dataType>& events) {
  gollumfit::fastmode::detail::dereference<dataType> dr;
  Event meta_e(events[0]);
  double w = meta_e.Get_MetaWeight();
  meta_e /= (1./meta_e.Get_MetaWeight());
  for(unsigned int i=1; i<events.size(); ++i) {
      meta_e += events[i];
      w += dr(events[i]).Get_MetaWeight();
  }
  if(w < 0)
    throw std::runtime_error("Spot a negative weight while combining events. Kill.");
  if(w > 0){
    meta_e /= w;
  }
  //std::cout << w << std::endl;
  //std::cout << meta_e << std::endl;
  // note that if w = 0; none of the events have weight and we just return the guy.
  return meta_e;
}

namespace H5Load {

namespace sterile {

/**
* @brief Reads data from a file and processes it using a callback function.
* 
* This function opens a file and reads several tables of data within it. Each table is
* processed to construct Event objects. An action, specified by a callback function, is performed 
* on each Event object. If a specific selectionStart is set, only a subset of the data
* according to some criterion is processed.
* 
* @tparam CallbackType The type of the callback function to be called for each Event.
* @param filePath The path to the file containing the data.
* @param energyName The name of the table within the file from which to read the energy data.
* @param selectionStart The starting value for selecting a subset of the data; if set to -1, 
*       all data is processed.
* @param action The callback function to be called on each Event. This function should take
*       two parameters: a RecordID and an Event.
* 
* @throw std::runtime_error If the file cannot be opened or if no tables are found within it.
* 
* Usage example:
* @code
*   void myAction(RecordID id, Event& e) {
*       // Do something with the event
*   }
*   
* readFile<MyActionType>("path/to/datafile.h5", "SomeEnergyTable", 0, myAction);
* @endcode
* 
* @note The tables read are determined by checking the presence of table names in the file.
*       The specific set of tables processed is hardcoded in the function.
*/
template<typename CallbackType>
void readFile(const std::string& filePath, const std::string& energyName, const double selectionStart, CallbackType action){
    using namespace phys_tools::cts;
    using namespace phys_tools::tableio;
    H5File h5file(filePath);
    if(!h5file)
        throw std::runtime_error("Unable to open "+filePath);
    std::set<std::string> tables;
    H5Giterate(h5file,"/",NULL,&collectTableNames,&tables);
    if(tables.empty())
        throw std::runtime_error(filePath+" contains no tables");
    std::map<RecordID,Event> intermediateData;

    using double_value_field = TableRow<field<double,CTS("value")>>;

    if(tables.count("PrimaryType")){
        readTable<double_value_field>(h5file,"PrimaryType", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryType=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    if(tables.count("FinalStateX")){
        readTable<double_value_field>(h5file,"FinalStateX", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.intX=p.get<CTS("value")>();
            });
    }

    if(tables.count("FinalStateY")){
        readTable<double_value_field>(h5file,"FinalStateY", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.intY=p.get<CTS("value")>();
            });
    }

    if(tables.count("FinalType0")){
        readTable<double_value_field>(h5file,"FinalType0", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.final_state_particle_0=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    if(tables.count("FinalType1")){
        readTable<double_value_field>(h5file,"FinalType1", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.final_state_particle_1=static_cast<LW::ParticleType>(p.get<CTS("value")>());
            });
    }

    if (selectionStart!=-1) {
        if(tables.count("DeepStart")){
            readTable<double_value_field>(h5file,"DeepStart", intermediateData,
                [&](const double_value_field& p, Event& e){
                    double aux = (p.get<CTS("value")>());
                    e.topology = (aux>selectionStart) ? 0 : 1;
                });
        }
    }
    else {
        readTable<double_value_field>(h5file,"MuExZenith", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.topology = 3;
            });
    }

    if(tables.count("MuExZenith")){
        readTable<double_value_field>(h5file,"MuExZenith", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.zenith=(p.get<CTS("value")>());
            });
    }

    if(tables.count(energyName)){
        readTable<double_value_field>(h5file,energyName, intermediateData,
            [&](const double_value_field& p, Event& e){
                e.energy=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuAzimuth")){
        readTable<double_value_field>(h5file,"NuAzimuth", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryAzimuth=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuZenith")){
        readTable<double_value_field>(h5file,"NuZenith", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryZenith=(p.get<CTS("value")>());
            });
    }

    if(tables.count("NuEnergy")){
        readTable<double_value_field>(h5file,"NuEnergy", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.primaryEnergy=(p.get<CTS("value")>());
            });
    }

    if(tables.count("TotalColumnDepth")){
        readTable<double_value_field>(h5file,"TotalColumnDepth", intermediateData,
            [&](const double_value_field& p, Event& e){
                e.totalColumnDepth=(p.get<CTS("value")>());
            });
    }

    for(std::map<RecordID,Event>::value_type& item : intermediateData)
        action(item.first,item.second);

    intermediateData.clear();
}

} // close io sterile namespace

} // close h5load input namespace

#endif /* EVENT_H_ */
