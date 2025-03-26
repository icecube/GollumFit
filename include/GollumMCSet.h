#ifndef GOLLUMMCSET_H_
#define GOLLUMMCSET_H_

#include <LeptonWeighter/Event.h>
#include <LeptonWeighter/ParticleType.h>
#include <LeptonWeighter/Generator.h>

#include "GollumEnumDefinitions.h"

/**
* @file GollumMCSet.h
* @brief MCSet struct definition
*/

namespace gollumfit {


/**
* @brief Monte Carlo configuration set for simulations.
*
* This struct represents a set of configurations for running Monte Carlo simulations, including
* paths, filenames, splitting options, and physical parameters related to the simulation.
*/
struct MCSet {

  /// Path to the directory containing simulation files.
	const std::string path;

  /// Name of the file to be used in the simulation.
  const std::string filename;

  /**
   * A pair indicating whether the simulation should be split and, if so, the number of parts to split into.
   * The boolean value indicates the presence of splitting (true for split, false for no split), while the
   * unsigned int represents the number of splits.
  */
  const std::pair<bool, unsigned int> split;

  /// The total number of files to be considered in the simulation.
  const double number_of_files;

  /// The fraction of unshadowed area, also known as dom efficiency.
  const double unshadowedFraction ;// aka domefficiency

  /// Forward scattering parameter for hole ice, also known as p2.
  const double holeiceForward;// aka p2

  /// List of generator pointers used for the simulation.
  std::vector<std::shared_ptr<LW::Generator>> generators;

  /**
  * @brief Construct a new MCSet object with specified parameters.
  *
  * @param p Path to the simulation files.
  * @param filename Name of the simulation file.
  * @param split Splitting configuration as a pair.
  * @param number_of_files Number of files in the simulation.
  * @param unshadowedFraction Dome efficiency or unshadowed fraction.
  * @param holeiceForward Hole ice forward scattering parameter.
  * @param g Vector of shared pointers to generators.
  */
	MCSet(const std::string& p,const std::string& filename, const std::pair<bool, unsigned int> split, double number_of_files, double unshadowedFraction, double holeiceForward, std::vector<std::shared_ptr<LW::Generator>> g):
	path(p),filename(filename),split(split),number_of_files(number_of_files),unshadowedFraction(unshadowedFraction),holeiceForward(holeiceForward),generators(g){}
  
  /**
  * @brief Copy constructor for MCSet object.
  *
  * @param other The MCSet object to copy from.
  */
	MCSet(MCSet& other):path(other.path),
  filename(other.filename),split(other.split),number_of_files(other.number_of_files),
  unshadowedFraction(other.unshadowedFraction),holeiceForward(other.holeiceForward),generators(other.generators){}
  
  // copy constructor
	MCSet(const MCSet& other):path(other.path),
  filename(other.filename),split(other.split),number_of_files(other.number_of_files),
  unshadowedFraction(other.unshadowedFraction),holeiceForward(other.holeiceForward),generators(other.generators){}
  
  /**
  * @brief Move constructor for MCSet object.
  *
  * @param other The MCSet object to move from.
  */
	MCSet(MCSet&& other):path(std::move(other.path)),
  filename(std::move(other.filename)),split(std::move(other.split)),number_of_files(std::move(other.number_of_files)),
  unshadowedFraction(std::move(other.unshadowedFraction)),holeiceForward(std::move(other.holeiceForward)),generators(std::move(other.generators)){}
};

} // close gollumfit namespace

#endif /* GOLLUMMCSET_H_ */
