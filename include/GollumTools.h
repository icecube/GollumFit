#ifndef _GOLLUMTOOLS_H_
#define _GOLLUMTOOLS_H_

#include <functional>
#include <iterator>
#include <set>
#include <string>
#include <chrono>
#include <queue>
#include <vector>
#include <memory>
#include <fstream>
#include <LeptonWeighter/ParticleType.h>
#include <nuSQuIDS/xsections.h>

/**
* @file GollumTools.cpp
* @brief Helper function declarations and ScaledNeutrinoCrossSections class.
*/

namespace gollumfit {
  namespace tools {
    using namespace nusquids;

    /**
    * @brief Check if the given particle type is a neutrino.
    *
    * @param p The particle type to check.
    * @return True if the particle type is an antineutrino, false otherwise.
    */
    bool isNeutrino(LW::ParticleType p);

    /**
    * @brief Check if the given particle type is an antineutrino.
    *
    * @param p The particle type to check.
    * @return True if the particle type is an antineutrino, false otherwise.
    */
    bool isAntineutrino(LW::ParticleType p);

    /**
    * @brief Checks if the file at the given file path exists and is readable.
    *
    * This function attempts to open a file at the specified file path to determine its existence
    * and readability. If the file cannot be opened or does not exist, the function throws a
    * runtime error. An optional message can be printed to standard output indicating the file
    * being read.
    *
    * @param FilePath The path to the file that needs to be checked.
    * @param quiet If set to false, the function prints the file reading action to standard output.
    *
    * @return The original file path is returned if the file exists and is readable.
    *
    * @throws std::runtime_error If the file at the given path does not exist or cannot be opened.
    */
    std::string CheckedFilePath(std::string FilePath, bool quiet = false);

    /**
    * @class ScaledNeutrinoCrossSections
    * @brief A class to represent scaled neutrino cross sections.
    *
    * Inherits from NeutrinoDISCrossSectionsFromTables and applies energy-dependent scale factors to 
    * the deep inelastic scattering (DIS) cross sections for neutrinos. These scale factors are 
    * applied to cross sections read from tables or provided vectors.
    */
    class ScaledNeutrinoCrossSections : public NeutrinoDISCrossSectionsFromTables {
    private:
      std::vector<double> edges;
      std::vector<double> scales;

    public:
      /// \brief Destructor.
      virtual ~ScaledNeutrinoCrossSections(){} 
      
      /// \brief Default constructor with built-in tables
      ScaledNeutrinoCrossSections()
        : edges(0), scales(0) {}
      /**
      * @brief Constructor to initialize from tables at a given path.
      * 
      * @param path The file path to the cross-section table data.
      */
      ScaledNeutrinoCrossSections(std::string path)
        : NeutrinoDISCrossSectionsFromTables(path), edges(0), scales(0) {}
      /**
      * @brief Constructor to initialize with provided energy bin edges and scale factors.
      * 
      * @param e Vector of energy bin edges.
      * @param s Vector of scale factors corresponding to the energy bins.
      */
      ScaledNeutrinoCrossSections(std::vector<double> e, std::vector<double> s)
        : edges(e), scales(s) {}
      /**
      * @brief Constructor to initialize from tables at a given path with provided energy bin edges and scale factors.
      * 
      * @param path The file path to the cross-section table data.
      * @param e Vector of energy bin edges.
      * @param s Vector of scale factors corresponding to the energy bins.
      */
      ScaledNeutrinoCrossSections(std::string path, std::vector<double> e, std::vector<double> s)
        : NeutrinoDISCrossSectionsFromTables(path), edges(e), scales(s) {}

      /**
      * @brief Determines the scale factor to apply to cross sections for a neutrino with a given energy.
      * 
      * @param Enu The energy of the neutrino.
      * @return The scale factor for the given neutrino energy.
      */
      double GetScale(double Enu) const {
        unsigned int nedges = edges.size();
        if (nedges == 0) return 1.;
        if (Enu < edges[0]) return 1.;
        if (Enu > edges[nedges-1]) return 1.;
        for (unsigned int i=1;i<nedges;i++) {
          if (Enu < edges[i]) {
            return scales[i-1];
          }
        }
        return std::numeric_limits<double>::quiet_NaN();
      }

      /**
      * @brief Returns the total neutrino cross section scaled by the energy-dependent scale factor.
      * 
      * @param Enu The energy of the neutrino.
      * @param flavor The flavor of the neutrino.
      * @param neutype The type of the neutrino (neutrino or antineutrino).
      * @param current The current type (charged or neutral).
      * @return The scaled total cross section.
      */
      double TotalCrossSection(double Enu, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override {
        return GetScale(Enu)*NeutrinoDISCrossSectionsFromTables::TotalCrossSection(Enu,flavor,neutype,current);
      }

      /**
        * @brief Returns the Differential cross section with respect to the outgoing lepton energy.
        * 
        * @param E1 Incident lepton energy.
        * @param E2 Outgoing lepton energy.
        * @param flavor Flavor index.
        * @param neutype Can be either neutrino or antineutrino.
        * @param current Can be either CC or NC.
        * 
        * @return The cross section in cm^2 GeV^-1.
      */
      double SingleDifferentialCrossSection(double E1, double E2, NeutrinoFlavor flavor, NeutrinoType neutype, Current current) const override {
        return GetScale(E1)*NeutrinoDISCrossSectionsFromTables::SingleDifferentialCrossSection(E1,E2,flavor,neutype,current);
      }
    };

  }// closing tools namespace
}// close gollumfit namespace

#endif /* GOLLUMTOOLS_H_ */
