/*
 * GollumParameters.h
 *
 *  Created on: Nov 21, 2017
 *      Author: aschneider
 */

#ifndef GOLLUMPARAMETERS_H_
#define GOLLUMPARAMETERS_H_

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <optional>
#include <sys/stat.h>
#include <boost/math/constants/constants.hpp>

#include "Event.h"
#include "GollumEnumDefinitions.h"
#include "utils.h"

/**
* @file GollumParameters.h
* @brief Definitions for useful data structures.
*/

namespace gollumfit{

/**
* 
* @struct FitParameters
* @brief Struct to store all the nuisance parameters for the fit.
*/
struct FitParameters {

  float convNorm;
  float promptNorm;
  float zenithCorrection;
  float kaonLosses;
  float hadronicHEkp;
  float hadronicHEkm;
  float hadronicVHE1pip;
  float hadronicVHE1pim;
  float hadronicVHE3kp;
  float hadronicVHE3km;
  float hadronicVHE3pip;
  float hadronicVHE3pim;
  float hadronicVHE3p;
  float hadronicVHE3n;
  float cosmicRay1;
  float cosmicRay2;
  float cosmicRay3;
  float cosmicRay4;
  float cosmicRay5;
  float cosmicRay6;
  float icegrad0;
  float icegrad1;
  float icegrad2;
  float icegrad3;
  float icegrad4;
  float icegrad5;
  float icegrad6;
  float icegrad7;
  float icegrad8;
  float domEfficiency;
  float holeiceForward;
  float astroNorm;
  float astroDeltaGamma;
  float astroDeltaGammaSec;
  float astroPivot;
  float NeutrinoAntineutrinoRatio;
  float nuxs;
  float nubarxs;

  /**
  * @brief Default constructor
  */
  FitParameters(){};

};

/**
* @brief Function to define the << operator to cout fit parameters.
*
* @param os std::ostream
* @param fp FitParameters object to be displayed
* 
*/
std::ostream& operator<<(std::ostream& os, const FitParameters& fp);

// true means its fix in the minimization. CA.
// defaults corresponds to the ClassicHESE paramters
/**
* 
* @struct FitParametersFlag
* @brief Struct to store a bool for each nuisance parameter
*
* True means it is fixed in the minimization (not minimized over)
*/
struct FitParametersFlag {

  bool convNorm = true;
  bool promptNorm = true;
  bool zenithCorrection = true;
  bool kaonLosses = true;
  bool hadronicHEkp = true;
  bool hadronicHEkm = true;
  bool hadronicVHE1pip = true;
  bool hadronicVHE1pim = true;
  bool hadronicVHE3kp = true;
  bool hadronicVHE3km = true;
  bool hadronicVHE3pip = true;
  bool hadronicVHE3pim = true;
  bool hadronicVHE3p = true;
  bool hadronicVHE3n = true;
  bool cosmicRay1 = true;
  bool cosmicRay2 = true;
  bool cosmicRay3 = true;
  bool cosmicRay4 = true;
  bool cosmicRay5 = true;
  bool cosmicRay6 = true;
  bool icegrad0 = true;
  bool icegrad1 = true;
  bool icegrad2 = true;
  bool icegrad3 = true;
  bool icegrad4 = true;
  bool icegrad5 = true;
  bool icegrad6 = true;
  bool icegrad7 = true;
  bool icegrad8 = true;
  bool domEfficiency = true;
  bool holeiceForward = true;
  bool astroNorm = true;
  bool astroDeltaGamma = true;
  bool astroDeltaGammaSec = true;
  bool astroPivot = true;
  bool NeutrinoAntineutrinoRatio = true;
  bool nuxs = true;
  bool nubarxs = true;

  /**
  * @brief Constructor
  */
  FitParametersFlag(){};

  /**
  * @brief Constructor
  *
  * @param flag bool assigned to each nuisance parameter
  */
  FitParametersFlag(bool flag) {
    convNorm = flag;
    promptNorm = flag;
    zenithCorrection = flag;
    kaonLosses = flag;
    hadronicHEkp = flag;
    hadronicHEkm = flag;
    hadronicVHE1pip = flag;
    hadronicVHE1pim = flag;
    hadronicVHE3kp = flag;
    hadronicVHE3km = flag;
    hadronicVHE3pip = flag;
    hadronicVHE3pim = flag;
    hadronicVHE3p = flag;
    hadronicVHE3n = flag;
    cosmicRay1 = flag;
    cosmicRay2 = flag;
    cosmicRay3 = flag;
    cosmicRay4 = flag;
    cosmicRay5 = flag;
    cosmicRay6 = flag;
    icegrad0 = flag;
    icegrad1 = flag;
    icegrad2 = flag;
    icegrad3 = flag;
    icegrad4 = flag;
    icegrad5 = flag;
    icegrad6 = flag;
    icegrad7 = flag;
    icegrad8 = flag;
    domEfficiency = flag;
    holeiceForward = flag;
    astroNorm = flag;
    astroDeltaGamma = flag;
    astroDeltaGammaSec = flag;
    astroPivot = flag;
    NeutrinoAntineutrinoRatio = flag;
    nuxs = flag;
    nubarxs = flag;    
  }

};

/**
* @struct FitParametersBound
*
* @brief Struct to store min and max bounds for each nuisance parameter
*/
struct FitParametersBound {

  float convNormMin;
  float convNormMax;
  float promptNormMin;
  float promptNormMax;
  float zenithCorrectionMin;
  float zenithCorrectionMax;
  float kaonLossesMin;
  float kaonLossesMax;
  float hadronicHEkpMin;
  float hadronicHEkpMax;
  float hadronicHEkmMin;
  float hadronicHEkmMax;
  float hadronicVHE1pipMin;
  float hadronicVHE1pipMax;
  float hadronicVHE1pimMin;
  float hadronicVHE1pimMax;
  float hadronicVHE3kpMin;
  float hadronicVHE3kpMax;
  float hadronicVHE3kmMin;
  float hadronicVHE3kmMax;
  float hadronicVHE3pipMin;
  float hadronicVHE3pipMax;
  float hadronicVHE3pimMin;
  float hadronicVHE3pimMax;
  float hadronicVHE3pMin;
  float hadronicVHE3pMax;
  float hadronicVHE3nMin;
  float hadronicVHE3nMax;
  float cosmicRay1Min;
  float cosmicRay1Max;
  float cosmicRay2Min;
  float cosmicRay2Max;
  float cosmicRay3Min;
  float cosmicRay3Max;
  float cosmicRay4Min;
  float cosmicRay4Max;
  float cosmicRay5Min;
  float cosmicRay5Max;
  float cosmicRay6Min;
  float cosmicRay6Max;
  float icegrad0Min;
  float icegrad0Max;
  float icegrad1Min;
  float icegrad1Max;
  float icegrad2Min;
  float icegrad2Max;
  float icegrad3Min;
  float icegrad3Max;
  float icegrad4Min;
  float icegrad4Max;
  float icegrad5Min;
  float icegrad5Max;
  float icegrad6Min;
  float icegrad6Max;
  float icegrad7Min;
  float icegrad7Max;
  float icegrad8Min;
  float icegrad8Max;
  float domEfficiencyMin;
  float domEfficiencyMax;
  float holeiceForwardMin;
  float holeiceForwardMax;
  float astroNormMin;
  float astroNormMax;
  float astroDeltaGammaMin;
  float astroDeltaGammaMax;
  float astroDeltaGammaSecMin;
  float astroDeltaGammaSecMax;
  float astroPivotMin;
  float astroPivotMax;
  float NeutrinoAntineutrinoRatioMin;
  float NeutrinoAntineutrinoRatioMax;
  float nuxsMin;
  float nuxsMax;
  float nubarxsMin;
  float nubarxsMax;

  /**
  * @brief Default constructor
  */
  FitParametersBound(){};

};


/**
* @struct Priors
* 
* @brief Struct to store prior for each nuisance parameter
*
* Each parameter is defined by a center and a width.
* Option to define a 16 by 16 correlation matrix for flux parameters and
* a 9 by 9 correlation matrix for ice parameters.
*/
struct Priors {

  float convNormWidth; 
  float convNormCenter;
  float promptNormWidth; 
  float promptNormCenter;
  float zenithCorrectionWidth; 
  float zenithCorrectionCenter;
  float kaonLossesWidth; 
  float kaonLossesCenter;
  float hadronicHEkpWidth; 
  float hadronicHEkpCenter;
  float hadronicHEkmWidth; 
  float hadronicHEkmCenter;
  float hadronicVHE1pipWidth; 
  float hadronicVHE1pipCenter;
  float hadronicVHE1pimWidth; 
  float hadronicVHE1pimCenter;
  float hadronicVHE3kpWidth; 
  float hadronicVHE3kpCenter;
  float hadronicVHE3kmWidth; 
  float hadronicVHE3kmCenter;
  float hadronicVHE3pipWidth; 
  float hadronicVHE3pipCenter;
  float hadronicVHE3pimWidth; 
  float hadronicVHE3pimCenter;
  float hadronicVHE3pWidth; 
  float hadronicVHE3pCenter;
  float hadronicVHE3nWidth; 
  float hadronicVHE3nCenter;
  float cosmicRay1Width; 
  float cosmicRay1Center;
  float cosmicRay2Width; 
  float cosmicRay2Center;
  float cosmicRay3Width; 
  float cosmicRay3Center;
  float cosmicRay4Width; 
  float cosmicRay4Center;
  float cosmicRay5Width; 
  float cosmicRay5Center;
  float cosmicRay6Width; 
  float cosmicRay6Center;
  float icegrad0Width; 
  float icegrad0Center;
  float icegrad1Width; 
  float icegrad1Center;
  float icegrad2Width; 
  float icegrad2Center;
  float icegrad3Width; 
  float icegrad3Center;
  float icegrad4Width; 
  float icegrad4Center;
  float icegrad5Width; 
  float icegrad5Center;
  float icegrad6Width; 
  float icegrad6Center;
  float icegrad7Width; 
  float icegrad7Center;
  float icegrad8Width; 
  float icegrad8Center;
  float domEfficiencyWidth; 
  float domEfficiencyCenter;
  float holeiceForwardWidth; 
  float holeiceForwardCenter;
  float astroNormWidth; 
  float astroNormCenter;
  float astroDeltaGammaWidth; 
  float astroDeltaGammaCenter;
  float astroDeltaGammaSecWidth; 
  float astroDeltaGammaSecCenter;
  float astroPivotMin,astroPivotMax; 
  float NeutrinoAntineutrinoRatioWidth; 
  float NeutrinoAntineutrinoRatioCenter;
  float nuxsWidth; 
  float nuxsCenter;
  float nubarxsWidth; 
  float nubarxsCenter;

  float flux_corr[16][16];///< Flux parameters correlation matrix
  float ice_corr[9][9];///< Ice parameters correlation matrix

  /**
  * @brief Constructor
  */
  Priors(){}

  /**
  * @brief Constructor
  * 
  * @param fit_params Centers of the nuisance parameters
  */
  Priors(const FitParameters& fit_params){
    convNormCenter = fit_params.convNorm;
    promptNormCenter = fit_params.promptNorm;
    zenithCorrectionCenter = fit_params.zenithCorrection;
    kaonLossesCenter = fit_params.kaonLosses;
    hadronicHEkpCenter = fit_params.hadronicHEkp;
    hadronicHEkmCenter = fit_params.hadronicHEkm;
    hadronicVHE1pipCenter = fit_params.hadronicVHE1pip;
    hadronicVHE1pimCenter = fit_params.hadronicVHE1pim;
    hadronicVHE3kpCenter = fit_params.hadronicVHE3kp;
    hadronicVHE3kmCenter = fit_params.hadronicVHE3km;
    hadronicVHE3pipCenter = fit_params.hadronicVHE3pip;
    hadronicVHE3pimCenter = fit_params.hadronicVHE3pim;
    hadronicVHE3pCenter = fit_params.hadronicVHE3p;
    hadronicVHE3nCenter = fit_params.hadronicVHE3n;
    cosmicRay1Center = fit_params.cosmicRay1;
    cosmicRay2Center = fit_params.cosmicRay2;
    cosmicRay3Center = fit_params.cosmicRay3;
    cosmicRay4Center = fit_params.cosmicRay4;
    cosmicRay5Center = fit_params.cosmicRay5;
    cosmicRay6Center = fit_params.cosmicRay6;
    icegrad0Center = fit_params.icegrad0;
    icegrad1Center = fit_params.icegrad1;
    icegrad2Center = fit_params.icegrad2;
    icegrad3Center = fit_params.icegrad3;
    icegrad4Center = fit_params.icegrad4;
    icegrad5Center = fit_params.icegrad5;
    icegrad6Center = fit_params.icegrad6;
    icegrad7Center = fit_params.icegrad7;
    icegrad8Center = fit_params.icegrad8;
    domEfficiencyCenter = fit_params.domEfficiency;
    holeiceForwardCenter = fit_params.holeiceForward;
    astroNormCenter = fit_params.astroNorm;
    astroDeltaGammaCenter = fit_params.astroDeltaGamma;
    astroDeltaGammaSecCenter = fit_params.astroDeltaGammaSec;
    NeutrinoAntineutrinoRatioCenter = fit_params.NeutrinoAntineutrinoRatio;
    nuxsCenter = fit_params.nuxs;
    nubarxsCenter = fit_params.nubarxs;
  }
  
  /**
  * @brief Set Flux correlation matrix elements
  * 
  * @param i first index
  * @param j second index
  * @param value value to set element to
  */
  void SetFluxCorr(int i, int j, float value) { flux_corr[i][j] = value; }

  /**
  * @brief Set Ice correlation matrix elements
  * 
  * @param i first index
  * @param j second index
  * @param value value to set element to
  */
  void SetIceGradientsCorr(int i, int j, float value) { ice_corr[i][j] = value; }

};

/**
  * @struct FitResult
  * @brief Struct to contain result of fit
  * 
  */
struct FitResult {
  FitParameters params;
  double likelihood;
  double aux_likelihood; ///<TO DO
  unsigned int nEval, nGrad;
  bool succeeded;

  /**
  * @brief Constructor
  */
  FitResult(){};
};

/**
* @struct DataPaths
* @brief Holds file and directory paths used throughout the application.
*
* The DataPaths struct stores various file and directory paths as strings. It includes
* methods to check the validity of these paths to ensure they point to accessible
* directories or files as expected. If a path is found to be invalid or points to a non-directory
* item, an exception is thrown.
*/
struct DataPaths {
  private:
   
   /**
    * @brief Checks the validity of a data path.
    *
    * This method verifies that a given path points to an accessible directory.
    * If the directory does not exist or cannot be accessed, a runtime error is thrown.
    * Warnings are issued for empty paths.
    *
    * @param p The file system path to be checked.
    * @return True if the path is a directory and is accessible, false otherwise.
    * @throws std::runtime_error if the path cannot be accessed or is not a directory.
    */
   bool CheckDataPath(std::string p) const {
     struct stat info;
     bool status=true;
     if(p!="")
       {
         if( stat(p.c_str(), &info) != 0 )
         {
             status=true;
             throw std::runtime_error("cannot access "+ p);
         }
         else if( !(info.st_mode & S_IFDIR) )
         {
             status=false;
             throw std::runtime_error("is not a directory: " +p);
         }
     } else {
       std::cout<<"Warning, there are unset paths in DataPaths. Check you want this."<<std::endl;
       return false;
     }
     return status;
    }
  public:
    std::string compact_file_path;///< path to the FastMC template
    std::string neutrino_cc_xs_spline_path;
    std::string antineutrino_cc_xs_spline_path;
    std::string neutrino_nc_xs_spline_path;
    std::string antineutrino_nc_xs_spline_path;
    std::string diff_neutrino_cc_xs_spline_path;
    std::string diff_antineutrino_cc_xs_spline_path;
    std::string diff_neutrino_nc_xs_spline_path;
    std::string diff_antineutrino_nc_xs_spline_path;
    std::string data_path;
    std::string mc_path;
    std::string conventional_nusquids_atmospheric_file;
    std::string prompt_nusquids_atmospheric_file;
    std::string astro_nusquids_file;
    std::string domeff_spline_path;
    std::string holeice_spline_path;
    std::string attenuation_spline_path;
    std::string ice_gradient_spline_path;
    std::string atmospheric_density_spline_path;
    std::string atmospheric_kaonlosses_spline_path;
    std::string hadronic_spline_path;
    std::string cosmic_ray_spline_path;

    /**
    * @brief Default constructor.
    */
    DataPaths(){};

};

/**
* @struct SteeringParams
* @brief Holds parameters used for steering simulations or analyses.
*
* This structure contains configuration parameters that affect how
* simulations or analyses are steered. It includes parameters such as
* energy range, angular range, bin widths, file names, multithreading options,
* and tolerances for numerical computations.
*
* @todo add detailed descriptions for the more technical parameters, i.e. change_tol, etc. 
*/
struct SteeringParams {
  string readCompact="";
  double minFitEnergy;///< minimum energy
  double maxFitEnergy;///< maximum energy
  double minCosth;///< minimum cosine zenith
  double maxCosth;///< maximum cosine zenith
  double logEbinEdge;///< starting bin edge in log energy
  double logEbinWidth;///< bin width in log energy
  double cosThbinEdge;///< starting bin edge in cosine zenith
  double cosThbinWidth;///< bin width in cosine zenith
  std::optional<double> astroOscAvgScale = std::nullopt;///< energy scale above which the astro component will be averaged out in the nuSQuIDSAtm class

  std::vector<std::string> ice_gradient_filename;
  std::vector<std::string> active_hadronic_parameters;
  std::vector<std::string> active_cosmicray_parameters;

  size_t evalThreads=4;

  double fullLivetime;

  double change_tol = 1.e-20;
  double grad_tol = 1.e-20;
  double uncertaintyModSigmaOverMu = 0.;

  string model_label;
  string simToLoad;
  string energyName;
  double selectionStart;

  bool enableTotalNorm = true;///< If true, convNorm is applied to all flux components

  /**
  * @brief Default constructor.
  */
  SteeringParams(){};

};


} // namespace gollumfit

#endif /* GOLLUMPARAMETERS_H_ */
