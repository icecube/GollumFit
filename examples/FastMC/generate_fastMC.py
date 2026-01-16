#!/usr/bin/env python3
"""
FastMC Generation Example

This script demonstrates how to generate FastMC, which is a compressed
representation of Monte Carlo that significantly improves fitting efficiency.

Example command:
    nohup python generate_fastMC.py > generate_fastMC.log 2>&1 &
"""

import GollumFitPy as gf
import numpy as np
import os

#####################################################################################
# Configure Data Paths - Set paths for cross section splines
#####################################################################################
datapaths = gf.DataPaths()
datapaths.neutrino_cc_xs_spline_path             = "../../resources/Splines/CrossSections/sigma_nu_CC_iso.fits"
datapaths.antineutrino_cc_xs_spline_path         = "../../resources/Splines/CrossSections/sigma_nubar_CC_iso.fits"
datapaths.neutrino_nc_xs_spline_path             = "../../resources/Splines/CrossSections/sigma_nu_NC_iso.fits"
datapaths.antineutrino_nc_xs_spline_path         = "../../resources/Splines/CrossSections/sigma_nubar_NC_iso.fits"
datapaths.diff_neutrino_cc_xs_spline_path        = "../../resources/Splines/CrossSections/dsdxdy_nu_CC_iso.fits"
datapaths.diff_antineutrino_cc_xs_spline_path    = "../../resources/Splines/CrossSections/dsdxdy_nubar_CC_iso.fits"
datapaths.diff_neutrino_nc_xs_spline_path        = "../../resources/Splines/CrossSections/dsdxdy_nu_NC_iso.fits"
datapaths.diff_antineutrino_nc_xs_spline_path    = "../../resources/Splines/CrossSections/dsdxdy_nubar_CC_iso.fits"
datapaths.mc_path                                = "../../monte_carlo/"
datapaths.domeff_spline_path                     = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path                    = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path                = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.ice_gradient_spline_path               = "../../resources/Splines/IceGradientsSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.atmospheric_density_spline_path        = "../../resources/Splines/AtmosphericZenithVariationSplines/atm_density_1s.fits"
datapaths.atmospheric_kaonlosses_spline_path     = "../../resources/Splines/AtmosphericKaonLossesSplines/kaon_loses_1s.fits"

#####################################################################################
# Configure Flux Files - Load atmospheric, prompt, and astrophysical flux files
#####################################################################################
datapaths.conventional_nusquids_atmospheric_file = "../fluxes/atmospheric.hdf5"
datapaths.prompt_nusquids_atmospheric_file       = "../fluxes/prompt_atmospheric.hdf5"
datapaths.astro_nusquids_file                    = "../fluxes/astro.hdf5"

# Hadronic and cosmic ray correction splines (necessary for flux nuisance parameters)
hadronlist = ["he_K+", "he_K-", "vhe1_pi+", "vhe1_pi-", "vhe3_K+", "vhe3_K-", 
              "vhe3_pi+", "vhe3_pi-", "vhe3_p", "vhe3_n"]
crlist = ["GSF_1", "GSF_2", "GSF_3", "GSF_4", "GSF_5", "GSF_6"]

datapaths.hadronic_spline_path   = "../fluxes"
datapaths.cosmic_ray_spline_path = "../fluxes"

#####################################################################################
# Set Steering Parameters - Configure analysis binning and settings
# NOTE: Binning choices affect FastMC compression and must match analysis configuration
#####################################################################################
steering_params = gf.SteeringParams()
steering_params.minFitEnergy                    = 300
steering_params.maxFitEnergy                    = 1e5
steering_params.logEbinEdge                     = np.log10(300)
steering_params.logEbinWidth                    = (np.log10(1e5) - np.log10(300)) / 24
steering_params.minCosth                        = -1.0
steering_params.maxCosth                        = 0.0
steering_params.cosThbinEdge                    = 0.0
steering_params.cosThbinWidth                   = 0.05
steering_params.selectionStart                  = 0.99
steering_params.ice_gradient_filename           = ["Amp_0", "Amp_1", "Amp_2", "Amp_3", "Amp_4", 
                                                   "Phs_1", "Phs_2", "Phs_3", "Phs_4"]
steering_params.active_hadronic_parameters      = hadronlist
steering_params.active_cosmicray_parameters     = crlist

# Livetime for the corresponding Monte Carlo
years = 10.669
steering_params.fullLivetime                    = years * 365 * 24 * 60 * 60.
steering_params.simToLoad                       = "BDT_Split_HE"
steering_params.energyName                      = "DnnEnergy"
steering_params.model_label                     = ""  # Can be used for uniquely-labelled flux files

#####################################################################################
# Construct and Write FastMC
#####################################################################################
gollumfit = gf.GollumFit(datapaths, steering_params)

# Compression parameter: smaller values = higher compression but potential accuracy loss
metascaling = 0.25
gollumfit.ConstructFastMode(metascaling)

# Write to file (output will be STERILE.meows in the specified directory)
gollumfit.WriteCompact("example.fastmc")

print("Done generating FastMC.")





