#!/usr/bin/env python3

# example command
# nohup python generate_fastMC.py  > generate_fastMC.log 2>&1 &

import GollumFitPy as gf
import numpy as np
import os

#####################################################################################
# Set the paths for necessary spline files
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
# Set the paths for the flux files, and for the hadronic and cosmic ray splines
#####################################################################################
datapaths.conventional_nusquids_atmospheric_file = "../fluxes/atmospheric.hdf5"
datapaths.prompt_nusquids_atmospheric_file       = "../fluxes/prompt_atmospheric.hdf5"
datapaths.astro_nusquids_file                    = "../fluxes/astro.hdf5"

hadronlist = ["he_K+","he_K-","vhe1_pi+","vhe1_pi-","vhe3_K+","vhe3_K-","vhe3_pi+","vhe3_pi-","vhe3_p","vhe3_n"]
crlist     = ["GSF_1","GSF_2","GSF_3","GSF_4","GSF_5","GSF_6"]

datapaths.hadronic_spline_path   = "../fluxes"
datapaths.cosmic_ray_spline_path = "../fluxes"

#####################################################################################
# steering params to set the binning 
#####################################################################################
steering_params = gf.SteeringParams()
steering_params.minFitEnergy                    = 300
steering_params.maxFitEnergy                    = 1e5
steering_params.logEbinEdge                     = np.log10(300)
steering_params.logEbinWidth                    = (np.log10(1e5)-np.log10(300))/24
steering_params.minCosth                        = -1
steering_params.maxCosth                        = 0.
steering_params.cosThbinEdge                    = 0.0
steering_params.cosThbinWidth                   = 0.05
steering_params.selectionStart                  = float("0.99")
steering_params.ice_gradient_filename           = ["Amp_0","Amp_1","Amp_2","Amp_3","Amp_4","Phs_1","Phs_2","Phs_3","Phs_4"]
steering_params.active_hadronic_parameters      = hadronlist
steering_params.active_cosmicray_parameters     = crlist
years = 10.669 #livetime for the corresponding MC
steering_params.fullLivetime                    = float(years)*365*24*60*60.
steering_params.simToLoad                       = "BDT_Split_HE"
steering_params.energyName                      = "DnnEnergy"
steering_params.model_label                     = "" # empty here, but can be used to call uniquely-labelled flux files. e.g. "GSF_1_model_label.hdf5"

#####################################################################################
# declare gollumfit object
#####################################################################################
gollumfit = gf.GollumFit(datapaths,steering_params)

#####################################################################################
# make the FastMC
#####################################################################################
# metascaling controls the FastMC compression. Large values can cause loss of accuracy
metascaling = "0.25"
gollumfit.ConstructFastMode(float(metascaling))

#####################################################################################
# Save it
#####################################################################################
gollumfit.WriteCompact(".")

print("Done generating FastMC.")





