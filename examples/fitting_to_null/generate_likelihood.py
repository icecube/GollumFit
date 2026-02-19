#!/usr/bin/env python3

"""
Likelihood Minimization Example - Fitting to Null Hypothesis

This script demonstrates GollumFit's core functionality: fitting Monte Carlo to
data using likelihood minimization. We fit to null pseudo-data starting from
random initial parameter values.

Example command:
    nohup time python ./generate_likelihood.py > generate_likelihood.log 2>&1 &
"""

import GollumFitPy as gf
import numpy as np
import os
import sys
import subprocess
import scipy.stats as stats

print('Starting generate_likelihood.py.')

#####################################################################################
# Define Nuisance Parameters (All set to vary in fit with False flag)
# Format: [vary_flag, prior_type, center, width, lower_bound, upper_bound]
#####################################################################################
syst_dict     = { 
    'convNorm'                  : [ False, 'Gaussian',      1.,   0.2,                   0.1,                   3. ], 
    'zenithCorrection'          : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'kaonLosses'                : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'hadronicHEkp'              : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicHEkm'              : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE1pip'           : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE1pim'           : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3kp'            : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3km'            : [ False, 'Gaussian',      0.,    1.,                  -1.5,                   2. ], 
    'hadronicVHE3pip'           : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3pim'           : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3p'             : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3n'             : [ False, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'cosmicRay1'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay2'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay3'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay4'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay5'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay6'                : [ False, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'icegrad0'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad1'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad2'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad3'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad4'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad5'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad6'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad7'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad8'                  : [ False, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'domEfficiency'             : [ False, 'Gaussian',    1.27, 0.123,                 1.234,                1.346 ], 
    'holeiceForward'            : [ False, 'Gaussian',     -1.,   10.,                 -5.35,                 1.85 ], 
    'astroNorm'                 : [ False, 'Gaussian', 4.72/6.,  0.36,                    0.,                   3. ], 
    'astroDeltaGamma'           : [ False, 'Gaussian',      0.,  0.36,                   -2.,                   2. ], 
    'astroDeltaGammaSec'        : [ False, 'Gaussian',      0.,  0.36,                   -2.,                   2. ], 
    'nuxs'                      : [ False, 'Gaussian',      1.,   0.1,                 0.824,                1.176 ], 
    'nubarxs'                   : [ False, 'Gaussian',      1.,   0.1,                 0.824,                1.176 ], 
    'astroPivot'                : [ False,  'Uniform',      5.,    1.,                    4.,                   6. ], 
    'promptNorm'                : [ False, 'Gaussian',      1.,    1.,                    0.,                   3. ],
    'NeutrinoAntineutrinoRatio' : [ False, 'Gaussian',      1.,    1.,                    0.,                   2. ],
}

#####################################################################################
# Define Random Sampling Function
# Helper function to sample random parameter values from prior distributions
#####################################################################################
def throw(syst):
    """Sample random value from prior distribution."""
    if syst[1] == 'Gaussian':
        # Truncated normal distribution
        val = stats.truncnorm(
            (syst[4] - syst[2]) / syst[3],  # Lower bound (standardized)
            (syst[5] - syst[2]) / syst[3],  # Upper bound (standardized)
            syst[2],  # Mean
            syst[3]   # Std dev
        ).rvs(1)
        return val[0]
    else:
        # Uniform distribution
        return np.random.uniform(syst[4], syst[5])

#####################################################################################
# Initialize Fit Parameter Objects
# Create objects to manage fit configuration
#####################################################################################

fitparams_flag  = gf.FitParametersFlag()  # Which parameters to vary
fitparams_bound = gf.FitParametersBound()  # Parameter bounds
priors          = gf.Priors()              # Prior distributions
seed_fitparams  = gf.FitParameters()       # Initial values

#####################################################################################
# Set Priors and Random Initial Values
# Configure priors and randomly initialize starting parameter values
#####################################################################################
np.random.seed(100)  # For reproducibility
print('Initializing with the following randomly-seeded nuisance params:')

for sname in syst_dict.keys():
    # Set flags and bounds
    exec(f'fitparams_flag.{sname} = syst_dict["{sname}"][0]')
    exec(f'fitparams_bound.{sname}Min = syst_dict["{sname}"][4]')
    exec(f'fitparams_bound.{sname}Max = syst_dict["{sname}"][5]')
    
    # Set priors
    if syst_dict[sname][1] == 'Gaussian':
        exec(f'priors.{sname}Center = syst_dict["{sname}"][2]')
        exec(f'priors.{sname}Width  = syst_dict["{sname}"][3]')
    else:
        exec(f'priors.{sname}Min = syst_dict["{sname}"][4]')
        exec(f'priors.{sname}Max = syst_dict["{sname}"][5]')
    
    # Randomly initialize
    thrown_val = throw(syst_dict[sname])
    exec(f'seed_fitparams.{sname} = thrown_val')
    print(f'{sname}: {thrown_val}')

#####################################################################################
# Set Correlations and Paths
# Load correlation matrices for ice gradients and flux parameters
#####################################################################################
gollumdir = "../../gollumfit-release"

# Set correlations (required for fitting/minimization, not for likelihood evaluation)
iceg_corr = np.load('../../resources/correlation_matrices/icegrad_correlations.npy')
flux_corr = np.load('../../resources/correlation_matrices/flux_correlations_new_ddmnodeis.npy')
for idx, val in np.ndenumerate(iceg_corr):
    priors.SetIceGradientsCorr(idx[0], idx[1], val)
for idx, val in np.ndenumerate(flux_corr):
    priors.SetFluxCorr(idx[0], idx[1], val)

datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC/compact.fastmc"

#####################################################################################
# Configure Steering Parameters
# Set binning and convergence criteria (must match FastMC binning)
#####################################################################################
edges = np.logspace(np.log10(300), np.log10(1e5), 25)
steering_params                = gf.SteeringParams()
steering_params.minFitEnergy   = edges[0]
steering_params.maxFitEnergy   = edges[-1]
steering_params.logEbinEdge    = np.log10(edges[0])
steering_params.logEbinWidth   = np.log10(edges[1]) - np.log10(edges[0])
steering_params.minCosth       = -1.0
steering_params.maxCosth       = 0.0
steering_params.cosThbinEdge   = 0.0
steering_params.cosThbinWidth  = 0.05
steering_params.selectionStart = float("DnnEnergy_0.99".split("_")[1])
steering_params.evalThreads    = 1

# Convergence criteria (tight tolerances for accurate minimization)
steering_params.change_tol     = 1.e-20
steering_params.grad_tol       = 1.e-20
steering_params.uncertaintyModSigmaOverMu = 0.0

#####################################################################################
# Load Data and Configure Fit
# Create GollumFit object and load pseudo-data
#####################################################################################
gollumfit = gf.GollumFit(datapaths, steering_params)

#####################################################################################
# declare the fake data location and load it
#####################################################################################
realization  = "nullexpectation.npz"
total_data = gollumfit.SetData(np.load(realization)["realization"])

#####################################################################################
# feed the flags, bounds, priors, on the nuisance parameters into gollumfit
#####################################################################################
gollumfit.SetFitParametersFlag(fitparams_flag)
gollumfit.SetFitParametersBound(fitparams_bound)
gollumfit.SetFitParametersPriors(priors)
gollumfit.SetFitParametersSeed([seed_fitparams])
gollumfit.ConstructLikelihoodProblem()

#####################################################################################
# perform the minimization
#####################################################################################
print("Starting minimization...")
min_llh = gollumfit.MinLLH()

#####################################################################################
# results: print the best fit nuisance parameters, likelihood, and the number of LLH evaluations
#####################################################################################
systematics = ""
for sname in syst_dict.keys() :
    exec('print(\"'+sname+'\",min_llh.params.'+sname+')')
    exec('systematics += str(min_llh.params.'+sname+')+\" \"')

print('llh:',min_llh.likelihood)
print('nEval: '+str(min_llh.nEval))

print("Completed successfully. Bye!")
