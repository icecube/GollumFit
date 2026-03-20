"""
Generate Null Pseudo-Data Example

This script generates null pseudo-data (expectation with nominal nuisance parameters)
for use in fitting examples.

Example command:
    nohup time python ./generate_fakedata.py > generate_fakedata.log 2>&1 &
"""

import GollumFitPy as gf
import numpy as np
import os
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats

#####################################################################################
# Define Nuisance Parameters
# Format: [vary_flag, prior_type, center, width, lower_bound, upper_bound]
#####################################################################################
syst_dict     = OrderedDict({ 
    'convNorm'                  : [ True, 'Gaussian',      1.,   0.2,                   0.1,                   3. ], 
    'zenithCorrection'          : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'kaonLosses'                : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'hadronicHEkp'              : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicHEkm'              : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE1pip'           : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE1pim'           : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3kp'            : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3km'            : [ True, 'Gaussian',      0.,    1.,                  -1.5,                   2. ], 
    'hadronicVHE3pip'           : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3pim'           : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3p'             : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'hadronicVHE3n'             : [ True, 'Gaussian',      0.,    1.,                   -2.,                   2. ], 
    'cosmicRay1'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay2'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay3'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay4'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay5'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'cosmicRay6'                : [ True, 'Gaussian',      0.,    1.,                   -4.,                   4. ], 
    'icegrad0'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad1'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad2'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad3'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad4'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad5'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad6'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad7'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'icegrad8'                  : [ True, 'Gaussian',      0.,    1.,                   -3.,                   3. ], 
    'domEfficiency'             : [ True, 'Gaussian',    1.27, 0.123,                 1.234,                1.346 ], 
    'holeiceForward'            : [ True, 'Gaussian',     -1.,   10.,                 -5.35,                 1.85 ], 
    'astroNorm'                 : [ True, 'Gaussian', 4.72/6.,  0.36,                    0.,                   3. ], 
    'astroDeltaGamma'           : [ True, 'Gaussian',      0.,  0.36,                   -2.,                   2. ], 
    'astroDeltaGammaSec'        : [ True, 'Gaussian',      0.,  0.36,                   -2.,                   2. ], 
    'nuxs'                      : [ True, 'Gaussian',      1.,   0.1,                 0.824,                1.176 ], 
    'nubarxs'                   : [ True, 'Gaussian',      1.,   0.1,                 0.824,                1.176 ], 
    'astroPivot'                : [ True,  'Uniform',      5.,    1.,                    4.,                   6. ], 
    'promptNorm'                : [ True, 'Gaussian',      1.,    1.,                    0.,                   3. ],
    'NeutrinoAntineutrinoRatio' : [ True, 'Gaussian',      1.,    1.,                    0.,                   2. ],
})

#####################################################################################
# set paths to relevant splines and fastMC
#####################################################################################
datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC/compact.fastmc"

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
steering_params.selectionStart                  = float("DnnEnergy_0.99".split("_")[1])
steering_params.evalThreads                     = 1

#####################################################################################
# declare gollumfit object
#####################################################################################
gollumfit = gf.GollumFit(datapaths,steering_params)

#####################################################################################
# Input the nuisance parameter values (at null) and generate and save the expectation
# events. We can then use this as "fake data"
# In general, we can use any combination of nuisance parameter values to see what the 
# expecation would be, and fit to it. This is useful for doing mismodelling tests, for
# example. 
# And when we are ready for the analysis, the real data can be input with the same
# format as the fake data, and the fit to real data can be seamlessly performed. 
#####################################################################################

fitparams = gf.FitParameters()

# set the nuisance parameter values
for sname in syst_dict.keys():
    if sname: 
        exec('fitparams.'+sname+' = syst_dict[\"'+sname+'\"][2]')

null_dist = gollumfit.GetExpectationEvents(fitparams)

print('saving to npz file.')
realization = "nullexpectation.npz"
np.savez(realization, realization=null_dist)


print("Done. Data saved to "+realization)