"""
Generate Expectation Histograms Example

This script demonstrates how to generate expectation histograms for different
nuisance parameter values using a FastMC file.

Example command:
    nohup time python ./generate_expectation.py example_label > generate_expectation.out 2>&1 &
"""

import GollumFitPy as gf
import numpy as np
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats

expectation_label = sys.argv[1]

#####################################################################################
# Define Nuisance Parameters
# Format: [vary_flag, prior_type, center, width, lower_bound, upper_bound]
# - vary_flag: False to vary in fit, True to fix
# - prior_type: 'Gaussian' or 'Uniform'
# - center: Central value
# - width: 1-sigma width for Gaussian
# - lower_bound, upper_bound: Parameter bounds
#####################################################################################
syst_dict = OrderedDict({ 
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
datapaths.compact_file_path       = "../FastMC/example.fastmc"

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

# Calculate bin counts for output file structure
energybincount = int((np.log10(steering_params.maxFitEnergy) - 
                      np.log10(steering_params.minFitEnergy)) / steering_params.logEbinWidth)
costhbincount = int((steering_params.maxCosth - steering_params.minCosth) / 
                    steering_params.cosThbinWidth)

#####################################################################################
# Initialize GollumFit
#####################################################################################
gollumfit = gf.GollumFit(datapaths, steering_params)

#####################################################################################
# Prepare Output File
# Create an HDF5 file to store expectations and corresponding parameter values
#####################################################################################
print("Setting up expectation generation.")
outfilename = f'{expectation_label}.h5'
print(f'Results will be saved to {outfilename}.')

with h5py.File(outfilename, 'a') as hf:
    #####################################################################################
    # Access the existing dataset within the file or create a new one if it doesn't exist
    #####################################################################################
    if ('hists' in hf) and ('params' in hf):
        histdataset = hf['hists']
        paramsdataset = hf['params']
    else:
        # Create a new datasets with an initial shape
        histdataset = hf.create_dataset('hists', 
                                    shape=(0, 2, costhbincount, energybincount), 
                                    maxshape=(None, 2, costhbincount, energybincount), 
                                    dtype='float64')
        paramsdataset = hf.create_dataset('params', 
                                    shape=(0, len(syst_dict)), 
                                    maxshape=(None, len(syst_dict)), 
                                    dtype='float64')
        paramsdataset.attrs['syst_dict_keys'] = list(syst_dict.keys())
    print("Successfully created or accessed datasets in output file.")    
    #####################################################################################
    # enter a loop to generate toy events
    #####################################################################################
    print("Entering loop to generate weights.")
    np.random.seed(100)
    
    # we vary the dom eff at two values, and generate the corresponding expectation at each
    dofEff_list = [1.27, 1.27+0.02]
    
    for de in dofEff_list: 
        #####################################################################################
        # make the histogram
        #####################################################################################
        fitparams = gf.FitParameters()
        # set the nuisance parameters to the central values stored in the dict we declared
        for sname in syst_dict.keys():
            exec('fitparams.'+sname+' = syst_dict[\"'+sname+'\"][2]')
        exec('fitparams.domEfficiency = de')
        
        hist = gollumfit.GetExpectation(fitparams)
        
        # recover the 38 nusiance param values too
        nuisanceparams = np.empty((len(syst_dict),), dtype=np.float64)
        for i in range(len(list(syst_dict.keys()))):
            exec('nuisanceparams[i] = fitparams.'+list(syst_dict.keys())[i])
        #####################################################################################
        # Save to HDF5 File
        # Histograms are separated by event type: starting and throughgoing
        #####################################################################################
        # Get current dataset shapes
        s_h = histdataset.shape
        s_p = paramsdataset.shape
        
        # Resize datasets to accommodate new data
        histdataset.resize((s_h[0] + 1, s_h[1], s_h[2], s_h[3]))
        paramsdataset.resize((s_p[0] + 1, s_p[1]))
        
        # Append new data
        histdataset[s_h[0]:] = hist
        paramsdataset[s_p[0]:] = nuisanceparams
        ###########################################################################

print("Done. Data saved to "+outfilename)    