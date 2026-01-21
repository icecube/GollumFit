# Examples {#examples}

## Overview

This section provides comprehensive, self-contained examples demonstrating GollumFit's capabilities through:
- Toy Monte Carlo datasets
- Nuisance parameter splines and gradients
- Python scripts for common analysis tasks

### Required Resources

**Monte Carlo Data:**

Example Monte Carlo (5 files generated from IceCube software) must be downloaded and extracted:
- **Download**: [https://doi.org/10.7910/DVN/5DSDTD](https://doi.org/10.7910/DVN/5DSDTD)
- **Extract to**: `GollumFit/monte_carlo`

**Nuisance Parameter Resources:**

Splines and gradients for likelihood minimization are stored in:
- **Location**: [GollumFit/resources](../../resources)
- **Contents**: Pre-generated splines and scripts to create custom ones

**Example Scripts:**

All example code and outputs are in:
- **Location**: [GollumFit/examples](../../examples)

### Example Structure

Each example is presented as:
1. Annotated code blocks with explanations
2. Complete scripts available in the `examples` directory
3. Expected outputs and visualizations

---

## FastMC Generation

### Overview

FastMC is a compressed representation of Monte Carlo that significantly improves fitting efficiency by reducing memory footprint while preserving analysis-relevant information. This example demonstrates the generation process.

**Purpose**: Generate compressed Monte Carlo for use in subsequent examples

### Implementation

#### Import Required Modules

```python
import GollumFitPy as gf
import numpy as np
import os
```

#### Configure Data Paths

Set paths for cross section splines and other required data files:

```python
datapaths = gf.DataPaths()
datapaths.neutrino_cc_xs_spline_path             = "../../resources/Splines/CrossSections/sigma_nu_CC_iso.fits"
datapaths.antineutrino_cc_xs_spline_path         = "../../resources/Splines/CrossSections/sigma_nubar_CC_iso.fits"
# ... (additional cross section paths)
```

We then load the fluxes, including the locations of the hadronic and cosmic ray splines (necessary for flux nuisance parameters).

```python
datapaths.conventional_nusquids_atmospheric_file = "../fluxes/atmospheric.hdf5"
datapaths.prompt_nusquids_atmospheric_file       = "../fluxes/prompt_atmospheric.hdf5"
datapaths.astro_nusquids_file                    = "../fluxes/astro.hdf5"

hadronlist = ["he_K+","he_K-","vhe1_pi+","vhe1_pi-","vhe3_K+","vhe3_K-","vhe3_pi+","vhe3_pi-","vhe3_p","vhe3_n"]
crlist     = ["GSF_1","GSF_2","GSF_3","GSF_4","GSF_5","GSF_6"]
for hadron in hadronlist : datapaths.hadronic_spline_path   = os.path.dirname("../fluxes/"+hadron+".hdf5")
for cr in crlist         : datapaths.cosmic_ray_spline_path = os.path.dirname("../fluxes/"+cr+".hdf5")
```

#### Set Steering Parameters

Configure analysis binning and other parameters. **Note**: Binning choices affect FastMC compression and must match the analysis configuration.

```python
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
```

#### Construct and Write FastMC

Create the GollumFit object, generate FastMC with specified compression, and save to file:

```python
gollumfit = gf.GollumFit(datapaths, steering_params)

# Compression parameter: smaller values = higher compression but potential accuracy loss
metascaling = 0.25
gollumfit.ConstructFastMode(metascaling)

# Write to current directory
gollumfit.WriteCompact("compact.fastmc")
```

**Output**: A `compact.fastmc` file containing the compressed FastMC. To use this in subsequent analyses, set `DataPaths.compact_file_path` to this file.

---

## Generating Expectations

### Overview

This example demonstrates how to generate expectation histograms for different nuisance parameter values using a FastMC file. During fitting, GollumFit evaluates expectations repeatedly; here we show how to obtain them explicitly for analysis and validation.

**Use Case**: Examine the effect of varying specific nuisance parameters on the expected event distribution.

**Example**: We'll generate expectations for:
1. Nominal nuisance parameter values
2. Shifted DOM efficiency parameter

### Implementation

#### Import Required Modules

```python
import GollumFitPy as gf
import numpy as np
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats
```

#### Define Nuisance Parameters

Create a dictionary containing all nuisance parameters. Each entry specifies:
- **Name**: Parameter identifier
- **[0]**: `False` to vary in fit, `True` to fix
- **[1]**: Prior distribution type (`'Gaussian'` or `'Uniform'`)
- **[2]**: Central value
- **[3]**: Width (1-sigma for Gaussian)
- **[4]**: Lower bound
- **[5]**: Upper bound

```python
syst_dict = OrderedDict({ 
    'convNorm':          [True, 'Gaussian', 1.0, 0.2, 0.1, 3.0], 
    'zenithCorrection':  [True, 'Gaussian', 0.0, 1.0, -3.0, 3.0], 
    'kaonLosses':        [True, 'Gaussian', 0.0, 1.0, -3.0, 3.0], 
    # ... (additional parameters)
})
```

#### Configure Data Paths

Set paths to required splines and the FastMC file:

```python
datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC/compact.fastmc"
```

#### Set Steering Parameters and Initialize

Configure binning parameters and create the GollumFit object:

```python
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

# Calculate bin counts
energybincount = int((np.log10(steering_params.maxFitEnergy) - np.log10(steering_params.minFitEnergy)) / steering_params.logEbinWidth)
costhbincount = int((steering_params.maxCosth - steering_params.minCosth) / steering_params.cosThbinWidth)

gollumfit = gf.GollumFit(datapaths, steering_params)
```

#### Prepare Output File

Create an HDF5 file to store expectations and corresponding parameter values:

```python
print("Entering bin weights generation loop.")
outfilename = 'example_label.h5'
print(f'Results will be saved to {outfilename}.')

with h5py.File(outfilename, 'a') as hf:
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
```

#### Generate Expectations

Iterate over DOM efficiency values, generate expectations, and save results:

```python
    print("Entering loop to generate weights.")
    np.random.seed(100)  # For reproducibility
    
    # Vary DOM efficiency at two values
    domeff_list = [1.27, 1.27 + 0.02]  # Nominal and shifted
    
    for de in domeff_list: 
        # Initialize parameters to central values
        fitparams = gf.FitParameters()
        for sname in syst_dict.keys():
            exec(f'fitparams.{sname} = syst_dict["{sname}"][2]')
        
        # Set specific DOM efficiency value
        fitparams.domEfficiency = de
        
        # Generate expectation for these parameters
        hist = gollumfit.GetExpectation(fitparams)
        
        # Store nuisance parameter values
        nuisanceparams = np.empty((len(syst_dict),), dtype=np.float64)
        for i, sname in enumerate(syst_dict.keys()):
            exec(f'nuisanceparams[i] = fitparams.{sname}')
        
        # Save to HDF5 file
        # Histograms are separated by event type: starting and throughgoing
        s_h = histdataset.shape
        s_p = paramsdataset.shape
        
        # Resize datasets
        histdataset.resize((s_h[0] + 1, s_h[1], s_h[2], s_h[3]))
        paramsdataset.resize((s_p[0] + 1, s_p[1]))
        
        # Append new data
        histdataset[s_h[0]:] = hist
        paramsdataset[s_p[0]:] = nuisanceparams
```

**Output**: Expectation histograms (2D: energy × zenith, for starting and throughgoing events) with corresponding parameter values saved to `example_label.h5`.

### Visualization

Plot the generated expectations:

```python
import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib import rc
rc('text', usetex=True)

file = h5py.File('example_label.h5', 'r')

nominal = file['hists'][0]
perturbed = file['hists'][1]

#hist = nominal
hist = perturbed

fig, axes = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
aspect = (0.4)
vmin=1
vmax=10000

from matplotlib.colors import LogNorm

# Create a log-scaled color map
norm = LogNorm(vmin=vmin, vmax=vmax)

# Plot the first subplot (left)
im1 = axes[0].imshow(hist[0].T, cmap='GnBu', origin='lower',extent=[-1, 0, np.log10(300), np.log10(1e5)],aspect=aspect, norm=norm)
axes[0].set_ylabel(r'$\log_{10}$($E$/GeV)')
axes[0].set_xlabel(r'$\cos(\theta)$')
axes[0].set_title('starting')

# Plot the second subplot (right)
im2 = axes[1].imshow(hist[1].T, cmap='GnBu', origin='lower', extent=[-1, 0, np.log10(300), np.log10(1e5)],aspect=aspect, norm=norm)
axes[1].set_ylabel(r'$\log_{10}$($E$/GeV)')
axes[1].set_xlabel(r'$\cos(\theta)$')
axes[1].set_title('throughgoing')

# Add a colorbar to both subplots
cbar1 = fig.colorbar(im1, ax=axes[0],fraction=0.046, pad=0.04, label=r'$N$')
cbar2 = fig.colorbar(im2, ax=axes[1],fraction=0.046, pad=0.04, label=r'$N$')

xedges = np.linspace(-1, 0, 20 + 1)
yedges = np.linspace(np.log10(300), np.log10(1e5), 24 + 1)

for i in range(len(xedges) - 1):
    for j in range(len(yedges) - 1):
        axes[0].text(xedges[i] + 0.5 * (xedges[i + 1] - xedges[i]),
                 yedges[j] + 0.5 * (yedges[j + 1] - yedges[j]),
                 f'{int(hist[0][i, j])}',
                 color='black',
                 ha='center',
                 va='center', fontsize=4)
        axes[1].text(xedges[i] + 0.5 * (xedges[i + 1] - xedges[i]),
                 yedges[j] + 0.5 * (yedges[j + 1] - yedges[j]),
                 f'{int(hist[1][i, j])}',
                 color='black',
                 ha='center',
                 va='center', fontsize=4)

plt.tight_layout()
plt.savefig('example_hist.pdf')
```

---

## Fitting to Null Hypothesis

### Overview

This example demonstrates GollumFit's core functionality: fitting Monte Carlo to data using likelihood minimization. We use **null pseudo-data** (expectation with nominal nuisance parameters) as our test dataset.

**Workflow**:
1. Generate null pseudo-data (expectation at nominal parameters)
2. Initialize fit with random parameter values
3. Minimize likelihood to recover nominal parameters

**Expected Result**: Successful fit should recover the nominal parameter values used to generate the pseudo-data.

### Part 1: Generate Pseudo-Data

#### Setup and Configuration

```python
import GollumFitPy as gf
import numpy as np
import os
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats

# Define nuisance parameters (as in previous example)
syst_dict = OrderedDict({ 
    'convNorm':          [True, 'Gaussian', 1.0, 0.2, 0.1, 3.0], 
    'zenithCorrection':  [True, 'Gaussian', 0.0, 1.0, -3.0, 3.0], 
    'kaonLosses':        [True, 'Gaussian', 0.0, 1.0, -3.0, 3.0], 
    # ... (additional parameters)
})
```

#### Configure Paths and Parameters

Set up data paths and steering parameters (similar to previous examples):

```python
datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC/compact.fastmc"

steering_params = gf.SteeringParams()
steering_params.minFitEnergy                    = 300
steering_params.maxFitEnergy                    = 1e5
steering_params.logEbinEdge                     = np.log10(300)
steering_params.logEbinWidth                    = (np.log10(1e5)-np.log10(300))/24
steering_params.minCosth                        = -1
steering_params.maxCosth                        = 0.0
steering_params.cosThbinEdge                    = 0.0
steering_params.cosThbinWidth                   = 0.05
steering_params.selectionStart                  = float("DnnEnergy_0.99".split("_")[1])
steering_params.evalThreads                     = 1
```

#### Generate and Save Null Expectation

Create GollumFit object, set parameters to nominal values, and generate pseudo-data:

```python
gollumfit = gf.GollumFit(datapaths, steering_params)

fitparams = gf.FitParameters()

# Set all parameters to nominal (central) values
for sname in syst_dict.keys():
    if sname: 
        exec(f'fitparams.{sname} = syst_dict["{sname}"][2]')
```

#### Export Pseudo-Data Events

Generate weighted events corresponding to the nominal parameters and save:

```python
# Get weighted event list for nominal parameters
null_dist = gollumfit.GetExpectationEvents(fitparams)

print('Saving to npz file.')
realization = "nullexpectation.npz"
np.savez(realization, realization=null_dist)

print(f"Done. Data saved to {realization}")
```

**Note**: This method generates pseudo-data with arbitrary parameter values, useful for sensitivity studies. For real analyses, substitute actual data at the unblinding stage using the same interface.

### Part 2: Fitting to Pseudo-Data

#### Import and Setup

Now we'll fit to the pseudo-data, starting from random initial parameter values. Success means recovering the nominal values.

**Import modules and define parameters** (note: `False` means vary in fit):

```python
import GollumFitPy as gf
import numpy as np
import os
import sys
import subprocess
import scipy.stats as stats

syst_dict = OrderedDict({ 
    'convNorm' :          [ False, 'Gaussian', 1.0, 0.2, 0.1, 3.0 ], 
    'zenithCorrection' :  [ False, 'Gaussian', 0.0, 1.0, -3.0, 3.0 ], 
    'kaonLosses' :        [ False, 'Gaussian', 0.0, 1.0, -3.0, 3.0 ], 
    # ...
})
```

#### Define Random Sampling Function

Create a helper function to sample random parameter values from their prior distributions:

```python
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
```

#### Initialize Fit Parameter Objects

Create objects to manage fit configuration:

```python
fitparams_flag  = gf.FitParametersFlag()  # Which parameters to vary
fitparams_bound = gf.FitParametersBound()  # Parameter bounds
priors          = gf.Priors()              # Prior distributions
seed_fitparams  = gf.FitParameters()       # Initial values
```

#### Set Priors and Random Initial Values

Configure priors and randomly initialize starting parameter values:

```python
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
```

#### Set Correlations and Paths

```python
gollumdir = "../../gollumfit-release"

# Set correlations, which are not required for evaluating the likelihood but are required for fitting/minimization
iceg_corr = np.load('../../resources/IceGradientsMaker/icegrad_correlations.npy')
flux_corr = np.load('../../resources/DDMFluxMaker/flux_correlations_new_ddmnodeis.npy')
for idx, val in np.ndenumerate(iceg_corr): priors.SetIceGradientsCorr(idx[0], idx[1], val)
for idx, val in np.ndenumerate(flux_corr): priors.SetFluxCorr(idx[0], idx[1], val)

datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC/compact.fastmc"
```

#### Configure Steering Parameters

Set binning and convergence criteria (must match FastMC binning):

```python
edges = np.logspace(np.log10(300), np.log10(1e5), 25)

steering_params = gf.SteeringParams()
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
```

#### Load Data and Configure Fit

```python
gollumfit = gf.GollumFit(datapaths, steering_params)
realization = "nullexpectation.npz"
total_data = gollumfit.SetData(np.load(realization)["realization"])
```

We use setter methods to set the priors on the nuisance parameters and also initialize the seed, which determines from which nuisance parameter values the fit begins.

```python
gollumfit.SetFitParametersFlag(fitparams_flag)
gollumfit.SetFitParametersBound(fitparams_bound)
gollumfit.SetFitParametersPriors(priors)
gollumfit.SetFitParametersSeed([seed_fitparams])
gollumfit.ConstructLikelihoodProblem()
```

#### Perform Likelihood Minimization

Run the LBFGSB minimization algorithm:

```python
print("Starting minimization...")
min_llh = gollumfit.MinLLH()

for sname in syst_dict.keys():
    exec('print("' + sname + '", min_llh.params.' + sname + ')')
    exec('systematics += str(min_llh.params.' + sname + ') + " "')

print('llh:', min_llh.likelihood)
print('nEval: ' + str(min_llh.nEval))

print("Completed successfully. Bye!")
```

In the figure, we show a pull plot of the results. We plot the randomly-initialized nuisance parameter value and the corresponding "best fit" value after the fit is performed; this is shown in terms of pulls with respect to the nuisance parameter's respective central value and width. We see that the fit has successfully converged on the null nuisance parameters that were used to generate the pseudo-data. 

(Figure temporarily omitted due to ongoing publication review process.)
![The result of minimization, comparing the best-fit and initial value of the nuisance parameters.](figures/inject_recover_null.pdf)

The parameter `astroPivot`, in our fit, has no influence on the shape of the flux and therefore no effect on the likelihood. Hence there is no preferred value to be recovered.

---

## Adding a new systematic parameter

In many analyses, we may want to introduce new physics or nuisance parameters, or even extend to other experiments that may be characterized differently. Hence, there is often a need to include additional parameters to fit over. In the current edition of GollumFit, this requires modifications to the source code. Here, we outline the necessary changes in the code for adding a new systematic parameter to fit over.

1. In the files `GollumParameters.*`, append the new parameter and its priors (width, flag, center, etc.) to the appropriate objects.

2. Determine how to weight the new parameter; i.e., how the event weight will shift with the value of the parameter. This can be done with a spline, a gradient, or some analytical expression. A weighter struct must be created in `analysisWeighting.h` that inherits from `phys_tools::GenericWeighter`. This struct should be written so that when it is constructed, the parameter value is one of the inputs, which can be any data type but must include the ability to accept a `phys_tools::autodiff::FD<Dim>` type (where `Dim` is the number of dimensions in the fitting). This struct must also have an overloaded `()` operator which accepts an event type as input and outputs a number representing the change in event weight. It must be able to output a `phys_tools::autodiff::FD<Dim>` type. Of course, this weighter will need to contain the means to calculate the weight, so input from a spline file, or a formula, or any other method, is necessary. This step will be the most sophisticated—the other weighters in the same file serve as good examples on how to write this new weighter.

3. In the same source file, within the `WeighterMaker` struct, assign the new parameter to a variable from the list `params`, and construct an instance of the weighter struct you wrote previously. Determine how the factor will shift the weight by adding it to the expression to calculate the return value.

4. If the weighting requires some external spline, it should be loaded via a function declared in `GollumFit.*`.

5. In the files `GollumFit.*`, the function `GollumFit::ConstructLikelihoodProblem` should have the priors of the new parameter added to the existing lists of priors. Furthermore, the functions `GollumFit::MinLLH()`, as well as the helpers `GollumFit::ConvertFitParametersFlag()`, `GollumFit::ConvertFitParameters()`, and so on, should all have the new parameter appended in a straightforward way.

6. There are some locations in the code, for example in the file `GollumFit.cpp` in the function `GollumFit::EvalLLHGradient`, and some `assert` checks in the same file, where the number of nuisance parameters (38 by default) is hard-coded. This must be adjusted to account for any new nuisance parameters.

---

## Adding a new fitting dimension

Another possible modification that may be desirable is the ability to change what is histogrammed over. By default, since we are working primarily with neutrino telescopes, we concern ourselves primarily with analyses that are binned in the space of energy and declination (or zenith) of the event. If we simply wish to keep the same two dimensions but adjust the bin size, this is straightforward as we just need to set the right options in the `steering_params` object that is input to the `GollumFit` constructor. However, if additional binning dimensions need to be added, such as the azimuth or inelasticity, this requires some changes to the source code. We will briefly outline the key steps below.

1. If new information needs to be introduced to the event structure, it must be added in the event class that is written in `Event.*`.

2. The new binning information must be added to variables in the `SteeringParams` structure written in `GollumParameters.*`.

3. Within the `GollumFit.*` files, some functions will loop over the bins - these loops must be updated so that the new dimension is looped over as well. In particular, `GollumFit::GetDataDistribution`, `GollumFit::GetWeightedExpectation`, and `GollumFit::GetRealization` must be updated to include the new binning.

4. Within the same files, the function `GollumFit::ConstructFastMode` must also be updated to reflect the new binning, especially within the `MetaHistType metaHist` object. It is important to use the same binning scheme to generate FastMC and to run the analysis, as the binning will affect the way in which FastMC will compress the events.