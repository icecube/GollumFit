# Examples

To comprehensively demonstrate the capabilities of GollumFit, we provide a set of toy Monte Carlo, nuisance parameter splines & gradients, and a series of Python scripts to give a set of self-contained examples on the usage of `GollumFit`. In particular, we provide some examples on generating expectations, minimizing the LLH, and generating FastMC.

The Monte Carlo is generated from some version of IceCube software and includes 5 files. This is stored in `GollumFit/monte_carlo`. The splines and gradients required for the minimization of the likelihood with respect to the nuisance parameters are stored in `GollumFit/resources`, which also include the scripts to generate some of them. Finally, the actual example scripts and their outputs are contained in `GollumFit/examples`.

The examples featured here are split into code blocks with accompanying explanation. 
The full example scripts are contained in the `examples` directory.

## FastMC generation

We show the procedure of generating FastMC, which is used to generate an iteration of the Monte Carlo. This is necessary for subsequent examples and is an important step in ensuring that the fitting is efficient.

We begin with necessary imports.

```python
import GollumFitPy as gf
import numpy as np
import os
```

We then set the paths of the splines and cross sections by loading them into a `DataPaths` object.

```python
datapaths = gf.DataPaths()
datapaths.neutrino_cc_xs_spline_path             = "../../resources/Splines/CrossSections/sigma_nu_CC_iso.fits"
datapaths.antineutrino_cc_xs_spline_path         = "../../resources/Splines/CrossSections/sigma_nubar_CC_iso.fits"
...
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

We then set the steering parameters - in particular, we set the binning for the histograms. This information is loaded into a `SteeringParams` object.

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

We declare a `GollumFit` object, construct the FastMC, and write it to file. The inputs are the `DataPaths` and `SteeringParams` objects that we declared and filled previously.

```python
gollumfit = gf.GollumFit(datapaths,steering_params)

# metascaling controls the FastMC compression. Large values can cause loss of accuracy
metascaling = "0.25"
gollumfit.ConstructFastMode(float(metascaling))

gollumfit.WriteCompact(".")
```

The output is a `STERILE.meows` file that contains the FastMC, which can be loaded in another instance by assigning its directory path to the `DataPaths.compact_file_path` variable.

## Generating Expectations

In this example, we explain how to set up the GollumFit object and output an expectation histogram given certain nuisance parameters and a FastMC file. While performing a fit, `GollumFit` will internally evaluate this expectation many times; here, we show how to obtain it. In particular, we consider the expectation given the nominal set of nuisance parameters, and likewise for a set of nuisance parameters with a shifted value of the DOM efficiency nuisance parameter. It is a good way to explicitly see the effect of varying a particular nuisance parameter and is useful for various tests that one may want to do surrounding a complex physics analysis. Of course, in general, we may vary any nuisance parameter in any combination, in order to see the desired effect on the expectation.

As usual, we begin with imports.

```python
import GollumFitPy as gf
import numpy as np
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats
```

We then declare a dictionary containing all the nuisance parameters. In each entry, the name is indicated, followed by a list containing a boolean, string, and four floats. They indicate, respectively, whether the parameter is ignored in the fit (`False` to enable in fitting), the prior shape (`Gaussian` or `Uniform`), the center, the width, the lower bound, and the upper bound of the prior.

```python
syst_dict     = OrderedDict({ 
    'convNorm' :          [ True, 'Gaussian', 1. 0.2, 0.1, 3. ], 
    'zenithCorrection' :  [ True, 'Gaussian',  0., 1., -3.,  3. ], 
    'kaonLosses' :        [ True, 'Gaussian', 0., 1., -3., 3. ], 
...
})
```

We then set the paths of the splines and FastMC to be used in a fitting by loading them into a `DataPaths` object.

```python
datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC"
```

Next, we set the steering parameters which primarily determine the binning for the histograms. Also, we declare an instance of the main `GollumFit` object again with the `datapaths` and `steering_params` as the input.

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

energybincount = (np.log10(steering_params.maxFitEnergy) - np.log10(steering_params.minFitEnergy))/steering_params.logEbinWidth
costhbincount = (steering_params.maxCosth - steering_params.minCosth)/steering_params.cosThbinWidth

gollumfit = gf.GollumFit(datapaths,steering_params)
```

Now, we set up a `.h5` file with the correct structure to contain the expectation histogram, so that we may save to it later. We will also a list where we can save the corresponding values of the nuisance parameters.

```python
print("Entering bin weights generation loop.")
outfilename = 'example_label'+'.h5'
print('Results will be saved to '+outfilename+".")

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

Continuing from the previous block, we set a random seed to have reproducible results. We then specify the two values of the `domEfficiency` that we want to evaluate the expectation at. We iterate over these values and for each, set the value of the `domEfficiency` parameter in the `fitparams` object. To evaluate the expectation, the `GetExpectation` method is called with the `fitparams` object as an input, so that the updated nuisance parameter value is being included. Finally, we save the nuisance parameter values that were used to a list.

```python
    print("Entering loop to generate weights.")
    np.random.seed(100)
    
    # we vary the dom eff at two values, and generate the corresponding expectation at each
    dofEff_list = [1.27, 1.27+0.02]
    
    for de in dofEff_list: 
        fitparams = gf.FitParameters()
        # set the nuisance parameters to the central values stored in the dict we declared
        for sname in syst_dict.keys():
            exec('fitparams.'+sname+' = syst_dict[\"'+sname+'\"][2]')
        exec('fitparams.domEfficiency = de')
        
        hist = gollumfit.GetExpectation(fitparams)
        
        # recover the 38 nusiance param values too
        nusianceparams = np.empty((len(syst_dict),), dtype=np.float64)
        for i in range(len(list(syst_dict.keys()))):
            exec('nusianceparams[i] = fitparams.'+list(syst_dict.keys())[i])
```

In the final part, we save the expectation output (a histogram) to the `.h5` file. For clarity, we also save the corresponding list of nuisance parameters that were input. Note that we have a separate histogram for starting and throughgoing events, so that we are working with two 2D arrays.

```python
        # Get the current shape of the dataset
        s_h = histdataset.shape
        s_p = paramsdataset.shape

        # Calculate the new shape after appending
        s_h_new = (s_h[0] + 1, s_h[1], s_h[2], s_h[3])
        s_p_new = (s_p[0] + 1, s_p[1])
    
        # Resize the dataset to accommodate the new data
        histdataset.resize(s_h_new)
        paramsdataset.resize(s_p_new)
        
        # Append the new data to the dataset
        histdataset[s_h[0]:] = hist
        paramsdataset[s_p[0]:] = nusianceparams
```

The output of this example is plotted in figure \ref{fig:expectation_variation}, where we show the nominal expectation, and beside it the pull plot given the variation in `domEfficiency` that we have inserted here.

The code that follows is a script to plot the expectation that we have just generated and saved to the `.h5` file.

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

## Fitting to null

In this example, we showcase the core functionality of `GollumFit`, in a fit of the Monte Carlo to data. In this case, our data is the null pseudo-data, which is simply the expectation generated with the nominal or central values of all nuisance parameters, which we are going to treat as the data to fit to. The example consists of two parts: generating the pseudo-data to fit to, and then performing the actual fit.

### Generating Pseudo-Data

To generate pseudo-data, we begin with the usual procedure of necessary imports and declaring the relevant nuisance parameters in a dictionary.

```python
import GollumFitPy as gf
import numpy as np
import os
import sys
import h5py
from collections import OrderedDict
import scipy.stats as stats

syst_dict = OrderedDict({ 
    'convNorm' :          [ True, 'Gaussian', 1.0, 0.2, 0.1, 3.0 ], 
    'zenithCorrection' :  [ True, 'Gaussian', 0.0, 1.0, -3.0, 3.0 ], 
    'kaonLosses' :        [ True, 'Gaussian', 0.0, 1.0, -3.0, 3.0 ], 
    # ...
})
```

As before, we also set the paths to the relevant spline and FastMC files, and also set the steering parameters.

```python
datapaths = gf.DataPaths()
datapaths.domeff_spline_path      = "../../resources/Splines/DOMEffSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.holeice_spline_path     = "../../resources/Splines/HoleIceSplines/new_ddmnodeis/BDT/DnnEnergy_0.99"
datapaths.attenuation_spline_path = "../../resources/Splines/AttenuationSplines/new_ddmnodeis"
datapaths.compact_file_path       = "../FastMC"

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

We then declare the main `GollumFit` object in the standard way, and set the `fitparameters` object variables to those defined in the nuisance parameter dictionary `syst_dict`.

```python
gollumfit = gf.GollumFit(datapaths,steering_params)

fitparams = gf.FitParameters()

# set the nuisance parameter values
for sname in syst_dict.keys():
    if sname: 
        exec('fitparams.' + sname + ' = syst_dict["' + sname + '"][2]')
```

Finally, we use the method `GetExpectationEvents` to obtain a list of weighted events that correspond to the input nuisance parameters. We then save it to a numpy file `nullexpectation.npz`, which can be used later on.

```python
null_dist = gollumfit.GetExpectationEvents(fitparams)

print('saving to npz file.')
realization = "nullexpectation.npz"
np.savez(realization, realization=null_dist)

print("Done. Data saved to " + realization)
```

We emphasize that this is a way to make pseudo-data with arbitrary values of nuisance parameters, which is a useful way to perform sensitivity studies or experiment with fitting; real data is fed in the same way, and can be swapped in at the unblinding stage.

### Fitting to the Pseudo-Data

Now we will use the pseudo-data output of the previous part to perform a fit. We will randomly initialize values of the nuisance parameters, and then fit to the pseudo-data that we created, treating it as real data. If the fit succeeds, we should recover the central (null) values of the nuisance parameters, which we used to generate the pseudo-data.

To start, we perform the necessary imports and declare the nuisance parameters in a dictionary as usual. Note here that in the dictionary, we have `False` for all parameters as we are now varying them for the fit.

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

We define a helper function that throws a random value of the nuisance parameter taking into account its range, center, and width.

```python
def throw(syst):
    if syst[1] == 'Gaussian': 
        val = stats.truncnorm((syst[4] - syst[2]) / syst[3], (syst[5] - syst[2]) / syst[3], syst[2], syst[3]).rvs(1)
        return val[0]
    else: 
        return np.random.uniform(syst[4], syst[5])
```

We declare various objects to characterize the nuisance parameters - their flag (whether they are fit over), bounds, prior value, and initialized value.

```python
fitparams_flag  = gf.FitParametersFlag()
fitparams_bound = gf.FitParametersBound()
priors          = gf.Priors()
seed_fitparams  = gf.FitParameters()
```

Now we declare a random seed for reproducibility, and iterate over the nuisance parameters. For each parameter, we set the prior centers and widths and ranges. We then use the `throw` function we defined earlier to set each nuisance parameter to a random initial value according to their prior distributions. The fit will commence from this value.

```python
np.random.seed(100)
print('Initializing with the following randomly-seeded nuisance params:')

for sname in syst_dict.keys():
    exec('fitparams_flag.' + sname + ' = syst_dict["' + sname + '"][0]')
    exec('fitparams_bound.' + sname + 'Min = syst_dict["' + sname + '"][4]')
    exec('fitparams_bound.' + sname + 'Max = syst_dict["' + sname + '"][5]')
    if syst_dict[sname][1] == 'Gaussian':
        exec('priors.' + sname + 'Center = syst_dict["' + sname + '"][2]')
        exec('priors.' + sname + 'Width  = syst_dict["' + sname + '"][3]')
    else:
        exec('priors.' + sname + 'Min = syst_dict["' + sname + '"][4]')
        exec('priors.' + sname + 'Max = syst_dict["' + sname + '"][5]')
        
    thrown_val = throw(syst_dict[sname])
    exec('seed_fitparams.' + sname + ' = thrown_val')
    print(sname + ' ' + str(thrown_val))
```

Here, we set the correlations between the ice gradient and flux parameters. We also set the paths to the necessary spline and FastMC files, as in the previous examples.

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
datapaths.compact_file_path       = "../FastMC"
```

We also set the steering parameters, as in previous examples. In particular, the binning must match the binning used to generate the FastMC.

```python
edges = np.logspace(np.log10(300), np.log10(1e5), 25)
steering_params = gf.SteeringParams()
steering_params.minFitEnergy   = edges[0]
steering_params.maxFitEnergy   = edges[-1]
steering_params.logEbinEdge    = np.log10(edges[0])
steering_params.logEbinWidth   = np.log10(edges[1])-np.log10(edges[0])
steering_params.minCosth       = -1
steering_params.maxCosth       = 0.0
steering_params.cosThbinEdge   = 0.0
steering_params.cosThbinWidth  = 0.05
steering_params.selectionStart = float("DnnEnergy_0.99".split("_")[1])
steering_params.evalThreads    = 1
steering_params.change_tol     = 1.e-20
steering_params.grad_tol       = 1.e-20
steering_params.uncertaintyModSigmaOverMu = 0.0
```

We construct a `GollumFit` object and use the method `SetData` to read in the pseudo-data that we generated previously. In a real analysis, we would simply feed in the real data instead.

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

Finally, we invoke the `MinLLH` method to perform the minimization. This may typically take several hours. Afterwards, we print the best-fit nuisance parameters, as well as the minimum likelihood value and the number of likelihood evaluations that were required.

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

## Adding a new systematic parameter

In many analyses, we may want to introduce new physics or nuisance parameters, or even extend to other experiments that may be characterized differently. Hence, there is often a need to include additional parameters to fit over. Unfortunately, this requires modifications to the source code. Here, we outline the necessary changes in the code for adding a new systematic parameter to fit over.

1. In the files `GollumParameters.*`, append the new parameter and its priors (width, flag, center, etc.) to the appropriate objects.

2. Determine how to weight the new parameter; i.e., how the event weight will shift with the value of the parameter. This can be done with a spline, a gradient, or some analytical expression. A weighter struct must be created in `analysisWeighting.h` that inherits from `phys_tools::GenericWeighter`. This struct should be written so that when it is constructed, the parameter value is one of the inputs, which can be any data type but must include the ability to accept a `phys_tools::autodiff::FD<Dim>` type (where `Dim` is the number of dimensions in the fitting). This struct must also have an overloaded `()` operator which accepts an event type as input and outputs a number representing the change in event weight. It must be able to output a `phys_tools::autodiff::FD<Dim>` type. Of course, this weighter will need to contain the means to calculate the weight, so input from a spline file, or a formula, or any other method, is necessary. This step will be the most sophisticated—the other weighters in the same file serve as good examples on how to write this new weighter.

3. In the same source file, within the `WeighterMaker` struct, assign the new parameter to a variable from the list `params`, and construct an instance of the weighter struct you wrote previously. Determine how the factor will shift the weight by adding it to the expression to calculate the return value.

4. If the weighting requires some external spline, it should be loaded via a function declared in `GollumFit.*`.

5. In the files `GollumFit.*`, the function `GollumFit::ConstructLikelihoodProblem` should have the priors of the new parameter added to the existing lists of priors. Furthermore, the functions `GollumFit::MinLLH()`, as well as the helpers `GollumFit::ConvertFitParametersFlag()`, `GollumFit::ConvertFitParameters()`, and so on, should all have the new parameter appended in a straightforward way.

6. There are some locations in the code, for example in the file `GollumFit.cpp` in the function `GollumFit::EvalLLHGradient`, and some `assert` checks in the same file, where the number of nuisance parameters (38 by default) is hard-coded. This must be adjusted to account for any new nuisance parameters.

## Adding a new fitting dimension

Another possible modification that may be desirable is the ability to change what is histogrammed over. By default, since we are working primarily with neutrino telescopes, we concern ourselves primarily with analyses that are binned in the space of energy and declination (or zenith) of the event. If we simply wish to keep the same two dimensions but adjust the bin size, this is straightforward as we just need to set the right options in the `steering_params` object that is input to the `GollumFit` constructor. However, if additional binning dimensions need to be added, such as the azimuth or inelasticity, this requires some changes to the source code. We will briefly outline the key steps below.

1. If new information needs to be introduced to the event structure, it must be added in the event class that is written in `Event.*`.

2. The new binning information must be added to variables in the `SteeringParams` structure written in `GollumParameters.*`.

3. Within the `GollumFit.*` files, some functions will loop over the bins - these loops must be updated so that the new dimension is looped over as well. In particular, `GollumFit::GetDataDistribution`, `GollumFit::GetWeightedExpectation`, and `GollumFit::GetRealization` must be updated to include the new binning.

4. Within the same files, the function `GollumFit::ConstructFastMode` must also be updated to reflect the new binning, especially within the `MetaHistType metaHist` object. It is important to use the same binning scheme to generate FastMC and to run the analysis, as the binning will affect the way in which FastMC will compress the events.