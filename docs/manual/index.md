# GollumFit               {#index}

GollumFit is an open-source framework developed by IceCube for performing 
binned-likelihood, high energy, neutrino telescope analyses. 
In particular, GollumFit handles the fitting of monte carlo
to data where the monte carlo may be parametrized by a large number of 
nuisance parameters such as detector efficiencies and neutrino fluxes. 

### Key Features

GollumFit: 
* can perform a fast, high-dimensional (i.e. O(40)), minimization of the binned
log-likelihood using the LBFGSB algorithm (from the PhyTools package)
* integrates the latest atmospheric neutrino flux model, DaemonFlux
* integrates correlations in these parameters as priors
* has built-in nuisance parameters common to neutrino telescopes and a parametrizeation of the likelihood in terms of these parameters

GollumFit can also secondarily: 
* reduce the size of the monte carlo depending on the analysis binning choice, 
a feature known as FastMC
* output the likelihood function, for 
performing Bayesian analyses which sample the likelihood

### Acknowledgements

GollumFit benefits from the different fitting techniques and tools
developed over the years by the sterile taskforce. Special acknoledments to
Ben Jones, Chris Weaver, Jordi Salvado, Shivesh Mandalia, Teppei Katori,
and Gabriel Collin.

GollumFit uses the Chris Weaver Likelihood.
The GollumFit main class comes from the Sterilizer by Ben Jones
and Carlos Arguelles. Python bindings were written by Carlos Arguelles.
MCMC connections were curated by Gabriel Collin and Shivesh Mandalia.
Teppei Katori wrote the analytic formulas for the astrophysical
neutrino oscillations. Jordi's handywork and spirit permeate the whole
code. Adaptation to the HESE framework was done by Carlos Arguelles.

Many amazing improvements and bug fixes have been introduced in this
code since the first version. Lots of great contributions by Austin
Schneider, Tianlu Yuan, Hrvoje Dujmovic, and Nancy Wandkowsky.
Implementation of dark matter decay and annihilation was done by Hrvoje.
Implementation of categories for MC statistical error handling was done by
Austin. Implementation of cross section physics by Tianlu.
Implementation of the fast mode was done by Carlos.
Recent clean-up and documentation writing performed by Alex Wen.