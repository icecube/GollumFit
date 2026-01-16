# GollumFit {#index}

## Overview

GollumFit is an open-source framework developed by IceCube for performing binned-likelihood, high-energy neutrino telescope analyses. The framework handles the fitting of Monte Carlo simulations to data where the Monte Carlo may be parametrized by a large number of nuisance parameters such as detector efficiencies and neutrino fluxes.

## Key Features

**Core Capabilities:**
- Fast O(40)-dimensional minimization of the binned log-likelihood using the LBFGSB algorithm (from the PhysTools package)
- Integration of the latest atmospheric neutrino flux model, DaemonFlux
- Support for correlations in nuisance parameters as priors
- Built-in nuisance parameters common to neutrino telescopes with a parametrization of the likelihood in terms of these parameters

**Additional Features:**
- **FastMC**: Monte Carlo compression based on analysis binning choice
- **Likelihood Export**: Output likelihood function for Bayesian analyses that sample the likelihood

## Acknowledgements

GollumFit benefits from the different fitting techniques and tools developed over the years by the sterile taskforce. Special acknowledgements to Ben Jones, Chris Weaver, Jordi Salvado, Shivesh Mandalia, Teppei Katori, and Gabriel Collin.

### Core Development

- **Likelihood Framework**: Chris Weaver Likelihood implementation
- **Main Class**: Derived from the Sterilizer by Ben Jones and Carlos Argüelles
- **Python Bindings**: Carlos Argüelles
- **MCMC Integration**: Gabriel Collin and Shivesh Mandalia
- **Astrophysical Oscillations**: Teppei Katori
- **HESE Adaptation**: Carlos Argüelles
- **FastMC Implementation**: Carlos Argüelles
- **Documentation & Maintenance**: Alex Wen

### Major Contributions

- **Dark Matter Physics**: Hrvoje Dujmovic (decay and annihilation)
- **MC Statistical Error Handling**: Austin Schneider (category implementation)
- **Cross Section Physics**: Tianlu Yuan
- **Containerization**: Lukas Hennig

### Contributor List

Spencer Axani, Gabriel Collin, Alejandro Diaz, Hrvoje Dujmovic, Nancy Wandkowsky, Ben J. P. Jones and the UofTexas-Arlington group, Teppei Katori, Shivesh Mandalia, Marjon Moulai, Jordi Salvado, Austin Schneider, Chris Weaver, Tianlu Yuan, Carlos Argüelles-Delgado, Alfonso Garcia Soto, and Alex Wen.

Special thanks to Lukas Henning, Archie Millsop, and Tanvi Krishnan for early testing of the public release.