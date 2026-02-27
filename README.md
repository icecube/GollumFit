# GollumFit

GollumFit is an open-source framework developed by IceCube for performing binned-likelihood, high energy, diffuse, neutrino telescope analyses. In particular, GollumFit handles the fitting of monte carlo to data where the monte carlo may be parametrized by a large number of model parameters such as detector efficiencies and neutrino fluxes. We include model parameters common to all neutrino telescopes.

For details on the code, installation, usage, and examples, please visit the official documentation at [https://docs.icecube.aq/gollumfit/main/index.html](https://docs.icecube.aq/gollumfit/main/index.html).


## Getting Started

To build and run GollumFit, follow the instructions in the [installation manual](https://docs.icecube.aq/gollumfit/main/md__home_runner_work_GollumFit_GollumFit_docs_manual_installation.html).

### Diver

This branch adds support for [Diver](https://github.com/diveropt/Diver), a differential evolution optimizer. Requires `gfortran` and MPI (`mpif90`). To install:

```bash
git clone https://github.com/diveropt/Diver.git
cd Diver
make libdiver.so
```

Then add Diver to your environment:
```bash
export DIVER_DIR=/path/to/Diver
export LD_LIBRARY_PATH=$DIVER_DIR/lib:$LD_LIBRARY_PATH
export CPLUS_INCLUDE_PATH=$DIVER_DIR/include:$CPLUS_INCLUDE_PATH
```


## Folder Descriptions

- **src/**  
  Contains the C++ source code for GollumFit, including core modules for analysis weighting, event handling, and parameter management.

- **include/**  
  Contains header files (.h) for the C++ source code, along with additional headers such as `json.hpp` needed by the project.

- **examples/**  
  Provides several example use cases:
  - **expectations/**: Scripts and files for generating and plotting expectation histograms.
  - **FastMC/**: Scripts for fast Monte Carlo simulations.
  - **fitting_to_null/**: Example scripts for generating fake data, computing likelihoods, and running related SLURM jobs.
  - **fluxes/**: A collection of HDF5 files containing various neutrino flux models (e.g., astrophysical, atmospheric, etc.).

- **monte_carlo/**  
  Contains a subdirectory (STERILE) with licensing files for Monte Carlo generation data used by GollumFit.

- **docs/**  
  Includes Doxygen configuration files and CSS for generating API documentation, as well as Markdown files for the user manual, code structure, and installation instructions.

- **python/**  
  Contains files to build the Python extension module (GollumFitPy), including source code, setup scripts, and build artifacts.

- **docker/**  
  Contains Docker configuration files including the `Dockerfile` for building a container image of GollumFit and a script to start a Jupyter server within the container.

- **resources/**  
  Provides additional resources:
  - **AttenuationSplineMaker/**: Tools and settings for creating attenuation splines.
  - **correlation_matrices/**: Numpy files with flux and ice gradient correlations.
  - **DDMFluxMaker/**: Scripts and data files for generating covariance matrices from flux models.
  - **FluxOscCalculator/**: Tools for computing neutrino oscillation fluxes.
  - **Logo/**: Official logos and images for GollumFit.
  - **Splines/**: Directories with pre-computed splines for various parameters such as atmospheric kaon losses, ice gradients, and more.


## Contributing

Thank you for your interest in contributing to GollumFit! Please follow these steps:

1. **Fork & Clone**  
   Fork the repository on GitHub, then clone your fork locally.

2. **Create a Branch**  
   Create a feature branch using `git checkout -b feature/your-feature-name`

3. **Make Changes**

4. **Commit & Push**  
    Push using `git push origin feature/your-feature-name`

5. **Open a Pull Request**  
    Submit a pull request against the `main` branch and describe your changes.

6. **Reporting Issues**  
    Please open an issue if you find bugs or have feature requests.






