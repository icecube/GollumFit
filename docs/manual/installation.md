# Installation {#installation}

## Building from Source

### Prerequisites

#### Included Dependencies

The following dependencies are included as git submodules in the GollumFit repository:
- [SQuIDS](https://github.com/jsalvado/SQuIDS) - Solve Quantum Systems
- [nuSQuIDS](https://github.com/arguelles/nuSQuIDS) - Neutrino oscillation solver
- [PhysTools](https://github.com/icecube/PhysTools) - Physics utilities including LBFGSB
- [LeptonWeighter](https://github.com/alexwenym/LeptonWeighter.git) - Lepton weighting utilities

#### Required System Dependencies

The following tools and libraries must be installed on your system (tested versions in parentheses):

**Compilers and Build Tools:**
- `g++` (9.2.0 or later)
- `make` (4.2.1 or later)

**Scientific Libraries:**
- <a href="https://www.gnu.org/software/gsl/"><code>GSL</code></a> - GNU Scientific Library
- `boost` (1.71.0 or later)
- `hdf5` - Hierarchical Data Format 5
- `pthread` - POSIX threads
- <a href="https://github.com/cnweaver/photospline"><code>photospline</code></a> - Photo-spline fitting
- <a href="https://heasarc.gsfc.nasa.gov/fitsio/"><code>cfitsio</code></a> - FITS file I/O
- `pybind11` - Python bindings (highly recommended)

**For IceCube Users:**

Access all required dependencies via CVMFS:
```bash
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
```

**For Other Users:**

All dependencies must be installed manually and available in your system paths.

### Installation Steps

#### 1. Environment Setup

First, set up the environment variables and create the necessary directories: 
```bash
export CXX=g++
export CC=gcc
export GOLLUMSPACE=/choose/your/favourite/dir
export GOLLUMBUILDPATH=$GOLLUMSPACE/local
export GOLLUMSOURCEPATH=$GOLLUMSPACE
export PREFIX=$GOLLUMBUILDPATH
export PATH=$GOLLUMBUILDPATH/bin:$PATH
export LD_LIBRARY_PATH=$GOLLUMBUILDPATH/lib/:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=$GOLLUMBUILDPATH/include/:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$GOLLUMBUILDPATH/include/:$CPLUS_INCLUDE_PATH
export CXX_INCLUDE_PATH=$GOLLUMBUILDPATH/include/:$CXX_INCLUDE_PATH
export PKG_CONFIG_PATH=$GOLLUMBUILDPATH/lib/pkgconfig:$PKG_CONFIG_PATH
export PYTHONPATH=$GOLLUMBUILDPATH/lib/python3.7/site-packages:$PYTHONPATH
export HDF5_DISABLE_VERSION_CHECK=1
mkdir $GOLLUMBUILDPATH
mkdir $PREFIX/bin $PREFIX/include $PREFIX/lib $PREFIX/lib64
```
**Note:** Ensure the Boost library and header directories are correctly set in `$BOOST_DIR`.

#### 2. Clone Repository

Check out the GollumFit repository and initialize submodules:

```bash
cd $GOLLUMSOURCEPATH
git clone https://github.com/icecube/GollumFit
cd GollumFit
git submodule update --init --recursive
```

#### 3. Build Dependencies

Compile each dependency in order:

**SQuIDS:**
```bash
# SQuIDS
cd $GOLLUMSOURCEPATH/GollumFit/vendor/SQuIDS
./configure --prefix=$PREFIX
make
make install
```

**nuSQuIDS:**
```bash
# nuSQuIDS
cd $GOLLUMSOURCEPATH/GollumFit/vendor/nuSQuIDS
./configure --prefix=$PREFIX --with-squids=$PREFIX --with-python-bindings --with-boost=$BOOST_DIR
make
make python
make install
```

**PhysTools:**
```bash
# PhysTools
cd $GOLLUMSOURCEPATH/GollumFit/vendor/PhysTools
mkdir build build/lbfgsb
./configure --with-boost=$BOOST_DIR --prefix=$PREFIX
make
make install
```

**LeptonWeighter:**
```bash
# LeptonWeighter
cd $GOLLUMSOURCEPATH/GollumFit/vendor/LeptonWeighter
./configure --prefix=$PREFIX --with-boost=$BOOST_DIR
make
make install
```

#### 4. Build GollumFit

Compile the main GollumFit library:

```bash
cd $GOLLUMSOURCEPATH/GollumFit
./configure --prefix=$PREFIX --with-boost=$BOOST_DIR --with-nusquids=$PREFIX --with-leptonweighter=$PREFIX
make
make install
```

#### 5. Build Python Bindings

Compile the Python bindings (recommended for convenient access to GollumFit functionalities):

```bash
cd $GOLLUMSOURCEPATH/GollumFit/python
python setup.py install --prefix=$PREFIX
```

#### 6. Optional Helper Tools

After the main build, you may optionally compile various utility tools in the resources folder. For example, to build the flux oscillation calculator:

```bash
cd $GOLLUMSOURCEPATH/GollumFit/resources/FluxOscCalculator
make
```

Other available tools include:
- `AttenuationSplineMaker` - Generate attenuation splines
- `DDMFluxMaker` - Create DaemonFlux fluxes and covariance matrices
- `FluxOscCalculator` - Calculate oscillated astrophysical neutrino fluxes

### Troubleshooting

(This section is reserved for common installation issues and their solutions.)

## Using Containers

As an alternative to building from source, GollumFit can be run in a containerized environment.

### Overview

The repository includes a [Dockerfile](../../docker/Dockerfile) that can be used to create both Docker and Singularity containers:
- **Docker**: For local machines where you have admin privileges
- **Singularity/Apptainer**: For HPC clusters and shared computing environments

The workflow involves:
1. Building a Docker image on your local machine
2. Uploading to Docker Hub
3. Pulling to create a Singularity container on the cluster
4. Running Jupyter notebooks inside the container

### Docker Container Setup

#### Building and Uploading to Docker Hub
1. **Install Docker**: Download and install [Docker](https://docs.docker.com/engine/install/) (Docker Desktop recommended but not required).

2. **Create Docker Hub Account**: Sign up at [Docker Hub](https://app.docker.com/signup) using your GitHub or Google account. Ensure you have available repository slots (free plan includes 1 private repo).

3. **Create Repository**: Navigate to your [Docker Hub repositories](https://hub.docker.com/repositories/) and create a new repository.

4. **Clone GollumFit**: Clone the repository to your local machine:
   ```bash
   git clone https://github.com/icecube/GollumFit.git
   cd GollumFit
   ```

5. **Create Builder**: Set up a multi-platform builder:
   ```bash
   docker buildx create --name gollumfit-builder --use
   ```

6. **Build and Push**: Build for `amd64` architecture (most clusters and non-Mac PCs):
   ```bash
   docker buildx build --platform linux/amd64 \
     -f docker/Dockerfile \
     -t <your_dockerhub_username>/<repo_name>:<tag> \
     --push .
   ```
   
   Example: `gollumfit` as repo name and `latest` as tag.
   
   **For Mac users**: To also support `arm64` architecture:
   ```bash
   docker buildx build --platform linux/amd64,linux/arm64 \
     -f docker/Dockerfile \
     -t <your_dockerhub_username>/<repo_name>:<tag> \
     --push .
   ```
   
   **Note**: Building for multiple platforms may take significantly longer on non-Mac systems.

### Singularity Container Setup

Build a Singularity/Apptainer container from your Docker image for use on HPC clusters.

1. **Verify Singularity Installation**: Check that Singularity or Apptainer is available:
   ```bash
   singularity --version
   # or
   apptainer --version
   ```
   
   **Note**: Singularity and Apptainer are equivalent.

2. **Create Access Token**: Generate a read-only Personal Access Token (PAT) at [Docker Hub Settings](https://app.docker.com/settings/personal-access-tokens).

3. **Pull and Convert Image**:
   ```bash
   singularity pull --docker-login docker://<your_dockerhub_username>/<repo_name>:<tag>
   ```
   
   When prompted:
   - Username: Your Docker Hub username
   - Password: The PAT from step 2
   
   **Note**: You may see warnings about undefined variables; these can be safely ignored.

### Running Jupyter in the Container

#### Start Jupyter Server

The containers include a convenience command to launch Jupyter:

```bash
singularity exec <sif_file> start-jupyter-server
```

This executes the script [docker/start-jupyter-server-in-container.sh](../../docker/start-jupyter-server-in-container.sh).

The output will include a URL similar to:
```
http://127.0.0.1:<port>/tree?token=<token>
```

#### Connect to Jupyter

**Option 1: VS Code (Recommended)**

1. Open a Jupyter notebook in VS Code
2. Click `Select Kernel`
3. Choose `Existing Jupyter Server`
4. Paste the URL from the server output

**Option 2: Browser with SSH Port Forwarding**

1. On your local machine, set up port forwarding:
   ```bash
   ssh -L <port>:localhost:<port> <username>@<remote_host>
   ```
   
   Where:
   - `<port>` is from the Jupyter URL
   - `<username>@<remote_host>` matches your SSH config

2. Open a browser and navigate to the URL from the server output
