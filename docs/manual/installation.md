# Installation

## From Source

### Software prerequisites

The following necessary prerequisites are included as subdirectories in GollumFit already:
* SQuIDS (https://github.com/jsalvado/SQuIDS)
* nuSQuIDS (https://github.com/arguelles/nuSQuIDS)
* PhysTools (https://github.com/alexwenym/PhysTools.git)
* LeptonWeighter (https://github.com/alexwenym/LeptonWeighter.git)

The following tools and libraries must be available on your system. Where relevant, the version that was most recently checked to work is indicated:
* `g++` (9.2.0)
* `make` (4.2.1)
* `GSL`
* `boost` (1.7.1)
* `hdf5` 
* `pthread`
* `photospline` (https://github.com/cnweaver/photospline)
* `cfitsio` (https://heasarc.gsfc.nasa.gov/fitsio/)

To access these tools and libraries, IceCube users should use the `cvmfs`:
```
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/setup.sh`
```
Otherwise, they must be installed and available.

### Main Installation Procedure

First, set up the environment and create the necessary directories: 
```
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

Check out the necessary code, including necessary submodules:
```
cd $GOLLUMSOURCEPATH
git clone https://github.com/icecube/GollumFit
cd GollumFit
git submodule update --init --recursive
```

Compile the required software, one-by-one: 
```
# SQuIDS
cd $GOLLUMSOURCEPATH/GollumFit/vendor/SQuIDS
./configure --prefix=$PREFIX
make
make install

# nuSQuIDS
cd $GOLLUMSOURCEPATH/GollumFit/vendor/nuSQuIDS
./configure --prefix=$PREFIX --with-squids=$PREFIX --with-python-bindings --with-boost=${SROOT}
make
make python
make install

# PhysTools
cd $GOLLUMSOURCEPATH/GollumFit/vendor/PhysTools
mkdir build build/lbfgsb
./configure --with-boost=${SROOT} --prefix=$PREFIX
make
make install

# LeptonWeighter
cd $GOLLUMSOURCEPATH/GollumFit/vendor/LeptonWeighter
./configure --prefix=$PREFIX --with-boost=${SROOT}
make
make install
```

Now, compile GollumFit itself:
```
cd $GOLLUMSOURCEPATH/GollumFit
make
make install
```

Finally, compile the pybindings so that we can call functions and handle I/O in 
a python script. This is the most convienient way to access GollumFit functionalities: 

```
cd $GOLLUMSOURCEPATH/GollumFit/python
python setup.py install --prefix=$PREFIX
```

### Helper Tools

After building these, you may optionally build the various functions in the `$GOLLUMSOURCEPATH/GollumFit/resources` folder. 
For example, `FluxOscCalculator` can be built as follows:

```
cd $GOLLUMSOURCEPATH/GollumFit/resources/FluxOscCalculator
make
```

### Troubleshooting

## Containers

In addition to building from source, you may also choose to use a container.

The repository includes a `Dockerfile` that can be used to create both a Docker and a Singularity container. The Docker container can be used in environments where one has Admin priviliges (like your local machine), while the Singularity container should be used on, for example, HPC clusters.

The following section explains how to build a Docker image on your local machine and upload it to your personal `Docker Hub` account (which you can simply link to your Github account). The Singularity container can then be simply built on the cluster by pulling the Docker image from `Docker Hub`. Both the Docker container and the Singularity container can be used to start a Jupyter server inside the container, to which you then can connect to use GollumFit in Jupyter notebooks.

### Building the Docker container and uploading it to Docker Hub
1. Install [Docker](https://docs.docker.com/engine/install/) (I also recommend the desktop app, but its not needed).
2. Create a Docker Hub account [here](https://app.docker.com/signup) if you do not have one. You can use your Github or Google account. Make sure you have enough private repo slots available (the free plan offers only 1 private repo).
3. Go to your [Docker Hub repositories](https://hub.docker.com/repositories/) and create a repo. We have to set it to private as long as GollumFit is not yet public.
4. Clone the [GollumFit repo](https://github.com/icecube/GollumFit.git) to your local machine and `cd` into it.
5. Create a builder with `docker buildx create --name gollumfit-builder --use`
6. We want to build the container for `amd64` architectures (most clusters and non-Mac PCs) and push it to your Docker Hub. Optionally, as described further below, the container can additionally be built for `arm64` (most Mac PCs) architectures. The `amd64` build can be achieved with

    ```docker buildx build --platform linux/amd64 -f docker/Dockerfile -t <your_dockerhub_username>/<repo_name>:<tag> --push .```

    You can choose, for example, `gollumfit` as repo name and `latest` as tag. After executing the command, the images for both architectures can be seen in your [Docker Hub repository](https://hub.docker.com/repositories/). If you also want to run this container on your Mac (i.e. arm64 architecture), execute aboves command with `--platform linux/amd64,linux/arm64`. Note that we have experienced that the latter takes very long if the container is built on a non-Mac PC.

### Building the Singularity container
The Singularity container can be built from the Docker image that you uploaded to your Docker Hub in the previous step.

1. Check that Singularity is installed on your cluster via `singularity --version` or `apptainer --version` (singularity and apptainer are basically synonyms)
2. Add a read-only personal access token (PAT) to your Docker Hub account [here](https://app.docker.com/settings/personal-access-tokens).
3. Pull the Docker image and convert it to a singularity image file with `singularity pull --docker-login docker://<your_dockerhub_username>/<repo_name>:<tag>`

    You will be asked for your Docker Hub username and a password, for which you can take the PAT that you created in step 2. There might be warnings about the usage of undefined variables, which can be ignored.

### Starting a Jupyter server and connecting to it
1. The containers include a command `start-jupyter-server`, which executes the script that you can find in `docker/start-jupyter-server-in-container.sh`. Start the server with `singularity exec <sif_file> start-jupyter-server`
2. The output should include lines similar to `http://127.0.0.1:<port>/tree?token=<token>`.
3. If you use VS code, open a Jupyter notebook and click on `Select Kernel`. Then choose `Existing Jupyter Server` and copy/paste the URL from step 2.
4.  If you do not use VS code, you have to do port forwarding from your local machine to the cluster. This can be done by executing `ssh -L <port>:localhost:<port> <username>@<remote_host>` on your local machine, where `<port>` is the port found in the URL in step 2, and `<username>`/`<remote host>` are the username and the hostname in your ssh config (the one you use to connect to the cluster). Then you can open a browser and copy/paste the URL in step 2 into the address line.
