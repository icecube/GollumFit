# Use the official Ubuntu base image
FROM ubuntu:20.04

# Set Dockerfile author
LABEL maintainer="Lukas Hennig <lukas.hennig@fau.de>"

# Set non-interactive frontend for automatic installation
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    gcc \
    wget \
    python3 \
    python3-pip \
    python3-dev \
    git \
    openssh-client \
    pkg-config \
    libgsl-dev \
    libboost-all-dev \
    libboost-python-dev \
    libboost-python1.71.0 \
    libhdf5-dev \
    hdf5-tools \
    cmake \
    libcfitsio-dev \
    libblas-dev \
    liblapack-dev \
    libsuitesparse-dev \
    python3-sphinx \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && ARCH=$(uname -m) \
    && ln -sf /usr/lib/${ARCH}-linux-gnu/libpython3.8.so /usr/lib/${ARCH}-linux-gnu/libpython3.so \
    && pip install --upgrade pip \
    && pip install numpy h5py scipy pandas uproot toml matplotlib \
    && pip install notebook

ENV CXX=g++
ENV CC=gcc
ENV GOLLUMSPACE=/gollumspace
ENV GOLLUMBUILDPATH=/usr/local
ENV GOLLUMSOURCEPATH=$GOLLUMSPACE/GollumFit/sources
ENV PREFIX=$GOLLUMBUILDPATH
ENV PATH=$GOLLUMBUILDPATH/bin:$PATH
ENV LD_LIBRARY_PATH=$GOLLUMBUILDPATH/lib/:$LD_LIBRARY_PATH
ENV C_INCLUDE_PATH=$GOLLUMBUILDPATH/include/:$C_INCLUDE_PATH
ENV CPLUS_INCLUDE_PATH=/usr/include/hdf5/serial/:$GOLLUMBUILDPATH/include/:$CPLUS_INCLUDE_PATH
ENV CXX_INCLUDE_PATH=$GOLLUMBUILDPATH/include/:$CXX_INCLUDE_PATH
ENV PKG_CONFIG_PATH=$GOLLUMBUILDPATH/lib/pkgconfig:$PKG_CONFIG_PATH
ENV PYTHONPATH=$GOLLUMBUILDPATH/lib/python3.8/site-packages:$PYTHONPATH
ENV HDF5_DISABLE_VERSION_CHECK=1

RUN mkdir -p $GOLLUMSPACE $GOLLUMBUILDPATH $GOLLUMSOURCEPATH && \
    mkdir -p $PREFIX/bin $PREFIX/include $PREFIX/lib $PREFIX/lib64

ADD ../docs $GOLLUMSPACE/GollumFit/docs/
ADD ../Fit $GOLLUMSPACE/GollumFit/Fit/
ADD ../python $GOLLUMSPACE/GollumFit/python/
ADD ../sources $GOLLUMSPACE/GollumFit/sources/

# For most cases, one does not want the following directories to be
# included. It is better to bind them in whenever it is necessary
# ADD data $GOLLUMSPACE/GollumFit/data/
# ADD monte_carlo $GOLLUMSPACE/GollumFit/monte_carlo/
# ADD resources $GOLLUMSPACE/GollumFit/resources/

WORKDIR $GOLLUMSOURCEPATH/SQuIDS
RUN ./configure --prefix=/usr && make && make install


WORKDIR $GOLLUMSOURCEPATH/nuSQuIDS
RUN ./configure --prefix=/usr --with-python-bindings \
    && make \
    && make test \
    && make python \
    && make install


WORKDIR $GOLLUMSOURCEPATH/PhysTools
RUN mkdir build build/lbfgsb \
    && ARCH=$(uname -m) \
    && ./configure --prefix=/usr --with-boost-incdir=/usr/include/ --with-boost-libdir=/usr/lib/${ARCH}-linux-gnu \
    && make \
    && make install


WORKDIR $GOLLUMSPACE
RUN git clone https://github.com/cnweaver/photospline.git photospline
RUN mkdir $GOLLUMSPACE/photospline/build
WORKDIR $GOLLUMSPACE/photospline/build
RUN cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr ../ \
    && make \
    && make install


WORKDIR $GOLLUMSOURCEPATH/LeptonWeighter
RUN ARCH=$(uname -m) \
    && ./configure --prefix=/usr --with-boost-incdir=/usr/include/ --with-boost-libdir=/usr/lib/${ARCH}-linux-gnu \
    && make \
    && make install


WORKDIR $GOLLUMSPACE/GollumFit/Fit
RUN make && make install
WORKDIR $GOLLUMSPACE/GollumFit/python
RUN python setup.py install --prefix=$PREFIX \
    && python3 -c "import GollumFitPy"

# Expose the Jupyter server port
EXPOSE 8888

# Add a custom script for starting the Jupyter server
ADD docker/start-jupyter-server-in-container.sh /usr/local/bin/start-jupyter-server
RUN chmod +x /usr/local/bin/start-jupyter-server
    
WORKDIR /workspace