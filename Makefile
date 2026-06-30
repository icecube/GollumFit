# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))
DYN_OPT_PY=-shared -Wl,-soname,$(shell basename $(PYTHON_LIB))

VERSION=
PREFIX=/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local


PATH_GF=$(shell pwd)

SOURCES = $(PATH_GF)/src/analysisWeighting.cpp \
          $(PATH_GF)/src/GollumParameters.cpp \
          $(PATH_GF)/src/compactIO.cpp \
          $(PATH_GF)/src/Event.cpp \
          $(PATH_GF)/src/GollumTools.cpp \
          $(PATH_GF)/src/GollumEnumDefinitions.cpp \
		  $(PATH_GF)/src/utils.cpp \
		  $(PATH_GF)/src/GollumFit.cpp

HEADERS = $(PATH_GF)/include/analysisWeighting.h \
          $(PATH_GF)/include/GollumParameters.h \
          $(PATH_GF)/include/compactIO.h \
          $(PATH_GF)/include/Event.h \
          $(PATH_GF)/include/FastMode.h \
          $(PATH_GF)/include/GollumEnumDefinitions.h \
          $(PATH_GF)/include/GollumFit.h \
          $(PATH_GF)/include/GollumMCSet.h \
          $(PATH_GF)/include/GollumMCSpecifications.h \
		  $(PATH_GF)/include/GollumTools.h \
		  $(PATH_GF)/include/json.hpp \
		  $(PATH_GF)/include/utils.h




CFLAGS=-O0 -fPIC

CXXFLAGS +=-std=c++17 -O3

# Directories

SQUIDS_CFLAGS=-I/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/include -Wno-abi
SQUIDS_LDFLAGS=-L/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/lib -L/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/gsl-2.8-ojequbhgei7gix257tufrlkwsdrq3pia/lib -lSQuIDS -lgsl -lgslcblas -lm
NUSQUIDS_CFLAGS=-I/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/include
NUSQUIDS_LDFLAGS=-L/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/lib -lnuSQuIDS
LEPTONWEIGHTER_CFLAGS=-I/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/include
LEPTONWEIGHTER_LDFLAGS=-L/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/jvillarreal/higher-octant/software/local/lib -lLeptonWeighter
PHYSTOOLS_CFLAGS=
PHYSTOOLS_LDFLAGS=
BOOST_CFLAGS=-I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/include
BOOST_LDFLAGS=-Wl,-rpath -Wl,/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/lib -L/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/lib
PHOTOSPLINE_CFLAGS= -DPHOTOSPLINE_VERSION=2.4.1 -DPHOTOSPLINE_INCLUDES_SPGLAM -I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/cfitsio-4.4.1-jbrb54gbb7yqt27iakwft4diqo44rqe6/include -I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/include/suitesparse -I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/photospline-2.4.1-6t5vdztgab6yiz2hkass7ho6mfquej54/include
PHOTOSPLINE_LDFLAGS=-L/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/photospline-2.4.1-6t5vdztgab6yiz2hkass7ho6mfquej54/lib64 -lphotospline /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/cfitsio-4.4.1-jbrb54gbb7yqt27iakwft4diqo44rqe6/lib/libcfitsio.so -lspglam /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libcamd.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libccolamd.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libspqr.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libcholmod.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libamd.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libcolamd.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/suite-sparse-7.8.0-tngsi2da4ohorljomm37xeo42kmwq633/lib64/libsuitesparseconfig.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/openblas-0.3.27-az2vcj6xiqvl5ynofpmynzjsgjgzfiqo/lib/libopenblas.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/openblas-0.3.27-az2vcj6xiqvl5ynofpmynzjsgjgzfiqo/lib/libopenblas.so /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/openblas-0.3.27-az2vcj6xiqvl5ynofpmynzjsgjgzfiqo/lib/libopenblas.so -lpthread -lrt /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/metis-5.1.0-iku7m2rpcxmbnqja66wvohrbrgldhqxl/lib/libmetis.so
CFITSIO_CFLAGS=-I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/cfitsio-4.4.1-jbrb54gbb7yqt27iakwft4diqo44rqe6/include
CFITSIO_LDFLAGS=-L/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/cfitsio-4.4.1-jbrb54gbb7yqt27iakwft4diqo44rqe6/lib -lcfitsio
HDF5_CFLAGS=-I/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/hdf5-1.14.3-fryxdzmuzvvwtf33xt4iwq4xp2xewwms/include
HDF5_LDFLAGS=-L/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_8_x86_64_v2/spack/opt/spack/linux-centos8-x86_64_v2/gcc-13.3.0/hdf5-1.14.3-fryxdzmuzvvwtf33xt4iwq4xp2xewwms/lib -lhdf5


INC_GF=$(PATH_GF)/include
SRC_GF=$(PATH_GF)/src
LIB_GF=$(PATH_GF)/lib

# FLAGS
CXXFLAGS+= -fPIC -I$(INC_GF) $(SQUIDS_CFLAGS) $(NUSQUIDS_CFLAGS) $(LEPTONWEIGHTER_CFLAGS) $(PHYSTOOLS_CFLAGS) $(PHOTOSPLINE_CFLAGS) $(CFITSIO_CFLAGS) $(BOOST_CFLAGS) $(HDF5_CFLAGS)

LDFLAGS= -Wl,-rpath -Wl,$(LIB_GF) -L$(LIB_GF)
LDFLAGS+= $(NUSQUIDS_LDFLAGS) $(SQUIDS_LDFLAGS) $(LEPTONWEIGHTER_LDFLAGS) $(PHYSTOOLS_LDFLAGS) $(PHOTOSPLINE_LDFLAGS) $(CFITSIO_LDFLAGS) $(BOOST_LDFLAGS) $(HDF5_LDFLAGS)

CXXFLAGS+=-I/usr/local/include/
CXXFLAGS+=-I${PREFIX}/include/

LDFLAGS+=-L/usr/local/lib/
LDFLAGS+=-L${PREFIX}/lib/
LDFLAGS+=-L/usr/local/lib64/
LDFLAGS+=-L${PREFIX}/lib64/
ARCH := $(shell uname -m)
LDFLAGS+=-L/usr/lib/${ARCH}-linux-gnu/hdf5/serial/

LDFLAGS+=-lgsl -lgslcblas -lm
LDFLAGS+=-lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LDFLAGS+=-lhdf5 -lhdf5_hl -pthread
LDFLAGS+=-lphotospline -lcfitsio

LDFLAGS+=-lPhysTools
LDFLAGS+=-lSQuIDS -lnuSQuIDS
LDFLAGS+=-lLeptonWeighter

GOLLUM_HELPER_OBJECTS = Event.o analysisWeighting.o utils.o compactIO.o GollumFit.o GollumParameters.o GollumTools.o GollumEnumDefinitions.o

OS_NAME=$(shell uname -s)

ifeq (${OS_NAME},Linux)
DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))
endif

ifeq (${OS_NAME},Darwin)
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -compatibility_version $(VERSION) -current_version $(VERSION)
endif

# Project files
NAME=GollumFit
VERSION=1.0.0
STAT_PRODUCT:=$(LIB_GF)/lib$(NAME).a
DYN_PRODUCT:=$(LIB_GF)/lib$(NAME)$(DYN_SUFFIX)
OBJECTS = $(addprefix $(LIB_GF)/, $(GOLLUM_HELPER_OBJECTS))

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

.PHONY: all clean

clean:
	rm -rf $(OBJECTS) $(DYN_PRODUCT) $(STAT_PRODUCT) *.o

$(LIB_GF)/analysisWeighting.o : $(SRC_GF)/analysisWeighting.cpp $(INC_GF)/analysisWeighting.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/analysisWeighting.cpp -c -o $(LIB_GF)/analysisWeighting.o

$(LIB_GF)/GollumParameters.o : $(SRC_GF)/GollumParameters.cpp $(INC_GF)/GollumParameters.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/GollumParameters.cpp -c -o $(LIB_GF)/GollumParameters.o

$(LIB_GF)/compactIO.o : $(INC_GF)/compactIO.h $(SRC_GF)/compactIO.cpp $(INC_GF)/Event.h $(INC_GF)/analysisWeighting.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/compactIO.cpp -c -o $(LIB_GF)/compactIO.o

$(LIB_GF)/Event.o : $(SRC_GF)/Event.cpp $(INC_GF)/Event.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/Event.cpp -c -o $(LIB_GF)/Event.o

$(LIB_GF)/GollumTools.o : $(SRC_GF)/GollumTools.cpp $(INC_GF)/GollumTools.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/GollumTools.cpp -c -o $(LIB_GF)/GollumTools.o

$(LIB_GF)/GollumEnumDefinitions.o : $(SRC_GF)/GollumEnumDefinitions.cpp $(INC_GF)/GollumEnumDefinitions.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/GollumEnumDefinitions.cpp -c -o $(LIB_GF)/GollumEnumDefinitions.o

$(LIB_GF)/utils.o : $(SRC_GF)/utils.cpp $(INC_GF)/utils.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/utils.cpp -c -o $(LIB_GF)/utils.o

$(LIB_GF)/GollumFit.o : $(SRC_GF)/GollumFit.cpp $(INC_GF)/Event.h $(INC_GF)/analysisWeighting.h $(INC_GF)/compactIO.h $(INC_GF)/GollumFit.h $(INC_GF)/GollumParameters.h $(INC_GF)/GollumEnumDefinitions.h $(INC_GF)/FastMode.h
	$(CXX) $(CXXFLAGS) $(SRC_GF)/GollumFit.cpp -c -o $(LIB_GF)/GollumFit.o

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking dynamic library $(DYN_PRODUCT)
	@mkdir -p $(LIB_GF)/lib/
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking static library $(STAT_PRODUCT)
	@mkdir -p $(LIB_GF)/lib/
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/GollumFit
	@mkdir -p $(PREFIX)/include/GollumFit
	@cp $(HEADERS) $(PREFIX)/include/GollumFit
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp lib/gollumfit.pc $(PREFIX)/lib/pkgconfig


uninstall:
	@echo Removing headers from $(PREFIX)/include/GollumFit
	@rm -rf $(PREFIX)/include/GollumFit
	@echo Removing libraries from $(PREFIX)/lib
	@rm -f $(PREFIX)/$(DYN_PRODUCT)
	@rm -f $(PREFIX)/$(STAT_PRODUCT)
	@echo Removing config information from $(PREFIX)/lib/pkgconfig
	@rm -f $(PREFIX)/lib/pkgconfig/gollumfit.pc

