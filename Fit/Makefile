CFLAGS=-O0 -fPIC

CXXFLAGS+=-std=c++17 -fPIC
CXXFLAGS+=-I. -O3 -Wall -Wunused-but-set-variable -Wunused-variable

ifeq (${PREFIX},)
PREFIX=/usr/local
endif

ifeq (${SROOT},)
SROOT=/usr/local
endif
CXXFLAGS+=-I${SROOT}/include
LDFLAGS+=-L${SROOT}/lib
LDFLAGS+=-L${SROOT}/lib64

# Determine the architecture dynamically
# Necessary to find hdf5 libs for some versions
ARCH := $(shell uname -m)

ifeq (${GOLLUMBUILDPATH},)
GOLLUMBUILDPATH=/usr/local
endif
CXXFLAGS+=-I${GOLLUMBUILDPATH}/include/
CXXFLAGS+=-I${PREFIX}/include/
LDFLAGS+=-L${GOLLUMBUILDPATH}/lib/
LDFLAGS+=-L${PREFIX}/lib/
LDFLAGS+=-L${GOLLUMBUILDPATH}/lib64/
LDFLAGS+=-L${PREFIX}/lib64/
LDFLAGS+=-L/usr/lib/${ARCH}-linux-gnu/hdf5/serial/


LDFLAGS+=-lgsl -lgslcblas -lm
LDFLAGS+=-lboost_system -lboost_iostreams -lboost_filesystem -lboost_regex
LDFLAGS+=-lhdf5 -lhdf5_hl
LDFLAGS+=-pthread

# the following should be found in the include and lib directories
# of the SNOTBUILDPATH
LDFLAGS+=-lPhysTools
LDFLAGS+=-lSQuIDS -lnuSQuIDS
LDFLAGS+=-lLeptonWeighter
LDFLAGS+=-lphotospline -lcfitsio

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
STAT_PRODUCT=../lib/lib$(NAME).a
DYN_PRODUCT=../lib/lib$(NAME)$(DYN_SUFFIX)
OBJECTS= ${GOLLUM_HELPER_OBJECTS}

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

.PHONY: all clean

clean:
	rm -rf $(GOLLUM_HELPER_OBJECTS) $(DYN_PRODUCT) $(STAT_PRODUCT) *.o

analysisWeighting.o : analysisWeighting.cpp analysisWeighting.h
	$(CXX) $(CXXFLAGS) analysisWeighting.cpp -c -o analysisWeighting.o

GollumParameters.o : GollumParameters.cpp GollumParameters.h
	$(CXX) $(CXXFLAGS) GollumParameters.cpp -c -o GollumParameters.o

compactIO.o : compactIO.h compactIO.cpp Event.h analysisWeighting.h
	$(CXX) $(CXXFLAGS) compactIO.cpp -c -o compactIO.o

Event.o : Event.cpp Event.h
	$(CXX) $(CXXFLAGS) Event.cpp -c -o Event.o

GollumTools.o : GollumTools.cpp GollumTools.h
	$(CXX) $(CXXFLAGS) GollumTools.cpp -c -o GollumTools.o

GollumEnumDefinitions.o : GollumEnumDefinitions.cpp GollumEnumDefinitions.h
	$(CXX) $(CXXFLAGS) GollumEnumDefinitions.cpp -c -o GollumEnumDefinitions.o

utils.o : utils.cpp utils.h
	$(CXX) $(CXXFLAGS) utils.cpp -c -o utils.o

GollumFit.o : GollumFit.cpp Event.h analysisWeighting.h compactIO.h GollumFit.h GollumParameters.h GollumEnumDefinitions.h FastMode.h
	$(CXX) $(CXXFLAGS) GollumFit.cpp -c -o GollumFit.o

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking dynamic library $(DYN_PRODUCT)
	@mkdir -p ../lib/
	@$(CXX) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking static library $(STAT_PRODUCT)
	@mkdir -p ../lib/
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)

install: $(DYN_PRODUCT) $(STAT_PRODUCT)
	@echo Installing headers in $(PREFIX)/include/GollumFit
	@mkdir -p $(PREFIX)/include/GollumFit
	@cp ./*.h $(PREFIX)/include/GollumFit
	@cp ./*.hpp $(PREFIX)/include/GollumFit
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp gollumfit.pc $(PREFIX)/lib/pkgconfig

