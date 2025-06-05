SRC_DIR = src
INC_DIR = include

CFLAGS=-O0 -fPIC

CXXFLAGS+=-std=c++17 -fPIC
CXXFLAGS+=-I$(INC_DIR) -O3 -Wall -Wunused-but-set-variable -Wunused-variable

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
OBJ_DIR = ../lib
OBJECTS = $(addprefix $(OBJ_DIR)/, $(GOLLUM_HELPER_OBJECTS))

# Compilation rules
all: $(STAT_PRODUCT) $(DYN_PRODUCT)

.PHONY: all clean

clean:
	rm -rf $(OBJECTS) $(DYN_PRODUCT) $(STAT_PRODUCT) *.o

$(OBJ_DIR)/analysisWeighting.o : $(SRC_DIR)/analysisWeighting.cpp $(INC_DIR)/analysisWeighting.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/analysisWeighting.cpp -c -o $(OBJ_DIR)/analysisWeighting.o

$(OBJ_DIR)/GollumParameters.o : $(SRC_DIR)/GollumParameters.cpp $(INC_DIR)/GollumParameters.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/GollumParameters.cpp -c -o $(OBJ_DIR)/GollumParameters.o

$(OBJ_DIR)/compactIO.o : $(INC_DIR)/compactIO.h $(SRC_DIR)/compactIO.cpp $(INC_DIR)/Event.h $(INC_DIR)/analysisWeighting.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/compactIO.cpp -c -o $(OBJ_DIR)/compactIO.o

$(OBJ_DIR)/Event.o : $(SRC_DIR)/Event.cpp $(INC_DIR)/Event.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/Event.cpp -c -o $(OBJ_DIR)/Event.o

$(OBJ_DIR)/GollumTools.o : $(SRC_DIR)/GollumTools.cpp $(INC_DIR)/GollumTools.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/GollumTools.cpp -c -o $(OBJ_DIR)/GollumTools.o

$(OBJ_DIR)/GollumEnumDefinitions.o : $(SRC_DIR)/GollumEnumDefinitions.cpp $(INC_DIR)/GollumEnumDefinitions.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/GollumEnumDefinitions.cpp -c -o $(OBJ_DIR)/GollumEnumDefinitions.o

$(OBJ_DIR)/utils.o : $(SRC_DIR)/utils.cpp $(INC_DIR)/utils.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/utils.cpp -c -o $(OBJ_DIR)/utils.o

$(OBJ_DIR)/GollumFit.o : $(SRC_DIR)/GollumFit.cpp $(INC_DIR)/Event.h $(INC_DIR)/analysisWeighting.h $(INC_DIR)/compactIO.h $(INC_DIR)/GollumFit.h $(INC_DIR)/GollumParameters.h $(INC_DIR)/GollumEnumDefinitions.h $(INC_DIR)/FastMode.h
	$(CXX) $(CXXFLAGS) $(SRC_DIR)/GollumFit.cpp -c -o $(OBJ_DIR)/GollumFit.o

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
	@cp $(INC_DIR)/*.h $(PREFIX)/include/GollumFit
	@cp $(INC_DIR)/*.hpp $(PREFIX)/include/GollumFit
	@echo Installing libraries in $(PREFIX)/lib
	@mkdir -p $(PREFIX)/lib
	@cp $(DYN_PRODUCT) $(STAT_PRODUCT) $(PREFIX)/lib
	@echo Installing config information in $(PREFIX)/lib/pkgconfig
	@mkdir -p $(PREFIX)/lib/pkgconfig
	@cp $(SRC_DIR)/gollumfit.pc $(PREFIX)/lib/pkgconfig

