CXXFLAGS+=-O3 -std=c++11 -fPIC -pthread -I.

ifeq (${PREFIX},)
PREFIX=/usr/local
endif

ifeq (${SROOT},)
SROOT=/usr/local
endif

HDF5_CFLAGS=-I${SROOT}/include
HDF5_LDFLAGS=-L${SROOT}/lib -lhdf5 -lhdf5_hl
LDFLAGS+= $(HDF5_LDFLAGS)

CXXFLAGS+=-I${SROOT}/include
CXXFLAGS+=-I${PREFIX}/include
CXXFLAGS+=`pkg-config nusquids --cflags`

LDFLAGS+=-L${SROOT}/lib
LDFLAGS+=-L${SROOT}/lib64

LDFLAGS+=-L${PREFIX}/lib
LDFLAGS+=-L${PREFIX}/lib64

LDFLAGS+=-L${GOLLUMBUILDPATH}/lib
LDFLAGS+=-L${GOLLUMBUILDPATH}/lib64

LDFLAGS+=`pkg-config nusquids --libs`

# Project files
SOURCES = $(wildcard ./*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: all clean

# Compilation rules
all: astrophysical_neutrino_flux_calculator_with_interactions conventional_atmospheric_flux_calculator_with_interactions_from_file prompt_atmospheric_flux_calculator_with_interactions_from_file errors_ddm_flux_calculator_with_interactions_from_file

astrophysical_neutrino_flux_calculator_with_interactions: astrophysical_neutrino_flux_calculator_with_interactions.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) astrophysical_neutrino_flux_calculator_with_interactions.o -o astrophysical_neutrino_flux_calculator_with_interactions

conventional_atmospheric_flux_calculator_with_interactions_from_file: conventional_atmospheric_flux_calculator_with_interactions_from_file.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) conventional_atmospheric_flux_calculator_with_interactions_from_file.o -o conventional_atmospheric_flux_calculator_with_interactions_from_file

prompt_atmospheric_flux_calculator_with_interactions_from_file: prompt_atmospheric_flux_calculator_with_interactions_from_file.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) prompt_atmospheric_flux_calculator_with_interactions_from_file.o -o prompt_atmospheric_flux_calculator_with_interactions_from_file

errors_ddm_flux_calculator_with_interactions_from_file: errors_ddm_flux_calculator_with_interactions_from_file.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) errors_ddm_flux_calculator_with_interactions_from_file.o -o errors_ddm_flux_calculator_with_interactions_from_file

install: astrophysical_neutrino_flux_calculator_with_interactions conventional_atmospheric_flux_calculator_with_interactions_from_file prompt_atmospheric_flux_calculator_with_interactions_from_file errors_ddm_flux_calculator_with_interactions_from_file
	@echo Installing binaries in $(PREFIX)/bin
	@mkdir -p $(PREFIX)/bin
	@cp conventional_atmospheric_flux_calculator_with_interactions_from_file $(PREFIX)/bin/conventional_atmospheric_flux_calculator_with_interactions_from_file
	@cp prompt_atmospheric_flux_calculator_with_interactions_from_file $(PREFIX)/bin/prompt_atmospheric_flux_calculator_with_interactions_from_file
	@cp astrophysical_neutrino_flux_calculator_with_interactions $(PREFIX)/bin/astrophysical_neutrino_flux_calculator_with_interactions
	@cp errors_ddm_flux_calculator_with_interactions_from_file $(PREFIX)/bin/errors_ddm_flux_calculator_with_interactions_from_file

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f *.o conventional_atmospheric_flux_calculator_with_interactions_from_file $(PREFIX)/bin/conventional_atmospheric_flux_calculator_with_interactions_from_file
	rm -f *.o prompt_atmospheric_flux_calculator_with_interactions_from_file $(PREFIX)/bin/prompt_atmospheric_flux_calculator_with_interactions_from_file
	rm -f *.o astrophysical_neutrino_flux_calculator_with_interactions $(PREFIX)/bin/astrophysical_neutrino_flux_calculator_with_interactions
	rm -f *.o errors_ddm_flux_calculator_with_interactions_from_file $(PREFIX)/bin/errors_ddm_flux_calculator_with_interactions_from_file

