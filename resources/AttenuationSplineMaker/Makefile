CXXFLAGS+=-O3 -std=c++11 -fPIC -pthread -I.

ifeq (${PREFIX},)
PREFIX=/usr/local
endif

ifeq (${SROOT},)
SROOT=/usr/local
endif

CXXFLAGS+=-I${SROOT}/include
CXXFLAGS+=-I${PREFIX}/include
CXXFLAGS+=`pkg-config nusquids --cflags`
CXXFLAGS+=`photospline-config --cflags`

LDFLAGS+=-L${SROOT}/lib
LDFLAGS+=-L${SROOT}/lib64

LDFLAGS+=-L${PREFIX}/lib
LDFLAGS+=-L${PREFIX}/lib64

LDFLAGS+=-L${GOLLUMBUILDPATH}/lib
LDFLAGS+=-L${GOLLUMBUILDPATH}/lib64

LDFLAGS+=`pkg-config nusquids --libs`
CXXFLAGS+=`photospline-config --libs`

# Project files
SOURCES = $(wildcard ./*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Compilation rules
all: conventional_attenuation_spit prompt_attenuation_spit astro_attenuation_spit

.PHONY: all clean

clean:
	rm -rf conventional_attenuation_spit prompt_attenuation_spit astro_attenuation_spit

conventional_attenuation_spit: conventional_attenuation_spit.cpp
	$(CXX) $(CXXFLAGS) conventional_attenuation_spit.cpp $(LDFLAGS) -o conventional_attenuation_spit

prompt_attenuation_spit: prompt_attenuation_spit.cpp
	$(CXX) $(CXXFLAGS) prompt_attenuation_spit.cpp $(LDFLAGS) -o prompt_attenuation_spit

astro_attenuation_spit: astro_attenuation_spit.cpp
	$(CXX) $(CXXFLAGS) astro_attenuation_spit.cpp $(LDFLAGS) -o astro_attenuation_spit
