# GNUmakefile for Beltrami 3D benchmark (benchmarks/beltrami3d/).
# Triply periodic, exercises all velocity components and D3Q19/27 isotropy.

TOP := ../../external/muphfasa

DIM ?= 3

EBASE = bench

CEXE_sources += ../bench.cpp
CEXE_sources += InitData.cpp
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa

INCLUDE_LOCATIONS += .
