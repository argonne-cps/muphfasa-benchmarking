# GNUmakefile for tgv3d benchmark (benchmarks/tgv3d/).
# Triply periodic 3-D Taylor-Green vortex, D3Q19/27.

TOP := ../../external/muphfasa

DIM ?= 3

EBASE = bench

CEXE_sources += ../bench.cpp
CEXE_sources += InitData.cpp
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa

INCLUDE_LOCATIONS += .
