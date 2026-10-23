# GNUmakefile for the TGV 3-D correctness test (benchmarks/tgv3d/validate/).
# Mirrors beltrami3d/validate/GNUmakefile.

TOP := ../../../external/muphfasa

DIM ?= 3

EBASE = test_tgv3d

CEXE_sources += test_tgv3d.cpp

vpath InitData.cpp ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa
