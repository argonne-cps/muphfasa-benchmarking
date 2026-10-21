# GNUmakefile for lid-driven cavity 3D benchmark (benchmarks/lid3d/).
# Wall-normal z (like couette3d): periodic in x,y, walls in z (NoSlip/MovingWall).

TOP := ../../external/muphfasa

DIM ?= 3

EBASE = bench

CEXE_sources += ../bench.cpp
# No custom InitData.cpp: uses default quiescent start, lid velocity via BC
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa
