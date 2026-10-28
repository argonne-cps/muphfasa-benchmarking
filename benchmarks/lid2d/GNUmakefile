# GNUmakefile for lid-driven cavity 2D benchmark (benchmarks/lid2d/).
# Closed square cavity: NoSlip walls on x-lo/x-hi/y-lo, MovingWall lid at y-hi.

TOP := ../../external/muphfasa

DIM ?= 2

EBASE = bench

CEXE_sources += ../bench.cpp
# No custom InitData.cpp: uses default quiescent start, lid velocity via BC
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa
