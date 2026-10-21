# GNUmakefile for tgv2d benchmark (benchmarks/tgv2d/).
# Uses Make.muphfasa from external/muphfasa, with local InitData.cpp
# (Taylor-Green) overriding the default quiescent start via VPATH.

TOP := ../../external/muphfasa

DIM ?= 2

EBASE = bench

CEXE_sources += ../bench.cpp
# Local InitData.cpp and TgvAnalytic.H provide TGV IC (VPATH precedence)
CEXE_sources += InitData.cpp
VPATH_LOCATIONS += ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa

# Ensure local TgvAnalytic.H is found
INCLUDE_LOCATIONS += .
