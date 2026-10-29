# GNUmakefile for the lid-driven cavity 3-D correctness test (benchmarks/lid3d/validate/).
# Closed cubic cavity: NoSlip on all faces except MovingWall lid at z-hi (x-direction).
# Uses the default quiescent InitData.cpp from src/Simulation/InitData.cpp (cavity
# starts from rest); no custom InitData.cpp. See README.md and inputs.

TOP := ../../../external/muphfasa

DIM ?= 3

EBASE = test_lid3d

CEXE_sources += test_lid3d.cpp

include $(TOP)/Make.muphfasa
