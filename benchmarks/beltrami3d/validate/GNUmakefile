# GNUmakefile for the Beltrami 3-D correctness test (benchmarks/beltrami3d/validate/).
# Mirrors external/muphfasa/tests/regression/tgv2d/GNUmakefile: builds
# test_beltrami3d directly, without going through CMake.
#
# Usage:
#   make                       # build test_beltrami3d (3-D, no MPI, by default)
#   ./test_beltrami3d3d.gnu.ex inputs
#
# InitData.cpp and BeltramiAnalytic.H live one directory up
# (benchmarks/beltrami3d/, shared with the bench.cpp perf driver); VPATH
# picks up InitData.cpp from there (overriding src/Simulation/InitData.cpp's
# default quiescent start) and INCLUDE_LOCATIONS lets this file #include
# "BeltramiAnalytic.H" from there too.

TOP := ../../../external/muphfasa

DIM ?= 3

EBASE = test_beltrami3d

CEXE_sources += test_beltrami3d.cpp

# Narrow vpath for InitData.cpp only (not a blanket VPATH_LOCATIONS += ..):
# AMReX's Make.rules generates tmp_build_dir/s/<dim>d.gnu.EXE/AMReX_Config.H
# as a relative-path target, and a blanket VPATH entry for ".." makes GNU
# Make match that same relative path in ../tmp_build_dir/s/<dim>d.gnu.EXE/
# (beltrami3d's own perf-bench build dir, one level up) and consider the
# header "up to date" there without ever generating/copying it locally,
# so the actual compile (which only has -Itmp_build_dir/...) fails to find
# it. vpath for the single named source file sidesteps that collision.
vpath InitData.cpp ..
INCLUDE_LOCATIONS += ..

include $(TOP)/Make.muphfasa
