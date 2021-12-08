cmake_minimum_required(VERSION $ENV{ACC_CMAKE_VERSION})

set(LIBNAME phase_advance_4py)

set(SHARED_DEPS 
  bmad
  sim_utils
  pgplot
  X11
  readline
  ncurses
  xsif
  recipes_f-90_LEPP
  forest
  fgsl
  gsl
  gslcblas
  xrlf03
  xrl
  lapack95
  lapack
  blas
)

set (SRC_DIRS
  modules
)

set(CREATE_SHARED true)
#
#SET (FFLAGS "-fPIC")
#SET (LINK_FLAGS "-fPIC")

#set (EXE_SPECS
# cmake.phase_advance
#)
##
include($ENV{ACC_BUILD_SYSTEM}/Master.cmake)
