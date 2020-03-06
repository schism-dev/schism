################################################################################
# Parallel SCHISM Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below as well as setting flags having to do with included algorithms (e.g. sediment)
# and the compiler configuration (debug, timing).
#
# The environment settings are based on the following options.
#
# Compiler name:
#   FCS: Serial compiler (for utilities)
#   FCP: Parallel compiler
#   FLD: Linker (in general same as parallel compiler)
#
# Compilation flags
#   FCSFLAGS: Flags for serial compilation
#   FCPFLAGS: Flags for parallel compilation (including all pre-processing flags)
#   FLDFLAGS: Flags for linker (e.g., -O2)
#
# Preprocessor flags:
#   DEBUG: Enable debugging code
#   ORDERED_SUM: Enable globally ordered sums & dot-products for bit reproducibility
#     of state quantities independent of number of processors (note: this can
#     significantly degrade performance);
#   INCLUDE_TIMING: Enable wallclock timing of code (note: this can have slight
#     effect on performance);
#   MPI_VERSION = 1 or 2: Version of MPI (try 2 first, if compile fails due to mpi
#     related errors then switch to version 1;
#
# Libraries (needed for parallel code)
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   ALTLIBS: Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
#
################################################################################

# Environment for Intel Compiler / ifort (Thor cluster)
ENV = thor.intel


################################################################################
# Fortran compilers
FCS = f90
FCP = mpiifort -f90=ifort
FLD = $(FCP)


################################################################################
# MPI version (1 or 2)
PPFLAGS := $(PPFLAGS) -DMPIVERSION=2


################################################################################
# Optimization and linking flags for debugging and development
#OPTFLAGS = -O2 -xHost -fPIC -g -fp-model precise -traceback -assume byterecl -check uninit -check bounds -check pointers -warn all,nounused -gen-interfaces
# Optimization flags for higher performances
OPTFLAGS = -O2 -xHost -fPIC -g -fp-model precise
FCPFLAGS = $(PPFLAGS) $(OPTFLAGS) -lnetcdff
FLDFLAGS = $(OPTFLAGS)


################################################################################
# Parmetis
METIS_PATH=../src/ParMetis-3.1-Sep2010/
MTSLIBS = -L$(METIS_PATH) -lparmetis -lmetis


################################################################################
# Netcdf libraries (Be consistent with your .bashrc !)
# Installed on Thor
#export NETCDF=/sw/lib/netcdf/4.5.0
#export NETCDF_FORTRAN=/sw/lib/netcdf/ifort/4.4.4_linked_with_4.5.0
#CDFLIBS =-L${NETCDF_FORTRAN}/lib -lnetcdff -L${NETCDF}/lib -lnetcdf
#CDFMOD  =-I${NETCDF}/include -I${NETCDF_FORTRAN}/include

# Installed locally on Xavier's space
export NETCDF=/home/xbertin/Partage_DPL/libraries/netcdfC
export NETCDF_FORTRAN=/home/xbertin/Partage_DPL/libraries/netcdfF
CDFLIBS =-L${NETCDF_FORTRAN}/lib -lnetcdff -L${NETCDF}/lib -lnetcdf
CDFMOD  =-I${NETCDF}/include -I${NETCDF_FORTRAN}/include


################################################################################
# Alternate executable name if you do not want the default.
#EXEC   := SCHISM_TRUNK_KM_STABLE_HD_WWM.exe


################################################################################
# Algorithm preference flags.
# Comment out unwanted modules and flags.
# -DSCHISM is always on and is defined elsewhere
include ../mk/include_modules

# Don't comment out the follow ifdef
# Note: currently GOTM4 may give reasonable results only with k-omega
ifdef USE_GOTM
  #Following for GOTM4
  #GTMMOD =  -I/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/modules/PGF90/ #modules
  #GTMLIBS = -L/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/lib/PGF90/ -lturbulence_prod -lutil_prod

  # GOTM4.0
  #GTMMOD =  -I/home/bmengu01/SCHISM_WWMIII_SED3D/trunk_share_WWM/src/GOTM4.0/modules/IFORT/
  #GTMLIBS = -L/home/bmengu01/SCHISM_WWMIII_SED3D/trunk_share_WWM/src/GOTM4.0/lib/IFORT/ -lturbulence_prod -lutil_prod

  # GOTM3.2.5
  GTMMOD =  -I/home/bmengu01/SCHISM_WWMIII_SED3D/trunk_share_WWM/src/GOTM3.2.5/modules/IFORT/
  GTMLIBS = -L/home/bmengu01/SCHISM_WWMIII_SED3D/trunk_share_WWM/src/GOTM3.2.5/lib/IFORT/ -lturbulence_prod -lutil_prod
  
  #GTMMOD =  -I/sciclone/home04/yinglong/gotm-3.2.5/modules/PGF90/ #modules
  #GTMLIBS = -L/sciclone/home04/yinglong/gotm-3.2.5/lib/PGF90/ -lturbulence_prod -lutil_prod
else
  GTMMOD =
  GTMLIBS =
endif

#############################################################################################
