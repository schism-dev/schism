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
################################################################################

ENV = LNEC_WS_GNU

################################################################################
# Alternate executable name if you do not want the default. 
################################################################################
EXEC   := pschism_$(ENV)

################################################################################
# Environment 
################################################################################

FCP = /usr/lib64/openmpi/bin/mpif90
FLD = $(FCP)
# MPI vserion (1 or 2)
PPFLAGS := $(PPFLAGS) -DMPIVERSION=2
#-CB is much slower to compile
FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -ffree-line-length-none #MPI code
FLDFLAGS = -O2 #for final linking of object files
#####Libraries
##MTSLIBS = -L/home/share/binec/AAzevedo/SELFE_DEV/ParMetis-3.1-Sep2010 -lparmetis -lmetis -L/usr/lib64/openmpi/lib -lmpi -lmpi_f90
CDFLIBS = -L/usr/lib64 -lnetcdf -lnetcdff
CDFMOD = -I/usr/include -I/usr/lib64/gfortran/modules # modules for netcdf

################################################################################
# Algorithm preference flags.
# Comment out unwanted modules and flags.
################################################################################

# -DSCHISM is always on and is defined elsewhere
include ../mk/include_modules

# Don't comment out the follow ifdef
ifdef USE_GOTM
   GTMMOD =  -I/home/yinglong/GOTM/gotm-3.2.5/TSUNAMI/modules/IFORT/ #modules
   GTMLIBS = -L/home/yinglong/GOTM/gotm-3.2.5/TSUNAMI/lib/IFORT/ -lturbulence_prod  -lutil_prod
else
   GTMMOD =
   GTMLIBS =
endif


######### Specialty compiler flags and workarounds

# Add -DNO_TR_15581 like below for allocatable array problem in sflux_subs.F90
# PPFLAGS := $(PPFLAGS) -DNO_TR_15581

# For openMPI compiler, search for "USE_OPEN64" below for compiler flags
USE_OPEN64 = no

# Obsolete flags: use USE_WRAP flag to avoid problems in ParMetis lib (calling C from FORTRAN)
# PPFLAGS := $(PPFLAGS) -DUSE_WRAP 



