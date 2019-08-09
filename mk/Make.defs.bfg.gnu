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

################################################################################
## Environment for the oss-hpc01 cluster of the BfG
#################################################################################
ENV = hpc01
VER = trunk

################################################################################
# Environment for BfG oss-hpc1 (HPC Xeon Cluster, Xeon E5-2670, Qlogic Infiniband (QDR))
# with GNU compilers
# requirement: gfortran, gcc, libnetcdf-dev, libnetcdff, mvapich2
#
#   before make check if following modules are loaded
#     module load parallel/mvapich2-2.1rc2
#     module load libraries/netcdf-fortran
#
# Note: use mpiexec to run compiled SCHISM executable   
################################################################################

# GNU compiler
FCP = mpif90 -ffree-line-length-none
FCS = gfortran
FLD = $(FCP)
# MPI vserion (1 or 2) 
PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 #-DUSE_WRAP
OPTFLAGS = -O2
FCPFLAGS = $(PPFLAGS) $(OPTFLAGS) -static -static-libgfortran -finit-local-zero
FLDFLAGS = $(OPTFLAGS) #for final linking of object files

#debugging
#OPTFLAGS = -O
#FCPFLAGS = $(PPFLAGS) $(OPTFLAGS) -static -static-libgfortran -Wuninitialized  
#FLDFLAGS = $(OPTFLAGS) #for final linking of object files

##### Libraries
# ParMETIS
MTSLIBS = -L./ParMetis-3.1-Sep2010  -lparmetis -lmetis  
CDFLIBS = -L/usr/lib64 -lnetcdf  -L/opt/netcdf-fortran/lib -lnetcdff #-L/opt/mvapich2-2.1rc2/lib -lfmpich -lmpich -lmpi # oss-hpc01
CDFMOD = -I/usr/include -I/opt/netcdf-fortran/include # modules for netcdf 



################################################################################
# Alternate executable name if you do not want the default. 
################################################################################
EXEC   := pschism_$(ENV)_$(VER)

################################################################################
# Algorithm preference flags.
# Comment out unwanted modules and flags.
################################################################################

# -DSCHISM is always on and is defined elsewhere
include ../mk/include_modules

# Don't comment out the follow ifdef
# Note: currently GOTM4 may give reasonable results only with k-omega
ifdef USE_GOTM
  #Following for GOTM4
  #GTMMOD =  -I/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/modules/PGF90/ #modules
  #GTMLIBS = -L/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/lib/PGF90/ -lturbulence_prod -lutil_prod

  #Following for GOTM3
  GTMMOD =  -I/sciclone/home04/yinglong/gotm-3.2.5/modules/PGF90/ #modules
  GTMLIBS = -L/sciclone/home04/yinglong/gotm-3.2.5/lib/PGF90/ -lturbulence_prod -lutil_prod
else
  GTMMOD =
  GTMLIBS =
endif


######### Specialty compiler flags and workarounds
# Add -DNO_TR_15581 like below for allocatable array problem in sflux_subs.F90
# PPFLAGS := $(PPFLAGS) -DNO_TR_15581

# Obsolete flags: use USE_WRAP flag to avoid problems in ParMetis lib (calling C from FORTRAN)
# PPFLAGS := $(PPFLAGS) -DUSE_WRAP 

