################################################################################
# Parallel SCHISM Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below. These settings are based on the following options.
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
# Libraries (needed for parallel code). Need to re-compile these when compilers are changed.
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   netcdf: v3.x has been tested
#   GOTM77051428:    v3.2.5
#   ALTLIBS: (obsolete) Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
################################################################################

##########

########### Ensure that the user's Make.defs.local exists and include it

make_defs_local_exists := $(wildcard ../mk/Make.defs.local)
ifeq ($(strip $(make_defs_local_exists)),)
   $(error Make.defs.local does not exist in the /mk directory (parallel to this /src directory). See readme)
endif


-include ../mk/Make.defs.local  #optional include
include ../mk/Make.defs.config

ifeq ($(findstring GNU Fortran,$(shell $(FCP) --version)),GNU Fortran)
  $(print Here)
  MODULEOUTDIR = -J
else
  ifeq ($(ENV),DKRZold) ##please edit here, Eli
    MODULEOUTDIR = -qmoddir=
  else
    MODULEOUTDIR = -module 
  endif
endif


##########
# Check for serial compiler and revert to parallel if not provided
#########
ifeq ($(FCS),)
    FCS := ${FCP}
endif

##PetSc config file that defines things like PETSC_FC_INCLUDES below
ifeq ($(PETSC_VERSION),v1) #v3.5
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DPETSCV=1
  include ${PETSC_DIR}/conf/variables
else
  ifeq ($(PETSC_VERSION),v2) ##v3.6, 3.7
    PPFLAGS := $(PPFLAGS) $(IBMCPP)-DPETSCV=2
    include ${PETSC_DIR}/lib/petsc/conf/variables
  endif
  ifeq ($(PETSC_VERSION),v3) #v3.10
    PPFLAGS := $(PPFLAGS) $(IBMCPP)-DPETSCV=3 -I${PETSC_DIR}/include/petsc
    include ${PETSC_DIR}/lib/petsc/conf/variables
  endif
endif

############
# Derived compiler flags
###########

ifeq ($(ENV),DKRZold)
  IBMCPP = -WF,
else
  IBMCPP =
endif

PPFLAGS := $(PPFLAGS) $(IBMCPP)-DSCHISM

ifneq ($(USE_OPEN64),)
    ifeq ($(USE_64), yes)
      AVOID_ADV_WRITE := $(USE_OPEN64)
      PPFLAGS := $(PPFLAGS) $(IBMCPP)-DAVOID_ADV_WRITE $(IBMCPP)-DUSE_OPEN64
    endif
endif

################################################################################
# Object files; modules compile first; modules that depend on other modules compile last.
################################################################################

# rules to build objects for executable programs and libraries

#Search path
VPATH = Core/ Core/PaHM/ Driver/ Hydro/ EcoSim/ COSINE/ Sediment/ Sed2d/ \
        WWMIII/ ICM/ TIMOR/ FIB/ Ice/ Multi_ice/ Fabm/

ifdef FUJITSU
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DFUJITSU
endif

include o/$(config)/autosrc.mk

#Cores
#OBJS =  $(AUTO_OBJ) PaHM_Sizes.o PaHM_Global.o PaHM_Messages.o TimeDateUtils.o PaHM_Utilities.o \
#        csv_module.o SortUtils.o csv_parameters.o PaHM_Vortex.o \
#        csv_utilities.o  parwind.o 
OBJS :=  $(AUTO_OBJ) schism_glbl.o schism_msgp.o scribe_io.o misc_modules.o schism_io.o schism_driver.o

ifdef USE_PETSC
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_PETSC
  PETSCLIBS = ${PETSC_LIB}
  PETSCMOD = ${PETSC_FC_INCLUDES}

  OBJS := $(OBJS) petsc_schism.o
else
  PETSCLIBS =
  PETSCMOD =
endif

ifdef NO_PARMETIS
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DNO_PARMETIS
endif

ifdef SH_MEM_COMM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DSH_MEM_COMM
endif

ifdef USE_OLDIO
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DOLDIO
endif

ifdef USE_ATMOS
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DATMOS
endif

ifdef USE_BMI
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_BMI
endif

ifdef USE_PREC_EVAP
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DPREC_EVAP
endif

ifdef USE_BULK_FAIRALL
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_BULK_FAIRALL
endif

ifdef USE_MM5
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DMM5
endif

ifeq ($(TVD_LIM),SB)
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DTVD_SB
else
  ifeq ($(TVD_LIM),VL)
    PPFLAGS := $(PPFLAGS) $(IBMCPP)-DTVD_VL
  else
    ifeq ($(TVD_LIM),OS)
      PPFLAGS := $(PPFLAGS) $(IBMCPP)-DTVD_OS
    else 
      ifeq ($(TVD_LIM),MM)
        PPFLAGS := $(PPFLAGS) $(IBMCPP)-DTVD_MM
      else
        $(error TVD_LIM must be defined!)
      endif
    endif
  endif
endif

ifdef USE_DEBUG
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DDEBUG 
endif

ifdef USE_TIMER2
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DTIMER2
endif

ifdef USE_TIMER
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DINCLUDE_TIMING
endif

ifdef USE_TIMOR_FLMUD
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_TIMOR_FLMUD
  OBJS := $(OBJS) flmud_input.o flmud_main.o flmud_pool.o
endif

ifdef USE_TIMOR
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_TIMOR
  OBJS := $(OBJS) timor_datapl.o\
       timor_3D_Boden.o timor_area.o timor_Bank.Avala.o  timor_bed.o \
       timor_Bilanzen.o timor_cr_w.o timor_dredge.o timor_fluc.o \
       timor_GravTrans.o timor_hydro_Anlagen.o\
       timor_init.o timor_Q3D.o timor_ReadWrite.o timor_schub.o\
       timor_sukzession.o timor_SusTrans.o timor_timeloop.o \
       timor_upw_ups.o timor_waves.o
endif

ifdef USE_GOTM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_GOTM
endif

ifdef USE_HA
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_HA
  OBJS := $(OBJS) harm.o
endif

ifdef USE_PAHM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_PAHM
  OBJS := $(OBJS) PaHM_Sizes.o PaHM_Global.o PaHM_Messages.o TimeDateUtils.o PaHM_Utilities.o \
        csv_module.o SortUtils.o csv_parameters.o PaHM_Vortex.o \
        csv_utilities.o  parwind.o
endif

ifdef USE_GEN
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_GEN
endif

ifdef USE_AGE
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_AGE
endif

ifdef USE_ECO
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ECO
  OBJS := $(OBJS) bio_param.o eclight.o biology.o bio_init.o ecosim.o read_ecoin.o specir.o
endif

ifdef USE_ICM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ICM
  OBJS := $(OBJS) icm_mod.o icm_misc.o icm_init.o icm_sfm.o icm.o
endif

ifdef USE_COSINE
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_COSINE
  OBJS := $(OBJS) cosine_misc.o cosine_mod.o cosine_init.o cosine.o 
endif

ifdef USE_FIB
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_FIB
###Add FIB file here
  OBJS := $(OBJS) fib_param.o fib.o
endif

ifdef USE_NAPZD
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_NAPZD
  OBJS := $(OBJS) biology_napzd.o napzd_spitz.o read_napzd_input.o
endif

ifdef USE_SED
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_SED
  OBJS := $(OBJS) sed_mod.o sed_misc_subs.o \
                  read_sed_input.o sed_init.o  \
                  sed_friction.o sed_filter.o \
                  sediment.o sed_bedload.o

endif

ifdef USE_FABM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_FABM
  FABMINSTALLDIR = ./Fabm/build
  ifneq ($(wildcard $(FABMINSTALLDIR)/lib64),)
    FABMLIBS = -L$(FABMINSTALLDIR)/lib64 -lfabm
  else
    FABMLIBS = -L$(FABMINSTALLDIR)/lib -lfabm
  endif
  OBJS := $(OBJS) fabm_schism.o
  FCPFLAGS += -I$(FABMINSTALLDIR)/include
endif

ifdef USE_DVD
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_DVD
endif

ifdef USE_OIL
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_OIL
endif

ifdef USE_WWM
# please take care we have now ST41 (ardhuin_old) and ST42 (ardhuin_new)
  PPFLAGS := $(PPFLAGS) -DUSE_WWM -DNCDF -DST42 -DWWM_SOLVER 
  OBJS := $(OBJS) wwm_datapl.o wwm_petscpool.o  wwm_petsc_seriell.o wwm_petsc_parallel.o \
       wwm_petsc_block.o wwm_petsc_controller.o wwm_aux.o wwm_aux_parall.o wwm_mjdv2.o wwm_blasaux.o wwm_sparskit.o \
       wwm_ardhuin_old.o wwm_wind.o wwm_ardhuin_new.o wwm_breaking.o wwm_friction.o wwm_cycle3.o  \
       wwm_gridcf.o wwm_hotfile.o wwm_parall_solver.o wwm_m_constants.o wwm_m_fileio.o wwm_serv_xnl4v5.o wwm_mod_xnl4v5.o wwm_snl4_wrt.o \
       wwm_wave_setup.o wwm_initio.o wwm_netcdf.o wwm_input.o wwm_grid.o \
       wwm_bdcons_init.o wwm_bdcons.o wwm_bdcons_wam.o wwm_nesting.o wwm_export_ww3.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_snl42.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_snonlin.o wwm_stress.o wwm_stresso.o wwm_sbottom.o wwm_sdiss_ardh_vec.o wwm_sinput.o wwm_sinput_ard.o wwm_wsigstar.o wwm_tauhf.o wwm_airsea.o wwm_jafu.o wwm_nlweigt.o wwm_buildstress.o wwm_inisnonlin.o wwm_fkmean.o wwm_frcutindex.o wwm_femeanws.o wwm_sdissip.o wwm_implsch.o wwm_implsch2.o \
       wwm_output.o wwm_compute.o wwm_compute_roller.o wwm_diffrac.o wwm_coupl_roms_pipe.o wwm_coupl_roms_pgmcl.o wwm_jacobi.o wwm_coupl_shyfem.o wwm_coupl_selfe.o wwm_coupl_timor.o wwm_vegdiss.o wwm_main.o wwm_lpp_filter.o
endif

ifdef USE_WW3
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_WW3
endif

ifdef USE_ICE
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ICE
  OBJS := $(OBJS) ice_module.o ice_therm_mod.o ice_init.o ice_step.o ice_evp.o ice_mevp.o \
                  ice_fct.o ice_thermodynamics.o
endif

ifdef USE_MICE
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_MICE
  OBJS := $(OBJS) icepack_kinds.o icepack_warnings.o icepack_parameters.o icepack_age.o icepack_tracers.o \
                  icepack_mushy_physics.o icepack_orbital.o icepack_ocean.o icepack_firstyear.o \
                  icepack_therm_shared.o icepack_therm_bl99.o icepack_meltpond_lvl.o \
                  icepack_meltpond_topo.o icepack_zbgc_shared.o icepack_fsd.o icepack_atmo.o icepack_therm_mushy.o \
                  icepack_isotope.o icepack_flux.o icepack_itd.o icepack_shortwave.o icepack_brine.o \
                  icepack_aerosol.o icepack_wavefracspec.o \
                  icepack_mechred.o icepack_algae.o icepack_therm_vertical.o icepack_snow.o icepack_zbgc.o icepack_therm_itd.o \
                  icepack_intfc.o icedrv_kinds.o icedrv_constants.o mice_module.o icedrv_system.o \
                  icedrv_main.o icedrv_allocate.o icedrv_set.o mice_init.o mice_evp.o icedrv_init.o\
                  mice_mevp.o icedrv_advection.o icedrv_transfer.o icedrv_io.o mice_step.o gen_modules_clock.o
endif

ifdef USE_MARSH
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_MARSH
endif

ifdef USE_ANALYSIS
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ANALYSIS
endif

ifdef USE_SED2D
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_SED2D
  OBJS := $(OBJS) sed2d_mod.o sed2d_read.o sed2d_filter.o sed2d_transport.o \
                  sed2d_friction.o sed2d_morpho.o sed2d_check.o sed2d_init.o \
                  sed2d_main.o
endif

#Remaining Hydro
OBJS := $(OBJS) grid_subs.o hydraulic_structures.o schism_init.o \
        schism_step.o schism_finalize.o bktrk_subs.o solver_subs.o \
        misc_subs.o transport_TVD_imp.o \
        sflux_9c.o lap.o

ifdef USE_SIMPLE_WIND
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_SIMPLE_WIND
   OBJS := $(OBJS) simpler_wind.o
endif

# Names of the object files with configuration intermediate file (.o,.mod) output directory prepended
CONFIG_OBJS :=  $(addprefix o/$(config)/, $(OBJS))

# Names of source files with relative path as found on $VPATH
FULLSRC := $(foreach f,$(subst .o,.F90,$(OBJS)),$(firstword $(wildcard $(VPATH:%=%$f) not-found:$f)))

PARMETISDIR := ./ParMetis-4.0.3
ifeq ($(ENV),EDISON)
    PARMETISDIR := $(PARMETISDIR).cray
endif
ifeq ($(ENV),CORI2)
    PARMETISDIR := $(PARMETISDIR).cray
endif
ifeq ($(ENV),BEBOP)
    PARMETISDIR := $(PARMETISDIR).mpiicc
endif
ifeq ($(ENV),ANVIL)
    PARMETISDIR := $(PARMETISDIR).mpiicc
endif
ifeq ($(ENV),COMET)
    PARMETISDIR := $(PARMETISDIR).mpiicc
endif
ifeq ($(ENV),STAMPEDE2)
    PARMETISDIR := $(PARMETISDIR).mpiicc
endif

################################################################################
# Library list
################################################################################
ifdef NO_PARMETIS
#Use offline partitioning and disable ParMETIS
  MTSLIBS= 
else
  MTSLIBS= -L$(PARMETISDIR) -lparmetis -lmetis
endif
LIBS = $(MTSLIBS) $(FABMLIBS) $(CDFLIBS) $(GTMLIBS) $(PETSCLIBS)

################################################################################
# Build commands for SCHISM
################################################################################

.PRECIOUS: $(SELFE_HOME)/mk/.check.$(config) o/$(config)/%.o o/$(config)/.dummy


# Allow for the override of executable name by the user
ifeq ($(EXEC),)
    EXEC := pschism.$(config).ex
endif




o/$(config)/autosrc.mk: o/$(config)/.dummy
	python Core/gen_version.py
	echo "AUTO_OBJ=`echo schism_version.o`"  > o/$(config)/autosrc.mk



.PHONY: pschism psc

pschism: $(EXEC) o/$(config)/.dummy 

psc: $(EXEC) o/$(config)/.dummy 

.PHONY: all
all: pschism utility


# This is the main executable target
ifdef NO_PARMETIS
$(EXEC): $(CONFIG_OBJS) 
	$(FLD) $(FLDFLAGS) $(CONFIG_OBJS) $(LIBS) ${MODULEOUTDIR}. -I. -o $@
	ar -ru libschism.a $(CONFIG_OBJS)
else
$(EXEC): $(CONFIG_OBJS) $(PARMETISDIR)/parmetis.a 
	$(FLD) $(FLDFLAGS) $(CONFIG_OBJS) $(LIBS) ${MODULEOUTDIR}. -I. -o $@
	ar -ru libschism.a $(CONFIG_OBJS)


$(PARMETISDIR)/parmetis.a:
##	+$(MAKE) distclean --directory=$(PARMETISDIR)
	+$(MAKE) config --directory=$(PARMETISDIR)
	+$(MAKE) --directory=$(PARMETISDIR)
endif

#.PHONY: parmetis
#parmetis:
#	+$(MAKE) clean --directory=$(PARMETISDIR)
#	+$(MAKE) --directory=$(PARMETISDIR)


.PHONY: depend
depend: o/$(config)/.dummy o/$(config)/.depend ../bin/.dummy
	@echo Dependencies generated



o/$(config)/.depend: o/$(config)/.dummy 
	@echo Building depend list
	@echo $(util_srcs_abs)
	@ ../mk/sfmakedepend.pl -d o/$(config) -I $(VPATH) -f o/$(config)/.depend_extended $(FULLSRC) > depend.out
	@echo $(FULL_CONFIG_OBJS)
	@ ../mk/cull_depends.py o/$(config)/.depend_extended o/$(config)/.depend "$(FULL_CONFIG_OBJS)"


o/$(config)/schism_version.mod: o/$(config)/schism_version.o

o/$(config)/schism_version.o: Core/schism_version.F90 

Core/schism_version.F90:

ifneq "$(MAKECMDGOALS)" "clean"
  -include o/$(config)/.depend
endif

o/$(config)/%.o: %.F90 o/$(config)/.dummy o/$(config)/.depend
	$(QUIET)$(FCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTMMOD) $(CDFMOD) $(PETSCMOD) $< -o $@

o/$(config)/fabm_schism.o: $(FABMINSTALLDIR)/lib/libfabm.a

$(FABMINSTALLDIR)/lib/libfabm.a:
	mkdir -p $(FABMINSTALLDIR)
	cd $(FABMINSTALLDIR) && cmake $(FABM_BASE)/src -DCMAKE_Fortran_COMPILER=$(FCP) -DFABM_HOST=schism -DCMAKE_INSTALL_PREFIX=. && make install

################# Utility
# The utility directory contains a lot of single-file executables plus a few that rely on shared subroutines
# for things like geometry calculations.
#
# The executables will be discovered automatically
# Files with re-usable routines should be maually added to the list util_lib_src a few lines below.
# Please consider moving the common files to UtilLib if you think the subroutines can be re-used 
# by utility scripts in several directories ... we'd like to start creating a modern API-driven library.


src_extension := .f90
util_lib_src := schism_geometry.f90 argparse.f90  extract_mod2.f90  interpolate_unstructured_module.f90  pt_in_poly_test.f90  stripesearch_unstr.f90 compute_zcor.f90 extract_mod.f90 julian_date.f90                     

util_lib_obj := $(subst $(src_extension),.o,$(util_lib_src))
util_lib_config_objs := $(addprefix o/$(config)/, $(util_lib_obj))


util_dirs :=  Combining_Scripts Grid_Scripts NCOM OneWayNestScripts Particle_Tracking Post-Processing-Fortran Sflux_nc SMS Stat3 UtilLib ArcGIS
util_dirs := $(addprefix Utility/,$(util_dirs))
VPATH := $(VPATH) $(util_dirs)


# list of source files in all the other source directories
ifneq ($(util_dirs),)
  util_srcs_abs := $(wildcard $(foreach ext,$(src_extension),$(addsuffix /*$(ext),$(util_dirs))))
  util_srcs_all := $(notdir $(util_srcs_abs) )
endif

# remove library items that don't compile as executables
util_srcs  := $(filter-out $(util_lib_src),$(util_srcs_all))
util_execs := $(addprefix ../bin/,$(subst $(src_extension),.ex,$(util_srcs)))
util_objs  := $(subst $(src_extension),.o,$(util_srcs))
util_config_objs := $(addprefix o/$(config)/, $(util_objs)) 

FULLSRC := $(FULLSRC) $(util_srcs_abs)
FULL_CONFIG_OBJS = $(CONFIG_OBJS) $(util_config_objs) $(util_lib_config_objs)

.PHONY: utility

utility: $(util_execs) o/$(config)/libutil_lib.a


#$(util_execs): $(util_config_objs) o/$(config)/libutil_lib.a
#	$(FCS) $(SCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Lo/$(config) -lutil_lib -Io/$(config)  $(CDFLIBS) $(GTTMOD) $(CDFMOD) $< -o $@
$(util_execs): ../bin/%.ex : o/$(config)/%.o o/$(config)/libutil_lib.a
	$(FCS) $(SCPFLAGS) $< ${MODULEOUTDIR}o/$(config) -Lo/$(config) -lutil_lib -Io/$(config)  $(CDFLIBS) $(GTTMOD) $(CDFMOD) -o $@

o/$(config)/libutil_lib.a: $(util_lib_config_objs)
	$(AR) $(ARFLAGS) $@ $^


#$(util_config_obs) : $(util_srcs) o/$(config)/.dummy  ../bin/.dummy
#	$(SCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@


o/$(config)/%.o : %.f90 o/$(config)/.dummy  ../bin/.dummy
#	$(FCS) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@
	$(FCS) -c -cpp $(FCPFLAGS) $< ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) -o $@


../bin/%.ex: o/$(config)/%.o  ../bin/.dummy
	$(FLD) $(FLDFLAGS) $< $(LIBS) ${MODULEOUTDIR}. -I. -o $@





#################################################################################

# This file exists solely to force creation of the subdirectories that
# hold the dependency, object and preprocessor files.
# NOTE: we cant use a dependency for this because .dummy will always
#       be out-of-date with respect to the subdirs because they are
#       modified after .dummy is created.  But since there are no
#       dependencies, this rule only runs when the file doesn't exist.
outdir_exists := $(wildcard o/$config)
make_defs_local_exists := $(wildcard ../mk/Make.defs.local)

o/$(config)/.dummy:
	@echo $(make_defs_local_exists)
        ifeq ( $(strip $(make_defs_local_exists)), )
             $(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template.  \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])
        endif
	$(ECHO)mkdir -p o/$(config)
	@touch o/$(config)/.dummy


../bin/.dummy:
	$(ECHO)mkdir -p ../bin
	@touch ../bin/.dummy

../mk/Make.defs.local :
	$(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template. \n \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])



################################################################################
# Clean up or Clobber
################################################################################
# Clean up objects and modules. Clean kills all configurations.

.PHONY: clean
clean:
	$(QUIET)rm -f $(EXEC)
	$(QUIET)rm -f libschism.a
	$(QUIET)rm -rf o/*
	$(QUIET)rm -rf ../bin/*
	$(QUIET)rm -rf WWMIII/*.mod
	$(QUIET)rm -rf WWMIII/*.o
	(cd $(PARMETISDIR); make distclean; rm -f *.a )
	##(cd $(PARMETISDIR)/ParMETISLib ; make realclean )
	$(QUIET)rm -rf $(FABMINSTALLDIR)

