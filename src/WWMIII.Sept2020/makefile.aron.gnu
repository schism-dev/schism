#include ${PETSC_DIR}/conf/variables
PROG = $(HOME)/bin/wwmadv
#PROG = wwmadv

FFLAGS = -DST42 -DWWM_SETUP

# Put MPI = on for compiling with mpi, nothing if no MPI is wished.
MPI = on
# Put PETSC = on for compiling with petsc, nothing if no petsc is wished.
PETSC =
# Put NETCDF = on for compiling with netcdf, nothing if no netcdf is wished.
NETCDF = on
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
OPENMP =
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
DARKO =
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
SNL4_TSA = 
# Put GRIB = on for compiling with grib library, nothing if GRIB is not wished
GRIB =  
# Put TIMINGS = on for timings with
TIME = on
# Put PDLIB = on for using the pdlib library.
PDLIB = on


#F90 = gfortran
F90 = mpif90
#F90 = mpif90-vt
#F90 = ifort
#F90 = $(FLINKER)
#F90 = pgf90
#F90 = pathf95

#INTEL
#debugging flags
#F90OPTI = -fpp -traceback -g
#F90OPTI = -O1 -axSSE4.2 -traceback -g  -check uninit -check bounds -check pointers -warn interfaces,nouncalled
#F90OPTI = -fpp -O2 -stand f03 -assume realloc_lhs -fstack-protector -assume protect_parens
#F90OPTI = -warn interfaces,nouncalled -fpp -gen-interface -g -traceback -check all

#normal optimization
#F90OPTI =  -O1 -g -traceback
#F90OPTI =  -O1 -axSSE4.2 -unroll-aggressive -assume byterecl 
#F90OPTI = -ip -O2
#F90OPTI = -warn interfaces,nouncalled -fpp -gen-interface -g -traceback -check all -fp-model precise
#F90OPTI =  -O2 -assume realloc_lhs -check all -traceback -warn all -fstack-protector -assume protect_parens -implicitnone
#F90OPTI =  -traceback -fstack-protector -assume protect_parens -g
#F90OPTI = -g -traceback -O5 -axSSE4.2 -unroll-aggressive -assume byterecl

#F90OPTI = -fpp -O2 -stand f03 -assume realloc_lhs -check all -traceback -fstack-protector -assume protect_parens -implicitnone
#extreme aggressive opti ...
#F90OPTI = -g -traceback -O1 -axSSE4.2 -unroll-aggressive -assume byterecl #-openmp#-check all -warn all

#aggressive opti
#F90OPTI = -O1 -fp-model fast -opt-multi-version-aggressive -axSSE4.2 -no-heap-arrays -fpp -stack_temps -unroll-aggressive -vec-guard-write -traceback -g
#aggressive opti no openmp
#F90OPTI = -O1 -axSSE4.2 -fpp
#debug 1
#F90OPTI = -O1 -mp1 -fpp -traceback -g -traceback -check pointers -check bound -check uninit 
#debug 2
#F90OPTI =  -fpp -g -traceback -check all
#aggresive debug openmp
#F90OPTI = -fpp -g -traceback -check all -openmp

#PGI
#F90OPTI = -g -O1 -Mprof=dwarf -Mprof=time -Minfo=ccff -mp
#F90OPTI = -fastsse -Munroll -Minline=reshape -Minline=levels:10 -Mvect -Mipa=fast Mipa=levels:10 -Mlarge_arrays -mp
#F90OPTI  = -g -Mbounds -Mchkfpstk -Mchkptr -Mchkstk

#GFORTRAN
F90OPTI = -g -O1 -ffree-line-length-0 -fbacktrace -Wall

OBJSELFE = elfe_glbl.o elfe_msgp.o solver_subs.o io_subs.o grid_subs.o misc_subs.o

OBJWWM = wwm_datapl.o wwm_petscpool.o  wwm_petsc_seriell.o wwm_petsc_parallel.o \
       wwm_petsc_block.o wwm_petsc_controller.o wwm_aux.o wwm_aux_parall.o wwm_mjdv2.o wwm_blasaux.o wwm_sparskit.o \
       wwm_ardhuin_old.o wwm_wind.o wwm_ardhuin_new.o wwm_breaking.o wwm_friction.o wwm_cycle3.o wwm_dislin.o wwm_diclass.o \
       wwm_gridcf.o wwm_hotfile.o wwm_parall_solver.o wwm_m_constants.o wwm_m_fileio.o wwm_serv_xnl4v5.o wwm_mod_xnl4v5.o wwm_snl4_wrt.o \
       wwm_wave_setup.o wwm_initio.o wwm_netcdf.o wwm_input.o wwm_grid.o \
       wwm_bdcons_init.o wwm_bdcons.o wwm_bdcons_wam.o wwm_nesting.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_snl42.o wwm_snl4_tsa.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_snonlin.o wwm_snonlin_local.o wwm_stress.o wwm_stresso.o wwm_stresso_local.o wwm_sbottom.o wwm_sdiss_ardh_vec.o wwm_sdiss_ardh_vec_local.o wwm_sinput.o wwm_sinput_local.o wwm_sinput_ard.o wwm_sinput_ard_local.o wwm_wsigstar.o wwm_wsigstar_local.o wwm_tauhf.o wwm_airsea.o wwm_airsea_local.o wwm_jafu.o wwm_nlweigt.o wwm_buildstress.o wwm_inisnonlin.o wwm_fkmean.o wwm_fkmean_local.o wwm_frcutindex.o wwm_frcutindex_local.o wwm_femeanws.o wwm_femeanws_local.o wwm_sdissip.o wwm_sdissip_local.o wwm_implsch.o wwm_implsch_local.o wwm_implsch2.o \
       wwm_output.o wwm_compute.o wwm_diffrac.o wwm_coupl_roms_pipe.o wwm_coupl_roms_pgmcl.o wwm_jacobi.o wwm_coupl_shyfem.o wwm_coupl_selfe.o wwm_coupl_timor.o wwm_main.o wwm_export_ww3.o
#       dislin.o diclass.o \

FFLAGS += -DWWM_SOLVER
ifdef NETCDF
  all: $(NETC) $(PROG)
  LIBS += $(NETCDF_FORTRAN_LINK)
  FFLAGS +=  -DNCDF -I$(NETCDF_INCDIR)
else
  all: $(PROG)
endif

ifdef TIME
  FFLAGS +=  -DTIMINGS
endif

ifdef PDLIB
  LIBS += $(PDLIB_PATH)/libpd.a
  FFLAGS += -I$(PDLIB_PATH)/modules/ -I$(PDLIB_PATH)/include -DPDLIB
endif

ifdef MPI
  LIBS += $(METIS_PATH)/lib/libparmetis.a $(METIS_PATH)/lib/libmetis.a
  MPIFLAG  = -DWWM_MPI -DMPIVERSION=2 -DUSE_WWM #-DSHYFEM_COUPLING
  ifdef PDLIB
    OBJS= wwm_pdlib.o $(OBJWWM)
  else
    OBJS= $(OBJSELFE) $(OBJWWM)
  endif
else
  MPIFLAG = 
  OBJS= $(OBJWWM)
endif

ifdef OPENMP 
  OPENMPFLAG  = -openmp 
endif

ifdef GRIB
  LIBS += $(GRIB_FORTRAN_LINK)
  FFLAGS += -DGRIB_API_ECMWF -I$(GRIB_PATH)/include
endif


ifdef PETSC
  MPIFLAG += -DPETSC ${PETSC_FC_INCLUDES} -DDIRECT_METHOD
  LIBS  += ${PETSC_LIB}  
endif

ifdef DARKO 
  DARKOFLAG = -DDARKO 
endif

ifdef SNL4_TSA 
  SNL4_TSAFLAG = -DSNL4_TSA
endif

F90FLAGS = ${MPIFLAG} ${FFLAGS} ${F90OPTI} ${OPENMPFLAG} ${DARKOFLAG} ${SNL4_TSAFLAG}

NETC = NETC

$(NETC):
	cp -f $(NETCDF_INCDIR)/netcdf.mod .
	cp -f $(NETCDF_INCDIR)/typesizes.mod .

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(F90FLAGS)  $(OBJS) $(LIBS) 

clean:
	rm -f *.o *.oo *.obj *.ipo *.mod *.map *__genmod.f90 *.ilk *.pdb $(PROG)

cleanall:
	rm -f *.*~ *.spg *.smb *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb *genmod* $(PROG)

.SUFFIXES: $(SUFFIXES) .F90 .f .ftn

.F90.o:
	$(F90) $(F90FLAGS)   -c $<

.ftn.o:
	$(F90) $(F90FLAGS) -c $<

