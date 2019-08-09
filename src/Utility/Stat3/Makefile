PROG = $(HOME)/bin/stat3

OBJS = stat_pool.f90 stat3.f90 \
#       dislin.o diclass.o \

#LIBS = /opt/dislin/lib/dislin-9.5.so
#LIBS = /opt/dislin/libdislin.so

F90 = ifort
#F90 = $(FLINKER)
#F90 = gfortran
#F90 = pgf90
#F90 = pathf95
#F90 = f77

#INTEL
#save opti
F90FLAGS = -O1 -axSSE4.2 -fpp -DWWMONLY -DNCDF -DST4 -assume byterecl -traceback -g 
#F90FLAGS = -O1 -mp1 -fpp -DWWMONLY -traceback -g -axSSE4.2 -g -traceback -check all -assume byterecl

#LDFLAGS = -Mprof=dwarf -Mprof=lines -Mprof=func -Mprof=time Minfo=ccff
#TCFLAG  = -tcheck
#LDFLAGS = -lpthread

#ITT_BASE = /home/aron/intel/itt

#LIBS = #libgoto_penrynp-r1.26.a #sparskit_pgi.a libgoto.a #sparskit.a

#ARCH  = $(shell uname -m | sed -e s/x86_64/32e/ -e s/i.86/32/ -e s/ia// )
#LPATH = $(ITT_BASE)/tcheck/bin/$(ARCH)

F90FLAGS += ${PETSC_INCLUDE}

NETC = NETC

all: $(NETC) $(PROG)

$(NETC):
#	cp -f $(NETCDF_INCDIR)/netcdf.mod .
#	cp -f $(NETCDF_INCDIR)/typesizes.mod .

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(F90FLAGS) $(TCFLAG) $(OBJS) $(LIBS) $(LDFLAGS) ${PETSC_LIB}

clean:
	rm *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb $(PROG)

cleanall:
	rm  *.*~ *.spg *.smb *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb $(PROG)

.SUFFIXES: $(SUFFIXES) .F90 .f .ftn

.F90.o:
	$(F90) $(F90FLAGS) $(TCFLAG)   -c $<

.ftn.o:
	$(F90) $(F90FLAGS) -c $<
