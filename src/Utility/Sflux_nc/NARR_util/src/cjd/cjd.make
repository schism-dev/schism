###################################################################
############### Overall definition of build target ################
###################################################################

PROGRAM = cjd

###################################################################
############### Options for this particular build #################
###################################################################

# do we need to enable F2KCLI (for command line options),
# and if so, what version
NEED_F2KCLI    = yes
F2KCLI_VERSION = unix_f90

###################################################################
############### Options for this specific arcitecture #############
###################################################################

FC            = ifc
LL            = $(FC)
FPP_CMD       = -D

FFLAGS        = -I/usr/local/netcdf/include -Vaxlib
LFLAGS        = -Vaxlib -static

LIBS          = -L/usr/local/netcdf/lib -lnetcdf

OBJS          = $(PROGRAM).o

###################################################################
#### Set the preprocessing flags based on configurations above ####
###################################################################

#### F2KCLI for command line options? ####
ifdef NEED_F2KCLI
  FPP_FLAGS := $(FPP_CMD)NEED_F2KCLI $(FPP_CMD)$(F2KCLI_VERSION) \
               $(FPP_FLAGS)
  OBJS := f2kcli.o $(OBJS)
endif

###################################################################
######################### Build commands ##########################
###################################################################

$(PROGRAM): $(OBJS)
	$(LL) -o $(PROGRAM) $(OBJS) $(LIBS) $(LFLAGS)

%.o: %.F 
	$(FC) $(FPP_FLAGS) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(PROGRAM)
