#$Id: Makefile,v 1.9 2005-08-16 14:38:37 kbk Exp $
#
# Master Makefile for making the GOTM executable.
#

include Rules.make

DOCSRC	= main.F90 gotm.F90

SUBDIRS	= util observations airsea meanflow turbulence output $(FEATURES)

CORE_LIBS	=	\
		-lairsea$(buildtype)		\
		-lmeanflow$(buildtype) 		\
		-lturbulence$(buildtype) 	\
		-lobservations$(buildtype)	\
		-loutput$(buildtype)		\
		-lutil$(buildtype)

ALL_LIBS	= $(FEATURE_LIBS) $(CORE_LIBS) $(EXTRA_LIBS)

EXEC		= gotm$(buildtype)_$(FORTRAN_COMPILER)

all: ../VERSION $(EXEC)

../VERSION: ../Makefile
	$(MAKE) -C ..

$(EXEC): makedirs subdirs features gotmlib 
	$(FC) -o $@ ./gotm/main.o $(LDFLAGS) -lgotm$(buildtype) $(ALL_LIBS)
	$(RM) ./gotm/main.o

makedirs:
	mkdir -p $(MODDIR) $(LIBDIR)

gotmlib:
	$(MAKE) -C gotm all

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

install: $(EXEC) 
	$(MV) $(EXEC) $(BINDIR)

$(SUBDIRS):
	$(MAKE) -C $@

observations airsea meanflow turbulence output: util

features:
ifdef FEATURES
	set -e; for i in $(FEATURES); do $(MAKE) -C $$i; done
endif

doc:
ifdef SUBDIRS
	set -e; for i in gotm $(SUBDIRS); do $(MAKE) -C $$i $@; done
endif
ifdef FEATURES
	set -e; for i in $(FEATURES); do $(MAKE) -C $$i $@; done
endif

bkup: distclean
	tar -czvf ../../backups/bkup.`date +%y%m%d%H`_`cat ../VERSION`.tar.gz *

clean:

realclean: clean
	$(RM) *~
	$(RM) ../modules/$(FORTRAN_COMPILER)/*.{m,mod} 
	$(RM) ../lib/$(FORTRAN_COMPILER)/lib*$(buildtype).a
ifdef SUBDIRS
	set -e; for i in gotm $(SUBDIRS); do $(MAKE) -C $$i $@; done
endif
ifdef FEATURES
	set -e; for i in $(FEATURES); do $(MAKE) -C $$i $@; done
endif

distclean: realclean
	$(RM) gotm_debug_* gotm_prof_* gotm_prod_*

