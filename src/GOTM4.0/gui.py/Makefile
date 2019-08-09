#$Id: Makefile,v 1.2 2007-03-02 14:11:15 jorn Exp $

# Makefile for creating the Python based GUI for GOTM

include ../src/Rules.make

# The flag -fexceptions is needed to generate FORTRAN code capable of 
# passing C++ exceptions
EXTRA_FFLAGS+=-fexceptions

NUMPYDIR = /usr/local/lib/python2.4/site-packages/numpy
NUMPYDIR = /usr/lib/python2.4/site-packages/numpy
NUMPYINC = $(NUMPYDIR)/core/include
F2PYDIR = $(NUMPYDIR)/f2py

PYTHONINC = /usr/include/python2.4
PYTHONLIB = /usr/lib

CORE_LIBS	=	\
		-lairsea$(buildtype)		\
		-lmeanflow$(buildtype) 		\
		-lturbulence$(buildtype) 	\
		-lobservations$(buildtype)	\
		-loutput$(buildtype)		\
		-lutil$(buildtype)

ALL_LIBS	= $(FEATURE_LIBS) $(CORE_LIBS) $(EXTRA_LIBS)

# Extra include directories for the C++ code (Python, NumPy and F2PY)
CPPFLAGS += -I$(PYTHONINC) -I$(F2PYDIR)/src -I$(NUMPYINC)

# Extra linker options for our Python-GOTM library. Respectively:
#   -shared: for building a shared library (*.so) rather than an executable
#   -lgotm$(buildtype): the GOTM library built from gotm.F90
#   $(ALL_LIBS): libaries built from GOTM modules (defined above)
#   -lpyhton2.4: the Python library
#   -L$(PYTHONLIB): directory that contains the Python library
LDFLAGS += -shared -lgotm$(buildtype) $(ALL_LIBS) -lpython2.4 -L$(PYTHONLIB) -lstdc++

.PHONY: clean gotm all

all: gotm gotmmodule.o fortranobject.o gotm-f2pywrappers2.o
	$(FC) ./gotm-f2pywrappers2.o ./gotmmodule.o ./fortranobject.o ./gui_util.o $(LDFLAGS) -o gotm.so

gotm:
	$(MAKE) EXTRA_FFLAGS=$(EXTRA_FFLAGS) -C ../src all

gotmmodule.o:
	$(CXX) -c gotmmodule.cpp $(CXXFLAGS) $(CPPFLAGS)

fortranobject.o:
	$(CXX) -c $(F2PYDIR)/src/fortranobject.c $(CXXFLAGS) $(CPPFLAGS)

gotm-f2pywrappers2.o: gui_util.o

clean: 
	$(MAKE) -C ../src $@
	$(RM) *~ *.o

realclean: clean 
	$(MAKE) -C ../src $@
	$(RM) *.pyc

distclean: realclean 
	$(MAKE) -C ../src $@
	$(RM) *.so

