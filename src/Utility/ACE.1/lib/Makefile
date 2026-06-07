# Makefile for graphics support libraries
#

INSTALLDIR = /usr/local/ace

all:
	cd graphics; make
	cd fli; make

install:
	cp -p graphics/libacegraph.a ${INSTALLDIR}/lib
	cp -p fli/libfpp.a ${INSTALLDIR}/lib
	cp -p graphics/*.h ${INSTALLDIR}/include
	cp -p fli/*.h ${INSTALLDIR}/include

clean:
	cd graphics; make clean
	cd fli; make clean
