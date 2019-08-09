# $Id: Makefile,v 1.7 2011/09/14 17:44:21 pturner Exp $
#
# makefile for xmgredit version 4.0
#

INSTALLDIR = /usr/local/ace

# NetCDF support
#################################
# Uncomment the following 3 defines for netCDF support
# adjust paths to suit local conditions
#
#NETCDF = -DHAVE_NETCDF
#NETCDFINCLUDES = -I/usr/local/ace/netcdf/include
#NETCDFLIBS = -L/usr/local/ace/netcdf/lib -lnetcdf

############################################
# LINUX
LIBS = -L/usr/X11R6/lib64 -L/opt/local/lib -lXm -lXp -lXt -lXext -lXpm -lX11 -lICE -lSM -lm
INCLUDES = -I/usr/X11R6/include  -I/opt/local/include
CC = gcc -g

#################
####### End of configuration, no changes should be required below #########
#################
CFLAGS = -DDO_TRIANGLE $(INCLUDES)

OBJS = main.o\
	vers.o\
	xmgredit.o\
	gridwin.o\
	fileswin.o\
	statuswin.o\
	propwin.o\
	editwin.o\
	bcwin.o\
	modwin.o\
	qualwin.o\
	timerwin.o\
	slicewin.o\
	graphwin.o\
	isolwin.o\
	buildwin.o\
	printwin.o\
	locatewin.o\
	wvwin.o\
	imagewin.o\
	monwin.o\
	helpwin.o\
	scalewin.o\
	cedit.o\
	events.o\
	strwin.o\
	objutils.o\
	find.o\
	allogrid.o\
	allobound.o\
	allobuild.o\
	gridio.o\
	gridutils.o\
	ibilinear.o\
	buildio.o\
	buildutils.o\
	boundio.o\
	boundutils.o\
	drawobjs.o\
	drawgrid.o\
	drawbound.o\
	isol.o\
	isolutils.o\
	stubs.o\
	stations.o\
	utils.o\
	io.o\
	params.o\
	malerts.o\
	motifutils.o\
	triangle.o\
	tritest.o\
	vedglist.o\
	vgeometry.o\
	vheap.o\
	vmemory.o\
	voronoi.o\
	voutput.o\
	vtriang.o\
	chersh.o\
	draw.o\
	image_f_io.o\
	ps.o\
	xvlib.o

SRCS = main.c\
	vers.c\
	xmgredit.c\
	gridwin.c\
	fileswin.c\
	statuswin.c\
	propwin.c\
	editwin.c\
	bcwin.c\
	modwin.c\
	qualwin.c\
	timerwin.c\
	slicewin.c\
	graphwin.c\
	isolwin.c\
	buildwin.c\
	printwin.c\
	locatewin.c\
	wvwin.c\
	imagewin.c\
	monwin.c\
	helpwin.c\
	scalewin.c\
	cedit.c\
	events.c\
	strwin.c\
	objutils.c\
	find.c\
	allogrid.c\
	allobound.c\
	allobuild.c\
	gridio.c\
	gridutils.c\
	ibilinear.c\
	buildio.c\
	buildutils.c\
	boundio.c\
	boundutils.c\
	drawobjs.c\
	drawgrid.c\
	drawbound.c\
	isol.c\
	isolutils.c\
	stubs.c\
	stations.c\
	utils.c\
	io.c\
	params.c\
	malerts.c\
	motifutils.c\
	triangle.c\
	tritest.c\
	vedglist.c\
	vgeometry.c\
	vheap.c\
	vmemory.c\
	voronoi.c\
	voutput.c\
	vtriang.c\
	chersh.c\
	draw.c\
	image_f_io.c\
	ps.c\
	xvlib.c

#
# Command parser
#
PARSOBJS = pars.o
PARSSRCS = pars.y

INCS = motifinc.h\
	f2c.h\
	defines.h\
	globals.h\
	gredit_icon.h\
	bitmaps.h\
	symdefs.h\
	externs.h\
	hersh.h\
	patterns.h\
	special.h\
	graphics.h\
	triangle.h\
	vdefines.h\
	vdefs.h

all: xmgredit5

xmgredit5: $(OBJS) $(PARSOBJS)
	$(CC) $(OBJS) $(PARSOBJS) -o xmgredit5 $(LIBS)

tx: tx.o
	$(CC) tx.o -o tx -L/usr/local/lib -lgd $(LIBS)

$(OBJS): defines.h globals.h
eventproc.c: defines.h globals.h

triangle.o:
	$(CC) -c -DTRILIBRARY triangle.c

vers.o: $(SRCS) pars.y
	sh newvers.sh
	$(CC) -c vers.c

pars.o: pars.c
pars.c: pars.y
	bison -y pars.y
	mv y.tab.c pars.c

print: Makefile $(SRCS) $(INCS) $(PARSSRCS)
	a2ps $? | lpr -h
	touch print

distrib: 
	-/bin/rm -rf xmgredit5-dist
	mkdir xmgredit5-dist
	/bin/cp -p Makefile vers.h newvers.sh version gredit_icon.h $(SRCS) $(INCS) $(PARSSRCS) xmgredit5-dist
	tar cvf xmgredit5.tar xmgredit5-dist
	gzip xmgredit5.tar

tarfile: 
	tar cvf xmgredit5.tar xmgredit5
	gzip xmgredit5.tar

rcs: Makefile $(SRCS) $(INCS) $(PARSSRCS)
	ci -l $?
	touch rcs

rel: Makefile $(SRCS) $(INCS) $(PARSSRCS)
	-ci -l $?
	touch rcs
	grep Id:  Makefile $(SRCS) $(INCS) $(PARSSRCS) > release.tmp

install: xmgredit5
	cp -p xmgredit5 $(INSTALLDIR)/bin/xmgredit5

lint:
	lint  -axv -wk -Nn10000 -Nd10000 $(SRCS) $(PARSSRCS)

getrcs:
	strings xmgredit | grep Id

clean:
	/bin/rm *.o xmgredit5

