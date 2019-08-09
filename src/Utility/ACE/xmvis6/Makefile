#
# ACE/vis - Visualization of Flow and Transport
#
# Paul J. Turner and Antonio M. Baptista
#
# Copyright 1990-2003 Oregon Health and Science University
# All Rights Reserved
#
# 

# $Id: Makefile,v 1.11 2007/03/13 15:34:15 pturner Exp $
#
# make for VIS
#

# xmvis - Motif with X

INSTALLDIR = /usr/local/ace

###################################################
# 'make rcs' to check everything in.
# 'make distrib' to build a clean tree containg source.
# 'make tarfile' to bundle everything up for shipment off-site.
###################################################

###################################################
# to compile a particular version set the following
# items - LIBS, CFLAGS, CC
###################################################

# LINUX and Cygwin
#LIBS = -L/usr/local/ace/lib -L/usr/local/lib -lgd -lfpp ../elio/libelio.a -L/usr/X11R6/lib64   -lXm -lXp -lXt -lXext -lXpm -lX11 -lICE -lSM -lm
LIBS = -L../lib -lgd -lfpp -lelio -L/usr/X11R6/lib64 -lXm -lXp -lXt -lXext -lXpm -lX11 -lICE -lSM -lm
#Chaged to local include
INCLUDES = -I/usr/X11R6/include -I../include #-I/usr/local/include

# for Linux on Alpha for reading binary date written on Intel
CFLAGS = -g -O $(INCLUDES) -mieee
# for Linux on Intel.
CFLAGS = -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -g $(INCLUDES)

CC = gcc

#graphics.h
#externs.h
#hersh.h
#image_f_io.c

# core routines
#
OBJS = gram.o\
	vers.o\
	main.o\
	cedit.o\
        defaults.o\
        tideutils.o\
	stations.o\
        graphutils.o\
        graphutils2.o\
	drawticks.o\
	gridutils.o\
	boundutils.o\
	allogrid.o\
	allobound.o\
	alloflow.o\
	allodrogs.o\
	gridio.o\
	gridtio.o\
	flowio.o\
	concio.o\
	drogsio.o\
        drawobjs.o\
	drawgrid.o\
	drawbound.o\
	drawflow.o\
	drawconc.o\
	drawelev.o\
        drawhist.o\
        drawtide.o\
        drawslice.o\
        drawzoom.o\
	drawdrogs.o\
	isol.o\
	isolutils.o\
        actions.o\
        animat.o\
        interp.o\
        objutils.o\
	find.o\
	getparms.o\
	putparms.o\
	stubs.o\
	utils.o\
	draw.o\
	gifd.o\
	ps.o\
	xvlib.o\
	chersh.o

SRCS = main.c\
	vers.c\
	cedit.c\
        tideutils.c\
	stations.c\
        graphutils.c\
        graphutils2.c\
	drawticks.c\
	defaults.c\
        gridio.c\
	gridtio.c\
	gridutils.c\
	drawgrid.c\
	allogrid.c\
	boundutils.c\
	drawbound.c\
	allobound.c\
	flowio.c\
	alloflow.c\
	drawflow.c\
	concio.c\
	drawconc.c\
	drawelev.c\
	drawhist.c\
	drawtide.c\
	drawslice.c\
	drawzoom.c\
	drogsio.c\
	allodrogs.c\
	drawdrogs.c\
	isol.c\
	isolutils.c\
	drawobjs.c\
	actions.c\
        animat.c\
        interp.c\
        objutils.c\
        find.c\
	getparms.c\
	putparms.c\
	stubs.c\
        gram.c\
        utils.c\
	draw.c\
	gifd.c\
	ps.c\
	xvlib.c\
	chersh.c

INCS = globals.h defines.h gram.h xgif.h bitmaps.h symdefs.h motifinc.h vis_icon.h

# for the motif interface
#
MOTIF_SRCS = newwin.c\
	timewin.c\
	configwin.c\
	strwin.c\
	pathlwin.c\
	trackwin.c\
	drogwin.c\
	isolwin.c\
	gridwin.c\
	gridtwin.c\
	boundwin.c\
	tidewin.c\
	teanlwin.c\
	sampflowwin.c\
	velhistwin.c\
	histwin.c\
	slicewin.c\
	zoomwin.c\
	adcircwin.c\
	adc3dwin.c\
	transwin.c\
	newtrans.c\
	sigmatrans.c\
	sigmaztrans.c\
	surfwin.c\
	annotwin.c\
	controlwin.c\
	placewin.c\
	comwin.c\
	graphwin.c\
	framewin.c\
	imagewin.c\
	worldwin.c\
	tickwin.c\
	locwin.c\
	printwin.c\
	locatewin.c\
	statuswin.c\
	parmsiowin.c\
	malerts.c\
	events.c\
	showdwin.c\
	monwin.c\
	helpwin.c\
	gifread.c\
	image_f_io.c\
	xgifload.c\
	motifutils.c

MOTIF_OBJS = newwin.o\
	tidewin.o\
	timewin.o\
	configwin.o\
	strwin.o\
	pathlwin.o\
	trackwin.o\
	drogwin.o\
	isolwin.o\
	gridwin.o\
	gridtwin.o\
	boundwin.o\
	teanlwin.o\
	sampflowwin.o\
	velhistwin.o\
	histwin.o\
	slicewin.o\
	zoomwin.o\
	adcircwin.o\
	adc3dwin.o\
	transwin.o\
	newtrans.o\
	sigmatrans.o\
	sigmaztrans.o\
	surfwin.o\
	annotwin.o\
	controlwin.o\
	placewin.o\
	comwin.o\
	graphwin.o\
	framewin.o\
	imagewin.o\
	worldwin.o\
	tickwin.o\
	locwin.o\
	printwin.o\
	statuswin.o\
	locatewin.o\
	parmsiowin.o\
	malerts.o\
	events.o\
	showdwin.o\
	monwin.o\
	helpwin.o\
	gifread.o\
	image_f_io.o\
	xgifload.o\
	motifutils.o

xmvis6: $(OBJS) xmvis.o $(MOTIF_OBJS)
	$(CC) $(MCFLAGS) xmvis.o  $(OBJS) $(MOTIF_OBJS) -o xmvis6 $(LIBS)

gram.o: gram.c
gram.c: gram.y
	bison -y -d gram.y
	mv y.tab.c gram.c
	mv y.tab.h gram.h


vers.o: $(SRCS) $(MOTIF_SRCS)
	sh newvers.sh
	$(CC) -c vers.c

main.o: $(INCS)
xgvis.o: $(INCS)
xmvis.o: $(INCS)
teanlwin.o: $(INCS)
configwin.o: $(INCS)
pathlwin.o: $(INCS)
drogwin.o: $(INCS)
draw.o: $(INCS)

isol.o: $(INCS)
gridio.o: $(INCS)
gridutils.o: $(INCS)
drawgrid.o: $(INCS)
allogrid.o: $(INCS)
boundio.o: $(INCS)
boundutils.o: $(INCS)
drawbound.o: $(INCS)
allobound.o: $(INCS)
flowio.o: $(INCS)
alloflow.o: $(INCS)
allovel.o: $(INCS)
alloelev.o: $(INCS)
drawflow.o: $(INCS)
drogsio.o: $(INCS)
allodrogs.o: $(INCS)
drawdrogs.o: $(INCS)
alerts.o: $(INCS)
find.o: $(INCS)
utils.o: $(INCS)


distrib:
	-/bin/rm -rf xmvis6-dist
	mkdir xmvis6-dist
	/bin/cp -p Makefile xmvis.c vers.h version newvers.sh externs.h gram.y $(INCS) $(SRCS) $(MOTIF_SRCS) xmvis6-dist
	tar cvf xmvis6.tar xmvis6-dist
	gzip xmvis6.tar

rcs: Makefile xmvis.c gram.y $(SRCS) $(INCS) $(MOTIF_SRCS)
	ci -l $?
	touch rcs

indent: $(SRCS) $(INCS) $(MOTIF_SRCS)
	indent -kr -l256 $?

rel: Makefile xmvis.c $(SRCS) $(INCS) $(MOTIF_SRCS)
	-ci -l $?
	touch rcs
	grep Id: Makefile xmvis.c $(SRCS) $(INCS) $(MOTIF_SRCS) > release.tmp

clean:
	/bin/rm *.o xmvis6

install:
	/bin/cp -p xmvis6 $(INSTALLDIR)/bin/xmvis6
