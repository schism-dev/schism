# $Id: Makefile,v 1.9 2008/10/06 15:45:46 pturner Exp $
#
# Makefile for xmgr Release 5
#

# Step 1.

# you may need to change the following
# BINDIR - location of the ACE/gr (xmgr) executable
# APPDEFAULTS - directory to put the ACE/gr application defaults file.
# HELPFILE - if not defined, then "ACEgr.html" (shouldn't need to change this).
# HELPVIEWER - if not defined, "Mosaic" - others should work.
##################################
INSTALLDIR = /usr/local/ace
BINDIR = /usr/local/ace/bin
APPDEFAULTS = /usr/lib/X11/app-defaults
HELPFILE = -DGR_HELPFILE='"ACEgr.html"'
HELPVIEWER = -DGR_HELPVIEWER='"netscape"'

# Step 2.

# You need to select the machine type
# uncomment one of the following, make sure the
# others are commented out
##################################

### Linux w X11R6/Motif 1.2
#LIBS = -L/usr/X11R6/lib -L/usr/local/ace/lib -lfpp -lgd -lXm -lXt -lXp -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
#LIBS = -L/usr/X11R6/lib -L/usr/local/lib -L/usr/local/ace/lib -lfpp -L/usr/local/ace/netcdf/lib -lgd -lXm -lXt -lXp -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
# netcdf, pgsql
#LIBS = -L/usr/X11R6/lib -L/usr/local/ace/lib -lfpp -L/usr/local/pgsql/lib -lpq -lgd -lXm -lXt -lXp -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
#LIBS = -L/usr/X11R6/lib64 -L/usr/local/ace/lib -lfpp -L/usr/local/pgsql/lib -lpq -L/usr/local/lib -lgd -lXm -lXt -lXp -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
LIBS = -L/usr/X11R6/lib64 -L../lib -L/usr/local/intel64/broadwell/gcc/pgsql-10.2-1/lib -L/usr/local/lib -lfpp -lgd -lXm -lXt -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
#LIBS = -L/usr/X11R6/lib -L/usr/local/ace/lib -lfpp -L/usr/local/pgsql/lib -lpq -L/usr/local/lib -lgd -lXm -lXt -lXp -lXext -lX11 -lICE -lSM -lICE -dynamic -lXmu -ldl -lm
#INCLUDES = -I/usr/X11R6/include -I/usr/local/ace/include
#INCLUDES = -I/usr/X11R6/include -I/usr/local/ace/netcdf/include -I/usr/local/ace/include
# netcdf, pgsql
INCLUDES = -I/usr/X11R6/include -I..//include -I/usr/local/intel64/broadwell/gcc/pgsql-10.2-1/include -I/usr/local/include
##INCLUDES = -I/usr/X11R6/include -I/usr/local/ace/include -I/usr/local/intel64/broadwell/gcc/pgsql-10.2-1/include -I/usr/local/include
CC = gcc -g
# Compile in support for postgres.
#PGSQL = -DPGSQL

# Step 3.

# Limited NetCDF/HDF support
#################################
# support for netCDF -> implemented
# Uncomment the following 3 defines for netCDF support
# adjust paths to suit local conditions
#
#NETCDF = -DHAVE_NETCDF
#NETCDFINCLUDES = -I/usr/local/netcdf/include
#NETCDFLIBS = -L/usr/local/netcdf/lib -lnetcdf

#################################

# Step 4.

# Install a printer
#################################
#
# The default is PostScript landscape using "lpr -h" to spool
#
#
# To install the default hardcopy device uncomment ONE of the
# following (HARDCOPYDEV lines), the hardcopy device can also 
# be set using environment variables, and on the command line:
#
# PostScript landscape [device number 1]
HARDCOPYDEV = -DHDEV=GR_PS_L 

# Step 5.

# To install the printer for the hardcopy devices above
# set the following to the appropriate value for your site.
# These can also be set by environment variables. The current
# definitions of these variables are the same as the definitions
# found in ps.c (PostScript), and mif.c (Maker Interchange format).
#
# NOTE: if you use lp (as on Solaris 2.x, IRIX, and other System V systems),
#       add the -c option or xmgr will erase the temporary file
#       before lp can spool it -> "lp -c "
# 
# For Postscript
# 
PS_PRSTR = -DPS_PRSTR='"lp -c "'

# Step 6.

# Xmgr uses 4 scratch arrays, designated in the function that
# interprets expressions as a, b, c, and d. MAXARR is the amount
# of memory that will be allocated for these arrays and at present
# there is no provision for dynamically extending the length of these 
# arrays. This sets a limit on the length of a set inasmuch
# as no operation between sets and the scratch arrays should exceed 
# this amount. MAXARR does NOT set a limit on the size of a set, on our 
# Sun4 (32Mb) we run MAXARR = 60000, I've plotted data sets in excess 
# of 700k points.
#
MAXARR = -DMAXARR=60000

###############################
# End installation parameters #
###############################

CFLAGS = $(XBAE_INCLUDES) $(INCLUDES) $(MAXARR) $(HARDCOPYDEV) \
	$(PGSQL) $(PS_PRSTR) $(MIF_PRSTR) $(HP_PRSTR1) \
	$(LEAF_PRSTR) $(NETCDF) $(MFHDF) \
	$(NETCDFINCLUDES) \
	$(HELPFILE) $(HELPVIEWER) $(HAS_XBAE)

SRCS =  main.c\
	xmgr.c\
	events.c\
	malerts.c\
	fileswin.c\
	blockwin.c\
	eblockwin.c\
	printwin.c\
	drawwin.c\
	comwin.c\
	compwin.c\
	setwin.c\
	hotwin.c\
	ptswin.c\
	editpwin.c\
	regionwin.c\
	statuswin.c\
	symwin.c\
	tickwin.c\
	graphwin.c\
	worldwin.c\
	strwin.c\
	labelwin.c\
	locatewin.c\
	framewin.c\
	pagewin.c\
	monwin.c\
	helpwin.c\
	miscwin.c\
	isolwin.c\
	dbwin.c\
	viswin.c\
	visutils.c\
	stubs.c\
	plotone.c\
	files.c\
	utils.c\
	drawticks.c\
	fit.c\
	fourier.c\
	io.c\
	monitors.c\
	isolines.c\
	graphutils.c\
	graphutils2.c\
	setutils.c\
	setprops.c\
	regionutils.c\
	objutils.c\
	computils.c\
	defaults.c\
	checkon.c\
	params.c\
	getparms.c\
	compute.c\
	xvlib.c\
	image_f_io.c\
	draw.c\
	ps.c\
	mif.c\
	gifd.c\
	gourad5.c\
	chersh.c\
	pars.c\
	motifutils.c

OBJS =  main.o\
	xmgr.o\
	events.o\
	malerts.o\
	fileswin.o\
	blockwin.o\
	eblockwin.o\
	printwin.o\
	drawwin.o\
	comwin.o\
	compwin.o\
	setwin.o\
	hotwin.o\
	ptswin.o\
	editpwin.o\
	regionwin.o\
	statuswin.o\
	symwin.o\
	tickwin.o\
	graphwin.o\
	worldwin.o\
	strwin.o\
	labelwin.o\
	locatewin.o\
	framewin.o\
	pagewin.o\
	monwin.o\
	helpwin.o\
	miscwin.o\
	isolwin.o\
	dbwin.o\
	viswin.o\
	visutils.o\
	stubs.o\
	plotone.o\
	files.o\
	utils.o\
	drawticks.o\
	fit.o\
	fourier.o\
	io.o\
	monitors.o\
	isolines.o\
	graphutils.o\
	graphutils2.o\
	setutils.o\
	setprops.o\
	regionutils.o\
	objutils.o\
	computils.o\
	defaults.o\
	checkon.o\
	params.o\
	getparms.o\
	compute.o\
	xvlib.o\
	image_f_io.o\
	draw.o\
	ps.o\
	mif.o\
	gifd.o\
	gourad5.o\
	chersh.o\
	pars.o\
	motifutils.o

D3DINCS = defs.h

# for printing and making a distribution
#
INCS = defines.h vis.h externs.h globals.h motifinc.h pars.h\
	 xmgr_icon.h symdefs.h isolines.h\
	 xprotos.h noxprotos.h draw.h bitmaps.h\
	patterns.h patchlevel.h elcirc.h
HERSHINCS = hersh.h symdef.h special.h
DISTDIR = xmgr5-dist
HELPFILES = doc/ACEgr.html\
	doc/README.html\
	doc/graphs.html\
	doc/points.html\
	doc/trans.html\
	doc/commands.html\
	doc/main.html\
	doc/regions.html\
	doc/view.html\
	doc/graphops.html\
	doc/page.html\
	doc/status.html\
	doc/sets.html

xmgr5: $(OBJS)
	$(CC) $(OBJS) -o xmgr5 $(NETCDFLIBS) $(XBAE_LIBS) $(LIBS)

local: $(OBJS) viswin.o visutils.o
	$(CC) $(OBJS) viswin.o visutils.o -o xmgr5 $(NETCDFLIBS) $(LIBS) 

pars.o: pars.c
pars.c: pars.y
	bison -y -d pars.y
	mv y.tab.c pars.c
	mv y.tab.h pars.h

install:
	install -s -m 751 xmgr5 $(INSTALLDIR)/bin
#	install -c -m 644 XMgr.ad $(APPDEFAULTS)/XMgr

lint:
	lint  $(INCLUDES) $(SRCS)

clean: 
	/bin/rm xmgr5 $(OBJS)

distrib:
	/bin/rm -rf $(DISTDIR)
	mkdir $(DISTDIR)
	cp -p Makefile $(SRCS) pars.y $(INCS) $(HERSHINCS) $(DISTDIR)

localdistrib:
	mkdir $(DISTDIR)
	mkdir $(DISTDIR)/src
	mkdir $(DISTDIR)/bin
	cp -p Makefile xmgr.mak $(SRCS) pars.y $(INCS) $(HERSHINCS) $(DISTDIR)/src
	cp -p Makefile $(DISTDIR)/src/Makefile.orig
	cp -p XMgr.ad $(DISTDIR)/src/XMgr.ad

tarfile:
	tar cvf $(DISTDIR).tar $(DISTDIR)
	gzip $(DISTDIR).tar

vi:
	vi Makefile README $(SRCS) $(INCS) hersh.h special.h symdef.h

wc:
	wc Makefile README $(SRCS) $(INCS) hersh.h special.h symdef.h

dorcs: Makefile patchlevel.h pars.y $(SRCS) $(INCS)
	ci -l $?
	touch dorcs

protos: $(SRCS)
	cextract $(SRCS) >> prototypes.h

