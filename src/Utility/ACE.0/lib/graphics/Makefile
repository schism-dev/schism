#
# Makefile for graphics routines
#
# 

LIBDIR = "/usr/local/ace/lib"

# LINUX
INCLUDES = -I/usr/X11R6/include
CC = gcc -D_NO_PROTO $(INCLUDES)
CFLAGS = -I/usr/local/include

SRCS = chersh.c draw.c image_f_io.c gifd.c ps.c xvlib.c
OBJS = chersh.o draw.o image_f_io.o gifd.o ps.o xvlib.o
INCS = externs.h hersh.h patterns.h special.h symdefs.h graphics.h

lib: $(OBJS)
	-/bin/rm libacegraph.a
	ar q libacegraph.a $(OBJS)

distrib:
	tar cvf acegraph.tar Makefile $(SRCS) $(INCS)

install:
	/bin/cp -p libacegraph.a $(LIBDIR)

clean:
	/bin/rm libacegraph.a $(OBJS)

help:
	echo 'make lib | distrib | help'
