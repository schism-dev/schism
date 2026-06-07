#
# makefile for FBM2FLI
#
CC = gcc #-64

CFLAGS = -O
CFLAGS = -O -I/usr/X11R6/include
.c.o:
	$(CC) $(CFLAGS) -c $*.c

OBJ = \
	fppmain.o \
	fppmake.o \
	fppframe.o \
	fppbrun.o \
	fppdelta.o \
	fppcolor.o \
	fppfile.o

#
# Library
#
lib: $(OBJ)
	-/bin/rm libfpp.a
	ar q libfpp.a $(OBJ)

clean:
	/bin/rm libfpp.a $(OBJ)

#
# Dependecies
#
$(OBJ) : fpfli.h
