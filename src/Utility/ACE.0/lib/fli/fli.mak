#
# makefile for FBM2FLI
#

#
# for WIN32
#

XDK=c:\home\pturner\exceed\xdk
INCLUDES = -I $(XDK)\include

CC = cl

.c.obj:
	$(CC) $(INCLUDES) -c $*.c

OBJ = \
	fppmain.obj\
	fppmake.obj\
	fppframe.obj\
	fppbrun.obj\
	fppdelta.obj\
	fppcolor.obj\
	fppfile.obj

#
# Library
#
lib: $(OBJ)
	lib /out:fli.lib $(OBJ)

#
# Dependecies
#
#$(OBJ) : fpfli.h
