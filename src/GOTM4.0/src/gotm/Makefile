#$Id: Makefile,v 1.4 2004-01-12 08:51:16 lars Exp $
#
# Makefile to build the the main library - libgotm.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libgotm$(buildtype).a

DOCSRC	= \
main.F90  \
gotm.F90

MOD	= \
${LIB}(gotm.o)


OBJ   = \
main.o

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/gotm.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/gotm.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2000 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
