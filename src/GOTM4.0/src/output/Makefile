#$Id: Makefile,v 1.4 2005-07-06 14:22:40 kbk Exp $
#
# Makefile to build the the output library - liboutput.a
#

include ../Rules.make

LIB	= $(LIBDIR)/liboutput$(buildtype).a

DOCSRC	= output.F90 asciiout.F90 ncdfout.F90


MOD	= \
${LIB}(asciiout.o)		\
${LIB}(ncdfout.o)		\
${LIB}(output.o)

#${LIB}(gradsout.o)		\

all: $(MOD)
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/output.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/{asciiout,ncdfout,output}.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
