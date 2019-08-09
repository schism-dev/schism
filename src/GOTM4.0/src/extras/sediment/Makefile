#$Id: Makefile,v 1.3 2004-01-12 08:51:01 lars Exp $
#
# Makefile to build the Sediment Module - libsediment.a
#

include ../../Rules.make

LIB	= $(LIBDIR)/libsediment$(buildtype).a

DOCSRC	= sediment.F90

OBJ   = \
${LIB}(sediment.o)

all: ${OBJ}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/sediment.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2001 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
