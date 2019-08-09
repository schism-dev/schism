#$Id: Makefile,v 1.3 2004-01-12 08:49:22 lars Exp $
#
# Makefile to build the the airsea interaction library - libairsea.a
#

include ../Rules.make

LIB	= ../../libairsea$(buildtype).a
LIB	:= $(LIBDIR)/libairsea$(buildtype).a

DOCSRC	= \
airsea.F90

MOD	= \
${LIB}(airsea.o)

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

$(OBJ): $(MOD)

tst_airsea: $(MOD)
	$(FC) -o tst_airsea tst_airsea.F90 $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/airsea.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/airsea.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
