#$Id: Makefile,v 1.8 2005-07-06 15:00:36 kbk Exp $
#
# Makefile to build the the Util library - libutil.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libutil$(buildtype).a

DOCSRC	=  util.F90                                                           \
           diff_center.F90 diff_face.F90  adv_center.F90                      \
           tridiagonal.F90 eqstate.F90  gridinterpol.F90  convert_fluxes.F90  \
           time.F90


MOD   = \
${LIB}(util.o)			\
${LIB}(time.o)			\
${LIB}(tridiagonal.o)		\
${LIB}(eqstate.o)

OBJ   = \
${LIB}(gridinterpol.o)		\
${LIB}(convert_fluxes.o)	\
${LIB}(diff_center.o)           \
${LIB}(diff_face.o)             \
${LIB}(adv_center.o)            

all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/util.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/{util,eqstate,time,tridiagonal}.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
