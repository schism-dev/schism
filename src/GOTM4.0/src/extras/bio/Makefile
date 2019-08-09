#$Id: Makefile,v 1.14 2007-04-18 07:24:51 kbk Exp $
#
# Makefile to build the Bio Module - libbio.a
#

include ../../Rules.make

LIB	= $(LIBDIR)/libbio$(buildtype).a

DOCSRC	=  bio.F90 \
bio_var.F90 bio_template.F90 bio_npzd.F90 bio_iow.F90 \
          bio_sed.F90 bio_fasham.F90 \
          bio_save.F90

OBJ   = \
${LIB}(bio_var.o)		\
${LIB}(bio_fluxes.o)		\
${LIB}(bio_template.o)		\
${LIB}(bio_npzd.o)		\
${LIB}(bio_iow.o)		\
${LIB}(bio_mab.o)		\
${LIB}(bio_fasham.o)		\
${LIB}(bio_sed.o)		\
${LIB}(bio_save.o)		\
${LIB}(bio.o)

all: ${OBJ}
	$(MOVE_MODULES_COMMAND)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../../doc/bio.tex 
	touch doc

clean:
	$(RM) ${LIB} doc

realclean: clean
	$(RM) *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 2003 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
