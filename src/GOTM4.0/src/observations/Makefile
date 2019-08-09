#$Id: Makefile,v 1.6 2007-01-04 12:08:12 kbk Exp $
#
# Makefile to build the the Observtions library - libobservations.a
#

include ../Rules.make

LIB	= $(LIBDIR)/libobservations$(buildtype).a

DOCSRC	= \
observations.F90 	\
get_s_profile.F90       \
get_t_profile.F90       \
get_o2_profile.F90      \
get_ext_pressure.F90    \
get_int_pressure.F90    \
read_extinction.F90     \
get_w_adv.F90           \
get_zeta.F90            \
get_wave.F90            \
get_vel_profile.F90     \
get_eps_profile.F90	\
analytical_profile.F90	\
const_NNT.F90		\
const_NNS.F90             


MOD	= \
${LIB}(observations.o)

OBJ	= \
${LIB}(get_s_profile.o)		\
${LIB}(get_t_profile.o)		\
${LIB}(get_o2_profile.o)	\
${LIB}(get_ext_pressure.o)	\
${LIB}(get_int_pressure.o)	\
${LIB}(read_extinction.o)	\
${LIB}(get_w_adv.o)		\
${LIB}(get_zeta.o)		\
${LIB}(get_wave.o)		\
${LIB}(get_vel_profile.o)	\
${LIB}(get_eps_profile.o)	\
${LIB}(analytical_profile.o)	\
${LIB}(const_NNT.o)             \
${LIB}(const_NNS.o)

#kbkall: $(MOD)
all: $(MOD) $(OBJ)
	$(MOVE_MODULES_COMMAND)

$(OBJ): $(MOD)

doc:    $(DOCSRC)
	$(PROTEX) $(DOCSRC) > ../../doc/observations.tex 
	touch doc

clean:
	$(RM) -f ${LIB} $(MODDIR)/observations.{m,mod} doc

realclean: clean
	$(RM) -f *.o *~

distclean: realclean

#-----------------------------------------------------------------------
# Copyright (C) 1999 - Hans Burchard and Karsten Bolding (BBH)         !
#-----------------------------------------------------------------------
