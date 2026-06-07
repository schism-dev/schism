# Makefile for libelio.a and supporting codes.
#
INSTALLDIR = ./

all: lib ${EXECS}

OBJS = elio.o

CFLAGS = -I./
LIB = .

CC = gcc -g
LDFLAGS = -L${LIB} -lelio -lm

lib: $(OBJS) elio.h elio.c
	$(CC) $(CFLAGS) -c -g elio.c
	-/bin/rm libelio.a
	ar q libelio.a $(OBJS)
