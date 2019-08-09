/*
 * The author of this software is Steven Fortune.  Copyright (c) 1994 by AT&T
 * Bell Laboratories.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.- */

#
#include "vdefs.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef lint
static char RCSid[] = "$Id: vmemory.c,v 1.4 2007/02/21 00:21:21 pturner Exp $";
#endif

/* vmemory.c */
void freeinit(struct Freelist *fl, int size);
char *getfree(struct Freelist *fl);
void makefree(struct Freenode *curr, struct Freelist *fl);
char *myalloc(unsigned int n);

void freeinit(struct Freelist * fl, int size)
{
    fl->head = (struct Freenode *) NULL;
    fl->nodesize = size;
}

char *getfree(struct Freelist * fl)
{
    int i;
    struct Freenode *t;

    if (fl->head == (struct Freenode *) NULL) {
	t = (struct Freenode *) myalloc(sqrt_nsites * fl->nodesize);
	for (i = 0; i < sqrt_nsites; i += 1)
	    makefree((struct Freenode *) ((char *) t + i * fl->nodesize), fl);
    };
    t = fl->head;
    fl->head = (fl->head)->nextfree;
    return ((char *) t);
}

void makefree(struct Freenode * curr, struct Freelist * fl)
{
    curr->nextfree = fl->head;
    fl->head = curr;
}

int total_alloc, called;
char *myalloc(unsigned int n)
{
    char *t;
    extern char *ptrs[];

    if ((t = (char *) malloc(n)) == (char *) 0) {
	fprintf(stderr, "Insufficient memory processing site %d (%d bytes in use)\n",
		siteidx, total_alloc);
	exit(1);
    };
    total_alloc += n;
    ptrs[called] = (char *) t;
    called++;
    return (t);
}
