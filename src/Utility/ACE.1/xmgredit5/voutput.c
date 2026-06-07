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

#ifndef lint
static char RCSid[] = "$Id: voutput.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include "vdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

void out_bisector(struct Edge * e)
{
}


void out_ep(struct Edge * e)
{
/*
        printf("e %d %d: %d", nedges, nels + 1, e->edgenbr + 1);
        printf(" %d ", e->ep[le] != (struct Site *) NULL ? e->ep[le]->sitenbr  + 1: -1);
        printf("%d\n", e->ep[re] != (struct Site *) NULL ? e->ep[re]->sitenbr  + 1: -1);
*/
}

void out_vertex(struct Site * v)
{
}


void out_site(struct Site * s)
{
}

#define BUFS 1000

void out_triple(struct Site * s1, struct Site * s2, struct Site * s3)
{
    if (tmptable == NULL) {
	tmptable = (Element *) calloc(BUFS, sizeof(Element));
    } else {
	if (nels % BUFS == 0) {
	    tmptable = (Element *) realloc(tmptable, (nels + BUFS) * sizeof(Element));
	}
    }

    tmptable[nels].nl[0] = s3->sitenbr;
    tmptable[nels].nl[1] = s2->sitenbr;
    tmptable[nels++].nl[2] = s1->sitenbr;
}
