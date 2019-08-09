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



#ifndef lint
static char RCSid[] = "$Id: vtriang.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#
#include <stdio.h>
#include "vdefines.h"

#include "defines.h"
#include "globals.h"

/* vtriang.c */
int triang(int bno);
int scomp(struct Point *s1, struct Point *s2);
struct Site *nextone(void);
void readsites1(double *bx, double *by, int n);
struct Site *readone(void);

static int nb;
static double *bx, *by;

char *ptrs[5000];

int triang(int bno)
{
    int c, i;
    struct Site *(*next) ();
    extern int total_alloc, called;
    sorted = 0;
    triangulate0 = 1;
    plot = 0;
    debug = 0;
    total_alloc = 0;
    called = 0;
    tmptable = NULL;
    nels = 0;
    freeinit(&sfl, sizeof *sites);
    nb = build[bno].nbuild;
    bx = build[bno].bx;
    by = build[bno].by;
    readsites1(bx, by, nb);
    next = nextone;
    siteidx = 0;
    geominit();
    voronoi0(triangulate0, next);
    for (i = 0; i < called; i++) {
	free(ptrs[i]);
    }
    return 0;
}

/* sort sites on y, then x, coord */
int scomp(struct Point * s1, struct Point * s2)
{
    if (s1->y < s2->y)
	return (-1);
    if (s1->y > s2->y)
	return (1);
    if (s1->x < s2->x)
	return (-1);
    if (s1->x > s2->x)
	return (1);
    return (0);
}

/* return a single in-storage site */
struct Site *nextone(void)
{
    struct Site *s;

    if (siteidx < nsites) {
	s = &sites[siteidx];
	siteidx += 1;
	return (s);
    } else
	return ((struct Site *) NULL);
}

/* read all sites, sort, and compute xmin, xmax, ymin, ymax */
void readsites1(double *bx, double *by, int n)
{
    int i;
    nsites = n;
    sites = (struct Site *) malloc(n * sizeof(*sites));
    ptrs[0] = (char *) sites;
    for (i = 0; i < n; i++) {
	sites[i].coord.x = bx[i];
	sites[i].coord.y = by[i];
	sites[i].sitenbr = i;
	sites[i].refcnt = 0;
    }
    qsort(sites, nsites, sizeof *sites, scomp);
    vxmin = sites[0].coord.x;
    vxmax = sites[0].coord.x;
    for (i = 1; i < nsites; i += 1) {
	if (sites[i].coord.x < vxmin)
	    vxmin = sites[i].coord.x;
	if (sites[i].coord.x > vxmax)
	    vxmax = sites[i].coord.x;
    };
    vymin = sites[0].coord.y;
    vymax = sites[nsites - 1].coord.y;
}

/* read one site */
struct Site *readone(void)
{
    struct Site *s;

    s = (struct Site *) getfree(&sfl);
    s->refcnt = 0;
    s->sitenbr = siteidx;
    siteidx += 1;
    if (scanf("%lf %lf", &(s->coord.x), &(s->coord.y)) == EOF)
	return ((struct Site *) NULL);
    return (s);
}
