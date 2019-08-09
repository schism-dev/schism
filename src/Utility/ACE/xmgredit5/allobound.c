/*
 * ACE/gredit - 2d finite element grid generation
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *                      All Rights Reserved.
 *
 */

/*
 *
 * allobound.c - allocate a boundary
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allobound.c,v 1.3 2007/02/21 00:21:20 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

/* allobound.c */
int Allocate_boundary(int bno, int npts, int btype);
int copy_boundary(int bno);
void copy_boundaries(int g1, int g2);
void Reallocate_boundary(int bno, int npts);
void Free_boundary(int bno);
void Free_boundaries(int gridno);
int nboundaries(int gridno);
int nextboundary(void);

/*
 * Allocate memory for boundary bno of length npts, 
 * btype == 0 => external boundary
 */
int Allocate_boundary(int bno, int npts, int btype)
{
    if (bno < 0) {
	return 0;
    }
    boundary[bno].nbpts = npts;
    boundary[bno].boundtype = btype;
    boundary[bno].bactive = 1;
    boundary[bno].boundx = (double *) calloc(npts, sizeof(double));
    boundary[bno].boundy = (double *) calloc(npts, sizeof(double));
    boundary[bno].nodes = (int *) calloc(npts, sizeof(int));
    boundary[bno].btype = (int *) calloc(npts, sizeof(int));
    return 1;
}

/*
 * Copy boundary bno and return a handle to the copy
 */
int copy_boundary(int bno)
{
    int i, ib;
    if (bno < 0) {
	return 0;
    }
    if ((ib = nextboundary()) == -1) {
	return -1;
    }
    Allocate_boundary(ib, boundary[bno].nbpts, boundary[bno].boundtype);
    for (i = 0; i < boundary[bno].nbpts; i++) {
	boundary[ib].boundx[i] = boundary[bno].boundx[i];
	boundary[ib].boundy[i] = boundary[bno].boundy[i];
    }
    return ib;
}

/*
 * Copy a collection of boundaries from one grid to another
 */
void copy_boundaries(int g1, int g2)
{
    int i, ib;
    void Free_boundaries(int gridno);
    if (grid[g2].nbounds) {
	Free_boundaries(g2);
    }
    grid[g2].nbounds = grid[g1].nbounds;
    for (i = 0; i < grid[g1].nbounds; i++) {
	grid[g2].boundaries[i] = copy_boundary(grid[g1].boundaries[i]);
    }
}

/*
 * Resize a boundary
 */
void Reallocate_boundary(int bno, int npts)
{
    if (bno < 0) {
	return;
    }
    boundary[bno].nbpts = npts;
    boundary[bno].boundx = (double *) realloc(boundary[bno].boundx, npts * sizeof(double));
    boundary[bno].boundy = (double *) realloc(boundary[bno].boundy, npts * sizeof(double));
    boundary[bno].nodes = (int *) realloc(boundary[bno].nodes, npts * sizeof(int));
    boundary[bno].btype = (int *) realloc(boundary[bno].nodes, npts * sizeof(int));
}

/*
 * Free a boundary
 */
void Free_boundary(int bno)
{
    if (bno < 0) {
	return;
    }
    boundary[bno].nbpts = 0;
    boundary[bno].boundtype = 0;
    boundary[bno].bactive = 0;
    if (boundary[bno].boundx != NULL) {
	free(boundary[bno].boundx);
    }
    if (boundary[bno].boundy != NULL) {
	free(boundary[bno].boundy);
    }
    if (boundary[bno].nodes != NULL) {
	free(boundary[bno].nodes);
    }
    if (boundary[bno].btype != NULL) {
	free(boundary[bno].btype);
    }
}

/*
 * Free a collection of boundaries
 */
void Free_boundaries(int gridno)
{
    int i;

    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(grid[gridno].boundaries[i]);
    }
}


/*
 * Return the number of boundaries attached to a grid
 */
int nboundaries(int gridno)
{
    return grid[gridno].nbounds;
}

/*
 * Return the next available boundary
 */
int nextboundary(void)
{
    int i;

    for (i = 0; i < MAXBOUNDS; i++) {
	if (boundary[i].bactive == 0) {
	    return i;
	}
    }
    return -1;
}
