/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/*
 *
 * allobound.c - allocate a boundary
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allobound.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

int Allocate_boundary(int gridno, int bno, int npts, int btype)
{
    Boundary *b = grid[gridno].boundaries;
    if (bno < 0) {
	return 1;
    }
    b[bno].active = ON;
    b[bno].nbpts = npts;
    b[bno].type = btype;
    b[bno].boundx = (double *) calloc(npts, sizeof(double));
    b[bno].boundy = (double *) calloc(npts, sizeof(double));
    b[bno].nodes = (int *) calloc(npts, sizeof(int));
    b[bno].openclose = (int *) calloc(npts, sizeof(int));
    return 0;
}

void Reallocate_boundary(int gridno, int bno, int npts)
{
    Boundary *b = grid[gridno].boundaries;
    if (bno < 0) {
	return;
    }
    b[bno].nbpts = npts;
    b[bno].boundx = (double *) realloc(b[bno].boundx, npts * sizeof(double));
    b[bno].boundy = (double *) realloc(b[bno].boundy, npts * sizeof(double));
    b[bno].nodes = (int *) realloc(b[bno].nodes, npts * sizeof(int));
    b[bno].openclose = (int *) realloc(b[bno].nodes, npts * sizeof(int));
}

void Allocate_boundaries(int gridno, int nb)
{
    if (nb <= 0 || nb == grid[gridno].nbounds) {
	return;
    }
    grid[gridno].nbounds = nb;
    if (grid[gridno].boundaries == NULL) {
	grid[gridno].boundaries = (Boundary *) calloc(nb, sizeof(Boundary));
    } else {
	grid[gridno].boundaries = (Boundary *) realloc(grid[gridno].boundaries, nb * sizeof(Boundary));
    }
}

void Free_boundary(int gridno, int bno)
{
    Boundary *b = grid[gridno].boundaries;
    if (bno < 0) {
	return;
    }
    b[bno].nbpts = 0;
    b[bno].type = 0;
    b[bno].active = OFF;
    if (b[bno].boundx != NULL) {
	free(b[bno].boundx);
    }
    if (b[bno].boundy != NULL) {
	free(b[bno].boundy);
    }
    if (b[bno].nodes != NULL) {
	free(b[bno].nodes);
    }
    if (b[bno].openclose != NULL) {
	free(b[bno].openclose);
    }
}

void Free_boundaries(int gridno)
{
    int i;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(gridno, i);
    }
    if (grid[gridno].boundaries != NULL) {
	free(grid[gridno].boundaries);
    }
    grid[gridno].nbounds = 0;
}

int Allocate_boundary2(int boundno, int bno, int npts, int btype)
{
    Boundary *b = bounds[boundno].boundaries;
    if (bno < 0) {
	return 1;
    }
    b[bno].active = ON;
    b[bno].nbpts = npts;
    b[bno].type = btype;
    b[bno].boundx = (double *) calloc(npts, sizeof(double));
    b[bno].boundy = (double *) calloc(npts, sizeof(double));
    b[bno].openclose = (int *) calloc(npts, sizeof(int));
    return 0;
}

void Reallocate_boundary2(int boundno, int bno, int npts)
{
    Boundary *b = bounds[boundno].boundaries;
    if (bno < 0) {
	return;
    }
    b[bno].nbpts = npts;
    b[bno].boundx = (double *) realloc(b[bno].boundx, npts * sizeof(double));
    b[bno].boundy = (double *) realloc(b[bno].boundy, npts * sizeof(double));
    b[bno].openclose = (int *) realloc(b[bno].nodes, npts * sizeof(int));
}

void Allocate_boundaries2(int boundno, int nb)
{
    if (nb <= 0) {
	return;
    }
    bounds[boundno].nbounds = nb;
    if (bounds[boundno].boundaries == NULL) {
	bounds[boundno].boundaries = (Boundary *) calloc(nb, sizeof(Boundary));
    } else {
	bounds[boundno].boundaries = (Boundary *) realloc(bounds[boundno].boundaries, nb * sizeof(Boundary));
    }
}

void Free_boundary2(int boundno, int bno)
{
    Boundary *b = grid[boundno].boundaries;
    if (bno < 0) {
	return;
    }
    b[bno].nbpts = 0;
    b[bno].type = 0;
    b[bno].active = OFF;
    if (b[bno].boundx != NULL) {
	free(b[bno].boundx);
    }
    if (b[bno].boundy != NULL) {
	free(b[bno].boundy);
    }
    if (b[bno].nodes != NULL) {
	free(b[bno].nodes);
    }
    if (b[bno].openclose != NULL) {
	free(b[bno].openclose);
    }
}

void Free_boundaries2(int boundno)
{
    int i;
    for (i = 0; i < bounds[boundno].nbounds; i++) {
	Free_boundary2(boundno, i);
    }
    if (bounds[boundno].boundaries != NULL) {
	free(bounds[boundno].boundaries);
    }
    bounds[boundno].nbounds = 0;
}
