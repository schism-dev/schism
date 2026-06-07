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
 * allogrid.c - allocate a grid
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allogrid.c,v 1.4 2003/10/17 02:38:47 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

int Allocate_grid(int gridno, int nmel, int nmnp)
{
    int i;
    grid[gridno].nmel = nmel;
    grid[gridno].nmnp = nmnp;
    grid[gridno].active = ON;
    grid[gridno].nbounds = 0;
    grid[gridno].boundaries = NULL;
    grid[gridno].ellist = NULL;
    grid[gridno].xord = (double *) malloc(nmnp * sizeof(double));
    grid[gridno].yord = (double *) malloc(nmnp * sizeof(double));
    grid[gridno].depth = (double *) malloc(nmnp * sizeof(double));
    grid[gridno].icon = (Element *) malloc(nmel * sizeof(Element));
    return 1;
}

void Free_grid(int gridno)
{
    int i;
    if (grid[gridno].boundaries != NULL) {
	free(grid[gridno].boundaries);
    }
    if (grid[gridno].xord != NULL) {
	free(grid[gridno].xord);
    }
    if (grid[gridno].yord != NULL) {
	free(grid[gridno].yord);
    }
    if (grid[gridno].icon != NULL) {
	free(grid[gridno].icon);
    }
    grid[gridno].ellist = NULL;
    grid[gridno].nmel = 0;
    grid[gridno].nmnp = 0;
    grid[gridno].active = OFF;
}

int isactive_grid(int gridno)
{
    return grid[gridno].active == ON;
}

int AllocateGrid(Grid * g, int nmel, int nmnp)
{
    int i;
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->active = ON;
    g->nbounds = 0;
    g->boundaries = NULL;
    g->ellist = NULL;
    g->xord = (double *) malloc(nmnp * sizeof(double));
    g->yord = (double *) malloc(nmnp * sizeof(double));
    g->depth = (double *) malloc(nmnp * sizeof(double));
    g->u = (double *) malloc(nmnp * sizeof(double));
    g->v = (double *) malloc(nmnp * sizeof(double));
    g->icon = (Element *) malloc(nmel * sizeof(Element));
    return 1;
}

int FreeGrid(Grid * g)
{
    int i;
    if (g->xord != NULL) {
	free(g->xord);
    }
    if (g->yord != NULL) {
	free(g->yord);
    }
    if (g->depth != NULL) {
	free(g->depth);
    }
    if (g->u != NULL) {
	free(g->u);
    }
    if (g->v != NULL) {
	free(g->v);
    }
    if (g->icon != NULL) {
	free(g->icon);
    }
    if (g->ellist != NULL) {
	free(g->ellist);
    }
    g->nmel = 0;
    g->nmnp = 0;
    return 1;
}
