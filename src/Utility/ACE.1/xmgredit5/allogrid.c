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
 * allogrid.c - allocate a grid
 *
 */

#ifndef lint
static char RCSid[] = "$Id: allogrid.c,v 1.4 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "globals.h"

/*
 * Allocate a grid and initialize the boundary
 */
int Allocate_grid(int gridno, int nmel, int nmnp)
{
    int i;
    grid[gridno].nmel = nmel;
    grid[gridno].nmnp = nmnp;
    grid[gridno].gactive = 1;
    grid[gridno].bactive = 0;
    grid[gridno].nbounds = 0;
    for (i = 0; i < MAXGRIDBOUNDS; i++) {
	grid[gridno].boundaries[i] = -1;
    }
    grid[gridno].xord = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].yord = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].depth = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    grid[gridno].ellist = (int *) calloc(nmel, sizeof(int));
    grid[gridno].nlist = (int *) calloc(nmnp, sizeof(int));
    grid[gridno].ntype = (int *) calloc(nmnp, sizeof(int));
    grid[gridno].icon = (Element *) calloc(nmel, sizeof(Element));
    return 1;
}

/*
 * Allocate a grid, but not the boundary
 */
int Allocate_grid_only(int gridno, int nmel, int nmnp)
{
    int i;
    grid[gridno].nmel = nmel;
    grid[gridno].nmnp = nmnp;
    grid[gridno].gactive = 1;
    if (grid[gridno].xord == NULL) {
    }
    grid[gridno].xord = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].yord = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].depth = (double *) calloc(nmnp, sizeof(double));
    grid[gridno].nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    grid[gridno].nlist = (int *) calloc(nmnp, sizeof(int));
    grid[gridno].ntype = (int *) calloc(nmnp, sizeof(int));
    grid[gridno].icon = (Element *) calloc(nmel, sizeof(Element));
    grid[gridno].ellist = (int *) calloc(nmel, sizeof(int));
    return 1;
}

int ReAllocate_grid(int gridno, int nmel, int nmnp)
{
    grid[gridno].nmel = nmel;
    grid[gridno].nmnp = nmnp;
    grid[gridno].gactive = 1;
    grid[gridno].xord = (double *) realloc(grid[gridno].xord, nmnp * sizeof(double));
    grid[gridno].yord = (double *) realloc(grid[gridno].yord, nmnp * sizeof(double));
    grid[gridno].depth = (double *) realloc(grid[gridno].depth, nmnp * sizeof(double));
    grid[gridno].nodecon = (conlist *) realloc(grid[gridno].nodecon, nmnp * sizeof(conlist));
    grid[gridno].ellist = (int *) realloc(grid[gridno].ellist, nmel * sizeof(int));
    grid[gridno].nlist = (int *) realloc(grid[gridno].nlist, nmnp * sizeof(int));
    grid[gridno].ntype = (int *) realloc(grid[gridno].ntype, nmnp * sizeof(int));
    grid[gridno].icon = (Element *) realloc(grid[gridno].icon, nmel * sizeof(Element));
    return 1;
}

int add_node(int gridno, double x, double y, double d)
{
    int nmnp = grid[gridno].nmnp + 1;

    grid[gridno].nmnp = nmnp;
    grid[gridno].xord = (double *) realloc(grid[gridno].xord, nmnp * sizeof(double));
    grid[gridno].yord = (double *) realloc(grid[gridno].yord, nmnp * sizeof(double));
    grid[gridno].depth = (double *) realloc(grid[gridno].depth, nmnp * sizeof(double));
    grid[gridno].nlist = (int *) realloc(grid[gridno].nlist, nmnp * sizeof(int));
    grid[gridno].ntype = (int *) realloc(grid[gridno].ntype, nmnp * sizeof(int));
    grid[gridno].xord[nmnp - 1] = x;
    grid[gridno].yord[nmnp - 1] = y;
    grid[gridno].depth[nmnp - 1] = d;
    return grid[gridno].nmnp - 1;
}

int add_node_type(int gridno, int type, double x, double y, double d)
{
    int nmnp = grid[gridno].nmnp + 1;

    grid[gridno].nmnp = nmnp;
    grid[gridno].xord = (double *) realloc(grid[gridno].xord, nmnp * sizeof(double));
    grid[gridno].yord = (double *) realloc(grid[gridno].yord, nmnp * sizeof(double));
    grid[gridno].depth = (double *) realloc(grid[gridno].depth, nmnp * sizeof(double));
    grid[gridno].nlist = (int *) realloc(grid[gridno].nlist, nmnp * sizeof(int));
    grid[gridno].ntype = (int *) realloc(grid[gridno].ntype, nmnp * sizeof(int));
    grid[gridno].xord[nmnp - 1] = x;
    grid[gridno].yord[nmnp - 1] = y;
    grid[gridno].depth[nmnp - 1] = d;
    grid[gridno].ntype[nmnp - 1] = type;
    return grid[gridno].nmnp - 1;
}

int add_nodes(int gridno, int n, double *x, double *y, double *d)
{
    int i;
    int nmnp = grid[gridno].nmnp;

    grid[gridno].nmnp += n;
    if (nmnp == 0) {
	grid[gridno].xord = (double *) calloc((nmnp + n), sizeof(double));
	grid[gridno].yord = (double *) calloc((nmnp + n), sizeof(double));
	grid[gridno].depth = (double *) calloc((nmnp + n), sizeof(double));
	grid[gridno].nlist = (int *) calloc((nmnp + n), sizeof(int));
	grid[gridno].ntype = (int *) calloc((nmnp + n), sizeof(int));
    } else {
	grid[gridno].xord = (double *) realloc(grid[gridno].xord, (nmnp + n) * sizeof(double));
	grid[gridno].yord = (double *) realloc(grid[gridno].yord, (nmnp + n) * sizeof(double));
	grid[gridno].depth = (double *) realloc(grid[gridno].depth, (nmnp + n) * sizeof(double));
	grid[gridno].nlist = (int *) realloc(grid[gridno].nlist, (nmnp + n) * sizeof(int));
	grid[gridno].ntype = (int *) realloc(grid[gridno].ntype, (nmnp + n) * sizeof(int));
    }
    for (i = nmnp; i < nmnp + n; i++) {
	grid[gridno].xord[i] = x[i - nmnp];
	grid[gridno].yord[i] = y[i - nmnp];
	grid[gridno].depth[i] = d[i - nmnp];
    }
    return 1;
}

/*
 * Add a triangle to the grid
 */
int add_element(int gridno, int n0, int n1, int n2)
{
    int nmel = grid[gridno].nmel + 1;

    grid[gridno].nmel = nmel;
    grid[gridno].ellist = (int *) realloc(grid[gridno].ellist, nmel * sizeof(int));
    grid[gridno].icon = (Element *) realloc(grid[gridno].icon, nmel * sizeof(Element));
    grid[gridno].icon[nmel - 1].type = 3;
    grid[gridno].icon[nmel - 1].nn = 3;
    grid[gridno].icon[nmel - 1].ngeom = 3;
    grid[gridno].ellist[nmel - 1] = 0;
    grid[gridno].icon[nmel - 1].nl[0] = n0;
    grid[gridno].icon[nmel - 1].nl[1] = n1;
    grid[gridno].icon[nmel - 1].nl[2] = n2;
    return nmel - 1;
}

/*
 * Add a quadrangle to the grid
 */
int add_quad_element(int gridno, int n0, int n1, int n2, int n3)
{
    int nmel = grid[gridno].nmel + 1;

    grid[gridno].nmel = nmel;
    grid[gridno].ellist = (int *) realloc(grid[gridno].ellist, nmel * sizeof(int));
    grid[gridno].icon = (Element *) realloc(grid[gridno].icon, nmel * sizeof(Element));
    grid[gridno].icon[nmel - 1].type = 4;
    grid[gridno].icon[nmel - 1].nn = 4;
    grid[gridno].icon[nmel - 1].ngeom = 4;
    grid[gridno].ellist[nmel - 1] = 0;
    grid[gridno].icon[nmel - 1].nl[0] = n0;
    grid[gridno].icon[nmel - 1].nl[1] = n1;
    grid[gridno].icon[nmel - 1].nl[2] = n2;
    grid[gridno].icon[nmel - 1].nl[3] = n3;
    return 1;
}

int add_element2(int gridno, int type, int n0, int n1, int n2, int n3, int n4, int n5)
{
    int nmel = grid[gridno].nmel + 1;

    grid[gridno].nmel = nmel;
    grid[gridno].ellist = (int *) realloc(grid[gridno].ellist, nmel * sizeof(int));
    grid[gridno].icon = (Element *) realloc(grid[gridno].icon, nmel * sizeof(Element));
    grid[gridno].icon[nmel - 1].type = type;
    grid[gridno].ellist[nmel - 1] = 0;
    switch (type) {
    case 3:
	grid[gridno].icon[nmel - 1].nn = 3;
	grid[gridno].icon[nmel - 1].ngeom = 3;
	grid[gridno].icon[nmel - 1].nl[0] = n0;
	grid[gridno].icon[nmel - 1].nl[1] = n1;
	grid[gridno].icon[nmel - 1].nl[2] = n2;
	break;
    case 4:
	grid[gridno].icon[nmel - 1].nn = 4;
	grid[gridno].icon[nmel - 1].ngeom = 4;
	grid[gridno].icon[nmel - 1].nl[0] = n0;
	grid[gridno].icon[nmel - 1].nl[1] = n1;
	grid[gridno].icon[nmel - 1].nl[2] = n2;
	grid[gridno].icon[nmel - 1].nl[3] = n3;
	break;
    case 6:
	grid[gridno].icon[nmel - 1].nn = 6;
	grid[gridno].icon[nmel - 1].ngeom = 3;
	grid[gridno].icon[nmel - 1].nl[0] = n0;
	grid[gridno].icon[nmel - 1].nl[1] = n1;
	grid[gridno].icon[nmel - 1].nl[2] = n2;
	grid[gridno].icon[nmel - 1].nl[3] = n3;
	grid[gridno].icon[nmel - 1].nl[4] = n4;
	grid[gridno].icon[nmel - 1].nl[5] = n5;
	break;
    }
    return nmel - 1;
}

int add_elements(int gridno, int n, int *n0, int *n1, int *n2)
{
    int i;
    int nmel = grid[gridno].nmel;
    grid[gridno].nmel += n;
    grid[gridno].ellist = (int *) realloc(grid[gridno].ellist, (nmel + n) * sizeof(int));
    grid[gridno].icon = (Element *) realloc(grid[gridno].icon, (nmel + n) * sizeof(Element));
    for (i = nmel; i < nmel + n; i++) {
	grid[gridno].icon[i].type = 3;
	grid[gridno].icon[i].nn = 3;
	grid[gridno].icon[i].ngeom = 3;
	grid[gridno].icon[i].nl[0] = n0[i - nmel];
	grid[gridno].icon[i].nl[1] = n1[i - nmel];
	grid[gridno].icon[i].nl[2] = n2[i - nmel];
    }
    return 1;
}

void Free_grid(int gridno)
{
    int i;
    grid[gridno].nmel = 0;
    grid[gridno].nmnp = 0;
    grid[gridno].gactive = 0;
    grid[gridno].bactive = 0;
    strcpy(grid[gridno].projection, "UNKNOWN");
    grid[gridno].proj = -1;
    grid[gridno].invproj = 0;
    grid[gridno].lat = 0.0;
    grid[gridno].lon = 0.0;
    for (i = 0; i < grid[gridno].nbounds; i++) {
	Free_boundary(grid[gridno].boundaries[i]);
    }
    if (grid[gridno].xord != NULL) {
	free(grid[gridno].xord);
    }
    if (grid[gridno].yord != NULL) {
	free(grid[gridno].yord);
    }
    if (grid[gridno].depth != NULL) {
	free(grid[gridno].depth);
    }
    if (grid[gridno].icon != NULL) {
	free(grid[gridno].icon);
    }
    if (grid[gridno].ellist != NULL) {
	free(grid[gridno].ellist);
    }
    if (grid[gridno].nlist != NULL) {
	free(grid[gridno].nlist);
    }
    if (grid[gridno].ntype != NULL) {
	free(grid[gridno].ntype);
    }
}

void Free_grid_only(int gridno)
{
    int i;
    grid[gridno].nmel = 0;
    grid[gridno].nmnp = 0;
    grid[gridno].gactive = 0;
    if (grid[gridno].xord != NULL) {
	free(grid[gridno].xord);
    }
    if (grid[gridno].yord != NULL) {
	free(grid[gridno].yord);
    }
    if (grid[gridno].depth != NULL) {
	free(grid[gridno].depth);
    }
    if (grid[gridno].icon != NULL) {
	free(grid[gridno].icon);
    }
    if (grid[gridno].ellist != NULL) {
	free(grid[gridno].ellist);
    }
    if (grid[gridno].nlist != NULL) {
	free(grid[gridno].nlist);
    }
    if (grid[gridno].ntype != NULL) {
	free(grid[gridno].ntype);
    }
}

void compact_grid(int gridno)
{
    int nn, i, j, itmp, n1 = 0, nodetmp = 0, eltmp = 0;
    double *xtmp, *ytmp, *dtmp;
    int *ntmp;

    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = -1;
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!grid[gridno].ellist[i]) {
	    nn = grid[gridno].icon[i].nn;
	    if (nn >= 0 && nn <= 8) {
		for (j = 0; j < nn; j++) {
		    itmp = grid[gridno].icon[i].nl[j];
		    if (grid[gridno].nlist[itmp] == -1) {
			grid[gridno].nlist[itmp] = nodetmp;
			nodetmp++;
		    }
		    grid[gridno].icon[eltmp].nl[j] = grid[gridno].nlist[itmp];
		}
		grid[gridno].icon[eltmp].type = grid[gridno].icon[i].type;
		grid[gridno].icon[eltmp].nn = grid[gridno].icon[i].nn;
		grid[gridno].icon[eltmp].ngeom = grid[gridno].icon[i].ngeom;
		eltmp++;
	    } else {
		printf("compact_grid(%d): Error in nn = %d, el = %d\n", gridno, nn, i);
	    }
	} else {		/* element deleted */
/*
	    printf("Deleting element %d type = %d\n", i, grid[gridno].icon[i].type);
*/
	}
    }
    xtmp = (double *) calloc(nodetmp, sizeof(double));
    ytmp = (double *) calloc(nodetmp, sizeof(double));
    dtmp = (double *) calloc(nodetmp, sizeof(double));
    ntmp = (int *) calloc(nodetmp, sizeof(int));
    for (i = 0; i < grid[gridno].nmnp; i++) {
	itmp = grid[gridno].nlist[i];
	if (itmp >= 0) {
	    xtmp[itmp] = grid[gridno].xord[i];
	    ytmp[itmp] = grid[gridno].yord[i];
	    dtmp[itmp] = grid[gridno].depth[i];
	    ntmp[itmp] = grid[gridno].ntype[i];
	}
    }
    for (i = 0; i < nodetmp; i++) {
	grid[gridno].xord[i] = xtmp[i];
	grid[gridno].yord[i] = ytmp[i];
	grid[gridno].depth[i] = dtmp[i];
	grid[gridno].ntype[i] = ntmp[i];
    }
    free(xtmp);
    free(ytmp);
    free(dtmp);
    free(ntmp);
/*
    printf("%d %d %d %d\n", grid[gridno].nmel, grid[gridno].nmnp, eltmp, nodetmp);
*/
    ReAllocate_grid(gridno, eltmp, nodetmp);
}

int delete_node(int gridno, int n)
{
    int i, j;
/*
    int nmnp = grid[gridno].nmnp - 1;
    for (i = n; i < nmnp; i++) {
	grid[gridno].xord[i] = grid[gridno].xord[i + 1];
	grid[gridno].yord[i] = grid[gridno].yord[i + 1];
	grid[gridno].depth[i] = grid[gridno].depth[i + 1];
	grid[gridno].nlist[i] = grid[gridno].nlist[i + 1];
	grid[gridno].ntype[i] = grid[gridno].ntype[i + 1];
    }
    grid[gridno].nmnp = nmnp;
    grid[gridno].xord = (double *) realloc(grid[gridno].xord, nmnp * sizeof(double));
    grid[gridno].yord = (double *) realloc(grid[gridno].yord, nmnp * sizeof(double));
    grid[gridno].depth = (double *) realloc(grid[gridno].depth, nmnp * sizeof(double));
    grid[gridno].nlist = (int *) realloc(grid[gridno].nlist, nmnp * sizeof(int));
    grid[gridno].ntype = (int *) realloc(grid[gridno].ntype, nmnp * sizeof(int));
*/
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
	for (j = 0; j < numnodes(gridno, i); j++) {
	    if (n == grid[gridno].icon[i].nl[j]) {
		grid[gridno].ellist[i] = 1;
		break;
	    }
	}
    }
    compact_grid(gridno);
    return 1;
}

/*
 * Pointer based versions of the above routines
 */ 
int AllocateGrid(Grid * g, int nmel, int nmnp)
{
    int i;
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->gactive = 1;
    g->bactive = 0;
    g->nbounds = 0;
    for (i = 0; i < MAXGRIDBOUNDS; i++) {
	g->boundaries[i] = -1;
    }
    g->xord = (double *) calloc(nmnp, sizeof(double));
    g->yord = (double *) calloc(nmnp, sizeof(double));
    g->depth = (double *) calloc(nmnp, sizeof(double));
    g->nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    g->ellist = (int *) calloc(nmel, sizeof(int));
    g->nlist = (int *) calloc(nmnp, sizeof(int));
    g->ntype = (int *) calloc(nmnp, sizeof(int));
    g->icon = (Element *) calloc(nmel, sizeof(Element));
    return 1;
}

int AllocateGridOnly(Grid * g, int nmel, int nmnp)
{
    int i;
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->gactive = 1;
    if (g->xord == NULL) {
    }
    g->xord = (double *) calloc(nmnp, sizeof(double));
    g->yord = (double *) calloc(nmnp, sizeof(double));
    g->depth = (double *) calloc(nmnp, sizeof(double));
    g->nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    g->nlist = (int *) calloc(nmnp, sizeof(int));
    g->ntype = (int *) calloc(nmnp, sizeof(int));
    g->icon = (Element *) calloc(nmel, sizeof(Element));
    g->ellist = (int *) calloc(nmel, sizeof(int));
    return 1;
}

int AllocateGridNodes(Grid * g, int nmnp)
{
    int i;
    g->nmnp = nmnp;
    g->xord = (double *) calloc(nmnp, sizeof(double));
    g->yord = (double *) calloc(nmnp, sizeof(double));
    g->depth = (double *) calloc(nmnp, sizeof(double));
    g->nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    g->nlist = (int *) calloc(nmnp, sizeof(int));
    g->ntype = (int *) calloc(nmnp, sizeof(int));
    return 1;
}

int AllocateGridElements(Grid * g, int nmel)
{
    int i;
    g->nmel = nmel;
    g->icon = (Element *) calloc(nmel, sizeof(Element));
    g->ellist = (int *) calloc(nmel, sizeof(int));
    return 1;
}

void FreeGrid(Grid * g)
{
    int i;
    g->nmel = 0;
    g->nmnp = 0;
    g->gactive = 0;
    g->bactive = 0;
    strcpy(g->projection, "UNKNOWN");
    g->proj = -1;
    g->invproj = 0;
    g->lat = 0.0;
    g->lon = 0.0;
    for (i = 0; i < g->nbounds; i++) {
	Free_boundary(g->boundaries[i]);
    }
    if (g->xord != NULL) {
	free(g->xord);
    }
    if (g->yord != NULL) {
	free(g->yord);
    }
    if (g->depth != NULL) {
	free(g->depth);
    }
    if (g->icon != NULL) {
	free(g->icon);
    }
    if (g->ellist != NULL) {
	free(g->ellist);
    }
    if (g->nlist != NULL) {
	free(g->nlist);
    }
    if (g->ntype != NULL) {
	free(g->ntype);
    }
}

void FreeGridOnly(Grid * g)
{
    int i;
    g->nmel = 0;
    g->nmnp = 0;
    g->gactive = 0;
    if (g->xord != NULL) {
	free(g->xord);
    }
    if (g->yord != NULL) {
	free(g->yord);
    }
    if (g->depth != NULL) {
	free(g->depth);
    }
    if (g->icon != NULL) {
	free(g->icon);
    }
    if (g->ellist != NULL) {
	free(g->ellist);
    }
    if (g->nlist != NULL) {
	free(g->nlist);
    }
    if (g->ntype != NULL) {
	free(g->ntype);
    }
}

int ReallocateGrid(Grid * g, int nmel, int nmnp)
{
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->gactive = 1;
    g->xord = (double *) realloc(g->xord, nmnp * sizeof(double));
    g->yord = (double *) realloc(g->yord, nmnp * sizeof(double));
    g->depth = (double *) realloc(g->depth, nmnp * sizeof(double));
    g->nodecon = (conlist *) realloc(g->nodecon, nmnp * sizeof(conlist));
    g->ellist = (int *) realloc(g->ellist, nmel * sizeof(int));
    g->nlist = (int *) realloc(g->nlist, nmnp * sizeof(int));
    g->ntype = (int *) realloc(g->ntype, nmnp * sizeof(int));
    g->icon = (Element *) realloc(g->icon, nmel * sizeof(Element));
    return 1;
}

int ReallocateGridElements(Grid * g, int nmel)
{
    g->nmel = nmel;
    g->gactive = 1;
    g->ellist = (int *) realloc(g->ellist, nmel * sizeof(int));
    g->icon = (Element *) realloc(g->icon, nmel * sizeof(Element));
    return 1;
}

int ReallocateGridNodes(Grid * g, int nmnp)
{
    g->nmnp = nmnp;
    g->gactive = 1;
    g->xord = (double *) realloc(g->xord, nmnp * sizeof(double));
    g->yord = (double *) realloc(g->yord, nmnp * sizeof(double));
    g->depth = (double *) realloc(g->depth, nmnp * sizeof(double));
    g->nodecon = (conlist *) realloc(g->nodecon, nmnp * sizeof(conlist));
    g->nlist = (int *) realloc(g->nlist, nmnp * sizeof(int));
    g->ntype = (int *) realloc(g->ntype, nmnp * sizeof(int));
    return 1;
}

/*
 * delete elements as indicated by ellist != 0
 */
void DeleteElements(Grid * g)
{
    int nn, i, j, itmp, n1 = 0, nodetmp = 0, eltmp = 0;
    double *xtmp, *ytmp, *dtmp;
    int *ntmp;

    for (i = 0; i < g->nmel; i++) {
	if (!g->ellist[i]) {
	    nn = g->icon[i].nn;
	    for (j = 0; j < nn; j++) {
		g->icon[eltmp].nl[j] = g->icon[i].nl[j];
	    }
	    g->icon[eltmp].type = g->icon[i].type;
	    g->icon[eltmp].nn = g->icon[i].nn;
	    g->icon[eltmp].ngeom = g->icon[i].ngeom;
	    eltmp++;
	} else {		/* element deleted */
	}
    }
    if (eltmp < g->nmel) {
	ReallocateGridElements(g, eltmp);
    }
}

int InsideGrid(Grid * g, double x, double y)
{
    int i, ib, itmp;
    for (i = 0; i < g->nbounds; i++) {
	ib = g->boundaries[i];
	if (boundary[ib].bactive) {
	    itmp = inbound(x, y, boundary[ib].boundx,
			   boundary[ib].boundy,
			   boundary[ib].nbpts);
	    if (boundary[ib].boundtype == 0) {
		if (itmp == 0) {
		    return 0;
		}
	    } else if (itmp == 1) {
		return 0;
	    }
	}
    }
    return 1;
}

void CleanGrid(Grid * g)
{
    int i;
    double x, y;
    for (i = 0; i < g->nmel; i++) {
	GetCenter(g, i, &x, &y);
	if (!InsideGrid(g, x, y)) {
	    g->ellist[i] = 1;
	} else {
	    g->ellist[i] = 0;
	}
    }
}
