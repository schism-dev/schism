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
 * gridtio.c - read/write a finite element grid
 *
*/

#ifndef lint
static char RCSid[] = "$Id: gridtio.c,v 1.3 2003/10/02 04:59:00 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

void Free_one_grid(Grid * g);
void Free_gridt(Gridt * gt);
void set_limits_grid(Grid * g);

int readgridt(int gridno, char *fname)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type;
    double t;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readgrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    Free_gridt(&gridt[gridno]);
    strcpy(gridt[gridno].fname, fname);
    fgets(buf, 255, fp);	/* first line is ignored */
    fgets(buf, 255, fp);
    sscanf(buf, "%d", &gridt[gridno].nsteps);
    gridt[gridno].grids = (Grid *) calloc(gridt[gridno].nsteps, sizeof(Grid));
    if (gridt[gridno].grids != NULL) {
	for (i = 0; i < gridt[gridno].nsteps; i++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%lf", &t);
	    if (read_one_grid(fp, &gridt[gridno].grids[i])) {
		gridt[gridno].grids[i].time = t;
	    } else {
		fclose(fp);
		return 0;
	    }
	}
	fclose(fp);
	gridt[gridno].active = ON;
	return 1;
    } else {
	return 0;
    }
}

int read_one_grid(FILE * fp, Grid * g)
{
    char buf[256];
    int i, j, itmp, type;
    Grid gtmp;
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &gtmp.nmel, &gtmp.nmnp);
    if (!Allocate_one_grid(&gtmp, gtmp.nmel, gtmp.nmnp)) {
	errwin("Can't allocate memory for grid");
	return 0;
    }
    for (i = 0; i < gtmp.nmnp; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("Error reading table of nodes");
	    Free_one_grid(&gtmp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %lf %lf %lf", &itmp, &gtmp.xord[i], &gtmp.yord[i], &gtmp.depth[i]);
    }
    getlimits(gtmp.xord, gtmp.nmnp, &gtmp.xmin, &gtmp.xmax);
    getlimits(gtmp.yord, gtmp.nmnp, &gtmp.ymin, &gtmp.ymax);
    getlimits(gtmp.depth, gtmp.nmnp, &gtmp.dmin, &gtmp.dmax);
    for (i = 0; i < gtmp.nmel; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("Error reading table of elements");
	    Free_one_grid(&gtmp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %d", &itmp, &type);
	gtmp.icon[i].type = type;
	switch (type) {
	case 0:		/* temporary TODO */
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2], &gtmp.icon[i].nl[3]);
	    for (j = 0; j < 4; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].nn = gtmp.icon[i].ngeom = 4;
	    break;
	case 3:
	    sscanf(buf, "%d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2]);
	    for (j = 0; j < 3; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].nn = gtmp.icon[i].ngeom = 3;
	    break;
	case 4:
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2], &gtmp.icon[i].nl[3]);
	    for (j = 0; j < 4; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].nn = gtmp.icon[i].ngeom = 4;
	    break;
	case 6:
	    sscanf(buf, "%d %d %d %d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2], &gtmp.icon[i].nl[3], &gtmp.icon[i].nl[4], &gtmp.icon[i].nl[5]);
	    for (j = 0; j < 6; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].ngeom = 3;
	    gtmp.icon[i].nn = 6;
	    break;
	case 5:		/* connecting 2d and 1d elements */
	    sscanf(buf, "%d %d %d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2], &gtmp.icon[i].nl[3], &gtmp.icon[i].nl[4]);
	    for (j = 0; j < 5; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].ngeom = 5;
	    gtmp.icon[i].nn = 5;
	    break;
	case 8:
	    sscanf(buf, "%d %d %d %d %d %d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2], &gtmp.icon[i].nl[3], &gtmp.icon[i].nl[4], &gtmp.icon[i].nl[5], &gtmp.icon[i].nl[6], &gtmp.icon[i].nl[7]);
	    for (j = 0; j < 8; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].ngeom = 4;
	    gtmp.icon[i].nn = 8;
	    break;
	case 1:
	    sscanf(buf, "%d %d %d %d %d", &itmp, &itmp, &gtmp.icon[i].nl[0], &gtmp.icon[i].nl[1], &gtmp.icon[i].nl[2]);
	    for (j = 0; j < 3; j++)
		gtmp.icon[i].nl[j]--;
	    gtmp.icon[i].ngeom = 2;
	    gtmp.icon[i].nn = 3;
	    break;
	}
    }
    set_limits_grid(&gtmp);
    *g = gtmp;
    return 1;
}

int Allocate_one_grid(Grid * g, int nmel, int nmnp)
{
    int i;
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->active = ON;
    g->nbounds = 0;
    g->boundaries = NULL;
    g->xord = (double *) calloc(nmnp, sizeof(double));
    g->yord = (double *) calloc(nmnp, sizeof(double));
    g->depth = (double *) calloc(nmnp, sizeof(double));
    g->icon = (Element *) calloc(nmel, sizeof(Element));
    return 1;
}

void Free_gridt(Gridt * gt)
{
    int i;
    for (i = 0; i < gt->nsteps; i++) {
	Free_one_grid(&(gt->grids[i]));
    }
}

void Free_one_grid(Grid * g)
{
    int i;
    if (g->boundaries != NULL) {
	free(g->boundaries);
    }
    if (g->xord != NULL) {
	free(g->xord);
    }
    if (g->yord != NULL) {
	free(g->yord);
    }
    if (g->icon != NULL) {
	free(g->icon);
    }
    g->nmel = 0;
    g->nmnp = 0;
    g->active = OFF;
}

void set_limits_grid(Grid * g)
{
    int i;
    double xmin, xmax, ymin, ymax, dmin, dmax;

    if (g->active == ON) {
	xmin = xmax = g->xord[0];
	ymin = ymax = g->yord[0];
	dmin = dmax = g->depth[0];
	for (i = 1; i < g->nmnp; i++) {
	    if (g->xord[i] < xmin) {
		xmin = g->xord[i];
	    }
	    if (g->xord[i] > xmax) {
		xmax = g->xord[i];
	    }
	    if (g->yord[i] < ymin) {
		ymin = g->yord[i];
	    }
	    if (g->yord[i] > ymax) {
		ymax = g->yord[i];
	    }
	    if (g->depth[i] < dmin) {
		dmin = g->depth[i];
	    }
	    if (g->depth[i] > dmax) {
		dmax = g->depth[i];
	    }
	}
	g->xmin = xmin;
	g->xmax = xmax;
	g->ymin = ymin;
	g->ymax = ymax;
	g->dmin = dmin;
	g->dmax = dmax;
    }
}

void drawgridt(int gno, int gridno, double red, int step, double t)
{
    int i, j, n0, n1, n2, n3, redflag = (red == 1.0);
    double xg, yg, x, y;
    Grid gtmp;
    gtmp = gridt[gridno].grids[step];
    for (i = 0; i < gtmp.nmel; i++) {
	switch (gtmp.icon[i].type) {
	case 0:
	    n0 = gtmp.icon[i].nl[0];
	    n1 = gtmp.icon[i].nl[1];
	    n2 = gtmp.icon[i].nl[2];
	    n3 = gtmp.icon[i].nl[3];
	    my_move2(gtmp.xord[n0], gtmp.yord[n0]);
	    my_draw2(gtmp.xord[n1], gtmp.yord[n1]);
	    my_draw2(gtmp.xord[n2], gtmp.yord[n2]);
	    my_draw2(gtmp.xord[n3], gtmp.yord[n3]);
	    break;
	case 1:
	    n0 = gtmp.icon[i].nl[0];
	    n1 = gtmp.icon[i].nl[1];
	    n2 = gtmp.icon[i].nl[2];
	    my_move2(gtmp.xord[n0], gtmp.yord[n0]);
	    my_draw2(gtmp.xord[n2], gtmp.yord[n2]);
	    my_draw2(gtmp.xord[n1], gtmp.yord[n1]);
	    break;
	case 8:
	    n0 = gtmp.icon[i].nl[0];
	    n1 = gtmp.icon[i].nl[1];
	    n2 = gtmp.icon[i].nl[2];
	    n3 = gtmp.icon[i].nl[3];
	    my_move2(gtmp.xord[n0], gtmp.yord[n0]);
	    my_draw2(gtmp.xord[n1], gtmp.yord[n1]);
	    my_draw2(gtmp.xord[n2], gtmp.yord[n2]);
	    my_draw2(gtmp.xord[n3], gtmp.yord[n3]);
	    my_draw2(gtmp.xord[n0], gtmp.yord[n0]);
	    break;
	case 5:
	    n0 = gtmp.icon[i].nl[0];
	    n1 = gtmp.icon[i].nl[1];
	    n2 = gtmp.icon[i].nl[2];
	    n3 = gtmp.icon[i].nl[3];
	    my_move2(gtmp.xord[n0], gtmp.yord[n0]);
	    my_draw2(gtmp.xord[n1], gtmp.yord[n1]);
	    my_draw2(gtmp.xord[n2], gtmp.yord[n2]);
	    my_move2(gtmp.xord[n1], gtmp.yord[n1]);
	    my_draw2(gtmp.xord[n3], gtmp.yord[n3]);
	    break;
	case 3:
	case 6:
	    n0 = gtmp.icon[i].nl[0];
	    n1 = gtmp.icon[i].nl[1];
	    n2 = gtmp.icon[i].nl[2];
	    my_move2(gtmp.xord[n0], gtmp.yord[n0]);
	    my_draw2(gtmp.xord[n1], gtmp.yord[n1]);
	    my_draw2(gtmp.xord[n2], gtmp.yord[n2]);
	    my_draw2(gtmp.xord[n0], gtmp.yord[n0]);
	    break;
	}
    }
}

void draw_gridtisol(int gno, int gridno, int step, double t)
{
    extern int mapisolconc[];
    Grid gtmp;
    gtmp = gridt[gridno].grids[step];
    do_grid_isol(&gtmp, 0, gtmp.depth, g[gno].gridt[gridno].ip, mapisolconc, g[gno].gridt[gridno].ip.nisol);
}
