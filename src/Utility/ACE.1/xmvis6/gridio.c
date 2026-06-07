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
 * gridio.c - read/write a finite element grid
 *
 * Contents:
 *
 * int readgrid() ..... read a linear grid
 * int writegrid() ..... write a linear grid
 *
 * both fuctions return 1 if successful, 0 if not
*/

#ifndef lint
static char RCSid[] = "$Id: gridio.c,v 1.11 2006/08/18 15:40:36 pturner Exp $";
#endif

extern int winsetwidth, winsetheight;

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"
#include "elio.h"

int ReadElcircGrid(char *fname, Grid * g)
{
    int i, err;
    ElcircHeader h;
/* Get the header */
    if (err = ElioGetHeader(fname, &h)) {
	fprintf(stderr, "ReadElcircNew(): Error in ElioGetHeader(): Error # %d\n", err);
	return 0;
    }
    /*ElioPrintHeader(&h); */
    AllocateGrid(g, h.ne, h.np);
    for (i = 0; i < h.np; i++) {
	g->depth[i] = h.d[i];
	if (h.v == 4) {
	    g->xord[i] = h.x[i];
	    g->yord[i] = h.y[i];
	} else {
	    g->xord[i] = h.x[i];
	    g->yord[i] = h.y[i];
	}
    }
    for (i = 0; i < h.ne; i++) {
	g->icon[i].nl[0] = h.icon[0][i];
	g->icon[i].nl[1] = h.icon[1][i];
	g->icon[i].nl[2] = h.icon[2][i];
	if (h.etype[i] == 4) {
	    g->icon[i].nl[3] = h.icon[3][i];
	}
	g->icon[i].type = g->icon[i].nn = g->icon[i].ngeom = h.etype[i];
    }
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);

    ElioFreeHeader(&h);
}

int GetElcircGrid(ElcircHeader * h, Grid * g)
{
    int i, err;
    AllocateGrid(g, h->ne, h->np);
    for (i = 0; i < h->np; i++) {
	g->depth[i] = h->d[i];
	if (h->v == 4) {
	    g->xord[i] = h->x[i];
	    g->yord[i] = h->y[i];
	} else {
	    g->xord[i] = h->x[i];
	    g->yord[i] = h->y[i];
	}
    }
    for (i = 0; i < h->ne; i++) {
	g->icon[i].nl[0] = h->icon[0][i];
	g->icon[i].nl[1] = h->icon[1][i];
	g->icon[i].nl[2] = h->icon[2][i];
	if (h->etype[i] == 4) {
	    g->icon[i].nl[3] = h->icon[3][i];
	}
	g->icon[i].type = g->icon[i].nn = g->icon[i].ngeom = h->etype[i];
    }
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);
    return 0;
}

int readgrid(int gridno, char *fname)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type;

    sprintf(statusstr, "Start readgrid: %d %s", gridno, fname);
    writelogfile(statusstr);
    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readgrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    if (ElioGetFileType(fp) > 0) {
	Free_grid(gridno);
	fclose(fp);
	ReadElcircGrid(fname, &grid[gridno]);
    } else {
	Free_grid(gridno);
	strcpy(grid[gridno].fname, fname);
	fgets(buf, 255, fp);	/* first line is ignored */
	fgets(buf, 255, fp);
	convertchar(buf);
	sscanf(buf, "%d %d", &grid[gridno].nmel, &grid[gridno].nmnp);
	if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	    errwin("Can't allocate memory for grid");
	    Free_grid(gridno);
	    return 0;
	}
	for (i = 0; i < grid[gridno].nmnp; i++) {
	    if (fgets(buf, 255, fp) == NULL) {
		errwin("Error reading table of nodes");
		Free_grid(gridno);
		fclose(fp);
		return 0;
	    }
	    convertchar(buf);
	    sscanf(buf, "%d %lf %lf %lf", &itmp, &grid[gridno].xord[i], &grid[gridno].yord[i], &grid[gridno].depth[i]);
	}
	dminmax(grid[gridno].nmnp, grid[gridno].xord, &grid[gridno].xmin, &grid[gridno].xmax);
	dminmax(grid[gridno].nmnp, grid[gridno].yord, &grid[gridno].ymin, &grid[gridno].ymax);
	dminmax(grid[gridno].nmnp, grid[gridno].depth, &grid[gridno].dmin, &grid[gridno].dmax);
	for (i = 0; i < grid[gridno].nmel; i++) {
	    if (fgets(buf, 255, fp) == NULL) {
		errwin("Error reading table of elements");
		Free_grid(gridno);
		fclose(fp);
		return 0;
	    }
	    convertchar(buf);
	    sscanf(buf, "%d %d", &itmp, &type);
	    grid[gridno].icon[i].type = type;
	    switch (type) {
	    case 3:
		sscanf(buf, "%d %d %d %d %d", &itmp, &itmp, &grid[gridno].icon[i].nl[0], &grid[gridno].icon[i].nl[1], &grid[gridno].icon[i].nl[2]);
		for (j = 0; j < 3; j++)
		    grid[gridno].icon[i].nl[j]--;
		grid[gridno].icon[i].nn = grid[gridno].icon[i].ngeom = 3;
		break;
	    case 4:
		sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp, &grid[gridno].icon[i].nl[0], &grid[gridno].icon[i].nl[1], &grid[gridno].icon[i].nl[2], &grid[gridno].icon[i].nl[3]);
		for (j = 0; j < 4; j++)
		    grid[gridno].icon[i].nl[j]--;
		grid[gridno].icon[i].nn = grid[gridno].icon[i].ngeom = 4;
		break;
	    }
	}
/* There might be some boundary information */
	if ((fgets(buf, 255, fp) != NULL) && strlen(buf) != 0) {
	    int nob, tnon, nop, non, nn, is;
	    int nlb, tnln, nlp, nln;
	    is = sscanf(buf, "%d", &nob);	/* number of open boundaries */

	    if (is > 0 && (fgets(buf, 255, fp) != NULL)) {

		is = sscanf(buf, "%d", &tnon);	/* number of open nodes (ignored) */
		for (i = 0; i < nob; i++) {
		    if (fgets(buf, 255, fp) != NULL) {
			is = sscanf(buf, "%d", &nop);
			for (j = 0; j < nop; j++) {
			    fgets(buf, 255, fp);
			    is = sscanf(buf, "%d", &nn);
			    //printf("%d\n", nn);
			}
		    }
		}
	    }
	}
	fclose(fp);
    }
    setlimits_grid(gridno);
    sprintf(statusstr, "End readgrid: %d %s", gridno, fname);
    writelogfile(statusstr);
    return 1;
}

int writegrid(int gridno, char *fname)
{
    FILE *fp;
    int i;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin(stderr, "In writegrid, unable to open file %s\n", fname);
	fclose(fp);
	return 0;
    }
    fprintf(fp, fname);
    fprintf(fp, "\n%d %d\n", grid[gridno].nmel, grid[gridno].nmnp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	fprintf(fp, "%d %14.6lf %14.6lf %12.3lf\n", i + 1, grid[gridno].xord[i], grid[gridno].yord[i], grid[gridno].depth[i]);
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	fprintf(fp, "%d 3 %d %d %d\n", i + 1, grid[gridno].icon[i].nl[0] + 1, grid[gridno].icon[i].nl[1] + 1, grid[gridno].icon[i].nl[2] + 1);
    }
    fclose(fp);
    return 1;
}

void autoscalegrid(int gno, int gridno)
{
    g[gno].w.xg1 = grid[gridno].xmin;
    g[gno].w.xg2 = grid[gridno].xmax;
    g[gno].w.yg1 = grid[gridno].ymin;
    g[gno].w.yg2 = grid[gridno].ymax;
}

int readboundary2(int bno, char *fname)
{
    int i, j, nb, npts, ib;
    FILE *fp;
    char s[255];

    if ((fp = fopen(fname, "r")) == NULL) {
	errwin("readboundary: can't open file for reading");
	return 0;
    }
    Free_boundaries2(bno);
    fgets(s, 255, fp);
    fgets(s, 255, fp);
    sscanf(s, "%d", &bounds[bno].nbounds);
    Allocate_boundaries2(bno, bounds[bno].nbounds);
    bounds[bno].active = ON;
    for (i = 0; i < bounds[bno].nbounds; i++) {
	fgets(s, 255, fp);
	sscanf(s, "%d", &npts);
	Allocate_boundary2(bno, i, npts, (i != 0));
	for (j = 0; j < npts; j++) {
	    fgets(s, 255, fp);
	    sscanf(s, "%lf %lf", &bounds[bno].boundaries[i].boundx[j], &bounds[bno].boundaries[i].boundy[j]);
	}
    }
    fclose(fp);
    return 1;
}

int readbinboundary(int bno, char *fname, int val)
{
    int i, j, n, nb, ib, npts, itmp, magic;
    FILE *fp;
    char s[255];
    float *xbtmp, *ybtmp;
    double *xb, *yb;
    if ((fp = fopen(fname, "r")) == NULL) {
	return 0;
    }
    fread(&magic, 1, sizeof(int), fp);
    if (magic != 32) {
	return 0;
    }
    Free_boundaries2(bno);
    fread(&nb, 1, sizeof(int), fp);
    bounds[bno].nbounds = nb;
    Allocate_boundaries2(bno, bounds[bno].nbounds);
    bounds[bno].active = ON;
    for (i = 0; i < nb; i++) {
	fread(&n, 1, sizeof(int), fp);
	Allocate_boundary2(bno, i, n, (i != 0));
	xbtmp = (float *) calloc(n, sizeof(float));
	ybtmp = (float *) calloc(n, sizeof(float));
	fread(xbtmp, sizeof(float), n, fp);
	fread(ybtmp, sizeof(float), n, fp);
	xb = bounds[bno].boundaries[i].boundx;
	yb = bounds[bno].boundaries[i].boundy;
	for (j = 0; j < n; j++) {
	    xb[j] = xbtmp[j];
	    yb[j] = ybtmp[j];
	}
	free(xbtmp);
	free(ybtmp);
    }
    fclose(fp);
    return 1;
}

int ReadGrid(char *fname, Grid * g)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type;
    if ((fp = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "In ReadGrid, unable to open file %s\n", fname);
	return 0;
    }
    strcpy(g->fname, fname);
    fgets(buf, 255, fp);	/* first line is ignored */
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &g->nmel, &g->nmnp);
    if (!AllocateGrid(g, g->nmel, g->nmnp)) {
	fprintf(stderr, "Can't allocate memory for grid");
	return 0;
    }
    for (i = 0; i < g->nmnp; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    fprintf(stderr, "Error reading table of nodes");
	    fclose(fp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %lf %lf %lf", &itmp, &g->xord[i], &g->yord[i], &g->depth[i]);
    }
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);
    for (i = 0; i < g->nmel; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    fprintf(stderr, "Error reading table of elements");
	    fclose(fp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %d", &itmp, &type);
	g->icon[i].type = type;
	switch (type) {
	case 3:
	    sscanf(buf, "%d %d %d %d %d", &itmp, &itmp, &g->icon[i].nl[0], &g->icon[i].nl[1], &g->icon[i].nl[2]);
	    for (j = 0; j < 3; j++)
		g->icon[i].nl[j]--;
	    g->icon[i].nn = g->icon[i].ngeom = 3;
	    break;
	case 4:
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp, &g->icon[i].nl[0], &g->icon[i].nl[1], &g->icon[i].nl[2], &g->icon[i].nl[3]);
	    for (j = 0; j < 4; j++)
		g->icon[i].nl[j]--;
	    g->icon[i].nn = g->icon[i].ngeom = 4;
	    break;
	}
    }

/* There might be some boundary information */
    if ((fgets(buf, 255, fp) != NULL) && strlen(buf) != 0) {
	int nob, tnon, nop, non, nn, is;
	int nlb, tnln, nlp, nln;
	is = sscanf(buf, "%d", &nob);	/* number of open boundaries */

	if (is > 0 && (fgets(buf, 255, fp) != NULL)) {

	    is = sscanf(buf, "%d", &tnon);	/* number of open nodes (ignored) */
	    if (is > 0 && (fgets(buf, 255, fp) != NULL)) {
		for (i = 0; i < nob; i++) {
		    is = sscanf(buf, "%d", &nop);
		    for (j = 0; j < nop; j++) {
			is = sscanf(buf, "%d", &nn);
			//printf("%d\n", nn);
		    }
		}
	    }
	}
    }
    fclose(fp);
    return 1;
}

int WriteGrid(char *fname, Grid * g)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type;
    if ((fp = fopen(fname, "w")) == NULL) {
	fprintf(stderr, "In WriteGrid, unable to open file %s\n", fname);
	return 0;
    }
    fprintf(fp, "Grid\n");	/* first line is ignored */
    fprintf(fp, "%d %d\n", g->nmel, g->nmnp);
    for (i = 0; i < g->nmnp; i++) {
	fprintf(fp, "%d %lf %lf %lf\n", i + 1, g->xord[i], g->yord[i], g->depth[i]);
    }
    for (i = 0; i < g->nmel; i++) {
        type = g->icon[i].type;
	switch (type) {
	case 3:
	    fprintf(fp, "%d 3 %d %d %d\n", i + 1, g->icon[i].nl[0] + 1, g->icon[i].nl[1] + 1, g->icon[i].nl[2] + 1);
	    break;
	case 4:
	    fprintf(fp, "%d 4 %d %d %d %d\n", i + 1, g->icon[i].nl[0] + 1, g->icon[i].nl[1] + 1, g->icon[i].nl[2] + 1, g->icon[i].nl[3] + 1);
	    break;
	}
    }
    fclose(fp);
    return 1;
}
