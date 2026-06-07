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

 * gridio.c - read/write a finite element grid
 *
 * Contents:
 *
 * int readgrid(int gridno, char *fname);
 * int writegrid(int gridno, char *fname);
 * int writegridbin(int gridno, char *fname);
 * int readgridbin(int gridno, char *fname);
 * void set_grid_limits(int gridno);
 * Grid *NewGrid(int nmel, int nmnp);
 * void DeleteGrid(Grid *g);
 * int readpomgrid(int gridno, char *fname);
 * void check_err(const int stat, const int line, const char *file);
 * void CreateGridnetcdf(char *fname, int gridno);
 * void WriteGridnetcdf(char *fname, int gridno);
 * void ReadGridnetcdf(char *fname, int gridno);
 *
 */

#ifndef lint
static char RCSid[] = "$Id: gridio.c,v 1.5 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "globals.h"

/* gridio.c */
int readgrid(int gridno, char *fname);
int mergegrid(int gridno, char *fname);
int writegrid(int gridno, char *fname);
int writegridbin(int gridno, char *fname);
int readgridbin(int gridno, char *fname);
void set_grid_limits(int gridno);
Grid *NewGrid(int nmel, int nmnp);
void DeleteGrid(Grid *g);
int readpomgrid(int gridno, char *fname);
void check_err(const int stat, const int line, const char *file);
void CreateGridnetcdf(char *fname, int gridno);
void WriteGridnetcdf(char *fname, int gridno);
void ReadGridnetcdf(char *fname, int gridno);

/*
 * Read a grid in ASCII format, set defaults for isolines
 */
int readgrid(int gridno, char *fname)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type, nn, nnodes;
    extern int doadcirc;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readgrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    Free_grid(gridno);

    fgets(buf, 255, fp);
    strcpy(grid[gridno].alphid, buf);
    i = strlen(grid[gridno].alphid);
    grid[gridno].alphid[i - 1] = '\0';
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &grid[gridno].nmel, &grid[gridno].nmnp);
    if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	fclose(fp);
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
	    sscanf(buf, "%d %d %d %d %d", &itmp, &itmp,
		   &grid[gridno].icon[i].nl[0],
		   &grid[gridno].icon[i].nl[1],
		   &grid[gridno].icon[i].nl[2]);
	    for (j = 0; j < 3; j++) {
		grid[gridno].icon[i].nl[j]--;
		nn = grid[gridno].icon[i].nl[j];
		grid[gridno].ntype[nn] = CORNER;
	    }
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	    break;
	case 4:
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp,
		   &grid[gridno].icon[i].nl[0],
		   &grid[gridno].icon[i].nl[1],
		   &grid[gridno].icon[i].nl[2],
		   &grid[gridno].icon[i].nl[3]);
	    for (j = 0; j < 4; j++) {
		grid[gridno].icon[i].nl[j]--;
		nn = grid[gridno].icon[i].nl[j];
		grid[gridno].ntype[nn] = CORNER;
	    }
	    grid[gridno].icon[i].nn = 4;
	    grid[gridno].icon[i].ngeom = 4;
	    break;
	case 6:
	    sscanf(buf, "%d %d %d %d %d %d %d %d", &itmp, &itmp,
		   &grid[gridno].icon[i].nl[0],
		   &grid[gridno].icon[i].nl[1],
		   &grid[gridno].icon[i].nl[2],
		   &grid[gridno].icon[i].nl[3],
		   &grid[gridno].icon[i].nl[4],
		   &grid[gridno].icon[i].nl[5]);
	    for (j = 0; j < 6; j++) {
		grid[gridno].icon[i].nl[j]--;
		nn = grid[gridno].icon[i].nl[j];
		grid[gridno].ntype[nn] = (j > 2) ? MIDDLE : CORNER;
	    }
	    grid[gridno].icon[i].nn = 6;
	    grid[gridno].icon[i].ngeom = 3;
	    break;
	}
    }
    grid[gridno].gridtype = 1;	/* set gridtype to linear. */
    fclose(fp);
    set_grid_limits(gridno);
    default_isolines(&grid[gridno].ip);
    autoscale_isolines(grid[gridno].dmin, grid[gridno].dmax, &grid[gridno].ip);
    mindist = (grid[gridno].xmax - grid[gridno].xmin) * 0.01;
    return 1;
}

/*
 * Merge a grid to the edit grid
 */
int mergegrid(int gridno, char *fname)
{
    FILE *fp;
    char buf[1024];
    int n1, n2, n3, n4, i, j, itmp, type, nn, nnodes, nmel, nmnp;
    double *x, *y, *d;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In mergegrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    nn = grid[curgrid].nmnp;
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &nmel, &nmnp);
    x = (double *) malloc(nmnp * sizeof(double));
    y = (double *) malloc(nmnp * sizeof(double));
    d = (double *) malloc(nmnp * sizeof(double));
    for (i = 0; i < nmnp; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("Error reading table of nodes");
	    fclose(fp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %lf %lf %lf", &itmp, &x[i], &y[i], &d[i]);
    }
    add_nodes(curgrid, nmnp, x, y, d);
    for (i = 0; i < nmel; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("Error reading table of elements");
	    fclose(fp);
	    return 0;
	}
	convertchar(buf);
	sscanf(buf, "%d %d", &itmp, &type);
	switch (type) {
	case 3:
	    sscanf(buf, "%d %d %d %d %d", &itmp, &itmp, &n1, &n2, &n3);
	    n1--;
	    n2--;
	    n3--;
	    add_element(curgrid, nn + n1, nn + n2, nn + n3);
	    break;
	case 4:
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &itmp, &n1, &n2, &n3, &n4);
	    n1--;
	    n2--;
	    n3--;
	    n4--;
	    add_quad_element(curgrid, nn + n1, nn + n2, nn + n3, nn + n4);
	    break;
	}
    }
    fclose(fp);
    free(x);
    free(y);
    free(d);
    return 1;
}

/*
 * Write a grid in ASCII
 */
int writegrid(int gridno, char *fname)
{
    FILE *fp;
    int i, j;

    if ((fp = fopen(fname, "w")) == NULL) {
	errwin(stderr, "In writegrid, unable to open file %s\n", fname);
	fclose(fp);
	return 0;
    }
    fprintf(fp, fname);
    fprintf(fp, "\n%d %d\n", grid[gridno].nmel, grid[gridno].nmnp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	fprintf(fp, "%d %14.6lf %14.6lf %.7le\n", i + 1, grid[gridno].xord[i], grid[gridno].yord[i], grid[gridno].depth[i]);
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	switch (grid[gridno].icon[i].type) {
	case 3:		/* 2d triangle */
	    fprintf(fp, "%d 3 %d %d %d\n", i + 1,
		    grid[gridno].icon[i].nl[0] + 1,
		    grid[gridno].icon[i].nl[1] + 1,
		    grid[gridno].icon[i].nl[2] + 1);
	    break;
	case 4:		/* 2d quadrangle */
	    fprintf(fp, "%d 4 %d %d %d %d\n", i + 1,
		    grid[gridno].icon[i].nl[0] + 1,
		    grid[gridno].icon[i].nl[1] + 1,
		    grid[gridno].icon[i].nl[2] + 1,
		    grid[gridno].icon[i].nl[3] + 1);
	    break;
	case 6:		/* quadratic triangles */
	    fprintf(fp, "%d 6 %d %d %d %d %d %d\n", i + 1,
		    grid[gridno].icon[i].nl[0] + 1,
		    grid[gridno].icon[i].nl[1] + 1,
		    grid[gridno].icon[i].nl[2] + 1,
		    grid[gridno].icon[i].nl[3] + 1,
		    grid[gridno].icon[i].nl[4] + 1,
		    grid[gridno].icon[i].nl[5] + 1);
	    break;
	}
    }
    fclose(fp);
    return 1;
}

/*
 * Write a grid in binary format
 */
int writegridbin(int gridno, char *fname)
{
    FILE *fp;
    int i, itmp, *intptr;
    float *floatptr;
    char buf[1024];

    if ((fp = fopen(fname, "w")) == NULL) {
	sprintf(buf, "In writegrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    itmp = 4;
    write_int(&itmp, 1, fp);
    itmp = 2;			/* magic number for real*4 based binary grid
				 * files */
    write_int(&itmp, 1, fp);
    itmp = 4;
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    write_int(&grid[gridno].nmel, 1, fp);
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    write_int(&grid[gridno].nmnp, 1, fp);
    write_int(&itmp, 1, fp);

    floatptr = (float *) calloc(grid[gridno].nmnp, sizeof(float));

    if (floatptr == NULL) {
	errwin("Can't allocate memory for floatptr");
	Free_grid(gridno);
	fclose(fp);
	return 0;
    }
    itmp = grid[gridno].nmnp * sizeof(float);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	floatptr[i] = grid[gridno].xord[i];
    }
    write_float(floatptr, grid[gridno].nmnp, fp);
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	floatptr[i] = grid[gridno].yord[i];
    }
    write_float(floatptr, grid[gridno].nmnp, fp);
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	floatptr[i] = grid[gridno].depth[i];
    }
    write_float(floatptr, grid[gridno].nmnp, fp);
    write_int(&itmp, 1, fp);

    free(floatptr);

    intptr = (int *) calloc(grid[gridno].nmel, sizeof(int));

    if (intptr == NULL) {
	errwin("Can't allocate memory for intptr");
	Free_grid(gridno);
	free(floatptr);
	fclose(fp);
	return 0;
    }
    itmp = grid[gridno].nmel * sizeof(int);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	intptr[i] = grid[gridno].icon[i].nl[0];
    }
    write_int(intptr, grid[gridno].nmel, fp);
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	intptr[i] = grid[gridno].icon[i].nl[1];
    }
    write_int(intptr, grid[gridno].nmel, fp);
    write_int(&itmp, 1, fp);

    write_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	intptr[i] = grid[gridno].icon[i].nl[2];
    }
    write_int(intptr, grid[gridno].nmel, fp);
    write_int(&itmp, 1, fp);

    free(intptr);

    fclose(fp);
    return 1;
}

/*
 * Read a grid in binary format
 */
int readgridbin(int gridno, char *fname)
{
    FILE *fp;
    char buf[1024];
    int i, j, nn, itmp;
    int *intptr;
    float *floatptr;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readgrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    Free_grid(gridno);

    read_int(&itmp, 1, fp);
    read_int(&itmp, 1, fp);
    if (itmp != 2) {
	sprintf(buf, "Improper magic number for binary grid file = %d", itmp);
	errwin(buf);
	fclose(fp);
	return 0;
    }
    read_int(&itmp, 1, fp);

    read_int(&itmp, 1, fp);
    read_int(&grid[gridno].nmel, 1, fp);
    read_int(&itmp, 1, fp);

    read_int(&itmp, 1, fp);
    read_int(&grid[gridno].nmnp, 1, fp);
    read_int(&itmp, 1, fp);

    if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	fclose(fp);
	return 0;
    }
    floatptr = (float *) calloc(grid[gridno].nmnp, sizeof(float));

    if (floatptr == NULL) {
	errwin("Can't allocate memory for floatptr");
	Free_grid(gridno);
	fclose(fp);
	return 0;
    }
    itmp = grid[gridno].nmnp * sizeof(float);

    read_int(&itmp, 1, fp);
    read_float(floatptr, grid[gridno].nmnp, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].xord[i] = floatptr[i];
    }

    read_int(&itmp, 1, fp);
    read_float(floatptr, grid[gridno].nmnp, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].yord[i] = floatptr[i];
    }

    read_int(&itmp, 1, fp);
    read_float(floatptr, grid[gridno].nmnp, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].depth[i] = floatptr[i];
    }

    free(floatptr);

    intptr = (int *) calloc(grid[gridno].nmel, sizeof(int));

    if (intptr == NULL) {
	errwin("Can't allocate memory for intptr");
	Free_grid(gridno);
	free(floatptr);
	fclose(fp);
	return 0;
    }
    itmp = grid[gridno].nmel * sizeof(int);

    read_int(&itmp, 1, fp);
    read_int(intptr, grid[gridno].nmel, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].icon[i].nl[0] = intptr[i];
    }

    read_int(&itmp, 1, fp);
    read_int(intptr, grid[gridno].nmel, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].icon[i].nl[1] = intptr[i];
    }

    read_int(&itmp, 1, fp);
    read_int(intptr, grid[gridno].nmel, fp);
    read_int(&itmp, 1, fp);
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].icon[i].nl[2] = intptr[i];
	grid[gridno].icon[i].type = 3;
	for (j = 0; j < 3; j++) {
	    nn = grid[gridno].icon[i].nl[j];
	    grid[gridno].ntype[nn] = CORNER;
	}
	grid[gridno].icon[i].nn = 3;
	grid[gridno].icon[i].ngeom = 3;
    }
    fclose(fp);

    free(intptr);

    set_grid_limits(gridno);
    default_isolines(&grid[gridno].ip);
    autoscale_isolines(grid[gridno].dmin, grid[gridno].dmax, &grid[gridno].ip);
    return 1;
}

/*
 * Setting the limits of the grid
 */
void set_grid_limits(int gridno)
{
    int i;
    grid[gridno].xmin = grid[gridno].xmax = grid[gridno].xord[0];
    grid[gridno].ymin = grid[gridno].ymax = grid[gridno].yord[0];
    grid[gridno].dmin = grid[gridno].dmax = grid[gridno].depth[0];
    for (i = 1; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].xmin > grid[gridno].xord[i])
	    grid[gridno].xmin = grid[gridno].xord[i];
	if (grid[gridno].ymin > grid[gridno].yord[i])
	    grid[gridno].ymin = grid[gridno].yord[i];
	if (grid[gridno].xmax < grid[gridno].xord[i])
	    grid[gridno].xmax = grid[gridno].xord[i];
	if (grid[gridno].ymax < grid[gridno].yord[i])
	    grid[gridno].ymax = grid[gridno].yord[i];
	if (grid[gridno].dmin > grid[gridno].depth[i])
	    grid[gridno].dmin = grid[gridno].depth[i];
	if (grid[gridno].dmax < grid[gridno].depth[i])
	    grid[gridno].dmax = grid[gridno].depth[i];
    }
    grid[gridno].ip.cmin = grid[gridno].dmin;
    grid[gridno].ip.cmax = grid[gridno].dmax;
}

int readdepths(int gridno, char *fname)
{
    FILE *fp;
    char buf[1024];
    int i, itmp;
    double d;

    if ((fp = fopen(fname, "r")) == NULL) {
        errwin(stderr, "In readdepths(), unable to open file %s\n", fname);
        fclose(fp);
        return 0;
    }
    grid[gridno].edepth = (double *) malloc(grid[gridno].nmel * sizeof(double));
    if (grid[gridno].edepth != NULL) {
      for (i = 0; i < grid[gridno].nmel; i++) {
        fgets(buf, 255, fp);
        sscanf(buf, "%d %lf", &itmp, &d);
        grid[gridno].edepth[i] = d;
      }
        grid[gridno].depthflag = 1;
    } else {
    }
    fclose(fp);
    return 1;
}

/*
 * Better allocate and delete Grid routines
 */
Grid *NewGrid(int nmel, int nmnp);
void DeleteGrid(Grid * g);

Grid *NewGrid(int nmel, int nmnp)
{
    int i;
    Grid *g = (Grid *) malloc(sizeof(Grid));
    if (g == NULL) {
	return NULL;
    }
    g->nmel = nmel;
    g->nmnp = nmnp;
    g->gactive = 1;
    g->bactive = 0;
    g->xord = (double *) calloc(nmnp, sizeof(double));
    g->yord = (double *) calloc(nmnp, sizeof(double));
    g->depth = (double *) calloc(nmnp, sizeof(double));
    g->nodecon = (conlist *) calloc(nmnp, sizeof(conlist));
    g->ellist = (int *) calloc(nmel, sizeof(int));
    g->nlist = (int *) calloc(nmnp, sizeof(int));
    g->ntype = (int *) calloc(nmnp, sizeof(int));
    g->icon = (Element *) calloc(nmel, sizeof(Element));
    return g;
}

void DeleteGrid(Grid * g)
{
    int i;
    if (g == NULL) {
	return;
    }
    free(g->xord);
    free(g->yord);
    free(g->depth);
    free(g->nodecon);
    free(g->ellist);
    free(g->nlist);
    free(g->ntype);
    free(g->icon);
    free(g);
}

/*
 * read a finite difference grid
 */
int readpomgrid(int gridno, char *fname)
{
    FILE *fp;
    char buf[256];
    int i, j, itmp, type, nn, nnodes;
    int cnt, nx, ny;
    int idx, idy;
    extern int doadcirc;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "In readgrid, unable to open file %s\n", fname);
	errwin(buf);
	return 0;
    }
    Free_grid(gridno);

    fgets(buf, 255, fp);
    strcpy(grid[gridno].alphid, buf);
    i = strlen(grid[gridno].alphid);
    grid[gridno].alphid[i - 1] = '\0';
    fgets(buf, 255, fp);
    convertchar(buf);
    sscanf(buf, "%d %d", &nx, &ny);
    grid[gridno].nmnp = nx * ny;
    grid[gridno].nmel = (nx - 1) * (ny - 1);
    if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	fclose(fp);
	return 0;
    }
    cnt = 0;
    grid[gridno].fdgrid = 1;
    grid[gridno].nx = nx;
    grid[gridno].ny = ny;
    for (j = 0; j < ny; j++) {
	for (i = 0; i < nx; i++) {
	    if (fgets(buf, 255, fp) == NULL) {
		errwin("Error reading table of nodes");
		Free_grid(gridno);
		fclose(fp);
		return 0;
	    }
	    convertchar(buf);
	    sscanf(buf, "%d %d %lf %lf %lf", &itmp, &itmp, &grid[gridno].xord[cnt], &grid[gridno].yord[cnt], &grid[gridno].depth[cnt]);
	    cnt++;
	}
    }
    cnt = 0;
    for (j = 0; j < ny - 1; j++) {
	for (i = 0; i < nx - 1; i++) {
	    grid[gridno].icon[cnt].type = 4;
	    grid[gridno].icon[cnt].wetdry = 0;
	    grid[gridno].icon[cnt].nl[0] = j * nx + i;
	    grid[gridno].icon[cnt].nl[1] = j * nx + i + 1;
	    grid[gridno].icon[cnt].nl[2] = (j + 1) * nx + i + 1;
	    grid[gridno].icon[cnt].nl[3] = (j + 1) * nx + i;
	    grid[gridno].icon[cnt].nn = 4;
	    grid[gridno].icon[cnt].ngeom = 4;
	    cnt++;
	}
    }
/*
 * read dry nodes
 */
    while (fgets(buf, 255, fp) != NULL) {
	cnt = sscanf(buf, "%d %d", &i, &j);
	if (cnt == 2) {
	    cnt = (j - 1) * (nx - 1) + i - 1;
	    grid[gridno].icon[cnt].wetdry = 1;
	}
    }
    grid[gridno].gridtype = 1;	/* set gridtype to linear. */
    fclose(fp);
    set_grid_limits(gridno);
    default_isolines(&grid[gridno].ip);
    autoscale_isolines(grid[gridno].dmin, grid[gridno].dmax, &grid[gridno].ip);
    mindist = (grid[gridno].xmax - grid[gridno].xmin) * 0.01;
    return 1;
}

#ifdef HAVE_NETCDF
/*
 * netcdf file definition
 * netcdf grid {
 * dimensions:
 *         nmnp = 30495 ; // number of nodes
 *         nmel = 54237 ; // number of elements
 *         nnodes = 54237 ; // number of nodes mentioned in table of elements
 * variables:
 *         double x(nmnp) ;
 *                 x:units = "meters" ;
 *                 x:description = "x-coordinate" ;
 *         double y(nmnp) ;
 *                 y:units = "meters" ;
 *                 y:description = "y-coordinate" ;
 *         double z(nmnp) ;
 *                 z:units = "meters" ;
 *                 z:description = "depth" ;
 *         int nodelist(nnodes) ;
 *                 nodelist:description = "Table of elements." ;
 *         int nnodes(nmel, 2) ;
 *                 nnodes:description = "Where in the list of nodes the element begins and
 *  the type of element (3 triangle, 4 quadrangle, 6 quadratic triangles)." ;
 * }
 */

#include "netcdf.h"

/* check errors on netcdf file operations */
void check_err(const int stat, const int line, const char *file)
{
    if (stat != NC_NOERR) {
	(void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    }
}

/* Create the netcdf file */
void CreateGridnetcdf(char *fname, int gridno)
{
    int i;
    int ncid;			/* netCDF id */

    /* dimension ids */
    int nmnp_dim;
    int nmel_dim;
    int elements_dim;

    /* dimension lengths */
    size_t nmnp_len = grid[gridno].nmnp;
    size_t nmel_len = grid[gridno].nmel;
    size_t elements_len;

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;
    int elements_id;

    /* variable shapes */
    int x_dims[1];
    int y_dims[1];
    int z_dims[1];
    int elements_dims[1];
   
    int stat;

    int sum = 0;
    /* 
     * find the total number of nodes mentioned in the table of 
     * elements + 1 to give the number of nodes per element
     */
    for (i = 0; i < grid[gridno].nmel; i++) {
	sum += grid[gridno].icon[i].nn + 1;
    }
    elements_len = sum;

    /* enter define mode */
    stat = nc_create(fname, NC_CLOBBER, &ncid);
    check_err(stat, __LINE__, __FILE__);

    /* define dimensions */
    stat = nc_def_dim(ncid, "nmnp", nmnp_len, &nmnp_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_dim(ncid, "nmel", nmel_len, &nmel_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_def_dim(ncid, "elements", elements_len, &elements_dim);
    check_err(stat, __LINE__, __FILE__);

    /* define variables */
    x_dims[0] = nmnp_dim;
    stat = nc_def_var(ncid, "x", NC_DOUBLE, 1, x_dims, &x_id);
    check_err(stat, __LINE__, __FILE__);

    y_dims[0] = nmnp_dim;
    stat = nc_def_var(ncid, "y", NC_DOUBLE, 1, y_dims, &y_id);
    check_err(stat, __LINE__, __FILE__);

    z_dims[0] = nmnp_dim;
    stat = nc_def_var(ncid, "z", NC_DOUBLE, 1, z_dims, &z_id);
    check_err(stat, __LINE__, __FILE__);

    elements_dims[0] = elements_dim;
    stat = nc_def_var(ncid, "elements", NC_INT, 1, elements_dims, &elements_id);
    check_err(stat, __LINE__, __FILE__);

    /* leave define mode */
    stat = nc_enddef(ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Write the grid to the netcdf file fname */
void WriteGridnetcdf(char *fname, int gridno)
{
    int i, j, cnt;
    int *ellist;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int nmnp_dim;
    int nmel_dim;
    int elements_dim;

    /* dimension lengths */
    size_t nmnp_len;
    size_t nmel_len;
    size_t elements_len;
    size_t start[2], count[2];

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;
    int elements_id;

    int stat;

    stat = nc_open(fname, NC_WRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "nmel", &nmel_dim);
    stat = nc_inq_dim(ncid, nmel_dim, "nmel", &nmel_len);
    stat = nc_inq_dimid(ncid, "nmnp", &nmnp_dim);
    stat = nc_inq_dim(ncid, nmnp_dim, "nmnp", &nmnp_len);
    stat = nc_inq_dimid(ncid, "elements", &elements_dim);
    stat = nc_inq_dim(ncid, elements_dim, "elements", &elements_len);
    stat = nc_inq_varid(ncid, "x", &x_id);
    stat = nc_inq_varid(ncid, "y", &y_id);
    stat = nc_inq_varid(ncid, "z", &z_id);
    stat = nc_inq_varid(ncid, "elements", &elements_id);

    ellist = (int *) malloc(elements_len * sizeof(int));
    if (ellist == NULL) {
	stat = nc_close(ncid);
	return;
    }
    start[0] = 0;
    count[0] = nmnp_len;
    stat = nc_put_vara_double(ncid, x_id, start, count, grid[gridno].xord);
    stat = nc_put_vara_double(ncid, y_id, start, count, grid[gridno].yord);
    stat = nc_put_vara_double(ncid, z_id, start, count, grid[gridno].depth);

    cnt = 0;
    for (i = 0; i < grid[gridno].nmel; i++) {
	ellist[cnt++] = grid[gridno].icon[i].nn;
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    ellist[cnt++] = grid[gridno].icon[i].nl[j];
	}
    }
    /*
     * Oops, number of items in elements differs between the grid and the
     * file, bail out.
     */
    if (cnt != elements_len) {
	free(ellist);
	stat = nc_close(ncid);
	return;
    }
    start[0] = 0;
    count[0] = cnt;
    stat = nc_put_vara_int(ncid, elements_id, start, count, ellist);
    free(ellist);
    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);
}

/* Read the grid from the netcdf file fname */
void ReadGridnetcdf(char *fname, int gridno)
{
    int i, j, cnt, nel, newelement;
    int *ellist;

    int ncid;			/* netCDF id */

    /* dimension ids */
    int nmnp_dim;
    int nmel_dim;
    int elements_dim;

    /* dimension lengths */
    size_t nmnp_len;
    size_t nmel_len;
    size_t elements_len;
    size_t start[2], count[2];

    /* variable ids */
    int x_id;
    int y_id;
    int z_id;
    int elements_id;

    int stat;

    stat = nc_open(fname, NC_NOWRITE, &ncid);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "nmel", &nmel_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, nmel_dim, "nmel", &nmel_len);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "nmnp", &nmnp_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, nmnp_dim, "nmnp", &nmnp_len);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dimid(ncid, "elements", &elements_dim);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_dim(ncid, elements_dim, "elements", &elements_len);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "x", &x_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "y", &y_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "z", &z_id);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_inq_varid(ncid, "elements", &elements_id);
    check_err(stat, __LINE__, __FILE__);

    ellist = (int *) malloc(elements_len * sizeof(int));
    if (ellist == NULL) {
	stat = nc_close(ncid);
	check_err(stat, __LINE__, __FILE__);
	return;
    }

    Free_grid(gridno);

    grid[gridno].alphid[0] = '\0';
    if (!Allocate_grid(gridno, nmel_len, nmnp_len)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	free(ellist);
	stat = nc_close(ncid);
	return;
    }
    start[0] = 0;
    count[0] = elements_len;
    stat = nc_get_vara_int(ncid, elements_id, start, count, ellist);
    check_err(stat, __LINE__, __FILE__);

/*
 * Read table of nodes
 */
    start[0] = 0;
    count[0] = nmnp_len;
    stat = nc_get_vara_double(ncid, x_id, start, count, grid[gridno].xord);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, y_id, start, count, grid[gridno].yord);
    check_err(stat, __LINE__, __FILE__);
    stat = nc_get_vara_double(ncid, z_id, start, count, grid[gridno].depth);
    check_err(stat, __LINE__, __FILE__);

/*
 * Read table of elements
 */
    cnt = 0;
    newelement = 1;
    nel = 0;
    for (i = 0; i < elements_len; i++) {
	if (newelement) {
	    grid[gridno].icon[nel].nn = ellist[i];
	    grid[gridno].icon[nel].type = ellist[i]; /* for now, type and nn are the same */
	    cnt = 0;
	    newelement = 0;
	} else {
	    grid[gridno].icon[nel].nl[cnt] = ellist[i];
	    cnt++;
	    if (cnt == grid[gridno].icon[nel].nn) {
		newelement = 1;
		nel++;
	    }
	}
    }

    free(ellist);

    stat = nc_close(ncid);
    check_err(stat, __LINE__, __FILE__);

    grid[gridno].gridtype = 1;  /* set gridtype to linear. */
    set_grid_limits(gridno);
    default_isolines(&grid[gridno].ip);
    autoscale_isolines(grid[gridno].dmin, grid[gridno].dmax, &grid[gridno].ip);
    mindist = (grid[gridno].xmax - grid[gridno].xmin) * 0.01;
}

#endif
