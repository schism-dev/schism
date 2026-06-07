/************************************************************************
*
* Copyright 1990-2003 Oregon Health and Science University
*
*************************************************************************
*
* tritest.c - interface to Triangle 
*
* Contents:
*
*
*************************************************************************
*/
/*****************************************************************************/
/*                                                                           */
/*  (tricall.c)                                                              */
/*                                                                           */
/*  Example program that demonstrates how to call Triangle.                  */
/*                                                                           */
/*  Accompanies Triangle Version 1.3                                         */
/*  July 19, 1996                                                            */
/*                                                                           */
/*  This file is placed in the public domain (but the file that it calls     */
/*  is still copyrighted!) by                                                */
/*  Jonathan Richard Shewchuk                                                */
/*  School of Computer Science                                               */
/*  Carnegie Mellon University                                               */
/*  5000 Forbes Avenue                                                       */
/*  Pittsburgh, Pennsylvania  15213-3891                                     */
/*  jrs@cs.cmu.edu                                                           */
/*                                                                           */
/*****************************************************************************/

#ifndef lint
static char RCSid[] = "$Id: tritest.c,v 1.3 2004/01/06 00:08:06 pturner Exp $";
#endif

#define REAL double

#include <stdio.h>
#include <math.h>
#include <time.h>

#include "triangle.h"

#include "defines.h"
#include "globals.h"

#ifndef _STDLIB_H_
extern void *malloc();
extern void free();
#endif				/* _STDLIB_H_ */

/* testing the voronoi computation */
/*
static double *vorx, *vory;
static int nvor, nvoredges, *vorn1, *vorn2;
void drawvoronoi(int n, int nedges, double *x, double *y, int *n1, int *n2);

void do_vor(void)
{
    if (nvor) {
	drawvoronoi(nvor, nvoredges, vorx, vory, vorn1, vorn2);
    }
}

void drawvoronoi(int n, int nedges, double *x, double *y, int *n1, int *n2)
{
    int i;
    setcolor(9);
    for (i=0;i<nedges;i++) {
	my_move2(x[n1[i]], y[n1[i]]);
	my_draw2(x[n2[i]], y[n2[i]]);
    }
    setcolor(10);
    for (i=0;i<n;i++) {
	my_circlefilled(x[i], y[i]);
    }
    setcolor(1);
}
*/

double find_nearest_buildpt(int bno, double x, double y, int *ind);

void tritest(int bno, int gridno, int include_boundary, char *args)
{
    char buf[512];
    int i, j, ne, bndno, cnt, ind, ind1, ind2;
    double x1, y1, x2, y2;
    struct triangulateio in, mid, out, vorout;
    int nb = build[bno].nbuild;
    double *bx = build[bno].bx;
    double *by = build[bno].by;

/*
    nvor = nvoredges = 0;
    vorx = vory = NULL;
    vorn1 = vorn2 = NULL;
*/

    in.pointlist = NULL;
    in.pointattributelist = NULL;
    in.pointmarkerlist = NULL;
    in.regionlist = NULL;
    mid.pointlist = NULL;
    mid.pointattributelist = NULL;
    mid.pointmarkerlist = NULL;
    mid.trianglelist = NULL;
    mid.triangleattributelist = NULL;
    mid.trianglearealist = NULL;
    mid.neighborlist = NULL;
    mid.segmentlist = NULL;
    mid.segmentmarkerlist = NULL;
    mid.edgelist = NULL;
    mid.edgemarkerlist = NULL;

    vorout.pointlist = NULL;
    vorout.pointattributelist = NULL;
    vorout.edgelist = NULL;
    vorout.normlist = NULL;

    in.numberofpoints = nb;
    in.numberofpointattributes = 0;
    in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));
    in.pointmarkerlist = (int *) NULL;
    in.numberofsegments = 0;
    in.numberofholes = 0;
    in.numberofregions = 0;
    in.segmentlist = NULL;
    in.segmentmarkerlist = NULL;

/* load points */
    for (i = 0; i < nb; i++) {
	in.pointlist[2 * i] = bx[i];
	in.pointlist[2 * i + 1] = by[i];
    }

/* load boundary edges */
    ne = 0;
    if (include_boundary) {
        for (i = 0; i < grid[gridno].nbounds; i++) {
            bndno = grid[gridno].boundaries[i];
            ne += boundary[bndno].nbpts;
        }
    }
    in.numberofsegments = ne;
    if (in.numberofsegments > 0) {
        in.segmentlist = (int *) malloc(in.numberofsegments * 2 * sizeof(int));
        in.segmentmarkerlist = (int *) malloc(in.numberofsegments * sizeof(int));
        if (include_boundary) {
	    cnt = 0;
            for (i = 0; i < grid[gridno].nbounds; i++) {
                bndno = grid[gridno].boundaries[i];
                for (j = 0; j < boundary[bndno].nbpts; j++) {
                    in.segmentlist[2 * cnt] = boundary[bndno].nodes[j] + 1;
                    in.segmentlist[2 * cnt + 1] = 
			boundary[bndno].nodes[(j + 1) % boundary[bndno].nbpts] + 1;
                    in.segmentmarkerlist[cnt] = 0;
                    cnt++;
                }
            }
        }
    }

    mid.pointlist = (REAL *) NULL;	/* Not needed if -N switch used. */
    mid.pointattributelist = (REAL *) NULL;
    mid.pointmarkerlist = (int *) NULL;		/* Not needed if -N or -B switch used. */
    mid.trianglelist = (int *) NULL;	/* Not needed if -E switch used. */

/*
    switch (method) {
    case 0:
	strcpy(buf, "pa0.1");
	strcpy(buf, "pq30X");
	break;
    case 1:
	strcpy(buf, "pi");
	break;
    case 2:
	strcpy(buf, "pF");
	break;
    }
*/

    triangulate(args, &in, &mid, NULL);
    /*triangulate(buf, &in, &mid, &vorout);*/
/*
    printf("Number of input points = %d, number of output points = %d\n", in.numberofpoints, mid.numberofpoints);
*/

    nels = mid.numberoftriangles;
    Free_grid_only(gridno);
    Allocate_grid_only(gridno, nels, mid.numberofpoints);
    for (i = 0; i < mid.numberofpoints; i++) {
        grid[gridno].xord[i] = mid.pointlist[2 * i];
        grid[gridno].yord[i] = mid.pointlist[2 * i + 1];
	find_nearest_buildpt(bno, 
			grid[gridno].xord[i], grid[gridno].yord[i], &ind);
        grid[gridno].depth[i] = build[bno].db[ind];
    }
    for (i = 0; i < nels; i++) {
        grid[gridno].icon[i].type = 3;
        grid[gridno].icon[i].nn = 3;
        grid[gridno].icon[i].ngeom = 3;
	for (j = 0; j < mid.numberofcorners; j++) {
            grid[gridno].icon[i].nl[j] = mid.trianglelist[i * mid.numberofcorners + j] - 1;
	}
    }
    CleanGrid(&grid[gridno]);
    DeleteElements(&grid[gridno]);
    setlimits_grid(gridno);

    free(in.pointlist);
    free(in.pointattributelist);
    free(in.pointmarkerlist);
    free(in.regionlist);
    free(mid.pointlist);
    free(mid.pointattributelist);
    free(mid.pointmarkerlist);
    free(mid.trianglelist);
    free(mid.triangleattributelist);
    free(mid.trianglearealist);
    free(mid.neighborlist);
    free(mid.segmentlist);
    free(mid.segmentmarkerlist);
    free(mid.edgelist);
    free(mid.edgemarkerlist);
}
