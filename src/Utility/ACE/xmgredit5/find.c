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
 * find.c - routines for finding nodes, elements
 *
 */

#ifndef lint
static char RCSid[] = "$Id: find.c,v 1.4 2003/08/15 00:37:57 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

extern double belel_tol;
int belel(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3);

void find_nearest_node(int gridno, double x, double y, int *ind)
{
    int i;
    double radius;
    double tmp, *xord = grid[gridno].xord, *yord = grid[gridno].yord;

    if (grid[gridno].nmnp > 0) {
	radius = hypot((x - xord[0]), (y - yord[0]));
	*ind = 0;
	for (i = 1; i < grid[gridno].nmnp; i++) {
	    tmp = hypot((x - xord[i]), (y - yord[i]));
	    if (tmp < radius) {
		radius = tmp;
		*ind = i;
	    }
	}
    } else {
	*ind = -1;
    }
}

void find_element(int gridno, double x, double y, int *elem)
{
    int i;
    int n0, n1, n2;
    double *xord = grid[gridno].xord, *yord = grid[gridno].yord;

    for (i = 0; i < grid[gridno].nmel; i++) {
	if (gridbelel(gridno, i, x, y)) {
	    *elem = i;
	    return;
	}
    }
    *elem = -1;
    return;
}

int inside_element(int gridno, int elem, double x, double y)
{
    return gridbelel(gridno, elem, x, y);
}

void find_nearest_element(int gridno, double x, double y, int *elem)
{
    int i;
    double xg, yg, tmp, radius = 1e307;

    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xg, &yg);
	tmp = hypot((x - xg), (y - yg));
	if (i == 0) {
	    radius = tmp;
	    *elem = 0;
	}
	if (tmp < radius) {
	    radius = tmp;
	    *elem = i;
	}
    }
    return;
}

int belel(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3)
{
/* printf("1 %lf %lf\n", (x - x1) * (y1 - y2) + (y - y1) * (x2 - x1), belel_tol); */
    if ((x - x1) * (y1 - y2) + (y - y1) * (x2 - x1) < belel_tol) {
	return 0;
    }
/* printf("2 %lf %lf\n", (x - x2) * (y2 - y3) + (y - y2) * (x3 - x2), belel_tol); */
    if ((x - x2) * (y2 - y3) + (y - y2) * (x3 - x2) < belel_tol) {
	return 0;
    }
/* printf("3 %lf %lf\n", (x - x3) * (y3 - y1) + (y - y3) * (x1 - x3), belel_tol); */
    if ((x - x3) * (y3 - y1) + (y - y3) * (x1 - x3) < belel_tol) {
	return 0;
    }
    return 1;
}

int belel4(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    if ((x - x1) * (y1 - y2) + (y - y1) * (x2 - x1) < belel_tol) {
	return 0;
    }
    if ((x - x2) * (y2 - y3) + (y - y2) * (x3 - x2) < belel_tol) {
	return 0;
    }
    if ((x - x3) * (y3 - y4) + (y - y3) * (x4 - x3) < belel_tol) {
	return 0;
    }
    if ((x - x4) * (y4 - y1) + (y - y4) * (x1 - x4) < belel_tol) {
	return 0;
    }
    return 1;
}

int gridbelel(int gridno, int el, double x, double y)
{
    double x1, y1, x2, y2, x3, y3, x4, y4;
    int n0, n1, n2, n3;
    double *xord = grid[gridno].xord, *yord = grid[gridno].yord;
    n0 = grid[gridno].icon[el].nl[0];
    n1 = grid[gridno].icon[el].nl[1];
    n2 = grid[gridno].icon[el].nl[2];
    x1 = xord[n0];
    x2 = xord[n1];
    x3 = xord[n2];
    y1 = yord[n0];
    y2 = yord[n1];
    y3 = yord[n2];
    if (grid[gridno].icon[el].nn == 4) {
	n3 = grid[gridno].icon[el].nl[3];
	x4 = xord[n3];
	y4 = yord[n3];
	return belel4(x, y, x1, y1, x2, y2, x3, y3, x4, y4);
    } else {
	return belel(x, y, x1, y1, x2, y2, x3, y3);
    }
}
