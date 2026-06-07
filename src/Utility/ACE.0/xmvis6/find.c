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
 * find.c - routines for finding nodes, elements
 *
 */

#ifndef lint
static char RCSid[] = "$Id: find.c,v 1.7 2003/09/25 04:35:15 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"
#include <math.h>

int belel4(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
int belel(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3);
double belel_tol = -1e-04;

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

void FindNearestNode(Grid * g, double x, double y, int *ind)
{
    int i;
    double radius;
    double tmp, *xord = g->xord, *yord = g->yord;

    if (g->nmnp > 0) {
	radius = hypot((x - xord[0]), (y - yord[0]));
	*ind = 0;
	for (i = 0; i < g->nmnp; i++) {
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

void FindElement(Grid * g, double x, double y, int *elem)
{
    int i;
    int n0, n1, n2, n3;
    double *xord = g->xord, *yord = g->yord;
    *elem = -1;

    for (i = 0; i < g->nmel; i++) {
	if (g->icon[i].ngeom == 3) {
	    n0 = g->icon[i].nl[0];
	    n1 = g->icon[i].nl[1];
	    n2 = g->icon[i].nl[2];
	    if (ElioInsideElement(x, y, xord[n0], yord[n0], xord[n1], yord[n1], xord[n2], yord[n2])) {
		*elem = i;
		return;
	    }
	} else if (g->icon[i].ngeom == 4) {
	    n0 = g->icon[i].nl[0];
	    n1 = g->icon[i].nl[1];
	    n2 = g->icon[i].nl[2];
	    n3 = g->icon[i].nl[3];
	    if (belel4(x, y, xord[n0], yord[n0], xord[n1], yord[n1], xord[n2], yord[n2], xord[n3], yord[n3])) {
		*elem = i;
		return;
	    }
	}
    }
    *elem = -1;
}

void find_element(int gridno, double x, double y, int *elem)
{
    int i;
    int n0, n1, n2, n3;
    double *xord = grid[gridno].xord, *yord = grid[gridno].yord;
    *elem = -1;

    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].ngeom == 3) {
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    if (belel(x, y, xord[n0], yord[n0], xord[n1], yord[n1], xord[n2], yord[n2])) {
		*elem = i;
		return;
	    }
	} else if (grid[gridno].icon[i].ngeom == 4) {
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    if (belel4(x, y, xord[n0], yord[n0], xord[n1], yord[n1], xord[n2], yord[n2], xord[n3], yord[n3])) {
		*elem = i;
		return;
	    }
	}
    }
    *elem = -1;
}

void find_nearest_element(int gridno, double x, double y, int *elem)
{
    int i;
    double radius;
    double xg, yg, tmp;

    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xg, &yg);
	if (i == 0) {
	    radius = tmp = hypot((x - xg), (y - yg));
	    *elem = 0;
	} else {
	    tmp = hypot((x - xg), (y - yg));
	}
	if (tmp < radius) {
	    radius = tmp;
	    *elem = i;
	}
    }
}

int belel(double xp, double yp, double xp1, double yp1, double xp2, double yp2, double xp3, double yp3)
{
    if ((xp - xp1) * (yp1 - yp2) + (yp - yp1) * (xp2 - xp1) < belel_tol)
	return 0;
    if ((xp - xp2) * (yp2 - yp3) + (yp - yp2) * (xp3 - xp2) < belel_tol)
	return 0;
    if ((xp - xp3) * (yp3 - yp1) + (yp - yp3) * (xp1 - xp3) < belel_tol)
	return 0;
    /*
    int iv1 = 0;
    int iv2 = 0;
    int iv3 = 0;
    iv1 =  ((xp - xp1) * (yp1 - yp2) + (yp - yp1) * (xp2 - xp1) < belel_tol);
    if (iv1 != 0) { return 0; }
    iv2 =  ((xp - xp2) * (yp2 - yp3) + (yp - yp2) * (xp3 - xp2) < belel_tol);
    if (iv2 != 0) { return 0; }
    iv3 =  ((xp - xp3) * (yp3 - yp1) + (yp - yp3) * (xp1 - xp3) < belel_tol);
    if (iv3 != 0) { return 0; }
    printf("In belel logical: %d %d %d\n", iv1, iv2, iv3);
    printf("In belel logical: %d\n", (xp - xp1) * (yp1 - yp2) + (yp - yp1) * (xp2 - xp1) < belel_tol);
    printf("In belel logical: %d\n", (xp - xp2) * (yp2 - yp3) + (yp - yp2) * (xp3 - xp2) < belel_tol);
    printf("In belel logical: %d\n", (xp - xp3) * (yp3 - yp1) + (yp - yp3) * (xp1 - xp3) < belel_tol);
    printf("In belel: %lf %lf %lf %lf %lf %lf %lf %lf\n",  xp,  yp,  xp1,  yp1,  xp2,  yp2,  xp3,  yp3);
    printf("In belel: %lf %lf %lf %lf %lf\n", (xp - xp1) , (yp1 - yp2) , (yp - yp1) , (xp2 - xp1), (xp - xp1) * (yp1 - yp2) + (yp - yp1) * (xp2 - xp1));
    printf("In belel: %lf %lf %lf %lf %lf\n", (xp - xp2) , (yp2 - yp3) , (yp - yp2) , (xp3 - xp2), (xp - xp2) * (yp2 - yp3) + (yp - yp2) * (xp3 - xp2));
    printf("In belel: %lf %lf %lf %lf %lf\n", (xp - xp3) , (yp3 - yp1) , (yp - yp3) , (xp1 - xp3), (xp - xp3) * (yp3 - yp1) + (yp - yp3) * (xp1 - xp3));
    */
    return 1;
}

int belel4(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    if ((x - x1) * (y1 - y2) + (y - y1) * (x2 - x1) < belel_tol)
	return 0;
    if ((x - x2) * (y2 - y3) + (y - y2) * (x3 - x2) < belel_tol)
	return 0;
    if ((x - x3) * (y3 - y4) + (y - y3) * (x4 - x3) < belel_tol)
	return 0;
    if ((x - x4) * (y4 - y1) + (y - y4) * (x1 - x4) < belel_tol)
	return 0;
    return 1;
}
