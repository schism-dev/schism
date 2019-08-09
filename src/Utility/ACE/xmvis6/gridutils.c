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
 * delete/split nodes/elements
 *
 */

#ifndef lint
static char RCSid[] = "$Id: gridutils.c,v 1.4 2003/08/05 22:36:11 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

static double sm, sb;

char item_buf[256];

double get_depth_element(int gno, int ind, double x, double y);

int comparea(int gridno, int n1, int n2, int n3)
{
    double a[2], b[2], ar;

    a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
    a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
    b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
    b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    if (ar <= 0) {
	return 0;
    }
    return 1;
}

double area(int gridno, int elem)
{
    int n1, n2, n3;
    double a[2], b[2], ar = 0.0;

    n1 = grid[gridno].icon[elem].nl[0];
    n2 = grid[gridno].icon[elem].nl[1];
    n3 = grid[gridno].icon[elem].nl[2];
    if (grid[gridno].icon[elem].ngeom == 3) {
	a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
	a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
	b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
	b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
	ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    } else if (grid[gridno].icon[elem].ngeom == 4) {
    }
    return ar;
}

double area_from_nodes(int gridno, int n1, int n2, int n3)
{
    double a[2], b[2], ar;

    a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
    a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
    b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
    b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    return ar;
}

int retelement(int gridno, double wx, double wy)
{
    int ind;
    find_element(gridno, wx, wy, &ind);
    return ind;
}

void getelement(int gridno, double wx, double wy)
{
    int elem, n0, n1, n2;

    find_element(gridno, wx, wy, &elem);
    if (elem) {
	n0 = grid[gridno].icon[elem].nl[0];
	n1 = grid[gridno].icon[elem].nl[1];
	n2 = grid[gridno].icon[elem].nl[2];
	sprintf(item_buf, "(x, y, element, [n1, n2, n3]) = (%7.2lf, %7.2lf, %d, [%d %d %d])", wx, wy, elem, n0, n1, n2);
    } else {
	sprintf(item_buf, "(x, y, No element found)", wx, wy);
    }
}

void getnodes(int gridno, double wx, double wy)
{
    int ind;
    double px, py, d;

    find_nearest_node(gridno, wx, wy, &ind);
    if (ind >= 0) {
	px = grid[gridno].xord[ind];
	py = grid[gridno].yord[ind];
	d = grid[gridno].depth[ind];
	sprintf(item_buf, "(x,y), depth, node) = (%7.2lf,%7.2lf), %9.3lf, %d)", px, py, d, ind + 1);
    }
}

void getall(int gridno, double wx, double wy)
{
    int ind, elem;
    double px, py;

    find_nearest_node(gridno, wx, wy, &ind);
    find_element(gridno, wx, wy, &elem);
    sprintf(item_buf, "(node,element) = (%1d, %1d)", ind, elem);
}

void get_center(int gridno, int elem, double *cx, double *cy)
{
    double x = 0.0, y = 0.0;
    int i, n;

    for (i = 0; i < grid[gridno].icon[elem].ngeom; i++) {
	n = grid[gridno].icon[elem].nl[i];
	x += grid[gridno].xord[n];
	y += grid[gridno].yord[n];
    }
    *cx = x / grid[gridno].icon[elem].ngeom;
    *cy = y / grid[gridno].icon[elem].ngeom;
}

void get_grid_center(Grid * g, int elem, double *cx, double *cy)
{
    double x = 0.0, y = 0.0;
    int i, n;

    for (i = 0; i < g->icon[elem].ngeom; i++) {
	n = g->icon[elem].nl[i];
	x += g->xord[n];
	y += g->yord[n];
    }
    *cx = x / g->icon[elem].ngeom;
    *cy = y / g->icon[elem].ngeom;
}

void get_average_depth(int gridno, int elem, double *depth)
{
    double d = 0.0;
    int i, n;

    for (i = 0; i < grid[gridno].icon[elem].ngeom; i++) {
	n = grid[gridno].icon[elem].nl[i];
	d += grid[gridno].depth[n];
    }
    *depth = d / grid[gridno].icon[elem].ngeom;
}

double get_depth_element(int gno, int ind, double x, double y)
{
    double x1, x2, x3, y1, y2, y3, x4, y4;
    double h, r, s, area2, a[3], b[3], g[3];
    int n0, n1, n2, n3;

    if (grid[gno].depthflag && grid[gno].edepth != NULL) {
	return grid[gno].edepth[ind];
    }
    n0 = grid[gno].icon[ind].nl[0];
    n1 = grid[gno].icon[ind].nl[1];
    n2 = grid[gno].icon[ind].nl[2];
    x1 = grid[gno].xord[n0];
    x2 = grid[gno].xord[n1];
    x3 = grid[gno].xord[n2];
    y1 = grid[gno].yord[n0];
    y2 = grid[gno].yord[n1];
    y3 = grid[gno].yord[n2];
    if (grid[gno].icon[ind].ngeom == 3) {
	a[0] = x2 * y3 - x3 * y2;
	a[1] = x3 * y1 - x1 * y3;
	a[2] = x1 * y2 - x2 * y1;
	b[0] = x2 - y3;
	b[1] = y3 - y1;
	b[2] = y1 - y2;
	g[0] = x3 - x2;
	g[1] = x1 - x3;
	g[2] = x2 - x1;
	area2 = a[0] + a[1] + a[2];
	s = (a[2] + b[2] * x + g[2] * y) / area2;
	r = (a[1] + b[1] * x + g[1] * y) / area2;
	h = 1.0 - r - s;
	return (grid[gno].depth[n0] * h + grid[gno].depth[n1] * r + grid[gno].depth[n2] * s);
    } else if (grid[gno].icon[ind].ngeom == 4) {
	double d = 0.0, xi, eta, w[4];
	int i;
	n3 = grid[gno].icon[ind].nl[3];
	x4 = grid[gno].xord[n3];
	y4 = grid[gno].yord[n3];
	ibilinear(x1, x2, x3, x4, y1, y2, y3, y4, x, y, &xi, &eta, w);
	/*
	   printf("       x1 = %lf\n       x2 = %lf\n       x3 = %lf\n       x4 = %lf\n       y1 = %lf\n       y2 = %lf\n        y3 = %lf\n       y4 = %lf\n       x = %lf\n       y = %lf\n", x1, x2, x3, x4, y1, y2, y3, y4, x, y);
	 */
	for (i = 0; i < 4; i++) {
	    d += w[i] * grid[gno].depth[grid[gno].icon[ind].nl[i]];
	    /*
	       printf("%lf %lf %lf\n", d, w[i], grid[gno].depth[grid[gno].icon[ind].nl[i]]);
	     */
	}
	return d;
    } else {
	return -999999.0;
    }
}

void xy2rs(double x, double y, double area2, double *a, double *b, double *g, double *r, double *s)
{
    *s = (a[2] + b[2] * x + g[2] * y) / area2;
    *r = (a[1] + b[1] * x + g[1] * y) / area2;
}

void compute_grad(int gridno, double *gx, double *gy)
{
    int i, j, n[3];
    double a, b, c, d;
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < 3; j++) {
	    n[j] = grid[gridno].icon[i].nl[j];
	}
	a = b = c = d = 0.0;
	for (j = 0; j < 3; j++) {
	    a += grid[gridno].yord[n[j]] * (grid[gridno].depth[n[(j + 1) % 3]]
					    - grid[gridno].depth[n[(j + 2) % 3]]);
	    b += grid[gridno].depth[n[j]] * (grid[gridno].xord[n[(j + 1) % 3]] - grid[gridno].xord[n[(j + 2) % 3]]);
	    c += grid[gridno].xord[n[j]] * (grid[gridno].yord[n[(j + 1) % 3]] - grid[gridno].yord[n[(j + 2) % 3]]);
	}
	d = -a * grid[gridno].xord[n[0]] - b * grid[gridno].yord[n[0]] - c * grid[gridno].depth[n[0]];
	if (c != 0.0) {
	    gx[i] = a / c;
	    gy[i] = b / c;
	} else {
	    printf("c = 0.0\n");
	}
    }
}
