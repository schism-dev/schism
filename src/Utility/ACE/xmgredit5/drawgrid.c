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
 * drawgrid.c - Draw grids
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawgrid.c,v 1.3 2003/07/29 22:55:32 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

extern int cancel_flag;
extern int restrict_to_region;

extern int drawquadflag;
extern int drawtriflag;
extern int quadcolor;
extern int tricolor;

/*
 * Draw a grid with a reduction factor
 */
void drawgrid(int gridno, double red)
{
    int i, n0, n1, n2, n3, n4, n5, n6, n7, n8, redflag = (red == 1.0);
    double xg, yg, x, y;
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].fdgrid && grid[gridno].icon[i].wetdry) {
	    continue;
	}
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	if (restrict_to_region) {
	    get_center(gridno, i, &xg, &yg);
	    if (region_flag && !inregion(regionx, regiony, nregion, xg, yg)) {
		continue;
	    }
	}
	switch (grid[gridno].icon[i].type) {
	case 4:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    if (redflag) {
		my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
		my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
		my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
		my_draw2(grid[gridno].xord[n3], grid[gridno].yord[n3]);
		my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    } else {
		get_center(gridno, i, &xg, &yg);
		x = red * (grid[gridno].xord[n0] - xg) + xg;
		y = red * (grid[gridno].yord[n0] - yg) + yg;
		my_move2(x, y);
		x = red * (grid[gridno].xord[n1] - xg) + xg;
		y = red * (grid[gridno].yord[n1] - yg) + yg;
		my_draw2(x, y);
		x = red * (grid[gridno].xord[n2] - xg) + xg;
		y = red * (grid[gridno].yord[n2] - yg) + yg;
		my_draw2(x, y);
		x = red * (grid[gridno].xord[n3] - xg) + xg;
		y = red * (grid[gridno].yord[n3] - yg) + yg;
		my_draw2(x, y);
		x = red * (grid[gridno].xord[n0] - xg) + xg;
		y = red * (grid[gridno].yord[n0] - yg) + yg;
		my_draw2(x, y);
	    }
	    break;
	case 3:
	case 6:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    if (redflag) {
		my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
		my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
		my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
		my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    } else {
		get_center(gridno, i, &xg, &yg);
		x = red * (grid[gridno].xord[n0] - xg) + xg;
		y = red * (grid[gridno].yord[n0] - yg) + yg;
		my_move2(x, y);
		x = red * (grid[gridno].xord[n1] - xg) + xg;
		y = red * (grid[gridno].yord[n1] - yg) + yg;
		my_draw2(x, y);
		x = red * (grid[gridno].xord[n2] - xg) + xg;
		y = red * (grid[gridno].yord[n2] - yg) + yg;
		my_draw2(x, y);
		x = red * (grid[gridno].xord[n0] - xg) + xg;
		y = red * (grid[gridno].yord[n0] - yg) + yg;
		my_draw2(x, y);
	    }
	    break;
	default:
	    printf("Element type unknown %d\n", grid[gridno].icon[i].type);
	}
    }
}

/*
 * Fill each element of a grid with color color
 */
void drawgrid_filled(int gridno, int color)
{
    int i, n0, n1, n2, n3;
    double xg, yg, x[4], y[4];
    setcolor(color);
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	switch (grid[gridno].icon[i].nn) {
	case 3:
	case 6:
	    if (drawtriflag) {
		setcolor(tricolor);
	    }
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    x[0] = grid[gridno].xord[n0];
	    y[0] = grid[gridno].yord[n0];
	    x[1] = grid[gridno].xord[n1];
	    y[1] = grid[gridno].yord[n1];
	    x[2] = grid[gridno].xord[n2];
	    y[2] = grid[gridno].yord[n2];
	    fillcolor(3, x, y);
	    break;
	case 4:
	    if (drawquadflag) {
		setcolor(quadcolor);
	    }
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    x[0] = grid[gridno].xord[n0];
	    y[0] = grid[gridno].yord[n0];
	    x[1] = grid[gridno].xord[n1];
	    y[1] = grid[gridno].yord[n1];
	    x[2] = grid[gridno].xord[n2];
	    y[2] = grid[gridno].yord[n2];
	    x[3] = grid[gridno].xord[n3];
	    y[3] = grid[gridno].yord[n3];
	    fillcolor(4, x, y);
	    break;
	}
        setcolor(color);
    }
    setcolor(1);
}

/*
 * Fill a single element with color color
 */
void fillelement(int gridno, int el, int color)
{
    int n0, n1, n2, n3;
    double x[4], y[4];
    switch (grid[gridno].icon[el].type) {
    case 3:
    case 6:
	n0 = grid[gridno].icon[el].nl[0];
	n1 = grid[gridno].icon[el].nl[1];
	n2 = grid[gridno].icon[el].nl[2];
	x[0] = grid[gridno].xord[n0];
	y[0] = grid[gridno].yord[n0];
	x[1] = grid[gridno].xord[n1];
	y[1] = grid[gridno].yord[n1];
	x[2] = grid[gridno].xord[n2];
	y[2] = grid[gridno].yord[n2];
	fillcolor(3, x, y);
	break;
    case 4:
	n0 = grid[gridno].icon[el].nl[0];
	n1 = grid[gridno].icon[el].nl[1];
	n2 = grid[gridno].icon[el].nl[2];
	n3 = grid[gridno].icon[el].nl[3];
	x[0] = grid[gridno].xord[n0];
	y[0] = grid[gridno].yord[n0];
	x[1] = grid[gridno].xord[n1];
	y[1] = grid[gridno].yord[n1];
	x[2] = grid[gridno].xord[n2];
	y[2] = grid[gridno].yord[n2];
	x[3] = grid[gridno].xord[n3];
	y[3] = grid[gridno].yord[n3];
	fillcolor(4, x, y);
	break;
    }
}

/*
 * Draw the gradient of each element as a vector at the cell center.
 */
void drawgrad(int gridno)
{
    FILE *fp;
    int i, n0, n1, n2;
    extern int mapisolconc[];
    extern double grad_scale, grad_units;
    double q, *xg, *yg, x[4], y[4], xc, yc;
    fp = fopen("grad.d", "w");
    xg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    yg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    compute_grad(gridno, xg, yg);
    setcolor(1);
    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xc, &yc);
	velplt(xc, yc, xg[i], yg[i], grad_units * grad_scale);
	fprintf(fp, "%lf %lf\n", xg[i], yg[i]);
    }
    fclose(fp);
    free(xg);
    free(yg);
}

/*
 * Draw the slopes of a grid as isolines
 */
void drawslopes(int gridno)
{
    int i, n0, n1, n2;
    extern int mapisolconc[];
    extern double grad_scale, grad_units, slope_cutoff;
    double q, *xg, *yg, x[4], y[4], xc, yc;
    double slope;
    xg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    yg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    compute_grad(gridno, xg, yg);
    setcolor(2);
    for (i = 0; i < grid[gridno].nmel; i++) {
	slope = hypot(xg[i], yg[i]);
	if (slope > slope_cutoff) {
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    x[0] = grid[gridno].xord[n0];
	    y[0] = grid[gridno].yord[n0];
	    x[1] = grid[gridno].xord[n1];
	    y[1] = grid[gridno].yord[n1];
	    x[2] = grid[gridno].xord[n2];
	    y[2] = grid[gridno].yord[n2];
	    fillcolor(3, x, y);
	}
    }
    setcolor(1);
    free(xg);
    free(yg);
}


/* 
 * Draw grid color filled based on Courant number or dimensionless
 * wavelength.
 */
extern double cour_dt, cour1_level, cour2_level;
extern int display_cour, display_courn;

void drawgrid_cour_filled(int gridno)
{
    int i, k, n0, n1, n2;
    double xg, yg, x[4], y[4], a, d, cu;
    double area(), dt, cumin = 0.0, cumax = 10.0;
    char buf[50];
    extern double cour_dt, cour1_level, cour2_level;
    dt = cour_dt;
    if (display_cour || display_courn) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    if (!(i % 100) && check_action()) {
		cancel_action(0);
		if (cancel_flag)
		    return;
	    }
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    d = grid[gridno].depth[n0] + grid[gridno].depth[n1] + grid[gridno].depth[n2];
	    d = d * 0.3333333333333333;
	    a = area(gridno, i);
	    cu = sqrt(M_PI * 9.8 * d / (4.0 * a)) * dt;
	    if (display_cour) {
		if (cu < cour1_level) {
		    setcolor(3);
		} else if (cu >= cour1_level && cu < cour2_level) {
		    setcolor(11);
		} else if (cu >= cour2_level) {
		    setcolor(2);
		}
		x[0] = grid[gridno].xord[n0];
		y[0] = grid[gridno].yord[n0];
		x[1] = grid[gridno].xord[n1];
		y[1] = grid[gridno].yord[n1];
		x[2] = grid[gridno].xord[n2];
		y[2] = grid[gridno].yord[n2];
		fillcolor(3, x, y);
	    }
	    if (display_courn) {
		setcolor(1);
		get_center(gridno, i, &xg, &yg);
		sprintf(buf, "%.2lf", cu);
		writestr(xg, yg, 0, 0, buf);
	    }
	}
	setcolor(1);
    }
}

extern double dimw_dt, dimw1_level, dimw2_level;
extern int display_dimw, display_dimwn;

void drawgrid_dimw_filled(int gridno)
{
    int i, k, n0, n1, n2;
    double xg, yg, x[4], y[4], a, d, dw;
    double area(), dt, dwmin = 0.0, dwmax = 10.0;
    char buf[50];
    extern double dimw_dt, dimw1_level, dimw2_level;
    dt = dimw_dt;
    if (display_dimw || display_dimwn) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    if (!(i % 100) && check_action()) {
		cancel_action(0);
		if (cancel_flag)
		    return;
	    }
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    d = grid[gridno].depth[n0] + grid[gridno].depth[n1] + grid[gridno].depth[n2];
	    d = d * 0.3333333333333333;
	    a = area(gridno, i);
	    dw = sqrt(M_PI * 9.8 * d / (4.0 * a)) * dt;
	    if (display_dimw) {
		if (dw > dimw1_level) {
		    setcolor(3);
		} else if (dw <= dimw1_level && dw > dimw2_level) {
		    setcolor(11);
		} else if (dw <= dimw2_level) {
		    setcolor(2);
		}
		x[0] = grid[gridno].xord[n0];
		y[0] = grid[gridno].yord[n0];
		x[1] = grid[gridno].xord[n1];
		y[1] = grid[gridno].yord[n1];
		x[2] = grid[gridno].xord[n2];
		y[2] = grid[gridno].yord[n2];
		fillcolor(3, x, y);
	    }
	    if (display_dimwn) {
		setcolor(1);
		get_center(gridno, i, &xg, &yg);
		sprintf(buf, "%.2lf", dw);
		writestr(xg, yg, 0, 0, buf);
	    }
	}
	setcolor(1);
    }
}

void draw_elements(int gridno, int nels, int *nelbuf, double red)
{
    int i, j, n0, n1, n2, redflag = (red == 1.0);
    double xg, yg, x, y;

    for (j = 0; j < nels; j++) {
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	i = nelbuf[j];
	n0 = grid[gridno].icon[i].nl[0];
	n1 = grid[gridno].icon[i].nl[1];
	n2 = grid[gridno].icon[i].nl[2];
	if (redflag) {
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	} else {
	    get_center(gridno, i, &xg, &yg);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_move2(x, y);
	    x = red * (grid[gridno].xord[n1] - xg) + xg;
	    y = red * (grid[gridno].yord[n1] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n2] - xg) + xg;
	    y = red * (grid[gridno].yord[n2] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_draw2(x, y);
	}
    }
}

void draw_element(int gridno, int nel, double red)
{
    int i, j, n0, n1, n2, n3, n4, n5, n6, n7, redflag = (red == 1.0);
    double xg, yg, x, y;

    i = nel;
    switch (grid[gridno].icon[i].type) {
    case 4:
	n0 = grid[gridno].icon[i].nl[0];
	n1 = grid[gridno].icon[i].nl[1];
	n2 = grid[gridno].icon[i].nl[2];
	n3 = grid[gridno].icon[i].nl[3];
	if (redflag) {
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n3], grid[gridno].yord[n3]);
	    my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	} else {
	    get_center(gridno, i, &xg, &yg);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_move2(x, y);
	    x = red * (grid[gridno].xord[n1] - xg) + xg;
	    y = red * (grid[gridno].yord[n1] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n2] - xg) + xg;
	    y = red * (grid[gridno].yord[n2] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n3] - xg) + xg;
	    y = red * (grid[gridno].yord[n3] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_draw2(x, y);
	}
	break;
    case 3:
    case 6:
	n0 = grid[gridno].icon[i].nl[0];
	n1 = grid[gridno].icon[i].nl[1];
	n2 = grid[gridno].icon[i].nl[2];
	if (redflag) {
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	} else {
	    get_center(gridno, i, &xg, &yg);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_move2(x, y);
	    x = red * (grid[gridno].xord[n1] - xg) + xg;
	    y = red * (grid[gridno].yord[n1] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n2] - xg) + xg;
	    y = red * (grid[gridno].yord[n2] - yg) + yg;
	    my_draw2(x, y);
	    x = red * (grid[gridno].xord[n0] - xg) + xg;
	    y = red * (grid[gridno].yord[n0] - yg) + yg;
	    my_draw2(x, y);
	}
	break;
    default:
	printf("Element type unknown %d\n", grid[gridno].icon[i].type);
    }
}

void draw_element_filled(int gridno, int nel, double red)
{
    int i, j, n0, n1, n2, n3, n4, n5, n6, n7, n8, n9, redflag = (red == 1.0);
    double xg, yg, xt, yt, x[10], y[10];

    i = nel;
    switch (grid[gridno].icon[i].type) {
    case 4:
	n0 = grid[gridno].icon[i].nl[0];
	n1 = grid[gridno].icon[i].nl[1];
	n2 = grid[gridno].icon[i].nl[2];
	n3 = grid[gridno].icon[i].nl[3];
	if (redflag) {
	    x[0] = grid[gridno].xord[n0];
	    y[0] = grid[gridno].yord[n0];
	    x[1] = grid[gridno].xord[n1];
	    y[1] = grid[gridno].yord[n1];
	    x[2] = grid[gridno].xord[n2];
	    y[2] = grid[gridno].yord[n2];
	    x[3] = grid[gridno].xord[n3];
	    y[3] = grid[gridno].yord[n3];
	} else {
	    get_center(gridno, i, &xg, &yg);
	    x[0] = red * (grid[gridno].xord[n0] - xg) + xg;
	    y[0] = red * (grid[gridno].yord[n0] - yg) + yg;
	    x[1] = red * (grid[gridno].xord[n1] - xg) + xg;
	    y[1] = red * (grid[gridno].yord[n1] - yg) + yg;
	    x[2] = red * (grid[gridno].xord[n2] - xg) + xg;
	    y[2] = red * (grid[gridno].yord[n2] - yg) + yg;
	    x[3] = red * (grid[gridno].xord[n3] - xg) + xg;
	    y[3] = red * (grid[gridno].yord[n3] - yg) + yg;
	    x[4] = red * (grid[gridno].xord[n0] - xg) + xg;
	    y[4] = red * (grid[gridno].yord[n0] - yg) + yg;
	}
	fillcolor(4, x, y);
	break;
    case 3:
    case 6:
	n0 = grid[gridno].icon[i].nl[0];
	n1 = grid[gridno].icon[i].nl[1];
	n2 = grid[gridno].icon[i].nl[2];
	if (redflag) {
	    x[0] = grid[gridno].xord[n0];
	    y[0] = grid[gridno].yord[n0];
	    x[1] = grid[gridno].xord[n1];
	    y[1] = grid[gridno].yord[n1];
	    x[2] = grid[gridno].xord[n2];
	    y[2] = grid[gridno].yord[n2];
	} else {
	    get_center(gridno, i, &xg, &yg);
	    x[0] = red * (grid[gridno].xord[n0] - xg) + xg;
	    y[0] = red * (grid[gridno].yord[n0] - yg) + yg;
	    x[1] = red * (grid[gridno].xord[n1] - xg) + xg;
	    y[1] = red * (grid[gridno].yord[n1] - yg) + yg;
	    x[2] = red * (grid[gridno].xord[n2] - xg) + xg;
	    y[2] = red * (grid[gridno].yord[n2] - yg) + yg;
	}
	fillcolor(3, x, y);
	break;
    default:
	printf("Element type unknown %d\n", grid[gridno].icon[i].type);
    }
}

void draw_nodes(int n0, int n1, int n2)
{
    setcolor(4);
    my_move2(build[0].bx[n0], build[0].by[n0]);
    my_draw2(build[0].bx[n1], build[0].by[n1]);
    my_draw2(build[0].bx[n2], build[0].by[n2]);
    my_draw2(build[0].bx[n0], build[0].by[n0]);
}

/*
 * draw either the node numbers n = 0 or mark the nodes n = 1
 */
void drawgridnodes(int gridno, int l, int n)
{
    int i;
    char buf[80];

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	if (n) {
	    if (symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
		my_circle(grid[gridno].xord[i], grid[gridno].yord[i]);
	    }
	}
	if (l) {
	    if (symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
		if (n) {
		    sprintf(buf, " %d", i + 1);
		} else {
		    sprintf(buf, "%d", i + 1);
		}
		setcharsize(1.2);
		writestr(grid[gridno].xord[i], grid[gridno].yord[i], 0, 0, buf);
		setcharsize(1.0);
	    }
	}
    }
}

void drawgridelems(int gridno)
{
    int i;
    double xg, yg;
    char buf[80];

    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	get_center(gridno, i, &xg, &yg);
	if (symok(xg, yg)) {
	    sprintf(buf, " %d", i + 1);
	    writestr(xg, yg, 0, 0, buf);
	}
    }
}

void drawnodedepths(int gridno)
{
    int i;
    char buf[80];

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (!(i % 100) && check_action()) {
	    cancel_action(0);
	    if (cancel_flag)
		return;
	}
	if (symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
	    sprintf(buf, "[%d,%.2lf]", i + 1, grid[gridno].depth[i]);
	    writestr(grid[gridno].xord[i], grid[gridno].yord[i], 0, 0, buf);
	}
    }
}

void set_defaults(void)
{
    int i;
    double xmin = 1e307, xmax = -1e307, ymin = 1e307, ymax = -1e307;

    for (i = 0; i < MAXGRIDS; i++) {
	if (grid[i].gactive) {
	    if (grid[i].xmin < xmin) {
		xmin = grid[i].xmin;
	    }
	    if (grid[i].xmax > xmax) {
		xmax = grid[i].xmax;
	    }
	    if (grid[i].ymin < ymin) {
		ymin = grid[i].ymin;
	    }
	    if (grid[i].ymax > ymax) {
		ymax = grid[i].ymax;
	    }
	}
    }
    xg1 = sxg1 = xmin - 0.025 * (xmax - xmin);
    yg1 = syg1 = ymin - 0.025 * (ymax - ymin);
    xg2 = sxg2 = xmax + 0.025 * (xmax - xmin);
    yg2 = syg2 = ymax + 0.025 * (ymax - ymin);
}

void set_default_scale(double xmin, double xmax, double ymin, double ymax)
{
    int i;
    xg1 = sxg1 = xmin - 0.025 * (xmax - xmin);
    yg1 = syg1 = ymin - 0.025 * (ymax - ymin);
    xg2 = sxg2 = xmax + 0.025 * (xmax - xmin);
    yg2 = syg2 = ymax + 0.025 * (ymax - ymin);
}

void getlimits_grid(int gridno, double *xmin, double *xmax, double *ymin, double *ymax, double *dmin, double *dmax)
{
    int i;

    if (grid[gridno].gactive) {
	*xmin = *xmax = grid[gridno].xord[0];
	*ymin = *ymax = grid[gridno].yord[0];
	*dmin = *dmax = grid[gridno].depth[0];
	for (i = 1; i < grid[gridno].nmnp; i++) {
	    if (grid[gridno].xord[i] < *xmin) {
		*xmin = grid[gridno].xord[i];
	    }
	    if (grid[gridno].xord[i] > *xmax) {
		*xmax = grid[gridno].xord[i];
	    }
	    if (grid[gridno].yord[i] < *ymin) {
		*ymin = grid[gridno].yord[i];
	    }
	    if (grid[gridno].yord[i] > *ymax) {
		*ymax = grid[gridno].yord[i];
	    }
	    if (grid[gridno].depth[i] < *dmin) {
		*dmin = grid[gridno].depth[i];
	    }
	    if (grid[gridno].depth[i] > *dmax) {
		*dmax = grid[gridno].depth[i];
	    }
	}
	grid[gridno].xmin = *xmin;
	grid[gridno].xmax = *xmax;
	grid[gridno].ymin = *ymin;
	grid[gridno].ymax = *ymax;
	grid[gridno].dmin = *dmin;
	grid[gridno].dmax = *dmax;
    }
}

void setlimits_grid(int gridno)
{
    int i;
    double xmin, xmax, ymin, ymax, dmin, dmax;

    if (grid[gridno].gactive) {
	xmin = xmax = grid[gridno].xord[0];
	ymin = ymax = grid[gridno].yord[0];
	dmin = dmax = grid[gridno].depth[0];
	for (i = 1; i < grid[gridno].nmnp; i++) {
	    if (grid[gridno].xord[i] < xmin) {
		xmin = grid[gridno].xord[i];
	    }
	    if (grid[gridno].xord[i] > xmax) {
		xmax = grid[gridno].xord[i];
	    }
	    if (grid[gridno].yord[i] < ymin) {
		ymin = grid[gridno].yord[i];
	    }
	    if (grid[gridno].yord[i] > ymax) {
		ymax = grid[gridno].yord[i];
	    }
	    if (grid[gridno].depth[i] < dmin) {
		dmin = grid[gridno].depth[i];
	    }
	    if (grid[gridno].depth[i] > dmax) {
		dmax = grid[gridno].depth[i];
	    }
	}
	grid[gridno].xmin = xmin;
	grid[gridno].xmax = xmax;
	grid[gridno].ymin = ymin;
	grid[gridno].ymax = ymax;
	grid[gridno].dmin = dmin;
	grid[gridno].dmax = dmax;
    }
}

void getlimits_boundary(int ib, double *xmin, double *xmax, double *ymin, double *ymax)
{
    int j;

    *xmin = *xmax = boundary[ib].boundx[0];
    *ymin = *ymax = boundary[ib].boundy[0];
    for (j = 1; j < boundary[ib].nbpts; j++) {
	if (boundary[ib].boundx[j] < *xmin) {
	    *xmin = boundary[ib].boundx[j];
	}
	if (boundary[ib].boundx[j] > *xmax) {
	    *xmax = boundary[ib].boundx[j];
	}
	if (boundary[ib].boundy[j] < *ymin) {
	    *ymin = boundary[ib].boundy[j];
	}
	if (boundary[ib].boundy[j] > *ymax) {
	    *ymax = boundary[ib].boundy[j];
	}
    }
}

void getlimits_build(int bno, double *xmin, double *xmax, double *ymin, double *ymax, double *dmin, double *dmax)
{
    int j;

    *xmin = *xmax = build[bno].bx[0];
    *ymin = *ymax = build[bno].by[0];
    *dmin = *dmax = build[bno].db[0];
    for (j = 1; j < build[bno].nbuild; j++) {
	if (build[bno].bx[j] < *xmin) {
	    *xmin = build[bno].bx[j];
	}
	if (build[bno].bx[j] > *xmax) {
	    *xmax = build[bno].bx[j];
	}
	if (build[bno].by[j] < *ymin) {
	    *ymin = build[bno].by[j];
	}
	if (build[bno].by[j] > *ymax) {
	    *ymax = build[bno].by[j];
	}
	if (build[bno].db[j] < *dmin) {
	    *dmin = build[bno].db[j];
	}
	if (build[bno].db[j] > *dmax) {
	    *dmax = build[bno].db[j];
	}
    }
}

void getlimits(double *x, int n, double *xmin, double *xmax)
{
    int j;

    *xmin = *xmax = x[0];
    for (j = 1; j < n; j++) {
	if (x[j] < *xmin) {
	    *xmin = x[j];
	}
	if (x[j] > *xmax) {
	    *xmax = x[j];
	}
    }
}

void do_showelems_region(void)
{
    int i;
    double xg, yg;
    setcolor(1);
    for (i = 0; i < grid[0].nmel; i++) {
	get_center(0, i, &xg, &yg);
	if (region_flag && inregion(regionx, regiony, nregion, xg, yg)) {
	    if (symok(xg, yg)) {
		solidbox(xg, yg);
/*
   sprintf(buf, " %d", i + 1);
   writestr(xg, yg, 0, 0, buf);
 */
	    }
	}
    }
}

void do_shownodes_region(void)
{
    int i;
    char buf[256];
    double x, y;
    setcolor(1);
    for (i = 0; i < grid[0].nmnp; i++) {
	x = grid[0].xord[i];
	y = grid[0].yord[i];
	if (region_flag && inregion(regionx, regiony, nregion, x, y)) {
	    if (symok(x, y)) {
		solidbox(x, y);
/*
   sprintf(buf, " %d", i + 1);
   writestr(x, y, 0, 0, buf);
 */
	    }
	}
    }
}

void DrawGrid(Grid * g)
{
    int i, n0, n1, n2;
    setcolor(1);
    for (i = 0; i < g->nmel; i++) {
	n0 = g->icon[i].nl[0];
	n1 = g->icon[i].nl[1];
	n2 = g->icon[i].nl[2];
	my_move2(g->xord[n0], g->yord[n0]);
	my_draw2(g->xord[n1], g->yord[n1]);
	my_draw2(g->xord[n2], g->yord[n2]);
	my_draw2(g->xord[n0], g->yord[n0]);
    }
}

/*
 * Experimental code to draw the circumcenters of elements
 */
int draw_circum = 0;

void draw_circumcenters(int gridno)
{
    int i;
    int n1, n2, n3;
    double x1, y1, x2, y2, x3, y3, xm, ym, xg, yg;
    if (draw_circum) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    x1 = grid[gridno].xord[n1];
	    y1 = grid[gridno].yord[n1];
	    x2 = grid[gridno].xord[n2];
	    y2 = grid[gridno].yord[n2];
	    x3 = grid[gridno].xord[n3];
	    y3 = grid[gridno].yord[n3];
	    get_center(gridno, i, &xg, &yg);
	    compute_circumcenter(x1, y1, x2, y2, x3, y3, &xm, &ym);
	    my_move2(xg, yg);
	    my_draw2(xm, ym);
	    my_circlefilled(xm, ym);
	}
    }
}
