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
 * dgrid.c - a routine to draw a grid
 *
 * Contents:
 *
 * drawgrid(int gridno, double red) ..... draw the grid designated by gridno
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawgrid.c,v 1.9 2006/11/12 05:09:12 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"
#include "externs.h"

static double dxv, dyv, dxi, dyi, scalex, scaley;
static double xg2s, yg2s;

double cour_dt, cour1_level, cour2_level;

double dimw_dt, dimw1_level, dimw2_level;
int display_dimw, display_dimwn;

double grad_scale, grad_units;

void set_scale(int gno);

void set_defaults(int gno)
{
    int i, didone = 0;
    double xmin = 1e307, xmax = -1e307, ymin = 1e307, ymax = -1e307;
    Grid *gd;

    for (i = 0; i < MAXGRIDS; i++) {
	if (object_isactive(GRID, i)) {
	    didone = 1;
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
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, ELEV, i) || object_isactive(ADCIRC, FLOW, i)) {
	    didone = 2;
	    gd = &flowt[i].g;
	    if (gd->xmin < xmin) {
		xmin = gd->xmin;
	    }
	    if (gd->xmax > xmax) {
		xmax = gd->xmax;
	    }
	    if (gd->ymin < ymin) {
		ymin = gd->ymin;
	    }
	    if (gd->ymax > ymax) {
		ymax = gd->ymax;
	    }
	}
    }
    for (i = 0; i < MAXTRANS; i++) {
	if (object_isactive(TRANSECT, i)) {
	    didone = 3;
	    gd = &(trans[i].g[0]);
	    if (gd == NULL)
		break;
	    if (gd->xmin < xmin) {
		xmin = gd->xmin;
	    }
	    if (gd->xmax > xmax) {
		xmax = gd->xmax;
	    }
	    if (gd->ymin < ymin) {
		ymin = gd->ymin;
	    }
	    if (gd->ymax > ymax) {
		ymax = gd->ymax;
	    }
	}
    }
    if (didone) {
	g[gno].w.xg1 = xmin - 0.025 * (xmax - xmin);
	g[gno].w.yg1 = ymin - 0.025 * (ymax - ymin);
	g[gno].w.xg2 = xmax + 0.025 * (xmax - xmin);
	g[gno].w.yg2 = ymax + 0.025 * (ymax - ymin);
    } else {
	xmin = ymin = -10.0;
	xmax = ymax = 10.0;
	g[gno].w.xg1 = xmin - 0.025 * (xmax - xmin);
	g[gno].w.yg1 = ymin - 0.025 * (ymax - ymin);
	g[gno].w.xg2 = xmax + 0.025 * (xmax - xmin);
	g[gno].w.yg2 = ymax + 0.025 * (ymax - ymin);
    }
}

void set_defaults_grid(int gno, int gridno)
{
    int i;
    double xmin = 1e307, xmax = -1e307, ymin = 1e307, ymax = -1e307;

    if (object_isactive(GRID, gridno)) {
	if (grid[gridno].xmin < xmin) {
	    xmin = grid[gridno].xmin;
	}
	if (grid[gridno].xmax > xmax) {
	    xmax = grid[gridno].xmax;
	}
	if (grid[gridno].ymin < ymin) {
	    ymin = grid[gridno].ymin;
	}
	if (grid[gridno].ymax > ymax) {
	    ymax = grid[gridno].ymax;
	}
	g[gno].w.xg1 = xmin - 0.025 * (xmax - xmin);
	g[gno].w.yg1 = ymin - 0.025 * (ymax - ymin);
	g[gno].w.xg2 = xmax + 0.025 * (xmax - xmin);
	g[gno].w.yg2 = ymax + 0.025 * (ymax - ymin);
    }
}

/*
	define world coordinate system and initialize
	the graphics drivers
*/
void set_up_world(int gno)
{
    set_scale(gno);
    defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, 0, 0);
    viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
}

void set_scale(int gno)
{
    extern int geoflag;
    if (debuglevel == 20) {
	/* printf("Before %lf %lf %lf %lf  %lf %lf %lf %lf\n", g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2, g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2); */
    }
    if (g[gno].type == XYFIXED) {
	setfixedscale(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2, &g[gno].w.xg1, &g[gno].w.yg1, &g[gno].w.xg2, &g[gno].w.yg2);
    }
    if (debuglevel == 20) {
// backward compatibility
	printf("%lf %lf %lf %lf\n", g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2);
    }
    if (geoflag) {
// make the .gfw file for georeferencing
	extern char geoprintstr[];
	FILE *gp;
	double px, py, dx, dy;
	dx = g[gno].w.xg2 - g[gno].w.xg1;
	dy = g[gno].w.yg2 - g[gno].w.yg1;
	px = dx / devwidth;
	py = -dy / devheight;
	if ((gp = fopen(geoprintstr, "w")) != NULL) {
	    fprintf(gp, "%.8lf\n0\n0\n%.8lf\n%.8lf\n%.8lf\n", px, py, g[gno].w.xg1, g[gno].w.yg2);
	    fclose(gp);
	} else {
	    fprintf(stderr, "Unable to open .gfw file = [%s]\n", geoprintstr);
	}
    }
}

void set_up_mapscale(int gno, double mapscale)
{
    setmapscale(mapscale, g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2, &g[gno].w.xg1, &g[gno].w.yg1, &g[gno].w.xg2, &g[gno].w.yg2);
}

void reset_world(int gno)
{
/*
    autoscalegrid(gno, 0);
*/
    set_defaults(gno);
    set_scale(gno);
    defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, 0, 0);
    viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
}

int elem_ok(Grid * g, int elno)
{
    double x1, y1, x2, y2, x3, y3;
    int n1, n2, n3;
    n1 = g->icon[elno].nl[0];
    n2 = g->icon[elno].nl[1];
    n3 = g->icon[elno].nl[2];
    x1 = g->xord[n1];
    y1 = g->yord[n1];
    x2 = g->xord[n2];
    y2 = g->yord[n2];
    x3 = g->xord[n3];
    y3 = g->yord[n3];
    if (symok(x1, y1) || symok(x2, y2) || symok(x3, y3)) {
	return 1;
    } else {
	return 0;
    }
}

void drawgrid(int gno, int gridno, double red)
{
    int i, n0, n1, n2, n3, redflag = (red == 1.0);
    double xg, yg, x, y;
    setcolor(g[gno].grid[gridno].p.color);
    setlinestyle(g[gno].grid[gridno].p.lines);
    setlinewidth(g[gno].grid[gridno].p.linew);
    for (i = 0; i < grid[gridno].nmel; i++) {
	switch (grid[gridno].icon[i].type) {
	case 0:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n3], grid[gridno].yord[n3]);
	    break;
	case 1:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    break;
	case 4:
	case 8:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_draw2(grid[gridno].xord[n3], grid[gridno].yord[n3]);
	    my_draw2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    break;
	case 5:
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[3];
	    my_move2(grid[gridno].xord[n0], grid[gridno].yord[n0]);
	    my_draw2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n2], grid[gridno].yord[n2]);
	    my_move2(grid[gridno].xord[n1], grid[gridno].yord[n1]);
	    my_draw2(grid[gridno].xord[n3], grid[gridno].yord[n3]);
	    break;
	case 3:
	case 6:
	    if (!elem_ok(&grid[gridno], i)) {
		break;
	    }
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
	}
    }
}

void drawgrid_filled(int gridno, int pat)
{
    int i, n0, n1, n2, n3;
    double xg, yg, x[4], y[4];
    setcolor(pat);
    for (i = 0; i < grid[gridno].nmel; i++) {
	switch (grid[gridno].icon[i].type) {
	case 3:
	case 6:
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
	case 8:
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
    }
    setcolor(1);
}

void drawelement_filled(int gridno, int elno, int pat)
{
    int n0, n1, n2;
    double x[4], y[4];
    setcolor(pat);
    n0 = grid[gridno].icon[elno].nl[0];
    n1 = grid[gridno].icon[elno].nl[1];
    n2 = grid[gridno].icon[elno].nl[2];
    x[0] = grid[gridno].xord[n0];
    y[0] = grid[gridno].yord[n0];
    x[1] = grid[gridno].xord[n1];
    y[1] = grid[gridno].yord[n1];
    x[2] = grid[gridno].xord[n2];
    y[2] = grid[gridno].yord[n2];
    fillcolor(3, x, y);
}

void drawgrad(int gridno)
{
    int i, n0, n1, n2;
    extern int mapisolconc[];
    double q, *xg, *yg, x[4], y[4], xc, yc;
    xg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    yg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    compute_grad(gridno, xg, yg);
    setcolor(1);
    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xc, &yc);
	velplt(xc, yc, xg[i], yg[i], grad_units * grad_scale);
    }
    free(xg);
    free(yg);
}

void drawslopes(int gridno)
{
    int i, n0, n1, n2;
    extern int mapisolconc[];
    extern double grad_scale, grad_units;
    double q, *xg, *yg, x[4], y[4], xc, yc;
    xg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    yg = (double *) calloc(grid[gridno].nmel, sizeof(double));
    compute_grad(gridno, xg, yg);
    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xc, &yc);
	velplt(xc, yc, xg[i], yg[i], grad_units * grad_scale);
    }
    free(xg);
    free(yg);
}


/* draw color filled elements with courant number of the courant number
 * values themselves
 */
void drawgrid_cour_filled(int gno, int gridno)
{
    int i, k, n0, n1, n2;
    double xg, yg, x[4], y[4], a, d, cu;
    double area(int gridno, int elem), dt;
    char buf[50];
    dt = cour_dt;
    if (g[gno].grid[gridno].display_courant == ON || g[gno].grid[gridno].display_courantn == ON) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n0 = grid[gridno].icon[i].nl[0];
	    n1 = grid[gridno].icon[i].nl[1];
	    n2 = grid[gridno].icon[i].nl[2];
	    d = grid[gridno].depth[n0] + grid[gridno].depth[n1] + grid[gridno].depth[n2];
	    d = d * 0.3333333333333333;
	    a = area(gridno, i);
	    cu = sqrt(M_PI * 9.8 * d / (4.0 * a)) * dt;
	    if (g[gno].grid[gridno].display_courant == ON) {
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
	    if (g[gno].grid[gridno].display_courantn == ON) {
		setcolor(1);
		get_center(gridno, i, &xg, &yg);
		sprintf(buf, "%.2lf", cu);
		writestr(xg, yg, 0, 0, buf);
	    }
	}
	setcolor(1);
    }
}

void drawgrid_dimw_filled(int gridno)
{
    int i, k, n0, n1, n2;
    double xg, yg, x[4], y[4], a, d, dw;
    double area(int gridno, int elem), dt, dwmin = 0.0, dwmax = 10.0;
    char buf[50];
    dt = dimw_dt;
    if (display_dimw || display_dimwn) {
	for (i = 0; i < grid[gridno].nmel; i++) {
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
    int i, j, n0, n1, n2, redflag = (red == 1.0);
    double xg, yg, x, y;

    i = nel;
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

void drawgridnodes(int gridno)
{
    int i;
    char buf[80];

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
	    sprintf(buf, "%d", i + 1);
	    writestr(grid[gridno].xord[i], grid[gridno].yord[i], 0, 0, buf);
	}
    }
}

void drawgridelems(int gridno)
{
    int i;
    double xg, yg;
    char buf[80];

    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xg, &yg);
	if (symok(xg, yg)) {
	    sprintf(buf, "%d", i + 1);
	    writestr(xg, yg, 0, 0, buf);
	}
    }
}

void drawnodedepths(int gridno)
{
    int i;
    char buf[80];

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
	    sprintf(buf, "[%d,%.2lf]", i + 1, grid[gridno].depth[i]);
	    writestr(grid[gridno].xord[i], grid[gridno].yord[i], 0, 0, buf);
	}
    }
}

void getlimits_grid(int gridno, double *xmin, double *xmax, double *ymin, double *ymax, double *dmin, double *dmax)
{
    int i;

    if (grid[gridno].active == ON) {
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

    if (grid[gridno].active == ON) {
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

void DrawGrid(int gno, Grid * g, double red)
{
    int i, n0, n1, n2, n3, redflag = (red == 1.0);
    double xg, yg, x, y;
    setcolor(1);
    setlinestyle(1);
    setlinewidth(1);
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
