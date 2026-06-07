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
 * draw the grid and other items if on
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawobjs.c,v 1.6 2011/09/14 17:44:21 pturner Exp $";
#endif

#include <stdio.h>

#include "defines.h"
#include "globals.h"
#include "graphics.h"
#include "externs.h"

/* drawobjs.c */
void set_main_canvas(void);
void do_drawgrid(void);
void do_drawobjects(void);
void do_hardcopy(void);
void set_up_world(void);
void set_scale(void);
void set_up_mapscale(double mapscale);
void reset_world(void);
void draw_region(double *rx, double *ry, int nr);
void dolegend(Isolparms ip, int *cmap, int nmap);

extern double xv1, xv2, yv1, yv2;

extern int mapisolconc[], mapisolbath[];

extern int tdevice;

extern int pointset;
extern double dsx, dsy;
extern int win_h, win_w;
extern int doclear, hdevice, hardcopyflag, ptofile;
extern char printstr[];
extern int list_flag;

int auto_redraw = 1;
int forced_draw = 0;
int draw_filled = 0;
int draw_prop = 0;
int draw_isolfilled = 0;
extern int display_maxangle, display_nnels, display_grad, display_grad_legend, display_slope, display_delh;
int display_mapscale = 0;
extern int displaynl;
extern int display_shapiro;

void do_drawobjects(void);
void dolegend(Isolparms ip, int *cmap, int nmap);

int grid_color = 1;
int bound_color = 2;
int backgrid_color = 11;
int backbound_color = 4;

/* cancel drawing */
extern int cancel_flag;

/*
   set_abort_button();
   set_redraw_button();
 */

void set_main_canvas(void)
{
    set_canvas();
    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(xv1, yv1, xv2, yv2);
}

void do_drawgrid(void)
{
    if (auto_redraw || forced_draw) {
	write_mode_str("Re-drawing, press Stop to terminate drawing");
	forced_draw = 0;
	set_canvas();
	setup_abort();
	initgraphics(tdevice);
	defineworld(xg1, yg1, xg2, yg2, 0, 0);
	//printf("%lf %lf %lf %lf %d %d\n", xg1, yg1, xg2, yg2, devwidth, devheight);
	viewport(xv1, yv1, xv2, yv2);
	set_wait_cursor(NULL);
	drawimage();
	do_drawobjects();
	draw_annotation(0);
	draw_annotation(-1);
	leavegraphics();
	unset_wait_cursor(NULL);
	unsetup_abort();
	write_mode_str(NULL);
	reset_action();
    }
}

void do_drawobjects(void)
{
    int i;
    extern int mapisolbath[];
    cancel_flag = 0;
    setcolor(1);
    setlinewidth(1);
    setcharsize(0.6);
    if (display_flags[EDIT_GRID_ISOLINES]) {
	if (grid[curgrid].depthflag) {
	    do_isol2(curgrid, 0, grid[curgrid].edepth, grid[curgrid].ip, mapisolbath, grid[curgrid].ip.nisol);
	} else {
	    do_isol(curgrid, 0, grid[curgrid].depth, grid[curgrid].ip, mapisolbath, grid[curgrid].ip.nisol);
	}
	if (grid[curgrid].ip.lactive == ON) {
	    dolegend(grid[curgrid].ip, mapisolbath, grid[curgrid].ip.nisol);
	}
    }
    if (display_flags[BACKGROUND_GRID_ISOLINES]) {
	do_isol(backgrid, 0, grid[backgrid].depth, grid[backgrid].ip, mapisolbath, grid[backgrid].ip.nisol);
	if (grid[backgrid].ip.lactive == ON) {
	    dolegend(grid[backgrid].ip, mapisolbath, grid[backgrid].ip.nisol);
	}
    }
    if (display_flags[EDIT_GRID_FILLED]) {
	drawgrid_filled(curgrid, 7);
    }
    drawgrid_cour_filled(curgrid);
    drawgrid_dimw_filled(curgrid);
    drawMaxAngle(curgrid);
    draw_circumcenters(curgrid);

    if (draw_prop) {
	drawgrid_prop_filled(curgrid);
    }
    if (display_delh) {
	drawdelh(curgrid);
    }
    if (display_grad) {
	drawgrad(curgrid);
    }
    if (display_grad_legend) {
	drawgradlegend(curgrid);
    }
    if (display_slope) {
	drawslopes(curgrid);
    }
    if (display_nnels) {
        extern int nnelsop, nnels;
        drawSelectedNodes(curgrid, nnelsop, nnels, 4);
    }
    if (pointset) {
	my_circle(dsx, dsy);
    }
    setcolor(1);
    if (display_flags[EDIT_GRID_NODE_NUMBERS] || display_flags[EDIT_GRID_NODES]) {
	drawgridnodes(curgrid,
		      display_flags[EDIT_GRID_NODE_NUMBERS],
		      display_flags[EDIT_GRID_NODES]);
    }
    if (display_flags[EDIT_GRID_DEPTHS]) {
	drawnodedepths(curgrid);
    }
    if (display_flags[EDIT_GRID_ELEMENT_NUMBERS]) {
	drawgridelems(curgrid);
    }
    setcolor(backgrid_color);
    if (display_flags[BACKGROUND_GRID]) {
	drawgrid(backgrid, 1.0);
    }
    setcolor(1);
    if (display_flags[BUILD_POINTS]) {
	drawbuild(curbuild);
    }
    if (display_flags[BACKGROUND_BOUNDARY]) {
	setcolor(backbound_color);
	draw_boundary(backgrid);
	setcolor(1);
    }
    if (region_flag) {
	setcolor(1);
	draw_region(regionx, regiony, nregion);
    }
    if (display_flags[EDIT_GRID]) {
	setcolor(grid_color);
	drawgrid(curgrid, redfact);
	setcolor(1);
    }
    if (display_flags[EDIT_BOUNDARY]) {
	setcolor(bound_color);
	draw_boundary(curgrid);
	setcolor(1);
    }
    if (displaynl) {
	DrawPNodeLists(curgrid);
    }
    draw_adcirc();
    setcolor(1);
    if (display_mapscale) {
	drawscale();
    }
    if (nsta) {
	DrawStations(nsta, sta);
    }
    if (display_shapiro) {
	find_sh(curgrid);
    }
/*
   drawisollegend(ip.isollegx, ip.isollegy, ip.cis, ip.nisol, 0, 0,
   ip.isolformat, ip.isolprec, legend_color,
   cmap, nmap, setvideo);
 */
}

/*
 * set hardcopy flag and if writing to a file, check
 * to see if it exists
 */
void do_hardcopy(void)
{
    FILE *fp;
    int i, j, k;
    char buf[256];
    extern double mapscale;
    double x1 = xg1, x2 = xg2, y1 = yg1, y2 = yg2;
    write_mode_str("Printing, press Cancel to abort");
    if (ptofile) {
	if (fexists(printstr)) {
	    hardcopyflag = FALSE;
	    return;
	}
	fp = fopen(printstr, "w");
	if (fp == NULL) {
	    sprintf(buf, "Can't open %s for write, hardcopy aborted", printstr);
	    errwin(buf);
	    hardcopyflag = FALSE;
	    return;
	}
	fclose(fp);
    }
    i = 0;
    hardcopyflag = TRUE;
    setup_abort();
    if (initgraphics(hdevice) != 1) {
	if (mapscale != 1.0) {
	    set_up_mapscale(mapscale);
	} else {
	    set_up_world();
	}
	defineworld(xg1, yg1, xg2, yg2, 0, 0);
	viewport(xv1, yv1, xv2, yv2);
	set_wait_cursor(NULL);
	do_drawobjects();
	draw_annotation(0);
	draw_annotation(-1);
	leavegraphics();
	unset_wait_cursor(NULL);
	write_mode_str(NULL);
    } else {
	errwin("Hardcopy failed");
    }
    unsetup_abort();
    hardcopyflag = FALSE;
    doclear = 0;
    xg1 = x1;
    xg2 = x2;
    yg1 = y1;
    yg2 = y2;
    initgraphics(tdevice);
    set_up_world();
    doclear = 1;
    write_mode_str(NULL);
}

/*
 * define world coordinate system and initialize
 * the graphics drivers
 */
void set_up_world(void)
{
    set_scale();
    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(xv1, yv1, xv2, yv2);
}

void set_scale(void)
{
    double scalex, scaley;
    if (fixed_scale) {
	setfixedscale(xv1, yv1, xv2, yv2, &xg1, &yg1, &xg2, &yg2);
    }
}

void set_up_mapscale(double mapscale)
{
    setmapscale(mapscale, xv1, yv1, xv2, yv2, &xg1, &yg1, &xg2, &yg2);
}

void reset_world(void)
{
    xg1 = sxg1;
    yg1 = syg1;
    xg2 = sxg2;
    yg2 = syg2;
    set_scale();
    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(xv1, yv1, xv2, yv2);
    wfact = 0.30;
    do_drawgrid();
}

void draw_region(double *rx, double *ry, int nr)
{
    int i;

    if (nr >= 3) {
	setlinewidth(3);
	my_move2(rx[0], ry[0]);
	for (i = 0; i < nr; i++) {
	    my_draw2(rx[i], ry[i]);
	}
	my_draw2(rx[0], ry[0]);
	setlinewidth(1);
    }
}

void dolegend(Isolparms ip, int *cmap, int nmap)
{
    double x, y;
    int f;
    switch (ip.p.format) {
    case DECIMAL:
	f = 0;
	break;
    case EXPONENTIAL:
	f = 1;
	break;
    case GENERAL:
	f = 2;
	break;
    }
    if (ip.loctype == WORLD) {
	x = ip.x;
	y = ip.y;
    } else {
	view2world(ip.x, ip.y, &x, &y);
    }
    if (ip.layout == HORIZONTAL) {
	drawisollegendh(x, y, ip.cis, ip.nisol, 0, 0,
			f, ip.p.prec, ip.p.color,
			cmap, nmap, 0);
    } else {
	drawisollegend(x, y, ip.cis, ip.nisol, 0, 0,
		       f, ip.p.prec, ip.p.color,
		       cmap, nmap, 0);
    }
}
