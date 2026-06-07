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
 * draw the clipping boundary
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawbound.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

static int bcolor = 1;
static int fcolor = 0;

void draw_boundary(int gno, int gridno)
{
    int i, j, bno;
    char s[128];
    Boundary boundary;

    if (use_colors <= 2) {
	bcolor = 1;
	fcolor = 1;
    }
    setcolor(g[gno].grid[gridno].bp.color);
    setlinestyle(g[gno].grid[gridno].bp.lines);
    setlinewidth(g[gno].grid[gridno].bp.linew);
    for (i = 0; i < grid[gridno].nbounds; i++) {
	boundary = grid[gridno].boundaries[i];
	if (boundary.active == ON) {
	    switch (g[gno].grid[gridno].display_boundary) {
	    case ON:
		if (use_colors > 2) {
		    if (i != 0) {
			if (g[gno].grid[gridno].bp.fill == ON) {
			    setcolor(g[gno].grid[gridno].bp.fillcol);
			    fillcolor(boundary.nbpts, boundary.boundx, boundary.boundy);
			    setcolor(g[gno].grid[gridno].bp.color);
			}
		    } else {
/*
			setcolor(g[gno].grid[gridno].bp.color);
			fillcolor(boundary.nbpts, boundary.boundx, boundary.boundy);
*/
		    }
		}
		my_move2(boundary.boundx[0], boundary.boundy[0]);
		for (j = 1; j < boundary.nbpts; j++) {
		    my_draw2(boundary.boundx[j], boundary.boundy[j]);
		}
		my_draw2(boundary.boundx[0], boundary.boundy[0]);
		break;
	    case SYMBOL:
		drawpolysym(boundary.boundx, boundary.boundy, boundary.nbpts, 2, 0, 1, 0.5);
		break;
	    case NODES:
		for (j = 0; j < boundary.nbpts; j++) {
		    my_move2(boundary.boundx[j], boundary.boundy[j]);
		    sprintf(s, "%d, %d", j + 1, boundary.nodes[j] + 1);
		}
		break;
	    }
	}
    }
}

void draw_boundary2(int gno, int boundno)
{
    int i, j, bno;
    char s[128];
    Boundary boundary;

    if (use_colors <= 2) {
	bcolor = 1;
	fcolor = 1;
    }
    setcolor(g[gno].bounds[boundno].p.color);
    setlinestyle(g[gno].bounds[boundno].p.lines);
    setlinewidth(g[gno].bounds[boundno].p.linew);
    for (i = 0; i < bounds[boundno].nbounds; i++) {
	boundary = bounds[boundno].boundaries[i];
	if (boundary.active == ON) {
	    switch (g[gno].bounds[boundno].display) {
	    case ON:
		if (use_colors > 2) {
		    if (i != 0) {
			if (g[gno].bounds[boundno].p.fill == ON) {
			    setcolor(g[gno].bounds[boundno].p.fillcol);
			    fillcolor(boundary.nbpts, boundary.boundx, boundary.boundy);
			    setcolor(g[gno].bounds[boundno].p.color);
			}
		    } else {
/*
			setcolor(g[gno].bounds[boundno].p.color);
			fillcolor(boundary.nbpts, boundary.boundx, boundary.boundy);
*/
		    }
		}
		my_move2(boundary.boundx[0], boundary.boundy[0]);
		for (j = 1; j < boundary.nbpts; j++) {
		    my_draw2(boundary.boundx[j], boundary.boundy[j]);
		}
		my_draw2(boundary.boundx[0], boundary.boundy[0]);
		break;
	    case SYMBOL:
		drawpolysym(boundary.boundx, boundary.boundy, boundary.nbpts, g[gno].bounds[boundno].p.symbol, 0, 1, g[gno].bounds[boundno].p.symsize * 0.01);
		break;
	    case NODES:
		for (j = 0; j < boundary.nbpts; j++) {
		    my_move2(boundary.boundx[j], boundary.boundy[j]);
		    sprintf(s, "%d, %d", j + 1, boundary.nodes[j] + 1);
		}
		break;
	    case XYSEG:
		for (j = 0; j < boundary.nbpts; j += 2) {
		    if (j + 1 < boundary.nbpts) {
			my_move2(boundary.boundx[j], boundary.boundy[j]);
			my_draw2(boundary.boundx[j + 1], boundary.boundy[j + 1]);
		    }
		}
		break;
	    }
	}
    }
}
