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
 * Draw elevations markers
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawelev.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

double eval_teanl_node(int flowno, int node, double t);

void draw_elevationmarkers(int gno, int n, int flowno, double time, int which)
{
    int i, j;
    double e;
    Elevmarker el;
    for (i = 0; i < MAXELEVMARKERS; i++) {
	switch (which) {
	case 0:		/* TEANL */
	    el = g[gno].flowf[flowno].em[i];
	    if (el.display == ON && el.active == ON) {
		e = eval_teanl_node(flowno, el.node, time);
		putelevmarker(el.x, el.y, el.locx, el.locy, el.emin, el.emax, e);
	    }
	    break;
	case 1:		/* ADCIRC */
	    el = g[gno].flowt[flowno].em[i];
	    if (el.display == ON && el.active == ON) {
		e = flowt[flowno].f[n].e[el.node];
		putelevmarker(el.x, el.y, el.locx, el.locy, el.emin, el.emax, e);
	    }
	    break;
	case 3:		/* ADCIRC 3D */
	    break;
	}
    }
}

void register_elevmarker(int gno, int flowno, int elno, int type, double wx1, double wy1, double wx2, double wy2)
{
    int ind, k, gridno;
    Elevmarker *el;
    switch (type) {
    case TEANL:		/* TEA-NL */
	find_nearest_node(flowf[flowno].grid, wx1, wy1, &ind);
	gridno = flowf[flowno].grid;
	if (ind >= 0) {
	    el = &g[gno].flowf[flowno].em[elno];
	    el->x = grid[gridno].xord[ind];
	    el->y = grid[gridno].yord[ind];
	}
	break;
    case ADCIRC:		/* ADCIRC */
	FindNearestNode(&flowt[flowno].g, wx1, wy1, &ind);
	if (ind >= 0) {
	    el = &g[gno].flowt[flowno].em[elno];
	    el->x = flowt[flowno].g.xord[ind];
	    el->y = flowt[flowno].g.yord[ind];
	}
	break;
    case ADCIRC3DFLOW:		/* ADCIRC 3D */
	break;
    }
    if (ind >= 0) {
	el->active = ON;
	el->node = ind;
	el->loctype = WORLD;
	el->locx = wx2;
	el->locy = wy2;
	my_frontbuffer(1);
	putelevmarker(el->x, el->y, el->locx, el->locy, el->emin, el->emax, 0.0);
	my_frontbuffer(0);
    }
}

int getnearest_elevmarker(int gno, int flowno, int type, double wx1, double wy1)
{
    int ind = -1, i;
    Elevmarker *el;
    double dist, mind = 1e38, lx, ly;
    for (i = 0; i < MAXELEVMARKERS; i++) {
	switch (type) {
	case TEANL:
	    el = &g[gno].flowf[flowno].em[i];
	    break;
	case ADCIRC:
	    el = &g[gno].flowt[flowno].em[i];
	    break;
	case ADCIRC3DFLOW:
	    break;
	}
	if (el->active == ON) {
	    lx = el->locx;
	    ly = el->locy;
	    dist = hypot(wx1 - lx, wy1 - ly);
	    if (mind > dist) {
		ind = i;
		mind = dist;
	    }
	}
    }
    return ind;
}
