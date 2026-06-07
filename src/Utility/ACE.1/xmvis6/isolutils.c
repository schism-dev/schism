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
 * isoline utilities
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isolutils.c,v 1.3 2006/12/27 01:02:59 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

double expt(double a, register int n), nicenum(double x, int round);

void default_range(int num, double *gmin, double *gmax, double *dd)
{
    double range, d = *dd, tmpmax = *gmax, tmpmin = *gmin;
    range = nicenum(tmpmax - tmpmin, 0);
    d = nicenum(range / (num - 1), 1);
    tmpmin = floor(tmpmin / d) * d;
    tmpmax = ceil(tmpmax / d) * d;
    *gmax = tmpmax;
    *gmin = tmpmin;
    *dd = d;
}

void default_isolines(int gno, Isolparms * ip)
{
    double d = (ip->cmax - ip->cmin);
    int i;
    if (ip->nisol == 0) {
	ip->nisol = MAXISOLINES;
    }
    if (ip->nisol > 1) {
/*
    default_range(gno, ip->nisol, &(ip->cmin), &(ip->cmax), &(ip->cint));
*/
	ip->cint = (ip->cmax + 0.02 * d - ip->cmin) / (ip->nisol - 1);
	for (i = 0; i < ip->nisol; i++) {
	    ip->cis[i] = ip->cmin + i * ip->cint;
	}
    }
}

void autoscale_grid(int gno, int gridno)
{
    g[gno].grid[gridno].ip.cmin = grid[gridno].dmin;
    g[gno].grid[gridno].ip.cmax = grid[gridno].dmax;
    default_isolines(gno, &g[gno].grid[gridno].ip);
}

void autoscale_isolines(double cmin, double cmax, Isolparms * ip)
{
    ip->cmin = cmin;
    ip->cmax = cmax;
    default_isolines(cg, ip);
}
