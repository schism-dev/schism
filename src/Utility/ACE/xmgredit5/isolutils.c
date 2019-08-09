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
 * graph utilities
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isolutils.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

static Props d_props = {1,	/* color */
    1,				/* line width */
    1,				/* line style */
    4,				/* font */
    DECIMAL,			/* format */
    1,				/* precision */
    14,				/* point size */
    1.0,			/* char size */
    0,				/* symbol */
    1.0,			/* symbol size */
    OFF,			/* fill */
    COLOR,			/* fill using */
    1,				/* fill color */
    1,				/* fill pattern */
};

double expt(), nicenum();

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

void default_isolines(Isolparms * ip)
{
    double d = (ip->cmax - ip->cmin);
    int i;
    if (ip->nisol == 0) {
	ip->nisol = 16;
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
    grid[gridno].ip.cmin = grid[gridno].dmin;
    grid[gridno].ip.cmax = grid[gridno].dmax;
    default_isolines(&grid[gridno].ip);
}

void autoscale_isolines(double cmin, double cmax, Isolparms * ip)
{
    ip->cmin = cmin;
    ip->cmax = cmax;
    default_isolines(ip);
}

void set_isolines_defaults(Isolparms * ip)
{
    int i;
    ip->nisol = 16;
    ip->type = 1;
    ip->isoltype = 0;
    ip->cmin = 0.0;
    ip->cmax = 1.0;
    ip->cint = 0.0;
    ip->p = d_props;
    for (i = 0; i < MAXISOLINES; i++) {
	ip->cis[i] = 0.0;
	ip->linew[i] = 1;
	ip->lines[i] = 1;
	ip->color[i] = i + 1;
    }
    ip->marker = 0;
    ip->markstep = 0;
    ip->p.prec = 2;
    ip->p.format = DECIMAL;
    ip->nsplits = 0;
    ip->layout = VERTICAL;
    ip->llabels = 1;
    ip->xlen = 1;
    ip->ylen = 2;
    ip->loctype = WORLD;
    ip->x = 0.0;
    ip->y = 0.0;
    ip->xgap = 0;
    ip->ygap = 2;
}

void set_program_defaults(void)
{
    int i;
    for (i = 0; i < MAXGRIDS + 3; i++) {
	set_isolines_defaults(&grid[i].ip);
    }
    for (i = 0; i < MAXBUILD; i++) {
	DefaultsBuild(&build[i]);
    }
}
