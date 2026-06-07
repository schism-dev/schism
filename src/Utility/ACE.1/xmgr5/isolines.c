/*
 * $Id: isolines.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 * $Source: /home/workspace/ccalmr/src/ace/xmgr5/isolines.c,v $
 * isoline utilities
 */
 
#ifndef lint
static char RCSid[] = "$Id: isolines.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $";
#endif
 
#include <stdio.h>
#include <stdlib.h>
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
        ip->nisol = 20;
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

void autoscale_xybox(int gno)
{ 
    int i, j, first = 1;
    double cmin, cmax, *z;
    plotarr p;
    for (i = 0; i < g[gno].maxplot; i++) {
         if (isactive(gno, i)) {
            if (dataset_type(gno, i) == XYBOX) {
                p = g[gno].p[i];
		z = p.ex[4];
		for (j = 0; j < p.len; j++) {
		    if (first) {
		 	first = 0;
			cmin = cmax = z[0];
		    } else {
			if (z[j] < cmin) cmin = z[j];
			if (z[j] > cmax) cmax = z[j];
		    }
		}
            }
        }
    }
    if (!first)  {
	g[gno].isol.cmin = cmin; 
	g[gno].isol.cmax = cmax; 
	default_isolines(gno, &g[gno].isol);
    }
}

void get_xybox_minmax(int gno)
{ 
    int i, j, first = 1;
    double cmin, cmax, *z;
    plotarr p;
    for (i = 0; i < g[gno].maxplot; i++) {
         if (isactive(gno, i)) {
            if (dataset_type(gno, i) == XYBOX) {
                p = g[gno].p[i];
		z = p.ex[4];
		for (j = 0; j < p.len; j++) {
		    if (first) {
		 	first = 0;
			cmin = cmax = z[0];
		    } else {
			if (z[j] < cmin) cmin = z[j];
			if (z[j] > cmax) cmax = z[j];
		    }
		}
            }
        }
    }
    if (!first)  {
	g[gno].isol.cmin = cmin; 
	g[gno].isol.cmax = cmax; 
    }
}


void autoscale_isolines(int gno, double cmin, double cmax, Isolparms * ip)
{
    ip->cmin = cmin;
    ip->cmax = cmax;
    default_isolines(gno, ip);
}
