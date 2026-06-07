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
 * draw drogues
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawdrogs.c,v 1.3 2003/12/11 16:36:20 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

int show_drognumbs = 0;		/* TODO */

void draw_init_drogues(void)
{
    extern double drogx[], drogy[];
    extern int ndrogs;
    int i;
    char buf[10];
    if (ndrogs > 0) {
	for (i = 0; i < ndrogs; i++) {
	    my_circle(drogx[i], drogy[i]);
	    if (show_drognumbs) {
		sprintf(buf, " %d", i + 1);
		writestr(drogx[i], drogy[i], 0, 0, buf);
	    }
	}
    }
}

void drawdrogues(int gno, int drno, int itime)
{
    int j, k, ix, iy, ixp, iyp, dnum;
    double size;
    int sym;
    char buf[80];
    int xconv(), yconv();
    if (g[gno].drogues[drno].display_connect == OFF) {
	for (k = 0; k < drogues[drno].p[itime].npts; k++) {
	    dnum = drogues[drno].p[itime].drnum[k] - 1;
	    setcolor(g[gno].drogues[drno].p.color);
	    sym = g[gno].drogues[drno].p.symbol;
	    size = g[gno].drogues[drno].p.symsize;
/*
	    setcolor(g[gno].drogues[drno].color[dnum]);
	    sym = g[gno].drogues[drno].sym[dnum];
	    size = g[gno].drogues[drno].size[dnum];
    if (g[gno].drogues[drno].display_connect == TIME) {
	setcolor(g[gno].drogues[drno].color[dnum]);
*/
	    if (g[gno].drogues[drno].display_streaml == OFF) {
		drawpolysym(&drogues[drno].p[itime].x[k], &drogues[drno].p[itime].y[k], 1, sym, 0, 1, size * 0.8);
		if (g[gno].drogues[drno].display_type == NUMBER) {
		    sprintf(buf, " %d", drogues[drno].p[itime].drnum[k]);
		    writestr(drogues[drno].p[itime].x[k], drogues[drno].p[itime].y[k], 0, 0, buf);
		} else if (g[gno].drogues[drno].display_type == DEPTH) {
		    sprintf(buf, " %.2lf", -drogues[drno].p[itime].z[k]);
		    writestr(drogues[drno].p[itime].x[k], drogues[drno].p[itime].y[k], 0, 0, buf);
		}
	    } else {
		drawpolysym(&drogues[drno].p[0].x[k], &drogues[drno].p[0].y[k], 1, 2, 0, 1, 0.5);
		if (g[gno].drogues[drno].display_type == NUMBER) {
		    sprintf(buf, " %d", drogues[drno].p[0].drnum[k]);
		    writestr(drogues[drno].p[0].x[k], drogues[drno].p[0].y[k], 0, 0, buf);
		} else if (g[gno].drogues[drno].display_type == DEPTH) {
		    sprintf(buf, " %.2lf", -drogues[drno].p[itime].z[k]);
		    writestr(drogues[drno].p[itime].x[k], drogues[drno].p[itime].y[k], 0, 0, buf);
		}
		if (itime > 0) {
		    my_move2(drogues[drno].p[0].x[k], drogues[drno].p[0].y[k]);
		    for (j = 0; j < itime; j++) {
			my_draw2(drogues[drno].p[j].x[k], drogues[drno].p[j].y[k]);
		    }
		}
	    }
	}
    } else {
	setcolor(g[gno].drogues[drno].p.color);
	setlinestyle(g[gno].drogues[drno].p.lines);
	setlinewidth(g[gno].drogues[drno].p.linew);
	my_move2(drogues[drno].p[itime].x[0], drogues[drno].p[itime].y[0]);
	for (j = 0; j < drogues[drno].p[itime].npts; j++) {
	    my_draw2(drogues[drno].p[itime].x[j], drogues[drno].p[itime].y[j]);
	}
	my_draw2(drogues[drno].p[itime].x[0], drogues[drno].p[itime].y[0]);
    }
}

void set_drogscolor(int gno, int drno, int n, int c)
{
    /*g[gno].drogues[drno].color[n - 1] = c; */
}
