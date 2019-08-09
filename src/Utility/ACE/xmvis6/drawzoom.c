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
 * Draw zoom markers
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawzoom.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

void draw_zoom(int gno, int n)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vxinc, vyinc, vx1, vx2, vy1, vy2;
    char buf[10];

    /* fill the region on the graph */
    setcolor(g[gno].zbox[n].rp.color);
    rect(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1, g[gno].zbox[n].wx2, g[gno].zbox[n].wy2);
/*
    setcolor(g[gno].zbox[n].rp.fillcol);
    fillrectcolor(g[gno].zbox[n].wx1,
		  g[gno].zbox[n].wy1,
		  g[gno].zbox[n].wx2,
		  g[gno].zbox[n].wy2);
*/
    /* connect the dots */
    setcolor(g[gno].zbox[n].p.color);
    my_move2(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1);
    my_draw2(g[gno].zbox[n].locx, g[gno].zbox[n].locy);
    vxinc = g[gno].zbox[n].vx * g[gno].zbox[n].expand;
    vyinc = g[gno].zbox[n].vy * g[gno].zbox[n].expand;
    switch (g[gno].zbox[n].attach) {
    case 0:
	world2view(g[gno].zbox[n].locx, g[gno].zbox[n].locy, &vx1, &vy1);
	vx2 = vx1 + vxinc;
	vy2 = vy1 + vyinc;
	break;
    case 1:
	world2view(g[gno].zbox[n].locx, g[gno].zbox[n].locy, &vx2, &vy1);
	vx1 = vx2 - vxinc;
	vy2 = vy1 + vyinc;
	break;
    case 2:
	world2view(g[gno].zbox[n].locx, g[gno].zbox[n].locy, &vx1, &vy2);
	vx2 = vx1 + vxinc;
	vy1 = vy2 - vyinc;
	break;
    case 3:
	world2view(g[gno].zbox[n].locx, g[gno].zbox[n].locy, &vx2, &vy2);
	vx1 = vx2 - vxinc;
	vy1 = vy2 - vyinc;
	break;
    }

/* redefine the world and viewport */
    defineworld(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1, g[gno].zbox[n].wx2, g[gno].zbox[n].wy2, 0, 0);
    viewport(vx1, vy1, vx2, vy2);

/* the zoom box region */
    setcolor(g[gno].zbox[n].p.fillcol);
    fillrectcolor(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1, g[gno].zbox[n].wx2, g[gno].zbox[n].wy2);

/* the zoom box outline and ticks */
    setcolor(g[gno].zbox[n].p.color);
    rect(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1, g[gno].zbox[n].wx2, g[gno].zbox[n].wy2);
    rectstr(g[gno].zbox[n].wx1, g[gno].zbox[n].wy1, g[gno].zbox[n].wx2, g[gno].zbox[n].wy2, g[gno].zbox[n].precx, g[gno].zbox[n].precy);
}

int read_zbox(int zboxno, char *fname)
{
    char buf[255];
    int i, npts;
    double x, y;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 0;
    }
    fgets(buf, 255, fp);
}
