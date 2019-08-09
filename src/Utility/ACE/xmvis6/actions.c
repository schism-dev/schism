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
 * action procs
 */

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "symdefs.h"

#ifndef lint
static char RCSid[] = "$Id: actions.c,v 1.3 2003/08/25 15:32:14 pturner Exp $";
#endif

/* TODO */
extern int hardcopyflag;
extern int cursortype;
extern int cursor_oldx;
extern int cursor_oldy;

void do_animat();
void freshen_display(int restart);
void display_image(void);
void drawgraph(void);
void boxplot(int gno, framep f);
void draw_string(int gno, int i);
void draw_annotation(int gno);
void domapscale(int gno);
void dovscale(int gno);
void dolegend(Isolparms ip, int *cmap, int nmap);
void drawtimeline(int gno, double time);
double get_current_time(void);

extern int slice_adcircelev;

extern int tdevice;

static int restart = 0;

double compteanlelevation();
extern int mapisolbath[];
extern int mapisolconc[];
extern int nslice;

#include <Xm/Xm.h>

void my_blowup(double x1, double y1, double x2, double y2)
{
    int restart = 0;
    void setirun(void);
    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    g[cg].w.xg1 = x1;
    g[cg].w.xg2 = x2;
    g[cg].w.yg1 = y1;
    g[cg].w.yg2 = y2;
    set_up_world(cg);
    doclear = 1;
    dobackground = 1;
    display_image();
    if (restart) {
	setirun();
    }
}

void go_blowup(int ox, int oy, int old_x, int old_y)
{
    double dx, dy, x1, y1, x2, y2;
    device2world(ox, oy, &x1, &y1);
    device2world(old_x, old_y, &x2, &y2);
    if (x1 > x2)
	fswap(&x1, &x2);
    if (y1 > y2)
	fswap(&y1, &y2);
    g[cg].w.xg1 = x1;
    g[cg].w.xg2 = x2;
    g[cg].w.yg1 = y1;
    g[cg].w.yg2 = y2;
    set_up_world(cg);
    freshen_display(restart);
}

/*
 * define a region
 */
void do_select_region(void)
{
    region_pts = 0;
    set_action(0);
    set_action(DEF_REGION);
}

void do_clear_region2(void)
{
    setcolor(0);
    draw_single_region(regionx, regiony, nregion);
    setcolor(1);
    nregion = 0;
    region_flag = 0;
}

void do_select_single_region(void)
{
    if (region_flag) {
	if (yesno("Region defined, clear?", "Press Yes or No", "Yes", "No")) {
	    do_clear_region2();
	} else {
	    errwin("Region not cleared");
	    return;
	}
    }
    set_action(0);
    set_action(DEFINE_REGION);
}

/*
 * routines to determine if a point lies in a polygon
*/
int intersect_to_left2(double x, double y, double x1, double y1, double x2, double y2)
{
    double xtmp, m, b;

    /* ignore horizontal lines */
    if (y1 == y2) {
	return 0;
    }
    /* not contained vertically */
    if (((y < y1) && (y < y2)) || ((y > y1) && (y > y2))) {
	return 0;
    }
    /* none of the above, compute the intersection */
    if ((xtmp = x2 - x1) != 0.0) {
	m = (y2 - y1) / xtmp;
	b = y1 - m * x1;
	xtmp = (y - b) / m;
    } else {
	xtmp = x1;
    }
    if (xtmp <= x) {
	/* check for intersections at a vertex */
	/* if this is the max ordinate then accept */
	if (y == y1) {
	    if (y1 > y2) {
		return 1;
	    } else {
		return 0;
	    }
	}
	/* check for intersections at a vertex */
	if (y == y2) {
	    if (y2 > y1) {
		return 1;
	    } else {
		return 0;
	    }
	}
	/* no vertices intersected */
	return 1;
    }
    return 0;
}

/*
 * determine if (x,y) is in the polygon xlist[], ylist[]
 */
int inbound2(double x, double y, double *xlist, double *ylist, int n)
{
    int i, l = 0, ll = 0;

    for (i = 0; i < n; i++) {
	if ((y < ylist[i] && y < ylist[(i + 1) % n]) || (y > ylist[i] && y > ylist[(i + 1) % n])) {
	    continue;
	}
	if ((x < xlist[i] && x < xlist[(i + 1) % n])) {
	    continue;
	}
	l += intersect_to_left2(x, y, xlist[i], ylist[i], xlist[(i + 1) % n], ylist[(i + 1) % n]);
    }
    return l % 2;
}

int inregion2(double *boundx, double *boundy, int nbpts, double x, double y)
{
    return (inbound2(x, y, boundx, boundy, nbpts));
}

void setzoom(void)
{
    set_restart();
    set_action(0);
    set_action(ZOOM_1ST);
}

void set_histloc(void)
{
    set_restart();
    set_action(0);
    set_action(PLACE_HIST1ST);
}

void set_sliceloc(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_BOX1ST);
}

void set_fluxloc(void)
{
    set_restart();
    set_action(0);
    set_action(FLUX_BOX);
}

void set_zoomloc(void)
{
    set_restart();
    set_action(0);
    set_action(ZOOM_BOX1ST);
}

void set_teanlelev(void)
{
    set_restart();
    set_action(0);
    set_action(PLACE_TEANL_ELEV1ST);
}

void set_adcircelev(void)
{
    set_restart();
    set_action(0);
    set_action(PLACE_ADCIRC_ELEV1ST);
}

void set_teanlsample(void)
{
    set_restart();
    set_action(0);
    set_action(SAMPLE_TEANL);
}

void del_teanlsample(void)
{
    set_restart();
    set_action(0);
    set_action(DEL_SAMPLE_TEANL);
}

void do_slicebox(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_BOX1ST);
}

void do_sliceline(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_LINE1ST);
}

void do_slicepoly(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_POLY);
}

void do_fluxbox(void)
{
    set_restart();
    set_action(0);
    set_action(FLUX_BOX);
}

void do_fluxline(void)
{
    set_restart();
    set_action(0);
    set_action(FLUX_LINE1ST);
}

void do_fluxpoly(void)
{
    set_restart();
    set_action(0);
    set_action(FLUX_POLY);
}

void slice_teanl_elev(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_TEANL_ELEV1ST);
}

void slice_teanl_flow(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_TEANL_FLOW1ST);
}

void slice_ela(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_ELA1ST);
}

void slice_adcirc_elev(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_ADCIRC_ELEV1ST);
}

void slice_adcirc_flow(void)
{
    set_restart();
    set_action(0);
    set_action(SLICE_ADCIRC_FLOW1ST);
}

void set_adc3dloc(void)
{
    set_restart();
    set_action(0);
    set_action(PLACE_ADCIRC3D);
}

void pick_transxy(void)
{
    set_restart();
    set_action(0);
    set_action(PICK_ADCIRC3D_TRANSECT);
}

void pick_adc3dnode(void)
{
    set_restart();
    set_action(0);
    set_action(PICK_ADCIRC3D_NODE);
}

void pick_adc3dxy(void)
{
    set_restart();
    set_action(0);
    set_action(PICK_ADCIRC3D_XY);
}

void query_adc3dxy(void)
{
    set_restart();
    set_action(0);
    set_action(QUERY_ADCIRC3D_XY);
}

void query_adc3dnode(void)
{
    set_restart();
    set_action(0);
    set_action(QUERY_ADCIRC3D_NODE);
}

void set_loc(int gno, int type, int which, double x, double y)
{
    double xconv(), yconv();
    extern Isolparms curip;
    switch (type) {
    case ISOLINES:
	if (curip.loctype == VIEW) {
	    curip.x = xconv(x);
	    curip.y = yconv(y);
	} else {
	    curip.x = x;
	    curip.y = y;
	}
	break;
    case VSCALE:
	if (g[gno].vl.loctype == VIEW) {
	    g[gno].vl.x = xconv(x);
	    g[gno].vl.y = yconv(y);
	} else {
	    g[gno].vl.x = x;
	    g[gno].vl.y = y;
	}
	break;
    case WSCALE:
	if (g[gno].wl.loctype == VIEW) {
	    g[gno].wl.x = xconv(x);
	    g[gno].wl.y = yconv(y);
	} else {
	    g[gno].wl.x = x;
	    g[gno].wl.y = y;
	}
	break;
    case FLUX:
	if (g[gno].fl.loctype == VIEW) {
	    g[gno].fl.x = xconv(x);
	    g[gno].fl.y = yconv(y);
	} else {
	    g[gno].fl.x = x;
	    g[gno].fl.y = y;
	}
	break;
    case MAPSCALE:
	if (g[gno].mapscale.loctype == VIEW) {
	    g[gno].mapscale.x = xconv(x);
	    g[gno].mapscale.y = yconv(y);
	} else {
	    g[gno].mapscale.x = x;
	    g[gno].mapscale.y = y;
	}
	break;
    case TIDALCLOCK:
	if (g[gno].tidalclock.loctype == VIEW) {
	    g[gno].tidalclock.x = xconv(x);
	    g[gno].tidalclock.y = yconv(y);
	} else {
	    g[gno].tidalclock.x = x;
	    g[gno].tidalclock.y = y;
	}
	break;
    case TIMELINE:
	if (g[gno].timeline.loctype == VIEW) {
	    g[gno].timeline.x = xconv(x);
	    g[gno].timeline.y = yconv(y);
	} else {
	    g[gno].timeline.x = x;
	    g[gno].timeline.y = y;
	}
	break;
    case TIMEINFO:
	if (g[gno].timeinfo.loctype == VIEW) {
	    g[gno].timeinfo.x = xconv(x);
	    g[gno].timeinfo.y = yconv(y);
	} else {
	    g[gno].timeinfo.x = x;
	    g[gno].timeinfo.y = y;
	}
	break;
    case NORTH:
	if (g[gno].north.loctype == VIEW) {
	    g[gno].north.x = xconv(x);
	    g[gno].north.y = yconv(y);
	} else {
	    g[gno].north.x = x;
	    g[gno].north.y = y;
	}
	break;
    }
}
