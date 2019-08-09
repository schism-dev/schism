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
 * set defaults - changes to the types in defines.h
 * will require changes in here also
 *
 */

#ifndef lint
static char RCSid[] = "$Id: defaults.c,v 1.13 2007/01/13 21:48:46 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

void set_program_defaults(void);
void set_region_defaults(int i);
void set_default_defaults();
void set_default_framep(framep * f);
void set_default_world(world * w);
void set_default_view(view * v);
void set_default_string(plotstr * s);
void set_default_line(linetype * l);
void set_default_box(boxtype * b);
void set_default_graph(int gno);
void set_default_annotation(void);
void set_default_ticks(tickmarks * t, int a);
void set_default_elevmarkers(Elevmarker * elev, int n);
void set_default_concmarkers(Conc_marker * conc, int n);
void set_default_histmarkers(int gno);
void set_default_sliceboxes(int gno);
void set_default_fluxboxes(int gno);
void set_default_zoomboxes(int gno);
void set_isolines_defaults(int gno, Isolparms * ip);
void set_default_transect(Transect * t);

static Props d_props = { 1,	/* color */
    1,				/* line width */
    1,				/* line style */
    4,				/* font */
    DECIMAL,			/* format */
    1,				/* precision */
    14,				/* point size */
    1.0,			/* char size */
    0,				/* symbol */
    100,			/* symbol size */
    OFF,			/* fill */
    COLOR,			/* fill using */
    1,				/* fill color */
    1,				/* fill pattern */
};

static linetype d_l = { OFF,	/* active */
    VIEW,			/* location type */
    -1,				/* if loctype == WORLD then graph number */
    0.0,			/* x1 */
    0.0,			/* y1 */
    0.0,			/* x2 */
    0.0,			/* y2 */
    1,				/* line style */
    1,				/* line width */
    1,				/* color */
    0,				/* arrow type */
    0,				/* arrow head type */
    1.0				/* arrow size */
};

static boxtype d_b = { OFF,	/* active */
    VIEW,			/* location type */
    -1,				/* if loctype == WORLD then graph number */
    0.0,			/* x1 */
    0.0,			/* y1 */
    0.0,			/* x2 */
    0.0,			/* y2 */
    1,				/* line style */
    1,				/* line width */
    1,				/* color */
    OFF,			/* fill */
    1,				/* fill color */
    1				/* fill pattern */
};

static plotstr d_s = { OFF,	/* active */
    VIEW,			/* location type */
    -1,				/* if loctype == WORLD then graph number */
    0.0,			/* x */
    0.0,			/* y */
    1,				/* line style */
    1,				/* line width */
    1,				/* color */
    0,				/* rotation (degrees) */
    4,				/* font */
    0,				/* justification */
    0,				/* where */
    0,				/* sym */
    RIGHT,			/* symloc */
    0,				/* symfill */
    1,				/* symcolor */
    1.0,			/* sym size */
    1.0,			/* character size */
    0				/* string */
};

static framep d_f = { OFF,	/* active */
    0,				/* type */
    1,				/* color */
    1,				/* line style */
    1,				/* line width */
    OFF,			/* fill background */
    0				/* background color */
};

static world d_w = { 0.0, 1.0, 0.0, 1.0 };

static view d_v = { 0.0, 1.0, 0.0, 1.0 };

TimeClock d_timeclock = {
    ON,				/* active */
    0,				/* type (time or steps ) */
    1,				/* direction of stepping (forward or reverse
				 * +-1 resp.) */
    1,				/* wrap flag */
    0,				/* is running */
    0.0,			/* start time of clock */
    0.0,			/* stop time of clock */
    0.0,			/* time step */
    0.0,			/* current time */
    0,				/* nsteps */
    0,				/* start step */
    0,				/* stop step */
    0,				/* current step */
    0,				/* skip */
    0,				/* interpolate inbetweens */
    NULL			/* array to hold time */
};

void set_program_defaults(void)
{
    int i;
    extern Isolparms elevip, concip;
    g = (graph *) calloc(maxgraph, sizeof(graph));
    for (i = 0; i < maxgraph; i++) {
	set_default_graph(i);
    }
    for (i = 0; i < MAXTRANS; i++) {
	set_default_transect(&trans[i]);
    }
    for (i = 0; i < MAXREGION; i++) {
	set_region_defaults(i);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	flowt[i].active = OFF;
    }
    for (i = 0; i < MAXTEANL; i++) {
	flowf[i].sign_of_phase = 1;
	flowf[i].active = OFF;
    }
    for (i = 0; i < MAXELA; i++) {
	elaconc[i].active = OFF;
    }
    for (i = 0; i < MAXBOUNDS; i++) {
	bounds[i].active = OFF;
    }
    for (i = 0; i < MAXTRACK; i++) {
	track[i].active = OFF;
    }
    for (i = 0; i < MAXADCIRC3D; i++) {
	adc3d[i].active = OFF;
	adc3d[i].sal = NULL;
	adc3d[i].u = NULL;
	adc3d[i].v = NULL;
	adc3d[i].w = NULL;
	adc3d[i].time = NULL;
    }
    set_default_annotation();
    set_default_string(&timestamp);
    sysstr = d_s;
    sysbox = d_b;
    sysline = d_l;
    timestamp.x = 0.03;
    timestamp.y = 0.03;
    scrollper = shexper = 0.1;	/* scroll and shrink/expand factors */
    timeclock = d_timeclock;
    timeclock.t = NULL;		/* paranoia */
    ndisplay = NDISPLAY;
}

void set_default_elevmarkers(Elevmarker * elev, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	elev[i].active = OFF;
	elev[i].type = -1;
	elev[i].display = ON;	/* is displayed */
	elev[i].use_node = -1;
	elev[i].node = -1;
	elev[i].p = d_props;
	elev[i].loctype = WORLD;
	elev[i].x = 0.0;
	elev[i].y = 0.0;
	elev[i].locx = 0.0;
	elev[i].locy = 0.0;
	elev[i].emin = -4.0;
	elev[i].emax = 4.0;
    }
}

void set_default_velmarkers(Velocity_marker * vel, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
	vel[i].display = OFF;
	vel[i].display_marker = OFF;
	vel[i].display_sz = DEPTH;
	vel[i].type = -1;
	vel[i].p = d_props;
	vel[i].p.fillcol = 0;	/* set fill to 0 */
	vel[i].x = 0.0;
	vel[i].y = 0.0;
	vel[i].use_node = -1;
	vel[i].node = -1;
	vel[i].locx = 0.0;
	vel[i].locy = 0.0;
	vel[i].attach = 0;
	vel[i].vx = 0.2;
	vel[i].vx = 0.1;
	vel[i].precx = -1;
	vel[i].precy = -1;
	vel[i].wx1 = 0.0;
	vel[i].wx2 = 1.0;
	vel[i].wy1 = 0.0;
	vel[i].wy2 = 1.0;
    }
}

void set_default_concmarkers(Conc_marker * conc, int n)
{
    int i, j;
    for (i = 0; i < n; i++) {
	conc[i].display = OFF;
	conc[i].display_marker = OFF;
	conc[i].type = -1;
	conc[i].p = d_props;
	conc[i].p.fillcol = 0;	/* set fill to 0 */
	conc[i].x = 0.0;
	conc[i].y = 0.0;
	conc[i].use_node = -1;
	conc[i].node = -1;
	conc[i].locx = 0.0;
	conc[i].locy = 0.0;
	conc[i].attach = 0;
	conc[i].vx = 0.2;
	conc[i].vx = 0.1;
	conc[i].precx = 1;
	conc[i].precy = 0;
	conc[i].wx1 = 0.0;
	conc[i].wx2 = 1.0;
	conc[i].wy1 = 0.0;
	conc[i].wy2 = 1.0;
    }
}

void InitDisplayTideStation(DisplayTideStation * ts)
{
    int j;
    ts->type = -1;
    ts->display = OFF;		/* history point is displayed if ON */
    ts->display_marker = OFF;	/* history marker is displayed if ON */
    ts->display_ampphase = 0;	/* */
    ts->p = d_props;
    ts->p.fillcol = 0;		/* set fill to 0 */
    ts->x = 0.0;
    ts->y = 0.0;
    ts->locx = 0.0;
    ts->locy = 0.0;
    ts->attach = 0;
    ts->vx = 0.2;
    ts->vy = 0.1;
    ts->wx1 = 0.0;
    ts->wx2 = 1.0;
    ts->wy1 = 0.0;
    ts->wy2 = 1.0;
    for (j = 0; j < MAXTEANL; j++) {
	ts->teanl[j] = 0;
	ts->tp[j] = d_props;
    }
    for (j = 0; j < MAXADCIRC; j++) {
	ts->adcirc[j] = 0;
	ts->ap[j] = d_props;
    }
    for (j = 0; j < MAXADCIRC; j++) {
	ts->adcircflow[j] = 0;
	ts->apf[j] = d_props;
    }
}

void set_default_histmarkers(int gno)
{
    int i, j;
    for (i = 0; i < MAXHISTMARKERS; i++) {
	g[gno].hbox[i].type = -1;
	g[gno].hbox[i].elem = -1;
	g[gno].hbox[i].display = OFF;	/* history point is displayed */
	g[gno].hbox[i].display_marker = OFF;	/* history marker is
						 * displayed */
	g[gno].hbox[i].p = d_props;
	g[gno].hbox[i].p.fillcol = 0;	/* set fill to 0 */
	g[gno].hbox[i].x = 0.0;
	g[gno].hbox[i].y = 0.0;
	g[gno].hbox[i].locx = 0.0;
	g[gno].hbox[i].locy = 0.0;
	g[gno].hbox[i].attach = 0;
	g[gno].hbox[i].vx = 0.2;
	g[gno].hbox[i].vy = 0.1;
	g[gno].hbox[i].wx1 = 0.0;
	g[gno].hbox[i].wx2 = 1.0;
	g[gno].hbox[i].wy1 = 0.0;
	g[gno].hbox[i].wy2 = 1.0;
	g[gno].hbox[i].thist = 0;
	g[gno].hbox[i].hp = d_props;
	for (j = 0; j < MAXTEANL; j++) {
	    g[gno].hbox[i].teanl[j] = 0;
	    g[gno].hbox[i].tp[j] = d_props;
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    g[gno].hbox[i].adcirc[j] = 0;
	    g[gno].hbox[i].ap[j] = d_props;
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    g[gno].hbox[i].adcircflow[j] = 0;
	    g[gno].hbox[i].apf[j] = d_props;
	}
	for (j = 0; j < MAXELA; j++) {
	    g[gno].hbox[i].ela[j] = 0;
	    g[gno].hbox[i].ep[j] = d_props;
	}
    }
}

void set_default_sliceboxes(int gno)
{
    int i, j;
    for (i = 0; i < MAXSLICEBOXES; i++) {
	g[gno].sbox[i].active = OFF;
	g[gno].sbox[i].type = -1;
	g[gno].sbox[i].display = OFF;	/* slice is displayed */
	g[gno].sbox[i].display_marker = OFF;	/* slice marker is displayed */
	g[gno].sbox[i].p = d_props;
	g[gno].sbox[i].p.fillcol = 0;	/* set fill to 0 */
	g[gno].sbox[i].x = 0.0;
	g[gno].sbox[i].y = 0.0;
	g[gno].sbox[i].npts = 0;
	g[gno].sbox[i].loctype = WORLD;
	g[gno].sbox[i].locx = 0.0;
	g[gno].sbox[i].locy = 0.0;
	g[gno].sbox[i].attach = 0;
	g[gno].sbox[i].sx = NULL;
	g[gno].sbox[i].sy = NULL;
	g[gno].sbox[i].vx = 0.2;
	g[gno].sbox[i].vy = 0.1;
	g[gno].sbox[i].wx1 = 0.0;
	g[gno].sbox[i].wx2 = 1.0;
	g[gno].sbox[i].wy1 = 0.0;
	g[gno].sbox[i].wy2 = 1.0;
	for (j = 0; j < 250; j++) {
	    g[gno].sbox[i].elist[j] = -1;
	}
	for (j = 0; j < MAXGRIDS; j++) {
	    g[gno].sbox[i].bath[j] = 0;
	    g[gno].sbox[i].bp[j] = d_props;
	}
	for (j = 0; j < MAXTEANL; j++) {
	    g[gno].sbox[i].teanl[j] = 0;
	    g[gno].sbox[i].tp[j] = d_props;
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    g[gno].sbox[i].adcirc[j] = 0;
	    g[gno].sbox[i].ap[j] = d_props;
	}
	for (j = 0; j < MAXELA; j++) {
	    g[gno].sbox[i].ela[j] = 0;
	    g[gno].sbox[i].ep[j] = d_props;
	}
    }
}

void set_default_fluxboxes(int gno)
{
    int i, j;
    for (i = 0; i < MAXFLUXBOXES; i++) {
	g[gno].fbox[i].active = OFF;
	g[gno].fbox[i].type = -1;
	g[gno].fbox[i].display = OFF;	/* slice is displayed */
	g[gno].fbox[i].display_marker = OFF;	/* slice marker is displayed */
	g[gno].fbox[i].p = d_props;
	g[gno].fbox[i].p.fillcol = 0;	/* set fill to 0 */
	g[gno].fbox[i].x = 0.0;
	g[gno].fbox[i].y = 0.0;
	g[gno].fbox[i].loctype = WORLD;
	g[gno].fbox[i].locx = 0.0;
	g[gno].fbox[i].locy = 0.0;
	g[gno].fbox[i].attach = 0;
	g[gno].fbox[i].vx = 0.2;
	g[gno].fbox[i].vy = 0.1;
	g[gno].fbox[i].sx = NULL;
	g[gno].fbox[i].sy = NULL;
	g[gno].fbox[i].wx1 = 0.0;
	g[gno].fbox[i].wx2 = 1.0;
	g[gno].fbox[i].wy1 = 0.0;
	g[gno].fbox[i].wy2 = 1.0;
	for (j = 0; j < MAXTEANL; j++) {
	    g[gno].fbox[i].teanl[j] = 0;
	    g[gno].fbox[i].tp[j] = d_props;
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    g[gno].fbox[i].adcirc[j] = 0;
	    g[gno].fbox[i].ap[j] = d_props;
	}
	for (j = 0; j < MAXELA; j++) {
	    g[gno].fbox[i].ela[j] = 0;
	    g[gno].fbox[i].ep[j] = d_props;
	}
    }
}


void set_default_zoomboxes(int gno)
{
    int i, j;
    for (i = 0; i < MAXZOOMBOXES; i++) {
	g[gno].zbox[i].active = OFF;
	g[gno].zbox[i].type = -1;
	g[gno].zbox[i].display = OFF;	/* zoom region is displayed */
	g[gno].zbox[i].rp = d_props;
	g[gno].zbox[i].display_marker = OFF;	/* zoom marker is displayed */
	g[gno].zbox[i].p = d_props;
	g[gno].zbox[i].p.fillcol = 0;	/* set fill to 0 */
	g[gno].zbox[i].expand = 2.0;	/* expansion factor */
	g[gno].zbox[i].x = 0.0;
	g[gno].zbox[i].y = 0.0;
	g[gno].zbox[i].loctype = WORLD;
	g[gno].zbox[i].locx = 0.0;
	g[gno].zbox[i].locy = 0.0;
	g[gno].zbox[i].attach = 0;
	g[gno].zbox[i].vx = 0.2;
	g[gno].zbox[i].vx = 0.1;
	g[gno].zbox[i].wx1 = 0.0;
	g[gno].zbox[i].wx2 = 1.0;
	g[gno].zbox[i].wy1 = 0.0;
	g[gno].zbox[i].wy2 = 1.0;
    }
}

void set_default_flow3d(int gno)
{
    int i, j;
    for (i = 0; i < MAXADCIRC3D; i++) {
	g[gno].flow3d[i].type = -1;
	g[gno].flow3d[i].display = OFF;	/* history point is displayed */
	g[gno].flow3d[i].display_marker = OFF;	/* history marker is
						 * displayed */
	g[gno].flow3d[i].display_lines = OFF;	/* lines separating levels
						 * displayed */
	g[gno].flow3d[i].display_labels = ON;	/* Display tick labels */
	g[gno].flow3d[i].nlabels = 2;	/* Number of labels */
	g[gno].flow3d[i].p = d_props;
	g[gno].flow3d[i].p.fillcol = 0;	/* set fill to 0 */
	g[gno].flow3d[i].x = 0.0;
	g[gno].flow3d[i].y = 0.0;
	g[gno].flow3d[i].locx = 0.0;
	g[gno].flow3d[i].locy = 0.0;
	g[gno].flow3d[i].attach = 0;
	g[gno].flow3d[i].vx = 0.2;
	g[gno].flow3d[i].vy = 0.3;
	g[gno].flow3d[i].wx1 = 0.0;
	g[gno].flow3d[i].wx2 = 1.0;
	g[gno].flow3d[i].wy1 = 0.0;
	g[gno].flow3d[i].wy2 = 1.0;
    }
}

void set_isolines_defaults(int gno, Isolparms * ip)
{
    int i;
    ip->nisol = 16;
    ip->type = 1;
    ip->isoltype = 0;
    ip->visflag = 0;
    ip->cmin = 0.0;
    ip->cmax = 1.0;
    ip->cint = 0.0;
    ip->p = d_props;
    for (i = 0; i < MAXISOLINES; i++) {
	ip->cis[i] = 0.0;
	ip->linew[i] = 1;
	ip->lines[i] = 1;
	ip->color[i] = 1;
    }
    ip->frame = OFF;
    ip->framecol = 0;
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

void set_vis_defaults(int gno, int type)
{
    tidal_clock *tc = &g[gno].tidalclock;
    time_line *tl = &g[gno].timeline;
    time_info *ts = &g[gno].timeinfo;
    velocity_scale *vl = &g[gno].vl;
    velocity_scale *wl = &g[gno].wl;
    velocity_scale *fl = &g[gno].fl;
    map_scale *ms = &g[gno].mapscale;
    north_indicator *no = &g[gno].north;
    switch (type) {
    case TIDALCLOCK:		/* tidal_clock */
	tc->active = OFF;
	tc->type = 0;
	tc->p = d_props;
	tc->p.fillcol = 0;
	tc->loctype = WORLD;
	tc->x = 0.;
	tc->y = 0.;
	tc->total_time = 12.4;
	break;
    case TIMELINE:		/* time_line */
	tl->active = OFF;
	tl->type = 0;
	tl->p = d_props;
	tl->c1 = 1;
	tl->c2 = 0;
	tl->c3 = 1;
	tl->len = 20;
	tl->loctype = WORLD;
	tl->x = 0.;
	tl->y = 0.;
	tl->start = 0.;
	tl->stop = 10.;
	tl->step = 1.;
	tl->units = 3;
	break;
    case TIMEINFO:		/* time_info */
	ts->active = OFF;
	ts->loctype = WORLD;
	strcpy(ts->format, "%Y-%b-%d %H:%M:%S");
	ts->x = 0.;
	ts->y = 0.;
	ts->lines = 1;
	ts->linew = 1;
	ts->color = 1;
	ts->rot = 0;
	ts->font = 4;
	ts->just = 0;
	ts->charsize = 1.0;
	strcpy(ts->start, "20011218000000");
	ts->s[0] = 0;
	break;
    case VSCALE:		/* velocity_legend */
	vl->active = OFF;
	vl->type = 0;
	vl->p = d_props;
	vl->loctype = WORLD;
	vl->x = 0.;
	vl->y = 0.;
	vl->units = M;
	vl->unitfac = 1.0;
	vl->len = 1.0;
	vl->scale = 10.0;
	break;
    case WSCALE:		/* wind_legend */
	wl->active = OFF;
	wl->type = 0;
	wl->p = d_props;
	wl->loctype = WORLD;
	wl->x = 0.;
	wl->y = 0.;
	wl->units = M;
	wl->unitfac = 1.0;
	wl->len = 1.0;
	wl->scale = 10.0;
	break;
    case FLUX:			/* flux_legend */
	fl->active = OFF;
	fl->type = 0;
	fl->p = d_props;
	fl->loctype = WORLD;
	fl->x = 0.;
	fl->y = 0.;
	fl->units = M;
	fl->unitfac = 1.0;
	fl->len = 1.0;
	fl->scale = 10.0;
	break;
    case MAPSCALE:		/* map_scale */
	ms->active = OFF;
	ms->type = 0;
	ms->p = d_props;
	ms->loctype = WORLD;
	ms->x = 0.;
	ms->y = 0.;
	ms->units = KM;
	ms->unitfac = 1000.0;
	ms->len = 10.0;
	break;
    case NORTH:		/* North indicator */
	no->active = OFF;
	no->type = 0;
	no->p = d_props;
	no->loctype = WORLD;
	no->x = 0.;
	no->y = 0.;
	break;
    }
}

void set_default_data(int type, int i)
{
    switch (type) {
    case GRID:
	grid[i].active = OFF;
	break;
    }
}

void set_region_defaults(int rno)
{
    int j;
    rg[rno].active = OFF;
    rg[rno].type = 0;
    rg[rno].color = 1;
    rg[rno].lines = 1;
    rg[rno].linew = 1;
    for (j = 0; j < MAXGRAPH; j++) {
	rg[rno].linkto[j] = -1;
    }
    rg[rno].n = 0;
    rg[rno].x = rg[rno].y = (double *) NULL;
    rg[rno].x1 = rg[rno].y1 = rg[rno].x2 = rg[rno].y2 = 0.0;
}

void set_default_framep(framep * f)
{
    memcpy(f, &d_f, sizeof(framep));
}

void set_default_world(world * w)
{
    memcpy(w, &d_w, sizeof(world));
}

void set_default_view(view * v)
{
    memcpy(v, &d_v, sizeof(view));
}

void set_default_string(plotstr * s)
{
    memcpy(s, &d_s, sizeof(plotstr));
}

void set_default_line(linetype * l)
{
    memcpy(l, &d_l, sizeof(linetype));
}

void set_default_box(boxtype * b)
{
    memcpy(b, &d_b, sizeof(boxtype));
}

void set_default_transect(Transect * t)
{
    t->display = OFF;
    t->start = -1;
    t->stop = 0;
    t->skip = 1;
}

void set_default_displaytrans(DisplayTransect * t)
{
    t->display = OFF;
    t->display_mag = OFF;
    t->p = d_props;
    set_isolines_defaults(0, &(t->ip));
}

void set_default_displaygrid(int gno, DisplayGrid * g)
{
    set_isolines_defaults(gno, &(g->ip));
    g->p = d_props;
    g->bp = d_props;
    g->display_boundary = OFF;
    g->display_nodes = OFF;
    g->display_elements = OFF;
    g->display_depths = OFF;
    g->display_courant = OFF;
    g->display_courantn = OFF;
    g->display_dimw = OFF;
    g->display_gridf = OFF;
    g->display_flags[EDIT_GRID] = 0;
    g->display_flags[EDIT_GRID_ISOLINES] = 0;
    g->display_flags[EDIT_BOUNDARY] = 1;
    g->display_flags[EDIT_GRID_FILLED] = 0;
    g->display_flags[EDIT_GRID_NODE_NUMBERS] = 0;
    g->display_flags[EDIT_GRID_ELEMENT_NUMBERS] = 0;
    g->display_flags[EDIT_GRID_DEPTHS] = 0;
    g->courantdt = 90;
}

void set_default_graph(int gno)
{
    int i;
    char buf[256];

    g[gno].active = OFF;
    g[gno].hidden = FALSE;
    g[gno].label = OFF;
    g[gno].type = XYFIXED;
    g[gno].auto_type = AUTO;
    g[gno].revx = FALSE;
    g[gno].revy = FALSE;
    g[gno].ws_top = 0;
    g[gno].dsx = g[gno].dsy = 0.0;	/* locator props */
    g[gno].pointset = FALSE;
    g[gno].pt_type = 0;
    g[gno].fx = GENERAL;
    g[gno].fy = GENERAL;
    g[gno].px = 6;
    g[gno].py = 6;
    set_default_ticks(&g[gno].t[0], X_AXIS);
    set_default_ticks(&g[gno].t[1], Y_AXIS);
    set_default_ticks(&g[gno].t[2], ZX_AXIS);
    set_default_ticks(&g[gno].t[3], ZY_AXIS);
    set_default_ticks(&g[gno].t[4], XA_AXIS);
    set_default_ticks(&g[gno].t[5], YA_AXIS);
    set_default_framep(&g[gno].f);
    set_default_world(&g[gno].w);
    set_default_view(&g[gno].v);
    set_default_string(&g[gno].labs.title);
    g[gno].labs.title.charsize = 1.5;
    set_default_string(&g[gno].labs.stitle);
    g[gno].labs.stitle.charsize = 1.0;
    for (i = 0; i < MAXGRIDTS; i++) {
	set_default_displaygrid(gno, (DisplayGrid *) & g[gno].gridt[i]);
    }
    for (i = 0; i < MAXTRANS; i++) {
	set_default_displaytrans((DisplayTransect *) &g[gno].trans[i]);
    }
    for (i = 0; i < MAXGRIDS; i++) {
	set_default_displaygrid(gno, &g[gno].grid[i]);
    }
    g[gno].grid[0].display_flags[EDIT_GRID_ISOLINES] = 0;
    for (i = 0; i < MAXTEANL; i++) {
	set_isolines_defaults(gno, &g[gno].flowf[i].elevip);
	set_isolines_defaults(gno, &g[gno].flowf[i].ampip);
	set_isolines_defaults(gno, &g[gno].flowf[i].phaseip);
	set_isolines_defaults(gno, &g[gno].flowf[i].magip);
	g[gno].flowf[i].display = OFF;
	g[gno].flowf[i].display_elevmarkers = ON;
	g[gno].flowf[i].display_elev = OFF;
	g[gno].flowf[i].display_elevdepth = OFF;
	g[gno].flowf[i].display_amp = OFF;
	g[gno].flowf[i].display_phase = OFF;
	g[gno].flowf[i].display_mag = OFF;
	g[gno].flowf[i].display_inun = OFF;
	g[gno].flowf[i].display_irestrict = OFF;
	g[gno].flowf[i].display_wind = OFF;
	g[gno].flowf[i].sample = OFF;
	g[gno].flowf[i].samptype = NODE;
	g[gno].flowf[i].nsamples = 0;
	g[gno].flowf[i].samples = NULL;
	g[gno].flowf[i].sampx = NULL;
	g[gno].flowf[i].sampy = NULL;
	g[gno].flowf[i].flowfreq = ALL;
	g[gno].flowf[i].p = d_props;
	g[gno].flowf[i].wet = d_props;
	g[gno].flowf[i].dry = d_props;
	g[gno].flowf[i].wetdry = d_props;
	set_default_elevmarkers(g[gno].flowf[i].em, MAXELEVMARKERS);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	set_isolines_defaults(gno, &g[gno].flowt[i].elevip);
	set_isolines_defaults(gno, &g[gno].flowt[i].maxelevip);
	set_isolines_defaults(gno, &g[gno].flowt[i].ampip);
	set_isolines_defaults(gno, &g[gno].flowt[i].phaseip);
	set_isolines_defaults(gno, &g[gno].flowt[i].magip);
	g[gno].flowt[i].display = OFF;
	g[gno].flowt[i].display_elev = OFF;
	g[gno].flowt[i].display_maxelev = OFF;
	g[gno].flowt[i].display_maxelevval = OFF;
	g[gno].flowt[i].display_elevdepth = OFF;
	g[gno].flowt[i].display_elevmarkers = ON;
	g[gno].flowt[i].display_amp = OFF;
	g[gno].flowt[i].display_phase = OFF;
	g[gno].flowt[i].display_mag = OFF;
	g[gno].flowt[i].display_wind = OFF;
	g[gno].flowt[i].display_inun = OFF;
	g[gno].flowt[i].display_irestrict = OFF;
	g[gno].flowt[i].sample = OFF;
	g[gno].flowt[i].nsamples = 0;
	g[gno].flowt[i].samples = NULL;
	g[gno].flowt[i].p = d_props;
	set_default_elevmarkers(g[gno].flowt[i].em, MAXELEVMARKERS);
    }
    for (i = 0; i < MAXELA; i++) {
	set_isolines_defaults(gno, &g[gno].elaconc[i].ip);
	g[gno].elaconc[i].display = OFF;
	g[gno].elaconc[i].display_isolines = 0;
	g[gno].elaconc[i].display_nodes = 0;
	g[gno].elaconc[i].display_max = 0;
	g[gno].elaconc[i].display_min = 0;
	g[gno].elaconc[i].p = d_props;
    }
    for (i = 0; i < MAXPATHLINES; i++) {
	g[gno].drogues[i].display = OFF;
	g[gno].drogues[i].display_id = OFF;
	g[gno].drogues[i].display_streaml = OFF;
	g[gno].drogues[i].display_type = OFF;
	g[gno].drogues[i].display_connect = OFF;
	g[gno].drogues[i].p = d_props;
	g[gno].drogues[i].p.color = 1;
	g[gno].drogues[i].p.lines = 1;
	g[gno].drogues[i].p.linew = 1;
	g[gno].drogues[i].p.symbol = 2;
	g[gno].drogues[i].p.symsize = 1;
    }
    for (i = 0; i < MAXTRACK; i++) {
	g[gno].track[i].display = OFF;
	g[gno].track[i].display_streaml = OFF;
	g[gno].track[i].display_type = OFF;
	g[gno].track[i].display_data = OFF;
	g[gno].track[i].display_isolines = OFF;
	g[gno].track[i].display_interp = OFF;
	g[gno].track[i].display_connect = OFF;
	g[gno].track[i].symsize = 1.0;
	g[gno].track[i].p = d_props;
    }
    for (i = 0; i < MAXVELHIST; i++) {
	g[gno].flowh[i].display = OFF;
	g[gno].flowh[i].p = d_props;
	g[gno].flowh[i].circle = 0;
	g[gno].flowh[i].cp = d_props;
    }
    for (i = 0; i < MAXBOUNDS; i++) {
	g[gno].bounds[i].display = OFF;
	g[gno].bounds[i].p = d_props;
    }
    set_default_sliceboxes(gno);
    set_default_fluxboxes(gno);
    set_default_zoomboxes(gno);
    set_default_histmarkers(gno);
    set_default_flow3d(gno);
    set_vis_defaults(gno, TIMEINFO);
    set_vis_defaults(gno, TIMELINE);
    set_vis_defaults(gno, TIDALCLOCK);
    set_vis_defaults(gno, VSCALE);
    set_vis_defaults(gno, WSCALE);
    set_vis_defaults(gno, FLUX);
    set_vis_defaults(gno, MAPSCALE);
    set_vis_defaults(gno, NORTH);
    g[gno].display_flags[EDIT_GRID] = 0;
    g[gno].display_flags[EDIT_GRID_ISOLINES] = 0;
    g[gno].display_flags[EDIT_BOUNDARY] = 1;
    g[gno].display_flags[EDIT_GRID_FILLED] = 0;
    g[gno].display_flags[EDIT_GRID_NODE_NUMBERS] = 0;
    g[gno].display_flags[EDIT_GRID_ELEMENT_NUMBERS] = 0;
    g[gno].display_flags[EDIT_GRID_DEPTHS] = 0;
/*
 * initialize tide station display
 */
    for (i = 0; i < 100; i++) {
	InitDisplayTideStation(&g[gno].tidestat[i]);
    }
/*
    g[gno].tidestat = NULL;
*/
    g[gno].sta = (Station *) NULL;
    nsta = 0;
    set_isolines_defaults(gno, &g[gno].velmagip);
    set_isolines_defaults(gno, &g[gno].salip);
    set_isolines_defaults(gno, &g[gno].curip);
}

void realloc_graphs(void)
{
    int i, j;

    g = (graph *) realloc(g, maxgraph * sizeof(graph));
    for (j = MAXGRAPH; j < maxgraph; j++) {
	set_default_graph(j);
    }
}

void set_default_annotation(void)
{
    int i;

    for (i = 0; i < MAXBOXES; i++) {
	set_default_box(&boxes[i]);
    }
    for (i = 0; i < MAXLINES; i++) {
	set_default_line(&lines[i]);
    }
    for (i = 0; i < MAXSTR; i++) {
	set_default_string(&pstr[i]);
    }
}

void set_default_ticks(tickmarks * t, int a)
{
    int i;

    t->axis = a;
    switch (a) {
    case X_AXIS:
    case Y_AXIS:
	t->active = ON;
	t->alt = OFF;
	t->tl_flag = OFF;
	t->mtl_flag = OFF;
	t->t_flag = OFF;
	break;
    case XA_AXIS:
    case YA_AXIS:
	t->active = ON;
	t->alt = OFF;
	t->tl_flag = OFF;
	t->mtl_flag = OFF;
	t->t_flag = OFF;
	break;
    case ZX_AXIS:
    case ZY_AXIS:
	t->active = ON;
	t->alt = OFF;
	t->tl_flag = OFF;
	t->mtl_flag = OFF;
	t->t_flag = OFF;
	break;
    }
    set_default_string(&t->label);
    t->tmin = 0.0;
    t->tmax = 1.0;
    t->tmajor = 0.5;
    t->tminor = 0.25;
    t->offsx = 0.0;
    t->offsy = 0.0;
    t->label_layout = PARA;
    t->label_place = AUTO;
    t->tl_type = AUTO;
    t->tl_layout = HORIZONTAL;
    t->tl_sign = NORMAL;
    t->tl_prec = 1;
    t->tl_format = DECIMAL;
    t->tl_angle = 0;
    t->tl_just = (a % 2) ? RIGHT : CENTER;
    t->tl_skip = 0;
    t->tl_staggered = 0;
    t->tl_starttype = AUTO;
    t->tl_stoptype = AUTO;
    t->tl_start = 0.0;
    t->tl_stop = 0.0;
    t->tl_op = (a % 2) ? LEFT : BOTTOM;
    t->tl_vgap = 1.0;
    t->tl_hgap = 1.0;
    t->tl_font = 4;
    t->tl_charsize = 1.0;
    t->tl_color = 1;
    t->tl_linew = 1;
    t->tl_appstr[0] = 0;
    t->tl_prestr[0] = 0;
    t->t_color = 1;
    t->t_linew = 1;
    t->t_type = AUTO;
    t->t_mflag = ON;
    t->t_integer = OFF;
    t->t_num = 6;
    t->t_inout = IN;
    t->t_log = OFF;
    t->t_op = BOTH;
    t->t_size = 1.0;
    t->t_msize = 0.5;
    t->t_drawbar = OFF;
    t->t_drawbarcolor = 1;
    t->t_drawbarlines = 1;
    t->t_drawbarlinew = 1;
    t->t_gridflag = OFF;
    t->t_mgridflag = OFF;
    t->t_color = 1;
    t->t_lines = 1;
    t->t_linew = 1;
    t->t_mcolor = 1;
    t->t_mlines = 1;
    t->t_mlinew = 1;
    t->t_spec = 0;
    for (i = 0; i < MAX_TICK_LABELS; i++) {
	t->t_specloc[i] = 0.0;
	t->t_speclab[i].s[0] = '\0';
    }
}
