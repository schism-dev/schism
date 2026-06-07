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
 * draw the teanl flow field
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawflow.c,v 1.3 2004/01/03 18:36:43 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"
#include <math.h>

double eval_teanl_node();

#define FAC 180.0/M_PI

extern int mapisolconc[];

extern double unitfac;

int is_elementwetdry(int gno, int elno, int flowno, double t);
int is_nodewetdry(int gno, int n, int flowno, double t);

void ReadSampleFlow(DisplayFlow * f, const char *s);
void WriteSampleFlow(DisplayFlow * f, const char *s);
void AddSampleFlowNode(DisplayFlow * f, int node);
void AddSampleFlowElem(DisplayFlow * f, int elem);
void AddSampleFlowXY(DisplayFlow * f, double x, double y);
void DeleteSampleFlowNode(DisplayFlow * f, int node);
void DeleteSampleFlowElem(DisplayFlow * f, int elem);
void DeleteSampleFlowXY(DisplayFlow * f, double x, double y);
int DisplaySample(DisplayFlow * f, int node);

void drawflowcenter(int gno, int flowno, double t)
{
    int i, j, ind, ii, n1, n2, n3, gridno = flowf[flowno].grid;
    double xx, yy, uu, vv, angx, angy, tang;
    int iangx, iangy, ir;
    double pcor = flowf[flowno].sign_of_phase;
    int jmin, jmax;
    double s;
    if (g[gno].flowf[flowno].display_wind == ON) {
	s = g[gno].wl.unitfac * g[gno].wl.scale;
    } else {
	s = g[gno].vl.unitfac * g[gno].vl.scale;
    }
    if (g[gno].flowf[flowno].flowfreq == ALL) {
	jmin = 0;
	jmax = flowf[flowno].nfreq - 1;
    } else {
	jmin = jmax = g[gno].flowf[flowno].flowfreq;
    }

    if (flowf[flowno].active == ON) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    if (g[gno].flowf[flowno].display_irestrict == ON) {
		int ir = is_elementwetdry(gno, i, flowno, t);
		if (ir == DRY) {
		    continue;
		}
	    }
	    uu = 0.0;
	    vv = 0.0;
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    get_center(gridno, i, &xx, &yy);
	    for (j = jmin; j <= jmax; j++) {
		tang = flowf[flowno].omega[j] * t;
		angx = tang - pcor * flowf[flowno].phax[j][n1];
		angy = tang - pcor * flowf[flowno].phay[j][n1];
		uu = uu + flowf[flowno].ampx[j][n1] * cos(angx);
		vv = vv + flowf[flowno].ampy[j][n1] * cos(angy);
		angx = tang - pcor * flowf[flowno].phax[j][n2];
		angy = tang - pcor * flowf[flowno].phay[j][n2];
		uu = uu + flowf[flowno].ampx[j][n2] * cos(angx);
		vv = vv + flowf[flowno].ampy[j][n2] * cos(angy);
		angx = tang - pcor * flowf[flowno].phax[j][n3];
		angy = tang - pcor * flowf[flowno].phay[j][n3];
		uu = uu + flowf[flowno].ampx[j][n3] * cos(angx);
		vv = vv + flowf[flowno].ampy[j][n3] * cos(angy);
	    }
	    velplt(xx, yy, uu * 0.333333333, vv * 0.333333333, s);
	}
    }
}

void drawflownodes(int gno, int flowno, double t)
{
    int i, j, n0, n1, n2, ind, ii;
    double xx, yy, uu, vv, angx, angy;
    int iangx, iangy, gridno;
    int jmin, jmax;
    double pcor = flowf[flowno].sign_of_phase;
    double s;
    if (g[gno].flowf[flowno].display_wind == ON) {
	s = g[gno].wl.unitfac * g[gno].wl.scale;
    } else {
	s = g[gno].vl.unitfac * g[gno].vl.scale;
    }
    if (g[gno].flowf[flowno].flowfreq == ALL) {
	jmin = 0;
	jmax = flowf[flowno].nfreq - 1;
    } else {
	jmin = jmax = g[gno].flowf[flowno].flowfreq;
    }

    if (flowf[flowno].active == ON) {
	gridno = flowf[flowno].grid;
	for (i = 0; i < grid[gridno].nmnp; i++) {
	    if (flowf[flowno].ampx[0][i] < 9999.0) {
		if (g[gno].flowf[flowno].sample == ON) {
		    if (!DisplaySample(&(g[gno].flowf[flowno]), i)) {
			continue;
		    }
		}
		if (g[gno].flowf[flowno].display_irestrict == ON) {
		    int ir = is_nodewetdry(gno, i, flowno, t);
		    if (ir == DRY) {
			continue;
		    }
		}
		uu = 0.0;
		vv = 0.0;
		for (j = jmin; j <= jmax; j++) {
		    angx = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phax[j][i];
		    angy = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phay[j][i];
		    uu = uu + flowf[flowno].ampx[j][i] * cos(angx);
		    vv = vv + flowf[flowno].ampy[j][i] * cos(angy);
		}
		velplt(grid[gridno].xord[i], grid[gridno].yord[i], uu, vv, s);
	    }
	}
    }
}

void drawflowmag(int gno, int flowno, double t)
{
    int i, j, n0, n1, n2, ind, ii;
    double xx, yy, uu, vv, angx, angy;
    int iangx, iangy, gridno;
    int jmin, jmax;
    double *mag;
    double pcor = flowf[flowno].sign_of_phase;
    if (g[gno].flowf[flowno].flowfreq == ALL) {
	jmin = 0;
	jmax = flowf[flowno].nfreq - 1;
    } else {
	jmin = jmax = g[gno].flowf[flowno].flowfreq;
    }

    if (flowf[flowno].active == ON) {
	gridno = flowf[flowno].grid;
	mag = (double *) calloc(grid[gridno].nmnp, sizeof(double));
	if (mag == NULL) {
	    return;
	}
	for (i = 0; i < grid[gridno].nmnp; i++) {
	    if (flowf[flowno].ampx[0][i] < 9999.00) {
		uu = 0.0;
		vv = 0.0;
		for (j = jmin; j <= jmax; j++) {
		    angx = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phax[j][i];
		    angy = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phay[j][i];
		    uu = uu + flowf[flowno].ampx[j][i] * cos(angx);
		    vv = vv + flowf[flowno].ampy[j][i] * cos(angy);
		}
		mag[i] = hypot(uu, vv);
	    }
	}
	do_isol(flowf[flowno].grid, 0, mag, g[gno].flowf[flowno].magip, mapisolconc, g[gno].flowf[flowno].magip.nisol);
	free(mag);
    }
}

void drawelevations(int gno, int flowno, double t)
{
    int j, gridno;
    int jmin, jmax;
    double *c;
    double pcor = flowf[flowno].sign_of_phase;
    if (g[gno].flowf[flowno].flowfreq == ALL) {
	jmin = 0;
	jmax = flowf[flowno].nfreq - 1;
    } else {
	jmin = jmax = g[gno].flowf[flowno].flowfreq;
    }

    if (flowf[flowno].active == ON) {
	gridno = flowf[flowno].grid;
	c = (double *) calloc(grid[gridno].nmnp, sizeof(double));
	if (c == NULL) {
	    return;
	}
	for (j = 0; j < grid[gridno].nmnp; j++) {
	    c[j] = eval_teanl_node(flowno, j, t);
	    if (g[gno].flowf[flowno].display_elevdepth == ON) {
		c[j] += grid[gridno].depth[j];
	    }
	}
	if (g[gno].flowf[flowno].display_irestrict == ON) {
	    for (j = 0; j < grid[gridno].nmel; j++) {
		int ir = is_elementwetdry(gno, j, flowno, t);
		if (ir == DRY) {
		    continue;
		}
		do_element_isol(&grid[gridno], j, 0, c, g[gno].flowf[flowno].elevip, mapisolconc, g[gno].flowf[flowno].elevip.nisol);
	    }
	} else {
	    do_isol(gridno, 0, c, g[gno].flowf[flowno].elevip, mapisolconc, g[gno].flowf[flowno].elevip.nisol);
	}

	free(c);
    }
}

void drawwetdry(int gno, int flowno, double t)
{
    int i, j, ind, ii, n1, n2, n3, gridno = flowf[flowno].grid;
    double e1, e2, e3, eval_teanl_node();
    if (flowf[flowno].active == ON) {
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    e1 = eval_teanl_node(flowno, n1, t) + grid[gridno].depth[n1] + flowf[flowno].dcor;
	    e2 = eval_teanl_node(flowno, n2, t) + grid[gridno].depth[n2] + flowf[flowno].dcor;
	    e3 = eval_teanl_node(flowno, n3, t) + grid[gridno].depth[n3] + flowf[flowno].dcor;
	    if ((e1 < 0.0) && (e2 < 0.0) && (e3 < 0.0)) {
		drawelement_filled(gridno, i, g[gno].flowf[flowno].dry.color);
	    } else if ((e1 > 0.0) && (e2 > 0.0) && (e3 > 0.0)) {
		drawelement_filled(gridno, i, g[gno].flowf[flowno].wet.color);
	    } else {
		drawelement_filled(gridno, i, g[gno].flowf[flowno].wetdry.color);
	    }
	}
    }
}

int is_elementwetdry(int gno, int elno, int flowno, double t)
{
    int ind, ii, n1, n2, n3, gridno = flowf[flowno].grid;
    double e1, e2, e3, eval_teanl_node();
    if (flowf[flowno].active == ON) {
	n1 = grid[gridno].icon[elno].nl[0];
	n2 = grid[gridno].icon[elno].nl[1];
	n3 = grid[gridno].icon[elno].nl[2];
	e1 = eval_teanl_node(flowno, n1, t) + grid[gridno].depth[n1] + flowf[flowno].dcor;
	e2 = eval_teanl_node(flowno, n2, t) + grid[gridno].depth[n2] + flowf[flowno].dcor;
	e3 = eval_teanl_node(flowno, n3, t) + grid[gridno].depth[n3] + flowf[flowno].dcor;
	if ((e1 < 0.0) && (e2 < 0.0) && (e3 < 0.0)) {
	    return DRY;
	} else if ((e1 > 0.0) && (e2 > 0.0) && (e3 > 0.0)) {
	    return WET;
	} else {
	    return WETDRY;
	}
    }
    return -1;
}

int is_pointwetdry(int gno, int elno, double x, double y, int flowno, double t)
{
    double eval_triangle3(int gridno, int el, double xp, double yp, double *c, int redo);
    int ind, ii, n1, n2, n3, gridno = flowf[flowno].grid;
    double dp, e[3], eval_teanl_node();
    if (flowf[flowno].active == ON) {
	n1 = grid[gridno].icon[elno].nl[0];
	n2 = grid[gridno].icon[elno].nl[1];
	n3 = grid[gridno].icon[elno].nl[2];
	e[0] = eval_teanl_node(flowno, n1, t) + grid[gridno].depth[n1] + flowf[flowno].dcor;
	e[1] = eval_teanl_node(flowno, n2, t) + grid[gridno].depth[n2] + flowf[flowno].dcor;
	e[2] = eval_teanl_node(flowno, n3, t) + grid[gridno].depth[n3] + flowf[flowno].dcor;
	dp = eval_triangle3(gridno, elno, x, y, e, 0);

	if (dp < 0.0) {
	    return DRY;
	} else if (dp > 0.0) {
	    return WET;
	} else {
	    return WETDRY;
	}
    }
    return -1;
}

int is_nodewetdry(int gno, int n, int flowno, double t)
{
    int gridno = flowf[flowno].grid;
    double e1, eval_teanl_node();
    if (flowf[flowno].active == ON) {
	e1 = eval_teanl_node(flowno, n, t) + grid[gridno].depth[n] + flowf[flowno].dcor;
	if (e1 < 0.0) {
	    return DRY;
	} else if (e1 > 0.0) {
	    return WET;
	} else {
	    return WETDRY;
	}
    }
    return 0;
}

void drawflowtmag(int gno, int flowno, int nstep)
{
    int i, j, ind, ii;
    double xx, yy, uu, vv;
    double *mag;
    if (flowt[flowno].active == ON) {
	mag = (double *) malloc(flowt[flowno].g.nmnp * sizeof(double));
	if (mag == NULL) {
	    return;
	}
	for (i = 0; i < flowt[flowno].g.nmnp; i++) {
	    uu = flowt[flowno].f[nstep].u[i];
	    vv = flowt[flowno].f[nstep].v[i];
	    mag[i] = hypot(uu, vv);
	}
	do_grid_isol(&flowt[flowno].g, 0, mag, g[gno].flowt[flowno].magip, mapisolconc, g[gno].flowt[flowno].magip.nisol);
	free(mag);
    }
}

void drawflowtcenter(int gno, int flowno, int n)
{
    int i, j, ind, ii, n1, n2, n3;
    double xx, yy, uu, vv;
    double s;
    if (g[gno].flowt[flowno].display_wind == ON) {
	s = g[gno].wl.unitfac * g[gno].wl.scale;
    } else {
	s = g[gno].vl.unitfac * g[gno].vl.scale;
    }

    if (flowt[flowno].active == ON) {
	for (i = 0; i < flowt[flowno].g.nmel; i++) {
	    n1 = flowt[flowno].g.icon[i].nl[0];
	    n2 = flowt[flowno].g.icon[i].nl[1];
	    n3 = flowt[flowno].g.icon[i].nl[2];
	    get_grid_center(&flowt[flowno].g, i, &xx, &yy);
	    uu = flowt[flowno].f[n].u[n1] + flowt[flowno].f[n].u[n2] + flowt[flowno].f[n].u[n3];
	    vv = flowt[flowno].f[n].v[n2] + flowt[flowno].f[n].v[n2] + flowt[flowno].f[n].v[n3];
	    velplt(xx, yy, uu * 0.333333333, vv * 0.333333333, s);
	}
    }
}

void drawflowtnodes(int gno, int flowno, int n)
{
    int i;
    double s;
    if (g[gno].flowf[flowno].display_wind == ON) {
	s = g[gno].wl.unitfac * g[gno].wl.scale;
    } else {
	s = g[gno].vl.unitfac * g[gno].vl.scale;
    }

    if (flowt[flowno].active == ON) {
	if (g[gno].flowt[flowno].sample == ON) {
	    if (g[gno].flowt[flowno].samples == NULL) {
		LoadSampleLocations(&(g[gno].flowt[flowno]), &(flowt[flowno].g));
	    }
	}
	for (i = 0; i < flowt[flowno].g.nmnp; i++) {
	    if (g[gno].flowt[flowno].sample == ON) {
		if (!DisplaySample(&(g[gno].flowt[flowno]), i)) {
		    continue;
		}
	    }
	    velplt(flowt[flowno].g.xord[i], flowt[flowno].g.yord[i], flowt[flowno].f[n].u[i], flowt[flowno].f[n].v[i], s);
	}
    }
}

int compute_adcirc_maxelev(int flowno, int recomp)
{
    double emax;
    int i, j;
    if (flowt[flowno].active == ON) {
	if (recomp) {
	    if (flowt[flowno].global_emax != NULL) {
		free(flowt[flowno].global_emax);
		flowt[flowno].global_emax = NULL;
	    }
	}
	if (flowt[flowno].global_emax == NULL) {
	    flowt[flowno].global_emax = (double *) calloc(flowt[flowno].g.nmnp, sizeof(double));
	    if (flowt[flowno].global_emax == NULL) {
		return 0;
	    }
	    for (i = 0; i < flowt[flowno].g.nmnp; i++) {
		emax = flowt[flowno].f[0].e[i];
		for (j = 0; j < flowt[flowno].nsteps; j++) {
		    if (flowt[flowno].f[j].e[i] > emax) {
			emax = flowt[flowno].f[j].e[i];
		    }
		}
		flowt[flowno].global_emax[i] = emax;
	    }
	}
    }
    return 1;
}

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

int compute_adcirc_stats(int flowno, double *emin, double *emax, int *imin, int *imax, double *temin, double *temax, int *itmin, int *itmax, double *gemin, double *gemax)
 /* min max at each time step */
 /* node at step for min/max */
 /* min max at each node */
 /* time step for min/max */
 /* global min max */
{
    int i, j;
    double *e1, *e2;
    if (flowt[flowno].active == ON) {
	for (j = 0; j < flowt[flowno].nsteps; i++) {
	    minmax2(flowt[flowno].f[j].e, flowt[flowno].g.nmnp, &emin[j], &emax[j], &imin[j], &imax[j]);
	    for (i = 0; i < flowt[flowno].g.nmnp; i++) {
		temax[i] = max(temax[i], flowt[flowno].f[j].e[i]);
		temin[i] = min(temin[i], flowt[flowno].f[j].e[i]);
	    }
	    *gemax = max(*gemax, emax[j]);
	    *gemin = min(*gemin, emin[j]);
	}
    }
    return 1;
}

void drawflowh(int gno, int flowno, int n)
{
    int i, j, n0, n1, n2, ind, ii;
    double xx, yy, uu, vv, angx, angy;
    int iangx, iangy;
    if (flowh[flowno].active == ON) {
	if (g[gno].flowh[flowno].circle) {
	    drawpolysym(flowh[flowno].x, flowh[flowno].y, flowh[flowno].npts, 2, 0, 2, g[gno].flowh[flowno].cp.symsize * 0.01);
	}
	for (i = 0; i < flowh[flowno].npts; i++) {
	    velplt(flowh[flowno].x[i], flowh[flowno].y[i], flowh[flowno].f[n].u[i], flowh[flowno].f[n].v[i], g[gno].vl.unitfac * g[gno].vl.scale);
	}
    }
}

void compute_teanl_node(int flowno, int node, double time, double *e, double *u, double *v)
{
    int i, j;
    double angx, angy, ange;
    double pcor = flowf[flowno].sign_of_phase;
    *e = *u = *v = 0.0;
    for (j = 0; j < flowf[flowno].nfreq; j++) {
	ange = flowf[flowno].omega[j] * time - pcor * flowf[flowno].elphase[j][node];
	angx = flowf[flowno].omega[j] * time - pcor * flowf[flowno].phax[j][node];
	angy = flowf[flowno].omega[j] * time - pcor * flowf[flowno].phay[j][node];
	*e = *e + flowf[flowno].elamp[j][node] * cos(ange);
	*u = *u + flowf[flowno].ampx[j][node] * cos(angx);
	*v = *v + flowf[flowno].ampy[j][node] * cos(angy);
    }
}

void drawflowt_ellipse(int gno, int flowno, int n)
{
    int i, j;
    double s;
    if (g[gno].flowf[flowno].display_wind == ON) {
	s = g[gno].wl.unitfac * g[gno].wl.scale;
    } else {
	s = g[gno].vl.unitfac * g[gno].vl.scale;
    }

    if (n <= 0) {		/* need 2 time steps */
	return;
    }
    if (flowt[flowno].active == ON) {
	for (i = 0; i < flowt[flowno].g.nmnp; i++) {
	    if (g[gno].flowt[flowno].sample == ON) {
		if (!DisplaySample(&g[gno].flowt[flowno], i)) {
		    continue;
		}
	    }
	    for (j = 0; j < n; j++) {
		vellineplt(flowt[flowno].g.xord[i], flowt[flowno].g.yord[i], flowt[flowno].f[j].u[i], flowt[flowno].f[j].v[i], flowt[flowno].f[j + 1].u[i], flowt[flowno].f[j + 1].v[i], s);
	    }
	}
    }
}

void drawflowf_ellipse(int gno, int flowno, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	drawflownodes(gno, flowno, timeclock.t[i]);
    }
}

int computeflow(int gno, int flowno, double t, int magflag, double *u, double *v)
{
    int i, j, n0, n1, n2, ind, ii;
    double xx, yy, uu, vv, angx, angy;
    int iangx, iangy, gridno;
    double *mag;
    double pcor = flowf[flowno].sign_of_phase;

    if (magflag && u == NULL) {
	return 1;
    }
    if (!magflag && u == NULL && v == NULL) {
	return 2;
    }
    if (flowf[flowno].active == ON) {
	gridno = flowf[flowno].grid;
	for (i = 0; i < grid[gridno].nmnp; i++) {
	    if (flowf[flowno].ampx[0][i] < 9999.00) {
		uu = 0.0;
		vv = 0.0;
		for (j = 0; j < flowf[flowno].nfreq; j++) {
		    angx = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phax[j][i];
		    angy = flowf[flowno].omega[j] * t - pcor * flowf[flowno].phay[j][i];
		    uu = uu + flowf[flowno].ampx[j][i] * cos(angx);
		    vv = vv + flowf[flowno].ampy[j][i] * cos(angy);
		}
		if (magflag) {
		    u[i] = hypot(uu, vv);
		} else {
		    u[i] = uu;
		    v[i] = vv;
		}
	    }
	}
    }
    return 0;
}
