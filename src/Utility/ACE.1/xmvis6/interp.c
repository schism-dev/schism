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
 * Interpolate various models and gauss integration
 *
 */

#ifndef lint
static char RCSid[] = "$Id: interp.c,v 1.3 2003/08/25 15:45:40 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

double get_depth_element(int gno, int ind, double x, double y);
double eval_teanl(int flowno, int el, double t, double xp, double yp, int redo);
int eval_teanl_flow(int flowno, int el, double t, double xp, double yp, double *uret, double *vret, int redo);
int eval_teanl_flownode(int flowno, int nn, double t, double *uret, double *vret);
double eval_teanl_node(int flowno, int node, double t);
double eval_adcirc(int flowno, int el, int step, double t, double xp, double yp, int redo);
int eval_adcirc_flow(int flowno, int el, int step, double t, double xp, double yp, double *uret, double *vret, int redo);
double eval_ela(int concno, int el, int step, double t, double xp, double yp, int redo);
double gaussi(int gridno, int nel, double *c, double *ht);
double ufun(int ftype, int flowno, int ind, double x, double y, int step, double t);
double vfun(int ftype, int flowno, int ind, double x, double y, int step, double t);
double wfun(void);
int belel2(void);
int get_grid_number(void);
double eval_model(int obj, int w1, int w2, int step, int ind, double t, double x, double y, int redo);
double eval_triangle3(int gridno, int el, double xp, double yp, double *c, int redo);
double eval_triangle6(int gridno, int el, double xp, double yp, double *c, int redo);
double eval_quadrangle4(int gridno, int el, double r, double s, double *c, int redo);
double eval_quadrangle8(int gridno, int el, double r, double s, double *c, int redo);

double eval_teanl(int flowno, int el, double t, double xp, double yp, int redo)
{
    int i;
    double e;
    static int n1, n2, n3;
    static double h1, h2, h3;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double atmp, btmp;
    static double phn[MAXFREQ], an[MAXFREQ];
    double *x, *y;
    double pcor = flowf[flowno].sign_of_phase;
    if (!redo) {
	x = grid[flowf[flowno].grid].xord;
	y = grid[flowf[flowno].grid].yord;
	n1 = grid[flowf[flowno].grid].icon[el].nl[0];
	n2 = grid[flowf[flowno].grid].icon[el].nl[1];
	n3 = grid[flowf[flowno].grid].icon[el].nl[2];

	aum = x[n2] * y[n3] - x[n3] * y[n2];
	bum = y[n2] - y[n3];
	cum = x[n3] - x[n2];
	ado = x[n3] * y[n1] - x[n1] * y[n3];
	bdo = y[n3] - y[n1];
	cdo = x[n1] - x[n3];
	atr = x[n1] * y[n2] - x[n2] * y[n1];
	btr = y[n1] - y[n2];
	ctr = x[n2] - x[n1];
	arei = 1.0 / (aum + ado + atr);
	h3 = (atr + btr * xp + ctr * yp) * arei;
	h2 = (ado + bdo * xp + cdo * yp) * arei;
	h1 = 1.0 - h2 - h3;
	for (i = 0; i < flowf[flowno].nfreq; i++) {
	    atmp = h1 * flowf[flowno].elamp[i][n1] * cos(flowf[flowno].elphase[i][n1]) + h2 * flowf[flowno].elamp[i][n2] * cos(flowf[flowno].elphase[i][n2]) + h3 * flowf[flowno].elamp[i][n3] * cos(flowf[flowno].elphase[i][n3]);
	    btmp = h1 * flowf[flowno].elamp[i][n1] * sin(flowf[flowno].elphase[i][n1]) + h2 * flowf[flowno].elamp[i][n2] * sin(flowf[flowno].elphase[i][n2]) + h3 * flowf[flowno].elamp[i][n3] * sin(flowf[flowno].elphase[i][n3]);
	    phn[i] = atan2(btmp, atmp);
	    an[i] = atmp / cos(phn[i]);
	}
    }
    e = 0.0;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	e += an[i] * cos(flowf[flowno].omega[i] * t - pcor * phn[i]);
    }
    return e;
}

int eval_teanl_flow(int flowno, int el, double t, double xp, double yp, double *uret, double *vret, int redo)
{
    int i;
    double e;
    static int n1, n2, n3;
    static double h1, h2, h3;
    static double phnx[MAXFREQ], anx[MAXFREQ];
    static double phny[MAXFREQ], any[MAXFREQ];
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double atmpx, btmpx;
    double atmpy, btmpy;
    double *x, *y;
    double pcor = flowf[flowno].sign_of_phase;
    if (!redo) {
	x = grid[flowf[flowno].grid].xord;
	y = grid[flowf[flowno].grid].yord;
	n1 = grid[flowf[flowno].grid].icon[el].nl[0];
	n2 = grid[flowf[flowno].grid].icon[el].nl[1];
	n3 = grid[flowf[flowno].grid].icon[el].nl[2];

	aum = x[n2] * y[n3] - x[n3] * y[n2];
	bum = y[n2] - y[n3];
	cum = x[n3] - x[n2];
	ado = x[n3] * y[n1] - x[n1] * y[n3];
	bdo = y[n3] - y[n1];
	cdo = x[n1] - x[n3];
	atr = x[n1] * y[n2] - x[n2] * y[n1];
	btr = y[n1] - y[n2];
	ctr = x[n2] - x[n1];
	arei = 1.0 / (aum + ado + atr);
	h3 = (atr + btr * xp + ctr * yp) * arei;
	h2 = (ado + bdo * xp + cdo * yp) * arei;
	h1 = 1.0 - h2 - h3;
	for (i = 0; i < flowf[flowno].nfreq; i++) {
	    atmpx = h1 * flowf[flowno].ampx[i][n1] * cos(flowf[flowno].phax[i][n1]) + h2 * flowf[flowno].ampx[i][n2] * cos(flowf[flowno].phax[i][n2]) + h3 * flowf[flowno].ampx[i][n3] * cos(flowf[flowno].phax[i][n3]);
	    btmpx = h1 * flowf[flowno].ampx[i][n1] * sin(flowf[flowno].phax[i][n1]) + h2 * flowf[flowno].ampx[i][n2] * sin(flowf[flowno].phax[i][n2]) + h3 * flowf[flowno].ampx[i][n3] * sin(flowf[flowno].phax[i][n3]);
	    phnx[i] = atan2(btmpx, atmpx);
	    anx[i] = atmpx / cos(phnx[i]);
	    atmpy = h1 * flowf[flowno].ampy[i][n1] * cos(flowf[flowno].phay[i][n1]) + h2 * flowf[flowno].ampy[i][n2] * cos(flowf[flowno].phay[i][n2]) + h3 * flowf[flowno].ampy[i][n3] * cos(flowf[flowno].phay[i][n3]);
	    btmpy = h1 * flowf[flowno].ampy[i][n1] * sin(flowf[flowno].phay[i][n1]) + h2 * flowf[flowno].ampy[i][n2] * sin(flowf[flowno].phay[i][n2]) + h3 * flowf[flowno].ampy[i][n3] * sin(flowf[flowno].phay[i][n3]);
	    phny[i] = atan2(btmpy, atmpy);
	    any[i] = atmpy / cos(phny[i]);
	}
    }
    *uret = 0.0;
    *vret = 0.0;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	*uret += anx[i] * cos(flowf[flowno].omega[i] * t - pcor * phnx[i]);
	*vret += any[i] * cos(flowf[flowno].omega[i] * t - pcor * phny[i]);
    }
    return 1;
}

int eval_teanl_flownode(int flowno, int nn, double t, double *uret, double *vret)
{
    int i;
    double e;
    double pcor = flowf[flowno].sign_of_phase;
    *uret = 0.0;
    *vret = 0.0;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	*uret += flowf[flowno].ampx[i][nn] * cos(flowf[flowno].omega[i] * t - pcor * flowf[flowno].phax[i][nn]);
	*vret += flowf[flowno].ampy[i][nn] * cos(flowf[flowno].omega[i] * t - pcor * flowf[flowno].phay[i][nn]);
    }
    return 1;
}

double eval_teanl_node(int flowno, int node, double t)
{
    int i;
    double e = 0.0;
    double pcor = flowf[flowno].sign_of_phase;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	e += flowf[flowno].elamp[i][node] * cos(flowf[flowno].omega[i] * t - pcor * flowf[flowno].elphase[i][node]);
    }
    return e;
}

double eval_adcirc(int flowno, int el, int step, double t, double xp, double yp, int redo)
{
    static double w[4];
    int i, n;
    double x[4], y[4], sum;
    if (!redo) {
	for (i=0;i<flowt[flowno].g.icon[el].nn;i++) {
	    n = flowt[flowno].g.icon[el].nl[i];
	    x[i] = flowt[flowno].g.xord[flowt[flowno].g.icon[el].nl[i]];
	    y[i] = flowt[flowno].g.yord[flowt[flowno].g.icon[el].nl[i]];
	    /* printf("%d,%d: %lf %lf\n", el, n, x[i], y[i]); */
	}
        ElioGetCoefficientsXY(flowt[flowno].g.icon[el].nn, x, y, xp, yp, w);
    }
    sum = 0.0;
    for (i=0;i<flowt[flowno].g.icon[el].nn;i++) {
	    n = flowt[flowno].g.icon[el].nl[i];
	    sum += w[i] * flowt[flowno].f[step].e[n];
	    /* printf("%d,%d: %f\n", el, n, flowt[flowno].f[step].e[n]); */
    }
    return sum;
}

int eval_adcirc_flow(int flowno, int el, int step, double t, double xp, double yp, double *uret, double *vret, int redo)
{
    static double w[4];
    int i, n;
    double x[4], y[4];
    if (!redo) {
	for (i=0;i<flowt[flowno].g.icon[el].nn;i++) {
	    n = flowt[flowno].g.icon[el].nl[i];
	    x[i] = flowt[flowno].g.xord[flowt[flowno].g.icon[el].nl[i]];
	    y[i] = flowt[flowno].g.yord[flowt[flowno].g.icon[el].nl[i]];
	}
        ElioGetCoefficientsXY(flowt[flowno].g.icon[el].nn, x, y, xp, yp, w);
    }
    *uret = 0.0;
    *vret = 0.0;
    for (i=0;i<flowt[flowno].g.icon[el].nn;i++) {
	    n = flowt[flowno].g.icon[el].nl[i];
	    *uret += w[i] * flowt[flowno].f[step].u[n];
	    *vret += w[i] * flowt[flowno].f[step].v[n];
    }
    return 1;
}

double eval_ela(int concno, int el, int step, double t, double xp, double yp, int redo)
{
    double e;
    static int i, n1, n2, n3, n4, n5, n6;
    static double h1, h2, h3, h4, h5, h6;
    static double c1, c2, c3, c4, c5, c6;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double *x, *y;
    if (!redo) {
	x = grid[elaconc[concno].grid].xord;
	y = grid[elaconc[concno].grid].yord;
	n1 = grid[elaconc[concno].grid].icon[el].nl[0];
	n2 = grid[elaconc[concno].grid].icon[el].nl[1];
	n3 = grid[elaconc[concno].grid].icon[el].nl[2];
	n4 = grid[elaconc[concno].grid].icon[el].nl[3];
	n5 = grid[elaconc[concno].grid].icon[el].nl[4];
	n6 = grid[elaconc[concno].grid].icon[el].nl[5];

	aum = x[n2] * y[n3] - x[n3] * y[n2];
	bum = y[n2] - y[n3];
	cum = x[n3] - x[n2];
	ado = x[n3] * y[n1] - x[n1] * y[n3];
	bdo = y[n3] - y[n1];
	cdo = x[n1] - x[n3];
	atr = x[n1] * y[n2] - x[n2] * y[n1];
	btr = y[n1] - y[n2];
	ctr = x[n2] - x[n1];
	arei = 1.0 / (aum + ado + atr);
	h3 = (atr + btr * xp + ctr * yp) * arei;
	h2 = (ado + bdo * xp + cdo * yp) * arei;
	h1 = 1.0 - h2 - h3;
	h4 = 4.0 * h2 * h1;
	h5 = 4.0 * h2 * h3;
	h6 = 4.0 * h3 * h1;
	h1 = h1 - 0.5 * (h4 + h6);
	h2 = h2 - 0.5 * (h4 + h5);
	h3 = h3 - 0.5 * (h5 + h6);
    }
    c1 = elaconc[concno].data[step].c[n1];
    c2 = elaconc[concno].data[step].c[n2];
    c3 = elaconc[concno].data[step].c[n3];
    c4 = elaconc[concno].data[step].c[n4];
    c5 = elaconc[concno].data[step].c[n5];
    c6 = elaconc[concno].data[step].c[n6];
    return h1 * c1 + h2 * c2 + h3 * c3 + h4 * c4 + h5 * c5 + h6 * c6;
}

/*
 * compute mass in an element (quadratic)
 *  using a 7 point gauss rule
 */

static double r3[7] = { .1012965073235, .7974269853531, .1012965073235,
    .4701420641051, .4701420641051, .0597158717898, .3333333333333,
};

static double s3[7] = {
    .1012965073235, .1012965073235, .7974269853531, .0597158717898,
    .4701420641051, .4701420641051, .3333333333333
};

static double w3[7] = {
    .1259391805448, .1259391805448, .1259391805448, .1323941527885,
    .1323941527885, .1323941527885, .225
};

double gaussi(int gridno, int nel, double *c, double *ht)
{
    double ret_val;

    int i, j, n[6];
    double p[6], r, s, h1, ch;
    ret_val = 0.0;
    for (i = 0; i < 7; ++i) {
	ch = 0.;
	r = r3[i];
	s = s3[i];
	h1 = 1.0 - r - s;
	p[3] = r * 4.0 * h1;
	p[4] = r * 4.0 * s;
	p[5] = s * 4.0 * h1;
	p[0] = h1 - (p[3] + p[5]) * 0.5;
	p[1] = r - (p[3] + p[4]) * 0.5;
	p[2] = s - (p[4] + p[5]) * 0.5;
	for (j = 0; j < 6; ++j) {
	    ch += c[j] * ht[j] * p[j];
	}
	ret_val += ch * w3[i];
    }
    return ret_val;
}

double ufun(int ftype, int flowno, int ind, double x, double y, int step, double t)
{
}

double vfun(int ftype, int flowno, int ind, double x, double y, int step, double t)
{
}

double wfun(void)
{
}

int belel2(void)
{
    return 0;
}

int get_grid_number(void)
{
    return 0;
}

double eval_model(int obj, int w1, int w2, int step, int ind, double t, double x, double y, int redo)
{
    double tmp, mag, u, v;
    switch (obj) {
    case BATH:
	return get_depth_element(w2, ind, x, y);
	break;
    case TEANL:
	if (w1 == ELEV)
	    return eval_teanl(w2, ind, t, x, y, redo);
	if (w1 == MAG) {
	    eval_teanl_flow(w2, ind, t, x, y, &u, &v, redo);
	    mag = hypot(u, v);
	    return mag;
	}
	break;
    case ADCIRC:
	if (w1 == ELEV)
	    return eval_adcirc(w2, ind, step, t, x, y, redo);
	if (w1 == MAG) {
	    eval_adcirc_flow(w2, ind, step, t, x, y, &u, &v, redo);
	    mag = hypot(u, v);
	    return mag;
	}
	break;
    case ADCIRC3DFLOW:
	break;
    case ELA:
	return eval_ela(w2, ind, step, t, x, y, redo);
	break;
    }
}

double eval_triangle3(int gridno, int el, double xp, double yp, double *c, int redo)
{
    static int n1, n2, n3;
    static double h1, h2, h3;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr;
    double c1, c2, c3, arei;
    double *x, *y;
    if (!redo) {
	x = grid[gridno].xord;
	y = grid[gridno].yord;
	n1 = grid[gridno].icon[el].nl[0];
	n2 = grid[gridno].icon[el].nl[1];
	n3 = grid[gridno].icon[el].nl[2];

	aum = x[n2] * y[n3] - x[n3] * y[n2];
	bum = y[n2] - y[n3];
	cum = x[n3] - x[n2];
	ado = x[n3] * y[n1] - x[n1] * y[n3];
	bdo = y[n3] - y[n1];
	cdo = x[n1] - x[n3];
	atr = x[n1] * y[n2] - x[n2] * y[n1];
	btr = y[n1] - y[n2];
	ctr = x[n2] - x[n1];
	arei = 1.0 / (aum + ado + atr);
	h3 = (atr + btr * xp + ctr * yp) * arei;
	h2 = (ado + bdo * xp + cdo * yp) * arei;
	h1 = 1.0 - h2 - h3;
    }
    return (h1 * c[0] + h2 * c[1] + h3 * c[2]);
}

double eval_triangle6(int gridno, int el, double xp, double yp, double *c, int redo)
{
    double e;
    static int i, n1, n2, n3, n4, n5, n6;
    static double h1, h2, h3, h4, h5, h6;
    static double c1, c2, c3, c4, c5, c6;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double *x, *y;
    if (!redo) {
	x = grid[gridno].xord;
	y = grid[gridno].yord;
	n1 = grid[gridno].icon[el].nl[0];
	n2 = grid[gridno].icon[el].nl[1];
	n3 = grid[gridno].icon[el].nl[2];
	n4 = grid[gridno].icon[el].nl[3];
	n5 = grid[gridno].icon[el].nl[4];
	n6 = grid[gridno].icon[el].nl[5];

	aum = x[n2] * y[n3] - x[n3] * y[n2];
	bum = y[n2] - y[n3];
	cum = x[n3] - x[n2];
	ado = x[n3] * y[n1] - x[n1] * y[n3];
	bdo = y[n3] - y[n1];
	cdo = x[n1] - x[n3];
	atr = x[n1] * y[n2] - x[n2] * y[n1];
	btr = y[n1] - y[n2];
	ctr = x[n2] - x[n1];
	arei = 1.0 / (aum + ado + atr);
	h3 = (atr + btr * xp + ctr * yp) * arei;
	h2 = (ado + bdo * xp + cdo * yp) * arei;
	h1 = 1.0 - h2 - h3;
	h4 = 4.0 * h2 * h1;
	h5 = 4.0 * h2 * h3;
	h6 = 4.0 * h3 * h1;
	h1 = h1 - 0.5 * (h4 + h6);
	h2 = h2 - 0.5 * (h4 + h5);
	h3 = h3 - 0.5 * (h5 + h6);
    }
    return h1 * c[0] + h2 * c[1] + h3 * c[2] + h4 * c[3] + h5 * c[4] + h6 * c[5];
}

/*
 * eval r, s in a 4 node quadrangle
 */
double eval_quadrangle4(int gridno, int el, double r, double s, double *c, int redo)
{
    double h1, h2, h3, h4;
    h1 = 0.25 * (1.0 + r) * (1.0 + s);
    h2 = 0.25 * (1.0 - r) * (1.0 + s);
    h3 = 0.25 * (1.0 - r) * (1.0 - s);
    h4 = 0.25 * (1.0 + r) * (1.0 - s);
    return (h1 * c[0] + h2 * c[1] + h3 * c[2] + h4 * c[3]);
}

/*
 * eval r, s in an 8 node quadrangle
 */
double eval_quadrangle8(int gridno, int el, double r, double s, double *c, int redo)
{
    int i;
    double h[8], sum = 0.0;
    h[4] = 0.5 * (1.0 + r * r) * (1.0 + s);
    h[5] = 0.5 * (1.0 - r) * (1.0 + s * s);
    h[6] = 0.5 * (1.0 - r * r) * (1.0 - s);
    h[7] = 0.5 * (1.0 + r) * (1.0 - s * s);
    h[0] = 0.25 * (1.0 + r) * (1.0 + s) - 0.5 * h[4] - 0.5 * h[7];
    h[1] = 0.25 * (1.0 - r) * (1.0 + s) - 0.5 * h[4] - 0.5 * h[5];
    h[2] = 0.25 * (1.0 - r) * (1.0 - s) - 0.5 * h[5] - 0.5 * h[6];
    h[3] = 0.25 * (1.0 + r) * (1.0 - s) - 0.5 * h[6] - 0.5 * h[7];
    for (i = 0; i < 8; i++) {
	sum += c[i] * h[i];
    }
    return sum;
}

void GetTeanlElevAmpPhase(int flowno, int el, double xp, double yp, double *amp, double *pha)
{
    int i;
    double e;
    static int n1, n2, n3;
    static double h1, h2, h3;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double atmp, btmp;
    double *x, *y;
    double pcor = flowf[flowno].sign_of_phase;
    x = grid[flowf[flowno].grid].xord;
    y = grid[flowf[flowno].grid].yord;
    n1 = grid[flowf[flowno].grid].icon[el].nl[0];
    n2 = grid[flowf[flowno].grid].icon[el].nl[1];
    n3 = grid[flowf[flowno].grid].icon[el].nl[2];

    aum = x[n2] * y[n3] - x[n3] * y[n2];
    bum = y[n2] - y[n3];
    cum = x[n3] - x[n2];
    ado = x[n3] * y[n1] - x[n1] * y[n3];
    bdo = y[n3] - y[n1];
    cdo = x[n1] - x[n3];
    atr = x[n1] * y[n2] - x[n2] * y[n1];
    btr = y[n1] - y[n2];
    ctr = x[n2] - x[n1];
    arei = 1.0 / (aum + ado + atr);
    h3 = (atr + btr * xp + ctr * yp) * arei;
    h2 = (ado + bdo * xp + cdo * yp) * arei;
    h1 = 1.0 - h2 - h3;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	atmp = h1 * flowf[flowno].elamp[i][n1] * cos(flowf[flowno].elphase[i][n1]) + h2 * flowf[flowno].elamp[i][n2] * cos(flowf[flowno].elphase[i][n2]) + h3 * flowf[flowno].elamp[i][n3] * cos(flowf[flowno].elphase[i][n3]);
	btmp = h1 * flowf[flowno].elamp[i][n1] * sin(flowf[flowno].elphase[i][n1]) + h2 * flowf[flowno].elamp[i][n2] * sin(flowf[flowno].elphase[i][n2]) + h3 * flowf[flowno].elamp[i][n3] * sin(flowf[flowno].elphase[i][n3]);
	pha[i] = atan2(btmp, atmp);
	amp[i] = atmp / cos(pha[i]);
    }
}

void GetTeanlFlowAmpPhase(int flowno, int el, double xp, double yp, double *ampx, double *phax, double *ampy, double *phay)
{
    int i;
    double e;
    static int n1, n2, n3;
    static double h1, h2, h3;
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, arei;
    double atmpx, btmpx;
    double atmpy, btmpy;
    double *x, *y;
    double pcor = flowf[flowno].sign_of_phase;
    x = grid[flowf[flowno].grid].xord;
    y = grid[flowf[flowno].grid].yord;
    n1 = grid[flowf[flowno].grid].icon[el].nl[0];
    n2 = grid[flowf[flowno].grid].icon[el].nl[1];
    n3 = grid[flowf[flowno].grid].icon[el].nl[2];

    aum = x[n2] * y[n3] - x[n3] * y[n2];
    bum = y[n2] - y[n3];
    cum = x[n3] - x[n2];
    ado = x[n3] * y[n1] - x[n1] * y[n3];
    bdo = y[n3] - y[n1];
    cdo = x[n1] - x[n3];
    atr = x[n1] * y[n2] - x[n2] * y[n1];
    btr = y[n1] - y[n2];
    ctr = x[n2] - x[n1];
    arei = 1.0 / (aum + ado + atr);
    h3 = (atr + btr * xp + ctr * yp) * arei;
    h2 = (ado + bdo * xp + cdo * yp) * arei;
    h1 = 1.0 - h2 - h3;
    for (i = 0; i < flowf[flowno].nfreq; i++) {
	atmpx = h1 * flowf[flowno].ampx[i][n1] * cos(flowf[flowno].phax[i][n1]) + h2 * flowf[flowno].ampx[i][n2] * cos(flowf[flowno].phax[i][n2]) + h3 * flowf[flowno].ampx[i][n3] * cos(flowf[flowno].phax[i][n3]);
	btmpx = h1 * flowf[flowno].ampx[i][n1] * sin(flowf[flowno].phax[i][n1]) + h2 * flowf[flowno].ampx[i][n2] * sin(flowf[flowno].phax[i][n2]) + h3 * flowf[flowno].ampx[i][n3] * sin(flowf[flowno].phax[i][n3]);
	phax[i] = atan2(btmpx, atmpx);
	ampx[i] = atmpx / cos(phax[i]);
	atmpy = h1 * flowf[flowno].ampy[i][n1] * cos(flowf[flowno].phay[i][n1]) + h2 * flowf[flowno].ampy[i][n2] * cos(flowf[flowno].phay[i][n2]) + h3 * flowf[flowno].ampy[i][n3] * cos(flowf[flowno].phay[i][n3]);
	btmpy = h1 * flowf[flowno].ampy[i][n1] * sin(flowf[flowno].phay[i][n1]) + h2 * flowf[flowno].ampy[i][n2] * sin(flowf[flowno].phay[i][n2]) + h3 * flowf[flowno].ampy[i][n3] * sin(flowf[flowno].phay[i][n3]);
	phay[i] = atan2(btmpy, atmpy);
	ampy[i] = atmpy / cos(phay[i]);
    }
}
