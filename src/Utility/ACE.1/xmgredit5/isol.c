/*
 * ACE/gredit - 2d finite element grid generation
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *                      All Rights Reserved.
 */

/*
 * draw isolines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isol.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include "defines.h"
#include "globals.h"

#define QUADRATIC 0
#define LINEAR 1
#define TEANL_AMP 3
#define TEANL_PHA 4

/* isol.c */
void fillel(int el, double e, struct isolparms p, int *cmap, int nmap);
void fillel2(int el, double e, struct isolparms p, int *cmap, int nmap);
void do_isol(int gridno, int itype, double *q, struct isolparms p, int *cmap, int nmap);
void do_grid_isol(Grid *gr, int itype, double *q, struct isolparms p, int *cmap, int nmap);
int checkonl(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int *imin, int *imax);
int checkonq(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, double c4, double c5, double c6, int *imin, int *imax);
int inrange(double c, double c1, double c2);
void drawisopat(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int imin, int imax);
void drawiso(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int imin, int imax);
void split(double x1, double y1, double x2, double y2, double x3, double y3, int ilev, int itype, int ipen, int imin, int imax);
void concna(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3);
void concna2(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3);
void comp_amp(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3);
void comp_pha(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3);

static double cc1, cc2, cc3, cc4, cc5, cc6;
static double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr, are2;
static struct isolparms ip;
static int static_cmap[16];

void do_grid_isol(Grid * gr, int itype, double *q, struct isolparms p, int *cmap, int nmap);

extern int mapisolbath[];
extern int cancel_flag;

FILE *sfp;			/* write out an isoline level */
int save_isoline = -1;
char saveisol_fname[256];

static int labcnt = 0;

void fillel(int el, double e, struct isolparms p, int *cmap, int nmap)
{
    int i, k, n1, n2, n3;
    double *x = grid[0].xord, *y = grid[0].yord, cc, ccp1, px[4], py[4];
    n1 = grid[0].icon[el].nl[0];
    n2 = grid[0].icon[el].nl[1];
    n3 = grid[0].icon[el].nl[2];
    px[0] = x[n1];
    py[0] = y[n1];
    px[1] = x[n2];
    py[1] = y[n2];
    px[2] = x[n3];
    py[2] = y[n3];
    for (k = 0; k < nmap - 1; k++) {
	cc = p.cis[k];
	ccp1 = p.cis[k + 1];
	setcolor(cmap[k]);
	if (e > cc && e <= ccp1) {
	    fillcolor(3, px, py);
	}
    }
}

void fillel2(int el, double e, struct isolparms p, int *cmap, int nmap)
{
    int i, k, n1, n2, n3, n4;
    double *x = grid[0].xord, *y = grid[0].yord, cc, ccp1, px[4], py[4];
    n1 = grid[0].icon[el].nl[0];
    n2 = grid[0].icon[el].nl[1];
    n3 = grid[0].icon[el].nl[2];
    n4 = grid[0].icon[el].nl[3];
    px[0] = x[n1];
    py[0] = y[n1];
    px[1] = x[n2];
    py[1] = y[n2];
    px[2] = x[n3];
    py[2] = y[n3];
    px[3] = x[n4];
    py[3] = y[n4];
    for (i = 0; i < nmap; i++) {
	static_cmap[i] = cmap[i];
    }
    for (k = 0; k < nmap - 1; k++) {
	cc = p.cis[k];
	ccp1 = p.cis[k + 1];
	setcolor(static_cmap[k]);
	if (e > cc && e <= ccp1) {
	    fillcolor(4, px, py);
	}
    }
}

void do_isol2(int gridno, int itype, double *q, struct isolparms p, int *cmap, int nmap)
{
    int i;
    for (i = 0; i < grid[gridno].nmel; i++) {
        fillel(i, q[i], p, cmap, nmap);
    }
    setcolor(1);
    setlinestyle(1);
    setlinewidth(1);
}

void do_isol(int gridno, int itype, double *q, struct isolparms p, int *cmap, int nmap)
{
    labcnt = 0;
    do_grid_isol(&grid[gridno], itype, q, p, cmap, nmap);
}

void do_grid_isol(Grid * gr, int itype, double *q, struct isolparms p, int *cmap, int nmap)
{
    int i, n1, n2, n3, n4, n5, n6, imin, imax;
    double x1, x2, y1, y2, x3, y3;
    double c1, c2, c3;
    double *x = gr->xord, *y = gr->yord;

    if (save_isoline >= 0) {
	sfp = fopen(saveisol_fname, "w");
    }
    if (p.nisol <= 0) {
	return;
    }
    for (i = 0; i < nmap; i++) {
	static_cmap[i] = cmap[i];
    }
    ip = p;
    if (ip.type == 1 || ip.type == 2) {
	for (i = 0; i < gr->nmel; i++) {
	    if (!(i % 100) && check_action()) {
		cancel_action(0);
		if (cancel_flag)
		    return;
	    }
	    if (gr->fdgrid && gr->icon[i].wetdry) {
		continue;
	    }
	    imin = 0;
	    imax = ip.nisol - 1;
	    switch (gr->icon[i].type) {
	    case 0:		/* cell centered data */
		fillel2(i, q[i], ip, cmap, nmap);
		break;
	    case 3:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		if (itype > 1) {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		} else {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		}
		break;
	    case 4:
	    case 8:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		if (itype > 1) {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		} else {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		}
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[2];
		n3 = gr->icon[i].nl[3];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		if (itype > 1) {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		} else {
		    drawisopat(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		}
		break;
	    case 6:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		n4 = gr->icon[i].nl[3];
		n5 = gr->icon[i].nl[4];
		n6 = gr->icon[i].nl[5];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		cc4 = q[n4];
		cc5 = q[n5];
		cc6 = q[n6];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonq(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, cc4, cc5, cc6, &imin, &imax);
		split(x1, y1, x2, y2, x3, y3, p.nsplits, QUADRATIC, 1, imin, imax);
		break;
	    default:
		break;
	    }
	}
    }
    if (ip.type == 0 || ip.type == 2) {
	setcolor(1);
	for (i = 0; i < gr->nmel; i++) {
	    if (!(i % 100) && check_action()) {
		cancel_action(0);
		if (cancel_flag)
		    return;
	    }
	    if (gr->fdgrid && gr->icon[i].wetdry) {
		continue;
	    }
	    switch (gr->icon[i].type) {
	    case 3:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		drawiso(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		break;
	    case 4:
	    case 8:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		drawiso(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[2];
		n3 = gr->icon[i].nl[3];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonl(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, &imin, &imax);
		drawiso(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, imin, imax);
		break;
	    case 6:
		n1 = gr->icon[i].nl[0];
		n2 = gr->icon[i].nl[1];
		n3 = gr->icon[i].nl[2];
		n4 = gr->icon[i].nl[3];
		n5 = gr->icon[i].nl[4];
		n6 = gr->icon[i].nl[5];
		cc1 = q[n1];
		cc2 = q[n2];
		cc3 = q[n3];
		cc4 = q[n4];
		cc5 = q[n5];
		cc6 = q[n6];
		aum = x[n2] * y[n3] - x[n3] * y[n2];
		bum = y[n2] - y[n3];
		cum = x[n3] - x[n2];
		ado = x[n3] * y[n1] - x[n1] * y[n3];
		bdo = y[n3] - y[n1];
		cdo = x[n1] - x[n3];
		atr = x[n1] * y[n2] - x[n2] * y[n1];
		btr = y[n1] - y[n2];
		ctr = x[n2] - x[n1];
		are2 = aum + ado + atr;
		x1 = x[n1];
		y1 = y[n1];
		x2 = x[n2];
		y2 = y[n2];
		x3 = x[n3];
		y3 = y[n3];
		if (!(symok(x1, y1) || symok(x2, y2) || symok(x3, y3))) {
		    break;
		}
		checkonq(x1, y1, x2, y2, x3, y3, cc1, cc2, cc3, cc4, cc5, cc6, &imin, &imax);
		split(x1, y1, x2, y2, x3, y3, p.nsplits, QUADRATIC, 0, imin, imax);
		break;
	    default:
		break;
	    }
	}
    }
    if (save_isoline >= 0) {
	fclose(sfp);
    }
}

int checkonl(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int *imin, int *imax)
{
    int i, k, l = 0, mincnt = 0, maxcnt = ip.nisol - 1;
    double cc1, cc2, cmin, cmax;
    cmin = c1;
    cmax = c1;
    if (c2 < cmin) {
	cmin = c2;
    }
    if (c3 < cmin) {
	cmin = c3;
    }
    if (c2 > cmax) {
	cmax = c2;
    }
    if (c3 > cmax) {
	cmax = c3;
    }
    if (cmax < ip.cis[0] || cmin > ip.cis[ip.nisol - 1]) {
	*imin = *imax = 0;
	return 0;
    }
    for (k = 1; k < ip.nisol; k++) {
	cc1 = ip.cis[k - 1];
	cc2 = ip.cis[k];
	if (cc1 <= cmin && cc2 >= cmin) {
	    mincnt = k - 1;
	}
	if (cc1 <= cmax && cc2 >= cmax) {
	    maxcnt = k;
	}
    }
    *imin = mincnt;
    *imax = maxcnt;
    return (maxcnt - (mincnt + 1));
}

int checkonq(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, double c4, double c5, double c6, int *imin, int *imax)
{
    int i, k, l = 0, mincnt = 0, maxcnt = ip.nisol - 1;
    double cc1, cc2, cmin, cmax;
    cmin = c1;
    cmax = c1;
    if (c2 < cmin) {
	cmin = c2;
    }
    if (c3 < cmin) {
	cmin = c3;
    }
    if (c2 > cmax) {
	cmax = c2;
    }
    if (c3 > cmax) {
	cmax = c3;
    }
    if (c4 < cmin) {
	cmin = c4;
    }
    if (c5 < cmin) {
	cmin = c5;
    }
    if (c6 < cmin) {
	cmin = c6;
    }
    if (c4 > cmax) {
	cmax = c4;
    }
    if (c5 > cmax) {
	cmax = c5;
    }
    if (c6 > cmax) {
	cmax = c6;
    }
    if (cmax < ip.cis[0] || cmin > ip.cis[ip.nisol - 1]) {
	*imin = *imax = 0;
	return 0;
    }
    for (k = 1; k < ip.nisol; k++) {
	cc1 = ip.cis[k - 1];
	cc2 = ip.cis[k];
	if (cc1 <= cmin && cc2 >= cmin) {
	    mincnt = k - 1;
	}
	if (cc1 <= cmax && cc2 >= cmax) {
	    maxcnt = k;
	}
    }
    *imin = mincnt;
    *imax = maxcnt;
    return (mincnt == maxcnt);
}

int inrange(double c, double c1, double c2)
{
    return ((c >= c1 && c <= c2) || (c >= c2 && c <= c1));
}

void drawisopat(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int imin, int imax)
{
    int i, j, k, l = 0, cnt = 0;
    int n1, n2, n3, nmin;
    double cc;
    double px[50], py[50];
    double x[3], y[3], c[3];
    char buf[80];
    double cc1, cc2, cmin, cmax;
    x[0] = x1;
    y[0] = y1;
    x[1] = x2;
    y[1] = y2;
    x[2] = x3;
    y[2] = y3;
    c[0] = c1;
    c[1] = c2;
    c[2] = c3;
    for (k = imin; k < imax; k++) {
	cc1 = ip.cis[k];
	cc2 = ip.cis[k + 1];
	l = 0;
	for (j = 0; j < 3; j++) {
	    if (inrange(c[j], cc1, cc2)) {
		px[l] = x[j];
		py[l++] = y[j];
	    }
	    if (c[j] < c[(j + 1) % 3]) {
		if (inrange(cc1, c[j], c[(j + 1) % 3])) {
		    if (c[(j + 1) % 3] != c[j]) {
			px[l] = (cc1 - c[j]) * (x[(j + 1) % 3] - x[j]) / (c[(j + 1) % 3] - c[j]) + x[j];
			py[l++] = (cc1 - c[j]) * (y[(j + 1) % 3] - y[j]) / (c[(j + 1) % 3] - c[j]) + y[j];
		    }
		}
		if (inrange(cc2, c[j], c[(j + 1) % 3])) {
		    if (c[(j + 1) % 3] != c[j]) {
			px[l] = (cc2 - c[j]) * (x[(j + 1) % 3] - x[j]) / (c[(j + 1) % 3] - c[j]) + x[j];
			py[l++] = (cc2 - c[j]) * (y[(j + 1) % 3] - y[j]) / (c[(j + 1) % 3] - c[j]) + y[j];
		    }
		}
	    } else if (c[j] > c[(j + 1) % 3]) {
		if (inrange(cc2, c[j], c[(j + 1) % 3])) {
		    if (c[(j + 1) % 3] != c[j]) {
			px[l] = (cc2 - c[j]) * (x[(j + 1) % 3] - x[j]) / (c[(j + 1) % 3] - c[j]) + x[j];
			py[l++] = (cc2 - c[j]) * (y[(j + 1) % 3] - y[j]) / (c[(j + 1) % 3] - c[j]) + y[j];
		    }
		}
		if (inrange(cc1, c[j], c[(j + 1) % 3])) {
		    if (c[(j + 1) % 3] != c[j]) {
			px[l] = (cc1 - c[j]) * (x[(j + 1) % 3] - x[j]) / (c[(j + 1) % 3] - c[j]) + x[j];
			py[l++] = (cc1 - c[j]) * (y[(j + 1) % 3] - y[j]) / (c[(j + 1) % 3] - c[j]) + y[j];
		    }
		}
	    }
	}
	if (l >= 2 && k < ip.nisol - 1) {
	    setcolor(static_cmap[k]);
	    fillcolor(l, px, py);
	}
    }
}

/*
   compute the isolines
 */
void drawiso(double x1, double y1, double x2, double y2, double x3, double y3, double c1, double c2, double c3, int imin, int imax)
{
    int i, k, l = 0;
    int n1, n2, n3;
    double cc;
    double px[3], py[3];
    char buf[80];

    for (k = 0; k < ip.nisol; k++) {
	cc = ip.cis[k];
	l = 0;
	/* crosses at n1 n2 */
	if ((c1 >= cc && c2 <= cc) || (c1 <= cc && c2 >= cc)) {
	    if (c1 != c2) {
		px[l] = (cc - c1) * (x2 - x1) / (c2 - c1) + x1;
		py[l++] = (cc - c1) * (y2 - y1) / (c2 - c1) + y1;
	    }
	}
	/* crosses at n1 n3 */
	if ((c1 >= cc && c3 <= cc) || (c1 <= cc && c3 >= cc)) {
	    if (c1 != c3) {
		px[l] = (cc - c1) * (x3 - x1) / (c3 - c1) + x1;
		py[l++] = (cc - c1) * (y3 - y1) / (c3 - c1) + y1;
	    }
	}
	/* crosses at n2 n3 */
	if ((c2 >= cc && c3 <= cc) || (c2 <= cc && c3 >= cc)) {
	    if (c2 != c3) {
		px[l] = (cc - c2) * (x3 - x2) / (c3 - c2) + x2;
		py[l++] = (cc - c2) * (y3 - y2) / (c3 - c2) + y2;
	    }
	}
	if (l == 2) {
	    my_move2(px[0], py[0]);
	    my_draw2(px[1], py[1]);
	    if (save_isoline == k) {
		fprintf(sfp, "%lf %lf %lf\n", px[0], py[0], cc);
		fprintf(sfp, "%lf %lf %lf\n", px[1], py[1], cc);
	    }
	    labcnt++;
	    if (ip.marker && ip.markstep) {
		if ((labcnt % ip.markstep) == 0) {
		    sprintf(buf, "%.*lf", ip.p.prec, cc);
		    writestr(px[1], py[1], 0, 0, buf);
		    labcnt = 0;
		}
	    }
	}
	if (l == 3) {
	    my_move2(px[0], py[0]);
	    my_draw2(px[1], py[1]);
	    my_draw2(px[2], py[2]);
	    if (save_isoline == k) {
		fprintf(sfp, "%lf %lf %lf 3\n", px[0], py[0], cc);
		fprintf(sfp, "%lf %lf %lf 3\n", px[1], py[1], cc);
		fprintf(sfp, "%lf %lf %lf 3\n", px[1], py[1], cc);
		fprintf(sfp, "%lf %lf %lf 3\n", px[2], py[2], cc);
	    }
	    labcnt++;
	    if (ip.marker && ip.markstep) {
		if ((labcnt % ip.markstep) == 0) {
		    sprintf(buf, "%.*lf", ip.p.prec, cc);
		    writestr(px[1], py[1], 0, 0, buf);
		    labcnt = 0;
		}
	    }
	}
    }
}

void split(double x1, double y1, double x2, double y2, double x3, double y3, int ilev, int itype, int ipen, int imin, int imax)
{
    double x4, y4, x5, y5, x6, y6;

    if (ilev == 0) {
	double c1, c2, c3;

	switch (itype) {
	case QUADRATIC:
	    concna(x1, y1, x2, y2, x3, y3, &c1, &c2, &c3);
	    break;
	case LINEAR:
	    concna2(x1, y1, x2, y2, x3, y3, &c1, &c2, &c3);
	    break;
	case TEANL_AMP:
	    comp_amp(x1, y1, x2, y2, x3, y3, &c1, &c2, &c3);
	    break;
	case TEANL_PHA:
	    comp_pha(x1, y1, x2, y2, x3, y3, &c1, &c2, &c3);
	    break;
	}

	if (ipen) {
	    checkonl(x1, y1, x2, y2, x3, y3, c1, c2, c3, &imin, &imax);
	    drawisopat(x1, y1, x2, y2, x3, y3, c1, c2, c3, imin, imax);
	} else {
	    checkonl(x1, y1, x2, y2, x3, y3, c1, c2, c3, &imin, &imax);
	    drawiso(x1, y1, x2, y2, x3, y3, c1, c2, c3, imin, imax);
	}
    } else {
	x4 = 0.5 * (x1 + x2);
	y4 = 0.5 * (y1 + y2);
	x6 = 0.5 * (x1 + x3);
	y6 = 0.5 * (y1 + y3);
	x5 = 0.5 * (x2 + x3);
	y5 = 0.5 * (y2 + y3);
	split(x1, y1, x4, y4, x6, y6, ilev - 1, itype, ipen, imin, imax);
	split(x3, y3, x6, y6, x5, y5, ilev - 1, itype, ipen, imin, imax);
	split(x4, y4, x5, y5, x6, y6, ilev - 1, itype, ipen, imin, imax);
	split(x2, y2, x5, y5, x4, y4, ilev - 1, itype, ipen, imin, imax);
    }
}

void concna(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3)
{
    double h1, h2, h3, h4, h5, h6, arei;

    arei = 1. / are2;
    h3 = (atr + btr * x1 + ctr * y1) * arei;
    h2 = (ado + bdo * x1 + cdo * y1) * arei;
    h1 = 1. - h2 - h3;
    h4 = 4. * h2 * h1;
    h5 = 4. * h2 * h3;
    h6 = 4. * h3 * h1;
    h1 = h1 - 0.5 * (h4 + h6);
    h2 = h2 - 0.5 * (h4 + h5);
    h3 = h3 - 0.5 * (h5 + h6);
    *c1 = h1 * cc1 + h2 * cc2 + h3 * cc3 + h4 * cc4 + h5 * cc5 + h6 * cc6;
    h3 = (atr + btr * x2 + ctr * y2) * arei;
    h2 = (ado + bdo * x2 + cdo * y2) * arei;
    h1 = 1. - h2 - h3;
    h4 = 4. * h2 * h1;
    h5 = 4. * h2 * h3;
    h6 = 4. * h3 * h1;
    h1 = h1 - 0.5 * (h4 + h6);
    h2 = h2 - 0.5 * (h4 + h5);
    h3 = h3 - 0.5 * (h5 + h6);
    *c2 = h1 * cc1 + h2 * cc2 + h3 * cc3 + h4 * cc4 + h5 * cc5 + h6 * cc6;
    h3 = (atr + btr * x3 + ctr * y3) * arei;
    h2 = (ado + bdo * x3 + cdo * y3) * arei;
    h1 = 1. - h2 - h3;
    h4 = 4. * h2 * h1;
    h5 = 4. * h2 * h3;
    h6 = 4. * h3 * h1;
    h1 = h1 - 0.5 * (h4 + h6);
    h2 = h2 - 0.5 * (h4 + h5);
    h3 = h3 - 0.5 * (h5 + h6);
    *c3 = h1 * cc1 + h2 * cc2 + h3 * cc3 + h4 * cc4 + h5 * cc5 + h6 * cc6;
}

void concna2(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3)
{
    double h1, h2, h3, arei;

    arei = 1. / are2;
    h3 = (atr + btr * x1 + ctr * y1) * arei;
    h2 = (ado + bdo * x1 + cdo * y1) * arei;
    h1 = 1. - h2 - h3;
    *c1 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x2 + ctr * y2) * arei;
    h2 = (ado + bdo * x2 + cdo * y2) * arei;
    h1 = 1. - h2 - h3;
    *c2 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x3 + ctr * y3) * arei;
    h2 = (ado + bdo * x3 + cdo * y3) * arei;
    h1 = 1. - h2 - h3;
    *c3 = h1 * cc1 + h2 * cc2 + h3 * cc3;
}

void comp_amp(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3)
{
    double h1, h2, h3, arei;

    arei = 1. / are2;
    h3 = (atr + btr * x1 + ctr * y1) * arei;
    h2 = (ado + bdo * x1 + cdo * y1) * arei;
    h1 = 1. - h2 - h3;
    *c1 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x2 + ctr * y2) * arei;
    h2 = (ado + bdo * x2 + cdo * y2) * arei;
    h1 = 1. - h2 - h3;
    *c2 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x3 + ctr * y3) * arei;
    h2 = (ado + bdo * x3 + cdo * y3) * arei;
    h1 = 1. - h2 - h3;
    *c3 = h1 * cc1 + h2 * cc2 + h3 * cc3;
}

void comp_pha(double x1, double y1, double x2, double y2, double x3, double y3, double *c1, double *c2, double *c3)
{
    double h1, h2, h3, arei;

    arei = 1. / are2;
    h3 = (atr + btr * x1 + ctr * y1) * arei;
    h2 = (ado + bdo * x1 + cdo * y1) * arei;
    h1 = 1. - h2 - h3;
    *c1 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x2 + ctr * y2) * arei;
    h2 = (ado + bdo * x2 + cdo * y2) * arei;
    h1 = 1. - h2 - h3;
    *c2 = h1 * cc1 + h2 * cc2 + h3 * cc3;
    h3 = (atr + btr * x3 + ctr * y3) * arei;
    h2 = (ado + bdo * x3 + cdo * y3) * arei;
    h1 = 1. - h2 - h3;
    *c3 = h1 * cc1 + h2 * cc2 + h3 * cc3;
}
