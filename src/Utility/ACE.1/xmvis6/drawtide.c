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
 * Draw tide station markers
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawtide.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"
#include "externs.h"

double eval_teanl(int flowno, int el, double t, double xp, double yp, int redo);
double eval_adcirc(int flowno, int el, int step, double t, double xp, double yp, int redo);
double eval_model(int obj, int w1, int w2, int step, int ind, double t, double x, double y, int redo);
void load_tide(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);
void write_tide(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);

void drawtidestation(int gno, int n, int step, double t)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2, x, y;
    double wx1, wx2, wy1, wy2;
    char buf[256];
    int fudge = stringextenty(1.0 * devcharsize, "Ny");	/* TODO replace 1.0 with
							 * charsize from prop
							 * struct */

    if (g[gno].tidestat[n].display == ON) {
	setcolor(g[gno].tidestat[n].p.color);
	my_circle(g[gno].tidestat[n].x, g[gno].tidestat[n].y);
	if (g[gno].tidestat[n].display_marker == OFF) {
	    sprintf(buf, "  %s", tidestat[n]->name);
	    writestr(g[gno].tidestat[n].x, g[gno].tidestat[n].y, 0, 0, buf);
	    switch (g[gno].tidestat[n].display_ampphase) {
	    case 0:
		break;
	    case 1:
		for (i = 0; i < tidestat[n]->nfreq; i++) {
		    sprintf(buf, "  %s %.4lf %.2lf", tidestat[n]->freqname[i], tidestat[n]->elamp[i], tidestat[n]->elphase[i] * 180.0 / M_PI);
		    x = g[gno].tidestat[n].x;
		    y = g[gno].tidestat[n].y;
		    (*devwritestr) ((*devconvx) (x), (*devconvy) (y) - (i + 1) * fudge, 0, buf, 0, 1);
		}
		break;
	    case 2:
		i = g[gno].tidestat[n].use_freq;
		sprintf(buf, "  %s %s %.2lf %.2lf", tidestat[n]->name[i], tidestat[n]->freqname[i], tidestat[n]->elamp[i], tidestat[n]->elphase[i]);
		writestr(g[gno].tidestat[n].x, g[gno].tidestat[n].y, 0, 0, buf);
		break;
	    }			/* switch */
	}			/* if */
    }				/* if */
    if (g[gno].tidestat[n].display_marker == ON) {
/* connecting line */
	my_move2(g[gno].tidestat[n].x, g[gno].tidestat[n].y);
	my_draw2(g[gno].tidestat[n].locx, g[gno].tidestat[n].locy);
/* redefine the world and viewport */
	switch (g[gno].tidestat[n].attach) {
	case 0:
	    world2view(g[gno].tidestat[n].locx, g[gno].tidestat[n].locy, &vx1, &vy1);
	    vx2 = vx1 + g[gno].tidestat[n].vx;
	    vy2 = vy1 + g[gno].tidestat[n].vy;
	    break;
	case 1:
	    world2view(g[gno].tidestat[n].locx, g[gno].tidestat[n].locy, &vx2, &vy1);
	    vx1 = vx2 - g[gno].tidestat[n].vx;
	    vy2 = vy1 + g[gno].tidestat[n].vy;
	    break;
	case 2:
	    world2view(g[gno].tidestat[n].locx, g[gno].tidestat[n].locy, &vx1, &vy2);
	    vx2 = vx1 + g[gno].tidestat[n].vx;
	    vy1 = vy2 - g[gno].tidestat[n].vy;
	    break;
	case 3:
	    world2view(g[gno].tidestat[n].locx, g[gno].tidestat[n].locy, &vx2, &vy2);
	    vx1 = vx2 - g[gno].tidestat[n].vx;
	    vy1 = vy2 - g[gno].tidestat[n].vy;
	    break;
	}
	setcolor(g[gno].tidestat[n].p.color);
	defineworld(g[gno].tidestat[n].wx1, g[gno].tidestat[n].wy1, g[gno].tidestat[n].wx2, g[gno].tidestat[n].wy2, 0, 0);
	viewport(vx1, vy1, vx2, vy2);

	setcolor(g[gno].tidestat[n].p.fillcol);
	fillrectcolor(g[gno].tidestat[n].wx1, g[gno].tidestat[n].wy1, g[gno].tidestat[n].wx2, g[gno].tidestat[n].wy2);

	setcolor(g[gno].tidestat[n].p.color);
	rect(g[gno].tidestat[n].wx1, g[gno].tidestat[n].wy1, g[gno].tidestat[n].wx2, g[gno].tidestat[n].wy2);
	rectstr(g[gno].tidestat[n].wx1, g[gno].tidestat[n].wy1, g[gno].tidestat[n].wx2, g[gno].tidestat[n].wy2, g[gno].tidestat[n].precx, g[gno].tidestat[n].precy);

	setlinewidth(g[gno].tidestat[n].p.linew);
	if (display_tide(gno, n, TIDESTATION, 0, 0)) {
	    setcolor(g[gno].tidestat[n].p.color);
	    load_tide(gno, n, 0, TIDESTATION, 0, 0, step, t);
	}
	for (j = 0; j < MAXTEANL; j++) {
	    if (display_tide(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
		setcolor(g[gno].tidestat[n].tp[j].color);
		load_tide(gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	    }
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    if (display_tide(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
		setcolor(1);
		setcolor(g[gno].tidestat[n].ap[j].color);
		load_tide(gno, n, -1, ADCIRC, ELEV, j, step, t);
	    }
	}
	defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, islogx(gno), islogy(gno));
	viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
    }
    setlinewidth(1);
    setcolor(1);
}

int display_tide(int gno, int sno, int obj, int feature, int w2)
{
    switch (obj) {
    case TIDESTATION:
	return (g[gno].tidestat[sno].display);
	break;
    case TEANL:
	if (feature == ELEV) {
	    return (g[gno].tidestat[sno].teanl[w2]);
	}
	break;
    case ADCIRC:
	if (feature == ELEV) {
	    return (g[gno].tidestat[sno].adcirc[w2]);
	}
	break;
    default:
	return 0;
    }
    return 0;
}

double EvalTideStationElev(TideStation * tides, double t, int constituent);

void load_tide(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double d, x, y;
    int i, ind;
    if (obj == TIDESTATION) {
	d = EvalTideStationElev(tidestat[n], timeclock.start, -1);
	my_move2(timeclock.start, d);
	for (i = 1; i <= step; i++) {
	    d = EvalTideStationElev(tidestat[n], timeclock.t[i], -1);
	    my_draw2(timeclock.t[i], d);
	}
    } else {
	find_element(gridno, x = g[gno].tidestat[n].x, y = g[gno].tidestat[n].y, &ind);
	if (ind >= 0) {
	    d = eval_model(obj, w1, w2, 0, ind, timeclock.start, x, y, 0);
	    my_move2(timeclock.start, d);
	    for (i = 1; i <= step; i++) {
		d = eval_model(obj, w1, w2, i, ind, timeclock.t[i], x, y, 1);
		my_draw2(timeclock.t[i], d);
	    }
	} else {
	    fprintf(stderr, "Element not found, station [%s]\n", tidestat[n]->name);
	}
    }
}

int viewtide(FILE * fp, int gno, int n, int step, double t, int cnt)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    int writecount;
    char buf[10];
    writecount = 0;
    if (display_tide(gno, n, TIDESTATION, 0, 0)) {
	fprintf(fp, "@legend string %d \"Tide station\"\n", cnt++);
	writecount++;
	write_tide(fp, gno, n, 0, TIDESTATION, 0, 0, step, t);
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_tide(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
	    fprintf(fp, "@legend string %d \"TEANL set %d\"\n", cnt++, j + 1);
	    writecount++;
	    write_tide(fp, gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_tide(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
	    fprintf(fp, "@legend string %d \"ADCIRC/TABS set %d\"\n", cnt++, j + 1);
	    writecount++;
	    write_tide(fp, gno, n, -1, ADCIRC, ELEV, j, step, t);
	}
    }
    return writecount;
}

void write_tide(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double d, x, y;
    int i, ind = -1, didone = 0;
    if (obj == TIDESTATION) {
	if (step >= 0)
	    didone = 1;
	for (i = 0; i <= step; i++) {
	    d = EvalTideStationElev(tidestat[n], timeclock.t[i], -1);
	    fprintf(fp, "%lf %lf\n", timeclock.t[i], d);
	}
    } else {
	find_element(gridno, x = g[gno].tidestat[n].x, y = g[gno].tidestat[n].y, &ind);
	if (ind >= 0) {
	    if (step >= 0)
		didone = 1;
	    for (i = 0; i <= step; i++) {
		d = eval_model(obj, w1, w2, i, ind, timeclock.t[i], x, y, i == 0 ? 0 : 1);
		fprintf(fp, "%lf %lf\n", timeclock.t[i], d);
	    }
	} else {
	    fprintf(stderr, "Element not found, station [%s]\n", tidestat[n]->name);
	}
    }
    if (didone) {
	fprintf(fp, "&\n");
    }
}

double rms_tide(int gno,	/* graph number */
		int n,		/* tide station number */
		int gridno,	/* reference grid */
		int obj,	/* which model */
		int w1, int w2, int step)
{
    double d, dm, x, y, sums = 0.0;
    int i, ind = -1, didone = 0;
    find_element(gridno, x = g[gno].tidestat[n].x, y = g[gno].tidestat[n].y, &ind);
    if (ind >= 0) {
	for (i = 0; i <= step; i++) {
	    d = EvalTideStationElev(tidestat[n], timeclock.t[i], -1);
	    dm = eval_model(obj, w1, w2, i, ind, timeclock.t[i], x, y, i == 0 ? 0 : 1);
	    sums += (dm - d) * (dm - d);
	}
	sums = sqrt(sums / (step + 1));
	return sums;
    } else {
	fprintf(stderr, "Element not found, station [%s]\n", tidestat[n]->name);
    }
    return -1.0;
}

void viewrms(FILE * fp, int gno, int n, int step, double t)
{
    int i, j, cnt;
    double e, rms[200], vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    double sum;
    char buf[10];

    cnt = 0;
    for (j = 0; j < MAXTEANL; j++) {
	if (object_isactive(TEANL, j)) {
	    sum = 0.0;
	    for (i = 0; i < n; i++) {
		rms[i] = rms_tide(gno, i, flowf[j].grid, TEANL, ELEV, j, step);
		fprintf(fp, "%d %lf %s\n", i + 1, rms[i], tidestat[i]->name);
		sum += rms[i];
	    }
	    fprintf(fp, "&\n1 %lf\n&\n@s%1d symbol 15\n", sum / n, cnt * 2 + 1);
	    cnt++;
	}
    }

    for (j = 0; j < MAXADCIRC; j++) {
	if (object_isactive(ADCIRC, ELEV, j)) {
	    sum = 0.0;
	    for (i = 0; i < n; i++) {
		rms[i] = rms_tide(gno, i, -1, ADCIRC, ELEV, j, step);
		fprintf(fp, "%d %lf %s\n", i + 1, rms[i], tidestat[i]->name);
		sum += rms[i];
	    }
	    fprintf(fp, "&\n1 %lf\n&\n@s%1d symbol 15\n", sum / n, cnt * 2 + 1);
	    cnt++;
	}
    }
}
