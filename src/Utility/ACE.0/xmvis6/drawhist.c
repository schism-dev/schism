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
 * Draw history markers
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawhist.c,v 1.11 2007/01/13 16:48:42 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

double eval_teanl(int flowno, int el, double t, double xp, double yp, int redo);
double eval_ela(int concno, int el, int step, double t, double xp, double yp, int redo);
double eval_adcirc(int flowno, int el, int step, double t, double xp, double yp, int redo);
double eval_adc3d(int flowno, int el, int step, double t, double xp, double yp, int redo);
double eval_model(int obj, int w1, int w2, int step, int ind, double t, double x, double y, int redo);
void load_hist(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);
void write_hist(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);

int read_hist(int histno, int type, char *fname)
{
    char buf[255];
    int i, npts;
    double x, y;
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 0;
    }
    //fgets(buf, 255, fp);
    //fgets(buf, 255, fp);
    //sscanf(buf, "%d %lf %lf", &npts, &x, &y);
    if (type == FREQ) {
    } else if (type == TIME) {
	npts = 0;
	while (fgets(buf, 255, fp) != NULL) {
	    npts++;
	}
	if (npts < 2) {
	    fclose(fp);
	    errwin("No points read from time history file, must have at least 2 points in the file.");
	    return 0;
	}
	rewind(fp);
	hist[histno].active = 1;
	//hist[histno].x = x;
	//hist[histno].y = y;
	hist[histno].nsteps = npts;
	if (hist[histno].t) {
	    free(hist[histno].t);
	}
	if (hist[histno].e) {
	    free(hist[histno].e);
	}
	hist[histno].t = (double *) malloc(npts * sizeof(double));
	hist[histno].e = (double *) malloc(npts * sizeof(double));
	for (i = 0; i < npts; i++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%lf %lf", &hist[histno].t[i], &hist[histno].e[i]);
	}
        dminmax(npts, hist[histno].t, &hist[histno].tmin, &hist[histno].tmax);
        dminmax(npts, hist[histno].e, &hist[histno].emin, &hist[histno].emax);
    }
    fclose(fp);
    return 1;
}

void drawhist(int gno, int n, int step, double t)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    char buf[10];

    if (g[gno].hbox[n].display == ON) {
	setcolor(g[gno].hbox[n].p.color);
	my_move2(g[gno].hbox[n].x, g[gno].hbox[n].y);
	if (g[gno].hbox[n].loctype == VIEW) {
	    double x, y;
	    view2world(g[gno].hbox[n].locx, g[gno].hbox[n].locy, &x, &y);
	    my_draw2(x, y);
	} else {
	    my_draw2(g[gno].hbox[n].locx, g[gno].hbox[n].locy);
	}
    }
    if (g[gno].hbox[n].display_marker == ON) {
/* redefine the world and viewport */
	switch (g[gno].hbox[n].attach) {
	case 0:
	    if (g[gno].hbox[n].loctype == WORLD) {
	        world2view(g[gno].hbox[n].locx, g[gno].hbox[n].locy, &vx1, &vy1);
	    } else {
		vx1 = g[gno].hbox[n].locx;
		vy1 = g[gno].hbox[n].locy;
	    }
	    vx2 = vx1 + g[gno].hbox[n].vx;
	    vy2 = vy1 + g[gno].hbox[n].vy;
	    break;
	case 1:
	    if (g[gno].hbox[n].loctype == WORLD) {
	        world2view(g[gno].hbox[n].locx, g[gno].hbox[n].locy, &vx2, &vy1);
	    } else {
		vx2 = g[gno].hbox[n].locx;
		vy1 = g[gno].hbox[n].locy;
	    }
	    vx1 = vx2 - g[gno].hbox[n].vx;
	    vy2 = vy1 + g[gno].hbox[n].vy;
	    break;
	case 2:
	    if (g[gno].hbox[n].loctype == WORLD) {
	        world2view(g[gno].hbox[n].locx, g[gno].hbox[n].locy, &vx1, &vy2);
	    } else {
		vx1 = g[gno].hbox[n].locx;
		vy2 = g[gno].hbox[n].locy;
	    }
	    vx2 = vx1 + g[gno].hbox[n].vx;
	    vy1 = vy2 - g[gno].hbox[n].vy;
	    break;
	case 3:
	    if (g[gno].hbox[n].loctype == WORLD) {
	        world2view(g[gno].hbox[n].locx, g[gno].hbox[n].locy, &vx2, &vy2);
	    } else {
		vx2 = g[gno].hbox[n].locx;
		vy2 = g[gno].hbox[n].locy;
	    }
	    vx1 = vx2 - g[gno].hbox[n].vx;
	    vy1 = vy2 - g[gno].hbox[n].vy;
	    break;
	}
	setcolor(g[gno].hbox[n].p.color);
	defineworld(g[gno].hbox[n].wx1, g[gno].hbox[n].wy1, g[gno].hbox[n].wx2, g[gno].hbox[n].wy2, 0, 0);
	viewport(vx1, vy1, vx2, vy2);

	setcolor(g[gno].hbox[n].p.fillcol);
	fillrectcolor(g[gno].hbox[n].wx1, g[gno].hbox[n].wy1, g[gno].hbox[n].wx2, g[gno].hbox[n].wy2);

	setcolor(g[gno].hbox[n].p.color);
	rect(g[gno].hbox[n].wx1, g[gno].hbox[n].wy1, g[gno].hbox[n].wx2, g[gno].hbox[n].wy2);
	//rectstr(g[gno].hbox[n].wx1, g[gno].hbox[n].wy1, g[gno].hbox[n].wx2, g[gno].hbox[n].wy2, g[gno].hbox[n].precx, g[gno].hbox[n].precy);

	setlinewidth(g[gno].hbox[n].p.linew);
	j = 0;
	if (display_hist(gno, n, HISTORY, j, 0)) {
	    setcolor(1);
	    load_hist(gno, n, flowf[j].grid, HISTORY, 0, 0, step, t);
	}
	for (j = 0; j < MAXTEANL; j++) {
	    if (display_hist(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
		setcolor(1);
		load_hist(gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	    }
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    if (display_hist(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
		setcolor(1);
		load_hist(gno, n, -1, ADCIRC, ELEV, j, step, t);
	    }
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    if (display_hist(gno, n, ADCIRC, MAG, j) && object_isactive(ADCIRC, FLOW, j)) {
		setcolor(1);
		load_hist(gno, n, -1, ADCIRC, MAG, j, step, t);
	    }
	}
	defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, islogx(gno), islogy(gno));
	viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
    }
    setlinewidth(1);
    setcolor(1);
}

int display_hist(int gno, int sno, int obj, int feature, int w2)
{
    switch (obj) {
    case HISTORY:
	return (g[gno].hbox[sno].thist);
    case TEANL:
	if (feature == ELEV) {
	    return (g[gno].hbox[sno].teanl[w2]);
	} else if (feature == MAG) {
	    return (g[gno].hbox[sno].teanl[w2]);
	} else {
	    fprintf(stderr, "No such feature in TEANL histories\n");
	    return 0;
	}
    case ADCIRC:
	if (feature == ELEV)
	    return (g[gno].hbox[sno].adcirc[w2]);
	if (feature == MAG)
	    return (g[gno].hbox[sno].adcircflow[w2]);
	break;
    default:
	return 0;
    }
    return 0;
}

void load_hist(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double d, x, y;
    double dsave;
    int i, ind;
    if (obj == HISTORY) {
	my_move2(hist[n].t[0], hist[n].e[0]);
	for (i = 1; i < hist[n].nsteps; i++) {
	    my_draw2(hist[n].t[i], hist[n].e[i]);
	}
	//my_move2(timeclock.t[step], -5.0);
	//my_draw2(timeclock.t[step], 5.0);
	setcolor(2);
	//my_move2(t, hist[n].emin);
	//my_draw2(t, hist[n].emax);
	my_move2(t, g[gno].hbox[n].wy1);
	my_draw2(t, g[gno].hbox[n].wy2);
        return;
    }
    if (g[gno].hbox[n].elem == -1) {
	if (obj == ADCIRC) {
	    FindElement(&flowt[w2].g, x = g[gno].hbox[n].x, y = g[gno].hbox[n].y, &ind);
	} else {
	    find_element(gridno, x = g[gno].hbox[n].x, y = g[gno].hbox[n].y, &ind);
	}
    }
    if (ind >= 0) {
	d = eval_model(obj, w1, w2, 0, ind, timeclock.start, x, y, 0);
	my_move2(timeclock.start, d);
	for (i = 1; i <= timeclock.nsteps; i++) {
	    if (obj == ADCIRC && i >= flowt[w2].nsteps)
		break;
	    d = eval_model(obj, w1, w2, i, ind, timeclock.t[i], x, y, 0);
	    /* printf("%d %d %lf %lf\n", obj == ADCIRC, ind, timeclock.t[i], d); */
	    my_draw2(timeclock.t[i], d);
	    if (i == step) {
		dsave = d;
	    }
	}
	my_move2(timeclock.t[step], -5.0);
	my_draw2(timeclock.t[step], dsave);
    } else {
	fprintf(stderr, "Element not found\n");
    }
}

int viewhist(FILE * fp, int gno, int n, int step, double t, int cnt)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    int writecount;
    char buf[10];
    writecount = 0;
    j = 0;
    if (display_hist(gno, n, HISTORY, j, 0)) {
	fprintf(fp, "@legend string %d \"Time history\"\n", cnt++, j + 1);
	writecount++;
	write_hist(fp, gno, n, flowf[j].grid, HISTORY, 0, 0, step, t);
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_hist(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
	    fprintf(fp, "@legend string %d \"TEANL set %d\"\n", cnt++, j + 1);
	    writecount++;
	    write_hist(fp, gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_hist(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
	    fprintf(fp, "@legend string %d \"ELCIRC set %d\"\n", cnt++, j + 1);
	    writecount++;
	    write_hist(fp, gno, n, -1, ADCIRC, ELEV, j, step, t);
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_hist(gno, n, ADCIRC, MAG, j) && object_isactive(ADCIRC, FLOW, j)) {
	    fprintf(fp, "@legend string %d \"ELCIRC set %d\"\n", cnt++, j + 1);
	    writecount++;
	    write_hist(fp, gno, n, -1, ADCIRC, MAG, j, step, t);
	}
    }
    return writecount;
}

void write_hist(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double d, x, y;
    int i, ind = -1;
    if (obj == HISTORY) {
    }
    if (g[gno].hbox[n].elem == -1) {
	if (obj == ADCIRC) {
	    FindElement(&flowt[w2].g, x = g[gno].hbox[n].x, y = g[gno].hbox[n].y, &ind);
	    /*
	       printf("Element = %d x,y %lf %lf stuff %d %d %d %d %d %d %d %lf\n", ind, x, y, gno, n, gridno, obj, w1, w2, step, t);
	     */
	} else {
	    find_element(gridno, x = g[gno].hbox[n].x, y = g[gno].hbox[n].y, &ind);
	}
    }
    if (ind >= 0) {
	d = eval_model(obj, w1, w2, 0, ind, timeclock.start, x, y, 0);
	fprintf(fp, "%lf %lf\n", timeclock.start, d);
	for (i = 1; i <= step; i++) {
	    if (obj == ADCIRC && i >= flowt[w2].nsteps)
		break;
	    d = eval_model(obj, w1, w2, i, ind, timeclock.t[i], x, y, 1);
	    fprintf(fp, "%lf %lf\n", timeclock.t[i], d);
	}
    } else {
	fprintf(stderr, "Element not found\n");
    }
    fprintf(fp, "&\n");
}
