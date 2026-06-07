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
 * Draw slice markers
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawslice.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "defines.h"
#include "globals.h"

double get_depth_element(int gno, int ind, double x, double y);
double eval_teanl(int flowno, int el, double t, double xp, double yp, int redo);
double eval_ela(int concno, int el, int step, double t, double xp, double yp, int redo);
double eval_adcirc(int flowno, int el, int step, double t, double xp, double yp, int redo);
double eval_adc3d(int flowno, int el, int step, double t, double xp, double yp, int redo);
void load_slice(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);
void write_slice(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t);
void get_slice(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t, double *x1, double *y1);

double eval_model(int obj, int w1, int w2, int step, int ind, double t, double x, double y, int redo);

void drawslice(int gno, int n, int step, double t)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    char buf[10];

    if (g[gno].sbox[n].display_slice == ON) {
	setcolor(g[gno].sbox[n].p.color);
	if (g[gno].sbox[n].type == LINE) {
	    my_move2(g[gno].sbox[n].x1, g[gno].sbox[n].y1);
	    my_draw2(g[gno].sbox[n].x2, g[gno].sbox[n].y2);
	} else {
	    my_move2(g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0]);
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		my_draw2(g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i]);
	    }
	}
    }
    if (g[gno].sbox[n].display_marker == ON) {
/* redefine the world and viewport */
	switch (g[gno].sbox[n].attach) {
	case 0:
	    world2view(g[gno].sbox[n].locx, g[gno].sbox[n].locy, &vx1, &vy1);
	    vx2 = vx1 + 0.2;
	    vy2 = vy1 + 0.1;
	    view2world(vx1, vy1, &wx1, &wy1);
	    view2world(vx2, vy1, &wx2, &wy2);
	    break;
	case 1:
	    world2view(g[gno].sbox[n].locx, g[gno].sbox[n].locy, &vx2, &vy1);
	    vx1 = vx2 - 0.2;
	    vy2 = vy1 + 0.1;
	    view2world(vx1, vy2, &wx1, &wy1);
	    view2world(vx2, vy2, &wx2, &wy2);
	    break;
	case 2:
	    world2view(g[gno].sbox[n].locx, g[gno].sbox[n].locy, &vx1, &vy2);
	    vx2 = vx1 + 0.2;
	    vy1 = vy2 - 0.1;
	    view2world(vx1, vy1, &wx1, &wy1);
	    view2world(vx2, vy1, &wx2, &wy2);
	    break;
	case 3:
	    world2view(g[gno].sbox[n].locx, g[gno].sbox[n].locy, &vx2, &vy2);
	    vx1 = vx2 - 0.2;
	    vy1 = vy2 - 0.1;
	    view2world(vx1, vy2, &wx1, &wy1);
	    view2world(vx2, vy2, &wx2, &wy2);
	    break;
	}
	setcolor(g[gno].sbox[n].p.color);
	if (g[gno].sbox[n].type == LINE) {
	    my_move2(g[gno].sbox[n].x1, g[gno].sbox[n].y1);
	    my_draw2(wx1, wy1);
	    my_move2(g[gno].sbox[n].x2, g[gno].sbox[n].y2);
	    my_draw2(wx2, wy2);
	} else if (g[gno].sbox[n].npts > 1) {
	    my_move2(g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0]);
	    my_draw2(wx1, wy1);
	    my_move2(g[gno].sbox[n].sx[g[gno].sbox[n].npts - 1], g[gno].sbox[n].sy[g[gno].sbox[n].npts - 1]);
	    my_draw2(wx2, wy2);
	}
	defineworld(g[gno].sbox[n].wx1, g[gno].sbox[n].wy1, g[gno].sbox[n].wx2, g[gno].sbox[n].wy2, 0, 0);
	viewport(vx1, vy1, vx2, vy2);

	setcolor(g[gno].sbox[n].p.fillcol);
	fillrectcolor(g[gno].sbox[n].wx1, g[gno].sbox[n].wy1, g[gno].sbox[n].wx2, g[gno].sbox[n].wy2);

	setcolor(g[gno].sbox[n].p.color);
	rect(g[gno].sbox[n].wx1, g[gno].sbox[n].wy1, g[gno].sbox[n].wx2, g[gno].sbox[n].wy2);
	rectstr(g[gno].sbox[n].wx1, g[gno].sbox[n].wy1, g[gno].sbox[n].wx2, g[gno].sbox[n].wy2, g[gno].sbox[n].precx, g[gno].sbox[n].precy);

	setlinewidth(g[gno].sbox[n].p.linew);
	for (j = 0; j < MAXGRIDS; j++) {
	    if (display_slice(gno, n, BATH, j, 0) && object_isactive(GRID, j)) {
		setcolor(1);
		load_slice(gno, n, j, BATH, 0, j, step, t);
	    }
	}
	for (j = 0; j < MAXTEANL; j++) {
	    if (display_slice(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
		setcolor(1);
		load_slice(gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	    }
	}
	for (j = 0; j < MAXTEANL; j++) {
	    if (display_slice(gno, n, TEANL, MAG, j) && object_isactive(TEANL, j)) {
		setcolor(1);
		load_slice(gno, n, flowf[j].grid, TEANL, MAG, j, step, t);
	    }
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    if (display_slice(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
		setcolor(1);
		load_slice(gno, n, -1, ADCIRC, ELEV, j, step, t);
	    }
	}
	for (j = 0; j < MAXADCIRC; j++) {
	    if (display_slice(gno, n, ADCIRC, MAG, j) && object_isactive(ADCIRC, FLOW, j)) {
		setcolor(1);
		load_slice(gno, n, -1, ADCIRC, MAG, j, step, t);
	    }
	}
/* reset original world */
	defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, islogx(gno), islogy(gno));
	viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
    }
}

int display_slice(int gno, int sno, int obj, int w1, int w2)
{
    switch (obj) {
    case BATH:
	return (g[gno].sbox[sno].bath[w1]);
	break;
    case TEANL:
	if (w1 == ELEV)
	    return (g[gno].sbox[sno].teanl[w2]);
	if (w1 == MAG)
	    return (g[gno].sbox[sno].teanlmag[w2]);
	break;
    case ADCIRC:
	if (w1 == ELEV)
	    return (g[gno].sbox[sno].adcirc[w2]);
	if (w1 == MAG)
	    return (g[gno].sbox[sno].adcircmag[w2]);
	break;
    case ADCIRC3DFLOW:
	if (w1 == ELEV)
	    return (g[gno].sbox[sno].adc3d[w2]);
	if (w1 == MAG)
	    return (g[gno].sbox[sno].adc3dflow[w2]);
	break;
    case ELA:
	return (g[gno].sbox[sno].ela[w1]);
	break;
    default:
	return 0;
    }
    return 0;
}

int read_slice(int sliceno, char *fname)
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

void load_slice(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double x, d;
    int i, ind;
    if (g[gno].sbox[n].npts > 1) {
	if (g[gno].sbox[n].type == LINE) {
	    double dx = (g[gno].sbox[n].x2 - g[gno].sbox[n].x1) / (g[gno].sbox[n].npts - 1);
	    double dy = (g[gno].sbox[n].y2 - g[gno].sbox[n].y1) / (g[gno].sbox[n].npts - 1);
	    if (g[gno].sbox[n].elist[0] == -1) {
		find_element(gridno, g[gno].sbox[n].x1, g[gno].sbox[n].y1, &ind);
	    } else {
		ind = g[gno].sbox[n].elist[0];
	    }
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1, g[gno].sbox[n].y1, 0);
		if (obj == BATH) {
		    d = -d;
		}
		my_move2(0.0, d);
	    } else {
		my_move2(0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x = hypot(i * dx, i * dy);
		    my_draw2(x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	} else {
	    x = 0.0;
	    if (g[gno].sbox[n].elist[0] == -1) {
		find_element(gridno, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], &ind);
	    } else {
		ind = g[gno].sbox[n].elist[0];
	    }
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], 0);
		if (obj == BATH) {
		    d = -d;
		}
		my_move2(0.0, d);
	    } else {
		my_move2(0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x += hypot(g[gno].sbox[n].sx[i] - g[gno].sbox[n].sx[i - 1], g[gno].sbox[n].sy[i] - g[gno].sbox[n].sy[i - 1]);
		    my_draw2(x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	}
    }
}

void load_flux_slice(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double x, d;
    int i, ind;
    if (g[gno].sbox[n].npts > 1) {
	if (g[gno].sbox[n].type == LINE) {
	    double dx = (g[gno].sbox[n].x2 - g[gno].sbox[n].x1) / (g[gno].sbox[n].npts - 1);
	    double dy = (g[gno].sbox[n].y2 - g[gno].sbox[n].y1) / (g[gno].sbox[n].npts - 1);
	    find_element(gridno, g[gno].sbox[n].x1, g[gno].sbox[n].y1, &ind);
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1, g[gno].sbox[n].y1, 0);
		if (obj == BATH) {
		    d = -d;
		}
		my_move2(0.0, d);
	    } else {
		my_move2(0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		find_element(gridno, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, &ind);
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x = hypot(i * dx, i * dy);
		    my_draw2(x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	} else {
	    x = 0.0;
	    find_element(gridno, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], &ind);
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], 0);
		if (obj == BATH) {
		    d = -d;
		}
		my_move2(0.0, d);
	    } else {
		my_move2(0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		find_element(gridno, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], &ind);
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x += hypot(g[gno].sbox[n].sx[i] - g[gno].sbox[n].sx[i - 1], g[gno].sbox[n].sy[i] - g[gno].sbox[n].sy[i - 1]);
		    my_draw2(x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	}
    }
}

int nslices(int gno, int n)
{
    int j, writecount = 0;
    for (j = 0; j < MAXGRIDS; j++) {
	if (display_slice(gno, n, BATH, j, 0) && object_isactive(GRID, j)) {
	    writecount++;
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
	    writecount++;
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, MAG, j) && object_isactive(TEANL, j)) {
	    writecount++;
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_slice(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
	    writecount++;
	}
    }
    for (j = 0; j < MAXADCIRC3D; j++) {
	if (display_slice(gno, n, ADCIRC3DFLOW, ELEV, j) && object_isactive(ADCIRC3DFLOW, j)) {
	    writecount++;
	}
    }
    for (j = 0; j < MAXELA; j++) {
	if (display_slice(gno, n, ELA, j, 0) && object_isactive(ELA, j)) {
	    writecount++;
	}
    }
    return writecount;
}

int write_slice_binary(FILE * fp, int gno, int n, int step, double t, int cnt)
{
    int i, j;
    double *x, *y;
    int writecount = 0;

    x = (double *) malloc(g[gno].sbox[n].npts * sizeof(double));
    y = (double *) malloc(g[gno].sbox[n].npts * sizeof(double));
    if (x == NULL || y == NULL) {
	return;
    }
    for (j = 0; j < MAXGRIDS; j++) {
	if (display_slice(gno, n, BATH, j, 0) && object_isactive(GRID, j)) {
	    writecount++;
	    get_slice(gno, n, j, BATH, 0, j, step, t, x, y);
	    fwrite(&g[gno].sbox[n].npts, sizeof(int), 1, fp);
	    fwrite(x, sizeof(double), g[gno].sbox[n].npts, fp);
	    fwrite(y, sizeof(double), g[gno].sbox[n].npts, fp);
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
	    writecount++;
	    get_slice(gno, n, flowf[j].grid, TEANL, ELEV, j, step, t, x, y);
	    fwrite(&g[gno].sbox[n].npts, sizeof(int), 1, fp);
	    fwrite(x, sizeof(double), g[gno].sbox[n].npts, fp);
	    fwrite(y, sizeof(double), g[gno].sbox[n].npts, fp);
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, MAG, j) && object_isactive(TEANL, j)) {
	    writecount++;
	    get_slice(gno, n, flowf[j].grid, TEANL, MAG, j, step, t, x, y);
	    fwrite(&g[gno].sbox[n].npts, sizeof(int), 1, fp);
	    fwrite(x, sizeof(double), g[gno].sbox[n].npts, fp);
	    fwrite(y, sizeof(double), g[gno].sbox[n].npts, fp);
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_slice(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
	    writecount++;
	    get_slice(gno, n, -1, ADCIRC, ELEV, j, step, t, x, y);
	    fwrite(&g[gno].sbox[n].npts, sizeof(int), 1, fp);
	    fwrite(x, sizeof(double), g[gno].sbox[n].npts, fp);
	    fwrite(y, sizeof(double), g[gno].sbox[n].npts, fp);
	}
    }
    free(x);
    free(y);
    return writecount;
}

int viewslice(FILE * fp, int gno, int n, int step, double t, int cnt)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    char buf[10];
    int writecount = 0;

    for (j = 0; j < MAXGRIDS; j++) {
	if (display_slice(gno, n, BATH, j, 0) && object_isactive(GRID, j)) {
	    writecount++;
	    write_slice(fp, gno, n, j, BATH, 0, j, step, t);
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, ELEV, j) && object_isactive(TEANL, j)) {
	    writecount++;
	    write_slice(fp, gno, n, flowf[j].grid, TEANL, ELEV, j, step, t);
	}
    }
    for (j = 0; j < MAXTEANL; j++) {
	if (display_slice(gno, n, TEANL, MAG, j) && object_isactive(TEANL, j)) {
	    writecount++;
	    write_slice(fp, gno, n, flowf[j].grid, TEANL, MAG, j, step, t);
	}
    }
    for (j = 0; j < MAXADCIRC; j++) {
	if (display_slice(gno, n, ADCIRC, ELEV, j) && object_isactive(ADCIRC, ELEV, j)) {
	    writecount++;
	    write_slice(fp, gno, n, -1, ADCIRC, ELEV, j, step, t);
	}
    }
    return writecount;
}

void write_slice(FILE * fp, int gno, int n, int gridno, int obj, int w1, int w2, int step, double t)
{
    double x, d;
    int i, ind;
    if (g[gno].sbox[n].npts > 1) {
	if (g[gno].sbox[n].type == LINE) {
	    double dx = (g[gno].sbox[n].x2 - g[gno].sbox[n].x1) / (g[gno].sbox[n].npts - 1);
	    double dy = (g[gno].sbox[n].y2 - g[gno].sbox[n].y1) / (g[gno].sbox[n].npts - 1);
	    if (g[gno].sbox[n].elist[0] == -1) {
		find_element(gridno, g[gno].sbox[n].x1, g[gno].sbox[n].y1, &ind);
	    } else {
		ind = g[gno].sbox[n].elist[0];
	    }
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1, g[gno].sbox[n].y1, 0);
		if (obj == BATH) {
		    d = -d;
		}
		fprintf(fp, "%lf %lf\n", 0.0, d);
	    } else {
		fprintf(fp, "%lf %lf\n", 0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x = hypot(i * dx, i * dy);
		    fprintf(fp, "%lf %lf\n", x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	    fprintf(fp, "&\n");
	} else {
	    x = 0.0;
	    if (g[gno].sbox[n].elist[0] == -1) {
		find_element(gridno, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], &ind);
	    } else {
		ind = g[gno].sbox[n].elist[0];
	    }
	    if (ind >= 0) {
		d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[0], g[gno].sbox[n].sy[0], 0);
		if (obj == BATH) {
		    d = -d;
		}
		fprintf(fp, "%lf %lf\n", 0.0, d);
	    } else {
		fprintf(fp, "%lf %lf\n", 0.0, 0.0);
	    }
	    for (i = 1; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x += hypot(g[gno].sbox[n].sx[i] - g[gno].sbox[n].sx[i - 1], g[gno].sbox[n].sy[i] - g[gno].sbox[n].sy[i - 1]);
		    fprintf(fp, "%lf %lf\n", x, d);
		} else {
		    printf("element not found\n");
		}
	    }
	}
    }
}

void get_slice(int gno, int n, int gridno, int obj, int w1, int w2, int step, double t, double *x1, double *y1)
{
    double x, d;
    int i, ind;
    if (g[gno].sbox[n].npts > 1) {
	if (g[gno].sbox[n].type == LINE) {
	    double dx = (g[gno].sbox[n].x2 - g[gno].sbox[n].x1) / (g[gno].sbox[n].npts - 1);
	    double dy = (g[gno].sbox[n].y2 - g[gno].sbox[n].y1) / (g[gno].sbox[n].npts - 1);
	    for (i = 0; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].x1 + i * dx, g[gno].sbox[n].y1 + i * dy, 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    x = hypot(i * dx, i * dy);
		    x1[i] = x;
		    y1[i] = d;
		} else {
		}
	    }
	} else {
	    x = 0.0;
	    for (i = 0; i < g[gno].sbox[n].npts; i++) {
		if (g[gno].sbox[n].elist[i] == -1) {
		    find_element(gridno, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], &ind);
		} else {
		    ind = g[gno].sbox[n].elist[i];
		}
		if (ind >= 0) {
		    d = eval_model(obj, w1, w2, step, ind, t, g[gno].sbox[n].sx[i], g[gno].sbox[n].sy[i], 0);
		    if (obj == BATH) {
			d = -d;
		    }
		    if (i > 0) {
			x += hypot(g[gno].sbox[n].sx[i] - g[gno].sbox[n].sx[i - 1], g[gno].sbox[n].sy[i] - g[gno].sbox[n].sy[i - 1]);
		    }
		    x1[i] = x;
		    y1[i] = d;
		} else {
		}
	    }
	}
    }
}

int display_slice2(DisplaySlice * s, int obj, int w1, int w2)
{
    switch (obj) {
    case BATH:
	return (s->bath[w1]);
	break;
    case TEANL:
	if (w1 == ELEV)
	    return (s->teanl[w2]);
	if (w1 == MAG)
	    return (s->teanlmag[w2]);
	break;
    case ADCIRC:
	if (w1 == ELEV)
	    return (s->adcirc[w2]);
	if (w1 == MAG)
	    return (s->adcircmag[w2]);
	break;
    case ADCIRC3DFLOW:
	if (w1 == ELEV)
	    return (s->adc3d[w2]);
	if (w1 == MAG)
	    return (s->adc3dflow[w2]);
	break;
    case ELA:
	return (s->ela[w1]);
	break;
    default:
	return 0;
    }
    return 0;
}
