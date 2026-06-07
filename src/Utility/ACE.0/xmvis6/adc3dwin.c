/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 *
 */

/*
 * adcirc3dwin.c - Sample data from points selected in the domain 
 * for display in a gadget. Specifically for SELFE.
 *
 */

#ifndef lint
static char RCSid[] = "$Id: adc3dwin.c,v 1.10 2006/05/17 22:38:27 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "elio.h"

int ReadNodeDataNew(int n, int startstep, int stopstep, int append);
int ReadXYDataNew(int n, int startstep, int stopstep, int append);
int QueryNodeData(char *fname, int node, int step);
int QueryXYData(char *fname, double x, double y, int step);

extern Widget app_shell;
extern XmStringCharSet charset;

typedef struct {
    Widget top;
    Widget *which;
    Widget datafile;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget append;
    Widget *loctype;
    Widget loc;
    Widget xy;
    Widget node;		/* end */
    Widget toggle;
    Widget togglemarker;
    Widget togglelines;
    Widget togglelabels;
    Widget nlabels;
    Widget tmin;
    Widget tmax;
    Widget cmin;
    Widget cmax;
    Widget w;
    Widget h;
    Widget *color;
    Widget *linew;
    Widget *fill;
    Widget *precx;
    Widget *precy;
    Widget *attach;
    Browser b1;			/* data */
    Browser b2;			/* Z-data */
} SetADC3DUI;

void create_setadc3d_frame(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc adc3d_accept_proc(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc adc3d_pickxy_proc(Widget w, XtPointer clid, XtPointer calld);
void clear_adc3d(int n);
XtCallbackProc adc3d_setflow3d_proc(Widget w, XtPointer clid, XtPointer calld);

static Widget xyloc = (Widget) NULL;

XtCallbackProc adc3d_pick_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *a = (SetADC3DUI *) clid;
    xyloc = a->loc;
    if (GetChoice(a->loctype)) {
	pick_adc3dxy();
    } else {
	pick_adc3dnode();
    }
}

XtCallbackProc adc3d_query_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *a = (SetADC3DUI *) clid;
    int n = GetChoice(a->which);
    int start = 0;
    int node = 0;
    char loc[1024], buf[1024];
    char datafile[2048];
    double x, y;
    strcpy(loc, (char *) xv_getstr(a->loc));
    strcpy(datafile, (char *) xv_getstr(a->datafile));
    if (!strlen(datafile)) {
	errwin("Need a data file name, operation canceled.");
	return;
    }
    strcpy(buf, (char *) xv_getstr(a->start));
    start = atoi(buf) - 1;
    if (start < 0) {
	errwin("Need to set start step >= 0 for Query, operation canceled.");
	return;
    }
    if (GetChoice(a->loctype)) {
	node = atoi(loc);
	QueryNodeData(datafile, node - 1, start);
    } else {
	sscanf(loc, "%lf %lf", &x, &y);
	QueryXYData(datafile, x, y, start);
    }
}

XtCallbackProc adc3d_view_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *a = (SetADC3DUI *) clid;
    xyloc = a->loc;
    if (GetChoice(a->loctype)) {
    } else {
    }
}

void set_adc3d_xy(double wx1, double wy1)
{
    char buf[256];
    if (xyloc != (Widget) NULL) {
/* should check here for element containment if possible
        int elem.
        if ((elem = ElioFindElementXY(&eh, x, y)) < 0) {
            sprintf(buf, "Unable to locate element\n");
            errwin(buf);
        }
*/
	sprintf(buf, "%.2lf %.2lf", wx1, wy1);
	xv_setstr(xyloc, buf);
    }
}

void set_adc3d_node(double wx1, double wy1)
{
    char buf[256];
    if (xyloc != (Widget) NULL) {
	sprintf(buf, "%.2lf %.2lf", wx1, wy1);
	xv_setstr(xyloc, buf);
    }
}

void query_adc3d_xy(double wx1, double wy1)
{
    char buf[256];
    if (xyloc != (Widget) NULL) {
	sprintf(buf, "%.2lf %.2lf", wx1, wy1);
	xv_setstr(xyloc, buf);
    }
}

void query_adc3d_node(double wx1, double wy1)
{
    char buf[256];
    if (xyloc != (Widget) NULL) {
	sprintf(buf, "%.2lf %.2lf", wx1, wy1);
	xv_setstr(xyloc, buf);
    }
}

XtCallbackProc adc3d_accept_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *a = (SetADC3DUI *) clid;
    char datafile[2048];
    char elevfile[2048];
    char loc[1024];
    double xp, yp;
    int append, sample, n;
    int start = -1, stop = 0;

    n = GetChoice(a->which);
    append = XmToggleButtonGetState(a->append);
    sample = XmToggleButtonGetState(a->sample);
    if (adc3d[n].active == OFF && append) {
	errwin("Can't append data to an inactive set.");
	return;
    }
    if (adc3d[n].active == ON && !append) {
	if (!yesno("This SELFE sample gadget is active, kill it?", " ", " YES ", " NO ")) {
	    return;
	}
	clear_adc3d(n);
    }
    adc3d[n].loctype = GetChoice(a->loctype);
    strcpy(loc, (char *) xv_getstr(a->loc));
    strcpy(datafile, (char *) xv_getstr(a->b1.file));
    strcpy(elevfile, (char *) xv_getstr(a->b2.file));

    if (strlen(datafile)) {
	strncpy(adc3d[n].datafile, datafile, 2047);
	strncpy(adc3d[n].elevfile, elevfile, 2047);
    } else {
	errwin("Need a data file name, operation canceled.");
	return;
    }
    if (sample) {
	strcpy(buf, (char *) xv_getstr(a->start));
	start = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(a->stop));
	stop = atoi(buf) - 1;
    } else {
	start = -1;
	stop = -1;
    }
    if (adc3d[n].loctype) {
	sscanf(loc, "%d", &adc3d[n].node);
	if (ReadNodeDataNew(n, start, stop, append) != 0) {
	    clear_adc3d(n);
	    errwin("Error reading data file, operation canceled.");
	    return;
	}
    } else {
	int ecode;
	sscanf(loc, "%lf %lf", &adc3d[n].x, &adc3d[n].y);
	if ((ecode = ReadXYDataNew(n, start, stop, append)) != 0) {
	    clear_adc3d(n);
	    if (ecode == 1) {
		errwin("Error reading data file, operation canceled.");
	    }
	    return;
	}
    }
    create_setadc3d_frame(NULL, NULL, NULL);
}

int ReadNodeDataNew(int n, int startstep, int stopstep, int append)
{
    FILE *fp;
    int i, j, k, ix, ictr;
    char buf[1024];
    int nsteps;
    int it, scnt = 0;
    float ttime, *dd;
    int node = adc3d[n].node - 1;
    int type;
    ElcircHeader h;
    int doelcirc = 0;

    if ((fp = fopen(adc3d[n].datafile, "r")) == NULL) {
	sprintf(buf, "ReadNodeDataNew(): Unable to open file %s for reading\n", adc3d[n].datafile);
	errwin(buf);
	return 1;
    }

    if (ElioGetHeader(adc3d[n].datafile, &h)) {
	fprintf(stderr, "ReadNodeDataNew(): Error reading header\n");
	return 1;
    }

    if (adc3d[n].node > h.np) {
	errwin("Node number greater than the number of nodes in the grid\n");
	return 1;
    }

    if (h.v == 3) {
	doelcirc = 1;
    } else {
	doelcirc = 0;
    }

    adc3d[n].nlevels = h.nvrt;

    nsteps = ElioGetNStepsInFile(adc3d[n].datafile, &h);

    adc3d[n].type = type = (h.ivs == 1 ? SCALAR : VECTOR);

    if (append) {
	if (type == VECTOR) {
	    adc3d[n].u = (double **) realloc(adc3d[n].u, (adc3d[n].nsteps + 1) * sizeof(double *));
	    adc3d[n].v = (double **) realloc(adc3d[n].v, (adc3d[n].nsteps + 1) * sizeof(double *));
	} else if (type == SCALAR) {
	    adc3d[n].sal = (double **) realloc(adc3d[n].sal, (adc3d[n].nsteps + 1) * sizeof(double *));
	} else {
	    errwin("ReadNodeDataNew(): Internal error, type not found");
	}
    } else {
	if (type == VECTOR) {
	    adc3d[n].u = (double **) malloc(sizeof(double *));
	    adc3d[n].v = (double **) malloc(sizeof(double *));
	} else if (type == SCALAR) {
	    adc3d[n].sal = (double **) malloc(sizeof(double *));
	} else {
	    errwin("ReadNodeDataNew(): Internal error, type not found");
	}
    }

    if (startstep >= 0 && stopstep >= startstep) {
	if (append) {
	    nsteps = adc3d[n].nsteps + stopstep - startstep + 1;
	    scnt = adc3d[n].nsteps;
	} else {
	    nsteps = adc3d[n].nsteps = stopstep - startstep + 1;
	    scnt = 0;
	}
    } else {
	startstep = 0;
	stopstep = nsteps - 1;
	if (append) {
	    scnt = adc3d[n].nsteps;
	} else {
	    adc3d[n].nsteps = nsteps;
	    scnt = 0;
	}
    }

    if (append) {
	adc3d[n].time = (double *) realloc(adc3d[n].time, (adc3d[n].nsteps + 1) * sizeof(double));
	adc3d[n].umin = adc3d[n].vmin = 10000.0;
	adc3d[n].umax = adc3d[n].vmax = -10000.0;
    } else {
	adc3d[n].time = (double *) malloc(sizeof(double));
	adc3d[n].umin = adc3d[n].vmin = 10000.0;
	adc3d[n].umax = adc3d[n].vmax = -10000.0;
    }

    dd = (float *) malloc(h.nvrt * 4 * h.ivs);

    for (j = startstep; j <= stopstep; j++) {
	for (i = 0; i < h.nvrt * h.ivs; i++) {
	    if (type == VECTOR) {
		dd[i] = 0.0;
	    } else {
		dd[i] = -99.0;
	    }
	}
	if (ElioGetNodeOld(fp, j, node, &h, &ttime, &it, dd)) {
	    break;
	}
/*
	printf("Time = %f it = %d, index = %d\n", ttime, it, j);
*/
	if (type == VECTOR) {
	    adc3d[n].u[scnt] = (double *) malloc(h.nvrt * sizeof(double));
	    adc3d[n].v[scnt] = (double *) malloc(h.nvrt * sizeof(double));
	    for (i = 0; i < h.nvrt; i++) {
		adc3d[n].u[scnt][i] = dd[2 * i];
		adc3d[n].v[scnt][i] = dd[2 * i + 1];
		if (dd[2 * i] >= -98.0) {	/* u */
		    adc3d[n].umin = (dd[2 * i] < adc3d[n].umin) ? dd[2 * i] : adc3d[n].umin;
		    adc3d[n].umax = (dd[2 * i] > adc3d[n].umax) ? dd[2 * i] : adc3d[n].umax;
		}
		if (dd[2 * i + 1] >= -98.0) {	/* v */
		    adc3d[n].vmin = (dd[2 * i + 1] < adc3d[n].vmin) ? dd[2 * i + 1] : adc3d[n].vmin;
		    adc3d[n].vmax = (dd[2 * i + 1] > adc3d[n].vmax) ? dd[2 * i + 1] : adc3d[n].vmax;
		}
	    }
	    adc3d[n].time[scnt] = ttime;
	    adc3d[n].time = (double *) realloc(adc3d[n].time, (scnt + 2) * sizeof(double));
	    adc3d[n].u = (double **) realloc(adc3d[n].u, (scnt + 2) * sizeof(double *));
	    adc3d[n].v = (double **) realloc(adc3d[n].v, (scnt + 2) * sizeof(double *));
	    scnt++;
	} else if (type == SCALAR) {
	    adc3d[n].sal[scnt] = (double *) malloc(h.nvrt * sizeof(double));
	    for (i = 0; i < h.nvrt; i++) {
		adc3d[n].sal[scnt][i] = dd[i];
		if (dd[i] >= -98.0) {
		    adc3d[n].umin = (dd[i] < adc3d[n].umin) ? dd[i] : adc3d[n].umin;
		    adc3d[n].umax = (dd[i] > adc3d[n].umax) ? dd[i] : adc3d[n].umax;
		}
	    }
	    adc3d[n].time[scnt] = ttime;
	    adc3d[n].time = (double *) realloc(adc3d[n].time, (scnt + 2) * sizeof(double));
	    adc3d[n].sal = (double **) realloc(adc3d[n].sal, (scnt + 2) * sizeof(double *));
	    scnt++;
	}
    }
    fclose(fp);
    free(dd);
    ElioFreeHeader(&h);
    adc3d[n].nsteps = scnt;
    adc3d[n].active = ON;
    set_clock(0, adc3d[n].time[0], adc3d[n].time[scnt - 1], adc3d[n].time[1] - adc3d[n].time[0], scnt);
    load_clock(ADCIRC3DFLOW, n, 0);
    return 0;
}

int find_elementxy(int nmel, double *xord, double *yord, int *icon[3], double x, double y);
int eval_flowxy(double h1, double h2, double h3, double u1, double v1, double u2, double v2, double u3, double v3, double *uret, double *vret);
int eval_scalarxy(double h1, double h2, double h3, double u1, double u2, double u3, double *uret);
void getcoefs(double x1, double y1, double x2, double y2, double x3, double y3, double xp, double yp, double *h1, double *h2, double *h3);

int ReadXYDataNew(int n, int startstep, int stopstep, int append)
{
    FILE *gridfp, *fp, *zfp;
    int i, j, k, ix, ictr, bind[4], sind[4];
    int nsteps;
    double hh[4];
    double x = adc3d[n].x;
    double y = adc3d[n].y;
    int elem, nn[3];
    char buf[1024];
    int it, scnt = 0;
    float ttime, *dd, *z;
    long offset;
    double uret, vret, zret;
    int type = adc3d[n].type;
    int slabsize, np, nels;
    int ftype;
    char id[4];
    ElcircHeader eh, zh;
    int doelcirc, retval = 0;;
//printf("ReadXYDataNew %d %d %d\n", startstep, stopstep, append);

    if ((fp = fopen(adc3d[n].datafile, "r")) == NULL) {
	sprintf(buf, "ReadXYDataNew(): Unable to open file %s for reading\n", adc3d[n].datafile);
	errwin(buf);
	return 1;
    }
    if ((zfp = fopen(adc3d[n].elevfile, "r")) == NULL) {
	sprintf(buf, "ReadXYDataNew(): Unable to open file %s for reading\n", adc3d[n].datafile);
	errwin(buf);
	return 1;
    }
    if (ElioGetHeader(adc3d[n].datafile, &eh)) {
	fprintf(stderr, "ReadXYDataNew(): Error reading header\n");
	return 1;
    }

    if (eh.v == 3) {
	doelcirc = 1;
    } else {
	doelcirc = 0;
	if (ElioGetHeader(adc3d[n].elevfile, &zh)) {
	    fprintf(stderr, "ReadXYDataNew(): Error reading header\n");
	    return 1;
	}
    }

/* find the element that contains the point */
    if ((elem = ElioFindElementXY(&eh, x, y)) < 0) {
	sprintf(buf, "Unable to locate element\n");
	errwin(buf);
	retval = 2;
	goto out;
    }
/* get the coefficients for interpolation */
    ElioGetCoefficients(&eh, elem, x, y, hh);
    adc3d[n].nlevels = eh.nvrt;
    nsteps = ElioGetNStepsInFile(adc3d[n].datafile, &eh);
    adc3d[n].type = type = (eh.ivs == 1 ? SCALAR : VECTOR);

    if (append) {
	if (type == VECTOR) {
	    adc3d[n].u = (double **) realloc(adc3d[n].u, (adc3d[n].nsteps + 1) * sizeof(double *));
	    adc3d[n].v = (double **) realloc(adc3d[n].v, (adc3d[n].nsteps + 1) * sizeof(double *));
	    adc3d[n].z = (double **) realloc(adc3d[n].z, (adc3d[n].nsteps + 1) * sizeof(double *));
	} else if (type == SCALAR) {
	    adc3d[n].sal = (double **) realloc(adc3d[n].sal, (adc3d[n].nsteps + 1) * sizeof(double *));
	    adc3d[n].z = (double **) realloc(adc3d[n].z, (adc3d[n].nsteps + 1) * sizeof(double *));
	} else {
	    errwin("ReadXYDataNew(): Internal error, type not found");
	    retval = 1;
	    goto out;
	}
    } else {
	if (type == VECTOR) {
	    adc3d[n].u = (double **) malloc(sizeof(double *));
	    adc3d[n].v = (double **) malloc(sizeof(double *));
	    adc3d[n].z = (double **) malloc(sizeof(double *));
	} else if (type == SCALAR) {
	    adc3d[n].sal = (double **) malloc(sizeof(double *));
	    adc3d[n].z = (double **) malloc(sizeof(double *));
	} else {
	    errwin("ReadXYDataNew(): Internal error, type not found");
	    retval = 1;
	    goto out;
	}
    }

    if (startstep >= 0 && stopstep >= startstep) {
	if (append) {
	    nsteps = adc3d[n].nsteps + stopstep - startstep + 1;
	    scnt = adc3d[n].nsteps;
	} else {
	    nsteps = adc3d[n].nsteps = stopstep - startstep + 1;
	    scnt = 0;
	}
    } else {
	startstep = 0;
	stopstep = nsteps - 1;
	if (append) {
	    scnt = adc3d[n].nsteps;
	} else {
	    adc3d[n].nsteps = nsteps;
	    scnt = 0;
	}
    }

    if (append) {
	adc3d[n].time = (double *) realloc(adc3d[n].time, (adc3d[n].nsteps + 1) * sizeof(double));
	adc3d[n].umin = adc3d[n].vmin = adc3d[n].zmin = 10000.0;
	adc3d[n].umax = adc3d[n].vmax = adc3d[n].zmax = -10000.0;
    } else {
	adc3d[n].time = (double *) malloc(sizeof(double));
	adc3d[n].umin = adc3d[n].vmin = adc3d[n].zmin = 10000.0;
	adc3d[n].umax = adc3d[n].vmax = adc3d[n].zmax = -10000.0;
    }

    nels = eh.ivs;		/* 2 if vector, 1 if scalar */

    dd = (float *) malloc(eh.ivs * eh.nvrt * 4);
    z = (float *) malloc(eh.nvrt * 4);
    for (j = startstep; j <= stopstep; j++) {
	for (i = 0; i < eh.nvrt * eh.ivs; i++) {
	    if (type == VECTOR) {
		dd[i] = 0.0;
	    } else {
		dd[i] = -99.0;
	    }
	}
	if (type == VECTOR) {
	    if (ElioGetXYData2(fp, j, elem, &eh, hh, &ttime, &it, bind, sind, dd)) {
		goto out;
	    }
	    if (!doelcirc && ElioGetXYData2(zfp, j, elem, &zh, hh, &ttime, &it, bind, sind, z)) {
		goto out;
	    }
	    adc3d[n].u[scnt] = malloc(eh.nvrt * sizeof(double));
	    adc3d[n].v[scnt] = malloc(eh.nvrt * sizeof(double));
	    adc3d[n].z[scnt] = malloc(eh.nvrt * sizeof(double));
	    for (i = 0; i < eh.nvrt; i++) {
		uret = adc3d[n].u[scnt][i] = dd[2 * i];
		vret = adc3d[n].v[scnt][i] = dd[2 * i + 1];
		if (doelcirc) {
		    zret = adc3d[n].z[scnt][i] = i + 1;
		} else {
		    zret = adc3d[n].z[scnt][i] = z[i];
		}
		adc3d[n].zmin = (zret < adc3d[n].zmin) ? zret : adc3d[n].zmin;
		adc3d[n].zmax = (zret > adc3d[n].zmax) ? zret : adc3d[n].zmax;
		if (uret >= -98.0) {	/* u */
		    adc3d[n].umin = (uret < adc3d[n].umin) ? uret : adc3d[n].umin;
		    adc3d[n].umax = (uret > adc3d[n].umax) ? uret : adc3d[n].umax;
		}
		if (vret >= -98.0) {	/* u */
		    adc3d[n].vmin = (vret < adc3d[n].vmin) ? vret : adc3d[n].vmin;
		    adc3d[n].vmax = (vret > adc3d[n].vmax) ? vret : adc3d[n].vmax;
		}
	    }
	    /* output results */
	    adc3d[n].time[scnt] = ttime;
	    adc3d[n].time = realloc(adc3d[n].time, (scnt + 2) * sizeof(double));
	    adc3d[n].u = (double **) realloc(adc3d[n].u, (scnt + 2) * sizeof(double *));
	    adc3d[n].v = (double **) realloc(adc3d[n].v, (scnt + 2) * sizeof(double *));
	    adc3d[n].z = (double **) realloc(adc3d[n].z, (scnt + 2) * sizeof(double *));
	    scnt++;
	} else if (type == SCALAR) {
	    if (ElioGetXYData2(fp, j, elem, &eh, hh, &ttime, &it, bind, sind, dd)) {
		goto out;
	    }
	    if (!doelcirc && ElioGetXYData2(zfp, j, elem, &zh, hh, &ttime, &it, bind, sind, z)) {
		goto out;
	    }
	    adc3d[n].sal[scnt] = malloc(eh.nvrt * sizeof(double));
	    adc3d[n].z[scnt] = malloc(eh.nvrt * sizeof(double));
	    for (i = 0; i < eh.nvrt; i++) {
		uret = adc3d[n].sal[scnt][i] = dd[i];
		if (doelcirc) {
		    zret = adc3d[n].z[scnt][i] = i + 1;
		} else {
		    zret = adc3d[n].z[scnt][i] = z[i];
		}
		adc3d[n].zmin = (zret < adc3d[n].zmin) ? zret : adc3d[n].zmin;
		adc3d[n].zmax = (zret > adc3d[n].zmax) ? zret : adc3d[n].zmax;
		if (uret >= -98.0) {	/* u */
		    adc3d[n].umin = (uret < adc3d[n].umin) ? uret : adc3d[n].umin;
		    adc3d[n].umax = (uret > adc3d[n].umax) ? uret : adc3d[n].umax;
		}
	    }
	    adc3d[n].time[scnt] = ttime;
	    adc3d[n].time = realloc(adc3d[n].time, (scnt + 2) * sizeof(double));
	    adc3d[n].sal = (double **) realloc(adc3d[n].sal, (scnt + 2) * sizeof(double *));
	    adc3d[n].z = (double **) realloc(adc3d[n].z, (scnt + 2) * sizeof(double *));
	    scnt++;
	}
    }
  out:;
    fclose(fp);
    free(dd);
    ElioFreeHeader(&eh);
    adc3d[n].nsteps = scnt;
    adc3d[n].active = ON;
    set_clock(0, adc3d[n].time[0], adc3d[n].time[scnt - 1], adc3d[n].time[1] - adc3d[n].time[0], scnt);
    load_clock(ADCIRC3DFLOW, n, 0);
    //printf("%lf %lf\n", adc3d[n].umin, adc3d[n].umax);
    return retval;
}

int find_elementxy(int nmel, double *xord, double *yord, int *icon[3], double x, double y)
{
    int i;
    int n0, n1, n2;
    for (i = 0; i < nmel; i++) {
	n0 = icon[0][i];
	n1 = icon[1][i];
	n2 = icon[2][i];
	if (belel(x, y, xord[n0], yord[n0], xord[n1], yord[n1], xord[n2], yord[n2])) {
	    return i;
	}
    }
    return -1;
}

int eval_flowxy(double h1, double h2, double h3, double u1, double v1, double u2, double v2, double u3, double v3, double *uret, double *vret)
{
    *uret = (h1 * u1 + h2 * u2 + h3 * u3);
    *vret = (h1 * v1 + h2 * v2 + h3 * v3);
    return 1;
}

int eval_scalarxy(double h1, double h2, double h3, double u1, double u2, double u3, double *uret)
{
    *uret = (h1 * u1 + h2 * u2 + h3 * u3);
    return 1;
}

void getcoefs(double x1, double y1, double x2, double y2, double x3, double y3, double xp, double yp, double *h1, double *h2, double *h3)
{
    double aum, ado, atr, bum, bdo, btr, cum, cdo, ctr;
    double c1, c2, c3, arei;
    aum = x2 * y3 - x3 * y2;
    bum = y2 - y3;
    cum = x3 - x2;
    ado = x3 * y1 - x1 * y3;
    bdo = y3 - y1;
    cdo = x1 - x3;
    atr = x1 * y2 - x2 * y1;
    btr = y1 - y2;
    ctr = x2 - x1;
    arei = 1.0 / (aum + ado + atr);
    *h3 = (atr + btr * xp + ctr * yp) * arei;
    *h2 = (ado + bdo * xp + cdo * yp) * arei;
    *h1 = 1.0 - *h2 - *h3;
}

void clear_adc3d(int n)
{
    int i, j;
    if (adc3d[n].active == ON && adc3d[n].nsteps > 0) {
	for (j = 0; j < adc3d[n].nsteps; j++) {
	    if (adc3d[n].sal != NULL) {
		cxfree(adc3d[n].sal[j]);
	    }
	    if (adc3d[n].u != NULL) {
		cxfree(adc3d[n].u[j]);
	    }
	    if (adc3d[n].v != NULL) {
		cxfree(adc3d[n].v[j]);
	    }
	    if (adc3d[n].w != NULL) {
		cxfree(adc3d[n].w[j]);
	    }
	}
	cxfree(adc3d[n].sal);
	cxfree(adc3d[n].u);
	cxfree(adc3d[n].v);
	cxfree(adc3d[n].w);
	cxfree(adc3d[n].time);
    }
    adc3d[n].active = OFF;
    adc3d[n].sal = NULL;
    adc3d[n].u = NULL;
    adc3d[n].v = NULL;
    adc3d[n].w = NULL;
    adc3d[n].nsteps = 0;
    adc3d[n].type = VECTOR;
    adc3d[n].loctype = 0;
    adc3d[n].node = -1;
    adc3d[n].x = 0;
    adc3d[n].y = 0;
    adc3d[n].datafile[0] = 0;
}

void drawflow3d(int gno, int n, int step, double t)
{
    int i, j, ix1, iy1, ix2, iy2, ie;
    double e, vx1, vx2, vy1, vy2;
    double wx1, wx2, wy1, wy2;
    char buf[10];
    int imax = adc3d[n].nsteps;
    int savestep = step;
    extern int mapisolconc[];
    strcpy(statusstr, "Start in drawflow3d()");

    if (step < 0)
	return;
/* don't allow for step > nsteps or crash */
    if (step >= adc3d[n].nsteps) {
	step = adc3d[n].nsteps - 1;
    }
    if (g[gno].flow3d[n].display == ON) {
	setcolor(g[gno].flow3d[n].p.color);
	my_move2(adc3d[n].x, adc3d[n].y);
	my_draw2(g[gno].flow3d[n].locx, g[gno].flow3d[n].locy);
    }
    if (g[gno].flow3d[n].display_marker == ON) {
/* redefine the world and viewport */
	switch (g[gno].flow3d[n].attach) {
	case 0:
	    world2view(g[gno].flow3d[n].locx, g[gno].flow3d[n].locy, &vx1, &vy1);
	    vx2 = vx1 + g[gno].flow3d[n].vx;
	    vy2 = vy1 + g[gno].flow3d[n].vy;
	    break;
	case 1:
	    world2view(g[gno].flow3d[n].locx, g[gno].flow3d[n].locy, &vx2, &vy1);
	    vx1 = vx2 - g[gno].flow3d[n].vx;
	    vy2 = vy1 + g[gno].flow3d[n].vy;
	    break;
	case 2:
	    world2view(g[gno].flow3d[n].locx, g[gno].flow3d[n].locy, &vx1, &vy2);
	    vx2 = vx1 + g[gno].flow3d[n].vx;
	    vy1 = vy2 - g[gno].flow3d[n].vy;
	    break;
	case 3:
	    world2view(g[gno].flow3d[n].locx, g[gno].flow3d[n].locy, &vx2, &vy2);
	    vx1 = vx2 - g[gno].flow3d[n].vx;
	    vy1 = vy2 - g[gno].flow3d[n].vy;
	    break;
	}
	setcolor(g[gno].flow3d[n].p.color);
	defineworld(g[gno].flow3d[n].wx1, g[gno].flow3d[n].wy1, g[gno].flow3d[n].wx2, g[gno].flow3d[n].wy2, 0, 0);
	viewport(vx1, vy1, vx2, vy2);

	setcolor(g[gno].flow3d[n].p.fillcol);
	fillrectcolor(g[gno].flow3d[n].wx1, g[gno].flow3d[n].wy1, g[gno].flow3d[n].wx2, g[gno].flow3d[n].wy2);

	setcolor(g[gno].flow3d[n].p.color);
	rect(g[gno].flow3d[n].wx1, g[gno].flow3d[n].wy1, g[gno].flow3d[n].wx2, g[gno].flow3d[n].wy2);
	if (g[gno].flow3d[n].nlabels == 2) {
	    rectstr(g[gno].flow3d[n].wx1, g[gno].flow3d[n].wy1, g[gno].flow3d[n].wx2, g[gno].flow3d[n].wy2, g[gno].flow3d[n].precx, g[gno].flow3d[n].precy);
	} else if (g[gno].flow3d[n].nlabels > 2) {
	    double dx, dy, loc;
	    int i;
	    char label[32], format[32];
	    strcpy(format, "%.*lf");
	    dy = (g[gno].flow3d[n].wy2 - g[gno].flow3d[n].wy1) / (g[gno].flow3d[n].nlabels - 1);
	    for (i = 0; i < g[gno].flow3d[n].nlabels; i++) {
		sprintf(label, format, 1, g[gno].flow3d[n].wy1 + i * dy);
		writestr(g[gno].flow3d[n].wx1, g[gno].flow3d[n].wy1 + i * dy, 0, 1, label);
	    }
	}

	setlinewidth(g[gno].flow3d[n].p.linew);
	if (object_isactive(ADCIRC3DFLOW, n)) {
	    setcolor(1);
	    if (step > imax) {
		step = imax;
	    }
	    for (i = 0; i < adc3d[n].nlevels; i++) {
		switch (adc3d[n].type) {
		case VECTOR:
		    if (adc3d[n].u[step][i] > -98.0 && adc3d[n].v[step][i] > -98.0) {
			velplt(adc3d[n].time[step], adc3d[n].z[step][i], adc3d[n].u[step][i], adc3d[n].v[step][i], g[gno].vl.unitfac * g[gno].vl.scale);
		    }
		    break;
		case MAGNITUDE:	/* Magnitudes of velocities */
		case SCALAR:	/* salinity */
		    {
			int j;
			double bx[4], by[4], mag;
			double delta = (adc3d[n].time[1] - adc3d[n].time[0]);
			for (j = 0; j < step; j++) {
			    if (i + 1 < adc3d[n].nlevels) {
				bx[0] = adc3d[n].time[j];
				by[0] = adc3d[n].z[j][i];
				bx[1] = adc3d[n].time[j + 1];
				by[1] = adc3d[n].z[j + 1][i];
				bx[2] = adc3d[n].time[j + 1];
				by[2] = adc3d[n].z[j + 1][i + 1];
				bx[3] = adc3d[n].time[j];
				by[3] = adc3d[n].z[j][i + 1];
				if (adc3d[n].type == MAGNITUDE) {
				    mag = hypot(adc3d[n].u[j][i], adc3d[n].v[j][i]);
				    fillsquare(mag, bx, by, g[cg].velmagip, mapisolconc, g[cg].velmagip.nisol);
				} else {
				    fillsquare(adc3d[n].sal[j][i], bx, by, g[cg].salip, mapisolconc, g[cg].salip.nisol);
				}
			    } else {
				bx[0] = adc3d[n].time[j];
				by[0] = adc3d[n].z[j][i];
				bx[1] = adc3d[n].time[j + 1];
				by[1] = adc3d[n].z[j + 1][i];
			    }
			    if (g[gno].flow3d[n].display_lines == ON) {
				setcolor(1);
				my_move2(bx[0], by[0]);
				my_draw2(bx[1], by[1]);
			    }
			}
		    }
		    break;
		}
	    }
	}
	setcolor(g[gno].flow3d[n].p.color);
	my_move2(t, g[gno].flow3d[n].wy1);
	my_draw2(t, g[gno].flow3d[n].wy2);
	defineworld(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2, islogx(gno), islogy(gno));
	viewport(g[gno].v.xv1, g[gno].v.yv1, g[gno].v.xv2, g[gno].v.yv2);
    }
    setlinewidth(1);
    setcolor(1);
    strcpy(statusstr, "Done in drawflow3d()");
}

void update_adc3d(SetADC3DUI * ui);
XtCallbackProc setadc3d_accept_proc(Widget w, XtPointer clid, XtPointer calld);
void do_autoscaleadc3d(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc adc3d_isolines_proc(Widget w, XtPointer clid, XtPointer calld);

int curadc3d = 0;

XtCallbackProc set_adc3dlocCB(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *ui = (SetADC3DUI *) clid;
    curadc3d = GetChoice(ui->which);
    set_adc3dloc();
}

XtCallbackProc set_curadc3d(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *ui = (SetADC3DUI *) clid;
    curadc3d = GetChoice(ui->which);
    update_adc3d(ui);
}

void create_setadc3d_frame(Widget w, XtPointer clid, XtPointer calld)
{
    static SetADC3DUI ui;
    Arg al[8];
    Widget panel, wbut, fr, bb, sep, rc, rc1, rc2;
    int ac = 0, i;
    setistop();
    strcpy(statusstr, "Start create_setadc3d_frame()");
    writelogfile(statusstr);
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "SELFE Samples setup", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "setadc3drc", NULL, 0);
	ui.which = CreatePanelChoice1(panel, "SELFE Sample marker: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < 10; i++) {
	    XtAddCallback(ui.which[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadc3d, (XtPointer) & ui);
	}
	//ui.datafile = CreateTextItem2(panel, 15, "Data file: ");

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b1.file = CreateTextItem2(rc, 20, "Data file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b1);
	XtManageChild(rc);

	XtVaCreateManagedWidget("Use _elev.61 for Z-level _zcor.63 for Sigma.", xmLabelWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b2.file = CreateTextItem2(rc, 20, "Elev/Zcor file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b2);
	XtManageChild(rc);

	ui.append = XtVaCreateManagedWidget("Append to current data set", xmToggleButtonWidgetClass, panel, NULL);

	ui.sample = XtVaCreateManagedWidget("Sample steps:", xmToggleButtonWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.start = CreateTextItem2(rc, 4, "Start: ");
	xv_setstr(ui.start, "1");
	ui.stop = CreateTextItem2(rc, 4, "Stop: ");
	xv_setstr(ui.stop, "1");
	ui.skip = CreateTextItem2(rc, 4, "Skip: ");
	xv_setstr(ui.skip, "1");
	XtManageChild(rc);
/*
	ui.toggle = XtVaCreateManagedWidget("Data from node", xmToggleButtonWidgetClass, panel, NULL);
	ui.toggle = XtVaCreateManagedWidget("Data from XY", xmToggleButtonWidgetClass, panel, NULL);
*/
	ui.loctype = CreatePanelChoice1(panel, "Location type: ", 3, "X Y", "Node", 0, 0);
	ui.loc = CreateTextItem2(panel, 15, "X Y or Node # = ");
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc1 = XmCreateRowColumn(panel, "rc1", NULL, 0);
	XtVaSetValues(rc1, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Read data", xmPushButtonWidgetClass, rc1, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_accept_proc, &ui);
	wbut = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc1, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_pick_proc, &ui);
	wbut = XtVaCreateManagedWidget("*Query...", xmPushButtonWidgetClass, rc1, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_query_proc, &ui);
	wbut = XtVaCreateManagedWidget("*View(ACE/gr)...", xmPushButtonWidgetClass, rc1, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_view_proc, &ui);
	XtManageChild(rc1);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	XtVaCreateManagedWidget("Display:", xmLabelWidgetClass, panel, NULL);
	rc1 = XmCreateRowColumn(panel, "rc1", NULL, 0);
	XtVaSetValues(rc1, XmNorientation, XmHORIZONTAL, NULL);
	ui.toggle = XtVaCreateManagedWidget("Location", xmToggleButtonWidgetClass, rc1, NULL);
	ui.togglemarker = XtVaCreateManagedWidget("Marker", xmToggleButtonWidgetClass, rc1, NULL);
	XtManageChild(rc1);
	rc1 = XmCreateRowColumn(panel, "rc1", NULL, 0);
	XtVaSetValues(rc1, XmNorientation, XmHORIZONTAL, NULL);
	ui.togglelines = XtVaCreateManagedWidget("Lines", xmToggleButtonWidgetClass, rc1, NULL);
	ui.togglelabels = XtVaCreateManagedWidget("Labels", xmToggleButtonWidgetClass, rc1, NULL);
	ui.nlabels = CreateTextItem2(rc1, 4, "# labels: ");
	XtManageChild(rc1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);


	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	ui.tmin = CreateTextItem2(rc2, 10, "X min: ");
	ui.tmax = CreateTextItem2(rc2, 10, "X max: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	ui.cmin = CreateTextItem2(rc2, 10, "Y min: ");
	ui.cmax = CreateTextItem2(rc2, 10, "Y max: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	ui.w = CreateTextItem2(rc2, 6, "Width: ");
	ui.h = CreateTextItem2(rc2, 6, "Height: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_autoscaleadc3d, (XtPointer) & ui);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "setadc3drc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	XtManageChild(rc);
	XtManageChild(bb);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) setadc3d_accept_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_adc3dlocCB, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Props...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_setflow3d_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adc3d_isolines_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
    update_adc3d(&ui);
    strcpy(statusstr, "End create_setadc3d_frame()");
    writelogfile(statusstr);
}

void do_autoscaleadc3d(Widget w, XtPointer clid, XtPointer calld)
{
    char buf[256];
    SetADC3DUI *ui = (SetADC3DUI *) clid;
    strcpy(statusstr, "Start do_autoscaleadc3d()");
    writelogfile(statusstr);
    if (timeclock.nsteps) {
	g[cg].flow3d[curadc3d].wx1 = timeclock.t[0];
	g[cg].flow3d[curadc3d].wx2 = timeclock.t[timeclock.nsteps - 1];
// Autoscale on depths.
	//g[cg].flow3d[curadc3d].wy1 = 0;
	//g[cg].flow3d[curadc3d].wy2 = adc3d[curadc3d].nlevels + 1;
	g[cg].flow3d[curadc3d].wy1 = adc3d[curadc3d].zmin;
	g[cg].flow3d[curadc3d].wy2 = adc3d[curadc3d].zmax;

	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wx1);
	XmTextSetString(ui->tmin, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wx2);
	XmTextSetString(ui->tmax, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wy1);
	XmTextSetString(ui->cmin, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wy2);
	XmTextSetString(ui->cmax, buf);
    }
    strcpy(statusstr, "End do_autoscaleadc3d()");
    writelogfile(statusstr);
}

XtCallbackProc adc3d_isolines_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *ui = (SetADC3DUI *) clid;
    int n = curadc3d;
    strcpy(statusstr, "Start adc3d_isolines_proc()");
    writelogfile(statusstr);
/* get min/max */
    switch (adc3d[n].type) {
    case VECTOR:
	g[cg].velmagip.cmin = adc3d[n].umin;
	g[cg].velmagip.cmax = adc3d[n].umax;
	create_isolines_popup("Isolines of magnitudes", &g[cg].velmagip, 1);
	break;
    case SCALAR:
	g[cg].salip.cmin = adc3d[n].umin;
	g[cg].salip.cmax = adc3d[n].umax;
	create_isolines_popup("Isolines of salinity", &g[cg].salip, 1);
	break;
    }
    strcpy(statusstr, "End adc3d_isolines_proc()");
    writelogfile(statusstr);
}

void update_adc3d(SetADC3DUI * ui)
{
    char buf[256];
    int n;
    if (ui->top) {
	SetChoice(ui->which, curadc3d);
	n = curadc3d;
	SetChoice(ui->loctype, adc3d[n].loctype);
	if (adc3d[n].loctype) {
	    sprintf(buf, "%d", adc3d[n].node);
	    xv_setstr(ui->loc, buf);
	} else {
	    sprintf(buf, "%.1lf %.1lf", adc3d[n].x, adc3d[n].y);
	    xv_setstr(ui->loc, buf);
	}
	XmToggleButtonSetState(ui->toggle, g[cg].flow3d[curadc3d].display == ON, False);
	XmToggleButtonSetState(ui->togglemarker, g[cg].flow3d[curadc3d].display_marker == ON, False);
	XmToggleButtonSetState(ui->togglelines, g[cg].flow3d[curadc3d].display_lines == ON, False);
	XmToggleButtonSetState(ui->togglelabels, g[cg].flow3d[curadc3d].display_labels == ON, False);
	if (g[cg].flow3d[curadc3d].nlabels < 1) {
	    g[cg].flow3d[curadc3d].nlabels = 2;
	}
	sprintf(buf, "%d", g[cg].flow3d[curadc3d].nlabels);
	xv_setstr(ui->nlabels, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wx1);
	XmTextSetString(ui->tmin, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wx2);
	XmTextSetString(ui->tmax, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wy1);
	XmTextSetString(ui->cmin, buf);
	sprintf(buf, "%.2lf", g[cg].flow3d[curadc3d].wy2);
	XmTextSetString(ui->cmax, buf);
	sprintf(buf, "%.3lf", g[cg].flow3d[curadc3d].vx);
	XmTextSetString(ui->w, buf);
	sprintf(buf, "%.3lf", g[cg].flow3d[curadc3d].vy);
	XmTextSetString(ui->h, buf);
    }
}

XtCallbackProc setadc3d_accept_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetADC3DUI *ui = (SetADC3DUI *) clid;
    int h;
    double x, y, locx, locy;
    x = y = locx = locy = 0.0;
    strcpy(statusstr, "Start setadc3d_accept_proc()");
    writelogfile(statusstr);

    curadc3d = h = GetChoice(ui->which);

    g[cg].flow3d[curadc3d].wx1 = atof((char *) xv_getstr(ui->tmin));
    g[cg].flow3d[curadc3d].wx2 = atof((char *) xv_getstr(ui->tmax));
    g[cg].flow3d[curadc3d].wy1 = atof((char *) xv_getstr(ui->cmin));
    g[cg].flow3d[curadc3d].wy2 = atof((char *) xv_getstr(ui->cmax));
    g[cg].flow3d[curadc3d].vx = atof((char *) xv_getstr(ui->w));
    g[cg].flow3d[curadc3d].vy = atof((char *) xv_getstr(ui->h));

    g[cg].flow3d[curadc3d].display = XmToggleButtonGetState(ui->toggle) ? ON : OFF;
    g[cg].flow3d[curadc3d].display_marker = XmToggleButtonGetState(ui->togglemarker) ? ON : OFF;
    g[cg].flow3d[curadc3d].display_lines = XmToggleButtonGetState(ui->togglelines) ? ON : OFF;
    g[cg].flow3d[curadc3d].display_labels = XmToggleButtonGetState(ui->togglelabels) ? ON : OFF;
    g[cg].flow3d[curadc3d].nlabels = atoi((char *) xv_getstr(ui->nlabels));
    if (g[cg].flow3d[curadc3d].nlabels < 2) {
	g[cg].flow3d[curadc3d].nlabels = 2;
    }
    strcpy(statusstr, "End setadc3d_accept_proc()");
    writelogfile(statusstr);
}

void register_adc3dvelmarker(int gno, int n, double wx1, double wy1)
{
    int ind, k;
    double get_current_time(void);
    ind = adc3d[n].node - 1;
    g[gno].flow3d[n].locx = wx1;
    g[gno].flow3d[n].locy = wy1;
    drawflow3d(gno, n, get_current_step(), get_current_time());
}

typedef struct {
    Widget top;
    Widget *color;
    Widget *linew;
    Widget *fill;
    Widget *precx;
    Widget *precy;
    Widget *attach;
} SetFlow3dPropsUI;

void create_flow3dprops_popup(char *name, Display3dFlow * f);

XtCallbackProc adc3d_setflow3d_proc(Widget w, XtPointer clid, XtPointer calld)
{
    int n = curadc3d;
    create_flow3dprops_popup("Properties for SELFE samples", &g[cg].flow3d[n]);
}

XtCallbackProc setflow3dprops_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SetFlow3dPropsUI *ui = (SetFlow3dPropsUI *) clid;
    strcpy(statusstr, "Start setflow3dprops_proc()");
    writelogfile(statusstr);
    g[cg].flow3d[curadc3d].p.color = GetChoice(ui->color);
    g[cg].flow3d[curadc3d].p.linew = GetChoice(ui->linew);
    g[cg].flow3d[curadc3d].p.fillcol = GetChoice(ui->fill);
    g[cg].flow3d[curadc3d].precx = GetChoice(ui->precx) - 1;
    g[cg].flow3d[curadc3d].precy = GetChoice(ui->precy) - 1;
    g[cg].flow3d[curadc3d].attach = GetChoice(ui->attach);
    strcpy(statusstr, "End setflow3dprops_proc()");
    writelogfile(statusstr);
}

void update_flow3dprops(SetFlow3dPropsUI * ui)
{
    char buf[256];
    int n;
    if (ui->top) {
	SetChoice(ui->color, g[cg].flow3d[curadc3d].p.color);
	SetChoice(ui->fill, g[cg].flow3d[curadc3d].p.fillcol);
	SetChoice(ui->linew, g[cg].flow3d[curadc3d].p.linew);
	SetChoice(ui->precx, g[cg].flow3d[curadc3d].precx + 1);
	SetChoice(ui->precy, g[cg].flow3d[curadc3d].precy + 1);
	SetChoice(ui->attach, g[cg].flow3d[curadc3d].attach);
    }
}

void create_flow3dprops_popup(char *name, Display3dFlow * f)
{
    static SetFlow3dPropsUI ui;
    Arg al[8];
    Widget panel, wbut, fr, bb, sep, rc, rc1, rc2;
    int ac = 0, i;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "SELFE Samples props", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "flow3dprops", NULL, 0);
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, NULL);
	ui.color = CreateColorChoice(bb, "Frame color:", 1);
	ui.linew = CreatePanelChoice2(bb, "Frame line width: ", 2, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

	ui.fill = CreateColorChoice(bb, "Fill color:", 1);
	ui.precx = CreatePanelChoice2(bb, "Precision X labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	ui.precy = CreatePanelChoice2(bb, "Precision Y labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	ui.attach = CreatePanelChoice2(bb, "Attach to: ", 1, 5, "SW", "SE", "NW", "NE", 0, 0);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) setflow3dprops_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
    update_flow3dprops(&ui);
}

int QueryNodeData(char *fname, int node, int step)
{
    FILE *fp;
    int i;
    int it, scnt = 0;
    float ttime, *dd = NULL;
    ElcircHeader h;
    int retval = 1, bi, ti;
    char buf[1024];

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "QueryNodeData(): Unable to open file %s for reading\n", fname);
	errwin(buf);
	return retval;
    }
    if (ElioGetHeader(fname, &h)) {
	errwin("QueryNodeData(): Error reading header\n");
	goto out;
    }
    dd = (float *) malloc(h.nvrt * 4 * h.ivs);
    if (node >= h.np) {
	errwin("QueryNodeData(): Node number greater than the number of nodes in the grid.\n");
	goto out;
    }
    scnt = ElioGetNStepsInFile(fname, &h);
    if (step > scnt) {
	errwin("QueryNodeData(): Step to read greater than the number of steps in file.\n");
	goto out;
    }
    bi = h.bi[node];
    ti = ElioGetNodeSurfaceIndex(fp, node, step, &h);
    if (ElioGetNode(fp, step, node, &h, &ttime, &it, &bi, &ti, dd)) {
	errwin("QueryNodeData(): Call to ElioGetNode() failed to read data.\n");
	goto out;
    }
    create_monitor_frame();
    sprintf(buf, "Time: %f, Node: %d, Bottom index: %d, Surface index: %d\n", ttime, node + 1, bi + 1, ti);
    stufftext(buf, 1);
    for (i = bi; i < h.nvrt; i++) {
	if (i == ti) {
	    sprintf(buf, "-- Surface --\n");
	    stufftext(buf, 0);
	}
	if (h.ivs == 2) {
	    sprintf(buf, "%d %f %f\n", i + 1, dd[2 * (i - bi)], dd[2 * (i - bi) + 1]);
	} else {
	    sprintf(buf, "%d %f\n", i + 1, dd[i - bi]);
	}
	stufftext(buf, 0);
    }
    retval = 1;
  out:;
    fclose(fp);
    ElioFreeHeader(&h);
    if (dd != NULL) {
	free(dd);
    }
    return retval;
}

int QueryXYData(char *fname, double x, double y, int step)
{
    FILE *fp;
    int i, j, k;
    int nsteps;
    double hh[4], c[4];
    int elem, nn[4];
    char buf[1024];
    int bi, ti, it, scnt = 0, bind[4], sind[4];
    float ttime, *dd;
    double uret, vret;
    ElcircHeader eh;

    if ((fp = fopen(fname, "r")) == NULL) {
	sprintf(buf, "QueryXYData(): Unable to open file %s for reading.\n", fname);
	errwin(buf);
	return 1;
    }
    if (ElioGetHeader(fname, &eh)) {
	fprintf(stderr, "QueryXYData(): Error reading header.\n");
	return 1;
    }

/* find the element that contains the point */
    if ((elem = ElioFindElementXY(&eh, x, y)) < 0) {
	printf("QueryXYData: Unable to locate element.\n");
	return 1;
    }
/* get the nodes for the element */
    ElioGetCoefficients(&eh, elem, x, y, hh);
    nsteps = ElioGetNStepsInFile(fname, &eh);

    dd = (float *) malloc(eh.ivs * eh.nvrt * 4);
    scnt = 0;
    if (!ElioGetXYData2(fp, step, elem, &eh, hh, &ttime, &it, bind, sind, dd)) {
	create_monitor_frame();
	sprintf(buf, "Element: %d, Time: %f, X,Y: %lf, %lf\n", elem + 1, ttime, x, y);
	stufftext(buf, 1);
	for (i = 0; i < eh.nvrt; i++) {
	    if (eh.ivs == 2) {
		sprintf(buf, "%d %.5lf %.5lf\n", i + 1, dd[2 * i], dd[2 * i + 1]);
	    } else {
		sprintf(buf, "%d %.5lf\n", i + 1, dd[i]);
	    }
	    stufftext(buf, 0);
	}
    }
    free(dd);
    return 0;
}

/*
int ViewNodeData(char *fname, int node, int step)
{
    char tbuf[1024];
    char tmpname[1024];
    strcpy(tbuf, "/tmp/ACEadcirc3dnodeXXXXXX");
    mkstemp(tbuf);
    tmpname = tbuf;
    fp = fopen(tmpname, "wb");
    if (fp != NULL) {
        fprintf(fp, "@title \"SELFE vertical values at nodes\"\n");
        fprintf(fp, "@subtitle \"Step %d, Time = %lf\"\n", c + 1, t);
        fprintf(fp, "@xaxis label \"Time\"\n");
        fprintf(fp, "@yaxis label \"Level\"\n");
        for (i = 0; i < flowt[curadcirc].grid.nmnp; i++) {
            fprintf(fp, "%d %lf\n", i + 1, flowt[curadcirc].f[c].e[i]);
        }
        fclose(fp);
        sprintf(combuf, "/usr/local/ace/bin/xmgr5 -remove %s &", tmpname);
        system(combuf);
    } else {
        errwin("ViewNodeData(): Unable to open file");
    }
}
*/
