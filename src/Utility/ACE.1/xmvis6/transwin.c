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
 * transwin.c - transects from elcirc 3d files.
 *
 */

#ifndef lint
static char RCSid[] = "$Id: transwin.c,v 1.17 2004/08/05 21:32:18 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "elio.h"

void reversearray(int n, float *tmpu);

extern Widget app_shell;
extern XmStringCharSet charset;

static int ntrans;
static double xt[1000], yt[1000];

/*
    char tbuf[256];
    strcpy(tbuf, "/tmp/ACEtransXXXXXX");
    mkstemp(tbuf);
*/

typedef struct {
    Widget top;
    Widget *trans;
    Widget *flow;
    Widget *graph;
    Widget data;
    Widget pick;
    Widget transname;
    Widget *type;
    Widget *transtype;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget missing;
    Widget missingval;
    Widget append;
    Browser b;
} TransectUI;

void create_settrans_frame(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc trans_accept_proc(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc trans_pickxy_proc(Widget w, XtPointer clid, XtPointer calld);
void clear_trans(int n);
XtCallbackProc trans_setflow3d_proc(Widget w, XtPointer clid, XtPointer calld);

static Widget xyloc = (Widget) NULL;

XtCallbackProc trans_pickxy_proc(Widget w, XtPointer clid, XtPointer calld)
{
    TransectUI *a = (TransectUI *) clid;
    ntrans = 0;
    pick_transxy();
}

XtCallbackProc trans_accept_proc(Widget w, XtPointer clid, XtPointer calld)
{
    TransectUI *a = (TransectUI *) clid;
    int i;
    int datatype;
    int flowno;
    int gno;
    int sample, start = -1, stop = 1, skip = 1;
    int missing;
    int append = 0;
    int pick = 0;
    double missingval;
    char tbuf[256];
    FILE *fp;

    curtrans = GetChoice(a->trans);
    trans[curtrans].flowno = GetChoice(a->flow);
/* here check for active grid and flow # */
    trans[curtrans].gno = GetChoice(a->graph);
    datatype = GetChoice(a->type);
    if (pick = XmToggleButtonGetState(a->pick)) {
	if (ntrans == 0) {
	    errwin("No points picked for transect.");
	    return;
	}
	strcpy(tbuf, "/tmp/ACEtransXXXXXX");
	mkstemp(tbuf);
	strcpy(trans[curtrans].transname, tbuf);
	if ((fp = fopen(tbuf, "w")) != NULL) {
	    fprintf(fp, "Temporary transect file\n%d\n", ntrans);
	    for (i = 0; i < ntrans; i++) {
		fprintf(fp, "%d %lf %lf 1.0\n", i + 1, xt[i], yt[i]);
	    }
	    fclose(fp);
	} else {
	    errwin("Unable to open temporary transect file.");
	    return;
	}
    } else {
	strcpy(trans[curtrans].transname, (char *) xv_getstr(a->transname));
    }
    append = XmToggleButtonGetState(a->append);
    sample = XmToggleButtonGetState(a->sample);
    if (sample) {
	strcpy(buf, (char *) xv_getstr(a->start));
	start = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(a->stop));
	stop = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(a->skip));
	skip = atoi(buf);
    } else {
	start = -1;
	skip = 1;
    }
    trans[curtrans].start = start;
    trans[curtrans].stop = stop;
    trans[curtrans].skip = skip;
    /* no slot for missing data in transect struct */
    missing = XmToggleButtonGetState(a->missing);
    if (missing) {
	strcpy(buf, (char *) xv_getstr(a->missingval));
	missingval = atof(buf);
    } else {
	missing = 0;
    }
    set_wait_cursor(a->top);
    switch (datatype) {
    case 0:
	strcpy(trans[curtrans].salname, (char *) xv_getstr(a->data));
	trans[curtrans].type = SALINITY;
	ReadTransect(&trans[curtrans], append);
	set_clock(0, flowt[trans[curtrans].flowno].start, flowt[trans[curtrans].flowno].stop, flowt[trans[curtrans].flowno].step, flowt[trans[curtrans].flowno].nsteps);
	load_clock(ADCIRC, trans[curtrans].flowno);
	break;
    case 1:
	strcpy(trans[curtrans].uvname, (char *) xv_getstr(a->data));
	trans[curtrans].type = FLOW;
	ReadTransect(&trans[curtrans], append);
	set_clock(0, flowt[trans[curtrans].flowno].start, flowt[trans[curtrans].flowno].stop, flowt[trans[curtrans].flowno].step, flowt[trans[curtrans].flowno].nsteps);
	load_clock(ADCIRC, trans[curtrans].flowno);

	break;
    }
    unset_wait_cursor(a->top);
    XtUnmanageChild(a->top);
    create_sadcirc_frame();
}

void create_settrans_frame(Widget w, XtPointer clid, XtPointer calld)
{
    static TransectUI ui;
    Arg al[8];
    Widget panel, wbut, fr, bb, sep, rc, rc1, rc2;
    int ac = 0, i;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "SELFE Samples setup", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "settransrc", NULL, 0);

	ui.trans = CreatePanelChoice1(panel, "Use transect: ", 6, "1", "2", "3", "4", "5", 0, 0);
	ui.flow = CreatePanelChoice1(panel, "Read to ADCIRC 2d flow/elevations: ", 6, "1", "2", "3", "4", "5", 0, 0);
	ui.graph = CreatePanelChoice1(panel, "Display in graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	ui.type = CreatePanelChoice1(panel, "Quantity to read: ", 3, "3d scalars", "3d vectors", 0, 0);
	ui.data = CreateTextItem2(panel, 10, "Data file: ");

	rc1 = XmCreateRowColumn(panel, "rc1", NULL, 0);
	XtVaSetValues(rc1, XmNorientation, XmHORIZONTAL, NULL);
	ui.pick = XtVaCreateManagedWidget("Pick transect interactively:", xmToggleButtonWidgetClass, rc1, NULL);
	wbut = XtVaCreateManagedWidget("Pick transect", xmPushButtonWidgetClass, rc1, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) trans_pickxy_proc, &ui);
	XtManageChild(rc1);
	ui.transname = CreateTextItem2(panel, 10, "Read transect data file: ");
	ui.append = XtVaCreateManagedWidget("Append to active data set", xmToggleButtonWidgetClass, panel, NULL);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
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

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.missing = XtVaCreateManagedWidget("Missing data:", xmToggleButtonWidgetClass, rc, NULL);
	ui.missingval = CreateTextItem2(rc, 4, "Value: ");
	xv_setstr(ui.missingval, "-99.0");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) trans_accept_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
}

void add_trans_xy(double x, double y)
{
    xt[ntrans] = x;
    yt[ntrans] = y;
    ntrans++;
}

int find_elementxy(int nmel, double *xord, double *yord, int *icon[3], double x, double y);
int eval_flowxy(double h1, double h2, double h3, double u1, double v1, double u2, double v2, double u3, double v3, double *uret, double *vret);
int eval_scalarxy(double h1, double h2, double h3, double u1, double u2, double u3, double *uret);
void getcoefs(double x1, double y1, double x2, double y2, double x3, double y3, double xp, double yp, double *h1, double *h2, double *h3);

int ReadTransect(Transect * t, int append)
{
    FILE *fp;
    int ftype;
    char *fname;
    ElcircHeader h;
    if (t->type == FLOW) {
	fname = t->uvname;
    } else if (t->type == SALINITY) {
	fname = t->salname;
    } else {
	fprintf(stderr, "ReadTransect(): Transect data type not found (%d).\n", t->type);
	return (1);
    }
    if ((fp = fopen(fname, "rb")) == NULL) {
	fprintf(stderr, "ReadTransect(): Unable to open flow file %s for reading\n", fname);
	return (1);
    }
    ftype = ElioGetFileType(fp);
    fclose(fp);
    if (ftype == -1) {
	errwin("ReadTransect(): Can't read old format with this version of ACE/vis...");
	return 1;
    } else {
	if (ElioGetHeader(fname, &h)) {
	    fprintf(stderr, "ReadTransect(): Error returned from ElioGetHeader().\n");
	    return 1;
	}
	if (ReadTransectNew(h, t, append)) {
	    fprintf(stderr, "ReadTransect(): Error returned from ReadTransectNew().\n");
	} else {
	    settime_info_start(cg, h.start_time);
	}
	ElioFreeHeader(&h);
    }
    return 0;
}

void PrintTransect(Transect t)
{
    int i;
    printf("Graph# %d, flow# %d\n", t.gno, t.flowno);
    printf("Transect data source %d\n", t.type);
    printf("Grid file name: %s\n", t.gname);
    printf("Salinity file name: %s\n", t.salname);
    printf("Fort64 file name: %s\n", t.uvname);
    printf("Transect file name: %s\n", t.transname);
    printf("Transect type: %d\n", t.transtype);
    printf("Number of points: %d\n", t.npts);
    printf("Line: (%lf, %lf) <-> (%lf, %lf)\n", t.x1, t.y1, t.x2, t.y2);
    for (i = 0; i < t.npts; i++) {
	printf("%d: (%lf, %lf)\n", i, t.x[i], t.y[i]);
    }
}

void AddTransNXY(Transect * t, int n, double x, double y)
{
    if (t->npts == 0) {
	t->x = malloc(sizeof(double));
	t->y = malloc(sizeof(double));
    } else {
	t->x = realloc(t->x, sizeof(double) * (t->npts + 1));
	t->y = realloc(t->y, sizeof(double) * (t->npts + 1));
    }
    t->x[t->npts] = x;
    t->y[t->npts] = y;
    t->npts++;
}

void AddTransNode(Transect * t, int gridno, int node)
{
    double x, y;
    x = grid[gridno].xord[node];
    y = grid[gridno].yord[node];
    if (t->npts == 0) {
	t->x = malloc(sizeof(double));
	t->y = malloc(sizeof(double));
    } else {
	t->x = realloc(t->x, sizeof(double) * (t->npts + 1));
	t->y = realloc(t->y, sizeof(double) * (t->npts + 1));
    }
    t->x[t->npts] = x;
    t->y[t->npts] = y;
    t->npts++;
}

int ReadTransectNew(ElcircHeader eh, Transect * t, int append)
{
    FILE *f64i, *fp;
    int i, ii, j, k, ix, ictr, node, nmel, nmnp, elem, cnt;
    char *aix;
    int nnodes, *elcheck, numtmp, ntype;
    int np, nsteps;
    char buf[1024], fname[1024];
    float tmpu[400], ttime;
    double *xord, *yord, *depth;
    int *icon[3];
    int it;
    float z[400];
    double x, y;
    double x1, y1, x2, y2, delta;
    int nn[3], n1, n2, n3;
    long offset;
    double uret, vret;
    double xtmp, ytmp;

    double umin, umax;
    double tumin, tumax;
    double vmin, vmax;
    double tvmin, tvmax;

    int flowno = t->flowno;

    int first = 1;

    FILE *fpreg;
    int doLine;
    double *xb, *yb;
    double deltax, deltay;

/* new grid stuff */
    int nmnp2, nmel2, *icon2[3];
    double *xt, *yt, *dt;

/* elements to use for interpolation */
    int nels, bind[4], sind[4];
    int *ellist;
    int *nlist[4];
    double *h[4], hh[4];
    float uu[4][400], vv[4][400];
    double sum = 0.0;

    int scnt, start, stop, skip;
    int *readdata;

    Grid *g;

/* get the points for the transect and create the transect grid */
    if (t->transtype == 0) {
	doLine = 0;
	if ((fpreg = fopen(t->transname, "r")) == NULL) {
	    fprintf(stderr, "Unable to open region file %s\n", t->transname);
	    return (1);
	}
	if (fgets(buf, 255, fpreg) == NULL) {
	    fprintf(stderr, "Error reading 1st line of region file %s\n", t->transname);
	    return (1);
	}
	if (fgets(buf, 255, fpreg) == NULL) {
	    fprintf(stderr, "Error reading 2nd line of region file %s\n", t->transname);
	    return (1);
	}
	sscanf(buf, "%d", &t->npts);
	if (t->npts < 2) {
	    sprintf(buf, "Insufficient number of points for transect, %d.\n", t->npts);
	    errwin(buf);
	}
	xb = (double *) malloc(t->npts * sizeof(double));
	yb = (double *) malloc(t->npts * sizeof(double));
	ellist = (int *) malloc(t->npts * sizeof(int));
	nlist[0] = (int *) malloc(t->npts * sizeof(int));
	nlist[1] = (int *) malloc(t->npts * sizeof(int));
	nlist[2] = (int *) malloc(t->npts * sizeof(int));
	nlist[3] = (int *) malloc(t->npts * sizeof(int));
	h[0] = (double *) malloc(t->npts * sizeof(double));
	h[1] = (double *) malloc(t->npts * sizeof(double));
	h[2] = (double *) malloc(t->npts * sizeof(double));
	h[3] = (double *) malloc(t->npts * sizeof(double));
	cnt = 0;		/* counter for valid stations */
	for (i = 0; i < t->npts; i++) {
	    if (fgets(buf, 255, fpreg) == NULL) {
		sprintf(buf, "Error reading line %d of region file %s\n", i + 4, t->transname);
		errwin(buf);
		free(xb);
		free(yb);
		free(ellist);
		free(nlist[0]);
		free(nlist[1]);
		free(nlist[2]);
		free(nlist[3]);
		free(h[0]);
		free(h[1]);
		free(h[2]);
		free(h[3]);
		return (1);
	    }
	    sscanf(buf, "%*d %lf %lf", &xtmp, &ytmp);
	    if ((elem = ElioFindElementXY(&eh, xtmp, ytmp)) < 0) {
		sprintf(buf, "Unable to locate element for station %d at (%lf, %lf)\n", i + 1, xtmp, ytmp);
		errwin(buf);
	    } else {
		ellist[cnt] = elem;
		xb[cnt] = xtmp;
		yb[cnt] = ytmp;
		cnt++;
	    }
	}
	fclose(fpreg);
	t->npts = cnt;		/* reset in case any stations were skipped */
	if (t->npts == 0) {
	    sprintf(buf, "No Points! Transect computation cancelled.");
	    errwin(buf);
	    free(xb);
	    free(yb);
	    free(ellist);
	    free(nlist[0]);
	    free(nlist[1]);
	    free(nlist[2]);
	    free(nlist[3]);
	    free(h[0]);
	    free(h[1]);
	    free(h[2]);
	    free(h[3]);
	    return 1;
	}
    } else if (t->transtype == 1 && t->npts > 2) {
	doLine = 1;
	deltax = (t->x2 - t->x1) / (t->npts - 1);
	deltay = (t->y2 - t->y1) / (t->npts - 1);
	xb = (double *) malloc(t->npts * sizeof(double));
	yb = (double *) malloc(t->npts * sizeof(double));
	ellist = (int *) malloc(t->npts * sizeof(int));
	nlist[0] = (int *) malloc(t->npts * sizeof(int));
	nlist[1] = (int *) malloc(t->npts * sizeof(int));
	nlist[2] = (int *) malloc(t->npts * sizeof(int));
	nlist[3] = (int *) malloc(t->npts * sizeof(int));
	h[0] = (double *) malloc(t->npts * sizeof(double));
	h[1] = (double *) malloc(t->npts * sizeof(double));
	h[2] = (double *) malloc(t->npts * sizeof(double));
	h[3] = (double *) malloc(t->npts * sizeof(double));
	cnt = 0;
	for (i = 0; i < t->npts; i++) {
	    xb[cnt] = t->x1 + i * deltax;
	    yb[cnt] = t->y1 + i * deltay;
	    if ((elem = ElioFindElementXY(&eh, xb[cnt], yb[cnt])) < 0) {
		sprintf(buf, "Unable to locate element for point %d: %lf, %lf\n", i + 1, xb[cnt], yb[cnt]);
		errwin(buf);
	    } else {
		ellist[cnt] = elem;
		cnt++;
	    }
	}
	t->npts = cnt;		/* reset in case any stations were skipped */
	if (t->npts == 0) {
	    sprintf(buf, "No Points! Transect computation cancelled.");
	    errwin(buf);
	    free(xb);
	    free(yb);
	    free(ellist);
	    free(nlist[0]);
	    free(nlist[1]);
	    free(nlist[2]);
	    free(nlist[3]);
	    free(h[0]);
	    free(h[1]);
	    free(h[2]);
	    free(h[3]);
	    return 1;
	}
    } else if (t->transtype == 2) {
	xb = (double *) malloc(t->npts * sizeof(double));
	yb = (double *) malloc(t->npts * sizeof(double));
	ellist = (int *) malloc(t->npts * sizeof(int));
	nlist[0] = (int *) malloc(t->npts * sizeof(int));
	nlist[1] = (int *) malloc(t->npts * sizeof(int));
	nlist[2] = (int *) malloc(t->npts * sizeof(int));
	nlist[3] = (int *) malloc(t->npts * sizeof(int));
	h[0] = (double *) malloc(t->npts * sizeof(double));
	h[1] = (double *) malloc(t->npts * sizeof(double));
	h[2] = (double *) malloc(t->npts * sizeof(double));
	h[3] = (double *) malloc(t->npts * sizeof(double));
	cnt = 0;
	for (i = 0; i < t->npts; i++) {
	    xb[cnt] = t->x[i];
	    yb[cnt] = t->y[i];
	    if ((elem = ElioFindElementXY(&eh, xb[cnt], yb[cnt])) < 0) {
		sprintf(buf, "Unable to locate element for point %d: %lf, %lf\n", i + 1, xb[cnt], yb[cnt]);
		errwin(buf);
	    } else {
		ellist[cnt] = elem;
		cnt++;
	    }
	}
	t->npts = cnt;		/* reset in case any stations were skipped */
	if (t->npts == 0) {
	    sprintf(buf, "No Points! Transect computation cancelled.");
	    errwin(buf);
	    free(xb);
	    free(yb);
	    free(ellist);
	    free(nlist[0]);
	    free(nlist[1]);
	    free(nlist[2]);
	    free(nlist[3]);
	    free(h[0]);
	    free(h[1]);
	    free(h[2]);
	    free(h[3]);
	    return 1;
	}
    }

    t->x = xb;
    t->y = yb;

    nmel = eh.ne;
    nmnp = eh.np;

    xord = (double *) malloc(nmnp * sizeof(double));
    yord = (double *) malloc(nmnp * sizeof(double));
    depth = (double *) malloc(nmnp * sizeof(double));
    for (i = 0; i < nmnp; i++) {
	depth[i] = eh.d[i];
	if (eh.v == 4) {
	    xord[i] = eh.x[i];
	    yord[i] = eh.y[i];
	} else {
	    xord[i] = eh.x[i];
	    yord[i] = eh.y[i];
	}
    }
    icon[0] = eh.icon[0];
    icon[1] = eh.icon[1];
    icon[2] = eh.icon[2];
    icon[3] = eh.icon[3];

/* create transect grid */
    nmel2 = 2 * (t->npts - 1) * (eh.nvrt - 1);
    nmnp2 = t->npts * eh.nvrt;
    xt = (double *) malloc(nmnp2 * sizeof(double));
    yt = (double *) malloc(nmnp2 * sizeof(double));
    dt = (double *) malloc(nmnp2 * sizeof(double));
    icon2[0] = (int *) malloc(nmel2 * sizeof(int));
    icon2[1] = (int *) malloc(nmel2 * sizeof(int));
    icon2[2] = (int *) malloc(nmel2 * sizeof(int));
    nmnp2 = 0;
    for (j = 0; j < eh.nvrt; j++) {
	for (i = 0; i < t->npts; i++) {
	    if (i == 0) {
		sum = xt[nmnp2] = 0.0;
	    } else {
		sum = sum + hypot((xb[i] - xb[i - 1]), (yb[i] - yb[i - 1]));
		xt[nmnp2] = sum;
	    }
	    yt[nmnp2] = eh.zcor[eh.nvrt - j - 1] - eh.zmsl;
	    dt[nmnp2] = j * t->npts + i;
	    nmnp2++;
	}
    }
    nmel2 = 0;
    for (j = 0; j < eh.nvrt - 1; j++) {
	for (i = 0; i < t->npts - 1; i++) {
	    icon2[0][nmel2] = j * t->npts + i;
	    icon2[1][nmel2] = j * t->npts + i + 1;
	    icon2[2][nmel2++] = (j + 1) * t->npts + i;
	    icon2[0][nmel2] = j * t->npts + i + 1;
	    icon2[1][nmel2] = (j + 1) * t->npts + i + 1;
	    icon2[2][nmel2++] = (j + 1) * t->npts + i;
	}
    }
    g = &flowt[flowno].g;
    AllocateGrid(g, nmel2, nmnp2);
    for (i = 0; i < nmnp2; i++) {
	g->xord[i] = xt[i];
	g->yord[i] = yt[i];
	g->depth[i] = dt[i];
    }
    for (i = 0; i < nmel2; i++) {
	g->icon[i].nl[0] = icon2[0][i], g->icon[i].nl[1] = icon2[1][i], g->icon[i].nl[2] = icon2[2][i];
	g->icon[i].type = g->icon[i].nn = g->icon[i].ngeom = 3;
    }
    dminmax(g->nmnp, g->xord, &g->xmin, &g->xmax);
    dminmax(g->nmnp, g->yord, &g->ymin, &g->ymax);
    dminmax(g->nmnp, g->depth, &g->dmin, &g->dmax);

/* free stuff used to create new grid */
    free(xt);
    free(yt);
    free(dt);
    free(icon2[0]);
    free(icon2[1]);
    free(icon2[2]);

/* locate elements and compute weights */
    for (i = 0; i < t->npts; i++) {
	/* printf("point %d: %lf, %lf\n", i, xb[i], yb[i]); */
	if ((elem = ElioFindElementXY(&eh, xb[i], yb[i])) < 0) {
	    fprintf(stderr, "Unable to locate element for point %d: %lf, %lf\n", i + 1, xb[i], yb[i]);
	    sprintf(buf, "Unable to locate element\n");
	    errwin(buf);
	    ellist[i] = -1;
	} else {

	    if (debuglevel == 13) {
		fprintf(stderr, "El %d point %d: %lf, %lf\n", elem, i, xb[i], yb[i]);
	    }

	    ellist[i] = elem;
	    nn[0] = icon[0][elem];
	    nn[1] = icon[1][elem];
	    nn[2] = icon[2][elem];
	    nlist[0][i] = nn[0];
	    nlist[1][i] = nn[1];
	    nlist[2][i] = nn[2];
	    if (eh.etype[elem] == 4) {
		nn[3] = icon[3][elem];
		nlist[3][i] = nn[3];
	    }
	    /*
	       avgd[i] = 0.0;
	       for (j=0;j<eh.etype[elem];j++) {
	       avgd[i] += eh.d[nn[j]];
	       }
	     */

	    /* get the coefficients for interpolation */
	    ElioGetCoefficients(&eh, elem, xb[i], yb[i], hh);
	    h[0][i] = hh[0];
	    h[1][i] = hh[1];
	    h[2][i] = hh[2];
	    h[3][i] = hh[3];
	}
    }

    if (!append) {
	for (j = 0; j < flowt[flowno].nsteps; j++) {
	    if (flowt[flowno].f) {
		cxfree(flowt[flowno].f[j].u);
		cxfree(flowt[flowno].f[j].v);
		cxfree(flowt[flowno].f[j].e);
	    }
	}
	flowt[flowno].nsteps = 0;
	if (flowt[flowno].f) {
	    cxfree(flowt[flowno].f);
	    flowt[flowno].f = NULL;
	}
	flowt[flowno].f = NULL;
	flowt[flowno].nsteps = 0;
	flowt[flowno].npts = nmnp2;
	flowt[flowno].f = (Flow_gen *) malloc(sizeof(Flow_gen));
    } else {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + 1) * sizeof(Flow_gen));
    }

/* free orginal grid */
    free(xord);
    free(yord);
    free(depth);

    if (t->type == FLOW) {
	strcpy(fname, t->uvname);
	if ((f64i = fopen(t->uvname, "r")) == NULL) {
	    fprintf(stderr, "Unable to open flow file %s for reading\n", t->uvname);
	    return (1);
	}
    } else if (t->type == SALINITY) {
	strcpy(fname, t->salname);
	if ((f64i = fopen(t->salname, "r")) == NULL) {
	    fprintf(stderr, "Unable to open salinity file %s for reading\n", t->salname);
	    return (1);
	}
    }

/* compute number of time steps in file */

    nsteps = ElioGetNStepsInFile(fname, &eh);

/* TODO start here by not using t->start, etc.
    if (t->start < 0) {
	t->start = 0;
	t->stop = nsteps - 1;
    } else {
	if (t->start >= nsteps) {
	}
	if (t->stop >= nsteps) {
	    t->stop = nsteps - 1;
	}
    }
*/

    start = t->start;
    stop = t->stop;
    skip = t->skip;
    if (start < 0) {
	start = 0;
	stop = nsteps - 1;
	skip = 1;
    } else {
	if (start >= nsteps) {
	}
	if (stop >= nsteps) {
	    stop = nsteps - 1;
	}
    }


    if (debuglevel == 13) {
	fprintf(stderr, "Start, stop, skip %d %d %d\n", t->start, t->stop, t->skip);
    }
    readdata = (int *) malloc(nsteps * sizeof(int));
    for (i = 0; i < nsteps; i++) {
	readdata[i] = 0;
    }
    scnt = 0;

    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    /* fprintf(stderr, "Start, stop, skip, scnt %d %d %d %d\n", t->start, t->stop, t->skip, scnt); &/

/*
    for (i = t->start; i <= t->stop; i += t->skip) {
	readdata[i] = 1;
	scnt++;
    }
    fprintf(stderr, "Start, stop, skip, nsteps, scnt %d %d %d %d %d\n", t->start, t->stop, t->skip, nsteps, scnt);
*/
    if (scnt == 0) {
	errwin("No complete time step available, cancelling read");
	flowt[flowno].active = OFF;
	free(readdata);
	return 1;
    }
    flowt[flowno].npts = nnodes;

    if (append) {
	flowt[flowno].f = (Flow_gen *) realloc(flowt[flowno].f, (flowt[flowno].nsteps + scnt) * sizeof(Flow_gen));
	scnt = flowt[flowno].nsteps;
    } else {
	flowt[flowno].f = (Flow_gen *) malloc(scnt * sizeof(Flow_gen));
	scnt = 0;
    }

    if (flowt[flowno].f == NULL) {
	fprintf(stderr, "Unable to allocate memory for reading timesteps\n");
	free(readdata);
	return 1;
    }

    offset = 0;
    first = 0;
/* replace j inside this loop with scnt */
    for (j = 0; j < nsteps; j++) {
	if (readdata[j] == 0) {
	    continue;
	}
/* start of timestep */

	if (t->type == FLOW) {
	    flowt[flowno].f[scnt].u = (float *) malloc(nmnp2 * sizeof(float));
	    flowt[flowno].f[scnt].v = (float *) malloc(nmnp2 * sizeof(float));
	    flowt[flowno].f[scnt].e = NULL;
	} else {
	    flowt[flowno].f[scnt].e = (float *) malloc(nmnp2 * sizeof(float));
	    flowt[flowno].f[scnt].u = NULL;
	    flowt[flowno].f[scnt].v = NULL;
	}
	for (ii = 0; ii < t->npts; ii++) {
	    for (i = 0; i < eh.nvrt * eh.ivs; i++) {
		if (t->type == FLOW) {
		    tmpu[i] = 0.0;
		} else {
		    tmpu[i] = -99.0;
		}
	    }
	    hh[0] = h[0][ii];
	    hh[1] = h[1][ii];
	    hh[2] = h[2][ii];
	    hh[3] = h[3][ii];

	    if (ElioGetXYData2(f64i, j, ellist[ii], &eh, hh, &ttime, &it, bind, sind, tmpu)) {
		fprintf(stderr, "GettransectNew(): Error from ElioGetNodeOld().\n");
		goto out;
	    }
	    reversearray(eh.nvrt * eh.ivs, tmpu);
	    /* output results */
	    if (t->type == FLOW) {
		for (i = 0; i < eh.nvrt; i++) {
		    flowt[flowno].f[scnt].u[i * t->npts + ii] = tmpu[2 * i];
		    flowt[flowno].f[scnt].v[i * t->npts + ii] = tmpu[2 * i + 1];
		}
	    } else {
		for (i = 0; i < eh.nvrt; i++) {
		    flowt[flowno].f[scnt].e[i * t->npts + ii] = tmpu[i];
		}
	    }
	}
	if (scnt == 0) {
	    if (t->type == FLOW) {
		tumin = tumax = flowt[flowno].f[scnt].u[0];
		tvmin = tvmax = flowt[flowno].f[scnt].v[0];
	    } else {
		tumin = tumax = flowt[flowno].f[scnt].e[0];
	    }
	}
	if (t->type == FLOW) {
	    minmax(nmnp2, flowt[flowno].f[scnt].u, &umin, &umax);
	    minmax(nmnp2, flowt[flowno].f[scnt].v, &vmin, &vmax);
	    flowt[flowno].f[scnt].umin = umin;
	    flowt[flowno].f[scnt].umax = umax;
	    flowt[flowno].f[scnt].vmin = vmin;
	    flowt[flowno].f[scnt].vmax = vmax;
	    tumin = (tumin > umin) ? umin : tumin;
	    tvmin = (tvmin > vmin) ? vmin : tvmin;
	    tumax = (tumax < umax) ? umax : tumax;
	    tvmax = (tvmax < vmax) ? vmax : tvmax;
	} else {
	    minmax(nmnp2, flowt[flowno].f[scnt].e, &umin, &umax);
	    flowt[flowno].f[scnt].emin = umin;
	    flowt[flowno].f[scnt].emax = umax;
	    tumin = (tumin > umin) ? umin : tumin;
	    tumax = (tumax < umax) ? umax : tumax;
	}
	flowt[flowno].f[scnt].time = ttime;
	scnt++;
    }
  out:;
    fclose(f64i);
    flowt[flowno].nsteps = scnt;
    if (t->type == FLOW) {
	flowt[flowno].umin = tumin;
	flowt[flowno].umax = tumax;
	flowt[flowno].vmin = tvmin;
	flowt[flowno].vmax = tvmax;
    } else {
	flowt[flowno].emin = tumin;
	flowt[flowno].emax = tumax;
    }
    flowt[flowno].type = 0;
    flowt[flowno].start = flowt[flowno].f[0].time;
    flowt[flowno].stop = flowt[flowno].f[flowt[flowno].nsteps - 1].time;
    flowt[flowno].active = ON;
/* final cleanup */
    free(ellist);
    free(nlist[0]);
    free(nlist[1]);
    free(nlist[2]);
    free(h[0]);
    free(h[1]);
    free(h[2]);
    free(readdata);
    return (0);
}

/* helper function to invert the order of the array so top is 0 and bottom is n */
void reversearray(int n, float *tmpu)
{
    int i;
    float t;
    for (i = 0; i < n / 2; i++) {
	t = tmpu[i];
	tmpu[i] = tmpu[n - i - 1];
	tmpu[n - i - 1] = t;
    }
}
