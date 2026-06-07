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
 * newtrans.c - transects from elcirc 3d files.
 *
 */

#ifndef lint
static char RCSid[] = "$Id: newtrans.c,v 1.21 2005/06/17 15:23:27 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "elio.h"

extern Widget app_shell;
extern XmStringCharSet charset;

typedef struct {
    Widget top;
    Widget *trans;
    Widget *flow;
    Widget *graph;
    Widget pick;
    Widget *type;
    Widget *transtype;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget append;
    Browser b1;
    Browser b2;
    Browser b3;
} TransectUI;

void NewTransFrame(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc NewTransAccept(Widget w, XtPointer clid, XtPointer calld);

typedef struct {
    Widget top;
    Widget *trans;
    Widget display;
    Widget displayg;
    Widget *displayv;
} TransDisplayUI;

void NewTransDisplayFrame(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc NewTransDisplayAccept(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc TransIsolinesAccept(Widget w, XtPointer clid, XtPointer calld);

/*
 * find out about data at a node at this time step
 *   
 */
static void query_adcirc(void)
{
    extern Widget locate_point_message;
    set_action(0);
    set_action(QUERY_TRANSECT);
    create_points_frame();
    SetLabel(locate_point_message, "TRANSECT:");
}

XtCallbackProc NewTransDisplayAccept(Widget w, XtPointer clid, XtPointer calld)
{
    int disp, trans;
    TransDisplayUI *a = (TransDisplayUI *) clid;
    disp = GetChoice(a->displayv);
    trans = GetChoice(a->trans);
    switch (disp) {
	case 0:
	break;
	case 1:
	break;
	case 2:
	break;
	case 3:
	g[cg].trans[trans].display_mag = ON;
	break;
    }
    if (XmToggleButtonGetState(a->display)) {
	g[cg].trans[trans].display = ON;
    } else {
	g[cg].trans[trans].display = OFF;
    }
    if (XmToggleButtonGetState(a->displayg)) {
	g[cg].trans[trans].display_grid = ON;
    } else {
	g[cg].trans[trans].display_grid = OFF;
    }
}

XtCallbackProc TransIsolinesAccept(Widget w, XtPointer clid, XtPointer calld)
{
    TransDisplayUI *a = (TransDisplayUI *) clid;
    int n = GetChoice(a->trans);
    if (trans[n].active == ON) {
	g[cg].trans[n].ip.cmin = trans[n].emin;
	g[cg].trans[n].ip.cmax = trans[n].emax;
	create_isolines_popup("Transect Isolines", &g[cg].trans[n].ip, 1);
    } else {
    }
}

void NewTransDisplayFrame(Widget w, XtPointer clid, XtPointer calld)
{
    static TransDisplayUI ui;
    Arg al[8];
    Widget panel, wbut, fr, bb, sep, rc, rc1, rc2;
    int ac = 0, i;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "SELFE Transect Display", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "settransrc", NULL, 0);

	ui.trans = CreatePanelChoice1(panel, "Use transect: ", 6, "1", "2", "3", "4", "5", 0, 0);
	ui.display = XtVaCreateManagedWidget("Display transect", xmToggleButtonWidgetClass, panel, NULL);
	ui.displayg = XtVaCreateManagedWidget("Display transect grid", xmToggleButtonWidgetClass, panel, NULL);
	ui.displayv = CreatePanelChoice1(panel, "*Vectors: ", 5, "Not displayed", "Horizontal", "Along transect", "Magnitudes", 0, 0);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) NewTransDisplayAccept, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Isolines of scalars...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) TransIsolinesAccept, (XtPointer) & ui);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
}

XtCallbackProc NewTransAccept(Widget w, XtPointer clid, XtPointer calld)
{
    TransectUI *a = (TransectUI *) clid;
    int i;
    int datatype;
    int flowno;
    int gno;
    int sample, start = -1, stop = 1, skip = 1;
    int append = 0;
    int pick = 0;
    char tbuf[256];
    FILE *fp;

    curtrans = GetChoice(a->trans);
    if (trans[curtrans].active == ON) {
	if (!yesno("Transect is active, kill it?", " ", " YES ", " NO ")) {
	    return;
	} else {
	    free(trans[curtrans].g);
	    trans[curtrans].active = OFF;
	}
    }
/* here check for active grid and flow # */
    trans[curtrans].gno = GetChoice(a->graph);
    datatype = GetChoice(a->type);
    strcpy(trans[curtrans].salname, (char *) xv_getstr(a->b1.file));
    strcpy(trans[curtrans].uvname, (char *) xv_getstr(a->b1.file));
    strcpy(trans[curtrans].elevname, (char *) xv_getstr(a->b2.file));
    strcpy(trans[curtrans].transname, (char *) xv_getstr(a->b3.file));
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
    set_wait_cursor(a->top);
    switch (datatype) {
    case 0:
	trans[curtrans].type = SALINITY;
        if (ElioGetFileVersion(trans[curtrans].salname) == 4) {
	    if (ReadNewTrans(&trans[curtrans], append)) {
	        fprintf(stderr, "Error reading transect\n");
	        return;
            }
        } else if (ElioGetFileVersion(trans[curtrans].salname) == 5) {
	    if (ReadSigmaZTrans(&trans[curtrans], append)) {
	        fprintf(stderr, "Error reading transect\n");
	        return;
	    }
        } else {
	    return;
        }
	trans[curtrans].active = ON;
	set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
	load_clock(TRANSECT, curtrans);
	break;
    case 1:
	trans[curtrans].type = FLOW;
	if (ReadNewTrans(&trans[curtrans], append)) {
	    fprintf(stderr, "Error reading transect\n");
	    return;
	} else {
	    trans[curtrans].active = ON;
	    set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
	    load_clock(TRANSECT, curtrans);
	}
	break;
    }
    unset_wait_cursor(a->top);
    XtUnmanageChild(a->top);
    NewTransDisplayFrame(NULL, NULL, NULL);
}

void NewTransFrame(Widget w, XtPointer clid, XtPointer calld)
{
    static TransectUI ui;
    Arg al[8];
    Widget panel, wbut, fr, bb, sep, rc, rc1, rc2;
    int ac = 0, i;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "SELFE Transect setup", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "settransrc", NULL, 0);

	ui.trans = CreatePanelChoice1(panel, "Use transect: ", 6, "1", "2", "3", "4", "5", 0, 0);
	ui.graph = CreatePanelChoice1(panel, "Display in graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	ui.type = CreatePanelChoice1(panel, "Quantity to read: ", 3, "3d scalars", "3d vectors", 0, 0);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b1.file = CreateTextItem2(rc, 20, "Data file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b1);
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b2.file = CreateTextItem2(rc, 20, "Elevation file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b2);
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b3.file = CreateTextItem2(rc, 20, "Transect file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b3);
	XtManageChild(rc);

	//ui.data = CreateTextItem2(panel, 10, "Data file: ");
	//ui.elev = CreateTextItem2(panel, 10, "Elevation file: ");
	//ui.transname = CreateTextItem2(panel, 10, "Read transect data file: ");
	ui.append = XtVaCreateManagedWidget("*Append to active data set", xmToggleButtonWidgetClass, panel, NULL);
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

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) NewTransAccept, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
}

int ReadNewTrans(Transect * tr, int append)
{
    char fname[1024];
    int first = 1;
    int scnt, start, stop, skip;
    int *readdata;
    int i, j, k, n, n1, n2, n3;
    ElcircHeader h1, h2;
    ElcircTimeStep t1, t2;
    FILE *fp1, *fp2, *fp3;
    double x, y, sum;
    Point *p;
    int elem, *ellist, nn[4];
    int ne = 0, np = 0;
    int *icon[3];
    double *xp, *yp, *dp, *up;
    double hh[4];
    double dd[4], d, e;
    char buf[2048];
    float t, data[400];
    int it, bind, sind;
    int binds[4], sinds[4];
    int ij, nmin, nmax;
    int cnt = 0;

    if ((fp1 = fopen(tr->elevname, "r")) == NULL) {
	errwin("Unable to open elevation file.");
	return 1;
    }
    if (ElioGetFileType(fp1) == 4) {
	fclose(fp1);
	return ReadNewSigmaTrans(tr, append);
    } else if (ElioGetFileType(fp1) == 5) {
	fclose(fp1);
	return ReadSigmaZTrans(tr, append);
    }
    if (tr->type == SCALAR) {
	if ((fp2 = fopen(tr->salname, "r")) == NULL) {
	    errwin("Unable to open scalar file.");
	    fclose(fp1);
	    return 1;
	}
	strcpy(fname, tr->salname);
    } else {
	if ((fp2 = fopen(tr->uvname, "r")) == NULL) {
	    errwin("Unable to open vector file.");
	    fclose(fp1);
	    return 1;
	}
	strcpy(fname, tr->uvname);
    }
    if ((fp3 = fopen(tr->transname, "r")) == NULL) {
	errwin("Unable to open transect file");
	fclose(fp1);
	fclose(fp2);
	return 1;
    }

    ElioGetHeader(tr->elevname, &h1);
    ElioGetHeader(fname, &h2);
    h2.start_time[19] = 0;
    /*printf("[%s]\n", h2.start_time); */
    settime_info_start(cg, h2.start_time);

    if (h2.ivs == 2) {
	tr->type = VECTOR;
    } else {
	tr->type = SCALAR;
    }

    n1 = ElioGetNStepsInFile(tr->elevname, &h1);
    n2 = ElioGetNStepsInFile(fname, &h2);

    if (n1 != n2) {
	errwin("Number of steps in data files is not the same in both data files.");
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	ElioFreeHeader(&h1);
	ElioFreeHeader(&h2);
	return 1;
    }

    /* Read transect points */
    fgets(buf, 255, fp3);
    fgets(buf, 255, fp3);
    sscanf(buf, "%d", &n3);
    if (n3 < 2) {
	sprintf(buf, "Too few points to compute transect: %d\n", n3);
	errwin(buf);
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	ElioFreeHeader(&h1);
	ElioFreeHeader(&h2);
        return 1;
    }
    p = (Point *) malloc(n3 * sizeof(Point));
    ellist = (int *) malloc(n3 * sizeof(int));
    /* for drawing the transect in a plan view */
    tr->npts = n3;
    tr->x = (double *) malloc(n3 * sizeof(double));
    tr->y = (double *) malloc(n3 * sizeof(double));
    sum = 0.0;

    for (i = 0; i < n3; i++) {
	if (fgets(buf, 255, fp3) == NULL) {
	    return (1);
	}
	sscanf(buf, "%*d %lf %lf", &x, &y);
	tr->x[i] = x;
	tr->y[i] = y;
	if ((elem = ElioFindElementXY(&h1, x, y)) < 0) {
	    sprintf(buf, "Unable to locate element for point %d: %lf, %lf\n", i + 1, x, y);
	    errwin(buf);
	} else {
	    p[cnt].x = x;
	    p[cnt].y = y;
	    if (cnt == 0) {
		sum = p[cnt].dist = 0.0;
	    } else {
		sum = sum + hypot(p[cnt].x - p[cnt - 1].x, p[cnt].y - p[cnt - 1].y);
	    }
	    p[cnt].dist = sum;
	    nn[0] = h1.icon[0][elem];
	    nn[1] = h1.icon[1][elem];
	    nn[2] = h1.icon[2][elem];
	    if (h1.etype[elem] == 4) {
		nn[3] = h1.icon[3][elem];
	    }
	    /* get the coefficients for interpolation */
	    ElioGetCoefficients(&h1, elem, p[cnt].x, p[cnt].y, hh);
	    p[cnt].elem = elem;
	    p[cnt].type = h1.etype[elem];
	    p[cnt].h[0] = hh[0];
	    p[cnt].h[1] = hh[1];
	    p[cnt].h[2] = hh[2];
	    p[cnt].nn[0] = nn[0];
	    p[cnt].nn[1] = nn[1];
	    p[cnt].nn[2] = nn[2];
	    dd[0] = h1.d[nn[0]];
	    dd[1] = h1.d[nn[1]];
	    dd[2] = h1.d[nn[2]];
	    if (h1.etype[elem] == 4) {
		p[cnt].h[3] = hh[3];
		p[cnt].nn[3] = nn[3];
		dd[3] = h1.d[nn[3]];
	    }
	    /* get depth, fixed for all time */
	    ElioEval(p[cnt].type, hh, dd, &p[cnt].d);
	    cnt++;
	}
    }
    n3 = cnt;			/* reset to actual number of good points */
    fclose(fp3);

    ne = np = 0;
    xp = (double *) malloc(n3 * h2.nvrt * sizeof(double));
    yp = (double *) malloc(n3 * h2.nvrt * sizeof(double));
    dp = (double *) malloc(n3 * h2.nvrt * sizeof(double));
    up = (double *) malloc(n3 * h2.nvrt * sizeof(double));
    icon[0] = (int *) malloc(2 * n3 * h2.nvrt * sizeof(int));
    icon[1] = (int *) malloc(2 * n3 * h2.nvrt * sizeof(int));
    icon[2] = (int *) malloc(2 * n3 * h2.nvrt * sizeof(int));

    start = tr->start;
    stop = tr->stop;
    skip = tr->skip;
    if (start < 0) {
	start = 0;
	stop = n1 - 1;
	skip = 1;
    } else {
	if (start >= n1) {
	}
	if (stop >= n1) {
	    stop = n1 - 1;
	}
    }
    readdata = (int *) malloc(n1 * sizeof(int));
    for (i = 0; i < n1; i++) {
	readdata[i] = 0;
    }
    scnt = 0;
    for (i = start; i <= stop; i += skip) {
	readdata[i] = 1;
	scnt++;
    }
    if (scnt == 0) {
	errwin("No complete time step available, cancelling read");
	free(readdata);
	return 1;
    }
    tr->nsteps = scnt;
    /*fprintf(stderr, "Start, stop, skip, scnt %d %d %d %d\n", tr->start, tr->stop, tr->skip, scnt); */

    ElioAllocateTimeStep(&h1, &t1);
    ElioAllocateTimeStep(&h2, &t2);

    if ((tr->g = (Grid *) malloc(scnt * sizeof(Grid))) == NULL) {
	errwin("Unable to allocate memory for grids.");
	return 1;
    }
    if ((tr->t = (double *) malloc(scnt * sizeof(double))) == NULL) {
	errwin("Unable to allocate memory for time.");
	return 1;
    }

    scnt = 0;
    for (i = 0; i < n1; i++) {
	if (readdata[i] == 0) {
	    continue;
	}
	np = ne = 0;
	if (ElioGetTimeStep(fp1, i, &h1, &t1)) {
	    goto out;
	}
	if (ElioGetTimeStep(fp2, i, &h2, &t2)) {
	    goto out;
	}
	/*
	   printf("Readdata %d: %d\n", i, readdata[i]);
	   printf("%d %f %d %f\n", t1.it, t1.t, t2.it, t2.t);
	 */
	tr->t[scnt] = t1.t;

	for (j = 0; j < n3; j++) {
	    ElioInterpTimeStep(&h2, p[j].elem, p[j].x, p[j].y, p[j].h, &t2, binds, sinds, data);
	    p[j].nnodes = 0;
	    if (p[j].elem == -1) {
		printf("Skipping %d\n", j);
		continue;
	    }

	    dd[0] = t1.d[p[j].nn[0]];
	    dd[1] = t1.d[p[j].nn[1]];
	    dd[2] = t1.d[p[j].nn[2]];
	    if (p[j].type == 4) {
		dd[3] = t1.d[p[j].nn[3]];
	    }
	    ElioEval(p[j].type, hh, dd, &p[j].e);
	    //printf("%d %lf %lf\n", j, -p[j].d, p[j].e);
	    ElioGetZPos(&h2, p[j].d, p[j].e, &p[j].bi, &p[j].si, p[j].zpos1);
	    //printf("%d %d %d\n", k, p[j].bi, p[j].si);
	    //for (k=p[j].bi;k<=p[j].si;k++) {
	    //  printf("%d %lf %lf %lf\n", k, p[j].dist, p[j].zpos1[k], (double) k);
	    //}
	    for (k = p[j].bi; k <= p[j].si; k++) {
		xp[np] = p[j].dist;
		yp[np] = p[j].zpos1[k];
		if (h2.ivs == 2) {
		    dp[np] = p[j].val[2 * k] = data[2 * k];
		    up[np] = p[j].val[2 * k + 1] = data[2 * k + 1];
		} else if (h2.ivs == 1) {
		    dp[np] = p[j].val[k] = data[k];
		} else {
		}
		p[j].nodes[k] = np;
		p[j].nnodes++;
		np++;
	    }
	    if (j > 0) {
		if (p[j].bi == p[j - 1].bi) {	// both equal
		    if (p[j].si == p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 1 bi[j] == bi[j-1] si[j] == si[j - 1]\n");
			for (k = p[j].bi; k < p[j].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
		    } else if (p[j].si > p[j - 1].si) {	// bi equal si j > si j-1
			if (debuglevel == 13)
			    printf("case 2 bi[j] == bi[j-1] si[j] > si[j - 1]\n");
			for (k = p[j].bi; k < p[j - 1].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j - 1].nodes[p[j - 1].si];
			icon[1][ne] = p[j].nodes[p[j].si - 1];
			icon[2][ne] = p[j].nodes[p[j].si];
			ne++;
		    } else if (p[j].si < p[j - 1].si) {	// bi equal si j < si j-1
			if (debuglevel == 13)
			    printf("case 3 bi[j] == bi[j-1] si[j] < si[j - 1]\n");
			for (k = p[j].bi; k < p[j].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j].nodes[p[j].si];
			icon[1][ne] = p[j - 1].nodes[p[j - 1].si];
			icon[2][ne] = p[j - 1].nodes[p[j - 1].si - 1];
			ne++;
		    }
		} else if (p[j].bi > p[j - 1].bi) {
		    if (p[j].si == p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 4 bi[j] > bi[j-1] si[j] == si[j - 1]\n");
			nmin = p[j - 1].bi;
			nmax = p[j].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j].nodes[p[j].bi];
			    icon[1][ne] = p[j - 1].nodes[ij + 1];
			    icon[2][ne] = p[j - 1].nodes[ij];
			    ne++;
			}
			for (k = p[j].bi; k < p[j].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
		    } else if (p[j].si > p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 5 bi[j] > bi[j-1] si[j] > si[j - 1]\n");
			nmin = p[j - 1].bi;
			nmax = p[j].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j].nodes[p[j].bi];
			    icon[1][ne] = p[j - 1].nodes[ij + 1];
			    icon[2][ne] = p[j - 1].nodes[ij];
			    ne++;
			}
			for (k = p[j].bi; k < p[j - 1].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j - 1].nodes[k];
			icon[1][ne] = p[j].nodes[k];
			icon[2][ne] = p[j].nodes[k + 1];
			ne++;
		    } else if (p[j].si < p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 6 bi[j] > bi[j-1] si[j] < si[j - 1]\n");
			nmin = p[j - 1].bi;
			nmax = p[j].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j].nodes[p[j].bi];
			    icon[1][ne] = p[j - 1].nodes[ij + 1];
			    icon[2][ne] = p[j - 1].nodes[ij];
			    ne++;
			}
			for (k = p[j].bi; k < p[j].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j - 1].nodes[k];
			icon[1][ne] = p[j].nodes[k];
			icon[2][ne] = p[j - 1].nodes[k + 1];
			ne++;
		    }
		} else if (p[j].bi < p[j - 1].bi) {
		    if (p[j].si == p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 7 bi[j] < bi[j-1] si[j] == si[j - 1]\n");
			nmin = p[j].bi;
			nmax = p[j - 1].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j - 1].nodes[p[j - 1].bi];
			    icon[1][ne] = p[j].nodes[ij];
			    icon[2][ne] = p[j].nodes[ij + 1];
			    ne++;
			}
			for (k = p[j - 1].bi; k < p[j - 1].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
		    } else if (p[j].si > p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 8 bi[j] < bi[j-1] si[j] > si[j - 1]\n");
			nmin = p[j].bi;
			nmax = p[j - 1].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j - 1].nodes[p[j - 1].bi];
			    icon[1][ne] = p[j].nodes[ij];
			    icon[2][ne] = p[j].nodes[ij + 1];
			    ne++;
			}
			for (k = p[j - 1].bi; k < p[j - 1].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j - 1].nodes[k];
			icon[1][ne] = p[j].nodes[k];
			icon[2][ne] = p[j].nodes[k + 1];
			ne++;
		    } else if (p[j].si < p[j - 1].si) {
			if (debuglevel == 13)
			    printf("case 9 bi[j] < bi[j-1] si[j] < si[j - 1]\n");
			nmin = p[j].bi;
			nmax = p[j - 1].bi;
			for (ij = nmin; ij < nmax; ij++) {
			    icon[0][ne] = p[j - 1].nodes[p[j - 1].bi];
			    icon[1][ne] = p[j].nodes[ij];
			    icon[2][ne] = p[j].nodes[ij + 1];
			    ne++;
			}
			for (k = p[j - 1].bi; k < p[j].si; k++) {
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k];
			    icon[2][ne] = p[j].nodes[k + 1];
			    ne++;
			    icon[0][ne] = p[j - 1].nodes[k];
			    icon[1][ne] = p[j].nodes[k + 1];
			    icon[2][ne] = p[j - 1].nodes[k + 1];
			    ne++;
			}
			icon[0][ne] = p[j - 1].nodes[k];
			icon[1][ne] = p[j].nodes[k];
			icon[2][ne] = p[j - 1].nodes[k + 1];
			ne++;
		    }
		}
	    }
	}			/* j */

	AllocateGrid(&(tr->g[scnt]), ne, np);
	for (k = 0; k < np; k++) {
	    tr->g[scnt].xord[k] = xp[k];
	    tr->g[scnt].yord[k] = yp[k];
	    tr->g[scnt].depth[k] = dp[k];
	    if (tr->type == VECTOR) {
		tr->g[scnt].u[k] = dp[k];
		tr->g[scnt].v[k] = up[k];
		tr->g[scnt].depth[k] = hypot(dp[k], up[k]);
	    }
	}

	for (k = 0; k < ne; k++) {
	    tr->g[scnt].icon[k].nl[0] = icon[0][k];
	    tr->g[scnt].icon[k].nl[1] = icon[1][k];
	    tr->g[scnt].icon[k].nl[2] = icon[2][k];
	    tr->g[scnt].icon[k].type = tr->g[scnt].icon[k].nn = tr->g[scnt].icon[k].ngeom = 3;
	}
        /* WriteGrid("test.gr3", &tr->g[scnt]); */

	dminmax(tr->g[scnt].nmnp, tr->g[scnt].xord, &tr->g[scnt].xmin, &tr->g[scnt].xmax);
	dminmax(tr->g[scnt].nmnp, tr->g[scnt].yord, &tr->g[scnt].ymin, &tr->g[scnt].ymax);
	dminmax(tr->g[scnt].nmnp, tr->g[scnt].depth, &tr->g[scnt].dmin, &tr->g[scnt].dmax);
	if (debuglevel == 13)
	    printf("%d %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf\n", scnt, tr->g[scnt].xmin, tr->g[scnt].xmax, tr->g[scnt].ymin, tr->g[scnt].ymax, tr->g[scnt].dmin, tr->g[scnt].dmax);

	if (scnt == 0) {
	    tr->emin = tr->g[scnt].dmin;
	    tr->emax = tr->g[scnt].dmax;
	} else {
	    tr->emin = (tr->g[scnt].dmin < tr->emin) ? tr->g[scnt].dmin : tr->emin;
	    tr->emax = (tr->g[scnt].dmax > tr->emax) ? tr->g[scnt].dmax : tr->emax;
	}
	scnt++;
    }				/* i */

  out:;
    fclose(fp1);
    fclose(fp2);
    tr->nsteps = scnt;
    tr->tstart = tr->t[0];
    tr->tstop = tr->t[scnt - 1];
    if (tr->nsteps > 1) {
	tr->tstep = tr->t[1] - tr->t[0];
    } else {
	tr->tstep = 0.0;
    }
    free(xp);
    free(yp);
    free(dp);
    free(up);
    free(icon[0]);
    free(icon[1]);
    free(icon[2]);
    free(readdata);
    ElioFreeHeader(&h1);
    ElioFreeHeader(&h2);
    ElioFreeTimeStep(&t1);
    ElioFreeTimeStep(&t2);
    return 0;
}

void draw_trans(int gno, int stepno, int tno, double curtime)
{
    int j;
    double *tmp;
    double s;
    if (stepno >= trans[tno].nsteps) {
	return;
    }
    if (trans[tno].type == SCALAR) {
	tmp = (double *) malloc(trans[tno].g[stepno].nmnp * sizeof(double));
	if (tmp != NULL) {
	    for (j = 0; j < trans[tno].g[stepno].nmnp; j++) {
		tmp[j] = trans[tno].g[stepno].depth[j];
	    }
	    do_grid_isol(&trans[tno].g[stepno], 0, tmp, g[gno].trans[tno].ip, mapisolconc, g[gno].trans[tno].ip.nisol);
	    free(tmp);
	}
    } else if (trans[tno].type == VECTOR) {
	if (g[gno].trans[tno].display_mag == ON) {
	    tmp = (double *) malloc(trans[tno].g[stepno].nmnp * sizeof(double));
	    if (tmp != NULL) {
		for (j = 0; j < trans[tno].g[stepno].nmnp; j++) {
		    tmp[j] = trans[tno].g[stepno].depth[j];
		}
		do_grid_isol(&trans[tno].g[stepno], 0, tmp, g[gno].trans[tno].ip, mapisolconc, g[gno].trans[tno].ip.nisol);
		free(tmp);
	    }
	} else {
            setcolor(g[gno].trans[tno].p.color);
	    s = g[gno].vl.unitfac * g[gno].vl.scale;
	    for (j = 0; j < trans[tno].g[stepno].nmnp; j++) {
		velplt(trans[tno].g[stepno].xord[j], trans[tno].g[stepno].yord[j], trans[tno].g[stepno].u[j], trans[tno].g[stepno].v[j], s);
	    }
            setcolor(1);
	}
    }
    if (g[gno].trans[tno].ip.lactive == ON) {
	dolegend(g[gno].trans[tno].ip, mapisolconc, g[gno].trans[tno].ip.nisol);
    }
    if (g[gno].trans[tno].display_grid == ON) {
	DrawGrid(gno, &(trans[tno].g[stepno]), 1.0);
    }
}
