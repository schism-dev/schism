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
 * Track data display
 */

#ifndef lint
static char RCSid[] = "$Id: trackwin.c,v 1.2 2003/07/24 15:44:07 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

typedef struct {
    Widget top;
    Widget file;
    Widget display;
    Widget displayi;
    Widget displayd;
    Widget *track;
    Widget *color;
    Widget *sym;
    Widget symsize;
    Widget *linew;
    Widget *lines;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget missing;
    Widget missingval;
    Browser b;
} trackUI;

static trackUI ui;

static void accept_track(Widget w, XtPointer clientd, XtPointer calld);
static void set_curtrack(Widget w, XtPointer clientd, XtPointer calld);
static XtCallbackProc track_isolines_proc(Widget w, XtPointer clid, XtPointer calld);
int ReadTrack(int n, char *file, int start, int stop, int step, int missing, double missingval);
int GetIsolColor(double val, Isolparms p, int *cmap, int nmap);

void create_track_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i;
    Widget wbut, panel, rc, rc2, wtmp;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "Track", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "panel", NULL, 0);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	ui.b.file = CreateTextItem2(rc, 20, "Data file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & ui.b);
	XtManageChild(rc);

	ui.track = CreatePanelChoice1(panel, "Read to Track #: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < 10; i++) {
	    XtAddCallback(ui.track[i + 2], XmNactivateCallback, (XtCallbackProc) set_curtrack, (XtPointer) i);
	}
	ui.display = XtVaCreateManagedWidget("Display marker:", xmToggleButtonWidgetClass, panel, NULL);
	ui.displayi = XtVaCreateManagedWidget("Display iso-marker:", xmToggleButtonWidgetClass, panel, NULL);
	ui.displayd = XtVaCreateManagedWidget("Display values:", xmToggleButtonWidgetClass, panel, NULL);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	ui.color = CreateColorChoice(panel, "Color:", 1);
	ui.linew = CreatePanelChoice1(panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
	ui.lines = CreatePanelChoice1(panel, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);

	ui.sym = CreatePanelChoice1(panel, "Symbol: ", 13, "None", "Dot", "Circle", "Square", "Diamond", "Triangle up", "Triangle left", "Triangle down", "Triangle right", "Plus", "X", "Star", 0, 0);
	ui.symsize = XtVaCreateManagedWidget("symsize", xmScaleWidgetClass, panel,
					     XmNwidth, 150, XmNminimum, 0, XmNmaximum, 400, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNtitleString, XmStringCreateLtoR("Marker size", charset), XmNorientation, XmHORIZONTAL, NULL);

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
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_track, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) track_isolines_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) ui.top);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    XtRaise(ui.top);
}

void accept_track(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256], file[2048];
    trackUI *q = (trackUI *) clientd;
    int sample, start = -1, stop = 1, skip = 1;
    int tno;
    int missing;
    double missingval;
    strcpy(file, (char *) xv_getstr(q->b.file));
    sample = XmToggleButtonGetState(q->sample);
    if (sample) {
	strcpy(buf, (char *) xv_getstr(q->start));
	start = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->stop));
	stop = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->skip));
	skip = atoi(buf);
    } else {
	start = -1;
	skip = 1;
    }
    missing = XmToggleButtonGetState(q->missing);
    if (missing) {
	strcpy(buf, (char *) xv_getstr(q->missingval));
	missingval = atof(buf);
    } else {
	missing = 0;
    }
    tno = GetChoice(q->track);
    g[cg].track[tno].display = XmToggleButtonGetState(q->display) ? ON : OFF;
    g[cg].track[tno].display_data = XmToggleButtonGetState(q->displayd) ? ON : OFF;
    g[cg].track[tno].display_isolines = XmToggleButtonGetState(q->displayi) ? ON : OFF;
    g[cg].track[tno].p.color = GetChoice(q->color);
    g[cg].track[tno].p.symbol = GetChoice(q->sym);
    g[cg].track[tno].symsize = 1.0;
    set_wait_cursor(q->top);
    if (ReadTrack(tno, file, start, stop, skip, missing, missingval)) {
	return;
    }
    set_clock(0, track[tno].start, track[tno].stop, track[tno].step, track[tno].nsteps);
    load_clock(TRACK, tno);

    unset_wait_cursor(q->top);
}

static XtCallbackProc track_isolines_proc(Widget w, XtPointer clid, XtPointer calld)
{
    trackUI *ui = (trackUI *) clid;
    int n = GetChoice(ui->track);
    g[cg].salip.cmin = track[n].dmin;
    g[cg].salip.cmax = track[n].dmax;
    create_isolines_popup("Isolines of salinity", &g[cg].salip, 1);
}

static void set_curtrack(Widget w, XtPointer clientd, XtPointer calld)
{
    int n = (int) calld;
}

int ReadTrack(int n, char *file, int start, int stop, int step, int missing, double missingval)
{
    FILE *fp;
    int nsteps;
    int cnt, scnt = 0, i, j, k;
    double dmin, dmax;
    double *x, *y, *d, *t;
    char buf[2048];
    track[n].active = OFF;
    fp = fopen(file, "rb");
    if (fp == NULL) {
	errwin("Unable to open track file\n");
	return 1;
    }
    if (fgets(buf, 255, fp) == NULL) {
	errwin("Unable to open track file\n");
	return 1;
    }
    if (fgets(buf, 255, fp) == NULL) {
	errwin("Unable to open track file\n");
	return 1;
    }
    sscanf(buf, "%d", &nsteps);
    x = (double *) malloc(nsteps * sizeof(double));
    y = (double *) malloc(nsteps * sizeof(double));
    t = (double *) malloc(nsteps * sizeof(double));
    d = (double *) malloc(nsteps * sizeof(double));
    for (i = 0; i < nsteps; i++) {
	if (fgets(buf, 255, fp) == NULL) {
	    errwin("Error reading track file\n");
	    fclose(fp);
	    free(x);
	    free(y);
	    free(t);
	    free(d);
	    return 1;
	}
	sscanf(buf, "%lf %lf %lf %lf", &t[i], &x[i], &y[i], &d[i]);
	if (i == 0) {
	    dmin = dmax = d[i];
	} else {
	    dmin = dmin > d[i] ? d[i] : dmin;
	    dmax = dmax < d[i] ? d[i] : dmax;
	}
    }
    fclose(fp);
    track[n].active = ON;
    track[n].x = x;
    track[n].y = y;
    track[n].time = t;
    track[n].d = d;
    track[n].dmin = dmin;
    track[n].dmax = dmax;
    track[n].nsteps = nsteps;
    track[n].start = track[n].time[0];
    track[n].stop = track[n].time[nsteps - 1];
    return 0;
}

void drawtrack(int gno, int n, int step)
{
    int ctmp, color, sym;
    double symsize;
    char buf[256];
    if (step >= track[n].nsteps || g[gno].track[n].display == OFF || track[n].active != ON) {
	return;
    }
    ctmp = color = g[gno].track[n].p.color;
    sym = g[gno].track[n].p.symbol;
    symsize = g[gno].track[n].symsize;
    setcolor(color);
    if (g[gno].track[n].display_isolines == ON) {
	color = GetIsolColor(track[n].d[step], g[gno].salip, mapisolconc, g[gno].salip.nisol);
	if (color == 0) {
	    color = 1;
	}
    }
    setcolor(color);
    drawpolysym(&track[n].x[step], &track[n].y[step], 1, sym, 0, 1, symsize);
    setcolor(ctmp);
    drawpolysym(&track[n].x[step], &track[n].y[step], 1, sym, 0, 0, symsize);
    if (g[gno].track[n].display_data == ON) {
	setcolor(1);
	sprintf(buf, "  %.2lf", track[n].d[step]);
	writestr(track[n].x[step], track[n].y[step], 0, 0, buf);
    }
}

int GetIsolColor(double val, Isolparms p, int *cmap, int nmap)
{
    int i, k, n1, n2, n3;
    double cc, ccp1;
    for (k = 0; k < nmap - 1; k++) {
	cc = p.cis[k];
	ccp1 = p.cis[k + 1];
	if (val > cc && val <= ccp1) {
	    return cmap[k];
	}
    }
    return 0;
}
