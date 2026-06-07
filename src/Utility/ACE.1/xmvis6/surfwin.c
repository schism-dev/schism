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
 * surfwin.c - Get top/bottom values from elcirc data files
 *
 */

#ifndef lint
static char RCSid[] = "$Id: surfwin.c,v 1.3 2004/03/04 04:22:26 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

/* read routine */
void create_surf_frame(Widget w, XtPointer clientd, XtPointer calld);

/*
 * create the file selection dialog for SELFE surface and bottom
 */

typedef struct {
    Widget top;
    Widget file;
    Widget *flow;
    Widget *surf;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget missing;
    Widget missingval;
    Widget append;
} surfUI;

static surfUI eui;

void accept_surf(Widget w, XtPointer clientd, XtPointer calld);

/*
 * set the current ADCIRC flow
 */
static void set_curadcirc(Widget w, int cd)
{
    if (cd != MAXADCIRC) {
	curadcirc = cd;
    }
}

void create_surf_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i;
    Widget wbut, panel, rc, rc2, wtmp;
    setistop();
    if (!eui.top) {
	eui.top = XmCreateDialogShell(app_shell, "SELFE Surface/Bottom", NULL, 0);
	panel = XmCreateRowColumn(eui.top, "panel", NULL, 0);

	eui.file = CreateTextItem2(panel, 10, "Data file: ");
	eui.flow = CreatePanelChoice1(panel, "Read to SELFE: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(eui.flow[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadcirc, (XtPointer) i);
	}
	eui.surf = CreatePanelChoice1(panel, "Surface or Bottom: ", 3, "Surface", "Bottom", 0, 0);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);

	eui.sample = XtVaCreateManagedWidget("Sample steps:", xmToggleButtonWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	eui.start = CreateTextItem2(rc, 4, "Start: ");
	xv_setstr(eui.start, "1");
	eui.stop = CreateTextItem2(rc, 4, "Stop: ");
	xv_setstr(eui.stop, "1");
	eui.skip = CreateTextItem2(rc, 4, "Skip: ");
	xv_setstr(eui.skip, "1");
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	eui.missing = XtVaCreateManagedWidget("Missing data:", xmToggleButtonWidgetClass, rc, NULL);
	eui.missingval = CreateTextItem2(rc, 4, "Value: ");
	xv_setstr(eui.missingval, "-99.0");
	XtManageChild(rc);

	eui.append = XtVaCreateManagedWidget("Append to flow:", xmToggleButtonWidgetClass, panel, NULL);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_surf, (XtPointer) & eui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) eui.top);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    XtRaise(eui.top);
}

void accept_surf(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256], file[2048];
    surfUI *q = (surfUI *) clientd;
    int flowno, gridno, surf, nlevels, filet = 0;
    int append = 0, sample, start = -1, stop = 1, skip = 1;
    int missing = 0;
    double missingval = 0.0;
    strcpy(file, (char *) xv_getstr(q->file));
    sample = XmToggleButtonGetState(q->sample);
    append = XmToggleButtonGetState(q->append);
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
	missingval = 0.0;
    }
    surf = GetChoice(q->surf);
    flowno = GetChoice(q->flow);
    if (!append && flowt[flowno].active == ON) {
        if (!yesno("Flow is active, kill it?", " ", " YES ", " NO ")) {
            return;
        }
    }
    set_wait_cursor(q->top);
    if (!ReadElcircSurf(flowno, file, surf, start, stop, skip, missing, missingval, append)) {
	sprintf(buf, "Error reading file %s", file);
	errwin(buf);
	return;
    } else {
	set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
	load_clock(ADCIRC, flowno);
	XtUnmanageChild(q->top);
	create_sadcirc_frame();
    }
    unset_wait_cursor(q->top);
}
