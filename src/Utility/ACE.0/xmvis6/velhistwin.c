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
 * Read ADCIRC
 */

#ifndef lint
static char RCSid[] = "$Id: velhistwin.c,v 1.2 2003/07/24 15:44:07 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
Widget CreateTextItem(Widget parent, int x, int y, int len, char *s);

static Widget vhist_frame;
static Widget vhist_panel;

static Widget svhist_frame;
static Widget svhist_panel;
void create_svhist_frame(void);
void update_vhist_flow(void);

int curvhist = 0;

/*
 * Panel item declarations
 */
static Widget vhist_files_list_item;
static Widget vhist_start_item;
static Widget vhist_stop_item;
static Widget vhist_nsteps_item;
static Widget *vhist_flow1_item;
static Widget *vhist_flow2_item;
static Widget *vhist_color_item;
static Widget vhist_toggle_item;
static Widget vhist_sym_item;
static Widget vhist_symsize_item;

/*
 * Event and Notify proc declarations
 */

void vhist_done_proc();
void vhist_accept_proc(void);
void svhist_doelev_proc();
void svhist_start_proc();
void svhist_stop_proc();
void svhist_nsteps_proc();
void svhist_accept_proc(void);
void svhist_done_proc(void);
void vhistel_place_notify_proc();
void vhistel_clear_notify_proc();
void toggle_vhistflow();

/*
 * Create the vhist Frame and the vhist Panel
 */

extern Widget app_shell;

void create_vhist_frame(void)
{
    Arg wargs[8];
    Widget wbut, rc;
    setistop();
    if (!vhist_frame) {
	vhist_frame = XmCreateFileSelectionDialog(app_shell, "vhist_frame", NULL, 0);
	XtVaSetValues(vhist_frame, XmNtitle, "Time histories of velocities", XmNdirMask, XmStringCreate("*", charset), NULL);
	rc = XmCreateRowColumn(vhist_frame, "rc", NULL, 0);
	vhist_flow1_item = CreatePanelChoice1(rc, "Read to flow: ", 6, "1", "2", "3", "4", "5", 0, 0);
	XtManageChild(rc);
	XtAddCallback(vhist_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, vhist_frame);
	XtAddCallback(vhist_frame, XmNokCallback, (XtCallbackProc) vhist_accept_proc, NULL);
    }
    XtRaise(vhist_frame);
}

void vhist_accept_proc(void)
{
    Arg args;
    XmString list_item;
    FILE *fp;
    char *s, buf[256], errbuf[256];
    int itmp;
    int flowno;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(vhist_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(buf, s);
    XtFree(s);

    flowno = curvhist;
    if (flowh[flowno].active) {
	if (!yesno("Flow is active, kill it?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    fp = fopen(buf, "r");
    if (fp == NULL) {
	errwin("Couldn't open file");
	return;
    }
    flowno = GetChoice(vhist_flow1_item);
    set_wait_cursor(vhist_frame);
    if (!readflowh2d(flowno, 0, buf)) {
	sprintf(errbuf, "Error reading file %s", buf);
	errwin(errbuf);
    } else {
	set_clock(0, flowh[flowno].start, flowh[flowno].stop, flowh[flowno].step, flowh[flowno].nsteps);
	load_clock(HISTORY, FLOW, flowno);
	XtUnmanageChild(vhist_frame);
	create_svhist_frame();
    }
    unset_wait_cursor(vhist_frame);
}

void create_svhist_frame(void)
{
    Widget wbut, fr, bb, sep, rc;
    int i, got_one = 0;
    setistop();
    for (i = 0; i < MAXVELHIST; i++) {
	if (flowh[i].active == ON) {
	    got_one = 1;
	}
    }
    if (!got_one) {
	create_vhist_frame();
	return;
    }
    if (!svhist_frame) {

	svhist_frame = XmCreateDialogShell(app_shell, "Time histories of flows setup", NULL, 0);

	svhist_panel = XmCreateRowColumn(svhist_frame, "svhistrc", NULL, 0);
	vhist_flow2_item = CreatePanelChoice1(svhist_panel, "Apply to flow: ", 7, "1", "2", "3", "4", "5", "All", 0, 0);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, svhist_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, NULL);
	vhist_start_item = CreateTextItem2(bb, 15, "Start time:");
	vhist_stop_item = CreateTextItem2(bb, 15, "Time step:");
	vhist_nsteps_item = CreateTextItem2(bb, 15, "Number of time steps:");
	XtManageChild(bb);
	XtManageChild(fr);

	vhist_toggle_item = XtVaCreateManagedWidget("Display flow", xmToggleButtonWidgetClass, svhist_panel, NULL);

	vhist_color_item = CreateColorChoice(svhist_panel, "Color:", 1);

	vhist_sym_item = XtVaCreateManagedWidget("Display marker", xmToggleButtonWidgetClass, svhist_panel, NULL);
	vhist_symsize_item = XtVaCreateManagedWidget("markerscale", xmScaleWidgetClass, svhist_panel,
						     XmNwidth, 150, XmNminimum, 0, XmNmaximum, 400, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNtitleString, XmStringCreateLtoR("Marker size", charset), XmNorientation, XmHORIZONTAL, NULL);

	rc = XmCreateRowColumn(svhist_panel, "svhistrc2", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) svhist_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Files...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_vhist_frame, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) svhist_done_proc, NULL);
	XtManageChild(rc);
	XtManageChild(svhist_panel);
    }
    update_vhist_flow();
    XtRaise(svhist_frame);
}

void svhist_done_proc(void)
{
    XtUnmanageChild(svhist_frame);
}

void update_vhist_flow(void)
{
    char buf[256];
    int flowno = curvhist;
    int value;
    if (svhist_frame) {
	SetChoice(vhist_flow2_item, curvhist);
	XmToggleButtonSetState(vhist_toggle_item, g[cg].flowh[flowno].display == ON, False);
	SetChoice(vhist_color_item, g[cg].flowh[flowno].p.color);
	if (flowh[flowno].nsteps > 0) {
	    sprintf(buf, "%.2lf", flowh[flowno].start);
	    XmTextSetString(vhist_start_item, buf);
	    sprintf(buf, "%.2lf", flowh[flowno].stop);
	    XmTextSetString(vhist_stop_item, buf);
	    sprintf(buf, "%d", flowh[flowno].nsteps);
	    XmTextSetString(vhist_nsteps_item, buf);
	}
	XmToggleButtonSetState(vhist_sym_item, g[cg].flowh[flowno].circle, False);
	value = g[cg].flowh[flowno].cp.symsize;
	XtVaSetValues(vhist_symsize_item, XmNvalue, value, NULL);
    }
}

void svhist_accept_proc(void)
{
    char buf[256];
    int i, flowno, start, stop, value;
    double tmp;
    Arg a;
    curvhist = flowno = GetChoice(vhist_flow2_item);
    g[cg].flowh[flowno].display = XmToggleButtonGetState(vhist_toggle_item) ? ON : OFF;
    g[cg].flowh[flowno].circle = XmToggleButtonGetState(vhist_sym_item);
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(vhist_symsize_item, &a, 1);
    g[cg].flowh[flowno].cp.symsize = value;
    g[cg].flowh[flowno].p.color = GetChoice(vhist_color_item);
/*
    curvhist = flowno = GetChoice(vhist_flow2_item);
    if (flowno == MAXVELHIST) {
	start = 0;
	stop = MAXVELHIST - 1;
    } else {
	start = stop = flowno;
    }
    for (i = start; i <= stop; i++) {
	g[cg].flowh[i].display = XmToggleButtonGetState(vhist_toggle_item) ? ON : OFF;
	g[cg].flowh[i].p.color = GetChoice(vhist_color_item);
    }
*/
    update_display();
}
