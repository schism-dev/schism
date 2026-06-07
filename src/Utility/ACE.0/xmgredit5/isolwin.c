/*
 * ACE/gredit - 2d finite element grid generation
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *                      All Rights Reserved.
 *
 */

/*
 *
 * set the parameters for isolines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isolwin.c,v 1.3 2006/07/31 23:41:47 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Isolparms *iptr;

static void autoscale_isolines_proc(void);
void create_saveisol_frame(void);

char isoltitle[256] = "Isolines";

static Widget isolines_frame;
static Widget isolines_panel;
static Widget *isol_n_item;
static Widget *isol_precision_item;
static Widget isol_start_item;
static Widget isol_step_item;
static Widget isol_minmax_item;
static Widget *isol_line_item;
static Widget *isol_layout_item;
static Widget *isol_linepat_item;
static Widget *isol_type_item;
static Widget isol_defs_item[16];
static Widget isol_defs_frame;
static Widget isol_ssdefs_frame;
static Widget isol_label_item;
static Widget isol_markstep_item;

static Widget legend_toggle;

char *xv_getstr();
void update_isolines(void);

void do_place_label(void)
{
    set_action(0);
    set_action(PLACE_ISOLINE_LABEL);
}

void isolleg_place_proc(void)
{
    set_action(0);
    set_action(PLACE_ISOLINES_LEGEND);
}

void set_isolleg_loc(double wx, double wy)
{
    curip.x = wx;
    curip.y = wy;
}

static void define_isolines_proc(void)
{
    char buf[80];
    int n, i;

    if (curip.nisol <= 0) {
	curip.nisol = 1;
    }
    if (curip.nisol > 16) {
	curip.nisol = 16;
    }
    curip.nisol = GetChoice(isol_n_item) + 1;
    if (curip.nisol <= 0) {
	curip.nisol = 1;
    }
    if (curip.nisol > 16) {
	curip.nisol = 16;
    }
    curip.isoltype = (int) GetChoice(isol_type_item);
    if (curip.isoltype == 0) {
	curip.cis[0] = atof(xv_getstr(isol_start_item));
	curip.cint = atof(xv_getstr(isol_step_item));
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = curip.cis[0] + i * curip.cint;
	}
    } else {
	if (curip.nisol > 16) {
	    curip.nisol = 16;
	}
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = atof(xv_getstr(isol_defs_item[i]));
	}
    }
    curip.type = (int) GetChoice(isol_line_item);
    curip.layout = (int) GetChoice(isol_layout_item) ? HORIZONTAL : VERTICAL;
    curip.p.prec = (int) GetChoice(isol_precision_item);
    curip.lactive = XmToggleButtonGetState(legend_toggle) ? ON : OFF;
    curip.marker = XmToggleButtonGetState(isol_label_item);
    curip.markstep = atoi(xv_getstr(isol_markstep_item));
    *iptr = curip;
    update_isolines();
    /* XtUnmanageChild(isolines_frame); */
}

static void set_isolines_type_proc(void)
{
    int i;
    char s[80];
    Arg a;

    curip.type = (int) GetChoice(isol_type_item);
    curip.p.prec = (int) GetChoice(isol_precision_item);
    if (curip.isoltype == 0) {
	sprintf(s, "%.*lf", curip.p.prec, curip.cmin);
	xv_setstr(isol_start_item, s);
	curip.cis[0] = curip.cmin;
	curip.cint = (curip.cmax - curip.cmin) / curip.nisol;
	sprintf(s, "%.*lf", curip.p.prec, curip.cint);
	xv_setstr(isol_step_item, s);
    } else {
	for (i = 0; i < 16; i++) {
	    sprintf(s, "%.*lf", curip.p.prec, curip.cis[i]);
	    xv_setstr(isol_defs_item[i], s);
	}
    }
}

void update_isolines(void)
{
    int i;
    char buf[256];

    if (isolines_frame) {
	sprintf(buf, "%lf %lf", curip.cmin, curip.cmax);
	xv_setstr(isol_minmax_item, buf);
	if (curip.nisol <= 0) {
	    curip.nisol = 1;
	}
	if (curip.nisol > 16) {
	    curip.nisol = 16;
	}
	SetChoice(isol_n_item, curip.nisol - 1);

	SetChoice(isol_type_item, curip.isoltype);
	if (curip.isoltype == 0) {
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cis[0]);
	    xv_setstr(isol_start_item, buf);
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cint);
	    xv_setstr(isol_step_item, buf);
	    for (i = 0; i < curip.nisol; i++) {
		sprintf(buf, "%.*lf", curip.p.prec, curip.cis[0] + i * curip.cint);
		xv_setstr(isol_defs_item[i], buf);
	    }
	    XtSetSensitive(isol_defs_frame, False);
	    XtSetSensitive(isol_ssdefs_frame, True);
	} else {
	    for (i = 0; i < curip.nisol; i++) {
		sprintf(buf, "%.*lf", curip.p.prec, curip.cis[i]);
		xv_setstr(isol_defs_item[i], buf);
	    }
	    XtSetSensitive(isol_defs_frame, True);
	    XtSetSensitive(isol_ssdefs_frame, False);
	}
	SetChoice(isol_line_item, curip.type);
	SetChoice(isol_precision_item, curip.p.prec);
	XmToggleButtonSetState(legend_toggle, curip.lactive == ON, False);
	SetChoice(isol_layout_item, curip.layout == HORIZONTAL ? 1 : 0);
	XmToggleButtonSetState(isol_label_item, curip.marker, False);
	sprintf(buf, "%d", curip.markstep);
	xv_setstr(isol_markstep_item, buf);
    }
}

static void set_curtype_proc(Widget w, int cd)
{
    curip.isoltype = cd;
    update_isolines();
}

static void isolines_done_proc(void)
{
    XtUnmanageChild(isolines_frame);
}

static void autoscale_isolines_proc(void)
{
    int i;
    if (isolines_frame) {
	curip.p.prec = (int) GetChoice(isol_precision_item);
    }
    default_isolines(&curip);
    update_isolines();
}

void create_isolines_popup(char *title, Isolparms * cd)
{
    int i;
    char s[80];
    Widget bt, rc, rc1, fr;
    XmString xms;

    curip = *cd;
    iptr = cd;
    if (title != NULL) {
	strcpy(isoltitle, title);
    } else {
	strcpy(isoltitle, "Isolines");
    }
    if (isolines_frame) {
	update_isolines();
	XtVaSetValues(isolines_frame,
		      XmNtitle, isoltitle,
		      NULL);
	XtRaise(isolines_frame);
	return;
    }
    isolines_frame = XmCreateDialogShell(app_shell, "Isolines", NULL, 0);
    handle_close(isolines_frame);
    isolines_panel = XmCreateRowColumn(isolines_frame, "rc", NULL, 0);
    isol_minmax_item = CreateTextItem2(isolines_panel,
				       20, "Limits (status only) ");
    rc = XmCreateRowColumn(isolines_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_line_item = CreatePanelChoice1(rc, "Draw:",
					4, "lines", "filled", "Both", (char *) NULL, 0);
    XtManageChild(rc);

    rc = XmCreateRowColumn(isolines_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_label_item = XtVaCreateManagedWidget("Auto place labels",
				       xmToggleButtonWidgetClass, rc, NULL);
    isol_markstep_item = CreateTextItem2(rc, 5, "Label every:");
    XtManageChild(rc);

    rc = XmCreateRowColumn(isolines_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_n_item = CreatePanelChoice2(rc, "# of isolines:",
				     4, 17,
				"1", "2", "3", "4", "5", "6", "7", "8", "9",
				   "10", "11", "12", "13", "14", "15", "16",
				     0, 0);

    isol_type_item = CreatePanelChoice1(rc, "Spacing by:",
					3,
					"start, step",
					"specified values",
					0, 0);
    XtAddCallback(isol_type_item[2], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, 0);
    XtAddCallback(isol_type_item[3], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, (XtPointer) 1);
    XtManageChild(rc);

    isol_ssdefs_frame = XmCreateFrame(isolines_panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_ssdefs_frame, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_start_item = CreateTextItem2(rc, 8, "Start:");
    isol_step_item = CreateTextItem2(rc, 8, "Step:");
    XtManageChild(rc);
    XtManageChild(isol_ssdefs_frame);

    isol_defs_frame = XmCreateFrame(isolines_panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_defs_frame, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = 0; i < 8; i++) {
	sprintf(s, "#%2d: ", i + 1);
	isol_defs_item[i] = CreateTextItem2(rc1, 12, s);
    }
    XtManageChild(rc1);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = 8; i < 16; i++) {
	sprintf(s, "#%2d: ", i + 1);
	isol_defs_item[i] = CreateTextItem2(rc1, 12, s);
    }
    XtManageChild(rc1);
    XtManageChild(rc);
    XtManageChild(isol_defs_frame);

    rc = XmCreateRowColumn(isolines_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    legend_toggle = XtVaCreateManagedWidget("Legend", xmToggleButtonWidgetClass, rc, NULL);
    isol_layout_item = CreatePanelChoice1(rc, "Layout:",
					  3,
					  "Vertical",
					  "Horizontal",
					  0, 0);
    isol_precision_item = CreatePanelChoice2(rc, "Precision:",
					     2, 11,
			"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0,
					     0);
    XtManageChild(rc);

    rc = XmCreateRowColumn(isolines_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) define_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Place legend", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolleg_place_proc, NULL);

    bt = XtVaCreateManagedWidget("Place label", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_place_label, NULL);

    bt = XtVaCreateManagedWidget("Save isoline...", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_saveisol_frame, NULL);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolines_done_proc, NULL);
    XtManageChild(rc);
    XtManageChild(isolines_panel);
    XtRaise(isolines_frame);
    XtVaSetValues(isolines_frame,
		  XmNtitle, isoltitle,
		  NULL);
    update_isolines();
    if (curip.nisol == 0) {
	autoscale_isolines_proc();
    }
}

static Widget saveisol_frame;
static Widget *saveisol_item;

void close_saveisol_proc(void)
{
    XtUnmanageChild(saveisol_frame);
}

void saveisol_proc(void)
{
    extern int save_isoline;
    extern char saveisol_fname[];
    Arg args;
    XmString list_item;
    char *s;
    int c;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(saveisol_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    if (!fexists(s)) {
	strcpy(saveisol_fname, s);
	XtFree(s);
    }
    save_isoline = GetChoice(saveisol_item);
    if (save_isoline == 17) {
	save_isoline = -1;
    }
}

void create_saveisol_frame(void)
{
    int i;
    Widget rc;
    Arg a;

    if (saveisol_frame) {
	XtRaise(saveisol_frame);
	return;
    }
    saveisol_frame = XmCreateFileSelectionDialog(app_shell, "Read data", NULL, 0);
    XtSetArg(a, XmNdialogTitle, XmStringCreateLtoR("Read data", charset));
    XtSetValues(saveisol_frame, &a, 1);
    XtAddCallback(saveisol_frame, XmNcancelCallback, (XtCallbackProc) close_saveisol_proc, 0);
    XtAddCallback(saveisol_frame, XmNokCallback, (XtCallbackProc) saveisol_proc, 0);

    rc = XmCreateRowColumn(saveisol_frame, "rc", NULL, 0);
    saveisol_item = CreatePanelChoice2(rc, "Save isoline #:",
				       4, 19,
				"1", "2", "3", "4", "5", "6", "7", "8", "9",
		     "10", "11", "12", "13", "14", "15", "16", "All", "OFF",
				       0, 0);
    XtManageChild(rc);
    XtManageChild(saveisol_frame);
    SetChoice(saveisol_item, 17);
}
