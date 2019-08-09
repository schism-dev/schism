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
 * set the parameters for isolines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isolwin.c,v 1.3 2006/12/27 01:02:59 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;

static void autoscale_isolines_proc(void);

char isoltitle[256] = "Isolines";

typedef struct _Isol_ui {
    int type;
    int nitem;
    int stype;
    char title[256];
    Widget top;
    Widget *n_item;
    Widget *precision_item;
    Widget start_item;
    Widget step_item;
    Widget minmax_item;
    Widget *line_item;
    Widget *layout_item;
    Widget *val_item;
    Widget *linepat_item;
    Widget *type_item;
    int ndefs;
    Widget *defs_item;
    Widget *defslabels_item;
    Widget *defs_frame;
    Widget *ssdefs_frame;
} Isol_ui;

static Widget isolines_frame;
static Widget *isol_n_item;
static Widget *isol_precision_item;
static Widget isol_start_item;
static Widget isol_step_item;
static Widget isol_minmax_item;
static Widget *isol_line_item;
static Widget *isol_layout_item;
static Widget *isol_val_item;
static Widget *isol_linepat_item;
static Widget *isol_type_item;
static Widget isol_defs_item[MAXISOLINES];
static Widget isol_defslabels_item[MAXISOLINES];
static Widget isol_defs_frame;
static Widget isol_ssdefs_frame;

static Widget legend_toggle;

char *xv_getstr(Widget w);
void update_isolines(void);
void create_cedit_popup(void);
Widget CreateTextItem3();

void linep_accept_proc(void);
void update_linep(void);
void create_linep_frame(void);

void savep_accept_proc(void);
void update_savep(void);
void create_savep_frame(void);

static Isolparms *iptr;

void isolleg_place_proc(void)
{
    curplaceitem = 0;
    set_action(0);
    set_action(PLACE_ISOLINES_LEGEND);
}

static void define_isolines_proc(void)
{
    char buf[80];
    int n, i;

    if (curip.nisol <= 0) {
	curip.nisol = 1;
    }
    if (curip.nisol > MAXISOLINES) {
	curip.nisol = MAXISOLINES;
    }
    curip.nisol = GetChoice(isol_n_item) + 1;
    if (curip.nisol <= 0) {
	curip.nisol = 1;
    }
    if (curip.nisol > MAXISOLINES) {
	curip.nisol = MAXISOLINES;
    }
    curip.isoltype = (int) GetChoice(isol_type_item);
    if (curip.isoltype == 0) {
	curip.cis[0] = atof(xv_getstr(isol_start_item));
	curip.cint = atof(xv_getstr(isol_step_item));
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = curip.cis[0] + i * curip.cint;
	}
    } else {
	if (curip.nisol > MAXISOLINES) {
	    curip.nisol = MAXISOLINES;
	}
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = atof(xv_getstr(isol_defs_item[i]));
	}
    }
    curip.type = (int) GetChoice(isol_line_item);
    curip.layout = (int) GetChoice(isol_layout_item) ? HORIZONTAL : VERTICAL;
/*
    curip.valflag = (int) GetChoice(isol_val_item) ? UP : DOWN;
*/
    curip.p.prec = (int) GetChoice(isol_precision_item);
    curip.lactive = XmToggleButtonGetState(legend_toggle) ? ON : OFF;
    *iptr = curip;
    update_isolines();
    update_linep();
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
	for (i = 0; i < MAXISOLINES; i++) {
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
	if (curip.nisol > MAXISOLINES) {
	    curip.nisol = MAXISOLINES;
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
/*
	SetChoice(isol_val_item, curip.valflag == UP ? 1 : 0);
*/
	update_linep();
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
    default_isolines(cg, &curip);
    update_isolines();
}

void create_isolines_popup(char *title, Isolparms * cd, int cm)
{
    int i;
    char s[80];
    Widget panel, bt, rc, rc1, fr, lab;
    XmString xms;
    extern long colors[];

    curip = *cd;
    iptr = cd;
    if (title != NULL) {
	strcpy(isoltitle, title);
    } else {
	strcpy(isoltitle, "Isolines");
    }
    if (isolines_frame) {
	update_isolines();
	XtVaSetValues(isolines_frame, XmNtitle, isoltitle, NULL);
	for (i = 0; i < MAXISOLINES; i++) {

	    if (cm == 0) {	/* bathymetry */
		XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + MAXISOLINES], NULL);
		sprintf(s, "#%2d (C%2d): ", i + 1, i + MAXISOLINES);
	    } else {
		XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
		sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	    }
	    xms = XmStringCreateLtoR(s, charset);
	    XtVaSetValues(isol_defslabels_item[i], XmNlabelString, xms, NULL);
	    XmStringFree(xms);
	}
	XtRaise(isolines_frame);
	return;
    }
    isolines_frame = XmCreateDialogShell(app_shell, "Isolines", NULL, 0);
    handle_close(isolines_frame);
    panel = XmCreateRowColumn(isolines_frame, "rc", NULL, 0);
    isol_minmax_item = CreateTextItem2(panel, 20, "Limits (status only) ");
    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    isol_line_item = CreatePanelChoice1(rc, "Draw isolines:", 4, "As lines", "Filled", "Both", 0, 0);
    isol_precision_item = CreatePanelChoice2(rc, "Precision:", 2, 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
    XtManageChild(rc);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    isol_n_item = CreatePanelChoice2(rc, "# of isolines:", 4, 17, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", 0, 0);

    isol_type_item = CreatePanelChoice1(rc, "Spacing by:", 3, "Start, step", "Specified", 0, 0);
    XtAddCallback(isol_type_item[2], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, 0);
    XtAddCallback(isol_type_item[3], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, (XtPointer) 1);
    XtManageChild(rc);

    isol_ssdefs_frame = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_ssdefs_frame, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    isol_start_item = CreateTextItem2(rc, 8, "Start:");
    isol_step_item = CreateTextItem2(rc, 8, "Step:");
    XtManageChild(rc);
    XtManageChild(isol_ssdefs_frame);

    isol_defs_frame = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_defs_frame, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = 0; i < MAXISOLINES / 2; i++) {
	if (cm == 0) {		/* bathymetry */
	    sprintf(s, "#%2d (C%2d): ", i + 1, i + MAXISOLINES);
	    isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + MAXISOLINES], NULL);
	} else {
	    sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	    isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
	}
    }
    XtManageChild(rc1);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = MAXISOLINES / 2; i < MAXISOLINES; i++) {
	sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	if (cm == 0) {		/* bathymetry */
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + MAXISOLINES], NULL);
	} else {
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
	}
    }
    XtManageChild(rc1);
    XtManageChild(rc);
    XtManageChild(isol_defs_frame);

    legend_toggle = XtVaCreateManagedWidget("Display legend", xmToggleButtonWidgetClass, panel, NULL);
    isol_layout_item = CreatePanelChoice1(panel, "Legend layout:", 3, "Vertical", "Horizontal", 0, 0);
    isol_val_item = CreatePanelChoice1(panel, "Legend values:", 3, "Ascending", "Descending", 0, 0);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) define_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Place legend", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolleg_place_proc, NULL);

    bt = XtVaCreateManagedWidget("Colors...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_cedit_popup, NULL);

    bt = XtVaCreateManagedWidget("Line props...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_linep_frame, NULL);

    bt = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_savep_frame, NULL);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolines_done_proc, NULL);
    XtManageChild(rc);
    XtManageChild(panel);
    XtRaise(isolines_frame);
    XtVaSetValues(isolines_frame, XmNtitle, isoltitle, NULL);
    update_isolines();
    update_linep();
    if (curip.nisol == 0) {
	autoscale_isolines_proc();
    }
}

static Widget *linep_color_item[MAXISOLINES];
static Widget *linep_linew_item[MAXISOLINES];
static Widget *linep_lines_item[MAXISOLINES];
static Widget linep_frame;

void create_linep_frame(void)
{
    Widget panel, wbut, fr, bb, sep, rc, rc2;
    int i;
    if (!linep_frame) {
	linep_frame = XmCreateDialogShell(app_shell, "Line properties", NULL, 0);
	panel = XmCreateRowColumn(linep_frame, "lineprc", NULL, 0);
	for (i = 0; i < MAXISOLINES; i++) {
	    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	    linep_color_item[i] = CreateColorChoice(rc, "Line color:", 1);
	    linep_linew_item[i] = CreatePanelChoice2(rc, "Width: ", 3, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	    linep_lines_item[i] = CreatePanelChoice1(rc, "Style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
	    XtManageChild(rc);
	}

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) linep_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, linep_frame);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    update_linep();
    XtRaise(linep_frame);
}

void update_linep(void)
{
    char buf[256];
    int i;
    if (linep_frame) {
	for (i = 0; i < MAXISOLINES; i++) {
	    SetChoice(linep_color_item[i], curip.color[i]);
	    SetChoice(linep_linew_item[i], curip.linew[i] - 1);
	    SetChoice(linep_lines_item[i], curip.lines[i] - 1);
	}
    }
}

void linep_accept_proc(void)
{
    char buf[256];
    int i;
    for (i = 0; i < MAXISOLINES; i++) {
	curip.color[i] = GetChoice(linep_color_item[i]);
	curip.linew[i] = GetChoice(linep_linew_item[i]) + 1;
	curip.lines[i] = GetChoice(linep_lines_item[i]) + 1;
    }
    *iptr = curip;
    XtUnmanageChild(linep_frame);
}

static Widget savep_frame;
static Widget savep_write;
static Widget savep_fname;
static Widget savep_toggle[MAXISOLINES];

void create_savep_frame(void)
{
    Widget panel, wbut, fr, bb, sep, rc, rc2;
    int i;
    char buf[128];
    if (!savep_frame) {
	savep_frame = XmCreateDialogShell(app_shell, "Save isolines", NULL, 0);
	panel = XmCreateRowColumn(savep_frame, "saveprc", NULL, 0);
	savep_write = XtVaCreateManagedWidget("Write isolines", xmToggleButtonWidgetClass, panel, NULL);
	savep_fname = CreateTextItem2(panel, 20, "Write to file:");
	for (i = 0; i < MAXISOLINES; i++) {
	    sprintf(buf, "Level #%1d", i);
	    savep_toggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, panel, NULL);
	}

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) savep_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, savep_frame);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    update_savep();
    XtRaise(savep_frame);
}

/*
 *  int writeflag;
 *  char wname[256];
 *  int writelevel[MAXISOLINES];
*/
void update_savep(void)
{
    char buf[256];
    int i;
    if (savep_frame) {
	XmToggleButtonSetState(savep_write, curip.writeflag, False);
	xv_setstr(savep_fname, curip.wname);
	for (i = 0; i < MAXISOLINES; i++) {
	    XmToggleButtonSetState(savep_toggle[i], curip.writelevel[i], False);
	}
    }
}

void savep_accept_proc(void)
{
    char buf[256];
    int i;
    strcpy(curip.wname, xv_getstr(savep_fname));
    curip.writeflag = XmToggleButtonGetState(savep_write);
    for (i = 0; i < MAXISOLINES; i++) {
	curip.writelevel[i] = XmToggleButtonGetState(savep_toggle[i]);
    }
    *iptr = curip;
    XtUnmanageChild(savep_frame);
}
