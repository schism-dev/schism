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
 * ticks / tick labels / axes labels
 *
 */

#ifndef lint
static char RCSid[] = "$Id: tickwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

Widget ticks_frame = (Widget) 0;
Widget ticks_panel;

static Widget *editaxis;	/* which axis to edit */
static Widget *axis_applyto;	/* ovverride */
static Widget offx;		/* x offset of axis in viewport coords */
static Widget offy;		/* y offset of axis in viewport coords */
static Widget altmap;		/* alternate mapping for axis */
static Widget altmin;		/* alternate mapping for axis */
static Widget altmax;		/* alternate mapping for axis */
static Widget tonoff;		/* toggle display of axis ticks */
static Widget tlonoff;		/* toggle display of tick labels */
static Widget mtlonoff;		/* toggle display of tick labels */
static Widget axislabel;	/* axis label */
static Widget axislabelop;	/* axis label on opposite side */
static Widget *axislabellayout;	/* axis label layout (perp or parallel) */
static Widget *axislabelfont;	/* axis label font */
static Widget axislabelcharsize;	/* axis label charsize */
static Widget *axislabelcolor;	/* axis label color */
static Widget *axislabellinew;	/* axis label linew */
static Widget tmajor;		/* major tick spacing */
static Widget tminor;		/* minor tick spacing */
static Widget *tickop;		/* ticks opposite */
static Widget *ticklop;		/* tick labels opposite */
static Widget *ticklabel_applyto;	/* ovverride */
static Widget *tlform;		/* format for labels */
static Widget *tlprec;		/* precision for labels */
static Widget *tlfont;		/* tick label font */
static Widget tlcharsize;	/* tick label charsize */
static Widget *tlcolor;		/* tick label color */
static Widget *tllinew;		/* tick label color */
static Widget tlvgap;		/* */
static Widget tlhgap;		/* */
static Widget *tlskip;		/* tick marks to skip */
static Widget tltype;		/* tick label type (auto or specified) */
static Widget ttype;		/* tick mark type (auto or specified) */
static Widget *tlstarttype;	/* use graph min or starting value */
static Widget tlstart;		/* value to start tick labels */
static Widget *tlstoptype;	/* use graph max or stop value */
static Widget tlstop;		/* value to stop tick labels */
static Widget *tllayout;	/* tick labels perp or horizontal or use the *
				 * angle */
static Widget tlangle;		/* angle */
static Widget *tlstagger;	/* stagger */
static Widget *tlsign;		/* sign of tick label (normal, negate, *
				 * absolute) */
static Widget tlspec;		/* tick labels specified */
static Widget *tick_applyto;	/* override */
static Widget tnum;		/* number of ticks for autoscaling */
static Widget tgrid;		/* major ticks grid */
static Widget *tgridcol;
static Widget *tgridlinew;
static Widget *tgridlines;
static Widget tmgrid;		/* minor ticks grid */
static Widget *tmgridcol;
static Widget *tmgridlinew;
static Widget *tmgridlines;
static Widget tlen;		/* tick length */
static Widget tmlen;
static Widget *tinout;		/* ticks in out or both */
static Widget tspec;		/* tick marks specified */
static Widget baronoff;		/* axis bar */
static Widget *barcolor;
static Widget *barlinew;
static Widget *barlines;

static Widget specticks;	/* special ticks and tick labels */
static Widget specticklabels;
static Widget nspec;
static Widget specnum[MAX_TICK_LABELS];	/* label denoting which tick/label */
static Widget specloc[MAX_TICK_LABELS];
static Widget speclabel[MAX_TICK_LABELS];

/*
 * Event and Notify proc declarations
 */
static void ticks_Done_notify_proc(void);
static void ticks_define_notify_proc(void);
static void do_axis_proc(void);
static void do_axislabel_proc(void);
static void do_ticklabels_proc(void);
static void do_tickmarks_proc(void);
static void do_axisbar_proc(void);
static void do_special_proc(void);
void autoticks_proc(void);
void update_ticks_items(int gno);
static void update_axis_items(int gno);
static void update_axislabel_items(int gno);
static void update_ticklabel_items(int gno);
static void update_tickmark_items(int gno);
static void update_axisbar_items(int gno);
static void update_special_items(int gno);
static void load_special(int gno, int a);

void update_ticks(int gno)
{
    update_ticks_items(gno);
    update_axis_items(gno);
    update_axislabel_items(gno);
    update_ticklabel_items(gno);
    update_tickmark_items(gno);
    update_axisbar_items(gno);
    load_special(gno, curaxis);
    update_special_items(gno);
}

void update_ticks_items(int gno)
{
    tickmarks t;

    if (ticks_frame) {
	SetChoice(editaxis, curaxis);
	get_graph_tickmarks(gno, &t, curaxis);
	XmToggleButtonGadgetSetState(tlonoff, t.tl_flag == ON, False);
	XmToggleButtonGadgetSetState(mtlonoff, t.mtl_flag == ON, False);
	XmToggleButtonGadgetSetState(tonoff, t.t_flag == ON, False);
	XmToggleButtonGadgetSetState(baronoff, t.t_drawbar == ON, False);
	XmTextSetString(axislabel, t.label.s);

	if (islogx(gno) && (curaxis % 2 == 0)) {
	    t.tmajor = (int) t.tmajor;
	    if (t.tmajor == 0) {
		t.tmajor = 1;
	    }
	    sprintf(buf, "%.0g", t.tmajor);
	} else if (islogy(gno) && (curaxis % 2 == 1)) {
	    t.tmajor = (int) t.tmajor;
	    if (t.tmajor == 0) {
		t.tmajor = 1;
	    }
	    sprintf(buf, "%.0g", t.tmajor);
	} else if (t.tmajor > 0) {
	    sprintf(buf, "%.5g", t.tmajor);
	} else {
	    strcpy(buf, "UNDEFINED");
	}
	XmTextSetString(tmajor, buf);
	if (islogx(gno) && (curaxis % 2 == 0)) {
	    t.tminor = (int) t.tminor;
	    if (t.tminor < 0 || t.tminor > 5) {
		t.tminor = 0;
	    }
	    sprintf(buf, "%.0g", t.tminor);
	} else if (islogy(gno) && (curaxis % 2 == 1)) {
	    t.tminor = (int) t.tminor;
	    if (t.tminor < 0 || t.tminor > 5) {
		t.tminor = 0;
	    }
	    sprintf(buf, "%.0g", t.tminor);
	} else if (t.tminor > 0) {
	    sprintf(buf, "%.5g", t.tminor);
	} else {
	    strcpy(buf, "UNDEFINED");
	}
	XmTextSetString(tminor, buf);
    }
}

static void set_axis_proc(Widget w, int cd)
{
    void update_ticks(int gno);

    curaxis = cd;
    update_ticks(cg);
}

/*
 * Create the ticks popup
 */
void create_ticks_frame(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc;
    int x, y;
    int i;

    if (ticks_frame) {
	update_ticks_items(cg);
	XtManageChild(ticks_panel);
	XtManageChild(ticks_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    ticks_frame = XmCreateDialogShell(app_shell, "Axes", NULL, 0);
    XtVaSetValues(ticks_frame, XmNx, x, XmNy, y, NULL);
    ticks_panel = XmCreateRowColumn(ticks_frame, "ticks_rc", NULL, 0);

    editaxis = (Widget *) CreatePanelChoice1(ticks_panel, "Edit:", 7, "X axis", "Y axis", "Zero X axis", "Zero Y axis", "Alternate X axis", "Alternate Y axis", NULL, NULL);
    for (i = 0; i < 6; i++) {
	XtAddCallback(editaxis[2 + i], XmNactivateCallback, (XtCallbackProc) set_axis_proc, (XtPointer) i);
    }

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    axislabel = (Widget) CreateTextItem2(ticks_panel, 30, "Axis label:");
    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    tmajor = (Widget) CreateTextItem2(ticks_panel, 10, "Major tick spacing:");
    tminor = (Widget) CreateTextItem2(ticks_panel, 10, "Minor tick spacing:");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    tlonoff = XtVaCreateManagedWidget("Tick labels", xmToggleButtonGadgetClass, ticks_panel, NULL);
    mtlonoff = XtVaCreateManagedWidget("Minor tick labels", xmToggleButtonGadgetClass, ticks_panel, NULL);
    tonoff = XtVaCreateManagedWidget("Tick marks", xmToggleButtonGadgetClass, ticks_panel, NULL);
    baronoff = XtVaCreateManagedWidget("Axis bar", xmToggleButtonGadgetClass, ticks_panel, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    rc = XmCreateRowColumn(ticks_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Axis props...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_axis_proc, 0);
    wbut = XtVaCreateManagedWidget("Axis label...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_axislabel_proc, 0);
    wbut = XtVaCreateManagedWidget("Tick labels...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_ticklabels_proc, 0);
    XtManageChild(rc);

    rc = XmCreateRowColumn(ticks_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Tick marks...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_tickmarks_proc, 0);
    wbut = XtVaCreateManagedWidget("Axis bar...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_axisbar_proc, 0);
    wbut = XtVaCreateManagedWidget("User ticks/tick labels...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_special_proc, 0);
    XtManageChild(rc);
    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    axis_applyto = (Widget *) CreatePanelChoice1(ticks_panel, "Apply to:", 5, "Current axis", "All axes, current graph", "Current axis, all graphs", "All axes, all graphs", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticks_panel, NULL);
    rc = XmCreateRowColumn(ticks_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) ticks_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Auto Ticks", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) autoticks_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) ticks_Done_notify_proc, 0);
    XtManageChild(rc);

    update_ticks_items(cg);
    XtManageChild(ticks_panel);
    XtManageChild(ticks_frame);
}

static void ticks_Done_notify_proc(void)
{
    XtUnmanageChild(ticks_frame);
}

/*
 * define tick marks
 */
static void ticks_define_notify_proc(void)
{
    char val[80];
    int i, j;
    int applyto;
    extern double result;
    double x = (g[cg].w.xg2 - g[cg].w.xg1), y = (g[cg].w.yg2 - g[cg].w.yg1), a = g[cg].w.xg1, b = g[cg].w.yg1, c = g[cg].w.xg2, d = g[cg].w.yg2;
    int errpos;
    tickmarks t;

    get_graph_tickmarks(cg, &t, curaxis);

    applyto = GetChoice(axis_applyto);
    strcpy(val, (char *) xv_getstr(tmajor));
    fixupstr(val);
    scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
    if (errpos) {
	return;
    }
    t.tmajor = result;
    if (islogx(cg) && (curaxis % 2 == 0)) {
	t.tmajor = (int) t.tmajor;
    } else if (islogy(cg) && (curaxis % 2 == 1)) {
	t.tmajor = (int) t.tmajor;
    }
    strcpy(val, (char *) xv_getstr(tminor));
    fixupstr(val);
    scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
    if (errpos) {
	return;
    }
    t.tminor = result;
    if (islogx(cg) && (curaxis % 2 == 0)) {
	t.tminor = (int) t.tminor;
	if (t.tminor < 0 || t.tminor > 5) {
	    t.tminor = 0;
	}
    } else if (islogy(cg) && (curaxis % 2 == 1)) {
	t.tminor = (int) t.tminor;
	if (t.tminor < 0 || t.tminor > 5) {
	    t.tminor = 0;
	}
    }
    t.tl_flag = XmToggleButtonGadgetGetState(tlonoff) ? ON : OFF;
    t.mtl_flag = XmToggleButtonGadgetGetState(mtlonoff) ? ON : OFF;
    t.t_flag = XmToggleButtonGadgetGetState(tonoff) ? ON : OFF;
    t.t_drawbar = XmToggleButtonGadgetGetState(baronoff) ? ON : OFF;
    strcpy(t.label.s, (char *) xv_getstr(axislabel));

    switch (applyto) {
    case 0:			/* current axis */
	set_graph_tickmarks(cg, &t, curaxis);
	break;
    case 1:			/* all axes, current graph */
	for (i = 0; i < 6; i++) {
	    g[cg].t[i].tl_flag = t.tl_flag;
	    g[cg].t[i].t_flag = t.t_flag;
	    g[cg].t[i].t_drawbar = t.t_drawbar;
	    strcpy(g[cg].t[i].label.s, t.label.s);
	    g[cg].t[i].tmajor = t.tmajor;
	    g[cg].t[i].tminor = t.tminor;
	}
	break;
    case 2:			/* current axis, all graphs */
	for (i = 0; i < MAXGRAPH; i++) {
	    g[i].t[curaxis].tl_flag = t.tl_flag;
	    g[i].t[curaxis].t_flag = t.t_flag;
	    g[i].t[curaxis].t_drawbar = t.t_drawbar;
	    strcpy(g[i].t[curaxis].label.s, t.label.s);
	    g[i].t[curaxis].tmajor = t.tmajor;
	    g[i].t[curaxis].tminor = t.tminor;
	}
	break;
    case 3:			/* all axes, all graphs */
	for (i = 0; i < MAXGRAPH; i++) {
	    for (j = 0; j < 6; j++) {
		g[i].t[i].tl_flag = t.tl_flag;
		g[i].t[i].t_flag = t.t_flag;
		g[i].t[i].t_drawbar = t.t_drawbar;
		strcpy(g[i].t[i].label.s, t.label.s);
		g[i].t[i].tmajor = t.tmajor;
		g[i].t[i].tminor = t.tminor;
	    }
	}
	break;
    }
    drawgraph();
}

static props_Done_notify_proc(Widget w)
{
    XtUnmanageChild(XtParent(XtParent(w)));
}

static Widget axis_frame;
static Widget axis_panel;

static void accept_axis_proc(Widget w)
{
    tickmarks t;

    get_graph_tickmarks(cg, &t, curaxis);
    t.alt = XmToggleButtonGadgetGetState(altmap) ? ON : OFF;
    t.tmin = atof((char *) xv_getstr(altmin));
    t.tmax = atof((char *) xv_getstr(altmax));
    t.offsx = atof((char *) xv_getstr(offx));
    t.offsy = atof((char *) xv_getstr(offy));
    set_graph_tickmarks(cg, &t, curaxis);
}

static void update_axis_items(int gno)
{
    tickmarks t;

    if (axis_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	XmToggleButtonGadgetSetState(altmap, t.alt == ON, False);
	sprintf(buf, "%.5g", t.tmin);
	XmTextSetString(altmin, buf);
	sprintf(buf, "%.5g", t.tmax);
	XmTextSetString(altmax, buf);
	sprintf(buf, "%.5g", t.offsx);
	XmTextSetString(offx, buf);
	sprintf(buf, "%.5g", t.offsy);
	XmTextSetString(offy, buf);
    }
}

static void do_axis_proc(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc;
    int x, y;

    if (axis_frame) {
	update_axis_items(cg);
	XtManageChild(axis_panel);
	XtManageChild(axis_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    axis_frame = XmCreateDialogShell(app_shell, "Axis props", NULL, 0);
    XtVaSetValues(axis_frame, XmNx, x, XmNy, y, NULL);
    axis_panel = XmCreateRowColumn(axis_frame, "axis_rc", NULL, 0);

    altmap = XtVaCreateManagedWidget("Use alternate map", xmToggleButtonGadgetClass, axis_panel, NULL);
    altmin = (Widget) CreateTextItem2(axis_panel, 10, "Alternate min:");
    altmax = (Widget) CreateTextItem2(axis_panel, 10, "Alternate max:");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, axis_panel, NULL);

    wlabel = XtVaCreateManagedWidget("Axis offset (viewport coordinates):", xmLabelGadgetClass, axis_panel, NULL);
    offx = (Widget) CreateTextItem2(axis_panel, 10, "Left or bottom:");
    offy = (Widget) CreateTextItem2(axis_panel, 10, "Right or top:");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, axis_panel, NULL);
    rc = XmCreateRowColumn(axis_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_axis_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);
    XtManageChild(rc);

    update_axis_items(cg);
    XtManageChild(axis_panel);
    XtManageChild(axis_frame);
}

static Widget axislabel_frame;
static Widget axislabel_panel;

static void accept_axislabel_proc(Widget w)
{
    Arg a;
    tickmarks t;
    int iv;

    get_graph_tickmarks(cg, &t, curaxis);
    t.label_layout = GetChoice(axislabellayout) ? PERP : PARA;
    t.label.font = GetChoice(axislabelfont);
    t.label.color = GetChoice(axislabelcolor);
    t.label.linew = GetChoice(axislabellinew) + 1;
    XtSetArg(a, XmNvalue, &iv);
    XtGetValues(axislabelcharsize, &a, 1);
    t.label.charsize = iv / 100.0;
    set_graph_tickmarks(cg, &t, curaxis);
    drawgraph();
}

static void update_axislabel_items(int gno)
{
    Arg a;
    tickmarks t;
    int iv;

    if (axislabel_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	SetChoice(axislabellayout, t.label_layout == PERP ? 1 : 0);
	SetChoice(axislabelfont, t.label.font);
	SetChoice(axislabelcolor, t.label.color);
	SetChoice(axislabellinew, t.label.linew - 1);
	iv = (int) (100 * t.label.charsize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(axislabelcharsize, &a, 1);
    }
}

static void do_axislabel_proc(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc;
    int x, y;

    if (axislabel_frame) {
	update_axislabel_items(cg);
	XtManageChild(axislabel_panel);
	XtManageChild(axislabel_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    axislabel_frame = XmCreateDialogShell(app_shell, "Axis label", NULL, 0);
    XtVaSetValues(axislabel_frame, XmNx, x, XmNy, y, NULL);
    axislabel_panel = XmCreateRowColumn(axislabel_frame, "axislabel_rc", NULL, 0);

    axislabellayout = (Widget *) CreatePanelChoice1(axislabel_panel, "Axis layout:", 3, "Parallel to axis", "Perpendicular to axis", NULL, NULL);

    axislabelfont = CreatePanelChoice1(axislabel_panel, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
    axislabelcolor = CreatePanelChoice2(axislabel_panel, "Color:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);

    axislabellinew = CreatePanelChoice1(axislabel_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    wlabel = XtVaCreateManagedWidget("Size:", xmLabelGadgetClass, axislabel_panel, NULL);
    axislabelcharsize = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, axislabel_panel, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, axislabel_panel, NULL);
    rc = XmCreateRowColumn(axislabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_axislabel_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);
    XtManageChild(rc);

    update_axislabel_items(cg);
    XtManageChild(axislabel_panel);
    XtManageChild(axislabel_frame);
}

static Widget ticklabel_frame;
static Widget ticklabel_panel;

static void accept_ticklabel_proc(Widget w)
{
    Arg a;
    tickmarks t;
    int iv;
    int i, j, applyto, gstart, gstop, astart, astop;;
    applyto = GetChoice(ticklabel_applyto);
    switch (applyto) {
    case 0:
	gstart = gstop = cg;
	astart = astop = curaxis;
	break;
    case 1:
	gstart = gstop = cg;
	astart = 0;
	astop = 5;
	break;
    case 2:
	gstart = 0;
	gstop = MAXGRAPH - 1;
	astart = astop = curaxis;
	break;
    case 3:
	gstart = 0;
	gstop = MAXGRAPH - 1;
	astart = 0;
	astop = 5;
	break;
    }
    for (i = gstart; i <= gstop; i++) {
	for (j = astart; j <= astop; j++) {
	    get_graph_tickmarks(i, &t, j);
	    t.tl_font = GetChoice(tlfont);
	    t.tl_color = GetChoice(tlcolor);
	    t.tl_linew = GetChoice(tllinew) + 1;
	    t.tl_skip = GetChoice(tlskip);
	    t.tl_prec = GetChoice(tlprec);
	    t.tl_staggered = (int) GetChoice(tlstagger);
	    t.tl_starttype = (int) GetChoice(tlstarttype) == 0 ? AUTO : SPEC;
	    if (t.tl_starttype == SPEC) {
		t.tl_start = atof((char *) xv_getstr(tlstart));
	    }
	    t.tl_stoptype = (int) GetChoice(tlstoptype) == 0 ? AUTO : SPEC;
	    if (t.tl_stoptype == SPEC) {
		t.tl_stop = atof((char *) xv_getstr(tlstop));
	    }
	    t.tl_format = format_types[(int) GetChoice(tlform)];
	    switch (GetChoice(ticklop)) {
	    case 0:
		if (j % 2) {
		    t.tl_op = LEFT;
		} else {
		    t.tl_op = BOTTOM;
		}
		break;
	    case 1:
		if (j % 2) {
		    t.tl_op = RIGHT;
		} else {
		    t.tl_op = TOP;
		}
		break;
	    case 2:
		t.tl_op = BOTH;
		break;
	    }
	    switch ((int) GetChoice(tlsign)) {
	    case 0:
		t.tl_sign = NORMAL;
		break;
	    case 1:
		t.tl_sign = ABSOLUTE;
		break;
	    case 2:
		t.tl_sign = NEGATE;
		break;
	    }
	    switch ((int) GetChoice(tllayout)) {
	    case 0:
		t.tl_layout = HORIZONTAL;
		break;
	    case 1:
		t.tl_layout = VERTICAL;
		break;
	    case 2:
		t.tl_layout = SPEC;
		XtSetArg(a, XmNvalue, &iv);
		XtGetValues(tlangle, &a, 1);
		t.tl_angle = iv;
		break;
	    }
	    XtSetArg(a, XmNvalue, &iv);
	    XtGetValues(tlcharsize, &a, 1);
	    t.tl_charsize = iv / 100.0;
	    set_graph_tickmarks(i, &t, j);
	}
    }
    drawgraph();
}

static void update_ticklabel_items(int gno)
{
    Arg a;
    tickmarks t;
    int iv;

    if (ticklabel_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	SetChoice(tlfont, t.tl_font);
	SetChoice(tlcolor, t.tl_color);
	SetChoice(tllinew, t.tl_linew - 1);
	SetChoice(tlskip, t.tl_skip);
	SetChoice(tlstagger, t.tl_staggered);
	SetChoice(tlstarttype, t.tl_starttype == SPEC);
	if (t.tl_starttype == SPEC) {
	    sprintf(buf, "%lf", t.tl_start);
	    xv_setstr(tlstart, buf);
	    sprintf(buf, "%lf", t.tl_stop);
	    xv_setstr(tlstop, buf);
	}
	SetChoice(tlstoptype, t.tl_stoptype == SPEC);
	if (t.tl_stoptype == SPEC) {
	    sprintf(buf, "%lf", t.tl_stop);
	    xv_setstr(tlstop, buf);
	}
	iv = get_format_index(t.tl_format);
	SetChoice(tlform, iv);
	switch (t.tl_op) {
	case LEFT:
	    SetChoice(ticklop, 0);
	    break;
	case RIGHT:
	    SetChoice(ticklop, 1);
	    break;
	case BOTTOM:
	    SetChoice(ticklop, 0);
	    break;
	case TOP:
	    SetChoice(ticklop, 1);
	    break;
	case BOTH:
	    SetChoice(ticklop, 2);
	    break;
	}
	switch (t.tl_sign) {
	case NORMAL:
	    SetChoice(tlsign, 0);
	    break;
	case ABSOLUTE:
	    SetChoice(tlsign, 1);
	    break;
	case NEGATE:
	    SetChoice(tlsign, 2);
	    break;
	}
	SetChoice(tlprec, t.tl_prec);
	iv = (int) (100 * t.tl_charsize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(tlcharsize, &a, 1);
	switch (t.tl_layout) {
	case HORIZONTAL:
	    SetChoice(tllayout, 0);
	    break;
	case VERTICAL:
	    SetChoice(tllayout, 1);
	    break;
	case SPEC:
	    SetChoice(tllayout, 2);
	    break;
	}
	iv = (int) t.tl_angle % 360;
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(tlangle, &a, 1);
    }
}

static void do_ticklabels_proc(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc;
    int x, y;

    if (ticklabel_frame) {
	update_ticklabel_items(cg);
	XtManageChild(ticklabel_panel);
	XtManageChild(ticklabel_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    ticklabel_frame = XmCreateDialogShell(app_shell, "Tick label", NULL, 0);
    XtVaSetValues(ticklabel_frame, XmNx, x, XmNy, y, NULL);
    ticklabel_panel = XmCreateRowColumn(ticklabel_frame, "ticklabel_rc", NULL, 0);

    tlfont = CreatePanelChoice1(ticklabel_panel, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);

    tlcolor = CreatePanelChoice2(ticklabel_panel, "Color:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);

    tllinew = CreatePanelChoice1(ticklabel_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    rc = XmCreateRowColumn(ticklabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Char size:", xmLabelGadgetClass, rc, NULL);
    tlcharsize = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 0, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    tlform = CreatePanelChoice2(ticklabel_panel,
				"Format:", 4,
				27,
				"Decimal",
				"Exponential",
				"Power",
				"General",
				"DD-MM-YY",
				"MM-DD-YY",
				"MM-YY",
				"MM-DD",
				"Month-DD",
				"DD-Month",
				"Month (abrev.)",
				"Month",
				"Day of week (abrev.)",
				"Day of week",
				"Day of year",
				"HH:MM:SS.s", "MM-DD HH:MM:SS.s", "MM-DD-YY HH:MM:SS.s", "Degrees (lon)", "DD MM' (lon)", "DD MM' SS.s\" (lon)", "MM' SS.s\" (lon)", "Degrees (lat)", "DD MM' (lat)", "DD MM' SS.s\" (lat)", "MM' SS.s\" (lat)", 0, 0);

    tlstagger = CreatePanelChoice1(ticklabel_panel, "Stagger labels:", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    tlprec = CreatePanelChoice1(ticklabel_panel, "Precision:", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    tlskip = CreatePanelChoice1(ticklabel_panel, "Skip every:", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    rc = XmCreateRowColumn(ticklabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    tlstarttype = CreatePanelChoice1(rc, "Start labels at:", 3, "Graph min", "Specified:", 0, 0);
    tlstart = XtVaCreateManagedWidget("tlstart", xmTextWidgetClass, rc, XmNtraversalOn, True, XmNcolumns, 10, NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(ticklabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    tlstoptype = CreatePanelChoice1(rc, "Stop labels at:", 3, "Graph max", "Specified:", 0, 0);
    tlstop = XtVaCreateManagedWidget("tlstop", xmTextWidgetClass, rc, XmNtraversalOn, True, XmNcolumns, 10, NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(ticklabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    tllayout = (Widget *) CreatePanelChoice1(rc, "Layout:", 4, "Horizontal", "Vertical", "Specified (degrees):", NULL, NULL);
    tlangle = XtVaCreateManagedWidget("ticklangle", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 360, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    ticklop = CreatePanelChoice1(ticklabel_panel, "Draw tick labels on:", 4, "Normal side", "Opposite side", "Both", 0, 0);

    tlsign = CreatePanelChoice1(ticklabel_panel, "Sign of label:", 4, "As is", "Absolute value", "Negate", NULL, 0);

    ticklabel_applyto = CreatePanelChoice1(ticklabel_panel, "Apply to:", 4, "Current axis", "All axes, current graph", "Current axis, all graphs", "All axes, all graphs", NULL, 0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, ticklabel_panel, NULL);
    rc = XmCreateRowColumn(ticklabel_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_ticklabel_proc, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);

    XtManageChild(rc);

    update_ticklabel_items(cg);
    XtManageChild(ticklabel_panel);
    XtManageChild(ticklabel_frame);
}

static Widget tickmark_frame;
static Widget tickmark_panel;

static void update_tickmark_items(int gno)
{
    Arg a;
    tickmarks t;
    int iv;

    if (tickmark_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	switch (t.t_inout) {
	case IN:
	    SetChoice(tinout, 0);
	    break;
	case OUT:
	    SetChoice(tinout, 1);
	    break;
	case BOTH:
	    SetChoice(tinout, 2);
	    break;
	}
	switch (t.t_op) {
	case LEFT:
	    SetChoice(tickop, 0);
	    break;
	case RIGHT:
	    SetChoice(tickop, 1);
	    break;
	case BOTTOM:
	    SetChoice(tickop, 0);
	    break;
	case TOP:
	    SetChoice(tickop, 1);
	    break;
	case BOTH:
	    SetChoice(tickop, 2);
	    break;
	}
	SetChoice(tgridcol, t.t_color);
	SetChoice(tgridlinew, t.t_linew - 1);
	SetChoice(tgridlines, t.t_lines - 1);
	SetChoice(tmgridcol, t.t_mcolor);
	SetChoice(tmgridlinew, t.t_mlinew - 1);
	SetChoice(tmgridlines, t.t_mlines - 1);
	iv = (int) (100 * t.t_size);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(tlen, &a, 1);
	iv = (int) (100 * t.t_msize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(tmlen, &a, 1);
	XmToggleButtonGadgetSetState(tgrid, t.t_gridflag == ON, False);
	XmToggleButtonGadgetSetState(tmgrid, t.t_mgridflag == ON, False);
    }
}

static void accept_tickmark_proc(Widget w)
{
    Arg a;
    tickmarks t;
    int iv;
    int i, j, applyto, gstart, gstop, astart, astop;;
    applyto = GetChoice(tick_applyto);
    switch (applyto) {
    case 0:
	gstart = gstop = cg;
	astart = astop = curaxis;
	break;
    case 1:
	gstart = gstop = cg;
	astart = 0;
	astop = 5;
	break;
    case 2:
	gstart = 0;
	gstop = MAXGRAPH - 1;
	astart = astop = curaxis;
	break;
    case 3:
	gstart = 0;
	gstop = MAXGRAPH - 1;
	astart = 0;
	astop = 5;
	break;
    }
    for (i = gstart; i <= gstop; i++) {
	for (j = astart; j <= astop; j++) {
	    get_graph_tickmarks(i, &t, j);
	    switch ((int) GetChoice(tinout)) {
	    case 0:
		t.t_inout = IN;
		break;
	    case 1:
		t.t_inout = OUT;
		break;
	    case 2:
		t.t_inout = BOTH;
		break;
	    }
	    switch (GetChoice(tickop)) {
	    case 0:
		if (j % 2) {
		    t.t_op = LEFT;
		} else {
		    t.t_op = BOTTOM;
		}
		break;
	    case 1:
		if (j % 2) {
		    t.t_op = RIGHT;
		} else {
		    t.t_op = TOP;
		}
		break;
	    case 2:
		t.t_op = BOTH;
		break;
	    }
	    t.t_color = GetChoice(tgridcol);
	    t.t_linew = GetChoice(tgridlinew) + 1;
	    t.t_lines = GetChoice(tgridlines) + 1;
	    t.t_mcolor = GetChoice(tmgridcol);
	    t.t_mlinew = GetChoice(tmgridlinew) + 1;
	    t.t_mlines = GetChoice(tmgridlines) + 1;
	    XtSetArg(a, XmNvalue, &iv);
	    XtGetValues(tlen, &a, 1);
	    t.t_size = iv / 100.0;
	    XtSetArg(a, XmNvalue, &iv);
	    XtGetValues(tmlen, &a, 1);
	    t.t_msize = iv / 100.0;
	    t.t_gridflag = XmToggleButtonGadgetGetState(tgrid) ? ON : OFF;
	    t.t_mgridflag = XmToggleButtonGadgetGetState(tmgrid) ? ON : OFF;
	    set_graph_tickmarks(i, &t, j);
	}
    }
    drawgraph();
}

static void do_tickmarks_proc(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc, rc2, rc3, fr;
    int x, y;

    if (tickmark_frame) {
	update_tickmark_items(cg);
	XtManageChild(tickmark_panel);
	XtManageChild(tickmark_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    tickmark_frame = XmCreateDialogShell(app_shell, "Tick marks", NULL, 0);
    XtVaSetValues(tickmark_frame, XmNx, x, XmNy, y, NULL);
    tickmark_panel = XmCreateRowColumn(tickmark_frame, "tickmark_rc", NULL, 0);

    tinout = CreatePanelChoice1(tickmark_panel, "Tick marks pointing:", 4, "In", "Out", "Both", 0, 0);

    tickop = CreatePanelChoice1(tickmark_panel, "Draw tick marks on:", 4, "Normal side", "Opposite side", "Both sides", 0, 0);

    rc2 = XmCreateRowColumn(tickmark_panel, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);

/* major tick marks */
    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);

    tgrid = XtVaCreateManagedWidget("Major ticks grid lines", xmToggleButtonGadgetClass, rc, NULL);

    rc3 = XmCreateRowColumn(rc, "rc3", NULL, 0);
    wlabel = XtVaCreateManagedWidget("Major tick length:", xmLabelGadgetClass, rc3, NULL);
    tlen = XtVaCreateManagedWidget("ticklength", xmScaleWidgetClass, rc3, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc3);

    tgridcol = CreatePanelChoice2(rc, "Color:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);

    tgridlinew = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
    tgridlines = (Widget *) CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);

    tmgrid = XtVaCreateManagedWidget("Minor ticks grid lines", xmToggleButtonGadgetClass, rc, NULL);
    rc3 = XmCreateRowColumn(rc, "rc", NULL, 0);
    wlabel = XtVaCreateManagedWidget("Minor tick length:", xmLabelGadgetClass, rc3, NULL);
    tmlen = XtVaCreateManagedWidget("mticklength", xmScaleWidgetClass, rc3, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc3);

    tmgridcol = CreatePanelChoice2(rc, "Color:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);
    tmgridlinew = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
    tmgridlines = (Widget *) CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    tick_applyto = CreatePanelChoice1(tickmark_panel, "Apply to:", 4, "Current axis", "All axes, current graph", "Current axis, all graphs", "All axes, all graphs", NULL, 0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, tickmark_panel, NULL);
    rc = XmCreateRowColumn(tickmark_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_tickmark_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);

    XtManageChild(rc);

    update_tickmark_items(cg);
    XtManageChild(tickmark_panel);
    XtManageChild(tickmark_frame);
}

static Widget axisbar_frame;
static Widget axisbar_panel;

static void accept_axisbar_proc(Widget w)
{
    tickmarks t;

    get_graph_tickmarks(cg, &t, curaxis);
    t.t_drawbarcolor = GetChoice(barcolor);
    t.t_drawbarlinew = GetChoice(barlinew) + 1;
    t.t_drawbarlines = GetChoice(barlines) + 1;
    set_graph_tickmarks(cg, &t, curaxis);
}

static void update_axisbar_items(int gno)
{
    tickmarks t;

    if (axisbar_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	SetChoice(barcolor, t.t_drawbarcolor);
	SetChoice(barlinew, t.t_drawbarlinew - 1);
	SetChoice(barlines, t.t_drawbarlines - 1);
    }
}

static void do_axisbar_proc(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    int x, y;

    if (axisbar_frame) {
	update_axisbar_items(cg);
	XtManageChild(axisbar_panel);
	XtManageChild(axisbar_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    axisbar_frame = XmCreateDialogShell(app_shell, "Axis bar", NULL, 0);
    XtVaSetValues(axisbar_frame, XmNx, x, XmNy, y, NULL);
    axisbar_panel = XmCreateRowColumn(axisbar_frame, "axisbar_rc", NULL, 0);

    barcolor = CreatePanelChoice2(axisbar_panel, "Color:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);

    barlinew = CreatePanelChoice1(axisbar_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

    barlines = (Widget *) CreatePanelChoice1(axisbar_panel, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, axisbar_panel, NULL);
    rc = XmCreateRowColumn(axisbar_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_axisbar_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);

    update_axisbar_items(cg);
    XtManageChild(axisbar_panel);
    XtManageChild(axisbar_frame);
}

static Widget special_frame;
static Widget special_panel;

#define TPAGESIZE 5
#define NPAGES (MAX_TICK_LABELS/TPAGESIZE)
static int tcurpage = 0;

static void accept_special_proc(Widget w)
{
    Arg a;
    tickmarks t;
    int iv, i;

    get_graph_tickmarks(cg, &t, curaxis);
    t.t_type = XmToggleButtonGadgetGetState(specticks) ? SPEC : AUTO;
    t.tl_type = XmToggleButtonGadgetGetState(specticklabels) ? SPEC : AUTO;
    iv = atoi((char *) xv_getstr(nspec));
    if (iv > MAX_TICK_LABELS) {
	sprintf(buf, "Number of ticks/tick labels exceeds %d", MAX_TICK_LABELS);
	errwin(buf);
	return;
    }
    t.t_spec = iv;
    for (i = 0; i < MAX_TICK_LABELS; i++) {
	t.t_specloc[i] = atof((char *) xv_getstr(specloc[i]));
	strcpy(t.t_speclab[i].s, (char *) xv_getstr(speclabel[i]));
    }
    set_graph_tickmarks(cg, &t, curaxis);
    drawgraph();
}

static void update_special_items(int gno)
{
    Arg a;
    tickmarks t;
    int iv, i, itmp;

    if (special_frame) {
	get_graph_tickmarks(gno, &t, curaxis);
	XmToggleButtonGadgetSetState(specticks, t.t_type == SPEC, False);
	XmToggleButtonGadgetSetState(specticklabels, t.tl_type == SPEC, False);
    }
}

static void load_special(int gno, int a)
{
    int i;
    char buf[128];
    tickmarks t;

    if (special_frame) {
	get_graph_tickmarks(gno, &t, a);
	sprintf(buf, "%d", t.t_spec);
	xv_setstr(nspec, buf);
	for (i = 0; i < t.t_spec; i++) {
	    sprintf(buf, "%lf", t.t_specloc[i]);
	    xv_setstr(specloc[i], buf);
	    xv_setstr(speclabel[i], t.t_speclab[i].s);
	}
    }
}

static void page_special_notify_proc(void)
{
    update_special_items(cg);
}

static void do_special_proc(void)
{
    extern Widget app_shell;
    Widget wbut, wlabel, rc, rc2, rc3, sw, sb, fr;
    int i, x, y;
    char buf[10];

    if (special_frame) {
	update_special_items(cg);
	XtManageChild(special_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    special_frame = XmCreateDialogShell(app_shell, "Specified ticks/ticklabels", NULL, 0);
    XtVaSetValues(special_frame, XmNx, x, XmNy, y, NULL);
    special_panel = XmCreateRowColumn(special_frame, "special_rc", NULL, 0);

    specticks = XtVaCreateManagedWidget("Use special tick locations", xmToggleButtonGadgetClass, special_panel, NULL);
    specticklabels = XtVaCreateManagedWidget("Use special tick labels", xmToggleButtonGadgetClass, special_panel, NULL);

    nspec = (Widget) CreateTextItem2(special_panel, 10, "# of user ticks/labels to use:");

    wlabel = XtVaCreateManagedWidget("Tick location - Label:", xmLabelGadgetClass, special_panel, NULL);

    sw = XtVaCreateManagedWidget("sw", xmScrolledWindowWidgetClass, special_panel, XmNscrollingPolicy, XmAUTOMATIC, NULL);
/*
    fr = XmCreateFrame(sw, "fr", NULL, 0);
    XtVaSetValues(sw,
		  XmNworkWindow, fr,
		  NULL);
*/
    rc = XmCreateRowColumn(sw, "rc", NULL, 0);
    XtVaSetValues(sw, XmNworkWindow, rc, NULL);

    for (i = 0; i < MAX_TICK_LABELS; i++) {
	rc3 = XmCreateRowColumn(rc, "rc3", NULL, 0);
	XtVaSetValues(rc3, XmNorientation, XmHORIZONTAL, NULL);
	sprintf(buf, "%2d", i + 1);
	specnum[i] = XtVaCreateManagedWidget(buf, xmLabelGadgetClass, rc3, NULL);
	specloc[i] = XtVaCreateManagedWidget("tickmark", xmTextFieldWidgetClass, rc3, XmNcolumns, 10, NULL);
	speclabel[i] = XtVaCreateManagedWidget("ticklabel", xmTextFieldWidgetClass, rc3, XmNcolumns, 35, NULL);
	XtManageChild(rc3);
    }
    XtManageChild(rc);
/*
    XtManageChild(fr);
*/
    XtManageChild(sw);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, special_panel, NULL);
    rc = XmCreateRowColumn(special_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_special_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) props_Done_notify_proc, 0);
    XtManageChild(rc);

    load_special(cg, curaxis);
    update_special_items(cg);
    XtManageChild(special_panel);
    XtManageChild(special_frame);
}
