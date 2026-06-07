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
 * graphics set routines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: graphwin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget display_frame;
static Widget display_panel;

static Widget auto_frame;
static Widget auto_panel;

static Widget gwind_frame;
static Widget gwind_panel;
static Widget graphwind_item;
static Widget gwindscroll_item;
static Widget shrink_item;
static Widget reduce_item;
static Widget gl_tog;
static Widget beleltol_item;

#define NITEMS 12
#define NAUTOITEMS 5

Widget toggle_display_item[NITEMS];
Widget toggle_auto_item[NAUTOITEMS];
extern Widget main_display_item[];

char *display_items[] = {
    "Grid",
    "Isolines",
    "Boundary",
    "Grid filled",
    "Node numbers",
    "Element numbers",
    "Depths at nodes",
    "Grid",
    "Isolines",
    "Coastal outline",
    "Build points"
};

char *auto_items[] = {
    "Edit Grid",
    "Edit boundary",
    "Build points",
    "Background grid",
    "Coastal outline"
};

extern int draw_filled;
extern int draw_isolfilled;

static double scrollper = 0.15;

void gwind_done_proc(void)
{
    XtUnmanageChild(gwind_frame);
}

void do_gwind_set(void)
{
}

void gwindleft_proc(void)
{
    double dx = scrollper * (xg2 - xg1);

    xg1 = xg1 - dx;
    xg2 = xg2 - dx;
    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(0.0, 0.0, 1.0, 1.0);
    do_drawgrid();
}

void gwindright_proc(void)
{
    double dx = scrollper * (xg2 - xg1);
    xg1 = xg1 + dx;
    xg2 = xg2 + dx;

    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(0.0, 0.0, 1.0, 1.0);
    do_drawgrid();
}

void gwinddown_proc(void)
{
    double dy = scrollper * (yg2 - yg1);
    yg1 = yg1 - dy;
    yg2 = yg2 - dy;

    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(0.0, 0.0, 1.0, 1.0);
    do_drawgrid();
}

void gwindup_proc(void)
{
    double dy = scrollper * (yg2 - yg1);
    yg1 = yg1 + dy;
    yg2 = yg2 + dy;

    defineworld(xg1, yg1, xg2, yg2, 0, 0);
    viewport(0.0, 0.0, 1.0, 1.0);
    do_drawgrid();
}

double shexper = 0.05;

void gwindshrink_proc(void)
{
    double dx = shexper * (xg2 - xg1);
    double dy = shexper * (yg2 - yg1);

    xg1 = xg1 + dx;
    xg2 = xg2 - dx;
    yg1 = yg1 + dy;
    yg2 = yg2 - dy;
    set_up_world();
    do_drawgrid();
}

void gwindexpand_proc(void)
{
    double dx = shexper * (xg2 - xg1);
    double dy = shexper * (yg2 - yg1);

    xg1 = xg1 - dx;
    xg2 = xg2 + dx;
    yg1 = yg1 - dy;
    yg2 = yg2 + dy;
    set_up_world();
    do_drawgrid();
}

static void gwinddefault_proc(void)
{
    sxg1 = xg1;
    syg1 = yg1;
    sxg2 = xg2;
    syg2 = yg2;
}

static void gwindrestoredefault_proc(void)
{

    set_defaults();
    xg1 = sxg1;
    yg1 = syg1;
    xg2 = sxg2;
    yg2 = syg2;
    reset_world();
}

static void gwind_autoscale_proc(void)
{
    set_defaults();
    set_up_world();
    do_drawgrid();
}

static void set_graph_scale(void)
{
}

/*
 * set scroll amount
 */
static void scroll_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(gwindscroll_item, &a, 1);
    scrollper = value / 100.0;
}

/*
 * set the reduction factor for drawing elements
 */
static void reduce_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(reduce_item, &a, 1);
    redfact = 1.0 - value / 100.0;
}

/*
 * set the reduction/expansion factor for shrink and expand
 */
static void shrink_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(shrink_item, &a, 1);
    shexper = value / 500.0;
}

/*
 * enlarge window
 */
static void enlarge_wind_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(shrink_item, &a, 1);
    write_mode_str("Expanding window");
    wfact = value / 100.0;
    xg1 = sxg1 - wfact * dx;
    yg1 = syg1 - wfact * dy;
    xg2 = sxg1 + (1 + wfact) * dx;
    yg2 = syg1 + (1 + wfact) * dy;
    set_up_world();
}

static void gwind_apply_notify_proc(void)
{
    Arg a;
    int value;
    extern int tdevice;
    char buf[30];
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(reduce_item, &a, 1);
    redfact = 1.0 - value / 100.0;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(shrink_item, &a, 1);
    shexper = value / 500.0;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(gwindscroll_item, &a, 1);
    scrollper = value / 100.0;
    belel_tol = -atof(xv_getstr(beleltol_item));

/*
    tdevice = XmToggleButtonGetState(gl_tog) ? 7 : 0;
    save_images = XmToggleButtonGetState(save_button) ? 1 : 0;
*/
}

void create_gwind_popup(void)
{
    int i;
    Widget bt, rc, rc2;
    XmString str;
    char buf[256];
    sprintf(buf, "%.9lf", -belel_tol);
    if (gwind_frame) {
	xv_setstr(beleltol_item, buf);
	XtRaise(gwind_frame);
	return;
    }
    gwind_frame = XmCreateDialogShell(app_shell, "Settings", NULL, 0);
    gwind_panel = XmCreateRowColumn(gwind_frame, "rc", NULL, 0);

    rc = XmCreateRowColumn(gwind_panel, "rc", NULL, 0);

    gwindscroll_item = XtVaCreateManagedWidget("Scroll %", xmScaleWidgetClass, rc,
					       XmNwidth, 150,
					       XmNminimum, 0,
					       XmNmaximum, 100,
					  XmNvalue, (int) (scrollper * 100),
					       XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
		    XmNtitleString, XmStringCreateLtoR("Scroll %", charset),
					       XmNorientation, XmHORIZONTAL,
					       NULL);
    XtAddCallback(gwindscroll_item, XmNvalueChangedCallback, (XtCallbackProc) scroll_proc, NULL);

    reduce_item = XtVaCreateManagedWidget("reduce", xmScaleWidgetClass, rc,
					  XmNwidth, 150,
					  XmNminimum, 0,
					  XmNmaximum, 90,
					  XmNvalue, 0,
					  XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
	 XmNtitleString, XmStringCreateLtoR("Element reduction %", charset),
					  XmNorientation, XmHORIZONTAL,
					  NULL);
    XtAddCallback(reduce_item, XmNvalueChangedCallback, (XtCallbackProc) reduce_proc, NULL);

    shrink_item = XtVaCreateManagedWidget("shrink", xmScaleWidgetClass, rc,
					  XmNwidth, 150,
					  XmNminimum, 0,
					  XmNmaximum, 500,
					  XmNvalue, (int) (500 * shexper),
					  XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
	     XmNtitleString, XmStringCreateLtoR("Expand/shrink %", charset),
					  XmNorientation, XmHORIZONTAL,
					  NULL);
    XtAddCallback(shrink_item, XmNvalueChangedCallback, (XtCallbackProc) shrink_proc, NULL);

    beleltol_item = CreateTextItem2(rc, 10, "Tolerance for element location: ");

    rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwind_apply_notify_proc, 0);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, gwind_frame);
    XtManageChild(rc2);
    XtManageChild(rc);

    xv_setstr(beleltol_item, buf);

    XtManageChild(gwind_panel);
    XtManageChild(gwind_frame);
}

void update_display_items(void)
{
    int i;
    if (display_frame) {
	for (i = 0; i < ndisplay; i++) {
	    XmToggleButtonSetState(toggle_display_item[i], display_flags[i], False);
	}
    }
    for (i = 0; i < ndisplay; i++) {
	XmToggleButtonSetState(main_display_item[i], display_flags[i], False);
    }
}

void get_display_toggles(void)
{
    int i;
    for (i = 0; i < ndisplay; i++) {
	display_flags[i] = XmToggleButtonGetState(toggle_display_item[i]);
    }
    for (i = 0; i < ndisplay; i++) {
	XmToggleButtonSetState(main_display_item[i], display_flags[i], False);
    }
    do_drawgrid();
}

void display_done_proc(void)
{
    XtUnmanageChild(display_frame);
}

Widget *grid_color_item;

void create_display_popup(void)
{
    int i, ac;
    Widget bt, lab, rb, fr, rc, rc2, sep;
    Dimension ww, wh;
    Position x, y;
    Arg a[2];
    XmString str;
    if (display_frame) {
	update_display_items();
	XtRaise(display_frame);
	return;
    }
    display_frame = XmCreateDialogShell(app_shell, "Features", NULL, 0);
    display_panel = XmCreateRowColumn(display_frame, "rc", NULL, 0);
    rc = display_panel;
/*
    grid_color_item = CreateColorChoice(display_panel, "Grid color", 0);
*/
    rc = XmCreateRowColumn(display_panel, "rc", NULL, 0);
    XtVaCreateManagedWidget("Grid:", xmLabelGadgetClass, rc, XmNy, 35 * 0, NULL);
    fr = XmCreateFrame(rc, "frame_1", NULL, 0);
    ac = 0;
    XtSetArg(a[ac], XmNradioBehavior, False);
    ac++;
    rb = XmCreateRadioBox(fr, "radio_box_1", a, ac);
/*
display_flags[EDIT_GRID]
display_flags[EDIT_GRID_ISOLINES]
display_flags[EDIT_BOUNDARY]
display_flags[EDIT_GRID_FILLED]
display_flags[EDIT_GRID_NODE_NUMBERS]
display_flags[EDIT_GRID_ELEMENT_NUMBERS]
display_flags[EDIT_GRID_DEPTHS]
display_flags[BACKGROUND_GRID]
display_flags[BACKGROUND_GRID_ISOLINES]
display_flags[BACKGROUND_BOUNDARY]
display_flags[BUILD_POINTS]
*/
    toggle_display_item[EDIT_GRID] = XtVaCreateManagedWidget("Grid",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_ISOLINES] = XtVaCreateManagedWidget("Isolines",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_BOUNDARY] = XtVaCreateManagedWidget("Boundary",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_FILLED] = XtVaCreateManagedWidget("Fill grid",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_NODE_NUMBERS] = XtVaCreateManagedWidget("Node numbers",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_NODES] = XtVaCreateManagedWidget("Nodes (markers)",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_ELEMENT_NUMBERS] = XtVaCreateManagedWidget("Element numbers",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[EDIT_GRID_DEPTHS] = XtVaCreateManagedWidget("Node depths",
				       xmToggleButtonWidgetClass, rb, NULL);
    XtManageChild(fr);
    XtManageChild(rb);
    XtVaCreateManagedWidget("Background grid:", xmLabelGadgetClass, rc, XmNy, 35 * 6, NULL);
    fr = XmCreateFrame(rc, "frame_2", NULL, 0);
    ac = 0;
    XtSetArg(a[ac], XmNradioBehavior, False);
    ac++;
    rb = XmCreateRadioBox(fr, "radio_box_2", a, ac);

    toggle_display_item[BACKGROUND_GRID] = XtVaCreateManagedWidget("Grid",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[BACKGROUND_GRID_ISOLINES] = XtVaCreateManagedWidget("Isolines",
				       xmToggleButtonWidgetClass, rb, NULL);
    toggle_display_item[BACKGROUND_BOUNDARY] = XtVaCreateManagedWidget("Coastal outline",
				       xmToggleButtonWidgetClass, rb, NULL);
    XtManageChild(rb);
    XtManageChild(fr);

    XtVaCreateManagedWidget("Build points:", xmLabelGadgetClass, rc, XmNy, 35 * 10, NULL);
    fr = XmCreateFrame(rc, "frame_3", NULL, 0);
    ac = 0;
    XtSetArg(a[ac], XmNradioBehavior, False);
    ac++;
    rb = XmCreateRadioBox(fr, "radio_box_3", a, ac);
    toggle_display_item[BUILD_POINTS] = XtVaCreateManagedWidget("Build points",
				       xmToggleButtonWidgetClass, rb, NULL);
    XtManageChild(rb);
    XtManageChild(fr);

    fr = XmCreateFrame(rc, "frame_3", NULL, 0);
    ac = 0;
    XtSetArg(a[ac], XmNpacking, XmPACK_NONE);
    ac++;
    rc2 = XmCreateRowColumn(fr, "rc2", a, ac);

    sep = XmCreateSeparatorGadget(rc, "sep", NULL, 0);
    XtManageChild(sep);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) display_done_proc, NULL);

    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) get_display_toggles, 0);
    XtManageChild(fr);
    XtManageChild(rc2);

    XtManageChild(rc);

    XtManageChild(display_panel);
    XtManageChild(display_frame);
    update_display_items();
}

static int auto_edit_grid = 0;
static int auto_edit_grid_boundary = 0;
static int auto_back_grid = 0;
static int auto_back_boundary = 0;
static int auto_build_points = 0;

void do_autoscale(void)
{
    double xmin1, xmax1, ymin1, ymax1, dmin1, dmax1;
    double xmin2, xmax2, ymin2, ymax2, dmin2, dmax2;
    int doit = 0;
    int i, ib;
    xmin2 = 1e307;
    xmax2 = -1e307;
    ymin2 = 1e307;
    ymax2 = -1e307;
    dmin2 = 1e307;
    dmax2 = -1e307;
    if (auto_edit_grid) {
	getlimits_grid(curgrid, &xmin1, &xmax1, &ymin1, &ymax1, &dmin1, &dmax1);
	xmin2 = xmin2 < xmin1 ? xmin2 : xmin1;
	xmax2 = xmax2 > xmax1 ? xmax2 : xmax1;
	ymin2 = ymin2 < ymin1 ? ymin2 : ymin1;
	ymax2 = ymax2 > ymax1 ? ymax2 : ymax1;
	dmin2 = dmin2 < dmin1 ? dmin2 : dmin1;
	dmax2 = dmax2 > dmax1 ? dmax2 : dmax1;
	doit = 1;
    }
    if (auto_edit_grid_boundary) {
	for (i = 0; i < grid[curgrid].nbounds; i++) {
	    ib = grid[curgrid].boundaries[i];
	    getlimits_boundary(ib, &xmin1, &xmax1, &ymin1, &ymax1);
	    xmin2 = xmin2 < xmin1 ? xmin2 : xmin1;
	    xmax2 = xmax2 > xmax1 ? xmax2 : xmax1;
	    ymin2 = ymin2 < ymin1 ? ymin2 : ymin1;
	    ymax2 = ymax2 > ymax1 ? ymax2 : ymax1;
	    doit = 1;
	}
    }
    if (auto_back_grid) {
	getlimits_grid(MAXGRIDS, &xmin1, &xmax1, &ymin1, &ymax1, &dmin1, &dmax1);
	xmin2 = xmin2 < xmin1 ? xmin2 : xmin1;
	xmax2 = xmax2 > xmax1 ? xmax2 : xmax1;
	ymin2 = ymin2 < ymin1 ? ymin2 : ymin1;
	ymax2 = ymax2 > ymax1 ? ymax2 : ymax1;
	dmin2 = dmin2 < dmin1 ? dmin2 : dmin1;
	dmax2 = dmax2 > dmax1 ? dmax2 : dmax1;
	doit = 1;
    }
    if (auto_back_boundary) {
	for (i = 0; i < grid[MAXGRIDS].nbounds; i++) {
	    ib = grid[MAXGRIDS].boundaries[i];
	    getlimits_boundary(ib, &xmin1, &xmax1, &ymin1, &ymax1);
	    xmin2 = xmin2 < xmin1 ? xmin2 : xmin1;
	    xmax2 = xmax2 > xmax1 ? xmax2 : xmax1;
	    ymin2 = ymin2 < ymin1 ? ymin2 : ymin1;
	    ymax2 = ymax2 > ymax1 ? ymax2 : ymax1;
	    doit = 1;
	}
    }
    if (auto_build_points) {
	getlimits_build(curbuild, &xmin1, &xmax1, &ymin1, &ymax1, &dmin1, &dmax1);
	xmin2 = xmin2 < xmin1 ? xmin2 : xmin1;
	xmax2 = xmax2 > xmax1 ? xmax2 : xmax1;
	ymin2 = ymin2 < ymin1 ? ymin2 : ymin1;
	ymax2 = ymax2 > ymax1 ? ymax2 : ymax1;
	dmin2 = dmin2 < dmin1 ? dmin2 : dmin1;
	dmax2 = dmax2 > dmax1 ? dmax2 : dmax1;
	doit = 1;
    }
    if (doit) {
	set_default_scale(xmin2, xmax2, ymin2, ymax2);
	set_scale();
    } else {
    }
}

void set_autoscale(void)
{
    auto_edit_grid = (int) XmToggleButtonGetState(toggle_auto_item[0]);
    auto_edit_grid_boundary = (int) XmToggleButtonGetState(toggle_auto_item[1]);
    auto_build_points = (int) XmToggleButtonGetState(toggle_auto_item[2]);
    auto_back_grid = (int) XmToggleButtonGetState(toggle_auto_item[3]);
    auto_back_boundary = (int) XmToggleButtonGetState(toggle_auto_item[4]);
    do_autoscale();
    do_drawgrid();
}

void set_auto_toggles(void)
{
    XmToggleButtonSetState(toggle_auto_item[0], auto_edit_grid, False);
    XmToggleButtonSetState(toggle_auto_item[1], auto_edit_grid_boundary, False);
    XmToggleButtonSetState(toggle_auto_item[2], auto_build_points, False);
    XmToggleButtonSetState(toggle_auto_item[3], auto_back_grid, False);
    XmToggleButtonSetState(toggle_auto_item[4], auto_back_boundary, False);
}

void auto_done_proc(void)
{
    XtUnmanageChild(auto_frame);
}

void create_auto_popup(void)
{
    int i, ac;
    Widget bt, lab, rb, fr, rc, rc2, sep;
    Arg a[2];
    XmString str;
    if (auto_frame) {
	XtRaise(auto_frame);
	set_auto_toggles();
	return;
    }
    auto_frame = XmCreateDialogShell(app_shell, "Scaling", NULL, 0);
    auto_panel = XmCreateRowColumn(auto_frame, "rc", NULL, 0);
    rc = XmCreateRowColumn(auto_panel, "rc", NULL, 0);
    XtVaCreateManagedWidget("Autoscale on:", xmLabelGadgetClass, rc, NULL);
    for (i = 0; i < NAUTOITEMS; i++) {
	toggle_auto_item[i] = XtVaCreateManagedWidget(auto_items[i],
				       xmToggleButtonWidgetClass, rc, NULL);
    }
    rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) set_autoscale, 0);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) auto_done_proc, NULL);
    XtManageChild(rc2);

    XtManageChild(rc);

    XtManageChild(auto_panel);
    XtManageChild(auto_frame);
    auto_edit_grid = 0;
    auto_edit_grid_boundary = 0;
    auto_back_grid = 0;
    auto_back_boundary = 0;
    auto_build_points = 0;
    set_auto_toggles();
}
