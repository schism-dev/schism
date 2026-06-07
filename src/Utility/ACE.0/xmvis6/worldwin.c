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
 * Set the world scale, tick marks, and viewport
 */

#ifndef lint
static char RCSid[] = "$Id: worldwin.c,v 1.2 2003/07/24 15:44:07 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

static Widget world_frame = (Widget) NULL;
static Widget world_panel;

static Widget view_frame = (Widget) NULL;
static Widget view_panel;

static Widget arrange_frame = (Widget) NULL;
static Widget arrange_panel;

static Widget autos_frame = (Widget) NULL;
static Widget autos_panel;

/*
 * Panel item declarations
 */
static Widget define_world_xg1;
static Widget define_world_xg2;
static Widget define_world_yg1;
static Widget define_world_yg2;
static Widget define_world_xminor;
static Widget define_world_yminor;
static Widget define_world_ymajor;
static Widget define_world_xmajor;
static Widget *world_applyto_item;

static Widget define_view_xv1;
static Widget define_view_xv2;
static Widget define_view_yv1;
static Widget define_view_yv2;
static Widget *view_applyto_item;

static Widget *arrange_rows_item;
static Widget *arrange_cols_item;
static Widget arrange_vgap_item;
static Widget arrange_hgap_item;
static Widget arrange_startx_item;
static Widget arrange_starty_item;
static Widget arrange_widthx_item;
static Widget arrange_widthy_item;
static Widget *arrange_packed_item;
static Widget arrange_applyto_item;

static Widget *autos_on_item;
static Widget *autos_using_item;
static Widget *autos_method_item;
static Widget *autos_nticksx_item;
static Widget *autos_nticksy_item;
static Widget *autos_applyto_item;

#define PANEL_VALUE 500

char *xv_getstr(Widget w);

int maxmajorticks = 500;
int maxminorticks = 1000;

/*
 * update the items in define world/view popup
 */
void update_world(int gno)
{
    if (world_frame) {
	sprintf(buf, "%.9lg", g[gno].w.xg1);
	xv_setstr(define_world_xg1, buf);
	sprintf(buf, "%.9lg", g[gno].w.xg2);
	xv_setstr(define_world_xg2, buf);
	sprintf(buf, "%.9lg", g[gno].w.yg1);
	xv_setstr(define_world_yg1, buf);
	sprintf(buf, "%.9lg", g[gno].w.yg2);
	xv_setstr(define_world_yg2, buf);
	sprintf(buf, "%.9lg", g[gno].t[0].tmajor);
	xv_setstr(define_world_xmajor, buf);
	sprintf(buf, "%.9lg", g[gno].t[0].tminor);
	xv_setstr(define_world_xminor, buf);
	sprintf(buf, "%.9lg", g[gno].t[1].tmajor);
	xv_setstr(define_world_ymajor, buf);
	sprintf(buf, "%.9lg", g[gno].t[1].tminor);
	xv_setstr(define_world_yminor, buf);
    }
}

static void define_world_proc(void)
{
    char val[80];
    int i, which, ming, maxg, errpos;
    double x, y, a, b, c, d;
    double tx1, tx2, ty1, ty2;
    double tm1, tm2;
    extern double result;	/* result passed from the expression
				 * interpreter */

    if (which = (int) GetChoice(world_applyto_item)) {
	ming = 0;
	maxg = maxgraph - 1;
    } else {
	ming = cg;
	maxg = cg;
    }
    if (ming == cg && maxg == cg) {
	if (!isactive_graph(cg)) {
	    errwin("Current graph is not active!");
	    return;
	}
    }
    for (i = ming; i <= maxg; i++) {
	if (isactive_graph(i)) {
	    x = (g[i].w.xg2 - g[i].w.xg1);	/* DX */
	    y = (g[i].w.yg2 - g[i].w.yg1);	/* DY */
	    a = g[i].w.xg1;	/* XMIN */
	    b = g[i].w.yg1;	/* XMAX */
	    c = g[i].w.xg2;	/* YMIN */
	    d = g[i].w.yg2;	/* YMAX */
	    if (which <= 1 || which == 2) {
		strcpy(val, (char *) xv_getstr(define_world_xg1));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tx1 = result;

		strcpy(val, (char *) xv_getstr(define_world_xg2));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tx2 = result;

		strcpy(val, (char *) xv_getstr(define_world_xmajor));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tm1 = result;

		strcpy(val, (char *) xv_getstr(define_world_xminor));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tm2 = result;
		if (tx1 < tx2) {
		    g[i].w.xg1 = tx1;
		    g[i].w.xg2 = tx2;
		} else {
		    errwin("World scaling along X improperly set, scale not changed");
		    return;
		}
		if (tm1 > 0.0) {
		    if (check_nticks(i, X_AXIS, g[i].w.xg1, g[i].w.xg2, tm1, maxmajorticks)) {
			g[i].t[X_AXIS].tmajor = tm1;
		    } else {
			sprintf(buf, "Too many major ticks along X-axis, %d", maxmajorticks);
			errwin(buf);
			return;
		    }
		    g[i].t[0].tmajor = tm1;
		} else {
		    errwin("Major ticks must be > 0.0");
		    return;
		}
		if (tm2 > 0.0) {
		    if (check_nticks(i, X_AXIS, g[i].w.xg1, g[i].w.xg2, tm2, maxminorticks)) {
			g[i].t[X_AXIS].tminor = tm2;
		    } else {
			sprintf(buf, "Too many minor ticks along X-axis, %d", maxminorticks);
			errwin(buf);
			return;
		    }
		    g[i].t[0].tminor = tm2;
		} else {
		    errwin("Minor ticks must be > 0.0");
		    return;
		}
	    }
	    if (which <= 1 || which == 3) {

		strcpy(val, (char *) xv_getstr(define_world_yg1));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		ty1 = result;
		strcpy(val, (char *) xv_getstr(define_world_yg2));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		ty2 = result;

		strcpy(val, (char *) xv_getstr(define_world_ymajor));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tm1 = result;

		strcpy(val, (char *) xv_getstr(define_world_yminor));
		fixupstr(val);
		scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
		if (errpos) {
		    return;
		}
		tm2 = result;

		if (ty1 < ty2) {
		    g[i].w.yg1 = ty1;
		    g[i].w.yg2 = ty2;
		} else {
		    errwin("World scaling along Y improperly set, scale not changed");
		    return;
		}
		if (tm1 > 0.0) {
		    if (check_nticks(i, Y_AXIS, g[i].w.yg1, g[i].w.yg2, tm1, maxmajorticks)) {
			g[i].t[Y_AXIS].tmajor = tm1;
		    } else {
			sprintf(buf, "Too many minor ticks along Y-axis, %d", maxmajorticks);
			errwin(buf);
			return;
		    }
		} else {
		    errwin("Major ticks must be > 0.0");
		    return;
		}
		if (tm2 > 0.0) {
		    if (check_nticks(i, Y_AXIS, g[i].w.yg1, g[i].w.yg2, tm2, maxminorticks)) {
			g[i].t[Y_AXIS].tminor = tm2;
		    } else {
			sprintf(buf, "Too many minor ticks along Y-axis, %d", maxminorticks);
			errwin(buf);
			return;
		    }
		} else {
		    errwin("Minor ticks must be > 0.0");
		    return;
		}
	    }
	}
    }
    drawgraph();
}

void update_world_proc(void)
{
    update_world(cg);
}

/*
static void autoscale_type_proc(w, cd)
    Widget w;
    int cd;
{
    if (cd == 0) {
	g[cg].auto_type = AUTO;
    } else {
	g[cg].auto_type = SPEC;
    }
}

static void autoscale_set_proc()
{
    int value;
    value = (int) GetChoice(autoscale_set_item);
    if (isactive(cg, value)) {
	defaultsetgraph(cg, value);
	default_axis(cg, g[cg].auto_type, X_AXIS);
	default_axis(cg, g[cg].auto_type, ZX_AXIS);
	default_axis(cg, g[cg].auto_type, Y_AXIS);
	default_axis(cg, g[cg].auto_type, ZY_AXIS);
	update_world(cg);
	drawgraph();
    } else {
	errwin("Set not active!");
    }
}
*/

static void do_zoom_proc(void)
{
    set_action(0);
    set_action(ZOOM_1ST);
}

static void do_view_proc(void)
{
    set_action(0);
    set_action(VIEW_1ST);
}

/*
 * Create the world Frame and the world Panel
 */
void create_world_frame(void)
{
    extern Widget app_shell;
    int x, y;
    Widget wbut, rc, rc2, fr;
    Widget wlabel;

    if (world_frame) {
	update_world(cg);
	XtRaise(world_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    world_frame = XmCreateDialogShell(app_shell, "World", NULL, 0);
    handle_close(world_frame);
    XtVaSetValues(world_frame, XmNx, x, XmNy, y, NULL);
    world_panel = XmCreateRowColumn(world_frame, "world_rc", NULL, 0);

    rc2 = XmCreateRowColumn(world_panel, "rc", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    wlabel = XtVaCreateManagedWidget("World (axis scaling)", xmLabelGadgetClass, rc, NULL);
    define_world_xg1 = (Widget) CreateTextItem2(rc, 10, "Xmin:");
    define_world_xg2 = (Widget) CreateTextItem2(rc, 10, "Xmax:");
    define_world_yg1 = (Widget) CreateTextItem2(rc, 10, "Ymin:");
    define_world_yg2 = (Widget) CreateTextItem2(rc, 10, "Ymax:");
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    wlabel = XtVaCreateManagedWidget("Tick spacing", xmLabelGadgetClass, rc, NULL);
    define_world_xmajor = (Widget) CreateTextItem2(rc, 10, "X-major:");
    define_world_xminor = (Widget) CreateTextItem2(rc, 10, "X-minor:");
    define_world_ymajor = (Widget) CreateTextItem2(rc, 10, "Y-major:");
    define_world_yminor = (Widget) CreateTextItem2(rc, 10, "Y-minor:");
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    world_applyto_item = (Widget *) CreatePanelChoice1(world_panel, "Apply to:", 5, "Current graph", "All active graphs", "X, all active graphs", "Y, all active graphs", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, world_panel, NULL);
    rc = XmCreateRowColumn(world_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_world_proc, 0);
    wbut = XtVaCreateManagedWidget("Update World", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) update_world_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, world_frame);
    XtManageChild(rc);

    update_world(cg);
    XtManageChild(world_panel);
    XtManageChild(world_frame);
}				/* end create_world_panel */

void update_view(int gno)
{
    if (view_frame) {
	sprintf(buf, "%.9lg", g[gno].v.xv1);
	xv_setstr(define_view_xv1, buf);
	sprintf(buf, "%.9lg", g[gno].v.xv2);
	xv_setstr(define_view_xv2, buf);
	sprintf(buf, "%.9lg", g[gno].v.yv1);
	xv_setstr(define_view_yv1, buf);
	sprintf(buf, "%.9lg", g[gno].v.yv2);
	xv_setstr(define_view_yv2, buf);
    }
}

double fround(double x, double r)
{
    x = x + r * 0.5;
    x = (long) (x / r);
    return x * r;
}

static void define_view_proc(void)
{
    char val[80];
    double tmpx1, tmpx2;
    double tmpy1, tmpy2;
    double r;
    int i, snap = 0, itmp, ming, maxg, which, ierr = 0;

/*
    snap = (int) xv_get(view_snap_item, PANEL_VALUE);
*/
    if (which = (int) GetChoice(view_applyto_item)) {
	ming = 0;
	maxg = maxgraph - 1;
    } else {
	ming = cg;
	maxg = cg;
    }
    if (ming == cg && maxg == cg) {
	if (!isactive_graph(cg)) {
	    errwin("Current graph is not active!");
	    return;
	}
    }
    strcpy(val, xv_getstr(define_view_xv1));
    tmpx1 = atof(val);
    strcpy(val, xv_getstr(define_view_xv2));
    tmpx2 = atof(val);
    strcpy(val, xv_getstr(define_view_yv1));
    tmpy1 = atof(val);
    strcpy(val, xv_getstr(define_view_yv2));
    tmpy2 = atof(val);
    if (!fbounds(tmpx1, 0.0, 1.0, "View xmin")) {
	ierr = 1;
    } else if (!fbounds(tmpx2, tmpx1, 1.0, "View xmax")) {
	ierr = 1;
    } else if (!fbounds(tmpy1, 0.0, 1.0, "View ymin")) {
	ierr = 1;
    } else if (!fbounds(tmpy2, tmpy1, 1.0, "View ymax")) {
	ierr = 1;
    } else if (tmpx2 - tmpx1 <= 0.0) {
	ierr = 1;
    } else if (tmpy2 - tmpy1 <= 0.0) {
	ierr = 1;
    }
    if (ierr) {
	errwin("Viewport not set");
	return;
    }
    switch (snap) {
    case 0:			/* no snap */
	break;
    case 1:			/* .1 */
	r = 0.1;
	break;
    case 2:			/* .05 */
	r = 0.05;
	break;
    case 3:			/* .01 */
	r = 0.01;
	break;
    case 4:			/* .005 */
	r = 0.005;
	break;
    case 5:			/* .001 */
	r = 0.001;
	break;
    }
    if (snap) {
	tmpx1 = fround(tmpx1, r);
	tmpx2 = fround(tmpx2, r);
	tmpy1 = fround(tmpy1, r);
	tmpy2 = fround(tmpy2, r);
	update_view(cg);
    }
    for (i = ming; i <= maxg; i++) {
	if (isactive_graph(i)) {
	    if (which <= 1 || which == 2) {
		g[i].v.xv1 = tmpx1;
		g[i].v.xv2 = tmpx2;
	    }
	    if (which <= 1 || which == 3) {
		g[i].v.yv1 = tmpy1;
		g[i].v.yv2 = tmpy2;
	    }
	}
    }
    drawgraph();
}

static void define_viewm_proc(void)
{
    set_action(0);
    set_action(VIEW_1ST);
}

static void snap_proc(void)
{
}

void create_view_frame(void)
{
    extern Widget app_shell;
    int x, y;
    Widget wbut, rc;

    if (view_frame) {
	update_view(cg);
	XtRaise(view_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    view_frame = XmCreateDialogShell(app_shell, "Viewports", NULL, 0);
    handle_close(view_frame);
    XtVaSetValues(view_frame, XmNx, x, XmNy, y, NULL);
    view_panel = XmCreateRowColumn(view_frame, "view_rc", NULL, 0);

    define_view_xv1 = (Widget) CreateTextItem2(view_panel, 10, "Xmin:");
    define_view_xv2 = (Widget) CreateTextItem2(view_panel, 10, "Xmax:");
    define_view_yv1 = (Widget) CreateTextItem2(view_panel, 10, "Ymin:");
    define_view_yv2 = (Widget) CreateTextItem2(view_panel, 10, "Ymax:");

    view_applyto_item = (Widget *) CreatePanelChoice1(view_panel, "Apply to:", 5, "Current graph", "All active graphs", "X, all active graphs", "Y, all active graphs", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, view_panel, NULL);
    rc = XmCreateRowColumn(view_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_view_proc, 0);
    wbut = XtVaCreateManagedWidget("Pick view", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_view_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, view_frame);
    XtManageChild(rc);

    update_view(cg);
    XtManageChild(view_panel);
    XtManageChild(view_frame);
}

static void define_arrange_proc(void)
{
    int nrows, ncols, pack;
    double vgap, hgap, sx, sy, wx, wy;

    nrows = (int) GetChoice(arrange_rows_item) + 1;
    ncols = (int) GetChoice(arrange_cols_item) + 1;
    if (nrows == 1 && ncols == 1) {
	errwin("For a single graph use View/Viewport...");
	return;
    }
    pack = (int) GetChoice(arrange_packed_item);
    vgap = atof((char *) xv_getstr(arrange_vgap_item));
    hgap = atof((char *) xv_getstr(arrange_hgap_item));
    sx = atof((char *) xv_getstr(arrange_startx_item));
    sy = atof((char *) xv_getstr(arrange_starty_item));
    wx = atof((char *) xv_getstr(arrange_widthx_item));
    wy = atof((char *) xv_getstr(arrange_widthy_item));
    if (wx <= 0.0) {
	errwin("Graph width must be > 0.0");
	return;
    }
    if (wy <= 0.0) {
	errwin("Graph height must be > 0.0");
	return;
    }
    define_arrange(nrows, ncols, pack, vgap, hgap, sx, sy, wx, wy);
}

update_arrange(void)
{
    xv_setstr(arrange_vgap_item, "0.05");
    xv_setstr(arrange_hgap_item, "0.05");
    xv_setstr(arrange_startx_item, "0.1");
    xv_setstr(arrange_starty_item, "0.1");
    xv_setstr(arrange_widthx_item, "0.8");
    xv_setstr(arrange_widthy_item, "0.8");
}

void create_arrange_frame(void)
{
    extern Widget app_shell;
    int x, y;
    Widget wbut, rc;
    Widget wlabel;

    if (arrange_frame) {
/*	update_arrange(); */
	XtRaise(arrange_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    arrange_frame = XmCreateDialogShell(app_shell, "Arrange graphs", NULL, 0);
    handle_close(arrange_frame);
    XtVaSetValues(arrange_frame, XmNx, x, XmNy, y, NULL);
    arrange_panel = XmCreateRowColumn(arrange_frame, "arrange_rc", NULL, 0);

    arrange_rows_item = (Widget *) CreatePanelChoice1(arrange_panel, "Rows:", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", NULL, NULL);
    arrange_cols_item = (Widget *) CreatePanelChoice1(arrange_panel, "Columns:", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", NULL, NULL);
    arrange_packed_item = (Widget *) CreatePanelChoice1(arrange_panel, "Packing:", 5, "None", "Horizontal", "Vertical", "Both", NULL, NULL);

    arrange_vgap_item = (Widget) CreateTextItem2(arrange_panel, 10, "Vertical gap:");
    arrange_hgap_item = (Widget) CreateTextItem2(arrange_panel, 10, "Horizontal gap:");
    arrange_startx_item = (Widget) CreateTextItem2(arrange_panel, 10, "Start at X =");
    arrange_starty_item = (Widget) CreateTextItem2(arrange_panel, 10, "Start at Y =");
    arrange_widthx_item = (Widget) CreateTextItem2(arrange_panel, 10, "Graph width:");
    arrange_widthy_item = (Widget) CreateTextItem2(arrange_panel, 10, "Graph height:");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, arrange_panel, NULL);
    rc = XmCreateRowColumn(arrange_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_arrange_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, arrange_frame);
    XtManageChild(rc);

    update_arrange();
    XtManageChild(arrange_panel);
    XtManageChild(arrange_frame);
}

/*
 * autoscale popup
 */

static void define_autos_proc(void)
{
    int i, aon, ameth, antx, anty, au, ap, ming, maxg;
    aon = (int) GetChoice(autos_on_item) - 4;
    au = (int) GetChoice(autos_using_item) - 1;
    ap = (int) GetChoice(autos_applyto_item);
    ameth = (int) GetChoice(autos_method_item);
    antx = (int) GetChoice(autos_nticksx_item);
    anty = (int) GetChoice(autos_nticksy_item);
    define_autos(aon, au, ap, ameth, antx, anty);
}

void update_autos(int gno)
{
    if (autos_frame) {
	SetChoice(autos_method_item, g[gno].auto_type == SPEC);
	SetChoice(autos_nticksx_item, g[gno].t[0].t_num - 2);
	SetChoice(autos_nticksy_item, g[gno].t[1].t_num - 2);
    }
}

void create_autos_frame(void)
{
    extern Widget app_shell;
    int x, y;
    Widget wbut, rc;

    if (autos_frame) {
	update_autos(cg);
	XtRaise(autos_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    autos_frame = XmCreateDialogShell(app_shell, "Autoscale graphs", NULL, 0);
    handle_close(autos_frame);
    XtVaSetValues(autos_frame, XmNx, x, XmNy, y, NULL);
    autos_panel = XmCreateRowColumn(autos_frame, "autos_rc", NULL, 0);

    autos_on_item = (Widget *) CreatePanelChoice1(autos_panel, "Autoscale axis:", 11, "None", "All", "All X-axes", "All Y-axes", "X-axis", "Y-axis", "Zero X-axis", "Zero Y-axis", "Alternate X-axis", "Alternate Y-axis", NULL, NULL);

    autos_using_item = CreatePanelChoice2(autos_panel,
					  "Using set:", 3, 32, "All", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", 0, 0);

    autos_method_item = (Widget *) CreatePanelChoice1(autos_panel, "Autoscale type:", 3, "Heckbert", "Fixed", NULL, NULL);
    autos_nticksx_item = (Widget *) CreatePanelChoice1(autos_panel, "Number of ticks in X:", 12, "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", NULL, NULL);
    autos_nticksy_item = (Widget *) CreatePanelChoice1(autos_panel, "Number of ticks in Y:", 12, "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", NULL, NULL);
    autos_applyto_item = (Widget *) CreatePanelChoice1(autos_panel, "Apply to:", 3, "Current graph", "All active graphs", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, autos_panel, NULL);
    rc = XmCreateRowColumn(autos_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_autos_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, autos_frame);
    XtManageChild(rc);

    update_autos(cg);
    XtManageChild(autos_panel);
    XtManageChild(autos_frame);
}
