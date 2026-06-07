/*
 * ACE/gredit - 2d finite element grid generation
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *                      All Rights Reserved.
 *
 */

/* $Id: wvwin.c,v 1.3 2007/02/21 00:21:21 pturner Exp $
 *
 * Set the world scale and viewport
 */

#ifndef lint
static char RCSid[] = "$Id: wvwin.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "motifinc.h"
#include "defines.h"
#include "globals.h"
 
extern Widget app_shell;
extern Widget canvas;
extern XmStringCharSet charset;

static Widget world_frame;
static Widget world_panel;

static Widget view_frame;
static Widget view_panel;

static Widget define_world_xg1;
static Widget define_world_xg2;
static Widget define_world_yg1;
static Widget define_world_yg2;

static Widget define_view_xv1;
static Widget define_view_xv2;
static Widget define_view_yv1;
static Widget define_view_yv2;

static Widget but1[2];

static void define_world_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void define_view_proc(Widget w, XtPointer client_data, XtPointer call_data);

double xv1 = 0.0, xv2 = 1.0, yv1 = 0.0, yv2 = 1.0;

void set_actioncb(Widget w, XtPointer client_data, XtPointer call_data)
{
    int func = (long) client_data;
    set_action(0);
    set_action(func);
}

/*
 * update the items in define world/view popup
 */
void update_world(void)
{
    char buf[1024];
    if (world_frame) {
	sprintf(buf, "%.9lg", xg1);
	xv_setstr(define_world_xg1, buf);
	sprintf(buf, "%.9lg", xg2);
	xv_setstr(define_world_xg2, buf);
	sprintf(buf, "%.9lg", yg1);
	xv_setstr(define_world_yg1, buf);
	sprintf(buf, "%.9lg", yg2);
	xv_setstr(define_world_yg2, buf);
    }
}

static void define_world_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char val[80];
    int i, which, ming, maxg, errpos;
    double x, y, a, b, c, d;
    double tx1, tx2, ty1, ty2;
    double tm1, tm2;
    strcpy(val, (char *) xv_getstr(define_world_xg1));
    xg1 = atof(val);
    strcpy(val, (char *) xv_getstr(define_world_xg2));
    xg2 = atof(val);
    strcpy(val, (char *) xv_getstr(define_world_yg1));
    yg1 = atof(val);
    strcpy(val, (char *) xv_getstr(define_world_yg2));
    yg2 = atof(val);
    do_drawgrid();
}

void update_world_proc(void)
{
    update_world();
}

/*
 * Create the world Frame and the world Panel
 */
void create_world_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc, rc2, fr;
    Widget wlabel;
    Widget buts[3];
    set_wait_cursor();
    if (world_frame == NULL) {
	char *blabels[3];
	blabels[0] = "Accept";
	blabels[1] = "Update";
	blabels[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	world_frame = XmCreateDialogShell(app_shell, "World", NULL, 0);
	handle_close(world_frame);
	XtVaSetValues(world_frame, XmNx, x, XmNy, y, NULL);
	world_panel = XmCreateRowColumn(world_frame, "world_rc", NULL, 0);

	rc2 = XmCreateRowColumn(world_panel, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	fr = XmCreateFrame(rc2, "fr", NULL, 0);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	wlabel = XtVaCreateManagedWidget("World (axis scaling)", xmLabelWidgetClass, rc, NULL);
	define_world_xg1 = CreateTextItem2(rc, 10, "Xmin:");
	define_world_xg2 = CreateTextItem2(rc, 10, "Xmax:");
	define_world_yg1 = CreateTextItem2(rc, 10, "Ymin:");
	define_world_yg2 = CreateTextItem2(rc, 10, "Ymax:");
	XtManageChild(rc);
	XtManageChild(fr);
	XtManageChild(rc2);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, world_panel, NULL);

	CreateCommandButtons(world_panel, 3, buts, blabels);
	XtAddCallback(buts[0], XmNactivateCallback, (XtCallbackProc) define_world_proc, (XtPointer) NULL);
	XtAddCallback(buts[1], XmNactivateCallback, (XtCallbackProc) update_world_proc, (XtPointer) NULL);
	XtAddCallback(buts[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) world_frame);

	XtManageChild(world_panel);
    }
    XtRaise(world_frame);
    update_world();
    unset_wait_cursor();
}				/* end create_world_panel */

void update_view()
{
    char buf[1024];
    if (view_frame) {
	sprintf(buf, "%.9lg", xv1);
	xv_setstr(define_view_xv1, buf);
	sprintf(buf, "%.9lg", xv2);
	xv_setstr(define_view_xv2, buf);
	sprintf(buf, "%.9lg", yv1);
	xv_setstr(define_view_yv1, buf);
	sprintf(buf, "%.9lg", yv2);
	xv_setstr(define_view_yv2, buf);
    }
}

double fround(double x, double r)
{
    x = x + r * 0.5;
    x = (long) (x / r);
    return x * r;
}

static void define_view_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char val[80];
    strcpy(val, xv_getstr(define_view_xv1));
    xv1 = atof(val);
    strcpy(val, xv_getstr(define_view_xv2));
    xv2 = atof(val);
    strcpy(val, xv_getstr(define_view_yv1));
    yv1 = atof(val);
    strcpy(val, xv_getstr(define_view_yv2));
    yv2 = atof(val);
    do_drawgrid();
}

static void define_viewm_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    set_action(0);
    set_action(VIEW_1ST);
}

void create_view_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc;
    Widget buts[3];
    set_wait_cursor();
    if (view_frame == NULL) {
	char *blabels[3];
	blabels[0] = "Accept";
	blabels[1] = "Pick";
	blabels[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	view_frame = XmCreateDialogShell(app_shell, "Viewports", NULL, 0);
	handle_close(view_frame);
	XtVaSetValues(view_frame, XmNx, x, XmNy, y, NULL);
	view_panel = XmCreateRowColumn(view_frame, "view_rc", NULL, 0);

	define_view_xv1 = CreateTextItem2(view_panel, 10, "Xmin:");
	define_view_xv2 = CreateTextItem2(view_panel, 10, "Xmax:");
	define_view_yv1 = CreateTextItem2(view_panel, 10, "Ymin:");
	define_view_yv2 = CreateTextItem2(view_panel, 10, "Ymax:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, view_panel, NULL);

	CreateCommandButtons(view_panel, 3, buts, blabels);
	XtAddCallback(buts[0], XmNactivateCallback, (XtCallbackProc) define_view_proc, (XtPointer) NULL);
	XtAddCallback(buts[1], XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) VIEW_1ST);
	XtAddCallback(buts[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) view_frame);

	XtManageChild(view_panel);
    }
    XtRaise(view_frame);
    update_view();
    unset_wait_cursor();
}
