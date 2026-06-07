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
 * zoom Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: zoomwin.c,v 1.2 2003/07/24 15:44:07 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget zoom_frame;
static Widget zoom_panel;

void create_zoom_frame(void);
void update_zoom(void);

void set_zoomloc(void);

void define_zoom(int gno, int zno, double wx1, double wy1, double wx2, double wy2, double wx, double wy)
{
    double vx1, vy1;
    double vx2, vy2;
    double vx, vy;
    world2view(wx1, wy1, &vx1, &vy1);
    world2view(wx2, wy2, &vx2, &vy2);
    vx = (vx2 - vx1);
    vy = (vy2 - vy1);
    g[gno].zbox[zno].vx = vx;
    g[gno].zbox[zno].vy = vy;
    g[gno].zbox[zno].wx1 = wx1;
    g[gno].zbox[zno].wx2 = wx2;
    g[gno].zbox[zno].wy1 = wy1;
    g[gno].zbox[zno].wy2 = wy2;
    g[gno].zbox[zno].locx = wx;
    g[gno].zbox[zno].locy = wy;
}

/*
 * Panel item declarations
 */
static Widget zoom_toggle_item;
static Widget zoom_togglemarker_item;
static Widget zoom_type_item;
static Widget *zoom_which_item;
static Widget *zoom_precx_item;
static Widget *zoom_precy_item;
static Widget *zoom_color_item;
static Widget *zoom_fill_item;
static Widget *zoom_linew_item;
static Widget *zoom_attach_item;
static Widget *zoom_expand_item;

/*
 * Event and Notify proc declarations
 */

void zoom_done_proc(void);
void zoom_accept_proc(void);

static void set_zoomlocCB(void)
{
    curzoom = GetChoice(zoom_which_item);
    set_zoomloc();
}

static void set_curzoom(Widget w, int cd)
{
    curzoom = cd;
    update_zoom();
}

void toggle_zoom(void)
{
/*
     g[cg].display_zoom = XmToggleButtonGetState(zoom_toggle_item);
*/
}

/*
 * Create the zoom Frame and the zoom Panel
 */

extern Widget app_shell;

void create_zoom_frame(void)
{
    Widget wbut, fr, bb, sep, rc, rc2;
    int i;
    setistop();
    if (!zoom_frame) {
	zoom_frame = XmCreateDialogShell(app_shell, "Zoom box", NULL, 0);
	handle_close(zoom_frame);
	zoom_panel = XmCreateRowColumn(zoom_frame, "zoomrc", NULL, 0);
	zoom_which_item = CreatePanelChoice2(zoom_panel, "Zoom box: ", 2, 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXZOOMBOXES; i++) {
	    XtAddCallback(zoom_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curzoom, (XtPointer) i);
	}
	zoom_toggle_item = XtVaCreateManagedWidget("Display zoom region", xmToggleButtonWidgetClass, zoom_panel, NULL);
	zoom_togglemarker_item = XtVaCreateManagedWidget("Display zoom marker", xmToggleButtonWidgetClass, zoom_panel, NULL);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, zoom_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);

	zoom_color_item = CreateColorChoice(bb, "Frame color:", 1);
	zoom_linew_item = CreatePanelChoice2(bb, "Frame line width: ", 2, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

	zoom_fill_item = CreateColorChoice(bb, "Fill color:", 1);
	zoom_precx_item = CreatePanelChoice2(bb, "Precision X labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	zoom_precy_item = CreatePanelChoice2(bb, "Precision Y labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	zoom_attach_item = CreatePanelChoice2(bb, "Attach to: ", 1, 5, "SW", "SE", "NW", "NE", 0, 0);
	zoom_expand_item = CreatePanelChoice2(bb, "Expand: ", 2, 8, "2x", "3x", "4x", "5x", "6x", "10x", "20x", 0, 0);

	rc = XmCreateRowColumn(zoom_panel, "zoomrc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) zoom_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_zoomlocCB, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) zoom_done_proc, NULL);
	XtManageChild(rc);

	XtManageChild(zoom_panel);
    }
    update_zoom();
    XtRaise(zoom_frame);
}

void update_zoom(void)
{
    char buf[256];
    int itmp;
    if (zoom_frame) {
	SetChoice(zoom_which_item, curzoom);

	XmToggleButtonSetState(zoom_toggle_item, g[cg].zbox[curzoom].display == ON, False);
	XmToggleButtonSetState(zoom_togglemarker_item, g[cg].zbox[curzoom].display_marker == ON, False);
	SetChoice(zoom_color_item, g[cg].zbox[curzoom].p.color);
	SetChoice(zoom_fill_item, g[cg].zbox[curzoom].p.fillcol);
	SetChoice(zoom_linew_item, g[cg].zbox[curzoom].p.linew);
	SetChoice(zoom_precx_item, g[cg].zbox[curzoom].precx + 1);
	SetChoice(zoom_precy_item, g[cg].zbox[curzoom].precy + 1);
	SetChoice(zoom_attach_item, g[cg].zbox[curzoom].attach);
	itmp = g[cg].zbox[curzoom].expand;
	if (itmp <= 1) {
	    itmp = 2;
	}
	if (itmp == 10) {	/* last 2 items are 10 & 20 */
	    itmp = 5;
	} else if (itmp == 20) {
	    itmp = 6;
	} else {
	    itmp = itmp - 2;
	}
	SetChoice(zoom_expand_item, itmp);
    }
}

void zoom_done_proc(void)
{
    XtUnmanageChild(zoom_frame);
}

void zoom_accept_proc(void)
{
    int h, itmp;

    curzoom = h = GetChoice(zoom_which_item);

    g[cg].zbox[curzoom].display = XmToggleButtonGetState(zoom_toggle_item) ? ON : OFF;
    g[cg].zbox[curzoom].display_marker = XmToggleButtonGetState(zoom_togglemarker_item) ? ON : OFF;
    g[cg].zbox[curzoom].p.color = GetChoice(zoom_color_item);
    g[cg].zbox[curzoom].p.fillcol = GetChoice(zoom_fill_item);
    g[cg].zbox[curzoom].p.linew = GetChoice(zoom_linew_item);
    g[cg].zbox[curzoom].precx = GetChoice(zoom_precx_item) - 1;
    g[cg].zbox[curzoom].precy = GetChoice(zoom_precy_item) - 1;
    g[cg].zbox[curzoom].attach = GetChoice(zoom_attach_item);
    itmp = GetChoice(zoom_expand_item);
    if (itmp == 5) {
	itmp = 10;
    } else if (itmp == 6) {
	itmp = 20;
    } else {
	itmp = itmp + 2;
    }
    g[cg].zbox[curzoom].expand = itmp;
}
