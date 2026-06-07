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
 * set bounds properties
 *
 */

#ifndef lint
static char RCSid[] = "$Id: boundwin.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget bounds_frame;
static Widget bounds_panel;

void create_bound_frame(void);

/*
 * Widget declarations
 */
static Widget *bounds_choice_item;
static Widget *bounds_display_item;
static Widget bounds_filled_item;
static Widget *bounds_color_choice_item;
static Widget *bounds_fill_choice_item;
static Widget *bounds_fillcolor_choice_item;
static Widget *bounds_lines_choice_item;
static Widget *bounds_linew_choice_item;
static Widget *bounds_applyto_choice_item;
static Widget *bounds_symbols_item;
static Widget *bounds_symfill_item;
static Widget bounds_symsize_item;

/*
 * Event and Notify proc declarations
 */
static Widget boundsio_frame;
static Widget *boundsio_item;

/*
 * Event and Notify proc declarations
 */
void boundsio_done_proc(void);
void boundsio_accept_proc(void);
static void bounds_Done_notify_proc(void);
static void bounds_isolines_notify_proc();
static void bounds_define_notify_proc(void);

void update_bounds_items(int gno);

/*
 * Create the bounds Frame and the bounds Panel
 */
extern Widget app_shell;

void set_curbound(Widget w, int cd)
{
    g[cg].curbound = cd;
    update_bounds_items(cg);
}

void create_boundio_frame(void)
{
    Widget wbut, rc;
    int i;
    setistop();
    if (!boundsio_frame) {
	boundsio_frame = XmCreateFileSelectionDialog(app_shell, "Read boundary", NULL, 0);
	XtVaSetValues(boundsio_frame, XmNdirMask, XmStringCreate("*", charset), NULL);

	rc = XmCreateRowColumn(boundsio_frame, "rc", NULL, 0);
	boundsio_item = CreatePanelChoice1(rc, "Read to bounds: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXBOUNDS; i++) {
	    XtAddCallback(boundsio_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curbound, (XtPointer) i);
	}
	XtManageChild(rc);
	XtAddCallback(boundsio_frame, XmNcancelCallback, (XtCallbackProc) boundsio_done_proc, NULL);
	XtAddCallback(boundsio_frame, XmNokCallback, (XtCallbackProc) boundsio_accept_proc, NULL);
    }
    XtRaise(boundsio_frame);
    update_bounds_items(cg);
}

void boundsio_done_proc(void)
{
    XtUnmanageChild(boundsio_frame);
}

void boundsio_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int boundno;
    Widget textw;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(boundsio_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    boundno = GetChoice(boundsio_item);
    if (bounds[boundno].active == ON) {
	if (!yesno("Boundary is active, kill it ?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(boundsio_frame);
    if (!readbinboundary(boundno, s)) {
	if (!readboundary2(boundno, s)) {
	    unset_wait_cursor(boundsio_frame);
	    sprintf(buf, "Error reading file %s", s);
	    errwin(buf);
	}
    }
    XtUnmanageChild(boundsio_frame);
    create_bound_frame();
    unset_wait_cursor(boundsio_frame);
}

void update_bounds_items(int gno)
{
    int gd, value;
    int c = g[gno].curbound;
    if (bounds_frame) {
	SetChoice(bounds_choice_item, c);
	switch (g[gno].bounds[c].display) {
	case OFF:
	    SetChoice(bounds_display_item, 0);
	    break;
	case ON:
	    SetChoice(bounds_display_item, 1);
	    break;
	case SYMBOL:
	    SetChoice(bounds_display_item, 2);
	    break;
	case NODES:
	    SetChoice(bounds_display_item, 3);
	    break;
	case XYSEG:
	    SetChoice(bounds_display_item, 4);
	    break;
	}
/*
	XmToggleButtonSetState(bounds_filled_item,
			       g[gno].bounds[c].display_boundsf == ON, False);
*/
	SetChoice(bounds_color_choice_item, g[gno].bounds[c].p.color);
	SetChoice(bounds_linew_choice_item, g[gno].bounds[c].p.linew - 1);
	SetChoice(bounds_lines_choice_item, g[gno].bounds[c].p.lines - 1);
	SetChoice(bounds_fill_choice_item, g[gno].bounds[c].p.fill == ON);
	SetChoice(bounds_fillcolor_choice_item, g[gno].bounds[c].p.fillcol == ON);
	value = g[gno].bounds[c].p.symsize;
	XtVaSetValues(bounds_symsize_item, XmNvalue, value, NULL);
    }
}

/*
 * Create the frame Widget and the frame Widget
 */
void create_bound_frame(void)
{
    extern Widget app_shell;
    Widget wbut, lab, rc, rc2, fr;
    int i;

    if (bounds_frame) {
	update_bounds_items(cg);
	XtRaise(bounds_frame);
	return;
    }
    bounds_frame = XmCreateDialogShell(app_shell, "Boundary setup", NULL, 0);
    handle_close(bounds_frame);
    bounds_panel = XmCreateRowColumn(bounds_frame, "panel", NULL, 0);
    XtVaSetValues(bounds_panel, XmNorientation, XmHORIZONTAL, NULL);

    rc2 = XmCreateRowColumn(bounds_panel, "rc2", NULL, 0);
    bounds_choice_item = CreatePanelChoice1(rc2, "Apply to boundary:", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", NULL, NULL);
    for (i = 0; i < MAXBOUNDS; i++) {
	XtAddCallback(bounds_choice_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curbound, (XtPointer) i);
    }

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    lab = XmCreateLabel(rc, "Display:", NULL, 0);
    XtManageChild(lab);
    bounds_display_item = CreatePanelChoice1(rc2, "Boundary:", 6, "Not displayed", "Lines", "Symbols", "Nodes", "Segments", NULL, NULL);
    bounds_filled_item = XmCreateToggleButton(rc, "Filled islands", NULL, 0);
    XtManageChild(bounds_filled_item);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(bounds_panel, "rc2", NULL, 0);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    bounds_color_choice_item = CreateColorChoice(rc, "Boundary color:", 1);
    bounds_linew_choice_item = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    bounds_lines_choice_item = CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    bounds_fill_choice_item = CreatePanelChoice1(rc, "Fill for islands:", 3, "None", "Filled", NULL, NULL);
    bounds_fillcolor_choice_item = CreateColorChoice(rc, "Island fill color:", 1);
    bounds_symbols_item = CreatePanelChoice1(rc, "Symbol:", 13, "No symbol", "Dot", "Circle", "Square", "Diamond", "Triangle up", "Triangle left", "Triangle down", "Triangle right", "Plus", "X", "Star", 0, 0);

    bounds_symfill_item = CreatePanelChoice1(rc, "Sym fill:", 4, "None", "Filled", "Opaque", NULL, NULL);

    XtVaCreateManagedWidget("Sym size:", xmLabelGadgetClass, rc, NULL);
    bounds_symsize_item = XtVaCreateManagedWidget("s", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);

    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(rc2, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) bounds_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Files...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_boundio_frame, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) bounds_Done_notify_proc, 0);
    XtManageChild(rc);
    XtManageChild(rc2);

    XtManageChild(bounds_panel);
    XtRaise(bounds_frame);
    update_bounds_items(cg);
}

/*
 * Notify and event procs
 */

static void bounds_Done_notify_proc(void)
{
    XtUnmanageChild(bounds_frame);
}

static void bounds_define_notify_proc(void)
{
    int c, bdisp, value;
    c = g[cg].curbound;
    bdisp = GetChoice(bounds_display_item);
    switch (bdisp) {
    case 0:
	g[cg].bounds[c].display = OFF;
	break;
    case 1:
	g[cg].bounds[c].display = ON;
	break;
    case 2:
	g[cg].bounds[c].display = SYMBOL;
	break;
    case 3:
	g[cg].bounds[c].display = NODES;
	break;
    case 4:
	g[cg].bounds[c].display = XYSEG;
	break;
    }
/*
	    g[cg].bounds[c].display_boundsf =
		    XmToggleButtonGetState(bounds_filled_item) ? ON : OFF;
*/
    g[cg].bounds[c].p.color = GetChoice(bounds_color_choice_item);
    g[cg].bounds[c].p.fill = GetChoice(bounds_fill_choice_item) ? ON : OFF;
    g[cg].bounds[c].p.fillcol = GetChoice(bounds_fillcolor_choice_item);
    g[cg].bounds[c].p.linew = GetChoice(bounds_linew_choice_item) + 1;
    g[cg].bounds[c].p.lines = GetChoice(bounds_lines_choice_item) + 1;
    XtVaGetValues(bounds_symsize_item, XmNvalue, &value, NULL);
    g[cg].bounds[c].p.symsize = value;
}
