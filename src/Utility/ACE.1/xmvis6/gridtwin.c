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
 * set grid properties
 *
 */

#ifndef lint
static char RCSid[] = "$Id: gridtwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget gridt_frame;
static Widget gridt_panel;

void create_gridt_frame(void);

/*
 * Widget item declarations
 */
static Widget *gridt_choice_item;
static Widget gridt_display_item;
static Widget gridt_bath_item;
static Widget gridt_nodes_item;
static Widget gridt_elements_item;
static Widget gridt_depths_item;
/*
static Widget gridt_cour_item;
static Widget gridt_dimw_item;
*/
static Widget gridt_filled_item;
static Widget *gridt_color_choice_item;
static Widget *gridt_fillcolor_choice_item;
static Widget *gridt_lines_choice_item;
static Widget *gridt_linew_choice_item;
static Widget *gridt_applyto_choice_item;

/*
 * Event and Notify proc declarations
 */
static Widget gridio_frame;
static Widget *gridio_item;

/*
 * Event and Notify proc declarations
 */
static void gridio_done_proc(void);
static void gridio_accept_proc(void);
static void gridt_Done_notify_proc(void);
static void gridt_isolines_notify_proc(void);
static void gridt_define_notify_proc(void);

void update_gridt_items(int gno);

/*
 * Create the grid Frame and the grid Panel
 */
extern Widget app_shell;

static void set_curgridt(Widget w, int cd)
{
    g[cg].curgridt = cd;
    update_gridt_items(cg);
}

void create_gridtio_frame(void)
{
    Widget wbut, rc;
    int i;
    setistop();
    if (!gridio_frame) {
	gridio_frame = XmCreateFileSelectionDialog(app_shell, "Read data", NULL, 0);
	XtVaSetValues(gridio_frame, XmNdirMask, XmStringCreate("*.gr*", charset), NULL);

	rc = XmCreateRowColumn(gridio_frame, "rc", NULL, 0);
	gridio_item = CreatePanelChoice1(rc, "Read to grid: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < MAXGRIDS; i++) {
	    XtAddCallback(gridio_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curgridt, (XtPointer) i);
	}
	XtManageChild(rc);
	XtAddCallback(gridio_frame, XmNcancelCallback, (XtCallbackProc) gridio_done_proc, NULL);
	XtAddCallback(gridio_frame, XmNokCallback, (XtCallbackProc) gridio_accept_proc, NULL);
    }
    XtRaise(gridio_frame);
    update_gridt_items(cg);
}

static void gridio_done_proc(void)
{
    XtUnmanageChild(gridio_frame);
}

static void gridio_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int gridno;
    Widget textw;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(gridio_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    gridno = GetChoice(gridio_item);
    if (gridt[gridno].active == ON) {
	if (!yesno("Grid is active, kill it ?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(gridio_frame);
    if (!readgridt(gridno, s)) {
	sprintf(buf, "Error reading file %s", s);
	errwin(buf);
    } else {
	set_clock(0, gridt[gridno].start, gridt[gridno].stop, gridt[gridno].step, gridt[gridno].nsteps);
	load_clock(GRID, gridno);
	XtUnmanageChild(gridio_frame);
	create_gridt_frame();
    }
    unset_wait_cursor(gridio_frame);
}

void update_gridt_items(int gno)
{
    int gd;
    int c = g[gno].curgridt;
    if (gridt_frame) {
	XmToggleButtonSetState(gridt_display_item, g[gno].gridt[c].display == ON, False);
	XmToggleButtonSetState(gridt_bath_item, g[gno].gridt[c].display_bath == ON, False);
	XmToggleButtonSetState(gridt_nodes_item, g[gno].gridt[c].display_nodes == ON, False);
	XmToggleButtonSetState(gridt_elements_item, g[gno].gridt[c].display_elements == ON, False);
	XmToggleButtonSetState(gridt_depths_item, g[gno].gridt[c].display_depths == ON, False);
	XmToggleButtonSetState(gridt_filled_item, g[gno].gridt[c].display_gridf == ON, False);
	SetChoice(gridt_color_choice_item, g[gno].gridt[c].p.color);
	SetChoice(gridt_fillcolor_choice_item, g[gno].gridt[c].p.fillcol);
	SetChoice(gridt_linew_choice_item, g[gno].gridt[c].p.linew - 1);
	SetChoice(gridt_lines_choice_item, g[gno].gridt[c].p.lines - 1);
    }
}

/*
 * Create the frame Widget and the frame Widget
 */
void create_gridt_frame(void)
{
    extern Widget app_shell;
    Widget wbut, lab, rc, rc2, fr;

    if (gridt_frame) {
	update_gridt_items(cg);
	XtRaise(gridt_frame);
	return;
    }
    gridt_frame = XmCreateDialogShell(app_shell, "Grid setup", NULL, 0);
    handle_close(gridt_frame);
    gridt_panel = XmCreateRowColumn(gridt_frame, "panel", NULL, 0);
    XtVaSetValues(gridt_panel, XmNorientation, XmHORIZONTAL, NULL);

    rc2 = XmCreateRowColumn(gridt_panel, "rc2", NULL, 0);
    gridt_choice_item = CreatePanelChoice1(rc2, "Apply to time dependent grid:", 6, "1", "2", "3", "4", "5", NULL, NULL);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    lab = XmCreateLabel(rc, "Display:", NULL, 0);
    XtManageChild(lab);
    gridt_display_item = XmCreateToggleButton(rc, "Grid", NULL, 0);
    gridt_bath_item = XmCreateToggleButton(rc, "Bathymetry", NULL, 0);
    gridt_nodes_item = XmCreateToggleButton(rc, "Grid node numbers", NULL, 0);
    gridt_elements_item = XmCreateToggleButton(rc, "Grid element numbers", NULL, 0);
    gridt_depths_item = XmCreateToggleButton(rc, "Grid nodal depths", NULL, 0);
    gridt_filled_item = XmCreateToggleButton(rc, "Grid filled", NULL, 0);
    XtManageChild(gridt_display_item);
    XtManageChild(gridt_bath_item);
    XtManageChild(gridt_nodes_item);
    XtManageChild(gridt_elements_item);
    XtManageChild(gridt_depths_item);
    XtManageChild(gridt_filled_item);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(gridt_panel, "rc2", NULL, 0);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    gridt_color_choice_item = CreateColorChoice(rc, "Grid color:", 1);
    gridt_fillcolor_choice_item = CreateColorChoice(rc, "Grid fill color:", 1);
    gridt_linew_choice_item = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    gridt_lines_choice_item = CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(rc2, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) gridt_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Isolines...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) gridt_isolines_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Files...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gridtio_frame, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) gridt_Done_notify_proc, 0);
    XtManageChild(rc);
    XtManageChild(rc2);

    update_gridt_items(cg);
    XtManageChild(gridt_panel);
    XtRaise(gridt_frame);
}

/*
 * Notify and event procs
 */

static void gridt_Done_notify_proc(void)
{
    XtUnmanageChild(gridt_frame);
}

static void gridt_isolines_notify_proc(void)
{
    int a = GetChoice(gridt_choice_item);
    create_isolines_popup("Time dependent grids bathymetry", &g[cg].gridt[a].ip, 0);
}

static void gridt_define_notify_proc(void)
{
    int i, ming, maxg;
    int a;
    int c;
    a = GetChoice(gridt_choice_item);
    if (a == 0) {
	ming = maxg = cg;
    } else {
	ming = 0;
	g[cg].curgridt = maxg = maxgraph - 1;
    }
    if (a < MAXGRIDS) {
	g[cg].curgridt = a;
	ming = maxg = cg;
    } else if (a == MAXGRIDS) {
	ming = maxg = cg;
    } else {
	ming = 0;
	maxg = maxgraph - 1;
    }
    for (i = ming; i <= maxg; i++) {
	if (isactive_graph(i)) {
	    c = g[i].curgridt;
	    g[i].gridt[c].display = XmToggleButtonGetState(gridt_display_item) ? ON : OFF;
	    g[i].gridt[c].display_bath = XmToggleButtonGetState(gridt_bath_item) ? ON : OFF;
	    g[i].gridt[c].display_nodes = XmToggleButtonGetState(gridt_nodes_item) ? ON : OFF;
	    g[i].gridt[c].display_elements = XmToggleButtonGetState(gridt_elements_item) ? ON : OFF;
	    g[i].gridt[c].display_depths = XmToggleButtonGetState(gridt_depths_item) ? ON : OFF;
	    g[i].gridt[c].display_gridf = XmToggleButtonGetState(gridt_filled_item) ? ON : OFF;
	    g[i].gridt[c].p.color = GetChoice(gridt_color_choice_item);
	    g[i].gridt[c].p.fillcol = GetChoice(gridt_fillcolor_choice_item);
	    g[i].gridt[c].p.linew = GetChoice(gridt_linew_choice_item) + 1;
	    g[i].gridt[c].p.lines = GetChoice(gridt_lines_choice_item) + 1;
	}
    }
}
