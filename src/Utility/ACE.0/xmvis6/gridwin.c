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
static char RCSid[] = "$Id: gridwin.c,v 1.3 2005/04/07 14:17:02 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget grid_frame;
static Widget grid_panel;

void create_grid_frame(void);

/*
 * Widget item declarations
 */
static Widget *grid_choice_item;
static Widget grid_display_item;
static Widget grid_bath_item;
static Widget grid_bound_item;
static Widget grid_nodes_item;
static Widget grid_elements_item;
static Widget grid_depths_item;
static Widget grid_cour_item;
static Widget grid_courn_item;
/*
static Widget grid_dimw_item;
*/
static Widget grid_filled_item;
static Widget *grid_color_choice_item;
static Widget *grid_fillcolor_choice_item;
static Widget *grid_lines_choice_item;
static Widget *grid_linew_choice_item;
static Widget *grid_bcolor_choice_item;
static Widget *grid_blines_choice_item;
static Widget *grid_blinew_choice_item;
static Widget *grid_fillbg_choice_item;
static Widget *grid_bgcolor_choice_item;
static Widget *grid_applyto_choice_item;

/*
 * Event and Notify proc declarations
 */
static Widget gridio_frame;
static Widget *gridio_item;
static Widget gridio_swap_item;

/*
 * Event and Notify proc declarations
 */
void gridio_done_proc(void);
void gridio_accept_proc(void);
static void grid_Done_notify_proc(void);
static void grid_isolines_notify_proc(void);
static void grid_define_notify_proc(void);

void update_grid_items(int gno);

/*
 * Create the grid Frame and the grid Panel
 */
extern Widget app_shell;

void set_curgrid(Widget w, int cd)
{
    g[cg].curgrid = cd;
    update_grid_items(cg);
}

void create_gridio_frame(void)
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
	    XtAddCallback(gridio_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curgrid, (XtPointer) i);
	}
	gridio_swap_item = XtVaCreateManagedWidget("Swap bytes", xmToggleButtonWidgetClass, rc, NULL);

	XtManageChild(rc);
	XtAddCallback(gridio_frame, XmNcancelCallback, (XtCallbackProc) gridio_done_proc, NULL);
	XtAddCallback(gridio_frame, XmNokCallback, (XtCallbackProc) gridio_accept_proc, NULL);
    }
    XtRaise(gridio_frame);
    XmToggleButtonSetState(gridio_swap_item, swapBytes, False);
    update_grid_items(cg);
}

void gridio_done_proc(void)
{
    XtUnmanageChild(gridio_frame);
}

void gridio_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int gridno;
    Widget textw;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(gridio_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    swapBytes = XmToggleButtonGetState(gridio_swap_item);

    gridno = GetChoice(gridio_item);
    if (grid[gridno].active == ON) {
	if (!yesno("Grid is active, kill it ?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(gridio_frame);
    if (!readgrid(gridno, s)) {
	sprintf(buf, "Error reading file %s", s);
	errwin(buf);
    } else {
	XtUnmanageChild(gridio_frame);
	compute_boundary(gridno);
	create_grid_frame();
    }
    unset_wait_cursor(gridio_frame);
}

void update_grid_items(int gno)
{
    int gd;
    int c = g[gno].curgrid;
    if (grid_frame) {
	SetChoice(grid_choice_item, c);
	XmToggleButtonSetState(grid_display_item, g[gno].grid[c].display == ON, False);
	XmToggleButtonSetState(grid_bath_item, g[gno].grid[c].display_bath == ON, False);
	XmToggleButtonSetState(grid_bound_item, g[gno].grid[c].display_boundary == ON, False);
	XmToggleButtonSetState(grid_nodes_item, g[gno].grid[c].display_nodes == ON, False);
	XmToggleButtonSetState(grid_elements_item, g[gno].grid[c].display_elements == ON, False);
	XmToggleButtonSetState(grid_depths_item, g[gno].grid[c].display_depths == ON, False);
	XmToggleButtonSetState(grid_filled_item, g[gno].grid[c].display_gridf == ON, False);
	XmToggleButtonSetState(grid_cour_item,
			       g[gno].grid[gno].display_courant == ON, False);
	XmToggleButtonSetState(grid_courn_item,
			       g[gno].grid[gno].display_courantn == ON, False);
/*
	XmToggleButtonSetState(grid_dimw_item,
			       g[gno].grid[curgrid].display_dimw == ON, False);
*/
	SetChoice(grid_color_choice_item, g[gno].grid[c].p.color);
	SetChoice(grid_fillcolor_choice_item, g[gno].grid[c].p.fillcol);
	SetChoice(grid_linew_choice_item, g[gno].grid[c].p.linew - 1);
	SetChoice(grid_lines_choice_item, g[gno].grid[c].p.lines - 1);
	SetChoice(grid_bcolor_choice_item, g[gno].grid[c].bp.color);
	SetChoice(grid_blinew_choice_item, g[gno].grid[c].bp.linew - 1);
	SetChoice(grid_blines_choice_item, g[gno].grid[c].bp.lines - 1);
	SetChoice(grid_fillbg_choice_item, g[gno].grid[c].bp.fill == ON);
	SetChoice(grid_bgcolor_choice_item, g[gno].grid[c].bp.fillcol == ON);
    }
}

void compute_grid_boundary(void)
{
    int gridno = GetChoice(grid_choice_item);
    if (grid[gridno].active == ON) {
	compute_boundary(gridno);
    }
}

/*
 * Create the frame Widget and the frame Widget
 */
void create_grid_frame(void)
{
    extern Widget app_shell;
    Widget wbut, lab, rc, rc2, fr;

    if (grid_frame) {
	update_grid_items(cg);
	XtRaise(grid_frame);
	return;
    }
    grid_frame = XmCreateDialogShell(app_shell, "Grid setup", NULL, 0);
    handle_close(grid_frame);
    grid_panel = XmCreateRowColumn(grid_frame, "panel", NULL, 0);
    XtVaSetValues(grid_panel, XmNorientation, XmHORIZONTAL, NULL);

    rc2 = XmCreateRowColumn(grid_panel, "rc2", NULL, 0);
    grid_choice_item = CreatePanelChoice1(rc2, "Apply to grid:", 6, "1", "2", "3", "4", "5", NULL, NULL);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    lab = XmCreateLabel(rc, "Display:", NULL, 0);
    XtManageChild(lab);
    grid_display_item = XmCreateToggleButton(rc, "Grid", NULL, 0);
    grid_bath_item = XmCreateToggleButton(rc, "Bathymetry", NULL, 0);
    grid_bound_item = XmCreateToggleButton(rc, "Grid boundary", NULL, 0);
    grid_nodes_item = XmCreateToggleButton(rc, "Grid node numbers", NULL, 0);
    grid_elements_item = XmCreateToggleButton(rc, "Grid element numbers", NULL, 0);
    grid_depths_item = XmCreateToggleButton(rc, "Grid nodal depths", NULL, 0);
    grid_filled_item = XmCreateToggleButton(rc, "Grid filled", NULL, 0);
    grid_cour_item = XmCreateToggleButton(rc, "Grid Courant graph", NULL, 0);
    grid_courn_item = XmCreateToggleButton(rc, "Grid Courant numbers", NULL, 0);
/*
    grid_dimw_item = XmCreateToggleButton(rc, "Grid dimensionless wavelength", NULL, 0);
*/
    XtManageChild(grid_display_item);
    XtManageChild(grid_bath_item);
    XtManageChild(grid_bound_item);
    XtManageChild(grid_nodes_item);
    XtManageChild(grid_elements_item);
    XtManageChild(grid_depths_item);
    XtManageChild(grid_cour_item);
    XtManageChild(grid_courn_item);
/*
    XtManageChild(grid_dimw_item);
*/
    XtManageChild(grid_filled_item);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(grid_panel, "rc2", NULL, 0);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    grid_color_choice_item = CreateColorChoice(rc, "Grid color:", 1);
    grid_fillcolor_choice_item = CreateColorChoice(rc, "Grid fill color:", 1);
    grid_linew_choice_item = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    grid_lines_choice_item = CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    grid_bcolor_choice_item = CreateColorChoice(rc, "Grid boundary color:", 1);
    grid_blinew_choice_item = CreatePanelChoice1(rc, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    grid_blines_choice_item = CreatePanelChoice1(rc, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    grid_fillbg_choice_item = CreatePanelChoice1(rc, "Fill for islands:", 3, "None", "Filled", NULL, NULL);

    grid_bgcolor_choice_item = CreateColorChoice(rc, "Island color:", 1);
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(rc2, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) grid_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Isolines of bathymetry...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) grid_isolines_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Compute boundary", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) compute_grid_boundary, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) grid_Done_notify_proc, 0);
    XtManageChild(rc);
    XtManageChild(rc2);

    update_grid_items(cg);
    XtManageChild(grid_panel);
    XtRaise(grid_frame);
}

/*
 * Notify and event procs
 */

static void grid_Done_notify_proc(void)
{
    XtUnmanageChild(grid_frame);
}

static void grid_isolines_notify_proc(void)
{
    int a = GetChoice(grid_choice_item);
    g[cg].grid[a].ip.cmin = grid[a].dmin;
    g[cg].grid[a].ip.cmax = grid[a].dmax;
    create_isolines_popup("Grid bathymetry", &g[cg].grid[a].ip, 0);
}

static void grid_define_notify_proc(void)
{
    int i, ming, maxg;
    int a;
    int c;
    a = GetChoice(grid_choice_item);
    if (a == 0) {
	ming = maxg = cg;
    } else {
	ming = 0;
	g[cg].curgrid = maxg = maxgraph - 1;
    }
    if (a < MAXGRIDS) {
	g[cg].curgrid = a;
	ming = maxg = cg;
    } else if (a == MAXGRIDS) {
	ming = maxg = cg;
    } else {
	ming = 0;
	maxg = maxgraph - 1;
    }
    for (i = ming; i <= maxg; i++) {
	if (isactive_graph(i)) {
	    c = g[i].curgrid;
	    g[i].grid[c].display = XmToggleButtonGetState(grid_display_item) ? ON : OFF;
	    g[i].grid[c].display_bath = XmToggleButtonGetState(grid_bath_item) ? ON : OFF;
	    g[i].grid[c].display_boundary = XmToggleButtonGetState(grid_bound_item) ? ON : OFF;
	    g[i].grid[c].display_nodes = XmToggleButtonGetState(grid_nodes_item) ? ON : OFF;
	    g[i].grid[c].display_elements = XmToggleButtonGetState(grid_elements_item) ? ON : OFF;
	    g[i].grid[c].display_depths = XmToggleButtonGetState(grid_depths_item) ? ON : OFF;
	    g[i].grid[c].display_courant =
		    XmToggleButtonGetState(grid_cour_item) ? ON : OFF;
	    g[i].grid[c].display_courantn =
		    XmToggleButtonGetState(grid_courn_item) ? ON : OFF;
/*
	    g[i].grid[c].display_dimw =
		    XmToggleButtonGetState(grid_dimw_item) ? ON : OFF;
*/
	    g[i].grid[c].display_gridf = XmToggleButtonGetState(grid_filled_item) ? ON : OFF;
	    g[i].grid[c].p.color = GetChoice(grid_color_choice_item);
	    g[i].grid[c].p.fillcol = GetChoice(grid_fillcolor_choice_item);
	    g[i].grid[c].p.linew = GetChoice(grid_linew_choice_item) + 1;
	    g[i].grid[c].p.lines = GetChoice(grid_lines_choice_item) + 1;
	    g[i].grid[c].bp.color = GetChoice(grid_bcolor_choice_item);
	    g[i].grid[c].bp.linew = GetChoice(grid_blinew_choice_item) + 1;
	    g[i].grid[c].bp.lines = GetChoice(grid_blines_choice_item) + 1;
	    g[i].grid[c].bp.fill = GetChoice(grid_fillbg_choice_item) ? ON : OFF;
	    g[i].grid[c].bp.fillcol = GetChoice(grid_bgcolor_choice_item);
	}
    }
}
