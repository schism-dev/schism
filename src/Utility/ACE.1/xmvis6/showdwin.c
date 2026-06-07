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
 * graphics set routines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: showdwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
extern Widget app_shell;

static Widget panel;
static Widget gridtoggle_button[MAXGRIDS];
static Widget gridbathtoggle_button[MAXGRIDS];
static Widget gridboundtoggle_button[MAXGRIDS];
static Widget drogstoggle_button[MAXPATHLINES];
static Widget streamstoggle_button[MAXPATHLINES];
static Widget numberstoggle_button[MAXPATHLINES];
static Widget elaconc_toggle_item[MAXELA];
static Widget flow_toggle_item[MAXTEANL];
static Widget flow_type_item[MAXTEANL];
static Widget elev_toggle_item[MAXPATHLINES];
static Widget elevmarker_toggle_item;
static Widget adcircflow_toggle_item[MAXADCIRC];
static Widget adcircflow_type_item[MAXADCIRC];
static Widget adcircelev_toggle_item[MAXADCIRC];
static Widget adcircelevmarker_toggle_item;
static Widget conc_toggle_item;
static Widget grid_toggle_item;

void setdrognumbs();

Widget display_frame;
Widget display_panel;

#define NITEMS 11
#define NAUTOITEMS 5

Widget toggle_display_item[NITEMS];
Widget toggle_auto_item[NAUTOITEMS];
extern Widget main_display_item[];

void update_display(void);

void display_done_proc(void)
{
    XtUnmanageChild(display_frame);
}

/*
 * If reading in a parameter file, set the display of the buttons on the left
 */
void set_display_items(void)
{
    int i;
    int c = g[cg].curgrid;
    g[cg].grid[c].display_flags[EDIT_GRID] = (g[cg].grid[c].display == ON);
    g[cg].grid[c].display_flags[EDIT_BOUNDARY] = (g[cg].grid[c].display_boundary == ON);
    g[cg].grid[c].display_flags[EDIT_GRID_ISOLINES] = (g[cg].grid[c].display_bath == ON);
    g[cg].grid[c].display_flags[EDIT_GRID_NODE_NUMBERS] = (g[cg].grid[c].display_nodes == ON);
    g[cg].grid[c].display_flags[EDIT_GRID_ELEMENT_NUMBERS] = (g[cg].grid[c].display_elements == ON);
    g[cg].grid[c].display_flags[EDIT_GRID_FILLED] = (g[cg].grid[c].display_gridf == ON);
    g[cg].grid[c].display_flags[EDIT_GRID_DEPTHS] = (g[cg].grid[c].display_depths == ON);
}

void update_display_items(void)
{
    int i, j, c = g[cg].curgrid;
    for (i = 0; i < ndisplay; i++) {
	if (main_display_item[i]) {
	    XmToggleButtonSetState(main_display_item[i], g[cg].grid[c].display_flags[i], False);
	}
	switch (i) {
	case EDIT_GRID:
	    g[cg].grid[c].display = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_BOUNDARY:
	    g[cg].grid[c].display_boundary = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_GRID_ISOLINES:
	    g[cg].grid[c].display_bath = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_GRID_NODE_NUMBERS:
	    g[cg].grid[c].display_nodes = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_GRID_ELEMENT_NUMBERS:
	    g[cg].grid[c].display_elements = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_GRID_FILLED:
	    g[cg].grid[c].display_gridf = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	case EDIT_GRID_DEPTHS:
	    g[cg].grid[c].display_depths = g[cg].grid[c].display_flags[i] ? ON : OFF;
	    break;
	}
    }
    update_display();
    update_grid_items(cg);
}

void get_display_toggles(void)
{
    int i;
    int c = g[cg].curgrid;
    for (i = 0; i < ndisplay; i++) {
	g[cg].grid[c].display_flags[i] = XmToggleButtonGetState(toggle_display_item[i]);
    }
    for (i = 0; i < ndisplay; i++) {
	XmToggleButtonSetState(main_display_item[i], g[cg].grid[c].display_flags[i], False);
    }
    update_display();
    update_grid_items(cg);
}

void update_display(void)
{
    int i, c = g[cg].curgrid;
    if (display_frame) {
	for (i = 0; i < MAXGRIDS; i++) {
	    XmToggleButtonSetState(gridtoggle_button[i], g[cg].grid[i].display == ON, False);
	    XmToggleButtonSetState(gridbathtoggle_button[i], g[cg].grid[i].display_bath == ON, False);
	    XmToggleButtonSetState(gridboundtoggle_button[i], g[cg].grid[i].display_boundary == ON, False);
	}
	for (i = 0; i < MAXPATHLINES; i++) {
	    XmToggleButtonSetState(drogstoggle_button[i], g[cg].drogues[i].display == ON, False);
	    XmToggleButtonSetState(drogstoggle_button[i], g[cg].drogues[i].display_streaml == ON, False);
	    XmToggleButtonSetState(drogstoggle_button[i], g[cg].drogues[i].display_type == ON, False);
	}
	for (i = 0; i < MAXTEANL; i++) {
	    XmToggleButtonSetState(flow_toggle_item[i], g[cg].flowf[i].display == NODES, False);
	    XmToggleButtonSetState(flow_type_item[i], g[cg].flowf[i].display == CENTER, False);
	    XmToggleButtonSetState(elev_toggle_item[i], g[cg].flowf[i].display_elev == ON, False);
	}
/*

	XmToggleButtonSetState(flow_toggle_item, show_flow, False);
	XmToggleButtonSetState(flow_type_item, show_flow_center, False);
	XmToggleButtonSetState(elev_toggle_item, show_teanlelev, False);
	XmToggleButtonSetState(elevmarker_toggle_item, show_elev, False);

	XmToggleButtonSetState(adcircflow_toggle_item, show_adcirc, False);
	XmToggleButtonSetState(adcircflow_type_item, show_adcircflow_center, False);
	XmToggleButtonSetState(adcircelev_toggle_item, show_adcircelev, False);
	XmToggleButtonSetState(adcircelevmarker_toggle_item, show_adcircelevmarkers, False);

	XmToggleButtonSetState(drogstoggle_button, show_drogues, False);
	XmToggleButtonSetState(numberstoggle_button, show_drognumbs, False);

	XmToggleButtonSetState(elaconc_toggle_item, show_ela, False);

*/
    }
}

void accept_display(void)
{
    int i;

    for (i = 0; i < MAXGRIDS; i++) {
	g[cg].grid[i].display = XmToggleButtonGetState(gridtoggle_button[i]) ? ON : OFF;
	g[cg].grid[i].display_bath = XmToggleButtonGetState(gridbathtoggle_button[i]) ? ON : OFF;
	g[cg].grid[i].display_boundary = XmToggleButtonGetState(gridboundtoggle_button[i]) ? ON : OFF;
    }
/*
    for (i = 0; i < MAXPATHLINES; i++) {
	g[cg].drogues[i].display =
		XmToggleButtonGetState(drogstoggle_button[i]) ? ON : OFF;
	g[cg].drogues[i].display_streaml =
		XmToggleButtonGetState(drogstoggle_button[i]) ? ON : OFF;
	g[cg].drogues[i].display_type =
		XmToggleButtonGetState(drogstoggle_button[i]) ? ON : OFF;
    }
    for (i = 0; i < MAXTEANL; i++) {
	g[cg].flowf[i].display =
		XmToggleButtonGetState(flow_toggle_item) ? NODES : CENTER;
	g[cg].flowf[i].display =
		XmToggleButtonGetState(flow_type_item) ? CENTER : NODES;
	g[cg].flowf[i].display_elev =
		XmToggleButtonGetState(elev_toggle_item) ? ON : OFF;
    }
    for (i = 0; i < MAXTEANL; i++) {
	g[cg].flowf[i].display_at = XmToggleButtonGetState(flow_type_item);
    }
    show_teanlelev = XmToggleButtonGetState(elev_toggle_item);
    show_elev = XmToggleButtonGetState(elevmarker_toggle_item);

    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].flowt[i].display = XmToggleButtonGetState(adcircflow_toggle_item);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].flowt[i].display_at = XmToggleButtonGetState(adcircflow_type_item);
    }
    show_adcircelev = XmToggleButtonGetState(adcircelev_toggle_item);
    show_adcircelevmarkers = XmToggleButtonGetState(adcircelevmarker_toggle_item);

    show_drogues = XmToggleButtonGetState(drogstoggle_button);
    show_drognumbs = XmToggleButtonGetState(numberstoggle_button);

    show_ela = XmToggleButtonGetState(elaconc_toggle_item);

*/
    update_teanl_flow();
    update_adcirc_flow();

    setredraw_world();
}

void create_display_popup(void)
{
    int i, ac;
    Widget bt, rc, display_panel, panel, sep, rc2, fr;
    Dimension ww, wh;
    Position x, y;
    Arg a[2];
    XmString str;
    char buf[128];
    if (display_frame) {
	update_display();
	XtRaise(display_frame);
	return;
    }
    display_frame = XmCreateDialogShell(app_shell, "Display", NULL, 0);
    display_panel = XtVaCreateManagedWidget("dispbb", xmRowColumnWidgetClass, display_frame, NULL);
    panel = XmCreateRowColumn(display_panel, "rc", NULL, 0);
    XtVaSetValues(panel, XmNorientation, XmHORIZONTAL, NULL);
    fr = XmCreateFrame(panel, "frame", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Grid:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXGRIDS; i++) {
	sprintf(buf, "%1d", i + 1);
	gridtoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Bathymetry:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXGRIDS; i++) {
	sprintf(buf, "%1d", i + 1);
	gridbathtoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Boundary:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXGRIDS; i++) {
	sprintf(buf, "%1d", i + 1);
	gridboundtoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);

    sep = XmCreateSeparator(rc, "sep", NULL, 0);
    XtManageChild(sep);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Drogues:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXPATHLINES; i++) {
	sprintf(buf, "%1d", i + 1);
	drogstoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Numbers:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXPATHLINES; i++) {
	sprintf(buf, "%1d", i + 1);
	numberstoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("Path lines:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXPATHLINES; i++) {
	sprintf(buf, "%1d", i + 1);
	streamstoggle_button[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    sep = XmCreateSeparator(rc, "sep", NULL, 0);
    XtManageChild(sep);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("TEA-NL flow at nodes:", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXTEANL; i++) {
	sprintf(buf, "%1d", i + 1);
	flow_toggle_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("TEA-NL flow at center", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXTEANL; i++) {
	sprintf(buf, "%1d", i + 1);
	flow_type_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("TEA-NL elevations", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXTEANL; i++) {
	sprintf(buf, "%1d", i + 1);
	elev_toggle_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    elevmarker_toggle_item = XtVaCreateManagedWidget("Elevation markers", xmToggleButtonWidgetClass, rc, NULL);
    sep = XmCreateSeparatorGadget(rc, "sep", NULL, 0);
    XtManageChild(sep);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("ADCIRC flows at nodes", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXADCIRC; i++) {
	sprintf(buf, "%1d", i + 1);
	adcircflow_toggle_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("ADCIRC flows at center", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXADCIRC; i++) {
	sprintf(buf, "%1d", i + 1);
	adcircflow_type_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("ADCIRC elevations", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXADCIRC; i++) {
	sprintf(buf, "%1d", i + 1);
	adcircelev_toggle_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    adcircelevmarker_toggle_item = XtVaCreateManagedWidget("Elevation markers", xmToggleButtonWidgetClass, rc, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(panel, "frame", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    XtVaCreateManagedWidget("ELA concentrations", xmLabelGadgetClass, rc2, NULL);
    for (i = 0; i < MAXELA; i++) {
	sprintf(buf, "%1d", i + 1);
	elaconc_toggle_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc2, NULL);
    }
    XtManageChild(rc2);
    XtManageChild(rc);
    XtManageChild(fr);

    sep = XmCreateSeparatorGadget(rc, "sep", NULL, 0);
    XtManageChild(sep);

    rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, XmNpacking, XmPACK_TIGHT, NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) accept_display, 0);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) display_done_proc, NULL);
    XtManageChild(rc2);

    XtManageChild(rc);
    XtManageChild(fr);

    update_display();
    XtManageChild(panel);
    XtManageChild(display_panel);
    XtRaise(display_frame);
}
