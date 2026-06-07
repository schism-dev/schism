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
 * Grid popups
 *
 */

#ifndef lint
static char RCSid[] = "$Id: gridwin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

/*
 * Interpolate bathymetry
 */

/* 0 = no warn,
 * 1 = report on all nodes,
 * 2 = report total number of nodes only
 */
int interp_warn_level = 0;
int interp_in_region = 0;
int noskip = 0;

void do_interp_depths(void);

Widget interp_frame;
Widget *interp_warnlevel_item;
Widget interp_skip_item;
Widget interp_region_item;

void update_interp_items(void)
{
    if (interp_frame) {
	XmToggleButtonSetState(interp_skip_item, noskip, False);
	XmToggleButtonSetState(interp_region_item, interp_in_region, False);
	SetChoice(interp_warnlevel_item, interp_warn_level);
    }
}

void accept_interp_items(void)
{
    noskip = XmToggleButtonGetState(interp_skip_item);
    interp_in_region = XmToggleButtonGetState(interp_region_item);
    interp_warn_level = GetChoice(interp_warnlevel_item);
    do_interp_depths();
}

void create_interp_frame(void)
{
    extern Widget app_shell;
    Widget panel, wbut, lab, rc, rc2, fr;

    if (interp_frame) {
	update_interp_items();
	XtRaise(interp_frame);
	return;
    }
    interp_frame = XmCreateDialogShell(app_shell, "Grid setup", NULL, 0);
    handle_close(interp_frame);
    panel = XmCreateRowColumn(interp_frame, "panel", NULL, 0);

    interp_skip_item = XmCreateToggleButton(panel, "Interpolate depths to all nodes", NULL, 0);
    XtManageChild(interp_skip_item);
    interp_region_item = XmCreateToggleButton(panel, "Restrict interpolation to region", NULL, 0);
    XtManageChild(interp_region_item);
    interp_warnlevel_item = CreatePanelChoice1(panel, "Warn level:",
					       4,
					       "No warnings",
				   "All nodes not found in background grid",
		       "Total number of nodes not found in background grid",
					       NULL,
					       NULL);
    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_interp_items, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, interp_frame);
    XtManageChild(rc);

    update_interp_items();
    XtManageChild(panel);
    XtRaise(interp_frame);
}
