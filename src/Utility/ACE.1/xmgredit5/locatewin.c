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
 * Locator Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: locatewin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;

static Widget locator_frame = (Widget) 0;
static Widget locator_panel;

/*
 * Panel item declarations
 */
static Widget *locator_onoff_item;
static Widget *delta_item;
static Widget *loc_formatx;
static Widget *loc_formaty;
static Widget *loc_precx;
static Widget *loc_precy;
static Widget locx_item;
static Widget locy_item;
static Widget *fixedp_item;

/*
 * Event and Notify proc declarations
 */
static void locator_define_notify_proc(void);
static void locator_reset_notify_proc(void);

extern int go_locateflag, deltaflag;
extern double dsx, dsy;
int pointset, locfx = 2, locfy = 2, locpx = 6, locpy = 6;

static char buf[128];

void update_locator_items(void)
{
    if (locator_frame) {
	SetChoice(locator_onoff_item, go_locateflag == FALSE);
	SetChoice(fixedp_item, pointset == TRUE);
	SetChoice(delta_item, deltaflag);
	SetChoice(loc_formatx, locfx);
	SetChoice(loc_formaty, locfy);
	SetChoice(loc_precx, locpx);
	SetChoice(loc_precy, locpy);
	if (pointset) {
	    sprintf(buf, "%lf", dsx);
	    panel_setstr_value(locx_item, buf);
	    sprintf(buf, "%lf", dsy);
	    panel_setstr_value(locy_item, buf);
	}
    }
}

/*
 * Create the locator Panel
 */
void create_locator_frame(void)
{
    Widget wbut, rc, fr, rc2;
    extern Widget app_shell;
    int x, y;

    if (locator_frame) {
	update_locator_items();
	XtRaise(locator_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    locator_frame = XmCreateDialogShell(app_shell, "locator", NULL, 0);
    XtVaSetValues(locator_frame, XmNx, x, XmNy, y, NULL);
    locator_panel = XmCreateRowColumn(locator_frame, "ticks_rc", NULL, 0);

    locator_onoff_item = (Widget *) CreatePanelChoice1(locator_panel, "Locator:",
						       3,
						       "ON",
						       "OFF",
						       NULL,
						       NULL);
    delta_item = (Widget *) CreatePanelChoice1(locator_panel, "Locator display type:",
					       7,
					       "[X, Y]",
					       "[DX, DY]",
					       "[DISTANCE]",
					       "[R, Theta]",
					       "[VX, VY]",
					       "[SX, SY]",
					       NULL,
					       NULL);
    fixedp_item = CreatePanelChoice1(locator_panel, "Fixed point:",
				     3, "OFF", "ON", NULL,
				     NULL);

    rc2 = XmCreateRowColumn(locator_panel, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);

    loc_formatx = CreatePanelChoice1(rc, "Format X:",
				     4,
				     "Decimal",
				     "Exponential",
				     "General",
				     NULL,
				     NULL);
    loc_precx = CreatePanelChoice2(rc, "Precision X:",
				   4, 12,
		     "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
				   NULL,
				   NULL);
    locx_item = (Widget) CreateTextItem2(rc, 10, "Fixed point X:");
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    loc_formaty = CreatePanelChoice1(rc, "Format Y:",
				     4,
				     "Decimal",
				     "Exponential",
				     "General",
				     NULL,
				     NULL);

    loc_precy = CreatePanelChoice2(rc, "Precision Y:",
				   4, 12,
		     "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
				   NULL,
				   NULL);
    locy_item = (Widget) CreateTextItem2(rc, 10, "Fixed point Y:");
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, locator_panel, NULL);
    rc = XmCreateRowColumn(locator_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) locator_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Reset", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) locator_reset_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, locator_frame);

    XtManageChild(rc);

    update_locator_items();
    XtManageChild(locator_panel);
    XtManageChild(locator_frame);
}				/* end create_locator_panel */

static void locator_define_notify_proc(void)
{
    double atof(const char *);
    int type;

    go_locateflag = (int) GetChoice(locator_onoff_item) == 0;
    deltaflag = (int) GetChoice(delta_item);
    locfx = (int) GetChoice(loc_formatx);
    locfy = (int) GetChoice(loc_formaty);
    locpx = (int) GetChoice(loc_precx);
    locpy = (int) GetChoice(loc_precy);
    pointset = (int) GetChoice(fixedp_item);
    if (pointset) {
	strcpy(buf, (char *) xv_getstr(locx_item));
	if (buf[0]) {
	    dsx = atof(buf);
	}
	strcpy(buf, (char *) xv_getstr(locy_item));
	if (buf[0]) {
	    dsy = atof(buf);
	}
    }
    XtUnmanageChild(locator_frame);
}

static void locator_reset_notify_proc(void)
{
    dsx = dsy = 0.0;
    pointset = FALSE;
    deltaflag = 0;
    locfx = 2;
    locfy = 2;
    locpx = 6;
    locpy = 6;
    update_locator_items();
}
