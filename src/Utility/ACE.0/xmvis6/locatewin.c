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
 * Locator Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: locatewin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Frame.h>
#include <Xm/LabelG.h>
#include <Xm/PushBG.h>
#include <Xm/RowColumn.h>
#include <Xm/SeparatoG.h>

#include "defines.h"
#include "globals.h"
#include "motifinc.h"

extern Widget app_shell;

static Widget locator_frame;
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
static int locator_Done_notify_proc(void);
static int locator_define_notify_proc(void);
static int locator_reset_notify_proc(void);

extern int go_locateflag;
static int locfx = 2, locfy = 2, locpx = 6, locpy = 6;

/*
 * set action to SEL_POINT for selecting the locator reference point
 */
void do_select_point(void)
{
    set_action(0);
    set_action(SEL_POINT);
    g[cg].pointset = TRUE;
}

/*
 * clear the locator reference point
 */
void do_clear_point(void)
{
    g[cg].pointset = FALSE;
    g[cg].pt_type = 0;
    g[cg].dsx = g[cg].dsy = 0.0;
}

void update_locator_items(int gno)
{
    if (locator_frame) {
	SetChoice(locator_onoff_item, go_locateflag == FALSE);
	SetChoice(fixedp_item, g[gno].pointset == TRUE);
	SetChoice(delta_item, g[gno].pt_type);
	SetChoice(loc_formatx, get_format_index(g[gno].fx));
	SetChoice(loc_formaty, get_format_index(g[gno].fy));
	SetChoice(loc_precx, g[gno].px);
	SetChoice(loc_precy, g[gno].py);
	if (g[gno].pointset) {
	    sprintf(buf, "%lf", g[gno].dsx);
	    xv_setstr(locx_item, buf);
	    sprintf(buf, "%lf", g[gno].dsy);
	    xv_setstr(locy_item, buf);
	}
    }
}

/*
 * Create the locator Panel
 */
void create_locator_frame(void)
{
    Widget wbut, rc, fr, rc2;
    int x, y;

    if (locator_frame) {
	update_locator_items(cg);
	XtRaise(locator_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    locator_frame = XmCreateDialogShell(app_shell, "locator", NULL, 0);
    XtVaSetValues(locator_frame, XmNx, x, XmNy, y, NULL);
    locator_panel = XmCreateRowColumn(locator_frame, "ticks_rc", NULL, 0);

    locator_onoff_item = (Widget *) CreatePanelChoice1(locator_panel, "Locator:", 3, "ON", "OFF", NULL, NULL);
    delta_item = (Widget *) CreatePanelChoice1(locator_panel, "Locator display type:", 7, "[X, Y]", "[DX, DY]", "[DISTANCE]", "[R, Theta]", "[VX, VY]", "[SX, SY]", NULL, NULL);
    fixedp_item = CreatePanelChoice1(locator_panel, "Fixed point:", 3, "OFF", "ON", NULL, NULL);

    rc2 = XmCreateRowColumn(locator_panel, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);

    loc_formatx = CreatePanelChoice2(rc,
				     "Format X:",
				     4, 27,
				     "Decimal",
				     "Exponential",
				     "Power (decimal)",
				     "General",
				     "DD-MM-YY",
				     "MM-DD-YY",
				     "MM-YY",
				     "MM-DD",
				     "Month-DD",
				     "DD-Month",
				     "Month (abrev.)",
				     "Month",
				     "Day of week (abrev.)",
				     "Day of week",
				     "Day of year",
				     "HH:MM:SS.s",
				     "MM-DD HH:MM:SS.s", "MM-DD-YY HH:MM:SS.s", "Degrees (lon)", "DD MM' (lon)", "DD MM' SS.s\" (lon)", "MM' SS.s\" (lon)", "Degrees (lat)", "DD MM' (lat)", "DD MM' SS.s\" (lat)", "MM' SS.s\" (lat)", NULL, NULL);
    loc_precx = CreatePanelChoice1(rc, "Precision X:", 12, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", NULL, NULL);
    locx_item = (Widget) CreateTextItem2(rc, 10, "Fixed point X:");
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc2, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    loc_formaty = CreatePanelChoice2(rc,
				     "Format Y:",
				     4, 27,
				     "Decimal",
				     "Exponential",
				     "Power (decimal)",
				     "General",
				     "DD-MM-YY",
				     "MM-DD-YY",
				     "MM-YY",
				     "MM-DD",
				     "Month-DD",
				     "DD-Month",
				     "Month (abrev.)",
				     "Month",
				     "Day of week (abrev.)",
				     "Day of week",
				     "Day of year",
				     "HH:MM:SS.s",
				     "MM-DD HH:MM:SS.s", "MM-DD-YY HH:MM:SS.s", "Degrees (lon)", "DD MM' (lon)", "DD MM' SS.s\" (lon)", "MM' SS.s\" (lon)", "Degrees (lat)", "DD MM' (lat)", "DD MM' SS.s\" (lat)", "MM' SS.s\" (lat)", NULL, NULL);

    loc_precy = CreatePanelChoice1(rc, "Precision Y:", 12, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", NULL, NULL);
    locy_item = (Widget) CreateTextItem2(rc, 10, "Fixed point Y:");
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(rc2);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, locator_panel, NULL);
    rc = XmCreateRowColumn(locator_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) locator_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Reset", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) locator_reset_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) locator_Done_notify_proc, 0);

    XtManageChild(rc);

    update_locator_items(cg);
    XtManageChild(locator_panel);
    XtRaise(locator_frame);
}				/* end create_locator_panel */

/*
 * Notify and event procs
 */

static int locator_Done_notify_proc(void)
{
    XtUnmanageChild(locator_frame);
}

static int locator_define_notify_proc(void)
{
    double atof(const char *);
    int type;

    go_locateflag = (int) GetChoice(locator_onoff_item) == 0;
    type = g[cg].pt_type = (int) GetChoice(delta_item);
    locfx = g[cg].fx = format_types[(int) GetChoice(loc_formatx)];
    locfy = g[cg].fy = format_types[(int) GetChoice(loc_formaty)];
    locpx = g[cg].px = (int) GetChoice(loc_precx);
    locpy = g[cg].py = (int) GetChoice(loc_precy);
    g[cg].pointset = (int) GetChoice(fixedp_item);
    if (g[cg].pointset) {
	strcpy(buf, (char *) xv_getstr(locx_item));
	if (buf[0]) {
	    g[cg].dsx = atof(buf);
	}
	strcpy(buf, (char *) xv_getstr(locy_item));
	if (buf[0]) {
	    g[cg].dsy = atof(buf);
	}
    }
    make_format(cg);
    XtUnmanageChild(locator_frame);
}

static int locator_reset_notify_proc(void)
{
    g[cg].dsx = g[cg].dsy = 0.0;	/* locator props */
    g[cg].pointset = FALSE;
    g[cg].pt_type = 0;
    g[cg].fx = GENERAL;
    g[cg].fy = GENERAL;
    g[cg].px = 6;
    g[cg].py = 6;
    update_locator_items(cg);
}

XmString astring, pstring;
Widget arealab, perimlab;

void do_select_area(void);
void do_select_peri(void);

void create_area_frame(void)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;
    extern Widget app_shell;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Area/perimeter", NULL, 0);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

    arealab = XtVaCreateManagedWidget("label Area", xmLabelGadgetClass, dialog, XmNlabelString, astring = XmStringCreateLtoR("[    Area    ]", charset), NULL);

    perimlab = XtVaCreateManagedWidget("label Perim", xmLabelGadgetClass, dialog, XmNlabelString, pstring = XmStringCreateLtoR("[    Perim    ]", charset), NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Area", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_select_area, 0);

    wbut = XtVaCreateManagedWidget("Perimeter", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_select_peri, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}
