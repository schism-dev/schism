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
 * set the parameters for isolines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: placewin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern int running;

Widget place_frame;
Widget place_panel;
Widget place_items[10];

static int curitem = 0;

void define_place_proc(void)
{
    set_action(0);
    switch (curitem) {
    case 0:
	set_action(PLACE_CLOCK);
	break;
    case 1:
	set_action(PLACE_INFO);
	break;
    case 2:
	set_action(PLACE_TIMELINE);
	break;
    case 3:
	set_action(PLACE_VSCALE);
	break;
    case 4:
	set_action(PLACE_MAPSCALE);
	break;
    case 5:
	set_action(PLACE_TITLE);
	break;
    case 6:
	set_action(PLACE_BATH_LEGEND);
	break;
    case 7:
	set_action(PLACE_CONC_LEGEND);
	break;
    }
}

void set_type_proc(Widget w, int data, XtPointer call_data)
{
    curitem = data;
}

static void place_done_proc(void)
{
    XtUnmanageChild(place_frame);
}

void create_place_popup(void)
{
    int i, ac;
    char s[80];
    Widget bt, rc, lab, fr, rb, w[12];
    Arg a[3];
    XmString str;
    if (place_frame) {
	XtRaise(place_frame);
	return;
    }
    place_frame = XmCreateDialogShell(app_shell, "Place", NULL, 0);
    place_panel = XtVaCreateManagedWidget("Place panel", xmRowColumnWidgetClass, place_frame, NULL);
    ac = 0;
    XtSetArg(a[ac], XmNx, 0);
    ac++;
    XtSetArg(a[ac], XmNy, 0);
    ac++;
    lab = XmCreateLabelGadget(place_panel, "Place:", a, ac);
    ac = 0;
    XtSetArg(a[ac], XmNx, 0);
    ac++;
    XtSetArg(a[ac], XmNy, 35);
    ac++;
    fr = XmCreateFrame(place_panel, "frame_1", NULL, 0);
    rb = XmCreateRadioBox(fr, "radio_box_1", NULL, 0);
    w[0] = XmCreateToggleButton(rb, "Tidal clock", NULL, 0);
    w[1] = XmCreateToggleButton(rb, "Run info", NULL, 0);
    w[2] = XmCreateToggleButton(rb, "Time line", NULL, 0);
    w[3] = XmCreateToggleButton(rb, "Velocity scale", NULL, 0);
    w[4] = XmCreateToggleButton(rb, "Map scale", NULL, 0);
    w[5] = XmCreateToggleButton(rb, "Title", NULL, 0);
    w[6] = XmCreateToggleButton(rb, "Legend for bathymetry", NULL, 0);
    w[7] = XmCreateToggleButton(rb, "Legend for isolines", NULL, 0);
    for (i = 0; i < 8; i++) {
	XtAddCallback(w[i], XmNvalueChangedCallback, (XtCallbackProc) set_type_proc, (XtPointer) i);
    }
    XtManageChild(lab);
    XtManageChild(fr);
    XtManageChild(rb);
    XtManageChildren(w, 8);
    XmToggleButtonSetState(w[0], True, False);

    rc = XmCreateRowColumn(place_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) place_done_proc, NULL);

    bt = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) define_place_proc, NULL);
    XtManageChild(rc);

    XtManageChild(place_panel);
    XtRaise(place_frame);
}
