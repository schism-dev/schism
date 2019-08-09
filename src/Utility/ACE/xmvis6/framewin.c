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
 * frame Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: framewin.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include <stdio.h>

#include "motifinc.h"

#include "defines.h"
#include "globals.h"

static Widget frame_frame = (Widget) 0;
static Widget frame_panel;

/*
 * Widget item declarations
 */
static Widget *frame_frameactive_choice_item;
static Widget *frame_framestyle_choice_item;
static Widget *frame_color_choice_item;
static Widget *frame_lines_choice_item;
static Widget *frame_linew_choice_item;
static Widget *frame_fillbg_choice_item;
static Widget *frame_bgcolor_choice_item;
static Widget *frame_applyto_choice_item;

/*
 * Event and Notify proc declarations
 */
static int frame_Done_notify_proc(void);
static int frame_define_notify_proc(void);

void update_frame_items(int gno)
{
    if (frame_frame) {
	SetChoice(frame_frameactive_choice_item, g[gno].f.active == OFF);
	SetChoice(frame_framestyle_choice_item, g[gno].f.type);
	SetChoice(frame_color_choice_item, g[gno].f.color);
	SetChoice(frame_linew_choice_item, g[gno].f.linew - 1);
	SetChoice(frame_lines_choice_item, g[gno].f.lines - 1);
	SetChoice(frame_fillbg_choice_item, g[gno].f.fillbg == ON);
	SetChoice(frame_bgcolor_choice_item, g[gno].f.bgcolor);
    }
}

/*
 * Create the frame Widget and the frame Widget
 */
void create_frame_frame(void)
{
    extern Widget app_shell;
    Widget wbut, rc;

    if (frame_frame) {
	update_frame_items(cg);
	XtRaise(frame_frame);
	return;
    }
    frame_frame = XmCreateDialogShell(app_shell, "Frame", NULL, 0);
    frame_panel = XtVaCreateWidget("frame panel", xmRowColumnWidgetClass, frame_frame, NULL);

    frame_frameactive_choice_item = CreatePanelChoice1(frame_panel, "Graph frame:", 3, "ON", "OFF", NULL, NULL);
    frame_framestyle_choice_item = CreatePanelChoice1(frame_panel, "Frame type:", 3, "Closed", "Half open", NULL, NULL);
    frame_color_choice_item = CreateColorChoice(frame_panel, "Line color:", 1);

    frame_linew_choice_item = CreatePanelChoice1(frame_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    frame_lines_choice_item = CreatePanelChoice1(frame_panel, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    frame_fillbg_choice_item = CreatePanelChoice1(frame_panel, "Background fill:", 3, "None", "Filled", NULL, NULL);

    frame_bgcolor_choice_item = CreateColorChoice(frame_panel, "Background color:", 1);

    frame_applyto_choice_item = CreatePanelChoice1(frame_panel, "Apply to:", 3, "Current graph", "All active graphs", NULL, NULL);
    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, frame_panel, NULL);

    rc = XmCreateRowColumn(frame_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) frame_define_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) frame_Done_notify_proc, 0);
    XtManageChild(rc);

    update_frame_items(cg);
    XtManageChild(frame_panel);
    XtRaise(frame_frame);
}

/*
 * Notify and event procs
 */

static frame_Done_notify_proc(void)
{
    XtUnmanageChild(frame_frame);
}

static frame_define_notify_proc(void)
{
    int i, ming, maxg;
    int a = (int) GetChoice(frame_applyto_choice_item);
    if (a == 0) {
	ming = maxg = cg;
    } else {
	ming = 0;
	maxg = maxgraph - 1;
    }
    for (i = ming; i <= maxg; i++) {
	if (isactive_graph(i)) {
	    g[i].f.active = (int) GetChoice(frame_frameactive_choice_item) ? OFF : ON;
	    g[i].f.type = (int) GetChoice(frame_framestyle_choice_item);
	    g[i].f.color = (int) GetChoice(frame_color_choice_item);
	    g[i].f.linew = (int) GetChoice(frame_linew_choice_item) + 1;
	    g[i].f.lines = (int) GetChoice(frame_lines_choice_item) + 1;
	    g[i].f.fillbg = (int) GetChoice(frame_fillbg_choice_item) ? ON : OFF;
	    g[i].f.bgcolor = (int) GetChoice(frame_bgcolor_choice_item);
	}
    }
}
