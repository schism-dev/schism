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
 * display a timer
 *
 */

#ifndef lint
static char RCSid[] = "$Id: timerwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
/* used to set up XmStrings */
extern XmStringCharSet charset;

void create_timer_popup(void);

static Widget timer_frame;
static Widget timer_panel;

Widget timer_item;

void timer_done_proc(void);

static Widget write_mode_item;
static XmString mstring;

void set_timer_item(char *buf, double min, double max, double cur)
{
    Arg wargs[10];
    char tmpb[256];
    double d = (max - min);
    int per;
    per = (int) (cur / d * 100.0);
    if (per < 0) {
	per = 0;
    } else if (per > 100) {
	per = 100;
    }
    create_timer_popup();

    XtSetArg(wargs[0], XmNvalue, per);
    XtSetValues(timer_item, wargs, 1);

    if (buf == NULL) {
	strcpy(tmpb, "Idle...");
    } else {
	strcpy(tmpb, buf);
    }
    XmStringFree(mstring);
    mstring = XmStringCreateLtoR(tmpb, charset);
    XtSetArg(wargs[0], XmNlabelString, mstring);
    XtSetValues(write_mode_item, wargs, 1);

    XmUpdateDisplay(write_mode_item);
    XmUpdateDisplay(timer_item);
    XmUpdateDisplay(timer_panel);
    XmUpdateDisplay(timer_frame);
}

void create_timer_popup(void)
{
    Widget bt, fr, rc;
    Arg a[2];
    int ac;
    if (timer_frame) {
	XtManageChild(timer_frame);
	XmUpdateDisplay(timer_frame);
	return;
    }
    timer_frame = XmCreateDialogShell(app_shell, "% Complete", NULL, 0);
    timer_panel = XmCreateRowColumn(timer_frame, "rc", NULL, 0);

    timer_item = XtVaCreateManagedWidget("timer_item", xmScaleWidgetClass, timer_panel,
					 XmNwidth, 250,
					 XmNminimum, 0,
					 XmNmaximum, 100,
					 XmNvalue, 90,
					 XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
		  XmNtitleString, XmStringCreateLtoR("% complete", charset),
					 XmNorientation, XmHORIZONTAL,
					 NULL);
    mstring = XmStringCreateLtoR("msg\ninit\nET", charset);
    write_mode_item = XmCreateLabel(timer_panel, "msg\ninit\nET", NULL, 0);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, timer_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) timer_done_proc, NULL);
    XtManageChild(write_mode_item);
    XtManageChild(timer_panel);
    XtManageChild(timer_frame);
    XmUpdateDisplay(timer_frame);
}

void timer_done_proc(void)
{
    XmUpdateDisplay(timer_frame);
    XtUnmanageChild(timer_frame);
    XmUpdateDisplay(timer_frame);
}
