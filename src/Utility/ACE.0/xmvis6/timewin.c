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
 * VIS clock routine
 */

#ifndef lint
static char RCSid[] = "$Id: timewin.c,v 1.4 2003/10/07 04:40:14 pturner Exp $";

#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget time_frame;
static Widget time_panel;

void create_time_frame(void);
static void update_time(void);
void set_clock(int type, double start, double stop, double step, int nsteps);

/*
 * Panel item declarations
 */
static Widget time_start_item;
static Widget time_step_item;
static Widget time_nsteps_item;
static Widget time_interp_item;
static Widget *time_skip_item;

/*
 * Event and Notify proc declarations
 */
void time_accept_proc(void);

void create_time_frame(void)
{
    Widget wbut, fr, bb, sep, rc, rc2;
    int i, got_one = 0;;
    setistop();
    if (!time_frame) {
	time_frame = XmCreateDialogShell(app_shell, "Time", NULL, 0);
	handle_close(time_frame);
	time_panel = XmCreateRowColumn(time_frame, "time_rc", NULL, 0);
	time_start_item = CreateTextItem2(time_panel, 10, "Start time:");
	time_step_item = CreateTextItem2(time_panel, 10, "Time step:");
	time_nsteps_item = CreateTextItem2(time_panel, 10, "Number of steps:");
	rc = XmCreateRowColumn(time_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) time_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, time_frame);
	XtManageChild(rc);
	XtManageChild(time_panel);
    }
    update_time();
    XtRaise(time_frame);
}

void print_clock(void)
{
    printf("Start = %lf\n", timeclock.start);
    printf("Stop = %lf\n", timeclock.stop);
    printf("Step = %lf\n", timeclock.step);
    printf("Nsteps = %d\n", timeclock.nsteps);
}

void time_accept_proc(void)
{
    char buf[256];

    timeclock.type = 1;
    timeclock.start = atof((char *) xv_getstr(time_start_item));
    timeclock.step = atof((char *) xv_getstr(time_step_item));
    timeclock.nsteps = atoi((char *) xv_getstr(time_nsteps_item));
    set_clock(timeclock.type, timeclock.start, timeclock.stop, timeclock.step, timeclock.nsteps);
}

static void update_time(void)
{
    char buf[256];

    if (time_frame) {
	sprintf(buf, "%.3lf", timeclock.start);
	xv_setstr(time_start_item, buf);
	sprintf(buf, "%.3lf", timeclock.step);
	xv_setstr(time_step_item, buf);
	sprintf(buf, "%d", timeclock.nsteps);
	xv_setstr(time_nsteps_item, buf);
    }
}

void set_clock(int type, double start, double stop, double step, int nsteps)
{
    char buf[256];
    int i;

    if (timeclock.nsteps > 0 && timeclock.t != NULL) {
	free(timeclock.t);
    }
    timeclock.start = start;
    timeclock.stop = stop;
    timeclock.step = step;
    timeclock.curstep = 0;
    timeclock.nsteps = nsteps;
    timeclock.t = (double *) malloc(nsteps * sizeof(double));
    if (type == 1) {		/* setting time from the UI */
	for (i = 0; i < nsteps; i++) {
	    timeclock.t[i] = start + i * step;
	}
    }
    timeclock.curtime = timeclock.t[0];
    update_time();
}
