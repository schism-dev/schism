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
 * annotwin.c - Set velocity scale, mapscale, and parameters
 * associated with the time line tidal clock, etc.
 *
 */

#ifndef lint
static char RCSid[] = "$Id: annotwin.c,v 1.5 2007/03/13 15:00:57 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget timeline_frame;
static Widget tidalclock_frame;
static Widget vel_frame;
static Widget map_frame;
static Widget flux_frame;

/*
 * Panel item declarations
 */
static Widget fr2, fr3, fr4, fr5, vfr1, vfr2;

static Widget mapscale_item;
static Widget mapscale_length_item;
static Widget *mapscale_color_item;
static Widget *mapscale_font_item;
static Widget *mapscale_units_item;
static Widget *mapscale_loctype_item;

static Widget vscale_item;
static Widget *scale_item;
static Widget *vscale_units_item;
static Widget vscale_val_item;
static Widget vscale_length_item;
static Widget vscale_llength_item;
static Widget *vscale_color_item;
static Widget *vscale_font_item;
static Widget *vscale_loctype_item;

static Widget fscale_item;
static Widget *fscale_units_item;
static Widget fscale_val_item;
static Widget fscale_length_item;
static Widget fscale_llength_item;
static Widget *fscale_color_item;
static Widget *fscale_font_item;
static Widget *fscale_loctype_item;

static Widget title_item;
static Widget title_text_item;
static Widget *title_font_item;
static Widget *title_color_item;

static Widget clock_item;
static Widget clock_time_item;
static Widget *clock_color_item;
static Widget *clock_fillcolor_item;
static Widget *clock_loctype_item;

static Widget timeline_item;
static Widget timeline_start_item;
static Widget timeline_stop_item;
static Widget timeline_step_item;
static Widget *timeline_font_item;
static Widget *timeline_color_item;
static Widget *timeline_fillcolor_item;
static Widget *timeline_loctype_item;
static Widget *timeline_prec_item;
static Widget *timeline_units_item;

/*
 * the following update... procs are called here and in update_all();
 */
void update_timeline(void);
void update_tidalclock(void);
void update_scaling(void);
void update_flux_scaling(void);
void update_map_scaling(void);

static void timeline_accept_proc(void);
static void tidalclock_accept_proc(void);

static void vel_accept_proc(void);
static void map_accept_proc(void);
static void flux_accept_proc(void);

static void define_vplace_proc(void);
static void define_fplace_proc(void);
static void define_mplace_proc(void);

/*
 * return the string equiv of the scaling parameter
 */
char *units_str(int u)
{
    static char buf[256];
    switch (u) {
    case MM:
	strcpy(buf, "mm");
	break;
    case CM:
	strcpy(buf, "cm");
	break;
    case M:
	strcpy(buf, "m");
	break;
    case KM:
	strcpy(buf, "km");
	break;
    }
    return buf;
}

/*
 * place the time line
 */
static void place_timeline_proc(void)
{
    set_action(0);
    set_action(PLACE_TIMELINE);
}

/*
 * place the tidal clock
 */
static void place_tidalclock_proc(void)
{
    set_action(0);
    set_action(PLACE_CLOCK);
}

/*
 * popup for time line properties
 */
void create_timeline_frame(void)
{
    Arg wargs[8];
    int i;
    Widget timeline_panel, wbut, fr, rc, rb, bb, sep, w[10];

    setistop();
    if (!timeline_frame) {
	timeline_frame = XmCreateDialogShell(app_shell, "Timeline setup", NULL, 0);
	handle_close(timeline_frame);
	timeline_panel = XmCreateRowColumn(timeline_frame, "rc", NULL, 0);

	rc = XmCreateRowColumn(timeline_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	fr3 = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, rc, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr3, NULL);
	timeline_start_item = CreateTextItem2(bb, 15, "Start time: ");
	timeline_stop_item = CreateTextItem2(bb, 15, "End time: ");
	timeline_step_item = CreateTextItem2(bb, 15, "Label every: ");

	timeline_color_item = CreateColorChoice(bb, "Outline color: ", 1);
	timeline_fillcolor_item = CreateColorChoice(bb, "Fill color: ", 1);
	timeline_font_item = CreatePanelChoice1(bb, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);

	timeline_prec_item = CreatePanelChoice1(bb, "Labels precision: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	timeline_units_item = CreatePanelChoice1(bb, "Labels units: ", 8, "seconds", "minutes", "hours", "days", "weeks", "months", "years", 0, 0);
	timeline_loctype_item = CreatePanelChoice1(bb, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);

	timeline_item = XtVaCreateManagedWidget("Display time line", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(bb);
	XtManageChild(rc);

	rc = XmCreateRowColumn(timeline_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) timeline_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) place_timeline_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, timeline_frame);
	XtManageChild(rc);

	XtManageChild(timeline_panel);
    }
    XtRaise(timeline_frame);
    update_timeline();
}

/*
 * update the time line popup
 */
void update_timeline(void)
{
    char buf[256];
    if (timeline_frame) {
	XmToggleButtonSetState(timeline_item, g[cg].timeline.active == ON, False);
	sprintf(buf, "%.2lf", g[cg].timeline.start);
	xv_setstr(timeline_start_item, buf);
	sprintf(buf, "%.2lf", g[cg].timeline.stop);
	xv_setstr(timeline_stop_item, buf);
	sprintf(buf, "%.2lf", g[cg].timeline.step);
	xv_setstr(timeline_step_item, buf);
	SetChoice(timeline_color_item, g[cg].timeline.c1);
	SetChoice(timeline_fillcolor_item, g[cg].timeline.c2);
	SetChoice(timeline_loctype_item, g[cg].timeline.loctype == VIEW);
	SetChoice(timeline_units_item, g[cg].timeline.units);
	SetChoice(timeline_prec_item, g[cg].timeline.p.prec);
/*
	SetChoice(timeline_font_item);
*/
    }
}

/*
 * register the current state of the time line popup
 */
static void timeline_accept_proc(void)
{
    g[cg].timeline.active = XmToggleButtonGetState(timeline_item) ? ON : OFF;
    g[cg].timeline.start = atof((char *) xv_getstr(timeline_start_item));
    g[cg].timeline.stop = atof((char *) xv_getstr(timeline_stop_item));
    g[cg].timeline.step = atof((char *) xv_getstr(timeline_step_item));
    g[cg].timeline.loctype = GetChoice(timeline_loctype_item) ? VIEW : WORLD;
    g[cg].timeline.units = GetChoice(timeline_units_item);
    g[cg].timeline.c1 = g[cg].timeline.c3 = GetChoice(timeline_color_item);
    g[cg].timeline.c2 = GetChoice(timeline_fillcolor_item);
    g[cg].timeline.p.prec = GetChoice(timeline_prec_item);
    GetChoice(timeline_font_item);
}

/*
 * create a property sheet for the tidal clock
 */
void create_tidalclock_frame(void)
{
    int i;
    Widget tidalclock_panel, wbut, fr, rc, rb, bb, sep;

    setistop();
    if (!tidalclock_frame) {
	tidalclock_frame = XmCreateDialogShell(app_shell, "Tidal clock setup", NULL, 0);
	handle_close(tidalclock_frame);
	tidalclock_panel = XmCreateRowColumn(tidalclock_frame, "rc", NULL, 0);

	rc = XmCreateRowColumn(tidalclock_panel, "rc", NULL, 0);

	fr2 = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, rc, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr2, NULL);
	clock_time_item = CreateTextItem2(bb, 10, "Total time: ");
	clock_color_item = CreateColorChoice(bb, "Clock outline color: ", 1);
	clock_fillcolor_item = CreateColorChoice(bb, "Clock fill color: ", 1);
	clock_loctype_item = CreatePanelChoice1(bb, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);
	clock_item = XtVaCreateManagedWidget("Display tidal clock", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(rc);

	rc = XmCreateRowColumn(tidalclock_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) tidalclock_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) place_tidalclock_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, tidalclock_frame);
	XtManageChild(rc);

	XtManageChild(tidalclock_panel);
    }
    XtRaise(tidalclock_frame);
    update_tidalclock();
}

/*
 * make the state of the tidal clock current
 */
void update_tidalclock(void)
{
    char buf[256];
    if (tidalclock_frame) {
	XmToggleButtonSetState(clock_item, g[cg].tidalclock.active == ON, False);
	SetChoice(clock_color_item, g[cg].tidalclock.p.color);
	SetChoice(clock_fillcolor_item, g[cg].tidalclock.p.fillcol);
	SetChoice(clock_loctype_item, g[cg].tidalclock.loctype == VIEW);
	sprintf(buf, "%.2lf", g[cg].tidalclock.total_time);
	xv_setstr(clock_time_item, buf);
    }
}

static void tidalclock_accept_proc(void)
{
    g[cg].tidalclock.active = XmToggleButtonGetState(clock_item) ? ON : OFF;
    g[cg].tidalclock.p.color = GetChoice(clock_color_item);
    g[cg].tidalclock.p.fillcol = GetChoice(clock_fillcolor_item);
    g[cg].tidalclock.loctype = GetChoice(clock_loctype_item) ? VIEW : WORLD;
    g[cg].tidalclock.total_time = atof((char *) xv_getstr(clock_time_item));
}

/*
 * popup to set the map scale legend
 */
void create_map_frame(void)
{
    int i;
    Widget map_panel, wbut, fr, rc, rb, bb, sep, w[10];

    setistop();
    if (!map_frame) {
	map_frame = XmCreateDialogShell(app_shell, "Scaling", NULL, 0);
	handle_close(map_frame);
	map_panel = XmCreateRowColumn(map_frame, "rc", NULL, 0);

	vfr1 = XmCreateFrame(map_panel, "frame", NULL, 0);
	bb = XmCreateRowColumn(vfr1, "bb", NULL, 0);
	mapscale_length_item = CreateTextItem2(bb, 15, "Mapscale legend length:");
	mapscale_units_item = CreatePanelChoice1(bb, "Mapscale legend units: ", 5, "mm", "cm", "m", "km", 0, 0);
	mapscale_font_item = CreatePanelChoice1(bb, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
	mapscale_color_item = CreateColorChoice(bb, "Mapscale legend color: ", 1);
	mapscale_item = XtVaCreateManagedWidget("Display mapscale legend", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(bb);
	XtManageChild(vfr1);

	rc = XmCreateRowColumn(map_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) map_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_mplace_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, map_frame);
	XtManageChild(rc);
	XtManageChild(map_panel);
    }
    XtRaise(map_frame);
    update_map_scaling();
}

/*
 * place the map scale legend for the current graph
 */
static void define_mplace_proc(void)
{
    set_action(0);
    set_action(PLACE_MAPSCALE);
}

/*
 * update the map scale popup
 */
void update_map_scaling(void)
{
    char buf[256];
    if (vel_frame) {
	sprintf(buf, "%.2lf\n", g[cg].mapscale.len);
	xv_setstr(mapscale_length_item, buf);
	XmToggleButtonSetState(mapscale_item, g[cg].mapscale.active == ON, False);
	SetChoice(mapscale_color_item, g[cg].mapscale.p.color);
	switch (g[cg].mapscale.units) {
	case MM:
	    SetChoice(mapscale_units_item, 0);
	    break;
	case CM:
	    SetChoice(mapscale_units_item, 1);
	    break;
	case M:
	    SetChoice(mapscale_units_item, 2);
	    break;
	case KM:
	    SetChoice(mapscale_units_item, 3);
	    break;
	}
    }
}

/*
 * accept the current settings
 */
static void map_accept_proc(void)
{
    char buf[256];
    int u;
    g[cg].mapscale.active = XmToggleButtonGetState(mapscale_item) == True ? ON : OFF;
    g[cg].mapscale.p.color = GetChoice(mapscale_color_item);
    u = GetChoice(mapscale_units_item);
    switch (u) {
    case 0:
	g[cg].mapscale.units = MM;
	g[cg].mapscale.unitfac = 0.001;
	break;
    case 1:
	g[cg].mapscale.units = CM;
	g[cg].mapscale.unitfac = 0.01;
	break;
    case 2:
	g[cg].mapscale.units = M;
	g[cg].mapscale.unitfac = 1.0;
	break;
    case 3:
	g[cg].mapscale.units = KM;
	g[cg].mapscale.unitfac = 1000.0;
	break;
    }
    g[cg].mapscale.len = atof(XmTextGetString(mapscale_length_item));
}

void update_scaling(void);

/*
 * define the velocity scale for this graph and properties of the
 * velocity scale legend
 */
void create_vel_frame(void)
{
    Arg wargs[8];
    int i;
    Widget vel_panel, wbut, fr, rc, rb, bb, sep, w[10];

    setistop();
    if (!vel_frame) {
	vel_frame = XmCreateDialogShell(app_shell, "Scaling", NULL, 0);
	handle_close(vel_frame);
	vel_panel = XmCreateRowColumn(vel_frame, "rc", NULL, 0);

	vfr2 = XmCreateFrame(vel_panel, "frame", NULL, 0);
	bb = XmCreateRowColumn(vfr2, "rc", NULL, 0);
	vscale_length_item = CreateTextItem2(bb, 15, "1 cm on screen = ");
	scale_item = CreatePanelChoice1(bb, "Set scale for: ", 3, "Velocities", "Wind", 0, 0);
	XtAddCallback(scale_item[2], XmNactivateCallback, (XtCallbackProc) update_scaling, (XtPointer) NULL);
	XtAddCallback(scale_item[3], XmNactivateCallback, (XtCallbackProc) update_scaling, (XtPointer) NULL);

	vscale_units_item = CreatePanelChoice1(bb, "Units of velocity: ", 5, "mm/sec", "cm/sec", "m/sec", "km/sec", 0, 0);

	vscale_llength_item = CreateTextItem2(bb, 15, "Velocity scale legend length: ");
	vscale_color_item = CreateColorChoice(bb, "Velocity scale legend color: ", 1);
	vscale_font_item = CreatePanelChoice1(bb, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
	vscale_item = XtVaCreateManagedWidget("Display velocity scale legend", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(bb);
	XtManageChild(vfr2);

	rc = XmCreateRowColumn(vel_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) vel_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_vplace_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, vel_frame);
	XtManageChild(rc);
	XtManageChild(vel_panel);
    }
    XtRaise(vel_frame);
    update_scaling();
}

/*
 * place the velocity scale
 */
static void define_vplace_proc(void)
{
    set_action(0);
    if (GetChoice(scale_item)) {
	set_action(PLACE_WSCALE);
    } else {
	set_action(PLACE_VSCALE);
    }
}

/*
 * update the velocity scale popup
 */
void update_scaling(void)
{
    char buf[256];
    velocity_scale vl;
    if (vel_frame) {
	int c = GetChoice(scale_item);
	if (c == 0) {
	    vl = g[cg].vl;
	} else {
	    vl = g[cg].wl;
	}
	sprintf(buf, "%lf\n", 1.0 / vl.scale * 10.0);
	xv_setstr(vscale_length_item, buf);
	sprintf(buf, "%.2lf\n", vl.len);
	xv_setstr(vscale_llength_item, buf);
	XmToggleButtonSetState(vscale_item, vl.active == ON, False);
	SetChoice(vscale_color_item, vl.p.color);
	switch (vl.units) {
	case MM:
	    SetChoice(vscale_units_item, 0);
	    break;
	case CM:
	    SetChoice(vscale_units_item, 1);
	    break;
	case M:
	    SetChoice(vscale_units_item, 2);
	    break;
	case KM:
	    SetChoice(vscale_units_item, 3);
	    break;
	}
    }
}

/*
 * accept the current state of the velocity popup
 */
static void vel_accept_proc(void)
{
    char buf[256];
    velocity_scale *vl;
    int u = GetChoice(vscale_units_item);
    int c = GetChoice(scale_item);
    if (c == 0) {
	vl = &g[cg].vl;
    } else {
	vl = &g[cg].wl;
    }
    vl->scale = 10.0 / atof(XmTextGetString(vscale_length_item));
    vl->len = atof(XmTextGetString(vscale_llength_item));
    vl->active = XmToggleButtonGetState(vscale_item) == True ? ON : OFF;
    vl->p.color = GetChoice(vscale_color_item);
    switch (u) {
    case 0:
	vl->units = MM;
	vl->unitfac = 0.001;
	break;
    case 1:
	vl->units = CM;
	vl->unitfac = 0.01;
	break;
    case 2:
	vl->units = M;
	vl->unitfac = 1.0;
	break;
    case 3:
	vl->units = KM;
	vl->unitfac = 1000.0;
	break;
    }
}

/*
 * create the flux scale popup
 */
void create_fluxscale_frame(void)
{
    Arg wargs[8];
    int i;
    Widget flux_panel, wbut, fr, rc, rb, bb, sep, w[10];

    setistop();
    if (!flux_frame) {
	flux_frame = XmCreateDialogShell(app_shell, "Scaling", NULL, 0);
	handle_close(flux_frame);
	flux_panel = XmCreateRowColumn(flux_frame, "rc", NULL, 0);

	vfr2 = XmCreateFrame(flux_panel, "frame", NULL, 0);
	bb = XmCreateRowColumn(vfr2, "rc", NULL, 0);
	fscale_length_item = CreateTextItem2(bb, 15, "1 cm on screen = ");
	fscale_units_item = CreatePanelChoice1(bb, "Units of flux: ", 5, "mm^3/sec", "cm^3/sec", "m^3/sec", "km^3/sec", 0, 0);

	fscale_llength_item = CreateTextItem2(bb, 15, "Flux scale legend length: ");
	fscale_color_item = CreateColorChoice(bb, "Flux scale legend color: ", 1);
	fscale_font_item = CreatePanelChoice1(bb, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
	fscale_item = XtVaCreateManagedWidget("Display flux scale legend", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(bb);
	XtManageChild(vfr2);

	rc = XmCreateRowColumn(flux_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) flux_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_fplace_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, flux_frame);
	XtManageChild(rc);
	XtManageChild(flux_panel);
    }
    XtRaise(flux_frame);
    update_flux_scaling();
}

/*
 * set the posistion of the flux scale for this graph
 */
static void define_fplace_proc(void)
{
    set_action(0);
    set_action(PLACE_FSCALE);
}

/*
 * update the flux scale popup
 */
void update_flux_scaling(void)
{
    char buf[256];
    if (flux_frame) {
	sprintf(buf, "%lf\n", 1.0 / g[cg].fl.scale * 10.0);
	xv_setstr(fscale_length_item, buf);
	sprintf(buf, "%.2lf\n", g[cg].fl.len);
	xv_setstr(fscale_llength_item, buf);
	XmToggleButtonSetState(fscale_item, g[cg].fl.active == ON, False);
	SetChoice(fscale_color_item, g[cg].fl.p.color);
	switch (g[cg].fl.units) {
	case MM:
	    SetChoice(fscale_units_item, 0);
	    break;
	case CM:
	    SetChoice(fscale_units_item, 1);
	    break;
	case M:
	    SetChoice(fscale_units_item, 2);
	    break;
	case KM:
	    SetChoice(fscale_units_item, 3);
	    break;
	}
    }
}

/*
 * accept the current state of the flux scale popup
 */
static void flux_accept_proc(void)
{
    char buf[256];
    int u = GetChoice(fscale_units_item);
    g[cg].fl.scale = 10.0 / atof(XmTextGetString(fscale_length_item));
    g[cg].fl.len = atof(XmTextGetString(fscale_llength_item));
    g[cg].fl.active = XmToggleButtonGetState(fscale_item) == True ? ON : OFF;
    g[cg].fl.p.color = GetChoice(fscale_color_item);
    switch (u) {
    case 0:
	g[cg].fl.units = MM;
	g[cg].fl.unitfac = 0.001;
	break;
    case 1:
	g[cg].fl.units = CM;
	g[cg].fl.unitfac = 0.01;
	break;
    case 2:
	g[cg].fl.units = M;
	g[cg].fl.unitfac = 1.0;
	break;
    case 3:
	g[cg].fl.units = KM;
	g[cg].fl.unitfac = 1000.0;
	break;
    }
}

typedef struct {
    Widget top;
    Widget active;
    Widget *display;
    Widget *font;
    Widget *color;
    Widget *loc;
    Widget size;
    Widget start;
    Widget format;
    Widget frame;
    Widget *framefill;
    Widget *framefillcolor;
} timeUI;

void update_timeinfo(timeUI * ui);
void time_info_start(int gno);
static XtCallbackProc timeinfo_accept_proc(Widget w, XtPointer clid, XtPointer calld);

/*
 * place the timeinfo string
 */
static void place_timeinfo_proc(void)
{
    set_action(0);
    set_action(TIMEINFO);
}

/*
 * Display the timeinfo string and set properties
 */
void create_timeinfo_frame(void)
{
    static timeUI ui;
    Widget panel, wbut, rc;
    Widget wlabel;
    int x, y;

    setistop();
    if (ui.top) {
	update_timeinfo(&ui);
	XtRaise(ui.top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    ui.top = XmCreateDialogShell(app_shell, "Time Information", NULL, 0);
    handle_close(ui.top);
    XtVaSetValues(ui.top, XmNx, x, XmNy, y, NULL);
    panel = XmCreateRowColumn(ui.top, "strings_rc", NULL, 0);

    ui.active = XtVaCreateManagedWidget("Toggle time info string", xmToggleButtonWidgetClass, panel, NULL);
    ui.display = CreatePanelChoice1(panel, "Display:", 3, "Formatted", "DD:HH:MM:SS", 0, 0);
    ui.start = CreateTextItem2(panel, 12, "Start time (YYYYMMDDHHMMSS): ");
    ui.format = CreateTextItem2(panel, 15, "Time format: ");
    ui.font = CreatePanelChoice1(panel, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
    ui.color = CreateColorChoice(panel, "Color:", 1);
    ui.loc = CreatePanelChoice1(panel, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Size:", xmLabelGadgetClass, rc, NULL);
    ui.size = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) timeinfo_accept_proc, (XtPointer) & ui);
    wbut = XtVaCreateManagedWidget("Pick Location", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) place_timeinfo_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
    XtManageChild(rc);

    update_timeinfo(&ui);
    XtManageChild(panel);
    XtManageChild(ui.top);
}

/*
 * make the state of the tidal clock current
 */
void update_timeinfo(timeUI * ui)
{
    Arg a;
    int iv;
    if (ui->top) {
	xv_setstr(ui->format, g[cg].timeinfo.format);
	xv_setstr(ui->start, g[cg].timeinfo.start);
	XmToggleButtonSetState(ui->active, g[cg].timeinfo.active == ON, False);
	SetChoice(ui->display, g[cg].timeinfo.display);
	SetChoice(ui->font, g[cg].timeinfo.font);
	SetChoice(ui->color, g[cg].timeinfo.color);
	iv = (int) (100 * g[cg].timeinfo.charsize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(ui->size, &a, 1);
	SetChoice(ui->loc, g[cg].timeinfo.loctype == VIEW);
    }
}

static XtCallbackProc timeinfo_accept_proc(Widget w, XtPointer clid, XtPointer calld)
{
    timeUI *ui = (timeUI *) clid;
    Arg a;
    int value;
    if (ui->top) {
	int m, d, y, h, mm, s;
	int display;
	char buf[256];
	struct tm t;
	g[cg].timeinfo.display = (int) GetChoice(ui->display);
	if (g[cg].timeinfo.display == 0) {
	    strncpy(buf, XmTextGetString(ui->start), 255);
	    if (strlen(buf) != 14) {
		errwin("Start date format may be incorrect, string does not have 14 characters");
	    }
	    strncpy(g[cg].timeinfo.start, buf, 14);
	    time_info_start(cg);
	    strncpy(buf, XmTextGetString(ui->format), 255);
	    strcpy(g[cg].timeinfo.format, buf);
	} else {		/* display other stuff like DD:HH:MM format */
	}
	g[cg].timeinfo.active = XmToggleButtonGetState(ui->active) ? ON : OFF;
	g[cg].timeinfo.font = (int) GetChoice(ui->font);
	g[cg].timeinfo.color = (int) GetChoice(ui->color);
	XtSetArg(a, XmNvalue, &value);
	XtGetValues(ui->size, &a, 1);
	g[cg].timeinfo.charsize = value / 100.0;
	g[cg].timeinfo.loctype = (int) GetChoice(ui->loc) ? VIEW : WORLD;
    }
}

void time_info_start(int gno)
{
    int m, d, y, h, mm, s;
    char buf[256];
    char ebuf[1024];
    struct tm t;
    strcpy(buf, g[gno].timeinfo.start);
    if (strlen(buf) != 14) {
	sprintf(ebuf, "Start date format may be incorrect, string does not have 14 characters: %s", buf);
	errwin(ebuf);
    }
    sscanf(buf, "%4d%2d%2d%2d%2d%2d", &y, &m, &d, &h, &mm, &s);
    t.tm_year = y - 1900;
    t.tm_mon = m - 1;
    t.tm_mday = d;
    t.tm_hour = h;
    t.tm_min = mm;
    t.tm_sec = s;
    t.tm_isdst = 0;
    g[gno].timeinfo.bdtime = t;
    g[gno].timeinfo.time = mktime(&t);
}

void settime_info_start(int gno, char *start)
{
    int yy, mo, dd, hh, mm, ss;
    char tz[10], buf[1024];
    sscanf(start, "%2d/%2d/%4d %2d:%2d:%2d %s", &mo, &dd, &yy, &hh, &mm, &ss, tz);
    sprintf(buf, "%4d%02d%02d%02d%02d%02d", yy, mo, dd, hh, mm, ss);
    strcpy(g[gno].timeinfo.start, buf);
    time_info_start(gno);
}
