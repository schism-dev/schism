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
 * Set velocity scale, mapscale, etc.
 */

#ifndef lint
static char RCSid[] = "$Id: scalewin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";

#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

Widget annot_frame;
Widget annot_panel;

Widget vel_frame;
Widget vel_panel;

#define M 362
#define CM 317
#define MM 368
#define KM 350

int munits = KM;
double munitfac = 1000.0;

double mapx;
double mapy;
double maplen = 10000.0;
int my_mapcolor = 1;
extern int display_mapscale;

static Widget mapscale_item;
static Widget mapscale_length_item;
static Widget *mapscale_color_item;
static Widget *mapscale_font_item;
static Widget *mapscale_units_item;

static void update_scaling(void);
static void vel_accept_proc(void);
static void define_vplace_proc(void);

void drawscale(void)
{
    char buf[30];
    if (maplen > 0.0) {
	sprintf(buf, "%lg", maplen / munitfac);
	drawmapscale(mapx, mapy, maplen, buf);
    }
}

void set_mapscale_loc(double wx, double wy)
{
    mapx = wx;
    mapy = wy;
}

/*
 * Create the annot Frame and the annot Panel
 */


void create_vel_frame(void)
{
    Widget wbut, rc;

    if (!vel_frame) {
	vel_frame = XmCreateDialogShell(app_shell, "Map Scale", NULL, 0);
	vel_panel = XmCreateRowColumn(vel_frame, "rc", NULL, 0);
	mapscale_length_item = CreateTextItem2(vel_panel, 15, "Mapscale legend length:");
	mapscale_units_item = CreatePanelChoice1(vel_panel, "Mapscale legend units: ",
						 5,
						 "mm", "cm", "m", "km",
						 0, 0);
	mapscale_font_item = CreatePanelChoice1(vel_panel, "Legend label font: ",
						7,
						"0", "1", "2", "3", "4", "5",
						0, 0);
	mapscale_color_item = CreateColorChoice(vel_panel, "Mapscale legend color: ", 1);
	mapscale_item = XtVaCreateManagedWidget("Display mapscale legend",
				       xmToggleButtonWidgetClass, vel_panel,
						NULL);

	rc = XmCreateRowColumn(vel_panel, "rc", NULL, 0);
	XtVaSetValues(rc,
		      XmNorientation, XmHORIZONTAL,
		      XmNpacking, XmPACK_TIGHT,
		      NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) vel_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_vplace_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, vel_frame);
	XtManageChild(rc);
	XtManageChild(vel_panel);
    }
    update_scaling();
    XtRaise(vel_frame);
}

static void define_vplace_proc(void)
{
    set_action(0);
    set_action(PLACE_MAPSCALE);
}

static void update_scaling(void)
{
    char buf[256];
    extern int units, munits;
    if (vel_frame) {
	sprintf(buf, "%.2lf\n", maplen / munitfac);
	panel_setstr_value(mapscale_length_item, buf);
	XmToggleButtonSetState(mapscale_item, display_mapscale, False);
	SetChoice(mapscale_color_item, my_mapcolor);
	switch (munits) {
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

static void vel_accept_proc(void)
{
    char buf[256];
    extern int units;
    extern double unitfac;
    int u;
    display_mapscale = XmToggleButtonGetState(mapscale_item) == True ? 1 : 0;
    my_mapcolor = GetChoice(mapscale_color_item);
    u = GetChoice(mapscale_units_item);
    switch (u) {
    case 0:
	munits = MM;
	munitfac = 0.001;
	break;
    case 1:
	munits = CM;
	munitfac = 0.01;
	break;
    case 2:
	munits = M;
	munitfac = 1.0;
	break;
    case 3:
	munits = KM;
	munitfac = 1000.0;
	break;
    }
    maplen = atof(XmTextGetString(mapscale_length_item)) * munitfac;
}
