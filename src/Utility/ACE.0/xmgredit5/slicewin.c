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
 * popup for slicing a grid
 *
 */

#ifndef lint
static char RCSid[] = "$Id: slicewin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

int npts_to_slice = 0;
int slice_mode = 0;
int slice_out = 0;
double slice_delta;
char slicefile[256];
char grid1dfile[256];
char slicebackfile[256];
int slice_back = 0;

extern Widget graph2d_frame;
static Widget slice_frame;
static Widget slice_panel;
static Widget *slice_method;
static Widget slice_npts_item;
static Widget slice_file;
static Widget slice_file2;
static Widget slice_back_item;
static Widget grid_to_read;
static Widget *slice_output;

char *panel_getstr_value();

void slice_bath_line(void)
{
    write_mode_str("Select start of line for slice");
    set_action(0);
    set_action(SLICE_BATH1ST);
}

void slice_bath_gridline(void)
{
    write_mode_str("Select start of line for slice");
    set_action(0);
    set_action(SLICE_GRIDBATH1ST);
}

double slicex[400], slicey[400];
int slice_npts = 0;

void slice_points(void)
{
    write_mode_str("Use left mouse button to pick points, middle button to register, right to cancel");
    slice_npts = 0;
    set_action(0);
    set_action(PICK_BATH);
}

void write_slice(int gno, int which)
{
    int i, ind, last = 0;
    double get_depth_element(), bathx, x, y, dist;
    FILE *fp;
    char buf[256];
    if (which) {
	fp = fopen(slicefile, "w");
    } else {
	fp = fopen(slicebackfile, "w");
	slice_back = 0;
	last = 1;
    }
    if (fp == NULL) {
	errwin("Error opening file for output");
	return;
    }
    fprintf(fp, "Slice from xmgredit\n");
    if (slice_out) {
	fprintf(fp, "%d\n", slice_npts);
	for (i = 0; i < slice_npts; i++) {
	    x = slicex[i];
	    y = slicey[i];
	    find_element(gno, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(gno, ind, x, y);
		fprintf(fp, "%d %lf %lf %lf\n", i + 1, x, y, bathx);
	    }
	}
    } else {
	fprintf(fp, "1\n");
	fprintf(fp, "%d %d\n", slice_npts - 1, slice_npts);
	fprintf(fp, "1\n");
	fprintf(fp, "2 2\n");
	fprintf(fp, "1 1 %d 1 %d\n", slice_npts, slice_npts - 1);
	dist = 0.0;
	for (i = 0; i < slice_npts; i++) {
	    find_element(gno, slicex[i], slicey[i], &ind);
	    bathx = get_depth_element(gno, ind, slicex[i], slicey[i]);
	    if (i > 0) {
		dist = dist + hypot(slicex[i] - slicex[i - 1], slicey[i] - slicey[i - 1]);
	    }
	    fprintf(fp, "%lf 0.0 %lf 1 2\n", dist, bathx);
	}
	for (i = 0; i < slice_npts - 1; i++) {
	    fprintf(fp, "10 2 2 %d %d 1 2\n", i + 1, i + 2);
	}
	fprintf(fp, "1.000000000  10000.00000\n");
	fprintf(fp, "5000.000000  50000000.00\n");
	fprintf(fp, "1.000000000  10002.00000\n");
	fprintf(fp, "5000.000000  20000.00000\n");
    }
    fclose(fp);
    if (last) {
	sprintf(buf, "xmgr %s %s", slicefile, slicebackfile);
    }
}

static void do_slice_proc(void)
{
    char buf[256];
    slice_mode = GetChoice(slice_method);
    strcpy(buf, panel_getstr_value(slice_npts_item));
    sscanf(buf, "%d", &npts_to_slice);
    strcpy(slicefile, panel_getstr_value(slice_file));
    if (XmToggleButtonGetState(slice_back_item)) {
	strcpy(slicebackfile, panel_getstr_value(slice_file2));
	if (slicebackfile[0]) {
	    if (!fexists(slicebackfile)) {
	    }
	} else {
	    errwin("Define file name for background slice first");
	    slice_back = 0;
	    return;
	}
	slice_back = 1;
    }
    strcpy(grid1dfile, panel_getstr_value(grid_to_read));
    slice_out = GetChoice(slice_output);
    if (slicefile[0]) {
	if (!fexists(slicefile)) {
	}
    } else {
	errwin("Define file name first");
	return;
    }

    switch (slice_mode) {
    case 0:
	slice_bath_line();
	break;
    case 1:
	slice_bath_gridline();
	break;
    case 2:
	slice_points();
	break;
    }
    XtUnmanageChild(slice_frame);
}

void slice_bygrid(int gno, double x1, double y1, double x2, double y2, double x3, double y3)
{
    int i, ind, nmel, nmnp;
    double get_depth_element(), bathx, xr, x, y, xp, yp, ang = atan2(y2 - y3, x2 - x3);
    char buf[256];
    FILE *fpr = fopen(grid1dfile, "r");
    FILE *fp = fopen(slicefile, "w");
    if (fpr == NULL) {
	errwin("Error opening 1d grid file for input");
	slice_back = 0;
	return;
    }
    if (fp == NULL) {
	errwin("Error opening file for output");
	slice_back = 0;
	return;
    }
    fgets(buf, 255, fpr);
    fgets(buf, 255, fpr);
    sscanf(buf, "%d %d", &nmel, &nmnp);
    fprintf(fp, "Slice from xmgredit\n");
    if (slice_out == 1) {
	fprintf(fp, "%d\n", nmnp);
	for (i = 0; i < nmnp; i++) {
	    fgets(buf, 255, fpr);
	    sscanf(buf, "%lf", &xr);
	    xp = xr * cos(ang);
	    yp = xr * sin(ang);
	    x = x3 + xp;
	    y = y3 + yp;
	    box(x, y);
	    find_element(curgrid, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(curgrid, ind, x, y);
		fprintf(fp, "%d %lf %lf %lf\n", i + 1, x, y, bathx);
	    }
	}
    } else if (slice_out == 0) {
	fprintf(fp, "%d %d\n", nmnp - 1, nmnp);
	for (i = 0; i < nmnp; i++) {
	    fgets(buf, 255, fpr);
	    sscanf(buf, "%lf", &xr);
	    xp = xr * cos(ang);
	    yp = xr * sin(ang);
	    x = x3 + xp;
	    y = y3 + yp;
	    box(x, y);
	    find_element(curgrid, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(curgrid, ind, x, y);
		fprintf(fp, "%lf %lf %lf\n", xr, bathx, 1.0);
	    }
	}
    } else {
    }
    for (i = 0; i < nmnp - 1; i++) {
	fprintf(fp, "2 %d %d\n", i + 1, i + 2);
    }

    fclose(fpr);
    fclose(fp);
    slice_back = 0;
}

void get_slice(int gno, double x1, double y1, double x2, double y2, int which)
{
    int i, ind, last = 0;
    double get_depth_element(), bathx, x, y, dx = x2 - x1, dy = y2 - y1;
    FILE *fp;
    char buf[256];
    if (which) {
	fp = fopen(slicefile, "w");
    } else {
	fp = fopen(slicebackfile, "w");
	slice_back = 0;
	last = 1;
    }
    if (fp == NULL) {
	errwin("Error opening file for output");
	return;
    }
    fprintf(fp, "# Slice from xmgredit\n");
    if (slice_out == 1) {
	fprintf(fp, "%d\n", npts_to_slice);
	for (i = 0; i < npts_to_slice; i++) {
	    x = x1 + i * dx / (npts_to_slice - 1);
	    y = y1 + i * dy / (npts_to_slice - 1);
	    find_element(gno, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(gno, ind, x, y);
		fprintf(fp, "%d %lf %lf %lf\n", i + 1, x, y, bathx);
		fprintf(fp, "%lf %lf\n", hypot(x1 - x, y1 - y), bathx);
	    }
	}
    } else if (slice_out == 0) {
	fprintf(fp, "%d %d\n", npts_to_slice - 1, npts_to_slice);
	for (i = 0; i < npts_to_slice; i++) {
	    x = x1 + i * dx / (npts_to_slice - 1);
	    y = y1 + i * dy / (npts_to_slice - 1);
	    find_element(gno, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(gno, ind, x, y);
		fprintf(fp, "%lf %lf %lf\n", x, bathx, 1.0);
	    }
	}
	for (i = 0; i < npts_to_slice - 1; i++) {
	    fprintf(fp, "2 %d %d\n", i + 1, i + 2);
	}
    } else {
	for (i = 0; i < npts_to_slice; i++) {
	    x = x1 + i * dx / (npts_to_slice - 1);
	    y = y1 + i * dy / (npts_to_slice - 1);
	    find_element(gno, x, y, &ind);
	    if (ind >= 0) {
		bathx = get_depth_element(gno, ind, x, y);
		fprintf(fp, "%lf %lf\n", hypot(x1 - x, y1 - y), bathx);
	    }
	}
    }
    fclose(fp);
    if (last) {
	sprintf(buf, "xmgr %s %s &", slicefile, slicebackfile);
	system(buf);
    }
}

void do_set_slice_grid(void)
{
    XtManageChild(slice_frame);
}

void slice_done_proc(void)
{
    XtUnmanageChild(slice_frame);
}

void create_slice_popup(void)
{
    Widget bt, rc;
    if (slice_frame) {
	XtRaise(slice_frame);
	return;
    }
    slice_frame = XmCreateDialogShell(app_shell, "Slice", NULL, 0);
    slice_panel = XmCreateRowColumn(slice_frame, "rc", NULL, 0);

    slice_method = (Widget *) CreatePanelChoice1(slice_panel, "Points by:",
				4, "Line", "By 1d grid", "Select", 0, NULL);

    grid_to_read = CreateTextItem2(slice_panel,
				   20, "Read 1d grid from: ");
    slice_npts_item = CreateTextItem2(slice_panel,
				      10, "Number of points: ");

    slice_back_item = XtVaCreateManagedWidget("Slice background also",
				     xmToggleButtonWidgetClass, slice_panel,
					      NULL);
    slice_output = (Widget *) CreatePanelChoice1(slice_panel, "Output as:",
						 4,
						 "1d grid",
						 "Bathymetry",
						 "ACE/gr XY format",
						 0,
						 NULL);

    slice_file = CreateTextItem2(slice_panel,
				 30, "File to write:");

    slice_file2 = CreateTextItem2(slice_panel,
				  30, "Background slice to: ");

    rc = XmCreateRowColumn(slice_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) slice_done_proc, NULL);

    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_slice_proc, NULL);
    XtManageChild(rc);

    XtManageChild(slice_panel);
    XtManageChild(slice_frame);
}
