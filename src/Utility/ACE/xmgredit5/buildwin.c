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

 * popup for spreading points on the grid
 *
 */

#ifndef lint
static char RCSid[] = "$Id: buildwin.c,v 1.4 2011/09/14 17:44:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

int npts_to_place = 0;
int spreadflag = 0;
int rseed = 100;
double spread_dx, spread_dy;
extern double mindist;

static Widget mindist_frame;
static Widget mindist_panel;
static Widget mindist_item;

static Widget triang_frame;
static Widget triang_panel;
static Widget *triang_method_item;
static Widget triang_mindist_item;
static Widget triangleframe;
static int doboundary = 0;
static Widget doboundarytoggle;
static Widget maxareatoggle;
static Widget maxareatext;
static Widget minangletoggle;
static Widget minangletext;
static Widget maxsteinertoggle;
static Widget maxsteinertext;
static Widget exactarithmetictoggle;

static Widget spreadr_frame;
static Widget spreadr_panel;

static Widget build_frame;
static Widget build_panel;
static Widget build_centx_item;
static Widget build_centy_item;
static Widget build_maxradius_item;
static Widget build_nradii_item;
static Widget build_initdiv_item;
static Widget build_divinc_item;
static Widget build_grid_item;
static Widget build_dx_item;
static Widget build_dy_item;
static Widget build_grid_label;
static Widget build_dx_label;
static Widget build_dy_label;
static Widget build_region_item;
static Widget *build_spread_item;

static Widget genfd_frame;
static Widget genfd_panel;
static void do_genfd_proc(void);

void do_select_build_region();
void do_clear_region();		/* defined in editwin.c */
void delete_build_outside_region();
void delete_build_inside_region();

void update_mindist(void);

char *panel_getstr_value();

/*
 * function for printing to the text item in the
 * triagulation popup
 */
static Widget triang_text_item;

#include <stdarg.h>

void my_printf(char *fmt,...)
{
    char buf[1024];
    va_list argp;
    va_start(argp, fmt);
    vsprintf(buf, fmt, argp);
    XmTextInsert(triang_text_item, XmTextGetLastPosition(triang_text_item), buf);
    va_end(argp);
}

/* callback for loading centers of elements to build points */
void do_load_centers(void)
{
    load_centers(0, 0);
}

static void do_circular_points_proc(void)
{
    char buf[256];
    double cx, cy, r;
    int nr, na, nainc;
    strcpy(buf, panel_getstr_value(build_centx_item));
    sscanf(buf, "%lf", &cx);
    strcpy(buf, panel_getstr_value(build_centy_item));
    sscanf(buf, "%lf", &cy);
    strcpy(buf, panel_getstr_value(build_maxradius_item));
    sscanf(buf, "%lf", &r);
    strcpy(buf, panel_getstr_value(build_nradii_item));
    sscanf(buf, "%d", &nr);
    strcpy(buf, panel_getstr_value(build_initdiv_item));
    sscanf(buf, "%d", &na);
    strcpy(buf, panel_getstr_value(build_divinc_item));
    sscanf(buf, "%d", &nainc);
    gen_circ(curbuild, cx, cy, r, nr, na, nainc);
}

void do_set_build_grid(void)
{
    XtRaise(build_frame);
}

void do_set_mindist(void)
{
    char buf[256];
    strcpy(buf, panel_getstr_value(mindist_item));
    sscanf(buf, "%lf", &mindist);
}

void create_build_popup(void)
{
    Widget bt, rc;
    XmString str;
    Arg a[2];
    int ac = 0;
    spreadflag = 0;
    if (build_frame) {
	XtRaise(build_frame);
	return;
    }
    build_frame = XmCreateDialogShell(app_shell,
			       "Spread build points in circle", NULL, 0);
    build_panel = XmCreateRowColumn(build_frame, "rc", NULL, 0);
    build_centx_item = CreateTextItem2(build_panel, 10, "Center X:");
    build_centy_item = CreateTextItem2(build_panel, 10, "Center Y:");
    build_maxradius_item = CreateTextItem2(build_panel, 10, "Maximum radius:");
    build_nradii_item = CreateTextItem2(build_panel, 10, "# of radii:");
    build_initdiv_item = CreateTextItem2(build_panel, 10,
					 "# of divisions in 1st circle:");
    build_divinc_item = CreateTextItem2(build_panel, 10,
				    "# of additional arcs\nper circle:");

    rc = XmCreateRowColumn(build_panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_circular_points_proc, NULL);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, build_frame);
    XtManageChild(rc);
    XtManageChild(build_panel);
    XtManageChild(build_frame);
}

void place_build(void)
{
    write_mode_str("Placing building points, use left button to place, right button to stop");
    set_action(0);
    set_action(PLACE_BUILD);
}

void delete_build(void)
{
    set_action(0);
    write_mode_str("Delete build points, use left button to mark, middle to register");
    set_action(DELETE_BUILD);
}

void move_build(void)
{
    set_action(0);
    write_mode_str("Move build points, use left button to mark, then move");
    set_action(MOVE_BUILD1ST);
}

/*
 * use Fortune's sweepline code
 */
void triangulate_build(void)
{
    write_mode_str("Triangularizing the current set of grid definition points");
    if (mindist != 0.0) {
	elim_build_dupes(curbuild, mindist);
    }
    triang(curbuild);
    write_mode_str("Loading triangularized points to editable grid, please wait");
    load_grid(curgrid, curbuild);
    do_drawgrid();
    write_mode_str(NULL);
}

/*
 * use Shewchuck's Triangle 
 */
#ifdef DO_TRIANGLE
void triangulate_triangle(char *args)
{
    write_mode_str("Triangularizing the current set of build points");
    tritest(curbuild, curgrid, doboundary, args);
    do_drawgrid();
    write_mode_str(NULL);
}

#endif

/*
 * Make the nodes of the grid build points
 */
void loadgrid_build(void)
{
    write_mode_str("Loading building points from editable grid");
    load_from_grid(curgrid, curbuild);
    drawbuild(curbuild);
    write_mode_str(NULL);
}

void mergegrid_build(void)
{
    int i;
    if (yesno("Merge? Are you sure?", "Press Yes or No", "Yes", "No")) {
	write_mode_str("Merging (building pts. + editable grid = building pts.)");
	init_qadd_build(curbuild);
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    qadd_build(curbuild, grid[curgrid].xord[i], grid[curgrid].yord[i], grid[curgrid].depth[i]);
	}
	flush_qadd_build(curbuild);
	do_drawgrid();
    }
}

void mergebound_build(void)
{
    int i, j, ib;
    if (yesno("Merge? Are you sure?", "Press Yes or No", "Yes", "No")) {
	write_mode_str("Adding boundary definition to grid definition points");
	init_qadd_build(curbuild);
	for (i = 0; i < grid[curgrid].nbounds; i++) {
	    ib = grid[curgrid].boundaries[i];
	    for (j = 0; j < boundary[ib].nbpts; j++) {
		qadd_build(curbuild, boundary[ib].boundx[j], boundary[ib].boundy[j], 1.0);
	    }
	}
	flush_qadd_build(curbuild);
	do_drawgrid();
    }
}

void toggle_build(void)
{
    write_mode_str("Toggling display of grid definition points");
    display_build_points = !display_build_points;
    do_drawgrid();
}

void clear_build(void)
{
    write_mode_str("Clearing grid definition points");
    Free_build(curbuild);
    do_drawgrid();
}

Widget autob_frame, autob_panel;
static Widget crit1_item, crit2_item;
static Widget label1_item, label2_item;
static Widget autob_dx_label;
static Widget autob_dx_item;
static Widget autob_loop_label;
static Widget autob_loop_item;
static Widget autob_x1_label;
static Widget autob_x1_item;
static Widget autob_x2_label;
static Widget autob_x2_item;
static Widget autob_y1_label;
static Widget autob_y1_item;
static Widget autob_y2_label;
static Widget autob_y2_item;
static Widget npoints_item;
static Widget showgrid_item;

static Widget prev[4];
static int nprev;
static int autotype = 2;

static double autob_tol = 100000.0;
static double courtime = 600.0, maxcour = 1.0;
static double dimwtime = 4.1 * 3600.0, mindimw = 25.0;

static int autob_loop = 50;
static double autob_x1, autob_x2, autob_y1, autob_y2;
static double autob_dx, autob_dy;

static Widget cbound_item;
static Widget cboundmin_item;
static Widget celevmin_item;
static Widget *method_item;
static Widget file_func_item;

double cboundmin = 2.0;
double celevmin = 2.0;

void set_autotype_proc(Widget w, int cd)
{
    XmString string;
    Arg al;
    char s[256];
    switch (autotype = cd) {
    case 0:
	string = XmStringCreateLtoR("Maximum area", charset);
	XtSetArg(al, XmNlabelString, string);
	XtSetValues(label1_item, &al, 1);
	XmStringFree(string);
	XtUnmanageChild(crit2_item);
	XtUnmanageChild(label2_item);
	autotype = 0;
	break;
    case 1:
	sprintf(s, "%lg", maxcour);
	XmTextSetString(crit1_item, s);
	sprintf(s, "%lg", courtime);
	XmTextSetString(crit2_item, s);
	string = XmStringCreateLtoR("Maximum Courant number:", charset);
	XtSetArg(al, XmNlabelString, string);
	XtSetValues(label1_item, &al, 1);
	XmStringFree(string);
	string = XmStringCreateLtoR("Time step:", charset);
	XtSetArg(al, XmNlabelString, string);
	XtSetValues(label2_item, &al, 1);
	XmStringFree(string);
	XtManageChild(crit1_item);
	XtManageChild(label1_item);
	XtManageChild(crit2_item);
	XtManageChild(label2_item);
	autotype = 1;

	break;
    case 2:
	sprintf(s, "%lg", mindimw);
	XmTextSetString(crit1_item, s);
	sprintf(s, "%lg", dimwtime / 3600.0);
	XmTextSetString(crit2_item, s);
	string = XmStringCreateLtoR("Minimum dim. wavelength:", charset);
	XtSetArg(al, XmNlabelString, string);
	XtSetValues(label1_item, &al, 1);
	XmStringFree(string);
	string = XmStringCreateLtoR("Shortest wave (hours):", charset);
	XtSetArg(al, XmNlabelString, string);
	XtSetValues(label2_item, &al, 1);
	XmStringFree(string);
	XtManageChild(crit1_item);
	XtManageChild(label1_item);
	XtUnmanageChild(crit2_item);
	XtUnmanageChild(label2_item);
	XtManageChild(crit2_item);
	XtManageChild(label2_item);
	autotype = 2;
	break;
    case 3:			/* external function */
	XtUnmanageChild(crit1_item);
	XtUnmanageChild(label1_item);
	XtUnmanageChild(crit2_item);
	XtUnmanageChild(label2_item);
	break;
    }
}

void update_autob(void)
{
    char s[256];
    XmString string;
    int nx, ny;
    double wx, wy;
    set_autotype_proc(NULL, autotype);
    switch (autotype) {
    case 0:
	break;
    case 1:
	sprintf(s, "%lg", maxcour);
	XmTextSetString(crit1_item, s);
	sprintf(s, "%lg", courtime);
	XmTextSetString(crit2_item, s);
	break;
    case 2:
	sprintf(s, "%lg", mindimw);
	XmTextSetString(crit1_item, s);
	sprintf(s, "%lg", dimwtime / 3600.0);
	XmTextSetString(crit2_item, s);
	break;
    case 3:
	break;
    }
    sprintf(s, "%lg", cboundmin);
    XmTextSetString(cboundmin_item, s);
    sprintf(s, "%lg", celevmin);
    XmTextSetString(celevmin_item, s);
    sprintf(s, "%d", autob_loop);
    XmTextSetString(autob_loop_item, s);
    sprintf(s, "%.2lf", autob_dx);
    XmTextSetString(autob_dx_item, s);
    sprintf(s, "%lg", autob_x1);
    XmTextSetString(autob_x1_item, s);
    sprintf(s, "%lg", autob_x2);
    XmTextSetString(autob_x2_item, s);
    sprintf(s, "%lg", autob_y1);
    XmTextSetString(autob_y1_item, s);
    sprintf(s, "%lg", autob_y2);
    XmTextSetString(autob_y2_item, s);
    wx = autob_x2 - autob_x1;
    wy = autob_y2 - autob_y1;
    nx = (int) (wx / autob_dx + 1);
    ny = (int) (wy / autob_dx + 1);
    sprintf(s, "%d points", nx * ny);
    string = XmStringCreateLtoR(s, charset);
    XtVaSetValues(npoints_item,
		  XmNlabelString, string,
		  NULL);
    XmStringFree(string);
}

void draw_auxgrid(double ax1, double ay1, double ax2, double ay2, double dx)
{
    int i, j, nx, ny;
    double wx, wy;
    wx = ax2 - ax1;
    wy = ay2 - ay1;
    nx = (int) (wx / dx + 0.5 + 1);
    ny = (int) (wy / dx + 0.5 + 1);
    setcolor(9);
    for (i = 0; i < ny; i++) {
	my_move2(ax1, ay1 + dx * i);
	my_draw2(ax2, ay1 + dx * i);
    }
    for (i = 0; i < nx; i++) {
	my_move2(ax1 + dx * i, ay1);
	my_draw2(ax1 + dx * i, ay2);
    }
    flush_pending();
/*
   ax1 + dx * i, ay1 + dy * j;
   for (j = 0; j < ny; j++) {
   if ((j % 20 == 0) && check_action()) {
   extern int cancel_flag;
   cancel_action(0);
   if (cancel_flag) {
   return;
   }
   }
   }
   }
 */
}

void autob_preview_proc(void)
{
    char *ss;
    char s[256];
    XmString string;
    double atof(const char *);
    int cbound, crit = 0;
    int al, dg, nx, ny;
    double dx, ax1, ax2, ay1, ay2, wx, wy;
    al = atoi(ss = XmTextGetString(autob_loop_item));
    XtFree(ss);
    dx = atof(ss = XmTextGetString(autob_dx_item));
    XtFree(ss);
    ax1 = atof(ss = XmTextGetString(autob_x1_item));
    XtFree(ss);
    ax2 = atof(ss = XmTextGetString(autob_x2_item));
    XtFree(ss);
    ay1 = atof(ss = XmTextGetString(autob_y1_item));
    XtFree(ss);
    ay2 = atof(ss = XmTextGetString(autob_y2_item));
    XtFree(ss);
    if (dx == 0.0) {
	errwin("Grid delta == 0");
	return;
    }
    wx = ax2 - ax1;
    wy = ay2 - ay1;
    nx = (int) (wx / dx + 1);
    ny = (int) (wy / dx + 1);
    sprintf(s, "%d points", nx * ny);
    string = XmStringCreateLtoR(s, charset);
    XtVaSetValues(npoints_item,
		  XmNlabelString, string,
		  NULL);
    XmStringFree(string);
    dg = XmToggleButtonGetState(showgrid_item);
    if (dg) {
	draw_auxgrid(ax1, ay1, ax2, ay2, dx);
    }
}

void autob_accept_proc(Widget w, int cd)
{
    char *s;
    double atof(const char *);
    int cbound, crit = 0, use_vel, build_region;
    use_vel = XmToggleButtonGetState(file_func_item);
    build_region = XmToggleButtonGetState(build_region_item);

    autob_loop = atoi(s = XmTextGetString(autob_loop_item));
    XtFree(s);
    autob_dy = autob_dx = atof(s = XmTextGetString(autob_dx_item));
    XtFree(s);
    autob_x1 = atof(s = XmTextGetString(autob_x1_item));
    XtFree(s);
    autob_x2 = atof(s = XmTextGetString(autob_x2_item));
    XtFree(s);
    autob_y1 = atof(s = XmTextGetString(autob_y1_item));
    XtFree(s);
    autob_y2 = atof(s = XmTextGetString(autob_y2_item));
    XtFree(s);
    if (autob_dx == 0.0) {
	errwin("Grid delta == 0, operation cancelled");
	return;
    }
    if (build[curbuild].nbuild) {
	if (yesno("Clear current set of build points?", "Press Yes or No", "Yes", "No")) {
	    Free_build(curbuild);
	    XmUpdateDisplay(app_shell);
	}
    }
    cbound = XmToggleButtonGetState(cbound_item);
    cboundmin = atof(panel_getstr_value(cboundmin_item));
    celevmin = atof(panel_getstr_value(celevmin_item));
    switch (autotype) {
    case 0:
	autob_tol = atof(s = XmTextGetString(crit1_item));
	XtFree(s);
	if (cbound) {
	    cbound = do_newbound(0, autotype, autob_tol, 1.0, 0.0, MAXGRIDS);
	    if (!cbound) {
		errwin("No boundary formed given current\nsettings, operation cancelled");
		return;
	    }
	}
	break;
    case 1:
	maxcour = atof(s = XmTextGetString(crit1_item));
	XtFree(s);
	courtime = atof(s = XmTextGetString(crit2_item));
	XtFree(s);
	if (!use_vel) {
	    autob_tol = M_PI * 9.8 * courtime * courtime / (maxcour * maxcour * 4);
	} else {
	    autob_tol = M_PI * courtime * courtime / (maxcour * maxcour * 4);
	}
	if (cbound) {
	    cbound = do_newbound(0, autotype, maxcour,
		    courtime, autob_tol, use_vel ? buildgrid : backgrid);
	    if (!cbound) {
		errwin("No boundary formed given current\nsettings, operation cancelled");
		return;
	    }
	}
	break;
    case 2:
	mindimw = atof(s = XmTextGetString(crit1_item));
	XtFree(s);
	dimwtime = atof(s = XmTextGetString(crit2_item)) * 3600.0;
	XtFree(s);
	if (!use_vel) {
	    autob_tol = M_PI * 9.8 * dimwtime * dimwtime / (mindimw * mindimw * 4);
	} else {
	    autob_tol = M_PI * dimwtime * dimwtime / (mindimw * mindimw * 4);
	}
	if (cbound) {
	    cbound = do_newbound(0, autotype, mindimw,
		    dimwtime, autob_tol, use_vel ? buildgrid : backgrid);
	    if (!cbound) {
		errwin("No boundary formed given current\nsettings, operation cancelled");
		return;
	    }
	}
	break;
    case 3:
	break;
    }
    auto_spread(curgrid, autob_x1, autob_y1, autob_x2, autob_y2, autob_dx, autob_dy,
		autotype, autob_loop, autob_tol, use_vel ? buildgrid : backgrid, build_region);
    do_drawgrid();
}

void autob_defaults_proc(void)
{
    setlimits_grid(MAXGRIDS);
    autob_loop = 50;
    autob_x1 = floor(grid[MAXGRIDS].xmin);
    autob_x2 = ceil(grid[MAXGRIDS].xmax);
    autob_y1 = floor(grid[MAXGRIDS].ymin);
    autob_y2 = ceil(grid[MAXGRIDS].ymax);
    autob_dx = autob_dy = (autob_x2 - autob_x1) / 100.0;
    update_autob();
}

static Widget rdata_dialog;
static void read_vel_proc(void);

void autob_file_proc(void)
{
    int i;
    Widget lab, rc, fr, rb, w[10];
    Arg a;

    if (rdata_dialog) {
	XtManageChild(rdata_dialog);
	return;
    }
    rdata_dialog = XmCreateFileSelectionDialog(app_shell, "Read velocity data", NULL, 0);
    XtSetArg(a, XmNdialogTitle, XmStringCreateLtoR("Read velocity data", charset));
    XtSetValues(rdata_dialog, &a, 1);
    XtAddCallback(rdata_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, rdata_dialog);
    XtAddCallback(rdata_dialog, XmNokCallback, (XtCallbackProc) read_vel_proc, 0);
    XtManageChild(rdata_dialog);
}

static void read_vel_proc(void)
{
    Arg args;
    XmString list_item;
    char *s;
    int i;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rdata_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    set_wait_cursor(rdata_dialog);
    if (readgrid(buildgrid, s)) {
	for (i = 0; i < grid[buildgrid].nmnp; i++) {
	    grid[buildgrid].depth[i] = grid[buildgrid].depth[i] * grid[buildgrid].depth[i];
	}
    }
    unset_wait_cursor(rdata_dialog);
}

void set_dimwtype_proc(void)
{
}

void set_courtype_proc(void)
{
}

void create_autob_popup(void)
{
    Widget bt, lab, rc, rc2, rc3, fr, fr1, fr2, rb, w[10], sep;
    XmString str;
    Arg a[5];
    int i, ac = 0;
    if (autob_frame) {
	XtRaise(autob_frame);
	update_autob();
	return;
    }
    autob_frame = XmCreateDialogShell(app_shell, "Automatic placement", NULL, 0);
    handle_close(autob_frame);
    autob_panel = XmCreateForm(autob_frame, "autob_panel", NULL, 0);
    XtVaSetValues(autob_panel, XmNorientation, XmHORIZONTAL, NULL);

    fr = XmCreateFrame(autob_panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    lab = XmCreateLabel(rc, "Distance between build points set by:", NULL, 0);
    XtManageChild(lab);
    method_item = CreatePanelChoice1(rc, " ",
				     5,
				     "Maximum area",
				     "Maximum Courant number",
				     "Min. dimensionless wavelength",
				     "*External function",
				     NULL, 0);
    for (i = 0; i < 4; i++) {
	XtAddCallback(method_item[i + 2], XmNactivateCallback,
		      (XtCallbackProc) set_autotype_proc, (XtPointer) ((long) i));
    }
    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    label1_item = XmCreateLabel(rc2, "Min. dimensionless\nwavelength:", NULL, 0);
    crit1_item = XmCreateText(rc2, "crit1_text", NULL, 0);
    XtVaSetValues(crit1_item, XmNcolumns, 10, NULL);
    XtManageChild(label1_item);
    XtManageChild(crit1_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    label2_item = XmCreateLabel(rc2, "Shortest wave (hours):", NULL, 0);
    crit2_item = XmCreateText(rc2, "crit2_text", NULL, 0);
    XtVaSetValues(crit2_item, XmNcolumns, 10, NULL);


    XtManageChild(label2_item);
    XtManageChild(crit2_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    file_func_item = XmCreateToggleButton(rc2, "Use velocity read from file", NULL, 0);
    XtManageChild(file_func_item);
    bt = XtVaCreateManagedWidget("File...", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autob_file_proc, NULL);

    XtManageChild(rc2);

    sep = XmCreateSeparatorGadget(rc, "sep", NULL, 0);
    XtManageChild(sep);
    cbound_item = XmCreateToggleButton(rc, "Create boundary", NULL, 0);
    XtManageChild(cbound_item);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Minimum depth at boundary", NULL, 0);
    cboundmin_item = XmCreateText(rc2, "mindepth", NULL, 0);
    XtVaSetValues(cboundmin_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(cboundmin_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Minimum depth in domain", NULL, 0);
    celevmin_item = XmCreateText(rc2, "mindepth", NULL, 0);
    XtVaSetValues(celevmin_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(celevmin_item);
    XtManageChild(rc2);

    XtManageChild(rc);
    XtManageChild(fr);


    fr1 = XmCreateFrame(autob_panel, "fr1", NULL, 0);
    rc = XmCreateRowColumn(fr1, "rc", NULL, 0);
    lab = XmCreateLabelGadget(rc, "Auxillary grid parameters:", NULL, 0);
    XtManageChild(lab);
    XtVaSetValues(rc, XmNorientation, XmVERTICAL,
		  XmNpacking, XmPACK_COLUMN,
		  NULL);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_dx_label = XmCreateLabel(rc2, "Grid delta:  ", NULL, 0);
    autob_dx_item = XmCreateText(rc2, "dx", NULL, 0);
    XtVaSetValues(autob_dx_item, XmNcolumns, 10, NULL);
    XtManageChild(autob_dx_label);
    XtManageChild(autob_dx_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_loop_label = XmCreateLabel(rc2, "Maximum loop:  ", NULL, 0);
    autob_loop_item = XmCreateText(rc2, "loop", NULL, 0);
    XtVaSetValues(autob_loop_item, XmNcolumns, 7, NULL);
    XtManageChild(autob_loop_label);
    XtManageChild(autob_loop_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_x1_label = XmCreateLabel(rc2, "Xmin:  ", NULL, 0);
    autob_x1_item = XmCreateText(rc2, "xmin", NULL, 0);
    XtVaSetValues(autob_x1_item, XmNcolumns, 10, NULL);
    XtManageChild(autob_x1_label);
    XtManageChild(autob_x1_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_x2_label = XmCreateLabel(rc2, "Xmax:  ", NULL, 0);
    autob_x2_item = XmCreateText(rc2, "xmax", NULL, 0);
    XtVaSetValues(autob_x2_item, XmNcolumns, 10, NULL);
    XtManageChild(autob_x2_label);
    XtManageChild(autob_x2_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_y1_label = XmCreateLabel(rc2, "Ymin:  ", NULL, 0);
    autob_y1_item = XmCreateText(rc2, "ymin", NULL, 0);
    XtVaSetValues(autob_y1_item, XmNcolumns, 10, NULL);
    XtManageChild(autob_y1_label);
    XtManageChild(autob_y1_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    autob_y2_label = XmCreateLabel(rc2, "Ymax:  ", NULL, 0);
    autob_y2_item = XmCreateText(rc2, "ymax", NULL, 0);
    XtVaSetValues(autob_y2_item, XmNcolumns, 10, NULL);
    XtManageChild(autob_y2_label);
    XtManageChild(autob_y2_item);
    XtManageChild(rc2);

    npoints_item = XmCreateLabel(rc, "npoints", NULL, 0);
    XtVaSetValues(npoints_item, XmNalignment, XmALIGNMENT_BEGINNING, NULL);
    XtManageChild(npoints_item);
    showgrid_item = XmCreateToggleButton(rc, "Display auxillary grid", NULL, 0);
    XtManageChild(showgrid_item);
    build_region_item = XmCreateToggleButton(rc, "Restrict to region", NULL, 0);
    XtManageChild(build_region_item);

    XtManageChild(rc);
    XtManageChild(fr1);

    fr2 = XmCreateFrame(autob_panel, "fr", NULL, 0);
    rc2 = XmCreateRowColumn(fr2, "rc", a, ac);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autob_accept_proc, NULL);
    bt = XtVaCreateManagedWidget("Default grid", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autob_defaults_proc, NULL);
    bt = XtVaCreateManagedWidget("Preview", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autob_preview_proc, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, autob_frame);
    XtManageChild(rc2);
    XtManageChild(fr2);

    XtVaSetValues(fr,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(fr1,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_WIDGET,
		  XmNleftWidget, fr,
		  NULL);
    XtVaSetValues(fr2,
		  XmNbottomAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);

    autob_loop = 50;
    autob_x1 = floor(grid[MAXGRIDS].xmin);
    autob_x2 = ceil(grid[MAXGRIDS].xmax);
    autob_y1 = floor(grid[MAXGRIDS].ymin);
    autob_y2 = ceil(grid[MAXGRIDS].ymax);
    autob_dx = autob_dy = (autob_x2 - autob_x1) / 100.0;
    SetChoice(method_item, 2);
    update_autob();
    XtManageChild(autob_panel);
    XtManageChild(autob_frame);
}

void update_mindist(void)
{
    char buf[256];
    if (mindist_frame) {
	sprintf(buf, "%lf", mindist);
	panel_setstr_value(mindist_item, buf);
    }
}

void create_mindist_popup(void)
{
    Widget bt, rc;
    spreadflag = 0;
    if (!mindist_frame) {
	mindist_frame = XmCreateDialogShell(app_shell, "Set minimum distance", NULL, 0);
	mindist_panel = XmCreateRowColumn(mindist_frame, "rc", NULL, 0);
	mindist_item = CreateTextItem2(mindist_panel, 20, "Minimum distance:");

	rc = XmCreateRowColumn(mindist_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL,
		      XmNpacking, XmPACK_TIGHT, NULL);
	bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_set_mindist, NULL);
	bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, mindist_frame);
	XtManageChild(rc);

	XtManageChild(mindist_panel);
    }
    update_mindist();
    XtRaise(mindist_frame);
}

void set_recttype(Widget w, int cd)
{
    if (cd == 2) {
	XtSetSensitive(build_grid_label, True);
	XtSetSensitive(build_grid_item, True);
	XtSetSensitive(build_dx_label, False);
	XtSetSensitive(build_dx_item, False);
	XtSetSensitive(build_dy_label, False);
	XtSetSensitive(build_dy_item, False);
    } else {
	XtSetSensitive(build_grid_label, False);
	XtSetSensitive(build_grid_item, False);
	XtSetSensitive(build_dx_label, True);
	XtSetSensitive(build_dx_item, True);
	XtSetSensitive(build_dy_label, True);
	XtSetSensitive(build_dy_item, True);
    }
}

void update_spreadr(void)
{
    char buf[256];
    if (spreadr_frame) {
	set_recttype(NULL, GetChoice(build_spread_item));
    }
}

static void do_spread_points_proc(void)
{
    char buf[256];

    spreadflag = GetChoice(build_spread_item);
    strcpy(buf, panel_getstr_value(build_grid_item));
    sscanf(buf, "%d", &npts_to_place);
    strcpy(buf, panel_getstr_value(build_dx_item));
    sscanf(buf, "%lf", &spread_dx);
    strcpy(buf, panel_getstr_value(build_dy_item));
    sscanf(buf, "%lf", &spread_dy);
    srand48((long) rseed);
    set_action(0);
    set_action(ADD_BUILD1);
}

static void do_spreadflag_proc(int value)
{
    extern int spreadflag;
    spreadflag = value;
}

void create_spreadr_popup(void)
{
    Widget bt, rc;
    int i;
    spreadflag = 0;
    if (!spreadr_frame) {
	spreadr_frame = XmCreateDialogShell(app_shell, "Spread points in rectangles", NULL, 0);
	spreadr_panel = XmCreateRowColumn(spreadr_frame, "rc", NULL, 0);
	build_spread_item = (Widget *) CreatePanelChoice1(spreadr_panel,
						      "Generation mode:",
							  6,
							  "Rectangular",
						    "Rectangular offset",
						      "Random (uniform)",
						   "Rotated Rectangular",
					    "Rotated Rectangular offset",
							  0,
							  NULL);
	for (i = 0; i < 5; i++) {
	    XtAddCallback(build_spread_item[i + 2], XmNactivateCallback,
			  (XtCallbackProc) set_recttype, (XtPointer)((long) i));
	}

	build_grid_item = CreateTextItem3(spreadr_panel,
		     10, "How many points to place:", &build_grid_label);
	build_dx_item = CreateTextItem3(spreadr_panel,
				   10, "Spacing in X:", &build_dx_label);
	build_dy_item = CreateTextItem3(spreadr_panel,
				   10, "Spacing in Y:", &build_dy_label);

	rc = XmCreateRowColumn(spreadr_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL,
		      XmNpacking, XmPACK_TIGHT, NULL);
	bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_spread_points_proc, NULL);

	bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, spreadr_frame);
	XtManageChild(rc);
	XtManageChild(spreadr_panel);
    }
    update_spreadr();
    XtRaise(spreadr_frame);
}

/*
 * Section for Del. triangulation
 */

/*
 * Set visibility for various items depeneding on tringulation method
 */
void set_triangtype_proc(void)
{
/*
   switch (GetChoice(triang_method_item)) {
   case 0:
   XtSetSensitive(triangleframe, False);
   break;
   case 1:
   XtSetSensitive(triangleframe, True);
   break;
   }
 */
}

void update_triang(void)
{
    char buf[256];
    if (triang_frame) {
	sprintf(buf, "%lg", mindist);
	panel_setstr_value(triang_mindist_item, buf);
	set_triangtype_proc();
    }
}

void do_triang_proc(void)
{
    char buf[256];
    char args[256];
    int doexacta, dominangle, domaxarea, domaxsteiner;
    double minangle;
    double maxarea;
    double maxsteiner;
    strcpy(buf, panel_getstr_value(triang_mindist_item));
    sscanf(buf, "%lf", &mindist);
#ifdef DO_TRIANGLE
    doboundary = XmToggleButtonGetState(doboundarytoggle);
    doexacta = XmToggleButtonGetState(exactarithmetictoggle);
    dominangle = XmToggleButtonGetState(minangletoggle);
    domaxarea = XmToggleButtonGetState(maxareatoggle);
    domaxsteiner = XmToggleButtonGetState(maxareatoggle);
#endif
    switch (GetChoice(triang_method_item)) {
    case 0:
	my_printf("Triangulating build points using Fortune's Sweepline with a minimum distance of %lg\n", mindist);
	set_wait_cursor(triang_frame);
	triangulate_build();
	unset_wait_cursor(triang_frame);
	break;
#ifdef DO_TRIANGLE
    case 1:
	set_wait_cursor(triang_frame);
	strcpy(args, "p");
	if (dominangle) {
	    minangle = atof(panel_getstr_value(minangletext));
	    sprintf(buf, "q%lf", minangle);
	    strcat(args, buf);
	}
	if (domaxarea) {
	    maxarea = atof(panel_getstr_value(maxareatext));
	    sprintf(buf, "a%lf", maxarea);
	    strcat(args, buf);
	}
	if (doexacta) {
	    strcat(args, "X");
	}
	if (domaxsteiner) {
	    maxsteiner = atof(panel_getstr_value(maxsteinertext));
	    sprintf(buf, "S%lf", maxsteiner);
	    strcat(args, buf);
	}
	triangulate_triangle(args);
	unset_wait_cursor(triang_frame);
	break;
#endif
    }
    my_printf("Done: %d nodes, %d elements\n", grid[curgrid].nmnp, grid[curgrid].nmel);
}

void create_triang_popup(void)
{
    Widget bt, fr, fr2, rc, rc2;
    XmString str;
    Arg a[2];
    int ac = 0, x, y, i;
    if (!triang_frame) {
	XmGetPos(app_shell, 0, &x, &y);
	ac = 0;
	XtSetArg(a[ac], XmNx, x + 300);
	ac++;
	XtSetArg(a[ac], XmNy, y + 300);
	ac++;
	triang_frame = XmCreateDialogShell(app_shell, "Triangulate build points", a, ac);
	triang_panel = XtVaCreateWidget("triang",
				    xmRowColumnWidgetClass, triang_frame,
					NULL);
#ifdef DO_TRIANGLE
	triang_method_item = (Widget *) CreatePanelChoice1(triang_panel,
					      "Triangulation engine:", 3,
						   "Fortune's sweepline",
						  "Shewchuck's Triangle",
							   0, NULL);

	for (i = 0; i < 2; i++) {
	    XtAddCallback(triang_method_item[i + 2], XmNactivateCallback,
		    (XtCallbackProc) set_triangtype_proc, (XtPointer) ((long) i));
	}
#else
	triang_method_item = (Widget *) CreatePanelChoice1(triang_panel,
					      "Triangulation engine:", 2,
						   "Fortune's sweepline",
							   0, NULL);
#endif

	rc2 = XmCreateRowColumn(triang_panel, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	bt = XmCreateLabel(rc2, "Minimum distance:", NULL, 0);
	triang_mindist_item = XmCreateText(rc2, "mindist", NULL, 0);
	XtVaSetValues(triang_mindist_item, XmNcolumns, 20, NULL);
	XtManageChild(bt);
	XtManageChild(triang_mindist_item);
	XtManageChild(rc2);

#ifdef DO_TRIANGLE
	triangleframe = XmCreateFrame(triang_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(triangleframe, "rc", NULL, 0);
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	doboundarytoggle = XmCreateToggleButton(rc2, "Add boundary constraints", NULL, 0);
	XtManageChild(doboundarytoggle);
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	maxareatoggle = XmCreateToggleButton(rc2, "Apply maximum area constraint", NULL, 0);
	XtManageChild(maxareatoggle);
	maxareatext = CreateTextItem2(rc2, 10, "Max area:");
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	minangletoggle = XmCreateToggleButton(rc2, "Apply minimum angle constraint", NULL, 0);
	XtManageChild(minangletoggle);
	minangletext = CreateTextItem2(rc2, 10, "Min angle:");
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	maxsteinertoggle = XmCreateToggleButton(rc2, "Apply max # of Steiner points", NULL, 0);
	XtManageChild(maxsteinertoggle);
	maxsteinertext = CreateTextItem2(rc2, 10, "Max # Steiner points:");
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	exactarithmetictoggle = XmCreateToggleButton(rc2, "Suppress exact arithmetic", NULL, 0);
	XtManageChild(exactarithmetictoggle);
	XtManageChild(rc2);
	XtManageChild(rc);
	XtManageChild(triangleframe);
#endif

	fr2 = XmCreateFrame(triang_panel, "fr", NULL, 0);
	triang_text_item = XmCreateScrolledText(fr2, "triang_text_item", NULL, 0);
	XtVaSetValues(triang_text_item,
		      XmNrows, 10,
		      XmNcolumns, 40,
		      XmNeditMode, XmMULTI_LINE_EDIT,
		      XmNwordWrap, True,
		      NULL);
	XtManageChild(triang_text_item);
	XtManageChild(fr2);

	fr2 = XmCreateFrame(triang_panel, "fr", NULL, 0);
	rc2 = XmCreateRowColumn(fr2, "rc", NULL, 0);
	XtVaSetValues(rc2,
		      XmNorientation, XmHORIZONTAL,
		      XmNpacking, XmPACK_TIGHT,
		      NULL);
	bt = XtVaCreateManagedWidget("Apply", xmPushButtonWidgetClass, rc2,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_triang_proc, NULL);
	bt = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc2,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, triang_frame);
	XtManageChild(rc2);
	XtManageChild(fr2);
	XtManageChild(triang_panel);
    }
    update_triang();
    XtRaise(triang_frame);
}

static Widget genfd_dx_item, genfd_dy_item;
static Widget genfd_nx_item, genfd_ny_item;
static Widget genfd_lx_item, genfd_ly_item;
static Widget genfd_sx_item, genfd_sy_item;
static Widget genfd_rot_item;
static Widget *genfd_type_item;
static Widget genfd_dodxdy_item;
static Widget genfd_donxny_item;
void create_genfd_grid(void);

static int genfd_doit = 1;
static int genfd_signx = 1, genfd_signy = 1;

void update_genfd(void)
{
    char buf[256];
    if (genfd_frame) {
	sprintf(buf, "%.2lf", genfd_dx);
	panel_setstr_value(genfd_dx_item, buf);
	sprintf(buf, "%.2lf", genfd_dy);
	panel_setstr_value(genfd_dy_item, buf);
	sprintf(buf, "%d", genfd_nx);
	panel_setstr_value(genfd_nx_item, buf);
	sprintf(buf, "%d", genfd_ny);
	panel_setstr_value(genfd_ny_item, buf);
	sprintf(buf, "%lf", genfd_sx);
	panel_setstr_value(genfd_sx_item, buf);
	sprintf(buf, "%lf", genfd_sy);
	panel_setstr_value(genfd_sy_item, buf);
	sprintf(buf, "%lf", genfd_lx);
	panel_setstr_value(genfd_lx_item, buf);
	sprintf(buf, "%lf", genfd_ly);
	panel_setstr_value(genfd_ly_item, buf);
	sprintf(buf, "%lf", genfd_rot);
	panel_setstr_value(genfd_rot_item, buf);
    }
}

static void set_genfd_type(Widget w, int cd)
{
    if (cd) {
	XtSetSensitive(genfd_dodxdy_item, True);
	XtSetSensitive(genfd_donxny_item, False);
    } else {
	XtSetSensitive(genfd_donxny_item, True);
	XtSetSensitive(genfd_dodxdy_item, False);
    }
}

static void do_genfd_accept_proc(Widget w, int cd)
{
    char buf[256];

    genfd_doit = cd;

    genfdflag = GetChoice(genfd_type_item);
    strcpy(buf, panel_getstr_value(genfd_dx_item));
    sscanf(buf, "%lf", &genfd_dx);
    strcpy(buf, panel_getstr_value(genfd_dy_item));
    sscanf(buf, "%lf", &genfd_dy);
    strcpy(buf, panel_getstr_value(genfd_nx_item));
    sscanf(buf, "%d", &genfd_nx);
    strcpy(buf, panel_getstr_value(genfd_ny_item));
    sscanf(buf, "%d", &genfd_ny);
    strcpy(buf, panel_getstr_value(genfd_sx_item));
    sscanf(buf, "%lf", &genfd_sx);
    strcpy(buf, panel_getstr_value(genfd_sy_item));
    sscanf(buf, "%lf", &genfd_sy);
    strcpy(buf, panel_getstr_value(genfd_lx_item));
    sscanf(buf, "%lf", &genfd_lx);
    strcpy(buf, panel_getstr_value(genfd_ly_item));
    sscanf(buf, "%lf", &genfd_ly);
    strcpy(buf, panel_getstr_value(genfd_rot_item));
    sscanf(buf, "%lf", &genfd_rot);
    create_genfd_grid();
}

static void do_genfd_proc(void)
{
    char buf[256];

    genfdflag = GetChoice(genfd_type_item);
    strcpy(buf, panel_getstr_value(genfd_dx_item));
    sscanf(buf, "%lf", &genfd_dx);
    strcpy(buf, panel_getstr_value(genfd_dy_item));
    sscanf(buf, "%lf", &genfd_dy);
    strcpy(buf, panel_getstr_value(genfd_nx_item));
    sscanf(buf, "%d", &genfd_nx);
    strcpy(buf, panel_getstr_value(genfd_ny_item));
    sscanf(buf, "%d", &genfd_ny);
    strcpy(buf, panel_getstr_value(genfd_sx_item));
    sscanf(buf, "%lf", &genfd_sx);
    strcpy(buf, panel_getstr_value(genfd_sy_item));
    sscanf(buf, "%lf", &genfd_sy);
    strcpy(buf, panel_getstr_value(genfd_lx_item));
    sscanf(buf, "%lf", &genfd_lx);
    strcpy(buf, panel_getstr_value(genfd_ly_item));
    sscanf(buf, "%lf", &genfd_ly);
    strcpy(buf, panel_getstr_value(genfd_rot_item));
    sscanf(buf, "%lf", &genfd_rot);
    set_action(0);
    set_action(GENFD_1ST);
}

static void do_genfdflag_proc(int value)
{
    extern int genfdflag;
    genfdflag = value;
}

void create_genfd_popup(void)
{
    Widget bt, fr, rc;
    int i;
    genfdflag = 0;
    if (!genfd_frame) {
	genfd_frame = XmCreateDialogShell(app_shell, "Finite difference grid", NULL, 0);
	genfd_panel = XmCreateRowColumn(genfd_frame, "rc", NULL, 0);
	genfd_type_item = (Widget *) CreatePanelChoice1(genfd_panel,
						      "Generation mode:",
							3,
							"NX, NY",
							"DX, DY",
							0,
							NULL);
	for (i = 0; i < 2; i++) {
	    XtAddCallback(genfd_type_item[i + 2], XmNactivateCallback,
			  (XtCallbackProc) set_genfd_type, (XtPointer) ((long) i));
	}

	genfd_dodxdy_item = XmCreateFrame(genfd_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(genfd_dodxdy_item, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	genfd_dx_item = CreateTextItem2(rc, 12, "Spacing X:");
	genfd_dy_item = CreateTextItem2(rc, 12, "Y:");
	XtManageChild(rc);
	XtManageChild(genfd_dodxdy_item);

	genfd_donxny_item = XmCreateFrame(genfd_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(genfd_donxny_item, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	genfd_nx_item = CreateTextItem2(rc, 12, "# of points X:");
	genfd_ny_item = CreateTextItem2(rc, 12, "Y:");
	XtManageChild(rc);
	XtManageChild(genfd_donxny_item);

	fr = XmCreateFrame(genfd_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	genfd_lx_item = CreateTextItem2(rc, 12, "Length X:");
	genfd_ly_item = CreateTextItem2(rc, 12, "Y:");
	XtManageChild(rc);
	XtManageChild(fr);

	fr = XmCreateFrame(genfd_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	genfd_sx_item = CreateTextItem2(rc, 12, "Anchor point X:");
	genfd_sy_item = CreateTextItem2(rc, 12, "Y:");
	XtManageChild(rc);
	XtManageChild(fr);

	fr = XmCreateFrame(genfd_panel, "fr", NULL, 0);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	genfd_rot_item = CreateTextItem2(rc, 12, "Rotation (degrees):");
	XtManageChild(rc);
	XtManageChild(fr);

	rc = XmCreateRowColumn(genfd_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL,
		      XmNpacking, XmPACK_TIGHT, NULL);
	bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_genfd_accept_proc, (XtPointer) 1);
	bt = XtVaCreateManagedWidget("Pick rectangle", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_genfd_proc, 0);
	bt = XtVaCreateManagedWidget("*Specify rectangle...", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_genfd_proc, 0);
	bt = XtVaCreateManagedWidget("Preview", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_genfd_accept_proc, 0);
	bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				     NULL);
	XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, genfd_frame);
	XtManageChild(rc);
	XtManageChild(genfd_panel);
	set_genfd_type(NULL, 0);
    }
    XtRaise(genfd_frame);
    update_genfd();
}

void genfd_grid(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double dx1, dy1, dx2, dy2, m1, m2, b1, b2;
    double ang1, ang2;
    double x[4], y[4];
    double xc[4], yc[4], xt, yt, xb, yb, d;
    int i, j, ind;
    double get_depth_element();
    genfd_sx = x1;
    genfd_sy = y1;
/* translate to the origin */
    dx1 = (x2 - x1);
    dy1 = (y2 - y1);
    dx2 = (x3 - x1);
    dy2 = (y3 - y1);
/* compute angle to rotate */
    genfd_rot = ang1 = atan2(dy1, dx1);
    genfd_rot *= 180.0 / M_PI;
/* rotate so that (x1, y1) - (x2, y2) lies on the X-axis */
    x[0] = 0.0;
    y[0] = 0.0;
    x[2] = dx2 * cos(-ang1) - dy2 * sin(-ang1);
    y[2] = dx2 * sin(-ang1) + dy2 * cos(-ang1);
    x[1] = hypot(x[2], y[2]) * cos(atan(y[2] / x[2]));
    y[1] = 0.0;
    x[3] = 0.0;
    y[3] = hypot(x[2], y[2]) * sin(atan(y[2] / x[2]));

/*
   setcolor(3);
   my_move2(x[0], y[0]);
   for (i = 1; i < 4; i++) {
   my_draw2(x[i], y[i]);
   }
 */

/* rotate back to original position */
    for (i = 0; i < 4; i++) {
	xc[i] = x[i] * cos(ang1) - y[i] * sin(ang1) + x1;
	yc[i] = x[i] * sin(ang1) + y[i] * cos(ang1) + y1;
    }
    genfd_x1 = xc[0];
    genfd_x2 = xc[1];
    genfd_x3 = xc[2];
    genfd_x4 = xc[3];
    genfd_y1 = yc[0];
    genfd_y2 = yc[1];
    genfd_y3 = yc[2];
    genfd_y4 = yc[3];

    setcolor(2);
    my_move2(xc[0], yc[0]);
    for (i = 0; i < 4; i++) {
	my_draw2(xc[(i + 1) % 4], yc[(i + 1) % 4]);
    }
    setcolor(1);
    genfd_signx = x[2] < 0.0 ? -1 : 1;
    genfd_signy = y[3] < 0.0 ? -1 : 1;
    genfd_lx = fabs(x[2]);
    genfd_ly = fabs(y[3]);
    if (genfdflag) {
	genfd_nx = (int) (fabs(x[2] / genfd_dx) + 1);
	genfd_ny = (int) (fabs(y[3] / genfd_dy) + 1);
	dx1 = genfd_dx * genfd_signx;
	dy1 = genfd_dy * genfd_signy;
    } else {
	dx1 = x[2] / (genfd_nx - 1);
	dy1 = y[3] / (genfd_ny - 1);
	genfd_dx = fabs(dx1);
	genfd_dy = fabs(dy1);
    }
    update_genfd();
/*
   if (genfd_doit) {
   init_qadd_build(curbuild);
   }
   for (i = 0; i < genfd_nx; i++) {
   xb = i * dx1;
   for (j = 0; j < genfd_ny; j++) {
   yb = j * dy1;
   xt = xb * cos(ang1) - yb * sin(ang1) + x1;
   yt = xb * sin(ang1) + yb * cos(ang1) + y1;
   find_element(backgrid, xt, yt, &ind);
   if (ind == -1) {
   printf("Element not found in background grid\n");
   } else {
   d = get_depth_element(backgrid, ind, xt, yt);
   if (genfd_doit) {
   qadd_build(curbuild, xt, yt, d);
   }
   }
   }
   }
   if (genfd_doit) {
   flush_qadd_build(curbuild);
   }
   genfd_doit = 1;
 */
}

void create_genfd_grid(void)
{
    double dx1, dy1, ang1;
    double xc[4], yc[4], xt, yt, xb, yb, d;
    int i, j, k, ind, errcnt, cnt = 0;
    char buf[256];
    double get_depth_element();
    double cs, ss;

    if (genfd_dx == 0.0 || genfd_dy == 0.0) {
	errwin("DX or DY == 0.0");
	return;
    }
    if (genfd_nx == 0 || genfd_ny == 0) {
	errwin("NX or NY == 0");
	return;
    }
    sprintf(buf, "Generating grid with %d points, OK?", genfd_nx * genfd_ny);
    if (!yesno(buf, "Press Yes or No", "Yes", "No")) {
	return;
    }
    setcolor(1);

    prep_find(backgrid);
    if (genfd_doit) {
	init_qadd_build(curbuild);
    }
    dx1 = genfd_dx * genfd_signx;
    dy1 = genfd_dy * genfd_signy;
    ang1 = genfd_rot * M_PI / 180.0;
    cs = cos(ang1);
    ss = sin(ang1);
    errcnt = 0;
    for (i = 0; i < genfd_nx; i++) {
	xb = i * dx1;
	for (j = 0; j < genfd_ny; j++) {
	    yb = j * dy1;
	    xt = xb * cs - yb * ss + genfd_sx;
	    yt = xb * ss + yb * cs + genfd_sy;
	    ind = find_element3(backgrid, xt, yt);
	    if (ind <= -1 || ind >= grid[backgrid].nmel) {
		errcnt++;
		qadd_build(curbuild, xt, yt, 0.0);
	    } else {
		d = get_depth_element(backgrid, ind, xt, yt);
		if (genfd_doit) {
		    qadd_build(curbuild, xt, yt, d);
		}
		box(xt, yt);
	    }
	    cnt++;
	}
    }
    flush_pending();
    if (genfd_doit) {
	flush_qadd_build(curbuild);
    }
    end_find();
    if (errcnt) {
	sprintf(buf, "Warning, %d points were not located in the background grid", buf);
	errwin(buf);
    }
    genfd_doit = 1;
}

typedef struct {
    Widget top;
    Widget *sym_item;
    Widget symsize_item;
    Widget *symcolor_item;
    Widget symsamp_item;
    Widget isol_item;
} BuildPropsUI;

void accept_buildprops_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    BuildPropsUI *ui = (BuildPropsUI *) clientd;
    build[curbuild].color = GetChoice(ui->symcolor_item);
    build[curbuild].sym = GetChoice(ui->sym_item);
    build[curbuild].isol = XmToggleButtonGetState(ui->isol_item);
}

void update_buildprops(BuildPropsUI * ui)
{
    SetChoice(ui->symcolor_item, build[curbuild].color);
    SetChoice(ui->sym_item, build[curbuild].sym);
    XmToggleButtonSetState(ui->isol_item, build[curbuild].isol, False);
}

void create_buildprops_frame(void)
{
    static BuildPropsUI ui;
    Widget panel, bt, rc, sep;
    if (ui.top) {
	XtRaise(ui.top);
	update_buildprops(&ui);
	return;
    }
    ui.top = XmCreateDialogShell(app_shell, "Display build points properties", NULL, 0);
    handle_close(ui.top);
    panel = XmCreateRowColumn(ui.top, "panel", NULL, 0);
    ui.sym_item = CreatePanelChoice1(panel, "Symbol",
				     6,
				     "None",
				     "Dot",
				     "Square",
				     "Circle",
				     "Number",
				     NULL, 0);
    ui.symsize_item = XtVaCreateManagedWidget("Symbol size", xmScaleWidgetClass, panel,
					      XmNwidth, 150,
					      XmNminimum, 0,
					      XmNmaximum, 100,
			    XmNvalue, (int) (build[curbuild].size * 100),
					      XmNshowValue, True,
				  XmNprocessingDirection, XmMAX_ON_RIGHT,
	      XmNtitleString, XmStringCreateLtoR("Symbol size", charset),
					    XmNorientation, XmHORIZONTAL,
					      NULL);
    ui.symcolor_item = CreateColorChoice(panel, "Color", 1);
    ui.isol_item = XmCreateToggleButton(panel, "Display using isolines", NULL, 0);
    XtManageChild(ui.isol_item);

    sep = XmCreateSeparatorGadget(panel, "sep", NULL, 0);
    XtManageChild(sep);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) accept_buildprops_proc,
		  (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);

    XtManageChild(rc);

    update_buildprops(&ui);
    XtManageChild(panel);
    XtManageChild(ui.top);
}
