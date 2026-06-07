/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 *
 */

/*
 * adcircwin.c - read flows and set display parameters for
 * ADCIRC flows (2d time dependent data)
 *
 */

#ifndef lint
static char RCSid[] = "$Id: adcircwin.c,v 1.8 2005/06/17 15:23:26 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "elio.h"

extern Widget app_shell;
extern XmStringCharSet charset;

void create_isolines_popup(char *title, Isolparms * cd, int cm);
void create_adcircinundate_frame(void);

double get_current_time(void);

/* contents */
static void set_curadcirc(Widget w, int cd);
static void set_curadcircelm(Widget w, int cd);
static void query_adcirc(void);
void create_adcirc_frame(Widget w, XtPointer clientd, XtPointer calld);
void accept_adcirc(Widget w, XtPointer clientd, XtPointer calld);
void create_adcircw_frame(void);
static void do_adcirc_isolines(Widget w, int cd);

/* 3d read routine */
void create_elcirc_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_scalar_frame(void);
void update_scalar_frame(Widget w, XtPointer clientd, XtPointer calld);
static void accept_scalar_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_vector_frame(void);
void update_vector_frame(Widget w, XtPointer clientd, XtPointer calld);
static void accept_vector_frame(Widget w, XtPointer clientd, XtPointer calld);

void create_sadcirc_frame(void);
void update_adcirc_flow(void);
static void sadcirc_accept_proc(void);
static void sadcirc_view_proc(void);
static void sadcirc_viewmag_proc(void);
static void sadcirc_viewmax_proc(void);
static void adcircel_place_notify_proc(void);
static void adcircel_accept_notify_proc(void);
static void adcircel_edit_notify_proc(void);

static Widget adcirc_frame;

static Widget sadcirc_frame;
static Widget sadcirc_panel;

/*
 * Panel item declarations
 */
static Widget adcirc_files_list_item;
static Widget adcirc_swap_item;
static Widget adcircel_emin_text_item;
static Widget adcircel_emax_text_item;
static Widget *adcirc_filetype_item;
static Widget *adcirc_flow_item;
static Widget *flow_toggle_item;
static Widget *flow_color_item;
static Widget elev_toggle_item;
static Widget wind_toggle_item;
static Widget mag_toggle_item;
static Widget mag_region_item;
static Widget mag_depth_item;
static Widget *sadcirc_elevm_item;
static Widget elevm_toggle_item;
static Widget elevm_node_item;
static Widget elevm_attach_item;
static Widget flowelev_toggle_item;
static Widget flowelev_depth_item;
static Widget maxelev_toggle_item;
static Widget maxelev_val_item;
static Widget maxelev_depth_item;
static Widget append_item;

/*
 * set the current ADCIRC flow
 */
static void set_curadcirc(Widget w, int cd)
{
    if (cd != MAXADCIRC) {
	curadcirc = cd;
	update_adcirc_flow();
    }
}

/*
 * set the current ADCIRC elevation marker
 */
static void set_curadcircelm(Widget w, int cd)
{
    curadcircem = cd;
    update_adcirc_flow();
}

/*
 * find out about data at a node at this time step
 */
static void query_adcirc(void)
{
    extern Widget locate_point_message;
    set_action(0);
    set_action(QUERY_ADCIRC);
    create_points_frame();
    SetLabel(locate_point_message, "ADCIRC:");
}

/*
 * find out about maximum elevations at a node
 */
static void query_adcirc_maxelev(void)
{
    extern Widget locate_point_message;
    set_action(0);
    set_action(QUERY_ADCIRC_MAXELEV);
    create_points_frame();
    SetLabel(locate_point_message, "SELFE max scalars:");
}

/*
 * set isolines of elevations and isolines of maximum elevations
 */
static void do_adcirc_isolines(Widget w, int cd)
{
    switch (cd) {
    case 0:
	if (g[cg].flowt[curadcirc].display_elevdepth == ON) {
	    g[cg].flowt[curadcirc].elevip.cmin = flowt[curadcirc].emin + flowt[curadcirc].g.dmin;
	    g[cg].flowt[curadcirc].elevip.cmax = flowt[curadcirc].emax + flowt[curadcirc].g.dmax;
	} else {
	    g[cg].flowt[curadcirc].elevip.cmin = flowt[curadcirc].emin;
	    g[cg].flowt[curadcirc].elevip.cmax = flowt[curadcirc].emax;
	}
	create_isolines_popup("Scalar isolines", &g[cg].flowt[curadcirc].elevip, 1);
	break;
    case 1:
	g[cg].flowt[curadcirc].maxelevip.cmin = 0.0;
	g[cg].flowt[curadcirc].maxelevip.cmax = flowt[curadcirc].emax + g[cg].flowt[curadcirc].d;
	create_isolines_popup("Scalar max isolines", &g[cg].flowt[curadcirc].maxelevip, 1);
	break;
    case 2:
	g[cg].flowt[curadcirc].magip.cmin = 0.0;
	g[cg].flowt[curadcirc].magip.cmax = 2.0;
	create_isolines_popup("Magnitudes of velocity", &g[cg].flowt[curadcirc].magip, 1);
	break;
    }
}

/*
 * create the popup used to set properties for the display of ADCIRC flows
 */
void create_sadcirc_frame(void)
{
    Widget wbut, fr, sep, rc, rc2;
    int i, got_one = 0;;
    setistop();
    for (i = 0; i < MAXADCIRC; i++) {
	if (flowt[i].active == ON) {
	    got_one = 1;
	}
    }
    if (!sadcirc_frame) {
	sadcirc_frame = XmCreateDialogShell(app_shell, "Vector/Scalar setup", NULL, 0);
	handle_close(sadcirc_frame);
	sadcirc_panel = XmCreateRowColumn(sadcirc_frame, "sadcirc_rc", NULL, 0);
	adcirc_flow_item = CreatePanelChoice1(sadcirc_panel, "Apply to data set: ", 7, "1", "2", "3", "4", "5", "All", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(adcirc_flow_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadcirc, (XtPointer) i);
	}
/* elevation markers */
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sadcirc_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	sadcirc_elevm_item = CreatePanelChoice2(rc, "Scalar marker: ", 4, 21, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 0, 0);
	for (i = 0; i < MAXELEVMARKERS; i++) {
	    XtAddCallback(sadcirc_elevm_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadcircelm, (XtPointer) i);
	}
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	elevm_attach_item = XtVaCreateManagedWidget("*Attach to node", xmToggleButtonWidgetClass, rc2, NULL);
	elevm_node_item = CreateTextItem2(rc2, 5, "*Node: ");
	XtManageChild(rc2);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	adcircel_emin_text_item = CreateTextItem2(rc2, 7, "Marker min: ");
	adcircel_emax_text_item = CreateTextItem2(rc2, 7, "max: ");
	XtManageChild(rc2);
	elevm_toggle_item = XtVaCreateManagedWidget("Display this elevation marker", xmToggleButtonWidgetClass, rc, NULL);
	elev_toggle_item = XtVaCreateManagedWidget("Display all elevation markers", xmToggleButtonWidgetClass, rc, NULL);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_place_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Edit", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_edit_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_accept_notify_proc, NULL);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

/* cut here for velocities */
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sadcirc_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	flow_toggle_item = CreatePanelChoice1(rc, "Velocity vectors: ", 5, "Not displayed", "Display at nodes", "Display at centers", "Display ellipse", 0, 0);

	flow_color_item = CreateColorChoice(rc, "Display vectors using color:", 1);
	mag_toggle_item = XtVaCreateManagedWidget("Display magnitudes of velocity", xmToggleButtonWidgetClass, rc, NULL);
	wind_toggle_item = XtVaCreateManagedWidget("Set wind flag", xmToggleButtonWidgetClass, rc, NULL);

	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of magnitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, (XtPointer) 2);
	wbut = XtVaCreateManagedWidget("View magnitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_viewmag_proc, 0);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc, NULL);
	XtManageChild(rc2);

	XtManageChild(fr);
	XtManageChild(rc);
/* end of velocities */

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sadcirc_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	flowelev_toggle_item = XtVaCreateManagedWidget("Display isolines of scalars", xmToggleButtonWidgetClass, rc, NULL);
	flowelev_depth_item = XtVaCreateManagedWidget("Add depth", xmToggleButtonWidgetClass, rc, NULL);
	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of scalars...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, 0);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_view_proc, 0);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sadcirc_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	maxelev_toggle_item = XtVaCreateManagedWidget("Display isolines of maximum scalar values", xmToggleButtonWidgetClass, rc, NULL);
	maxelev_val_item = XtVaCreateManagedWidget("Display max scalars at nodes", xmToggleButtonWidgetClass, rc, NULL);
	maxelev_depth_item = CreateTextItem2(rc, 7, "Add depth: ");
	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, (XtPointer) 1);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc_maxelev, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_viewmax_proc, 0);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

	rc = XmCreateRowColumn(sadcirc_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_elcirc_frame, NULL);
	wbut = XtVaCreateManagedWidget("Inundation...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_adcircinundate_frame, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, sadcirc_frame);
	XtManageChild(rc);
	XtManageChild(sadcirc_panel);
	xv_setstr(maxelev_depth_item, "0.0");
    }
    update_adcirc_flow();
    XtRaise(sadcirc_frame);
}

static void sadcirc_view_proc(void)
{
    int i;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEadcircelevXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"SELFE scalars at nodes\"\n");
	fprintf(fp, "@subtitle \"Step %d, Time = %lf\"\n", c + 1, t);
	fprintf(fp, "@xaxis label \"Node number\"\n");
	fprintf(fp, "@yaxis label \"Elevation\"\n");
	for (i = 0; i < flowt[curadcirc].g.nmnp; i++) {
	    fprintf(fp, "%d %lf\n", i + 1, flowt[curadcirc].f[c].e[i]);
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    } else {
	errwin("Unable to open file");
    }
}

static void sadcirc_viewmag_proc(void)
{
    int i;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEadcircmagXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"SELFE magnitudes of velocity at nodes\"\n");
	fprintf(fp, "@subtitle \"Step %d, Time = %lf\"\n", c + 1, t);
	fprintf(fp, "@xaxis label \"Node number\"\n");
	fprintf(fp, "@yaxis label \"Magnitude\"\n");
	for (i = 0; i < flowt[curadcirc].g.nmnp; i++) {
	    fprintf(fp, "%d %lf\n", i + 1, hypot(flowt[curadcirc].f[c].u[i], flowt[curadcirc].f[c].v[i]));
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    } else {
	errwin("Unable to open file");
    }
}

static void sadcirc_viewmax_proc(void)
{
    int i;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    double d = atof(xv_getstr(maxelev_depth_item));
    strcpy(tbuf, "/tmp/ACEadcircelevmaxXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");

    if (fp != NULL) {
	fprintf(fp, "@title \"SELFE maximum scalars at nodes\"\n");
	fprintf(fp, "@xaxis label \"Depth\"\n");
	fprintf(fp, "@yaxis label \"Elevation\"\n");
	if (compute_adcirc_maxelev(curadcirc, 1)) {
	    for (i = 0; i < flowt[curadcirc].g.nmnp; i++) {
		if (nregion) {
		    double x, y;
		    x = flowt[curadcirc].g.xord[i];
		    y = flowt[curadcirc].g.yord[i];
		    if (inregion2(regionx, regiony, nregion, x, y)) {
			fprintf(fp, "%lf %lf %d\n", flowt[curadcirc].g.depth[i] + d, flowt[curadcirc].global_emax[i], i + 1);
		    }
		} else {
/*
		    fprintf(fp, "%d %lf\n", i + 1, flowt[curadcirc].global_emax[i]);
*/
		    fprintf(fp, "%lf %lf %d\n", flowt[curadcirc].g.depth[i] + d, flowt[curadcirc].global_emax[i], i + 1);
		}
	    }
	    fclose(fp);
	    sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	    system(combuf);
	}
    } else {
	errwin("Unable to open file");
    }
}

/*
 * freshen up the ADCIRC setup panel
 */
void update_adcirc_flow(void)
{
    int flowno;
    char buf[256];
    if (sadcirc_frame) {
	SetChoice(adcirc_flow_item, curadcirc);
	flowno = curadcirc;
	SetChoice(flow_color_item, g[cg].flowt[flowno].p.color);
	switch (g[cg].flowt[flowno].display) {
	case OFF:
	    SetChoice(flow_toggle_item, 0);
	    break;
	case NODES:
	    SetChoice(flow_toggle_item, 1);
	    break;
	case CENTER:
	    SetChoice(flow_toggle_item, 2);
	    break;
	case ELLIPSE:
	    SetChoice(flow_toggle_item, 3);
	    break;
	}
	XmToggleButtonSetState(elev_toggle_item, g[cg].flowt[flowno].display_elevmarkers == ON, False);
	XmToggleButtonSetState(flowelev_toggle_item, g[cg].flowt[flowno].display_elev == ON, False);
	XmToggleButtonSetState(flowelev_depth_item, g[cg].flowt[flowno].display_elevdepth == ON, False);
	XmToggleButtonSetState(maxelev_toggle_item, g[cg].flowt[flowno].display_maxelev == ON, False);
	XmToggleButtonSetState(maxelev_val_item, g[cg].flowt[flowno].display_maxelevval == ON, False);
	XmToggleButtonSetState(elevm_toggle_item, g[cg].flowt[flowno].em[curadcircem].display == ON, False);
	XmToggleButtonSetState(mag_toggle_item, g[cg].flowt[flowno].display_mag == ON, False);
	XmToggleButtonSetState(wind_toggle_item, g[cg].flowt[flowno].display_wind == ON, False);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].d);
	xv_setstr(maxelev_depth_item, buf);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].em[curadcircem].emin);
	xv_setstr(adcircel_emin_text_item, buf);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].em[curadcircem].emax);
	xv_setstr(adcircel_emax_text_item, buf);
	SetChoice(flow_color_item, g[cg].flowt[flowno].p.color);
	SetChoice(sadcirc_elevm_item, curadcircem);
    }
}

/*
 * Accept the current state of the ADCIRC popup
 */
static void sadcirc_accept_proc(void)
{
    int i, flowno, type, start, stop;
    char buf[256];
    flowno = GetChoice(adcirc_flow_item);
    if (flowno == MAXADCIRC) {
	start = 0;
	stop = MAXADCIRC - 1;
    } else {
	curadcirc = start = stop = flowno;
	if (flowt[curadcirc].active == OFF) {
	    sprintf(buf, "Warning, no data for ADCIRC flow %d", curadcirc + 1);
	    errwin(buf);
	}
    }
    type = GetChoice(flow_toggle_item);
    for (i = start; i <= stop; i++) {
	switch (type) {
	case 0:
	    g[cg].flowt[i].display = OFF;
	    break;
	case 1:
	    g[cg].flowt[i].display = NODES;
	    break;
	case 2:
	    g[cg].flowt[i].display = CENTER;
	    break;
	case 3:
	    g[cg].flowt[i].display = ELLIPSE;
	    break;
	}
	g[cg].flowt[i].display_elev = XmToggleButtonGetState(flowelev_toggle_item) ? ON : OFF;
	g[cg].flowt[flowno].display_elevdepth = XmToggleButtonGetState(flowelev_depth_item) ? ON : OFF;
	g[cg].flowt[i].display_maxelev = XmToggleButtonGetState(maxelev_toggle_item) ? ON : OFF;
	g[cg].flowt[i].display_maxelevval = XmToggleButtonGetState(maxelev_val_item) ? ON : OFF;
	g[cg].flowt[i].d = atof(XmTextGetString(maxelev_depth_item));
	g[cg].flowt[i].display_elevmarkers = XmToggleButtonGetState(elev_toggle_item) ? ON : OFF;
	g[cg].flowt[i].p.color = GetChoice(flow_color_item);
	g[cg].flowt[i].em[curadcircem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
	g[cg].flowt[i].em[curadcircem].emin = atof(XmTextGetString(adcircel_emin_text_item));
	g[cg].flowt[i].em[curadcircem].emax = atof(XmTextGetString(adcircel_emax_text_item));
	g[cg].flowt[i].display_mag = XmToggleButtonGetState(mag_toggle_item) ? ON : OFF;
	g[cg].flowt[i].display_wind = XmToggleButtonGetState(wind_toggle_item) ? ON : OFF;
    }
    update_display();
    if (flowt[curadcirc].active == ON) {
	display_image();
    }
}

/*
 * place the currently selected elevation marker
 */
static void adcircel_place_notify_proc(void)
{
    g[cg].flowt[curadcirc].em[curadcircem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
    g[cg].flowt[curadcirc].em[curadcircem].emin = atof(XmTextGetString(adcircel_emin_text_item));
    g[cg].flowt[curadcirc].em[curadcircem].emax = atof(XmTextGetString(adcircel_emax_text_item));
    set_adcircelev();
}

/*
 * edit an elevation marker
 */
static void adcircel_edit_notify_proc(void)
{
    set_action(0);
    set_action(EDIT_ADCIRC_ELEV);
}

static void adcircel_accept_notify_proc(void)
{
    g[cg].flowt[curadcirc].em[curadcircem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
    g[cg].flowt[curadcirc].em[curadcircem].emin = atof(XmTextGetString(adcircel_emin_text_item));
    g[cg].flowt[curadcirc].em[curadcircem].emax = atof(XmTextGetString(adcircel_emax_text_item));
}

static Widget inundate_frame;
static Widget inundate_toggle_item;
static Widget restrict_toggle_item;
static Widget *wet_color_item;
static Widget *dry_color_item;
static Widget *wetdry_color_item;

/*
    int display_inun;
    int display_irestrict;
    Props wet, dry, wetdry;
*/

void inundate_adcircaccept_proc()
{
    g[cg].flowt[curadcirc].wet.color = GetChoice(wet_color_item);
    g[cg].flowt[curadcirc].dry.color = GetChoice(dry_color_item);
    g[cg].flowt[curadcirc].wetdry.color = GetChoice(wetdry_color_item);
    g[cg].flowt[curadcirc].display_inun = XmToggleButtonGetState(inundate_toggle_item) ? ON : OFF;
    g[cg].flowt[curadcirc].display_irestrict = XmToggleButtonGetState(restrict_toggle_item) ? ON : OFF;
}

void update_adcircinundate_flow()
{
    if (inundate_frame) {
	SetChoice(wet_color_item, g[cg].flowt[curadcirc].wet.color);
	SetChoice(dry_color_item, g[cg].flowt[curadcirc].dry.color);
	SetChoice(wetdry_color_item, g[cg].flowt[curadcirc].wetdry.color);
	XmToggleButtonSetState(inundate_toggle_item, g[cg].flowt[curadcirc].display_inun == ON, False);
	XmToggleButtonSetState(restrict_toggle_item, g[cg].flowt[curadcirc].display_irestrict == ON, False);
    }
}

void create_adcircinundate_frame(void)
{
    Widget inundate_panel, lab, wbut, fr, bb, sep, rc, rc2;
    int i, got_one = 0;
    Arg al[10];
    setistop();
    if (!inundate_frame) {
	inundate_frame = XmCreateDialogShell(app_shell, "Inundation", NULL, 0);
	handle_close(inundate_frame);
	inundate_panel = XmCreateRowColumn(inundate_frame, "inundate_rc", NULL, 0);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, inundate_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	inundate_toggle_item = XtVaCreateManagedWidget("Display inundatations", xmToggleButtonWidgetClass, rc, NULL);
	restrict_toggle_item = XtVaCreateManagedWidget("Restrict drawing to wet elements", xmToggleButtonWidgetClass, rc, NULL);

	wet_color_item = CreateColorChoice(rc, "Display wet areas using color:", 1);
	dry_color_item = CreateColorChoice(rc, "Display dry areas using color:", 1);
	wetdry_color_item = CreateColorChoice(rc, "Display partially dry areas using color:", 1);

	XtManageChild(rc);
	XtManageChild(fr);

	rc = XmCreateRowColumn(inundate_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) inundate_adcircaccept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, inundate_frame);
	XtManageChild(rc);
	XtManageChild(inundate_panel);
    }
    update_adcircinundate_flow();
    XtRaise(inundate_frame);
}

/*
 * create the file selection dialog for SELFE 3d data
 */
typedef struct {
    Widget top;
    Widget file;
    Widget *flow;
    Widget *type;
    Widget typeflag;
    Widget level;
    Widget depth;
    Widget nlevels;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget missing;
    Widget missingval;
    Widget append;
    Browser b;
    Browser bz;
} elcircUI;

static elcircUI eui;

void accept_elcirc(Widget w, XtPointer clientd, XtPointer calld);

void create_elcirc_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i;
    Widget wbut, panel, rc, rc2, wtmp;
    setistop();
    if (!eui.top) {
	eui.top = XmCreateDialogShell(app_shell, "SELFE", NULL, 0);
	panel = XmCreateRowColumn(eui.top, "panel", NULL, 0);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	eui.b.file = CreateTextItem2(rc, 20, "Data file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & eui.b);
	XtManageChild(rc);

	eui.flow = CreatePanelChoice1(panel, "Read to data set: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(eui.flow[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadcirc, (XtPointer) i);
	}
	eui.typeflag = XtVaCreateManagedWidget("Read depth:", xmToggleButtonWidgetClass, panel, NULL);
        XtVaCreateManagedWidget("Need _zcor.63 for Sigma.", xmLabelWidgetClass, panel, NULL);
        rc = XmCreateRowColumn(panel, "rc", NULL, 0);
        XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
        eui.bz.file = CreateTextItem2(rc, 20, "Elev/Zcor file:");
        wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
        XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & eui.bz);
        XtManageChild(rc);

	eui.level = CreateTextItem2(panel, 7, "Read level (from 1) or depth: ");
	eui.append = XtVaCreateManagedWidget("Append to active data set", xmToggleButtonWidgetClass, panel, NULL);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	eui.sample = XtVaCreateManagedWidget("Sample steps:", xmToggleButtonWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	eui.start = CreateTextItem2(rc, 4, "Start: ");
	xv_setstr(eui.start, "1");
	eui.stop = CreateTextItem2(rc, 4, "Stop: ");
	xv_setstr(eui.stop, "1");
	eui.skip = CreateTextItem2(rc, 4, "Skip: ");
	xv_setstr(eui.skip, "1");
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	eui.missing = XtVaCreateManagedWidget("Missing data:", xmToggleButtonWidgetClass, rc, NULL);
	eui.missingval = CreateTextItem2(rc, 4, "Value: ");
	xv_setstr(eui.missingval, "-99.0");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_elcirc, (XtPointer) & eui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) eui.top);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    XtRaise(eui.top);
}

void accept_elcirc(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256], file[2048], zname[2048];
    elcircUI *q = (elcircUI *) clientd;
    int flowno, level, nlevels, filet = 0;
    int sample, start = -1, stop = 1, skip = 1;
    int missing;
    int append;
    int err;
    int typeflag;
    double depth;
    double missingval;
    ElcircHeader h;


    strcpy(file, (char *) xv_getstr(q->b.file));
    strcpy(zname, (char *) xv_getstr(q->bz.file));

    strcpy(buf, (char *) xv_getstr(q->level));
    typeflag = XmToggleButtonGetState(q->typeflag);
    if (typeflag) {
        depth = atof(buf);
    } else {
        level = atoi(buf) - 1;
    	if (level < 0) {
	    errwin("Level must start from 1, operation canceled.");
	    return;
        }
    }

    err = ElioGetHeader(file, &h);
    ElioFreeHeader(&h); /* don't need allocated members */
    if (err) {
        fprintf(stderr, "Error in ElioGetHeader(): Error # %d\n", err);
        return;
    }
    if (!typeflag && (level > h.nvrt)) {
	errwin("Level is greater that the total number of levels, operation canceled.");
	return;
    }

    append = XmToggleButtonGetState(q->append);
    sample = XmToggleButtonGetState(q->sample);
    if (sample) {
	strcpy(buf, (char *) xv_getstr(q->start));
	start = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->stop));
	stop = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->skip));
	skip = atoi(buf);
    } else {
	start = -1;
	skip = 1;
    }
    missing = XmToggleButtonGetState(q->missing);
    if (missing) {
	strcpy(buf, (char *) xv_getstr(q->missingval));
	missingval = atof(buf);
    } else {
	missing = 0;
    }
    flowno = GetChoice(q->flow);
    if (!append && flowt[flowno].active == ON) {
	if (!yesno("Flow is active, kill it?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(q->top);
    err = 0;
    if (typeflag) {
        if (!ReadElcircDepth(flowno, file, zname, depth, start, stop, skip, missing, missingval, append)) {
            sprintf(buf, "Error reading file %s", file);
            errwin(buf);
            return;
	}
        err = 0;
     } else {
        if (!ReadElcirc(flowno, file, level, start, stop, skip, missing, missingval, append)) {
	    sprintf(buf, "Error reading file %s", file);
	    errwin(buf);
	    return;
        }
        err = 0;
    }
    if (!err) {
	set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
	load_clock(ADCIRC, flowno);
	XtUnmanageChild(q->top);
	create_sadcirc_frame();
/* TODO activate when vectors and scalrs are split into their own objects.
	if (flowt[flowno].type == 1) {
	    create_scalar_frame();
	} else if (flowt[flowno].type == 2) {
	    create_vector_frame();
	} else {
	    create_sadcirc_frame();
	}
*/
    }
    unset_wait_cursor(q->top);
}

typedef struct {
    Widget top;
    Widget *flow;
    Widget *em;
    Widget emin;
    Widget emax;
    Widget etoggle;
    Widget etoggleall;
    Widget isol;
    Widget depth;
    Widget maxisol;
    Widget maxval;
    Widget adddepth;
} scalarUI;

static scalarUI sui;

static void set_scalar(Widget w, int cd)
{
    if (cd != MAXADCIRC) {
	curadcirc = cd;
	update_scalar_frame(NULL, NULL, NULL);
    }
}

static void set_vector(Widget w, int cd)
{
    if (cd != MAXADCIRC) {
	curadcirc = cd;
	update_vector_frame(NULL, NULL, NULL);
    }
}

static void set_scalar_em(Widget w, int cd)
{
    curadcircem = cd;
    update_scalar_frame(NULL, NULL, NULL);
}

/*
 * Scalar data set up frame
 */
void create_scalar_frame(void)
{
    Widget wbut, panel, fr, sep, rc, rc2;
    int i, got_one = 0;
    setistop();
    for (i = 0; i < MAXADCIRC; i++) {
	if (flowt[i].active == ON) {
	    got_one = 1;
	}
    }
    if (!sui.top) {
	sui.top = XmCreateDialogShell(app_shell, "Scalars setup", NULL, 0);
	handle_close(sui.top);
	panel = XmCreateRowColumn(sui.top, "sadcirc_rc", NULL, 0);
	sui.flow = CreatePanelChoice1(panel, "Apply to flow: ", 7, "1", "2", "3", "4", "5", "All", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(sui.flow[i + 2], XmNactivateCallback, (XtCallbackProc) set_scalar, (XtPointer) i);
	}
/* elevation markers */
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	sui.em = CreatePanelChoice2(rc, "Scalar marker: ", 4, 21, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 0, 0);
	for (i = 0; i < MAXELEVMARKERS; i++) {
	    XtAddCallback(sui.em[i + 2], XmNactivateCallback, (XtCallbackProc) set_scalar_em, (XtPointer) i);
	}
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	sui.emin = CreateTextItem2(rc2, 7, "Marker min: ");
	sui.emax = CreateTextItem2(rc2, 7, "max: ");
	XtManageChild(rc2);
	sui.etoggle = XtVaCreateManagedWidget("Display this elevation marker", xmToggleButtonWidgetClass, rc, NULL);
	sui.etoggleall = XtVaCreateManagedWidget("Display all elevation markers", xmToggleButtonWidgetClass, rc, NULL);

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_place_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Edit", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_edit_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) adcircel_accept_notify_proc, NULL);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

/* isolines */
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	sui.isol = XtVaCreateManagedWidget("Display isolines of scalars", xmToggleButtonWidgetClass, rc, NULL);
	sui.depth = XtVaCreateManagedWidget("Add depth", xmToggleButtonWidgetClass, rc, NULL);
	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of scalars...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, 0);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_view_proc, 0);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

/* isolines of max */
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	sui.maxisol = XtVaCreateManagedWidget("Display isolines of maximum scalar values", xmToggleButtonWidgetClass, rc, NULL);
	sui.maxval = XtVaCreateManagedWidget("Display maximum scalar values at nodes", xmToggleButtonWidgetClass, rc, NULL);
	sui.adddepth = CreateTextItem2(rc, 7, "Add depth: ");
	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, (XtPointer) 1);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc_maxelev, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_viewmax_proc, 0);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_scalar_frame, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_elcirc_frame, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, sui.top);
	XtManageChild(rc);
	XtManageChild(panel);
	xv_setstr(sui.adddepth, "0.0");
    }
    update_scalar_frame(NULL, NULL, NULL);
    XtRaise(sui.top);
}

typedef struct {
    Widget top;
    Widget *flow;
    Widget *display;
    Widget *color;
    Widget mag;
    Widget wind;
} vectorUI;

static vectorUI vui;

/*
 * Vector data set up frame
 */
void create_vector_frame(void)
{
    Widget wbut, panel, fr, sep, rc, rc2;
    int i, got_one = 0;;
    setistop();
    for (i = 0; i < MAXADCIRC; i++) {
	if (flowt[i].active == ON) {
	    got_one = 1;
	}
    }
    if (!vui.top) {
	vui.top = XmCreateDialogShell(app_shell, "ADCIRC setup", NULL, 0);
	handle_close(vui.top);
	panel = XmCreateRowColumn(vui.top, "panel", NULL, 0);
	vui.flow = CreatePanelChoice1(panel, "Apply to flow: ", 7, "1", "2", "3", "4", "5", "All", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(vui.flow[i + 2], XmNactivateCallback, (XtCallbackProc) set_vector, (XtPointer) i);
	}
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	vui.display = CreatePanelChoice1(rc, "Velocity vectors: ", 5, "Not displayed", "Display at nodes", "Display at centers", "Display ellipse", 0, 0);

	vui.color = CreateColorChoice(rc, "Display vectors using color:", 1);
	vui.mag = XtVaCreateManagedWidget("Display magnitudes of velocity", xmToggleButtonWidgetClass, rc, NULL);
	vui.wind = XtVaCreateManagedWidget("Set wind flag", xmToggleButtonWidgetClass, rc, NULL);

	rc2 = XmCreateRowColumn(rc, "rc", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of magnitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_adcirc_isolines, (XtPointer) 2);
	wbut = XtVaCreateManagedWidget("View magnitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sadcirc_viewmag_proc, 0);
	wbut = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) query_adcirc, NULL);
	XtManageChild(rc2);

	XtManageChild(fr);
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_vector_frame, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_elcirc_frame, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, vui.top);
	XtManageChild(rc);
	XtManageChild(panel);
    }
    update_vector_frame(NULL, NULL, NULL);
    XtRaise(vui.top);
}

/*
 * freshen up the ADCIRC setup panel
 */
void update_scalar_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int flowno;
    char buf[256];
    if (sui.top) {
	SetChoice(sui.flow, curadcirc);
	flowno = curadcirc;
	XmToggleButtonSetState(sui.etoggleall, g[cg].flowt[flowno].display_elevmarkers == ON, False);
	XmToggleButtonSetState(sui.isol, g[cg].flowt[flowno].display_elev == ON, False);
	XmToggleButtonSetState(sui.depth, g[cg].flowt[flowno].display_elevdepth == ON, False);
	XmToggleButtonSetState(sui.maxisol, g[cg].flowt[flowno].display_maxelev == ON, False);
	XmToggleButtonSetState(sui.maxval, g[cg].flowt[flowno].display_maxelevval == ON, False);
	XmToggleButtonSetState(sui.etoggle, g[cg].flowt[flowno].em[curadcircem].display == ON, False);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].d);
	xv_setstr(sui.adddepth, buf);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].em[curadcircem].emin);
	xv_setstr(sui.emin, buf);
	sprintf(buf, "%.3lf", g[cg].flowt[flowno].em[curadcircem].emax);
	xv_setstr(sui.emax, buf);
	SetChoice(sui.em, curadcircem);
    }
}

static void accept_scalar_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i, flowno, type, start, stop;
    char buf[256];
    flowno = GetChoice(sui.flow);
    if (flowno == MAXADCIRC) {
	start = 0;
	stop = MAXADCIRC - 1;
    } else {
	curadcirc = start = stop = flowno;
	if (flowt[curadcirc].active == OFF) {
	    sprintf(buf, "Warning, no data for ADCIRC flow %d", curadcirc + 1);
	    errwin(buf);
	}
    }
    for (i = start; i <= stop; i++) {
	g[cg].flowt[i].display_elev = XmToggleButtonGetState(sui.isol) ? ON : OFF;
	g[cg].flowt[flowno].display_elevdepth = XmToggleButtonGetState(sui.depth) ? ON : OFF;
	g[cg].flowt[i].display_maxelev = XmToggleButtonGetState(sui.maxisol) ? ON : OFF;
	g[cg].flowt[i].display_maxelevval = XmToggleButtonGetState(sui.maxval) ? ON : OFF;
	g[cg].flowt[i].d = atof(XmTextGetString(sui.adddepth));
	g[cg].flowt[i].display_elevmarkers = XmToggleButtonGetState(sui.etoggleall) ? ON : OFF;
	g[cg].flowt[i].em[curadcircem].display = XmToggleButtonGetState(sui.etoggle) ? ON : OFF;
	g[cg].flowt[i].em[curadcircem].emin = atof(XmTextGetString(sui.emin));
	g[cg].flowt[i].em[curadcircem].emax = atof(XmTextGetString(sui.emax));
    }
    update_display();
    if (flowt[curadcirc].active == ON) {
	display_image();
    }
}

/*
 * freshen up the ADCIRC setup panel
 */
void update_vector_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int flowno;
    char buf[256];
    if (vui.top) {
	SetChoice(vui.flow, curadcirc);
	flowno = curadcirc;
	SetChoice(vui.color, g[cg].flowt[flowno].p.color);
	switch (g[cg].flowt[flowno].display) {
	case OFF:
	    SetChoice(vui.display, 0);
	    break;
	case NODES:
	    SetChoice(vui.display, 1);
	    break;
	case CENTER:
	    SetChoice(vui.display, 2);
	    break;
	case ELLIPSE:
	    SetChoice(vui.display, 3);
	    break;
	}
	XmToggleButtonSetState(vui.mag, g[cg].flowt[flowno].display_mag == ON, False);
	XmToggleButtonSetState(vui.wind, g[cg].flowt[flowno].display_wind == ON, False);
    }
}

static void accept_vector_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i, flowno, type, start, stop;
    char buf[256];
    flowno = GetChoice(vui.flow);
    if (flowno == MAXADCIRC) {
	start = 0;
	stop = MAXADCIRC - 1;
    } else {
	curadcirc = start = stop = flowno;
	if (flowt[curadcirc].active == OFF) {
	    sprintf(buf, "Warning, no data for ADCIRC flow %d", curadcirc + 1);
	    errwin(buf);
	}
    }
    type = GetChoice(vui.display);
    for (i = start; i <= stop; i++) {
	switch (type) {
	case 0:
	    g[cg].flowt[i].display = OFF;
	    break;
	case 1:
	    g[cg].flowt[i].display = NODES;
	    break;
	case 2:
	    g[cg].flowt[i].display = CENTER;
	    break;
	case 3:
	    g[cg].flowt[i].display = ELLIPSE;
	    break;
	}
	g[cg].flowt[i].p.color = GetChoice(vui.color);
	g[cg].flowt[i].display_mag = XmToggleButtonGetState(vui.mag) ? ON : OFF;
	g[cg].flowt[i].display_wind = XmToggleButtonGetState(vui.wind) ? ON : OFF;
    }
    update_display();
    if (flowt[curadcirc].active == ON) {
	display_image();
    }
}

/*
 * create the file selection dialog for ADCIRC 2d data
 */
typedef struct {
    Widget top;
    Widget file;
    Widget grid;
    Widget *flow;
    Widget *type;
    Widget format;
    Widget level;
    Widget nlevels;
    Widget sample;
    Widget start;
    Widget stop;
    Widget skip;
    Widget missing;
    Widget missingval;
    Widget append;
    Browser b;
    Browser bg;
} adcircUI;

static adcircUI aui;

void accept_adcirc(Widget w, XtPointer clientd, XtPointer calld);

void create_adcirc_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int i;
    Widget wbut, panel, rc, rc2, wtmp;
    setistop();
    if (!aui.top) {
	aui.top = XmCreateDialogShell(app_shell, "ADCIRC", NULL, 0);
	panel = XmCreateRowColumn(aui.top, "panel", NULL, 0);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	aui.b.file = CreateTextItem2(rc, 20, "Data file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & aui.b);
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	aui.bg.file = CreateTextItem2(rc, 20, "Grid file:");
	wbut = XtVaCreateManagedWidget("Browse...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_browser, (XtPointer) & aui.bg);
	XtManageChild(rc);

	aui.flow = CreatePanelChoice1(panel, "Read to data set: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < MAXADCIRC; i++) {
	    XtAddCallback(aui.flow[i + 2], XmNactivateCallback, (XtCallbackProc) set_curadcirc, (XtPointer) i);
	}
	aui.type = CreatePanelChoice1(panel, "File type: ", 4, "fort.63 (scalars)", "fort.64 (vectors)", "fort.63 and fort.64", 0, 0);
	aui.format = XtVaCreateManagedWidget("File is ASCII", xmToggleButtonWidgetClass, panel, NULL);
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	aui.sample = XtVaCreateManagedWidget("*Sample steps:", xmToggleButtonWidgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	aui.start = CreateTextItem2(rc, 4, "Start: ");
	xv_setstr(aui.start, "1");
	aui.stop = CreateTextItem2(rc, 4, "Stop: ");
	xv_setstr(aui.stop, "1");
	aui.skip = CreateTextItem2(rc, 4, "Skip: ");
	xv_setstr(aui.skip, "1");
	XtManageChild(rc);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	aui.missing = XtVaCreateManagedWidget("*Missing data:", xmToggleButtonWidgetClass, rc, NULL);
	aui.missingval = CreateTextItem2(rc, 4, "Value: ");
	xv_setstr(aui.missingval, "-99.0");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_adcirc, (XtPointer) & aui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) aui.top);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    XtRaise(aui.top);
}

void accept_adcirc(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256], file[2048], gfile[2048];
    adcircUI *q = (adcircUI *) clientd;
    int flowno, level, nlevels, filet = 0, format = 0;
    int sample, start = -1, stop = 1, skip = 1;
    int missing;
    int append;
    double missingval;

    strcpy(file, (char *) xv_getstr(q->b.file));
    strcpy(gfile, (char *) xv_getstr(q->bg.file));
    filet = GetChoice(q->type);
    format = XmToggleButtonGetState(q->format);
    sample = XmToggleButtonGetState(q->sample);
    if (sample) {
	strcpy(buf, (char *) xv_getstr(q->start));
	start = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->stop));
	stop = atoi(buf) - 1;
	strcpy(buf, (char *) xv_getstr(q->skip));
	skip = atoi(buf);
    } else {
	start = -1;
	skip = 1;
    }
    missing = XmToggleButtonGetState(q->missing);
    if (missing) {
	strcpy(buf, (char *) xv_getstr(q->missingval));
	missingval = atof(buf);
    } else {
	missing = 0;
    }
    flowno = GetChoice(q->flow);
    if (!append && flowt[flowno].active == ON) {
	if (!yesno("Flow is active, kill it?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(q->top);
    switch (filet) {
    case 0:
	if (!format) {
	    if (!readbin_adcirc_elev(flowno, gfile, file)) {
		sprintf(buf, "Accept_adcirc(): error reading file %s", file);
		errwin(buf);
	    } else {
		set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		load_clock(ADCIRC, flowno);
		XtUnmanageChild(q->top);
		create_sadcirc_frame();
	    }
	} else {
	    if (!readflowt2d(flowno, filet, gfile, file)) {
		sprintf(buf, "Accept_adcirc(): error reading file %s", file);
		errwin(buf);
	    } else {
		set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		load_clock(ADCIRC, flowno);
		XtUnmanageChild(q->top);
		create_sadcirc_frame();
	    }
	}

	break;
    case 1:
	if (!format) {
	    if (!readbin_adcirc_flow(flowno, gfile, file, 0)) {
		sprintf(buf, "Error reading file %s", file);
		errwin(buf);
	    } else {
		set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		load_clock(ADCIRC, flowno);
		XtUnmanageChild(q->top);
		create_sadcirc_frame();
	    }
	} else {
	    if (!readflowt2d(flowno, filet, gfile, file)) {
		sprintf(buf, "Accept_adcirc(): error reading file %s", file);
		errwin(buf);
	    } else {
		set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		load_clock(ADCIRC, flowno);
		XtUnmanageChild(q->top);
		create_sadcirc_frame();
	    }
	}
	break;
    case 2:
	if (!format) {
	    if (!readbin_adcirc_elev(flowno, gfile, "fort.63")) {
		sprintf(buf, "Error reading file fort.63");
		errwin(buf);
	    } else {
		if (!readbin_adcirc_flow(flowno, gfile, "fort.64", 1)) {
		    sprintf(buf, "Error reading file fort.64");
		    errwin(buf);
		} else {
		    set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		    load_clock(ADCIRC, flowno);
		    XtUnmanageChild(q->top);
		    create_sadcirc_frame();
		}
	    }
	} else {
	    if (!readflowt2d(flowno, filet, gfile, file)) {
		sprintf(buf, "Accept_adcirc(): error reading file %s", file);
		errwin(buf);
	    } else {
		set_clock(0, flowt[flowno].start, flowt[flowno].stop, flowt[flowno].step, flowt[flowno].nsteps);
		load_clock(ADCIRC, flowno);
		XtUnmanageChild(q->top);
		create_sadcirc_frame();
	    }
	}
	break;
    default:
	break;
    }
    unset_wait_cursor(q->top);
}
