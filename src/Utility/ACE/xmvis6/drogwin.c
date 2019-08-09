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
 * Run drogues or read data
 */

#ifndef lint
static char RCSid[] = "$Id";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget drogues_frame;
static Widget drogues_panel;
void create_drogues_frame(void);
void create_pathl_frame(void);
void update_drogues(void);

extern int drog_curcolor;

/*
 * Panel item declarations
 */
static Widget drogues_start_text_item;
static Widget drogues_step_text_item;
static Widget drogues_flow_text_item;
static Widget *drogues_which_item;
static Widget *drogues_color_item;
static Widget drogues_grid_text_item;
static Widget drogues_pathl_text_item;
static Widget drogues_nsteps_text_item;
static Widget drogues_toggle_item;
static Widget drogues_num_item;
static Widget drogues_streams_item;
static Widget drogues_model_item;
static Widget drogues_backtrack_item;
static Widget *drogues_track_item;
static Widget *drogues_pos_item;

/*
 * Event and Notify proc declarations
 */

void drogues_doelev_proc();
void drogues_place_proc(void);
void drogues_delete_proc(void);
void drogues_pick_proc(void);
void drogues_load_proc(void);
void drogues_run_proc(void);
void drogues_accept_proc(void);
void drogues_done_proc(void);
void drogues_reset_proc(void);

void toggle_drogues();
void set_curdrog(Widget w, int cd);

void create_drogues_frame(void)
{
    Arg al[8];
    Widget wbut, fr, bb, sep, rc;
    int ac, i;
    setistop();
    if (!drogues_frame) {
	drogues_frame = XmCreateDialogShell(app_shell, "Drogues setup", NULL, 0);
	drogues_panel = XmCreateRowColumn(drogues_frame, "droguesrc", NULL, 0);
	drogues_which_item = CreatePanelChoice1(drogues_panel, "Use Drogues: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < 5; i++) {
	    XtAddCallback(drogues_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curdrog, (XtPointer) i);
	}

	fr = XmCreateFrame(drogues_panel, "fr", NULL, 0);
	bb = XmCreateRowColumn(fr, "droguesrc", NULL, 0);

	drogues_start_text_item = CreateTextItem2(bb, 5, "Start time (days):");
	drogues_step_text_item = CreateTextItem2(bb, 5, "Length of run (days):");
	drogues_nsteps_text_item = CreateTextItem2(bb, 5, "Number of time steps:");
	drogues_grid_text_item = CreateTextItem2(bb, 10, "Study: ");
	drogues_flow_text_item = CreateTextItem2(bb, 10, "Flow run: ");
	drogues_pathl_text_item = CreateTextItem2(bb, 10, "Pathline run: ");

	drogues_pos_item = CreatePanelChoice1(bb, "Positioning: ", 4, "All nodes", "Individual", "Distribution", 0, 0);
	drogues_track_item = CreatePanelChoice1(bb, "Tracking by: ", 5, "Euler", "RK2", "RK4", "RK5", 0, 0);
	drogues_backtrack_item = XtVaCreateManagedWidget("Back track", xmToggleButtonWidgetClass, bb, NULL);
	rc = XmCreateRowColumn(bb, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_place_proc, NULL);

	wbut = XtVaCreateManagedWidget("Delete", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_delete_proc, NULL);

	wbut = XtVaCreateManagedWidget("*Pick", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_pick_proc, NULL);

	wbut = XtVaCreateManagedWidget("Load", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_load_proc, NULL);

	wbut = XtVaCreateManagedWidget("Reset", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_reset_proc, NULL);

	wbut = XtVaCreateManagedWidget("Run", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_run_proc, NULL);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(fr);

	fr = XmCreateFrame(drogues_panel, "fr", NULL, 0);
	bb = XmCreateRowColumn(fr, "droguesrc", NULL, 0);
	drogues_toggle_item = XtVaCreateManagedWidget("Display drogues", xmToggleButtonWidgetClass, bb, NULL);
	drogues_num_item = XtVaCreateManagedWidget("Display drogue numbers", xmToggleButtonWidgetClass, bb, NULL);
	drogues_streams_item = XtVaCreateManagedWidget("Streamlines", xmToggleButtonWidgetClass, bb, NULL);
	drogues_color_item = CreateColorChoice(drogues_panel, "Color:", 1);
	XtManageChild(bb);
	XtManageChild(fr);

	rc = XmCreateRowColumn(drogues_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Files...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_pathl_frame, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) drogues_done_proc, NULL);
	XtManageChild(rc);

	XtManageChild(drogues_panel);
    }
    XtRaise(drogues_frame);
    update_drogues();
}

double drogx[1000], drogy[1000];
int ndrogs;

int find_nearest_drog(int drno, double x, double y)
{
    int i, ind;
    double radius = 1e307;
    double tmp;

    for (i = 0; i < ndrogs; i++) {
	tmp = hypot((x - drogx[i]), (y - drogy[i]));
	if (tmp < radius) {
	    radius = tmp;
	    ind = i;
	}
    }
    return ind;
}

int find_nearest_drogue(int drno, double x, double y)
{
    int i, ind, dnum = -1;
    double radius = 1e307;
    double tmp;
    int cstep = get_current_step();
    for (i = 0; i < drogues[drno].p[cstep].npts; i++) {
	tmp = hypot((x - drogues[drno].p[cstep].x[i]), (y - drogues[drno].p[cstep].y[i]));
	if (tmp < radius) {
	    radius = tmp;
	    ind = i;
	    dnum = drogues[drno].p[cstep].drnum[i] - 1;
	}
    }
    return dnum;
}


void delete_drog(int drno, int dr)
{
    int i, cnt = 0;
    if (dr != ndrogs - 1) {
	for (i = drno; i < ndrogs - 1; i++) {
	    drogx[i] = drogx[i + 1];
	    drogy[i] = drogy[i + 1];
	}
    }
    ndrogs--;
}

void pick_drog(int drno, int dr)
{
    /*g[cg].drogues[drno].color[dr] = GetChoice(drogues_color_item); */
}

void set_drogue_color(int drno, int dr, int c)
{
    /*g[cg].drogues[drno].color[dr] = c; */
}

void drogues_reset_proc(void)
{
    ndrogs = 0;
}

void drogues_run_proc(void)
{
    char buf[256];
    char study[128];
    char flowrun[128];
    char pathline[128];
    char fname[256];
    double start, days, t1, t2;
    int i, steps;
    FILE *fp;
    if (ndrogs == 0) {
	errwin("No drogues, run canceled");
	return;
    }
    fp = fopen("tmp.drogs", "w");
    if (fp == NULL) {
	errwin("Couldn't open file for run");
	return;
    }
    strcpy(study, (char *) xv_getstr(drogues_grid_text_item));
    fprintf(fp, "%s\n", study);
    strcpy(pathline, (char *) xv_getstr(drogues_pathl_text_item));
    fprintf(fp, "%s\n", pathline);
    strcpy(flowrun, (char *) xv_getstr(drogues_flow_text_item));
    fprintf(fp, "%s\n", flowrun);
    fclose(fp);

    strcpy(flowrun, study);
    strcat(flowrun, (char *) xv_getstr(drogues_flow_text_item));
    strcat(flowrun, ".tct");

    strcpy(pathline, study);
    strcat(pathline, (char *) xv_getstr(drogues_pathl_text_item));
    strcat(pathline, ".din");

    strcat(study, ".gr3");

    printf(">%s< >%s< >%s<\n", study, flowrun, pathline);
    fp = fopen(pathline, "w");
    if (fp == NULL) {
	errwin("Couldn't open file for run");
	return;
    }
    fprintf(fp, "From xmvis\n");
    fprintf(fp, "Run 00\n");
    fprintf(fp, "0.001 0.00004\n");
    start = atof((char *) xv_getstr(drogues_start_text_item));
    steps = atoi((char *) xv_getstr(drogues_nsteps_text_item));
    days = atof((char *) xv_getstr(drogues_step_text_item));
    fprintf(fp, "%lf %d %lf\n", days, steps, start);
    fprintf(fp, "%d\n", ndrogs);
    for (i = 0; i < ndrogs; i++) {
	fprintf(fp, "%lf %lf %d\n", drogx[i], drogy[i], i);
    }
    fclose(fp);
    sprintf(buf, "drogues < tmp.drogs");
    system(buf);
}

void drogues_done_proc(void)
{
    XtUnmanageChild(drogues_frame);
}

void update_drogues(void)
{
    char buf[256];
    if (drogues_frame) {
	sprintf(buf, "%lf", drogues[curdrog].start);
	xv_setstr(drogues_start_text_item, buf);
/* TODO
    sprintf(buf, "%lf", drogues[curdrog].stop);
    xv_setstr(drogues_stop_text_item, buf);
*/
	sprintf(buf, "%d", drogues[curdrog].nsteps);
	xv_setstr(drogues_nsteps_text_item, buf);
	XmToggleButtonSetState(drogues_toggle_item, g[cg].drogues[curdrog].display == ON, False);
	XmToggleButtonSetState(drogues_num_item, g[cg].drogues[curdrog].display_type == NUMBER, False);
	XmToggleButtonSetState(drogues_streams_item, g[cg].drogues[curdrog].display_streaml == ON, False);
    }
}

void drogues_accept_proc(void)
{
    char buf[256];
    g[cg].drogues[curdrog].display = XmToggleButtonGetState(drogues_toggle_item) ? ON : OFF;
    g[cg].drogues[curdrog].display_type = XmToggleButtonGetState(drogues_num_item) ? NUMBER : OFF;
    g[cg].drogues[curdrog].display_streaml = XmToggleButtonGetState(drogues_streams_item) ? ON : OFF;
    g[cg].drogues[curdrog].p.color = GetChoice(drogues_color_item);
}

void drogues_load_proc(void)
{
    int i, c = get_current_step();
    int drind;
    for (i = 0; i < drogues[curdrog].p[c].npts; i++) {
	drind = drogues[curdrog].p[c].drnum[i] - 1;
	drogx[i] = drogues[curdrog].p[0].x[drind];
	drogy[i] = drogues[curdrog].p[0].y[drind];
    }
    ndrogs = drogues[curdrog].p[c].npts;
}

void drogues_place_proc(void)
{
    set_restart();
    set_action(0);
    set_action(PLACE_DROGUE);
}

void drogues_delete_proc(void)
{
    set_restart();
    set_action(0);
    set_action(DELETE_DROGUE);
}

void drogues_pick_proc(void)
{
    set_restart();
    set_action(0);
    drog_curcolor = GetChoice(drogues_color_item);
    set_action(PICK_DROGUE_COLOR);
}
