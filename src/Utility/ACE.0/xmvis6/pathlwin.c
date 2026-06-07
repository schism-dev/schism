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
 * Set the properties of pathlines and drogues
 *
 */

#ifndef lint
static char RCSid[] = "$Id: pathlwin.c,v 1.3 2003/12/11 16:36:21 pturner Exp $";

#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget pathl_frame;
static Widget pathl_panel;
static Widget pathls_space_item;

static Widget *pathl_which_item;

void pathl_Done_notify_proc(void);
void pathl_file_notify_proc(void);
void pathl_files_notify_proc();
void create_pathls_frame(void);

void update_pathls(void);
void define_color_inregion(void);

void set_curdrog(Widget w, int cd)
{
    curdrog = cd;
    update_pathls();
}

void create_pathl_frame(void)
{
    Arg wargs[5];
    int i, n = 0;
    Widget rc;

    setistop();
    if (!pathl_frame) {
	pathl_frame = XmCreateFileSelectionDialog(app_shell, "pathl_frame", NULL, 0);
	XtVaSetValues(pathl_frame, XmNtitle, "Read Pathlines", XmNdirMask, XmStringCreate("*.pth", charset), NULL);
	rc = XmCreateRowColumn(pathl_frame, "rc", NULL, 0);
	pathl_which_item = CreatePanelChoice1(rc, "Read to Pathline: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXPATHLINES; i++) {
	    XtAddCallback(pathl_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curdrog, (XtPointer) i);
	}

	XtManageChild(rc);

	XtAddCallback(pathl_frame, XmNcancelCallback, (XtCallbackProc) pathl_Done_notify_proc, NULL);
	XtAddCallback(pathl_frame, XmNokCallback, (XtCallbackProc) pathl_file_notify_proc, NULL);
    }
    XtRaise(pathl_frame);
}

void pathl_Done_notify_proc(void)
{
    XtUnmanageChild(pathl_frame);
}

void pathl_file_notify_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(pathl_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    set_wait_cursor(pathl_frame);
    if (!readdrogues(curdrog, s, -1, 0, 0)) {
	sprintf(buf, "Error reading file %s", s);
	errwin(buf);
    } else {
	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop, drogues[curdrog].step, drogues[curdrog].nsteps);
	load_clock(DROGUES, curdrog);
    }
    unset_wait_cursor(pathl_frame);
    XtUnmanageChild(pathl_frame);
    create_pathls_frame();
}

static Widget pathls_frame;
static Widget pathls_panel;
void update_pathls(void);

int drog_curcolor = 1;
int drog_cursym = 2;

/*
 * Panel item declarations
 */
static Widget *pathls_which_item;
static Widget *pathls_color_item;
static Widget *pathls_symbol_item;
static Widget *pathls_colorapplyto_item;
static Widget pathls_toggle_item;
static Widget pathls_num_item;
static Widget pathls_depth_item;
static Widget pathls_streams_item;
static Widget pathls_id_item;
static Widget *pathls_linew_item;
static Widget *pathls_lines_item;

/*
 * Event and Notify proc declarations
 */

void pathls_pick_proc(void);
void pathls_region_proc(void);
void pathls_accept_proc(void);

void toggle_drogues();
void set_curdrog(Widget w, int cd);

void create_pathls_frame(void)
{
    Arg al[8];
    Widget wbut, fr, bb, sep, rc;
    int ac, i, j;

    setistop();
    if (!pathls_frame) {
	pathls_frame = XmCreateDialogShell(app_shell, "Drogues setup", NULL, 0);
	handle_close(pathls_frame);
	pathls_panel = XmCreateRowColumn(pathls_frame, "droguesrc", NULL, 0);
	pathls_which_item = CreatePanelChoice1(pathls_panel, "Use Drogues: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXPATHLINES; i++) {
	    XtAddCallback(pathls_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curdrog, (XtPointer) i);
	}

	fr = XmCreateFrame(pathls_panel, "fr", NULL, 0);
	bb = XmCreateRowColumn(fr, "rc", NULL, 0);
	pathls_toggle_item = XtVaCreateManagedWidget("Display drogues", xmToggleButtonWidgetClass, bb, NULL);
	pathls_num_item = XtVaCreateManagedWidget("Display drogue numbers", xmToggleButtonWidgetClass, bb, NULL);
	pathls_depth_item = XtVaCreateManagedWidget("Display drogue depths", xmToggleButtonWidgetClass, bb, NULL);
	pathls_id_item = XtVaCreateManagedWidget("*Display drogue IDs", xmToggleButtonWidgetClass, bb, NULL);
	pathls_streams_item = XtVaCreateManagedWidget("Display path lines", xmToggleButtonWidgetClass, bb, NULL);
	pathls_space_item = XtVaCreateManagedWidget("Display connection in space", xmToggleButtonWidgetClass, bb, NULL);
	XtManageChild(bb);
	XtManageChild(fr);

	fr = XmCreateFrame(pathls_panel, "fr", NULL, 0);
	bb = XmCreateRowColumn(fr, "rc", NULL, 0);
	pathls_color_item = CreateColorChoice(bb, "Color:", 1);
	pathls_linew_item = CreatePanelChoice1(bb, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
	pathls_lines_item = CreatePanelChoice1(bb, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
	pathls_symbol_item = CreatePanelChoice1(bb, "Symbol: ", 13, "None", "Dot", "Circle", "Square", "Diamond", "Triangle up", "Triangle left", "Triangle down", "Triangle right", "Plus", "X", "Star", 0, 0);
	pathls_colorapplyto_item = CreatePanelChoice1(bb, "*Apply color to: ", 4, "All drogues", "Drogues in region", "Pick drogues", 0, 0);
	rc = XmCreateRowColumn(bb, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Region", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) pathls_region_proc, NULL);
	wbut = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) pathls_pick_proc, NULL);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(fr);

	rc = XmCreateRowColumn(pathls_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) pathls_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Files...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_pathl_frame, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, pathls_frame);
	XtManageChild(rc);

	XtManageChild(pathls_panel);
/*
	for (i = 0; i < MAXPATHLINES; i++) {
	    g[cg].drogues[i].color = (int *) malloc(MAXDROGS * sizeof(int));
	    g[cg].drogues[i].sym = (int *) malloc(MAXDROGS * sizeof(int));
	    g[cg].drogues[i].size = (double *) malloc(MAXDROGS * sizeof(double));
	    for (j = 0; j < MAXDROGS; j++) {
		g[cg].drogues[i].color[j] = 1;
		g[cg].drogues[i].sym[j] = 3;
		g[cg].drogues[i].size[j] = 1.0;
	    }
	}
*/
    }
    update_pathls();
    XtRaise(pathls_frame);
}

void update_pathls(void)
{
    char buf[256];

    if (pathls_frame) {
	XmToggleButtonSetState(pathls_toggle_item, g[cg].drogues[curdrog].display == ON, False);
	XmToggleButtonSetState(pathls_num_item, g[cg].drogues[curdrog].display_type == NUMBER, False);
	XmToggleButtonSetState(pathls_depth_item, g[cg].drogues[curdrog].display_type == DEPTH, False);
	XmToggleButtonSetState(pathls_id_item, g[cg].drogues[curdrog].display_id == NUMBER, False);
	XmToggleButtonSetState(pathls_streams_item, g[cg].drogues[curdrog].display_streaml == ON, False);
	XmToggleButtonSetState(pathls_space_item, g[cg].drogues[curdrog].display_connect == ON, False);
	SetChoice(pathls_symbol_item, drog_cursym);

	SetChoice(pathls_color_item, g[cg].drogues[curdrog].p.color);
	SetChoice(pathls_linew_item, g[cg].drogues[curdrog].p.linew - 1);
	SetChoice(pathls_lines_item, g[cg].drogues[curdrog].p.lines - 1);
    }
}

void pathls_accept_proc(void)
{
    char buf[256];
    int j, sym, applyto;

    g[cg].drogues[curdrog].display = XmToggleButtonGetState(pathls_toggle_item) ? ON : OFF;
    g[cg].drogues[curdrog].display_type = XmToggleButtonGetState(pathls_num_item) ? NUMBER : OFF;
    g[cg].drogues[curdrog].display_type = XmToggleButtonGetState(pathls_depth_item) ? DEPTH : OFF;
    g[cg].drogues[curdrog].display_id = XmToggleButtonGetState(pathls_id_item) ? ON : OFF;
    g[cg].drogues[curdrog].display_streaml = XmToggleButtonGetState(pathls_streams_item) ? ON : OFF;
    g[cg].drogues[curdrog].display_connect = XmToggleButtonGetState(pathls_space_item) ? ON : OFF;
    g[cg].drogues[curdrog].p.color = GetChoice(pathls_color_item);
    g[cg].drogues[curdrog].p.linew = GetChoice(pathls_linew_item) + 1;
    g[cg].drogues[curdrog].p.lines = GetChoice(pathls_lines_item) + 1;
    sym = GetChoice(pathls_symbol_item);
    applyto = GetChoice(pathls_colorapplyto_item);
    switch (applyto) {
    case 0:
/*
	for (j = 0; j < MAXDROGS; j++) {
	    g[cg].drogues[curdrog].color[j] = g[cg].drogues[curdrog].p.color;
	    g[cg].drogues[curdrog].sym[j] = sym;
	}
*/
	break;
    case 1:
	define_color_inregion();
	break;
    case 2:
	break;
    }
}

void pathls_pick_proc(void)
{
    set_restart();
    set_action(0);
    drog_curcolor = GetChoice(pathls_color_item);
    set_action(PICK_DROGUE_COLOR);
}

void pathls_region_proc(void)
{
    set_restart();
    set_action(0);
    drog_curcolor = GetChoice(pathls_color_item);
    do_select_single_region();
/*
   set_action(PICK_DROGUE_REGION);
 */
}

void define_color_inregion(void)
{
    int i;
    int n;
    double x, y;
    if (region_flag) {
	for (i = 0; i < drogues[0].p[0].npts; i++) {
	    x = drogues[0].p[0].x[i];
	    y = drogues[0].p[0].y[i];
	    if (inregion2(regionx, regiony, nregion, x, y)) {
/*
printf("x, y == %lf %lf drog_curcolor = %d\n");
		g[cg].drogues[curdrog].color[i] = drog_curcolor;
*/
	    }
	}
    }
}
