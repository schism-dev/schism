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
 * Read Time histories
 */

#ifndef lint
static char RCSid[] = "$Id: histwin.c,v 1.3 2003/08/31 16:49:44 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget hist_frame;
static Widget hist_panel;

static Widget save_frame;
static Widget shist_frame;
static Widget shist_panel;
void create_shist_frame(void);
void create_selhist_frame(void);
void shist_save_proc(void);
void hist_save_accept_proc(void);
static int save_hist(int ch, char *fname);
void update_hist(void);
void update_selhist(void);
void do_autoscalehistbox(void);

/*
 * Panel item declarations
 */
static Widget *shist_which_item;
static Widget shist_tmin_text_item;
static Widget shist_tmax_text_item;
static Widget shist_cmin_text_item;
static Widget shist_cmax_text_item;
static Widget shist_w_text_item;
static Widget shist_h_text_item;
static Widget *shist_precx_item;
static Widget *shist_precy_item;
static Widget *shist_color_item;
static Widget *shist_fill_item;
static Widget *shist_linew_item;
static Widget *shist_attach_item;

static Widget hist_files_list_item;
static Widget hist_toggle_item;
static Widget hist_togglemarker_item;
static Widget hist_teanl_item;
static Widget *hist_item;
static Widget *hist_filetype_item;
static Widget *hist_teanlcolor_item;
static Widget hist_adcirc_item;
static Widget *hist_adcirccolor_item;
static Widget hist_file_item;
static Widget *hist_filecolor_item;

static Widget hist_type_item;

void hist_done_proc();
void hist_accept_proc(void);
void shist_accept_proc(void);
void shist_view_proc(void);

void set_histloc(void);

void set_histlocCB(void)
{
    curhist = GetChoice(shist_which_item);
    set_histloc();
}

void set_curhist(int w, int cd)
{
    curhist = cd;
    update_hist();
    update_selhist();
}

void toggle_hist(void)
{
    g[cg].hbox[curhist].display = XmToggleButtonGetState(hist_toggle_item);
    g[cg].hbox[curhist].display_marker = XmToggleButtonGetState(hist_togglemarker_item);
}

/*
 * Create the hist Frame and the hist Panel
 */

extern Widget app_shell;

void create_hist_frame(void)
{
    int i;
    Arg wargs[8];
    Widget wbut, rc;
    setistop();
    if (!hist_frame) {
	hist_frame = XmCreateFileSelectionDialog(app_shell, "Time histories", NULL, 0);
	XtVaSetValues(hist_frame, XmNdirMask, XmStringCreate("*", charset), NULL);
	rc = XmCreateRowColumn(hist_frame, "rc", NULL, 0);
	hist_item = CreatePanelChoice1(rc, "Use time history marker: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXHISTMARKERS; i++) {
	    XtAddCallback(hist_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curhist, (XtPointer) i);
	}
	hist_filetype_item = CreatePanelChoice1(rc, "File type: ", 3, "Time domain", "Frequency domain", 0, 0);
	XtManageChild(rc);
	XtAddCallback(hist_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, hist_frame);
	XtAddCallback(hist_frame, XmNokCallback, (XtCallbackProc) hist_accept_proc, NULL);
    }
    XtRaise(hist_frame);
}

static char hist_fname[128];

void hist_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int type = GetChoice(hist_filetype_item);

    if (type) {
	type = FREQ;
    } else {
	type = TIME;
    }
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(hist_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(hist_fname, s);
    if (read_hist(curhist, type, hist_fname)) {
    } else {
    }
}

void create_shist_frame(void)
{
    Arg al[8];
    Widget wbut, fr, bb, sep, rc, rc2;
    int ac = 0, i;
    setistop();
    if (!shist_frame) {
	shist_frame = XmCreateDialogShell(app_shell, "Histories setup", NULL, 0);
	shist_panel = XmCreateRowColumn(shist_frame, "shistrc", NULL, 0);

	shist_which_item = CreatePanelChoice1(shist_panel, "History marker: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXHISTMARKERS; i++) {
	    XtAddCallback(shist_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curhist, (XtPointer) i);
	}
	hist_toggle_item = XtVaCreateManagedWidget("Display location", xmToggleButtonWidgetClass, shist_panel, NULL);
	hist_togglemarker_item = XtVaCreateManagedWidget("Display history marker", xmToggleButtonWidgetClass, shist_panel, NULL);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, shist_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);


	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	shist_tmin_text_item = CreateTextItem2(rc2, 10, "X min: ");
	shist_tmax_text_item = CreateTextItem2(rc2, 10, "X max: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	shist_cmin_text_item = CreateTextItem2(rc2, 10, "Y min: ");
	shist_cmax_text_item = CreateTextItem2(rc2, 10, "Y max: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	shist_w_text_item = CreateTextItem2(rc2, 6, "Width: ");
	shist_h_text_item = CreateTextItem2(rc2, 6, "Height: ");

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_autoscalehistbox, NULL);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, shist_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, NULL);
	shist_color_item = CreateColorChoice(bb, "Frame color:", 1);
	shist_linew_item = CreatePanelChoice2(bb, "Frame line width: ", 2, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

	shist_fill_item = CreateColorChoice(bb, "Fill color:", 1);
	shist_precx_item = CreatePanelChoice2(bb, "Precision X labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	shist_precy_item = CreatePanelChoice2(bb, "Precision Y labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	shist_attach_item = CreatePanelChoice2(bb, "Attach to: ", 1, 5, "SW", "SE", "NW", "NE", 0, 0);

	rc = XmCreateRowColumn(shist_panel, "shistrc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_histlocCB, NULL);
	wbut = XtVaCreateManagedWidget("*Specify location...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_histlocCB, NULL);
	wbut = XtVaCreateManagedWidget("Display...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_selhist_frame, NULL);
	XtManageChild(rc);
	XtManageChild(bb);

	rc = XmCreateRowColumn(shist_panel, "shistrc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) shist_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_hist_frame, NULL);
	wbut = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) shist_save_proc, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) shist_view_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, shist_frame);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(shist_panel);
    }
    XtRaise(shist_frame);
    update_hist();
}

void do_autoscalehistbox(void)
{
    char buf[256];
    if (timeclock.nsteps) {
	g[cg].hbox[curhist].wx1 = timeclock.t[0];
	g[cg].hbox[curhist].wx2 = timeclock.t[timeclock.nsteps - 1];
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wx1);
	XmTextSetString(shist_tmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wx2);
	XmTextSetString(shist_tmax_text_item, buf);
    }
}

void shist_view_proc(void)
{
    int i, itmp, cnt = 0;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEtimehistxXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"Time histories\"\n");
	fprintf(fp, "@xaxis label \"Time (seconds)\"\n");
	for (i = 0; i < MAXHISTMARKERS; i++) {
	    if (itmp = viewhist(fp, cg, i, timeclock.nsteps - 1, 0.0, cnt)) {
		cnt += itmp;
	    }
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    }
}

void shist_save_proc(void)
{
    int i;
    Arg wargs[8];
    Widget wbut, rc;
    setistop();
    if (!save_frame) {
	save_frame = XmCreateFileSelectionDialog(app_shell, "Save time histories", NULL, 0);
	XtVaSetValues(save_frame, XmNdirMask, XmStringCreate("*", charset), NULL);
	rc = XmCreateRowColumn(save_frame, "rc", NULL, 0);
	XtManageChild(rc);
	XtAddCallback(save_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, save_frame);
	XtAddCallback(save_frame, XmNokCallback, (XtCallbackProc) hist_save_accept_proc, NULL);
    }
    XtManageChild(save_frame);
}

void hist_save_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256], save_fname[256];
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(save_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(save_fname, s);
    if (save_hist(curhist, save_fname)) {
    } else {
    }
}

static int save_hist(int ch, char *fname)
{
    int i, cnt = 0;
    int c = get_current_step();
    double t = get_current_time();
    FILE *fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"Time histories\"\n");
	fprintf(fp, "@xaxis label \"Time (seconds)\"\n");
/*
	for (i = 0; i < MAXHISTMARKERS; i++) {
	    viewhist(fp, cg, i, timeclock.nsteps - 1, 0.0);
	}
*/
	viewhist(fp, cg, ch, timeclock.nsteps - 1, 0.0, cnt);
	fclose(fp);
	return 1;
    } else {
	return 0;
    }
}

void update_hist(void)
{
    char buf[256];
    if (shist_frame) {
	SetChoice(shist_which_item, curhist);

	XmToggleButtonSetState(hist_toggle_item, g[cg].hbox[curhist].display == ON, False);
	XmToggleButtonSetState(hist_togglemarker_item, g[cg].hbox[curhist].display_marker == ON, False);
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wx1);
	XmTextSetString(shist_tmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wx2);
	XmTextSetString(shist_tmax_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wy1);
	XmTextSetString(shist_cmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].hbox[curhist].wy2);
	XmTextSetString(shist_cmax_text_item, buf);
	sprintf(buf, "%.3lf", g[cg].hbox[curhist].vx);
	XmTextSetString(shist_w_text_item, buf);
	sprintf(buf, "%.3lf", g[cg].hbox[curhist].vy);
	XmTextSetString(shist_h_text_item, buf);

	SetChoice(shist_color_item, g[cg].hbox[curhist].p.color);
	SetChoice(shist_fill_item, g[cg].hbox[curhist].p.fillcol);
	SetChoice(shist_linew_item, g[cg].hbox[curhist].p.linew);
	SetChoice(shist_precx_item, g[cg].hbox[curhist].precx + 1);
	SetChoice(shist_precy_item, g[cg].hbox[curhist].precy + 1);
	SetChoice(shist_attach_item, g[cg].hbox[curhist].attach);
    }
}

void shist_accept_proc(void)
{
    int h;
    double x, y, locx, locy;
    x = y = locx = locy = 0.0;

    curhist = h = GetChoice(shist_which_item);

    g[cg].hbox[curhist].wx1 = atof((char *) xv_getstr(shist_tmin_text_item));
    g[cg].hbox[curhist].wx2 = atof((char *) xv_getstr(shist_tmax_text_item));
    g[cg].hbox[curhist].wy1 = atof((char *) xv_getstr(shist_cmin_text_item));
    g[cg].hbox[curhist].wy2 = atof((char *) xv_getstr(shist_cmax_text_item));
    g[cg].hbox[curhist].vx = atof((char *) xv_getstr(shist_w_text_item));
    g[cg].hbox[curhist].vy = atof((char *) xv_getstr(shist_h_text_item));

    g[cg].hbox[curhist].display = XmToggleButtonGetState(hist_toggle_item) ? ON : OFF;
    g[cg].hbox[curhist].display_marker = XmToggleButtonGetState(hist_togglemarker_item) ? ON : OFF;

    g[cg].hbox[curhist].p.color = GetChoice(shist_color_item);
    g[cg].hbox[curhist].p.fillcol = GetChoice(shist_fill_item);
    g[cg].hbox[curhist].p.linew = GetChoice(shist_linew_item);
    g[cg].hbox[curhist].precx = GetChoice(shist_precx_item) - 1;
    g[cg].hbox[curhist].precy = GetChoice(shist_precy_item) - 1;
    g[cg].hbox[curhist].attach = GetChoice(shist_attach_item);
}

static Widget selhist_frame, selhist_panel;
static Widget *teanl_color;
static Widget *time_color;
static Widget *adcirc_color;
static Widget *adcircflow_color;
static Widget timetoggle;
static Widget teanltoggle[MAXTEANL];
static Widget adcirctoggle[MAXADCIRC];
static Widget adcircflowtoggle[MAXADCIRC];

void selhist_accept_proc(void)
{
    int i;
    g[cg].hbox[curhist].thist = XmToggleButtonGetState(timetoggle);
    for (i = 0; i < MAXTEANL; i++) {
	g[cg].hbox[curhist].tp[i].color = GetChoice(teanl_color);
	g[cg].hbox[curhist].teanl[i] = XmToggleButtonGetState(teanltoggle[i]);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].hbox[curhist].ap[i].color = GetChoice(adcirc_color);
	g[cg].hbox[curhist].adcirc[i] = XmToggleButtonGetState(adcirctoggle[i]);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].hbox[curhist].apf[i].color = GetChoice(adcircflow_color);
	g[cg].hbox[curhist].adcircflow[i] = XmToggleButtonGetState(adcircflowtoggle[i]);
    }
}

void update_selhist(void)
{
    int i;
    if (selhist_frame) {
	XmToggleButtonSetState(timetoggle, g[cg].hbox[curhist].thist, False);
	for (i = 0; i < MAXTEANL; i++) {
	    XmToggleButtonSetState(teanltoggle[i], g[cg].hbox[curhist].teanl[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcirctoggle[i], g[cg].hbox[curhist].adcirc[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcircflowtoggle[i], g[cg].hbox[curhist].adcircflow[i], False);
	}
    }
}

void create_selhist_frame(void)
{
    Widget wbut, fr, rc3, sep, rc, rc2;
    int i;
    setistop();
    if (!selhist_frame) {
	selhist_frame = XmCreateDialogShell(app_shell, "Time history selection", NULL, 0);
	selhist_panel = XmCreateRowColumn(selhist_frame, "selhistrc", NULL, 0);


	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selhist_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("Time history:", xmLabelGadgetClass, rc3, NULL);
	timetoggle = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	time_color = CreateColorChoice(rc3, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selhist_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("TEANL elevations:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXTEANL; i++) {
	    sprintf(buf, "%1d", i + 1);
	    teanltoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
	teanl_color = CreateColorChoice(rc3, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selhist_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ADCIRC elevations:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcirctoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
	adcirc_color = CreateColorChoice(rc3, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selhist_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);

	XtVaCreateManagedWidget("ADCIRC magnitudes of velocity:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcircflowtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
	adcircflow_color = CreateColorChoice(rc3, "Color:", 1);

	rc = XmCreateRowColumn(selhist_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) selhist_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, selhist_frame);
	XtManageChild(rc);
	XtManageChild(rc3);
	XtManageChild(selhist_panel);
    }
    XtRaise(selhist_frame);
    update_selhist();
}
