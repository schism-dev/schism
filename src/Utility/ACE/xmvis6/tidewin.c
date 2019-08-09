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
 * Tide stations
 */

#ifndef lint
static char RCSid[] = "$Id: tidewin.c,v 1.2 2003/07/24 15:44:07 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

static Widget tide_frame;
static Widget tide_panel;

static Widget save_frame;
static Widget stide_frame;
static Widget stide_panel;

int ReadTideStations(char *tide_fname);

void create_stide_frame(void);
void create_seltide_frame(void);
void stide_save_proc(void);
static void tide_save_accept_proc(void);
static void stide_rms_proc(void);
static void stide_view2_proc(void);
static void set_curfreq(int w, int cd, XtPointer cb);
static int save_tide(int ch, char *fname);
void update_tide(int which);
void update_seltide(void);
void do_autoscaletidebox(void);
void create_tideprops_frame(void);

/*
 * Panel item declarations
 */
static Widget stide_list_item;
static Widget *stide_applyto_item;
static Widget amppha_list_item;
static Widget stide_tmin_text_item;
static Widget stide_tmax_text_item;
static Widget stide_cmin_text_item;
static Widget stide_cmax_text_item;
static Widget stide_vx_text_item;
static Widget stide_vy_text_item;
static Widget *stide_precx_item;
static Widget *stide_precy_item;
static Widget *stide_color_item;
static Widget *stide_fill_item;
static Widget *stide_linew_item;
static Widget *stide_attach_item;
static Widget *stide_props_applyto;

static Widget tide_files_list_item;
static Widget tide_toggle_item;
static Widget tide_togglemarker_item;
static Widget *tide_displayamppha_item;
static Widget *tide_usefreq_item;
static Widget tide_teanl_item;
static Widget *tide_item;
static Widget *tide_filetype_item;
static Widget *tide_teanlcolor_item;
static Widget tide_adcirc_item;
static Widget *tide_adcirccolor_item;
static Widget tide_file_item;
static Widget *tide_filecolor_item;

static Widget tide_type_item;

static void tide_done_proc();
static void tide_accept_proc(void);
static void stide_accept_proc(void);
static void stide_view_proc(void);

void find_nearest_tidestation(double wx, double wy, int *ind);

static int curtide = 0;
static int curfreq = 0;
static int curamppha = 0;

static void set_tidepick(void)
{
    set_action(0);
    set_action(QUERY_TIDESTATION);
}

static void set_tidelocCB(void)
{
    set_action(0);
    set_action(PLACE_TIDESTATION);
}

void set_current_tidestation(int ind)
{
    if (ind >= 0) {
	curtide = ind;
	XmListSetPos(stide_list_item, curtide + 1);
	XmListSelectPos(stide_list_item, curtide + 1, False);
	update_tide(0);
    }
}

void set_current_tidestationloc(double wx, double wy)
{
    g[cg].tidestat[curtide].locx = wx;
    g[cg].tidestat[curtide].locy = wy;
}

static void set_curtide(int w, int cd, XtPointer cb)
{
    int pos;
    char buf[256], name[512];
    XmString xms;
    XmString *s, cs;
    int *pos_list;
    int i, j, pos_cnt, cnt;
    char *cstr;
    XmListCallbackStruct *cbs = (XmListCallbackStruct *) cb;

    curtide = cbs->item_position - 1;

    if (pos = XmListGetSelectedPos(stide_list_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(stide_list_item, XmNselectedItemCount, &cnt, XmNselectedItems, &s, NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(name, cstr);
	    XtFree(cstr);
	    XmStringFree(cs);
	    update_tide(0);
	}
    }
}

static void toggle_tide(void)
{
    g[cg].tidestat[curtide].display = XmToggleButtonGetState(tide_toggle_item) ? ON : OFF;
    g[cg].tidestat[curtide].display_marker = XmToggleButtonGetState(tide_togglemarker_item) ? ON : OFF;
    g[cg].tidestat[curtide].display_ampphase = GetChoice(tide_displayamppha_item);
    g[cg].tidestat[curtide].use_freq = GetChoice(tide_usefreq_item);
}

/*
 * Create the tide Frame and the tide Panel
 */

extern Widget app_shell;
static Widget freq_list_item;

void create_tide_frame(void)
{
    int i;
    Arg wargs[8];
    Widget wbut, rc;
    setistop();
    if (!tide_frame) {
	tide_frame = XmCreateFileSelectionDialog(app_shell, "Tide stations", NULL, 0);
	XtVaSetValues(tide_frame, XmNdirMask, XmStringCreate("*", charset), NULL);
	rc = XmCreateRowColumn(tide_frame, "rc", NULL, 0);
	tide_filetype_item = CreatePanelChoice1(rc, "*File type: ", 3, "Time domain", "Frequency domain", 0, 0);
	XtManageChild(rc);
	XtAddCallback(tide_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, tide_frame);
	XtAddCallback(tide_frame, XmNokCallback, (XtCallbackProc) tide_accept_proc, NULL);
    }
    XtRaise(tide_frame);
}

static char tide_fname[128];

static void tide_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int type = GetChoice(tide_filetype_item);

    if (type) {
	type = FREQ;
    } else {
	type = TIME;
    }
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(tide_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(tide_fname, s);
    if (ReadTideStations(tide_fname)) {
    } else {
	XtUnmanageChild(tide_frame);
	create_stide_frame();
    }
}

void create_stide_frame(void)
{
    Arg al[8];
    Widget wbut, fr, bb, sep, rc, lab, rc1, rc2;
    int ac = 0, i;
    int got_one = 0;
    setistop();
    if (ntidestat == 0) {
	create_tide_frame();
	return;
    }
/*
    for (i = 0; i < MAXTEANL; i++) {
        if (flowf[i].active == ON) {
            got_one = 1;
        }
    }
    if (!got_one) {
        create_tide_frame();
        return;
    }
*/

    if (!stide_frame) {
	stide_frame = XmCreateDialogShell(app_shell, "Tide stations", NULL, 0);
	stide_panel = XmCreateRowColumn(stide_frame, "stiderc", NULL, 0);

	XtSetArg(al[0], XmNlistSizePolicy, XmRESIZE_IF_POSSIBLE);
	XtSetArg(al[1], XmNvisibleItemCount, 5);
	XtSetArg(al[2], XmNselectionPolicy, XmSINGLE_SELECT);

	lab = XmCreateLabel(stide_panel, "Tide stations:", NULL, 0);
	XtManageChild(lab);
	stide_list_item = XmCreateScrolledList(stide_panel, "list", al, 3);
	XtAddCallback(stide_list_item, XmNsingleSelectionCallback, (XtCallbackProc) set_curtide, 0);
	XtManageChild(stide_list_item);

	tide_toggle_item = XtVaCreateManagedWidget("Display location", xmToggleButtonWidgetClass, stide_panel, NULL);
	tide_togglemarker_item = XtVaCreateManagedWidget("Display tide station marker", xmToggleButtonWidgetClass, stide_panel, NULL);
	tide_displayamppha_item = CreatePanelChoice1(stide_panel, "Display: ", 4, "None", "All frequencies", "*Selected", 0, 0);
	tide_usefreq_item = CreatePanelChoice1(stide_panel, "Use: ", 4, "All frequencies", "*Same as simulation", "*Selected", 0, 0);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, stide_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);

	stide_tmin_text_item = CreateTextItem2(bb, 15, "X min: ");
	stide_tmax_text_item = CreateTextItem2(bb, 15, "X max: ");
	stide_cmin_text_item = CreateTextItem2(bb, 15, "Y min: ");
	stide_cmax_text_item = CreateTextItem2(bb, 15, "Y max: ");

	XtSetArg(al[0], XmNlistSizePolicy, XmRESIZE_IF_POSSIBLE);
	XtSetArg(al[1], XmNvisibleItemCount, 5);
	XtSetArg(al[2], XmNselectionPolicy, XmMULTIPLE_SELECT);

	lab = XmCreateLabel(bb, "Use frequency:", NULL, 0);
	XtManageChild(lab);
	freq_list_item = XmCreateScrolledList(bb, "list", al, 3);
	XtAddCallback(freq_list_item, XmNsingleSelectionCallback, (XtCallbackProc) set_curfreq, 0);
	XtManageChild(freq_list_item);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, stide_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);

	rc2 = XtVaCreateManagedWidget("rc2", xmRowColumnWidgetClass, bb, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_autoscaletidebox, NULL);
	wbut = XtVaCreateManagedWidget("Props...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_tideprops_frame, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_tidelocCB, NULL);
	wbut = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_tidepick, NULL);
	wbut = XtVaCreateManagedWidget("Display...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_seltide_frame, NULL);

	rc = XmCreateRowColumn(bb, "stiderc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	stide_applyto_item = CreatePanelChoice1(rc, "Apply to: ", 4, "Current selection", "All", "All in region", 0, 0);
	XtManageChild(rc);

	rc = XmCreateRowColumn(bb, "stiderc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stide_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_tide_frame, NULL);
	wbut = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stide_save_proc, NULL);
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stide_view_proc, NULL);
	wbut = XtVaCreateManagedWidget("Cmp...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stide_view2_proc, NULL);
	wbut = XtVaCreateManagedWidget("RMS...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stide_rms_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, stide_frame);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(stide_panel);
    }
    XtRaise(stide_frame);
    update_tide(1);
}

void do_autoscaletidebox(void)
{
    double amin, amax;
    char buf[256];
    if (timeclock.nsteps) {
	switch (GetChoice(stide_applyto_item)) {
	case 0:
	    break;
	case 1:
	    break;
	case 2:
	    break;
	}
	g[cg].tidestat[curtide].wx1 = timeclock.t[0];
	g[cg].tidestat[curtide].wx2 = timeclock.t[timeclock.nsteps - 1];

	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].wx1);
	XmTextSetString(stide_tmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].wx2);
	XmTextSetString(stide_tmax_text_item, buf);

	TideStationElevationMinMax(tidestat[curtide], &amin, &amax);

	g[cg].tidestat[curtide].wy1 = amin;
	g[cg].tidestat[curtide].wy2 = amax;

	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].wy1);
	XmTextSetString(stide_cmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].wy2);
	XmTextSetString(stide_cmax_text_item, buf);

    } else {
	errwin("Time clock not set, can't autoscale");
    }
}

static void stide_view_proc(void)
{
    int i, itmp, cnt = 0;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEtidesXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"%s\"\n", tidestat[curtide]->name);
	fprintf(fp, "@xaxis label \"Time (seconds)\"\n");
	if (itmp = viewtide(fp, cg, curtide, timeclock.nsteps - 1, 0.0, cnt)) {
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    }
}

static void stide_rms_proc(void)
{
    int i, itmp, cnt = 0;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACErmsXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"RMS errors\"\n");
	fprintf(fp, "@xaxis label \"Station\"\n");
	fprintf(fp, "@yaxis label \"RMS Error\"\n");
	fprintf(fp, "@xaxis  ticklabel type spec\n");
	fprintf(fp, "@xaxis  ticklabel layout spec\n");
	fprintf(fp, "@xaxis  ticklabel angle 333\n");
	fprintf(fp, "@xaxis  tick type spec\n");
	fprintf(fp, "@xaxis  tick spec %d\n", ntidestat);
	fprintf(fp, "@s0 symbol 12\n");
	fprintf(fp, "@s0 linewidth 0\n");
	for (i = 0; i < ntidestat; i++) {
	    fprintf(fp, "@xaxis  tick %d, %d\n", i, i + 1);
	    fprintf(fp, "@xaxis  ticklabel %d, \"%s\"\n", i, tidestat[i]->name);
	}

	if (itmp = viewrms(fp, cg, ntidestat, timeclock.nsteps - 1, 0.0, cnt)) {
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    }
}

static void stide_view2_proc(void)
{
    int i, itmp, cnt = 0;
    char tbuf[256];
    char combuf[256];
    int ind = -1;
    double x, y, amp[120], pha[120];
    char *fname;
    FILE *fp;
    find_element(0, x = g[0].tidestat[curtide].x, y = g[0].tidestat[curtide].y, &ind);
    if (ind >= 0) {
	strcpy(tbuf, "/tmp/ACEcmptideXXXXXX");
	mkstemp(tbuf);
	fname = tbuf;
	fp = fopen(fname, "w");
	if (fp != NULL) {
	    fprintf(fp, "@xaxis label \"Frequency\"\n");
	    fprintf(fp, "@yaxis label \"Amplitude\"\n");
	    GetTeanlElevAmpPhase(0, ind, x, y, amp, pha);
	    for (i = 0; i < tidestat[curtide]->nfreq; i++) {
		fprintf(fp, "%.12lf %.12lf\n", tidestat[curtide]->omega[i], tidestat[curtide]->elamp[i]);
	    }
	    fprintf(fp, "&\n");
	    for (i = 0; i < flowf[0].nfreq; i++) {
		fprintf(fp, "%.12lf %.12lf\n", flowf[0].omega[i], amp[i]);
	    }
	    fprintf(fp, "&\n");
	    fprintf(fp, "@with g1\n");
	    for (i = 0; i < tidestat[curtide]->nfreq; i++) {
		fprintf(fp, "%.12lf %.12lf\n", tidestat[curtide]->omega[i], fmod(tidestat[curtide]->elphase[i] * 180.0 / M_PI + 360.0, 360.0));
	    }
	    fprintf(fp, "&\n");
	    for (i = 0; i < flowf[0].nfreq; i++) {
		fprintf(fp, "%.12lf %.12lf\n", flowf[0].omega[i], fmod(pha[i] * 180.0 / M_PI + 360.0, 360.0));
	    }
	    fprintf(fp, "&\n");
	    fprintf(fp, "@title \"Station %s\"\n", tidestat[curtide]->name);
	    fprintf(fp, "@yaxis label \"Phase\"\n");
	    fprintf(fp, "@autoscale\n");
	    fclose(fp);
	    sprintf(combuf, " /usr/local/ace/bin/xmgr5 -arrange 2 1 -remove %s &", fname);
	    system(combuf);
	}
    }
}

void stide_save_proc(void)
{
    int i;
    Arg wargs[8];
    Widget wbut, rc;
    setistop();
    if (!save_frame) {
	save_frame = XmCreateFileSelectionDialog(app_shell, "Save tide stations", NULL, 0);
	XtVaSetValues(save_frame, XmNdirMask, XmStringCreate("*", charset), NULL);
	rc = XmCreateRowColumn(save_frame, "rc", NULL, 0);
	XtManageChild(rc);
	XtAddCallback(save_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, save_frame);
	XtAddCallback(save_frame, XmNokCallback, (XtCallbackProc) tide_save_accept_proc, NULL);
    }
    XtManageChild(save_frame);
}

static void tide_save_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256], save_fname[256];
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(save_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(save_fname, s);
    if (save_tide(curtide, save_fname)) {
    } else {
    }
}

static int save_tide(int ch, char *fname)
{
    int i, cnt = 0;
    int c = get_current_step();
    double t = get_current_time();
    FILE *fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"%s\"\n", tidestat[curtide]->name);
	fprintf(fp, "@xaxis label \"Time (seconds)\"\n");
	viewtide(fp, cg, curtide, timeclock.nsteps - 1, 0.0, cnt);
	fclose(fp);
	return 1;
    } else {
	return 0;
    }
}

void update_tide(int which)
{
    char buf[256];
    XmString xms;
    int i;
    double period;
    if (stide_frame) {
	if (g[cg].tidestat == NULL && ntidestat > 0) {
/*
	    g[cg].tidestat = (DisplayTideStation *) malloc(ntidestat * sizeof(DisplayTideStation));
*/
	    for (i = 0; i < ntidestat; i++) {
		InitDisplayTideStation(&g[cg].tidestat[i]);
		g[cg].tidestat[i].x = tidestat[i]->x;
		g[cg].tidestat[i].y = tidestat[i]->y;
		g[cg].tidestat[i].locx = tidestat[i]->x;
		g[cg].tidestat[i].locy = tidestat[i]->y;
	    }
	}
	if (which) {
	    XmListDeleteAllItems(stide_list_item);
	    for (i = 0; i < ntidestat; i++) {
		if (strlen(tidestat[i]->name) > 0) {
		    sprintf(buf, "%s", tidestat[i]->name);
		} else {
		    strcpy(buf, "Unknown station");
		}
		xms = XmStringCreateLtoR(buf, charset);
		XmListAddItemUnselected(stide_list_item, xms, 0);
		XmStringFree(xms);
	    }
	}
	XmListSelectPos(stide_list_item, curtide + 1, False);
	update_seltide();

	XmListDeleteAllItems(freq_list_item);
/*
	xms = XmStringCreateLtoR("ALL frequencies", charset);
	XmListAddItemUnselected(freq_list_item, xms, 0);
	XmStringFree(xms);
*/
	for (i = 0; i < tidestat[curtide]->nfreq; i++) {
	    if (strlen(tidestat[curtide]->freqname[i]) > 0) {
		if (tidestat[curtide]->omega[i] != 0.0) {
		    period = 2.0 * M_PI / (tidestat[curtide]->omega[i] * 3600.0);
		} else {
		    period = 999999.0;
		}
		sprintf(buf, "%s (%lf, %lf)", tidestat[curtide]->freqname[i], tidestat[curtide]->omega[i], period);
	    } else {
		strcpy(buf, "UNKNOWN");
	    }
	    xms = XmStringCreateLtoR(buf, charset);
	    XmListAddItemUnselected(freq_list_item, xms, 0);
	    XmStringFree(xms);
	}
	XmListSelectPos(freq_list_item, curfreq + 1, False);
    }
}

static void stide_accept_proc(void)
{
    int h, i, imin, imax;
    double x, y, locx, locy;
    x = y = locx = locy = 0.0;
    switch (GetChoice(stide_applyto_item)) {
    case 0:
	imin = imax = curtide;
	break;
    case 1:
	imin = 0;
	imax = ntidestat - 1;
	break;
    case 2:
	if (nregion) {
	    double xtmp, ytmp;
	    for (i = 0; i < ntidestat; i++) {
		xtmp = tidestat[i]->x;
		ytmp = tidestat[i]->y;
		if (inregion2(regionx, regiony, nregion, xtmp, ytmp)) {
		    g[cg].tidestat[i].wx1 = atof((char *) xv_getstr(stide_tmin_text_item));
		    g[cg].tidestat[i].wx2 = atof((char *) xv_getstr(stide_tmax_text_item));
		    g[cg].tidestat[i].wy1 = atof((char *) xv_getstr(stide_cmin_text_item));
		    g[cg].tidestat[i].wy2 = atof((char *) xv_getstr(stide_cmax_text_item));

		    g[cg].tidestat[i].display = XmToggleButtonGetState(tide_toggle_item) ? ON : OFF;
		    g[cg].tidestat[i].display_marker = XmToggleButtonGetState(tide_togglemarker_item) ? ON : OFF;
		    g[cg].tidestat[i].display_ampphase = GetChoice(tide_displayamppha_item);
		    g[cg].tidestat[i].use_freq = GetChoice(tide_usefreq_item);
		    switch (g[cg].tidestat[i].use_freq) {
		    case 0:
/* All frequencies */
			break;
		    case 1:
/* same as simulation */
			break;
		    case 2:
/* here, allocate memory for g[cg].tidestat[i].freqs */
/* count the number of selected frequencies */
			break;
		    }
		}
	    }
	}
	return;
	break;
    }
    for (i = imin; i <= imax; i++) {
	g[cg].tidestat[i].wx1 = atof((char *) xv_getstr(stide_tmin_text_item));
	g[cg].tidestat[i].wx2 = atof((char *) xv_getstr(stide_tmax_text_item));
	g[cg].tidestat[i].wy1 = atof((char *) xv_getstr(stide_cmin_text_item));
	g[cg].tidestat[i].wy2 = atof((char *) xv_getstr(stide_cmax_text_item));

	g[cg].tidestat[i].display = XmToggleButtonGetState(tide_toggle_item) ? ON : OFF;
	g[cg].tidestat[i].display_marker = XmToggleButtonGetState(tide_togglemarker_item) ? ON : OFF;
	g[cg].tidestat[i].display_ampphase = GetChoice(tide_displayamppha_item);
	g[cg].tidestat[i].use_freq = GetChoice(tide_usefreq_item);
	switch (g[cg].tidestat[i].use_freq) {
	case 0:
/* All frequencies */
	    break;
	case 1:
/* same as simulation */
	    break;
	case 2:
/* here, allocate memory for g[cg].tidestat[i].freqs */
/* count the number of selected frequencies */
	    break;
	}
    }
}

static Widget seltide_frame, seltide_panel;
static Widget *seltide_display_applyto;
static Widget teanltoggle[MAXTEANL];
static Widget adcirctoggle[MAXADCIRC];
static Widget adcircflowtoggle[MAXADCIRC];

void seltide_accept_proc(void)
{
    int h, i, j, imin, imax;
    switch (GetChoice(seltide_display_applyto)) {
    case 0:
	imin = imax = curtide;
	break;
    case 1:
	imin = 0;
	imax = ntidestat - 1;
	break;
    case 2:
	if (nregion) {
	    double xtmp, ytmp;
	    for (i = 0; i < ntidestat; i++) {
		xtmp = tidestat[i]->x;
		ytmp = tidestat[i]->y;
		if (inregion2(regionx, regiony, nregion, xtmp, ytmp)) {
		    for (j = 0; j < MAXTEANL; j++) {
			g[cg].tidestat[i].teanl[j] = XmToggleButtonGetState(teanltoggle[j]);
			g[cg].tidestat[i].tp[j].color = j + 2;
		    }
		    for (j = 0; j < MAXADCIRC; j++) {
			g[cg].tidestat[i].adcirc[j] = XmToggleButtonGetState(adcirctoggle[j]);
			g[cg].tidestat[i].ap[j].color = j + 2;
		    }
		    for (j = 0; j < MAXADCIRC; j++) {
			g[cg].tidestat[i].adcircflow[j] = XmToggleButtonGetState(adcircflowtoggle[j]);
		    }
		}
	    }
	}
	return;
	break;
    }

    for (j = imin; j <= imax; j++) {
	for (i = 0; i < MAXTEANL; i++) {
	    g[cg].tidestat[j].teanl[i] = XmToggleButtonGetState(teanltoggle[i]);
	    g[cg].tidestat[j].tp[i].color = i + 2;
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    g[cg].tidestat[j].adcirc[i] = XmToggleButtonGetState(adcirctoggle[i]);
	    g[cg].tidestat[j].ap[i].color = i + 2;
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    g[cg].tidestat[j].adcircflow[i] = XmToggleButtonGetState(adcircflowtoggle[i]);
	}
    }
}

void update_seltide(void)
{
    int i;
    if (seltide_frame) {
	for (i = 0; i < MAXTEANL; i++) {
	    XmToggleButtonSetState(teanltoggle[i], g[cg].tidestat[curtide].teanl[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcirctoggle[i], g[cg].tidestat[curtide].adcirc[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcircflowtoggle[i], g[cg].tidestat[curtide].adcircflow[i], False);
	}
    }
}

void create_seltide_frame(void)
{
    Widget wbut, fr, rc3, sep, rc, rc2;
    int i;
    setistop();
    if (!seltide_frame) {
	seltide_frame = XmCreateDialogShell(app_shell, "Tide station vs model", NULL, 0);
	seltide_panel = XmCreateRowColumn(seltide_frame, "seltiderc", NULL, 0);


	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, seltide_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("Tide station:", xmLabelGadgetClass, rc3, NULL);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, seltide_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("TEANL elevations:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXTEANL; i++) {
	    sprintf(buf, "%1d", i + 1);
	    teanltoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
/*
	teanl_color = CreateColorChoice(rc3, "Color:", 1);
*/

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, seltide_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ADCIRC elevations:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcirctoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
/*
	adcirc_color = CreateColorChoice(rc3, "Color:", 1);
*/

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, seltide_panel, NULL);
	rc3 = XtVaCreateManagedWidget("rc3", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);

	XtVaCreateManagedWidget("ADCIRC magnitudes of velocity:", xmLabelGadgetClass, rc3, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcircflowtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, rc3, NULL);
	}
/*
	adcircflow_color = CreateColorChoice(rc3, "Color:", 1);
*/
	rc = XmCreateRowColumn(seltide_panel, "rc", NULL, 0);
	seltide_display_applyto = CreatePanelChoice1(rc, "Apply to: ", 4, "Current selection", "All", "All in region", 0, 0);
	XtManageChild(rc);

	rc = XmCreateRowColumn(seltide_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) seltide_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, seltide_frame);
	XtManageChild(rc);
	XtManageChild(rc3);
	XtManageChild(seltide_panel);
    }
    XtRaise(seltide_frame);
    update_seltide();
}

static Widget tideprops_frame;

static void tideprops_accept_proc(void)
{
    int h;
    int i, j, imin, imax;
    double x, y, locx, locy;
    x = y = locx = locy = 0.0;
    switch (GetChoice(seltide_display_applyto)) {
    case 0:
	imin = imax = curtide;
	break;
    case 1:
	imin = 0;
	imax = ntidestat - 1;
	break;
    case 2:
	if (nregion) {
	    double xtmp, ytmp;
	    for (i = 0; i < ntidestat; i++) {
		xtmp = tidestat[i]->x;
		ytmp = tidestat[i]->y;
		if (inregion2(regionx, regiony, nregion, xtmp, ytmp)) {
		    g[cg].tidestat[i].p.color = GetChoice(stide_color_item);
		    g[cg].tidestat[i].p.fillcol = GetChoice(stide_fill_item);
		    g[cg].tidestat[i].p.linew = GetChoice(stide_linew_item);
		    g[cg].tidestat[i].precx = GetChoice(stide_precx_item) - 1;
		    g[cg].tidestat[i].precy = GetChoice(stide_precy_item) - 1;
		    g[cg].tidestat[i].attach = GetChoice(stide_attach_item);
		    g[cg].tidestat[i].vx = atof((char *) xv_getstr(stide_vx_text_item));
		    g[cg].tidestat[i].vy = atof((char *) xv_getstr(stide_vy_text_item));
		}
	    }
	}
	return;
	break;
    }

    for (i = imin; i <= imax; i++) {
	g[cg].tidestat[i].p.color = GetChoice(stide_color_item);
	g[cg].tidestat[i].p.fillcol = GetChoice(stide_fill_item);
	g[cg].tidestat[i].p.linew = GetChoice(stide_linew_item);
	g[cg].tidestat[i].precx = GetChoice(stide_precx_item) - 1;
	g[cg].tidestat[i].precy = GetChoice(stide_precy_item) - 1;
	g[cg].tidestat[i].attach = GetChoice(stide_attach_item);
	g[cg].tidestat[i].vx = atof((char *) xv_getstr(stide_vx_text_item));
	g[cg].tidestat[i].vy = atof((char *) xv_getstr(stide_vy_text_item));
    }
}

void update_tideprops(void)
{
    char buf[256];
    if (tideprops_frame) {
	SetChoice(stide_color_item, g[cg].tidestat[curtide].p.color);
	SetChoice(stide_fill_item, g[cg].tidestat[curtide].p.fillcol);
	SetChoice(stide_linew_item, g[cg].tidestat[curtide].p.linew);
	SetChoice(stide_precx_item, g[cg].tidestat[curtide].precx + 1);
	SetChoice(stide_precy_item, g[cg].tidestat[curtide].precy + 1);
	SetChoice(stide_attach_item, g[cg].tidestat[curtide].attach);
	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].vx);
	XmTextSetString(stide_vx_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].tidestat[curtide].vy);
	XmTextSetString(stide_vy_text_item, buf);
    }
}

void create_tideprops_frame(void)
{
    Widget bb, wbut, tideprops_panel, fr, rc3, sep, rc, rc2;
    int i;
    setistop();
    if (!tideprops_frame) {
	tideprops_frame = XmCreateDialogShell(app_shell, "Time tideory selection", NULL, 0);
	tideprops_panel = XmCreateRowColumn(tideprops_frame, "tidepropsrc", NULL, 0);
	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, tideprops_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, NULL);
	stide_color_item = CreateColorChoice(bb, "Frame color:", 1);
	stide_linew_item = CreatePanelChoice2(bb, "Frame line width: ", 2, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

	stide_fill_item = CreateColorChoice(bb, "Fill color:", 1);
	stide_precx_item = CreatePanelChoice2(bb, "Precision X labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	stide_precy_item = CreatePanelChoice2(bb, "Precision Y labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	stide_attach_item = CreatePanelChoice2(bb, "Attach to: ", 1, 5, "SW", "SE", "NW", "NE", 0, 0);

	stide_vx_text_item = CreateTextItem2(bb, 15, "View DX: ");
	stide_vy_text_item = CreateTextItem2(bb, 15, "View DY: ");
	stide_props_applyto = CreatePanelChoice1(bb, "Apply to: ", 4, "Current selection", "All", "All in region", 0, 0);

	rc = XmCreateRowColumn(tideprops_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) tideprops_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, tideprops_frame);
	XtManageChild(rc);
	XtManageChild(tideprops_panel);
    }
    XtRaise(tideprops_frame);
    update_tideprops();
}


static void set_current_frequency(int ind)
{
    if (ind >= 0) {
	curfreq = ind;
	XmListSetPos(freq_list_item, curfreq + 1);
	XmListSelectPos(freq_list_item, curfreq + 1, False);
    }
}

static void set_current_amppha(int ind)
{
    if (ind >= 0) {
	curamppha = ind;
	XmListSetPos(amppha_list_item, curamppha + 1);
	XmListSelectPos(amppha_list_item, curamppha + 1, False);
    }
}

static void set_curfreq(int w, int cd, XtPointer cb)
{
    int pos;
    char buf[256], name[512];
    XmString xms;
    XmString *s, cs;
    int *pos_list;
    int i, j, pos_cnt, cnt;
    char *cstr;
    XmListCallbackStruct *cbs = (XmListCallbackStruct *) cb;

    curfreq = cbs->item_position - 1;

    if (pos = XmListGetSelectedPos(freq_list_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(freq_list_item, XmNselectedItemCount, &cnt, XmNselectedItems, &s, NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(name, cstr);
	    XtFree(cstr);
	    XmStringFree(cs);
	}
    }
}

static void set_curamppha(int w, int cd, XtPointer cb)
{
    int pos;
    char buf[256], name[512];
    XmString xms;
    XmString *s, cs;
    int *pos_list;
    int i, j, pos_cnt, cnt;
    char *cstr;
    XmListCallbackStruct *cbs = (XmListCallbackStruct *) cb;

    curamppha = cbs->item_position - 1;

    if (pos = XmListGetSelectedPos(amppha_list_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(amppha_list_item, XmNselectedItemCount, &cnt, XmNselectedItems, &s, NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(name, cstr);
	    XtFree(cstr);
	    XmStringFree(cs);
	}
    }
}
