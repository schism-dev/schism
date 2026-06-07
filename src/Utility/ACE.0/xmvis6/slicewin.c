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
 *
 * slice Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: slicewin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

extern Widget app_shell;

static Widget slice_frame;
static Widget slice_panel;
static Widget slicesave_frame;
static Widget slicesave_panel;

static Widget sslice_frame;
static Widget sslice_panel;
void create_sslice_frame(void);
void update_slice(void);
void do_sliceline(void);
void do_slicepoly(void);
void do_slicebox(void);

/*
 * Panel item declarations
 */
static Widget *sslice_which_item;
static Widget sslice_tmin_text_item;
static Widget sslice_tmax_text_item;
static Widget sslice_cmin_text_item;
static Widget sslice_cmax_text_item;
static Widget sslice_npts_text_item;
static Widget *sslice_precx_item;
static Widget *sslice_precy_item;
static Widget *sslice_color_item;
static Widget *sslice_fill_item;
static Widget *sslice_linew_item;
static Widget *sslice_attach_item;
static Widget slice_fluxtoggle_item;
static Widget slice_displaypts_item;

static Widget slice_files_list_item;
static Widget *slice_toggle_item;
static Widget slice_markertoggle_item;
static Widget slice_teanl_item;
static Widget *slice_teanlcolor_item;
static Widget slice_adcirc_item;
static Widget *slice_adcirccolor_item;
static Widget slice_adc3d_item;
static Widget *slice_adc3dccolor_item;
static Widget slice_file_item;
static Widget *slice_filecolor_item;

static Widget slice_type_item;

/*
 * Event and Notify proc declarations
 */

void do_viewslice(void);
void slice_done_proc(void);
void slice_accept_proc(void);
void slicesave_done_proc(void);
void slicesave_accept_proc(void);
void sslice_doelev_proc();
void sslice_start_proc();
void sslice_step_proc();
void sslice_nsteps_proc();
void sslice_accept_proc(void);
void sslice_done_proc(void);

void create_selslice_frame(void);
void update_slice(void);
void update_selslice(void);

void set_sliceloc(void);

void define_sliceline(int gno, int sno, double wx1, double wy1, double wx2, double wy2)
{
    int i;
    g[gno].sbox[sno].active = ON;
    g[gno].sbox[sno].x1 = wx1;
    g[gno].sbox[sno].y1 = wy1;
    g[gno].sbox[sno].x2 = wx2;
    g[gno].sbox[sno].y2 = wy2;
    /* g[gno].sbox[sno].npts = 100; */
    g[gno].sbox[sno].type = LINE;
    for (i = 0; i < 250; i++) {
	g[gno].sbox[sno].elist[i] = -1;
    }
}

void define_slicepoly(int gno, int sno, int n, double *wx, double *wy)
{
    int i;
    if (g[gno].sbox[sno].sx != NULL) {
	free(g[gno].sbox[sno].sx);
	free(g[gno].sbox[sno].sy);
    }
    g[gno].sbox[sno].type = POLY;
    g[gno].sbox[sno].active = ON;
    g[gno].sbox[sno].sx = (double *) calloc(n, sizeof(double));
    g[gno].sbox[sno].sy = (double *) calloc(n, sizeof(double));
    g[gno].sbox[sno].npts = n;
    for (i = 0; i < n; i++) {
	g[gno].sbox[sno].sx[i] = wx[i];
	g[gno].sbox[sno].sy[i] = wy[i];
    }
    for (i = 0; i < 250; i++) {
	g[gno].sbox[sno].elist[i] = -1;
    }
}

void define_slicebox(int gno, int sno, double wx, double wy)
{
    g[gno].sbox[sno].locx = wx;
    g[gno].sbox[sno].locy = wy;
}

void set_slicelocCB(void)
{
    curslice = GetChoice(sslice_which_item);
    set_sliceloc();
}

void set_curslice(void)
{
    curslice = GetChoice(sslice_which_item);
    update_slice();
    update_selslice();
}

void toggle_slice(void)
{
    switch (GetChoice(slice_toggle_item)) {
    case 0:
	g[cg].sbox[curslice].display_slice = OFF;
	break;
    case 1:
	g[cg].sbox[curslice].display_slice = ON;
	break;
    case 2:
	g[cg].sbox[curslice].display_slice = POINT;
	break;
    case 3:
	g[cg].sbox[curslice].display_slice = BOTH;
	break;
    }
    g[cg].sbox[curslice].display_marker = XmToggleButtonGetState(slice_markertoggle_item);
}

/*
 * Create the slice Frame and the slice Panel
 */

void create_slice_frame(void)
{
    Arg wargs[8];
    Widget wbut;
    setistop();
    if (!slice_frame) {
	slice_frame = XmCreateFileSelectionDialog(app_shell, "Read slice", NULL, 0);
	XtAddCallback(slice_frame, XmNcancelCallback, (XtCallbackProc) slice_done_proc, NULL);
	XtAddCallback(slice_frame, XmNokCallback, (XtCallbackProc) slice_accept_proc, NULL);
    }
    XtRaise(slice_frame);
}

void slice_done_proc(void)
{
    XtUnmanageChild(slice_frame);
}

static char slice_fname[128];

void slice_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(slice_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(slice_fname, s);
}

void create_slicesave_frame(void)
{
    setistop();
    if (!slicesave_frame) {
	slicesave_frame = XmCreateFileSelectionDialog(app_shell, "Save slice", NULL, 0);
	XtAddCallback(slicesave_frame, XmNcancelCallback, (XtCallbackProc) slicesave_done_proc, NULL);
	XtAddCallback(slicesave_frame, XmNokCallback, (XtCallbackProc) slicesave_accept_proc, NULL);
    }
    XtRaise(slicesave_frame);
}

void slicesave_done_proc(void)
{
    XtUnmanageChild(slicesave_frame);
}

static char slicesave_fname[256];

void slicesave_accept_proc(void)
{
    int i, j, nsets, itmp, cnt = 0, magic = 20;
    char combuf[256];
    FILE *fp;
    double t;

    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(slicesave_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(slicesave_fname, s);
    fp = fopen(slicesave_fname, "w");
    if (fp != NULL) {
	fwrite(&magic, sizeof(int), 1, fp);
	fwrite(&timeclock.nsteps, sizeof(int), 1, fp);
	nsets = nslices(cg, 0);
	fwrite(&nsets, sizeof(int), 1, fp);
	for (j = 0; j < timeclock.nsteps; j++) {

	    t = timeclock.t[j];

	    if (j > 0) {
	    } else {
	    }

	    if (itmp = write_slice_binary(fp, cg, 0, j, t, cnt)) {
		cnt += itmp;
	    }
/*
	    for (i = 0; i < MAXSLICEBOXES; i++) {
		if (itmp = write_slice_binary(fp, cg, i, j, t, cnt)) {
		    cnt += itmp;
		}
	    }
*/
	}
	fclose(fp);
/*
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 %s &", slicesave_fname);
	system(combuf);
*/
    } else {
	errwin("Unable to open file for slice view");
    }
}

void create_sslice_frame(void)
{
    Widget wbut, fr, bb, sep, rc, rc2;
    int i;
    setistop();
    if (!sslice_frame) {
	sslice_frame = XmCreateDialogShell(app_shell, "Slice setup", NULL, 0);
	sslice_panel = XmCreateRowColumn(sslice_frame, "sslicerc", NULL, 0);
	sslice_which_item = CreatePanelChoice1(sslice_panel, "Use slice: ", 11, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 0, 0);
	for (i = 0; i < MAXSLICES; i++) {
	    XtAddCallback(sslice_which_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curslice, (XtPointer) i);
	}
	slice_markertoggle_item = XtVaCreateManagedWidget("Display slice marker", xmToggleButtonWidgetClass, sslice_panel, NULL);
	slice_toggle_item = CreatePanelChoice1(sslice_panel, "Slice", 5, "Not displayed", "Display as line", "Display as points", "Display as points and line", NULL, NULL);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);

	sslice_tmin_text_item = CreateTextItem2(bb, 15, "X min: ");
	sslice_tmax_text_item = CreateTextItem2(bb, 15, "X max: ");
	sslice_cmin_text_item = CreateTextItem2(bb, 15, "Y min: ");
	sslice_cmax_text_item = CreateTextItem2(bb, 15, "Y max: ");

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmVERTICAL, NULL);
	sslice_color_item = CreateColorChoice(bb, "Frame color:", 1);
	sslice_linew_item = CreatePanelChoice2(bb, "Frame line width: ", 3, 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);

	sslice_fill_item = CreateColorChoice(bb, "Fill color:", 1);
	sslice_precx_item = CreatePanelChoice2(bb, "Precision X labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	sslice_precy_item = CreatePanelChoice2(bb, "Precision Y labels: ", 4, 12, "No labels", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, 0);
	sslice_attach_item = CreatePanelChoice2(bb, "Attach to: ", 1, 5, "SW", "SE", "NW", "NE", 0, 0);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sslice_panel, NULL);
	rc2 = XmCreateRowColumn(fr, "rc2", NULL, 0);
	sslice_npts_text_item = CreateTextItem2(rc2, 7, "Number of points: ");
	rc = XmCreateRowColumn(rc2, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place slice marker", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_slicebox, NULL);
	wbut = XtVaCreateManagedWidget("Slice by line", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_sliceline, NULL);
	wbut = XtVaCreateManagedWidget("*Slice by polyline", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_slicepoly, NULL);
	XtManageChild(rc);
	XtManageChild(rc2);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, sslice_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Apply slice to...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_selslice_frame, NULL);
	wbut = XtVaCreateManagedWidget("*Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_slice_frame, NULL);
	wbut = XtVaCreateManagedWidget("*Save...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_slicesave_frame, NULL);
#ifndef CORPS
	wbut = XtVaCreateManagedWidget("View...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_viewslice, NULL);
#endif
	XtManageChild(rc);

	rc = XmCreateRowColumn(sslice_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sslice_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sslice_done_proc, NULL);
	XtManageChild(rc);
	XtManageChild(sslice_panel);
    }
    update_slice();
    update_selslice();
    XtRaise(sslice_frame);
}

void update_slice(void)
{
    char buf[256];
    if (sslice_frame) {
	SetChoice(sslice_which_item, curslice);

	XmToggleButtonSetState(slice_markertoggle_item, g[cg].sbox[curslice].display_marker == ON, False);
	switch (g[cg].sbox[curslice].display_slice) {
	case OFF:
	    SetChoice(slice_toggle_item, 0);
	    break;
	case 1:
	    SetChoice(slice_toggle_item, 1);
	    break;
	case 2:
	    SetChoice(slice_toggle_item, 2);
	    break;
	case 3:
	    SetChoice(slice_toggle_item, 3);
	    break;
	}
	sprintf(buf, "%.2lf", g[cg].sbox[curslice].wx1);
	XmTextSetString(sslice_tmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].sbox[curslice].wx2);
	XmTextSetString(sslice_tmax_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].sbox[curslice].wy1);
	XmTextSetString(sslice_cmin_text_item, buf);
	sprintf(buf, "%.2lf", g[cg].sbox[curslice].wy2);
	XmTextSetString(sslice_cmax_text_item, buf);
	sprintf(buf, "%d", g[cg].sbox[curslice].npts);
	XmTextSetString(sslice_npts_text_item, buf);

/*
	XmToggleButtonSetState(slice_teanl_item, sbox[curslice].models[0], False);
	SetChoice(slice_teanlcolor_item, sbox[curslice].set_colors[0]);
	XmToggleButtonSetState(slice_adcirc_item, sbox[curslice].models[1], False);
	SetChoice(slice_adcirccolor_item, sbox[curslice].set_colors[1]);
	XmToggleButtonSetState(slice_file_item, sbox[curslice].models[2], False);
	SetChoice(slice_filecolor_item, sbox[curslice].set_colors[2]);
*/

	SetChoice(sslice_color_item, g[cg].sbox[curslice].p.color);
	SetChoice(sslice_fill_item, g[cg].sbox[curslice].p.fillcol);
	SetChoice(sslice_linew_item, g[cg].sbox[curslice].p.linew);
	SetChoice(sslice_precx_item, g[cg].sbox[curslice].precx + 1);
	SetChoice(sslice_precy_item, g[cg].sbox[curslice].precy + 1);
	SetChoice(sslice_attach_item, g[cg].sbox[curslice].attach);
    }
}

void sslice_done_proc(void)
{
    XtUnmanageChild(sslice_frame);
}

void do_viewslice(void)
{
    int i, itmp, cnt = 0;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEsliceXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	fprintf(fp, "@title \"Slice\"\n");
	fprintf(fp, "@xaxis label \"Distance\"\n");
	for (i = 0; i < MAXSLICEBOXES; i++) {
	    if (itmp = viewslice(fp, cg, i, c, t, cnt)) {
		cnt += itmp;
	    }
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    } else {
	errwin("Unable to open file for slice view");
    }
}

load_teanl_slice(void)
{
}

load_adcirc_slice(void)
{
}

void sslice_accept_proc(void)
{
    int i, h;
    double x, y, locx, locy;
    x = y = locx = locy = 0.0;

    curslice = h = GetChoice(sslice_which_item);

    g[cg].sbox[curslice].wx1 = atof((char *) xv_getstr(sslice_tmin_text_item));
    g[cg].sbox[curslice].wx2 = atof((char *) xv_getstr(sslice_tmax_text_item));
    g[cg].sbox[curslice].wy1 = atof((char *) xv_getstr(sslice_cmin_text_item));
    g[cg].sbox[curslice].wy2 = atof((char *) xv_getstr(sslice_cmax_text_item));
    g[cg].sbox[curslice].npts = atoi((char *) xv_getstr(sslice_npts_text_item));

    g[cg].sbox[curslice].display_marker = XmToggleButtonGetState(slice_markertoggle_item) ? ON : OFF;
    switch (GetChoice(slice_toggle_item)) {
    case 0:
	g[cg].sbox[curslice].display_slice = OFF;
	break;
    case 1:
	g[cg].sbox[curslice].display_slice = ON;
	break;
    case 2:
	g[cg].sbox[curslice].display_slice = POINT;
	break;
    case 3:
	g[cg].sbox[curslice].display_slice = BOTH;
	break;
    }

/*
    sbox[curslice].models[0] = XmToggleButtonGetState(slice_teanl_item);
    sbox[curslice].set_colors[0] = GetChoice(slice_teanlcolor_item);
    sbox[curslice].models[1] = XmToggleButtonGetState(slice_adcirc_item);
    sbox[curslice].set_colors[1] = GetChoice(slice_adcirccolor_item);
    sbox[curslice].models[2] = XmToggleButtonGetState(slice_file_item);
    sbox[curslice].set_colors[2] = GetChoice(slice_filecolor_item);
*/

    g[cg].sbox[curslice].p.color = GetChoice(sslice_color_item);
    g[cg].sbox[curslice].p.fillcol = GetChoice(sslice_fill_item);
    g[cg].sbox[curslice].p.linew = GetChoice(sslice_linew_item);
    g[cg].sbox[curslice].precx = GetChoice(sslice_precx_item) - 1;
    g[cg].sbox[curslice].precy = GetChoice(sslice_precy_item) - 1;
    g[cg].sbox[curslice].attach = GetChoice(sslice_attach_item);
    for (i = 0; i < 250; i++) {
	g[cg].sbox[curslice].elist[i] = -1;
    }

/*
    if (h >= 0) {
	if (sbox[curslice].models[2]) {
	    set_wait_cursor(sslice_frame);
	    if (!read_slice(curslice, slice_fname)) {
		sprintf(buf, "Error reading file %s", slice_fname);
		errwin(buf);
		sbox[curslice].models[2] = 0;
	    } else {
	    }
	    unset_wait_cursor(sslice_frame);
	}
	if (sbox[curslice].models[1]) {
	    load_adcirc_slice(curslice, 0);
	}
	if (sbox[curslice].models[0]) {
	    load_teanl_slice(curslice, 0);
	}
    }
*/
}

static Widget selslice_frame, selslice_panel;
static Widget *teanl_color;
static Widget *teanlmag_color;
static Widget *adcircmag_color;
static Widget *bath_color;
static Widget *ela_color;
static Widget *adcirc_color;
static Widget *adc3d_color;
static Widget bathtoggle[MAXGRIDS];
static Widget teanltoggle[MAXTEANL];
static Widget teanlmagtoggle[MAXTEANL];
static Widget adcirctoggle[MAXADCIRC];
static Widget adcircmagtoggle[MAXADCIRC];
static Widget adc3dtoggle[MAXADCIRC];
static Widget elatoggle[MAXELA];

void selslice_accept_proc(void)
{
    int i;
    for (i = 0; i < MAXGRIDS; i++) {
	g[cg].sbox[curslice].bath[i] = XmToggleButtonGetState(bathtoggle[i]);
    }
    for (i = 0; i < MAXTEANL; i++) {
	g[cg].sbox[curslice].teanl[i] = XmToggleButtonGetState(teanltoggle[i]);
    }
    for (i = 0; i < MAXTEANL; i++) {
	g[cg].sbox[curslice].teanlmag[i] = XmToggleButtonGetState(teanlmagtoggle[i]);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].sbox[curslice].adcirc[i] = XmToggleButtonGetState(adcirctoggle[i]);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	g[cg].sbox[curslice].adcircmag[i] = XmToggleButtonGetState(adcircmagtoggle[i]);
    }
    for (i = 0; i < MAXADCIRC3D; i++) {
	g[cg].sbox[curslice].adc3d[i] = XmToggleButtonGetState(adc3dtoggle[i]);
    }
    for (i = 0; i < MAXELA; i++) {
	g[cg].sbox[curslice].ela[i] = XmToggleButtonGetState(elatoggle[i]);
    }
}

void update_selslice(void)
{
    int i;
    if (selslice_frame) {
	for (i = 0; i < MAXGRIDS; i++) {
	    XmToggleButtonSetState(bathtoggle[i], g[cg].sbox[curslice].bath[i], False);
	}
	for (i = 0; i < MAXTEANL; i++) {
	    XmToggleButtonSetState(teanltoggle[i], g[cg].sbox[curslice].teanl[i], False);
	}
	for (i = 0; i < MAXTEANL; i++) {
	    XmToggleButtonSetState(teanlmagtoggle[i], g[cg].sbox[curslice].teanlmag[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcirctoggle[i], g[cg].sbox[curslice].adcirc[i], False);
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    XmToggleButtonSetState(adcircmagtoggle[i], g[cg].sbox[curslice].adcircmag[i], False);
	}
	for (i = 0; i < MAXADCIRC3D; i++) {
	    XmToggleButtonSetState(adc3dtoggle[i], g[cg].sbox[curslice].adc3d[i], False);
	}
	for (i = 0; i < MAXELA; i++) {
	    XmToggleButtonSetState(elatoggle[i], g[cg].sbox[curslice].ela[i], False);
	}
    }
}

void create_selslice_frame(void)
{
    Widget wbut, fr, bb, sep, rc, rc2;
    int i;
    setistop();
    if (!selslice_frame) {
	selslice_frame = XmCreateDialogShell(app_shell, "Slice setup", NULL, 0);
	selslice_panel = XmCreateRowColumn(selslice_frame, "selslicerc", NULL, 0);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("Grid bathymetry:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXGRIDS; i++) {
	    sprintf(buf, "%1d", i + 1);
	    bathtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	bath_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("TEANL elevations:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXTEANL; i++) {
	    sprintf(buf, "%1d", i + 1);
	    teanltoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	teanl_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("TEANL magnitudes of velocity:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXTEANL; i++) {
	    sprintf(buf, "%1d", i + 1);
	    teanlmagtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	teanlmag_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ADCIRC elevations:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcirctoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	adcirc_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ADCIRC magnitudes of velocity:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXADCIRC; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adcircmagtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	adcircmag_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ADCIRC3D elevations:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXADCIRC3D; i++) {
	    sprintf(buf, "%1d", i + 1);
	    adc3dtoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	adc3d_color = CreateColorChoice(bb, "Color:", 1);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, selslice_panel, NULL);
	bb = XtVaCreateManagedWidget("bb", xmRowColumnWidgetClass, fr, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("ELA concentrations:", xmLabelGadgetClass, bb, NULL);
	for (i = 0; i < MAXELA; i++) {
	    sprintf(buf, "%1d", i + 1);
	    elatoggle[i] = XtVaCreateManagedWidget(buf, xmToggleButtonWidgetClass, bb, NULL);
	}
	ela_color = CreateColorChoice(bb, "Color:", 1);

	rc = XmCreateRowColumn(selslice_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) selslice_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, selslice_frame);
	XtManageChild(rc);
	XtManageChild(bb);
	XtManageChild(selslice_panel);
    }
    XtRaise(selslice_frame);
    update_selslice();
}
