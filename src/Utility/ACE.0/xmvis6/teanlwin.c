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
 * Read TEA-NL files popup and set TEA-NL parameters popup
 */

#ifndef lint
static char RCSid[] = "$Id: teanlwin.c,v 1.3 2004/02/26 19:32:49 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

static Widget teanl_frame;
static Widget teanl_panel;

static Widget steanl_frame;
static Widget steanl_panel;
void create_steanl_frame(void);
void create_teanlprops_frame(void);
void create_isolines_popup(char *title, Isolparms * cd, int cm);
void create_inundate_frame(void);
void update_teanl_flow(void);
static void teanl_view_proc(Widget w, int cd);

/*
 * Panel item declarations
 */
static Widget steanl_start_text_item;
static Widget steanl_step_text_item;
static Widget steanl_nsteps_text_item;
static Widget teanlel_emin_text_item;
static Widget teanlel_emax_text_item;
static Widget *flow_toggle_item;
static Widget *flow_color_item;
static Widget flow_type_item;
static Widget elev_toggle_item;
static Widget elevm_toggle_item;
static Widget elevm_attach_item;
static Widget elevm_node_item;
static Widget flowelev_toggle_item;
static Widget flowelev_depth_item;
static Widget teanl_phase_item;
static Widget teanl_amp_item;
static Widget teanl_mag_item;
static Widget freq_list_item;
static Widget amppha_list_item;
static Widget *teanl_frequnits_item;
static Widget *teanl_flow_item;
static Widget *teanl_grid_item;
static Widget teanl_dcor_item;
static Widget teanlsign_toggle_item;
static Widget *steanl_flow_item;
static Widget *steanl_elevm_item;

/*
 * Event and Notify proc declarations
 */
void teanl_accept_proc(void);
void steanl_doelev_proc();
void steanl_start_proc();
void steanl_step_proc();
void steanl_nsteps_proc();
void steanl_accept_proc(void);
void teanlel_place_notify_proc(void);
void teanlel_edit_notify_proc(void);
void teanlel_accept_notify_proc(void);
void teanlel_clear_notify_proc();

static int curfreq = 0;
static int curamppha = 0;

/*
 * Create the teanl Frame and the teanl Panel
 */

void set_curteanl(Widget w, int cd)
{
    if (cd != curteanl) {
	curfreq = 0;
	curamppha = 0;
    }
    curteanl = cd;
    update_teanl_flow();
}

void set_curteanlem(Widget w, int cd)
{
    curteanlem = cd;
    update_teanl_flow();
}

void create_teanl_frame(void)
{
    Widget wbut, rc;
    int i;
    setistop();
    if (!teanl_frame) {
	teanl_frame = XmCreateFileSelectionDialog(app_shell, "Read TEA-NL", NULL, 0);
    	XtVaSetValues(teanl_frame,
       	XmNdialogTitle, XmStringCreate("Read TEA-NL file", charset),
       	XmNtitle, "Read TEA-NL file",
       	XmNcancelLabelString, XmStringCreate("Done", charset),
       	XmNdirMask, XmStringCreate("*.t[co][tu]", charset), 
       	NULL);

	rc = XmCreateRowColumn(teanl_frame, "rc", NULL, 0);
	teanl_flow_item = CreatePanelChoice1(rc, "Read to flow: ", 6, "1", "2", "3", "4", "5", 0, 0);
	teanl_grid_item = CreatePanelChoice1(rc, "Using reference grid: ", 6, "1", "2", "3", "4", "5", 0, 0);
	for (i = 0; i < 5; i++) {
	    XtAddCallback(teanl_flow_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curteanl, (XtPointer) i);
	}
	teanlsign_toggle_item = XtVaCreateManagedWidget("Phases are added (for .tou files)", xmToggleButtonWidgetClass, rc, NULL);
	teanl_dcor_item = CreateTextItem2(rc, 5, "Depth correction: ");
	XtManageChild(rc);
	XtAddCallback(teanl_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, teanl_frame);
	XtAddCallback(teanl_frame, XmNokCallback, (XtCallbackProc) teanl_accept_proc, NULL);
    }
    XtRaise(teanl_frame);
    update_teanl_flow();
    sprintf(buf, "%lf", flowf[curteanl].dcor);
    xv_setstr(teanl_dcor_item, buf);
    XmToggleButtonSetState(teanlsign_toggle_item, flowf[curteanl].sign_of_phase == 1 ? False : True, False);
}

void teanl_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int flowno, gridno, fbin;
    Widget textw;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(teanl_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    flowno = GetChoice(teanl_flow_item);
    gridno = GetChoice(teanl_grid_item);
    flowf[curteanl].sign_of_phase = XmToggleButtonGetState(teanlsign_toggle_item) ? -1 : 1;
    if (grid[gridno].active != ON) {
	errwin("Grid not active, read grid first");
	return;
    }
    if (flowf[flowno].active == ON) {
	if (!yesno("Flow is active, kill it ?", " ", " YES ", " NO ")) {
	    return;
	}
    }
    set_wait_cursor(teanl_frame);
    if (!readflow(flowno, gridno, s)) {
	unset_wait_cursor(teanl_frame);
	sprintf(buf, "Error reading file %s", s);
	errwin(buf);
	return;
    }
    flowf[flowno].dcor = atof((char *) xv_getstr(teanl_dcor_item));
    unset_wait_cursor(teanl_frame);
    XtUnmanageChild(teanl_frame);
    create_steanl_frame();
}

static void do_teanl_isolines(Widget w, int cd)
{
    int i;
    double cmax;
    switch (cd) {
    case 0:
	cmax = 0.0;
	for (i = 0; i < flowf[curteanl].nfreq; i++) {
	    cmax += flowf[curteanl].ampmax[i];
	}
	if (g[cg].flowf[curteanl].display_elevdepth == ON) {
	    g[cg].flowf[curteanl].elevip.cmin = -cmax + grid[flowf[curteanl].grid].dmin;
	    g[cg].flowf[curteanl].elevip.cmax = cmax + grid[flowf[curteanl].grid].dmax;
	    create_isolines_popup("TEANL elevations+depth", &g[cg].flowf[curteanl].elevip, 1);
	} else {
	    g[cg].flowf[curteanl].elevip.cmin = -cmax;
	    g[cg].flowf[curteanl].elevip.cmax = cmax;
	    create_isolines_popup("TEANL elevations", &g[cg].flowf[curteanl].elevip, 1);
	}
	break;
    case 1:
	g[cg].flowf[curteanl].ampip.cmin = flowf[curteanl].ampmin[curamppha];
	g[cg].flowf[curteanl].ampip.cmax = flowf[curteanl].ampmax[curamppha];
	create_isolines_popup("TEANL amplitudes", &g[cg].flowf[curteanl].ampip, 1);
	break;
    case 2:
	g[cg].flowf[curteanl].phaseip.cmin = flowf[curteanl].phamin[curamppha];
	g[cg].flowf[curteanl].phaseip.cmax = flowf[curteanl].phamax[curamppha];
	create_isolines_popup("TEANL phases", &g[cg].flowf[curteanl].phaseip, 1);
	break;
    case 3:
	create_isolines_popup("TEANL magnitudes of velocity", &g[cg].flowf[curteanl].magip, 1);
	break;
    }
}

void set_current_frequency(int ind)
{
    if (ind >= 0) {
	curfreq = ind;
	XmListSetPos(freq_list_item, curfreq + 1);
	XmListSelectPos(freq_list_item, curfreq + 1, False);
    }
}

void set_current_amppha(int ind)
{
    if (ind >= 0) {
	curamppha = ind;
	XmListSetPos(amppha_list_item, curamppha + 1);
	XmListSelectPos(amppha_list_item, curamppha + 1, False);
    }
}

void set_curfreq(int w, int cd, XtPointer cb)
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

void set_curamppha(int w, int cd, XtPointer cb)
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

void create_steanl_frame(void)
{
    Widget lab, wbut, fr, bb, sep, rc, rc2;
    int i, got_one = 0;
    Arg al[10];
    setistop();
    for (i = 0; i < MAXTEANL; i++) {
	if (flowf[i].active == ON) {
	    got_one = 1;
	}
    }
    if (!got_one) {
	create_teanl_frame();
	return;
    }
    if (!steanl_frame) {
	steanl_frame = XmCreateDialogShell(app_shell, "TEA-NL setup", NULL, 0);
	handle_close(steanl_frame);
	steanl_panel = XmCreateRowColumn(steanl_frame, "steanl_rc", NULL, 0);

	steanl_flow_item = CreatePanelChoice1(steanl_panel, "Apply to flow: ", 7, "1", "2", "3", "4", "5", "All", 0, 0);
	for (i = 0; i < 5; i++) {
	    XtAddCallback(steanl_flow_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curteanl, (XtPointer) i);
	}

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, steanl_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	steanl_elevm_item = CreatePanelChoice2(rc, "Elevation marker: ", 4, 21, "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 0, 0);
	for (i = 0; i < MAXELEVMARKERS; i++) {
	    XtAddCallback(steanl_elevm_item[i + 2], XmNactivateCallback, (XtCallbackProc) set_curteanlem, (XtPointer) i);
	}

/*
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2,
		XmNorientation, XmHORIZONTAL,
		NULL);
	elevm_attach_item = XtVaCreateManagedWidget("Attach to node",
					 xmToggleButtonWidgetClass, rc2,
					 NULL);
	elevm_node_item = CreateTextItem2(rc2, 5, "Node: ");
	XtManageChild(rc2);
*/

	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNpacking, XmPACK_COLUMN, XmNorientation, XmHORIZONTAL, XmNnumColumns, 2, XmNisAligned, True, XmNentryAlignment, XmALIGNMENT_BEGINNING, NULL);
	sep = XmCreateLabel(rc2, "Marker min: ", NULL, 0);
	XtManageChild(sep);
	teanlel_emin_text_item = XmCreateText(rc2, "t1", NULL, 0);
	XtVaSetValues(teanlel_emin_text_item, XmNcolumns, 7, NULL);
	XtManageChild(teanlel_emin_text_item);
	sep = XmCreateLabel(rc2, "Max: ", NULL, 0);
	XtManageChild(sep);
	teanlel_emax_text_item = XmCreateText(rc2, "t2", NULL, 0);
	XtVaSetValues(teanlel_emax_text_item, XmNcolumns, 7, NULL);
	XtManageChild(teanlel_emax_text_item);
	XtManageChild(rc2);

	elevm_toggle_item = XtVaCreateManagedWidget("Display this elevation marker", xmToggleButtonWidgetClass, rc, NULL);
	elev_toggle_item = XtVaCreateManagedWidget("Display all active elevation markers", xmToggleButtonWidgetClass, rc, NULL);
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) teanlel_place_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Edit", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) teanlel_edit_notify_proc, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) teanlel_accept_notify_proc, NULL);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, steanl_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	flow_toggle_item = CreatePanelChoice1(rc, "Flow: ", 5, "Not displayed", "Display at nodes", "Display at centers", "Display ellipse", 0, 0);

	XtSetArg(al[0], XmNlistSizePolicy, XmRESIZE_IF_POSSIBLE);
	XtSetArg(al[1], XmNvisibleItemCount, 5);
	XtSetArg(al[2], XmNselectionPolicy, XmSINGLE_SELECT);

	lab = XmCreateLabel(rc, "Use frequency:", NULL, 0);
	XtManageChild(lab);
	freq_list_item = XmCreateScrolledList(rc, "list", al, 3);
	XtAddCallback(freq_list_item, XmNsingleSelectionCallback, (XtCallbackProc) set_curfreq, 0);
	XtManageChild(freq_list_item);

	flow_color_item = CreateColorChoice(rc, "Display vectors using color:", 1);

	teanl_mag_item = XtVaCreateManagedWidget("Display isolines of magnitudes", xmToggleButtonWidgetClass, rc, NULL);
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of magnitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_teanl_isolines, (XtPointer) 3);
	wbut = XtVaCreateManagedWidget("1D View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) teanl_view_proc, (XtPointer) 1);
	XtManageChild(rc2);

	XtManageChild(rc);
	XtManageChild(fr);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, steanl_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	teanl_amp_item = XtVaCreateManagedWidget("Display isolines of amplitudes", xmToggleButtonWidgetClass, rc, NULL);
	teanl_phase_item = XtVaCreateManagedWidget("Display isolines of phases", xmToggleButtonWidgetClass, rc, NULL);

	lab = XmCreateLabel(rc, "Use frequency:", NULL, 0);
	XtManageChild(lab);
	amppha_list_item = XmCreateScrolledList(rc, "list", al, 3);
	XtAddCallback(amppha_list_item, XmNsingleSelectionCallback, (XtCallbackProc) set_curamppha, 0);
	XtManageChild(amppha_list_item);

/*
	teanl_frequnits_item = CreatePanelChoice2(rc,
					     "Units of phase: ",
					     4,
					     "Radians", "Degrees", "Minutes",
					     0, 0);
*/
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of amplitudes...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_teanl_isolines, (XtPointer) 1);
	wbut = XtVaCreateManagedWidget("Isolines of phases...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_teanl_isolines, (XtPointer) 2);
	XtManageChild(rc2);
	XtManageChild(rc);
	XtManageChild(fr);

	fr = XtVaCreateManagedWidget("frame", xmFrameWidgetClass, steanl_panel, NULL);
	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	flowelev_toggle_item = XtVaCreateManagedWidget("Display isolines of elevation", xmToggleButtonWidgetClass, rc, NULL);
	flowelev_depth_item = XtVaCreateManagedWidget("Add depth", xmToggleButtonWidgetClass, rc, NULL);
	rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Isolines of elevations...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_teanl_isolines, 0);
	wbut = XtVaCreateManagedWidget("1D View...", xmPushButtonWidgetClass, rc2, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) teanl_view_proc, 0);
	XtManageChild(rc2);
	XtManageChild(rc);
	XtManageChild(fr);

/*
	rc = XmCreateRowColumn(steanl_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("*Props...", xmPushButtonWidgetClass, rc,
			    NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_teanlprops_frame, NULL);

	XtManageChild(rc);
*/

	rc = XmCreateRowColumn(steanl_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) steanl_accept_proc, NULL);

	wbut = XtCreateManagedWidget("Files...", xmPushButtonWidgetClass, rc, NULL, 0);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_teanl_frame, NULL);

	wbut = XtCreateManagedWidget("Inundation...", xmPushButtonWidgetClass, rc, NULL, 0);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_inundate_frame, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, steanl_frame);
	XtManageChild(rc);
    }
    update_teanl_flow();
    XtManageChild(steanl_panel);
    XtRaise(steanl_frame);
}

void update_teanl_flow(void)
{
    char buf[256];
    int flowno, itmp, i;
    double period;
    XmString xms;
    if (steanl_frame) {
	SetChoice(steanl_flow_item, curteanl);
	flowno = curteanl;
	SetChoice(flow_color_item, g[cg].flowf[flowno].p.color);
	if (g[cg].flowf[flowno].flowfreq == ALL) {
	    itmp = 0;
	} else {
	    itmp = g[cg].flowf[flowno].flowfreq + 1;
	}
/*
	SetChoice(teanl_freq2_item, itmp);
	SetChoice(teanl_freq_item, g[cg].flowf[flowno].freq);
*/
	switch (g[cg].flowf[flowno].display) {
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
	XmToggleButtonSetState(flowelev_toggle_item, g[cg].flowf[flowno].display_elev == ON, False);
	XmToggleButtonSetState(flowelev_depth_item, g[cg].flowf[flowno].display_elevdepth == ON, False);
	XmToggleButtonSetState(elev_toggle_item, g[cg].flowf[flowno].display_elevmarkers == ON, False);
	XmToggleButtonSetState(teanl_phase_item, g[cg].flowf[flowno].display_phase == ON, False);
	XmToggleButtonSetState(teanl_amp_item, g[cg].flowf[flowno].display_amp == ON, False);
	XmToggleButtonSetState(teanl_mag_item, g[cg].flowf[flowno].display_mag == ON, False);
	XmToggleButtonSetState(elevm_toggle_item, g[cg].flowf[flowno].em[curteanlem].display == ON, False);
	sprintf(buf, "%.3lf", g[cg].flowf[curteanl].em[curteanlem].emin);
	xv_setstr(teanlel_emin_text_item, buf);
	sprintf(buf, "%.3lf", g[cg].flowf[curteanl].em[curteanlem].emax);
	xv_setstr(teanlel_emax_text_item, buf);
	SetChoice(steanl_elevm_item, curteanlem);

	XmListDeleteAllItems(freq_list_item);
	xms = XmStringCreateLtoR("ALL frequencies", charset);
	XmListAddItemUnselected(freq_list_item, xms, 0);
	XmStringFree(xms);
	for (i = 0; i < flowf[flowno].nfreq; i++) {
	    if (strlen(flowf[flowno].freqname[i]) > 0) {
		if (flowf[flowno].omega[i] != 0.0) {
		    period = 2.0 * M_PI / (flowf[flowno].omega[i] * 3600.0);
		} else {
		    period = 999999.0;
		}
		sprintf(buf, "%s (%lf, %lf)", flowf[flowno].freqname[i], flowf[flowno].omega[i], period);
	    } else {
		strcpy(buf, "UNKNOWN");
	    }
	    xms = XmStringCreateLtoR(buf, charset);
	    XmListAddItemUnselected(freq_list_item, xms, 0);
	    XmStringFree(xms);
	}
	XmListSelectPos(freq_list_item, curfreq + 1, False);

	XmListDeleteAllItems(amppha_list_item);
	for (i = 0; i < flowf[flowno].nfreq; i++) {
	    if (strlen(flowf[flowno].freqname[i]) > 0) {
		if (flowf[flowno].omega[i] != 0.0) {
		    period = 2.0 * M_PI / (flowf[flowno].omega[i] * 3600.0);
		} else {
		    period = 999999.0;
		}
		sprintf(buf, "%s (%lf, %lf)", flowf[flowno].freqname[i], flowf[flowno].omega[i], period);
	    } else {
		strcpy(buf, "UNKNOWN");
	    }
	    xms = XmStringCreateLtoR(buf, charset);
	    XmListAddItemUnselected(amppha_list_item, xms, 0);
	    XmStringFree(xms);
	}
	XmListSelectPos(amppha_list_item, curfreq + 1, False);
    }
}

void steanl_accept_proc(void)
{
    char buf[256];
    int i, itmp, flowno, type, start, stop;
    curteanl = flowno = GetChoice(steanl_flow_item);
    if (flowno == MAXTEANL) {
	start = 0;
	stop = MAXTEANL - 1;
    } else {
	start = stop = flowno;
    }
    type = GetChoice(flow_toggle_item);
    for (i = start; i <= stop; i++) {
	switch (type) {
	case 0:
	    g[cg].flowf[i].display = OFF;
	    break;
	case 1:
	    g[cg].flowf[i].display = NODES;
	    break;
	case 2:
	    g[cg].flowf[i].display = CENTER;
	    break;
	case 3:
	    g[cg].flowf[i].display = ELLIPSE;
	    break;
	}
	g[cg].flowf[i].em[curteanlem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
	g[cg].flowf[i].em[curteanlem].emin = atof(XmTextGetString(teanlel_emin_text_item));
	g[cg].flowf[i].em[curteanlem].emax = atof(XmTextGetString(teanlel_emax_text_item));
	g[cg].flowf[i].display_elev = XmToggleButtonGetState(flowelev_toggle_item) ? ON : OFF;
	if (g[cg].flowf[i].display_elev == ON && g[cg].flowf[i].elevip.cint == 0.0) {
	    g[cg].flowf[i].elevip.cmin = -g[cg].flowf[i].elevip.cmax;
	    default_isolines(cg, &g[cg].flowf[i].elevip);
	}
	g[cg].flowf[i].display_elevdepth = XmToggleButtonGetState(flowelev_depth_item) ? ON : OFF;
	g[cg].flowf[i].display_elevmarkers = XmToggleButtonGetState(elev_toggle_item) ? ON : OFF;
	g[cg].flowf[i].display_amp = XmToggleButtonGetState(teanl_amp_item) ? ON : OFF;
	g[cg].flowf[i].display_phase = XmToggleButtonGetState(teanl_phase_item) ? ON : OFF;
	g[cg].flowf[i].display_mag = XmToggleButtonGetState(teanl_mag_item) ? ON : OFF;
	g[cg].flowf[i].p.color = GetChoice(flow_color_item);
/*
	g[cg].flowf[i].freq = GetChoice(teanl_freq_item);
	itmp = GetChoice(teanl_freq2_item);
*/
	g[cg].flowf[i].freq = curamppha;
	itmp = curfreq;
	if (itmp == 0) {
	    g[cg].flowf[flowno].flowfreq = ALL;
	} else {
	    g[cg].flowf[flowno].flowfreq = itmp - 1;
	}
	if (g[cg].flowf[i].freq < flowf[i].nfreq) {
	} else {
	}
    }
    if (flowf[curteanl].active == ON) {
	display_image();
    }
    update_display();
}

void teanlel_place_notify_proc(void)
{
    g[cg].flowf[curteanl].em[curteanlem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
    g[cg].flowf[curteanl].em[curteanlem].emin = atof(XmTextGetString(teanlel_emin_text_item));
    g[cg].flowf[curteanl].em[curteanlem].emax = atof(XmTextGetString(teanlel_emax_text_item));
    set_teanlelev();
}

void teanlel_edit_notify_proc(void)
{
    set_action(0);
    set_action(EDIT_TEANL_ELEV);
}

void teanlel_accept_notify_proc(void)
{
    g[cg].flowf[curteanl].em[curteanlem].display = XmToggleButtonGetState(elevm_toggle_item) ? ON : OFF;
    g[cg].flowf[curteanl].em[curteanlem].emin = atof(XmTextGetString(teanlel_emin_text_item));
    g[cg].flowf[curteanl].em[curteanlem].emax = atof(XmTextGetString(teanlel_emax_text_item));
}

static void teanl_view_proc(Widget w, int cd)
{
    int i;
    double eval_teanl_node(int flowno, int node, double t);
    double get_current_time(void);
    double u, v;
    char tbuf[256];
    char combuf[256];
    char *fname;
    FILE *fp;
    int c = get_current_step();
    double t = get_current_time();
    strcpy(tbuf, "/tmp/ACEteanlelevxXXXXXX");
    mkstemp(tbuf);
    fname = tbuf;
    fp = fopen(fname, "w");
    if (fp != NULL) {
	switch (cd) {
	case 0:
	    fprintf(fp, "@title \"TEANL elevations at nodes\"\n");
	    fprintf(fp, "@subtitle \"Step %d, Time = %lf\"\n", c + 1, t);
	    fprintf(fp, "@xaxis label \"Node number\"\n");
	    fprintf(fp, "@yaxis label \"Elevation\"\n");
	    for (i = 0; i < grid[flowf[curteanl].grid].nmnp; i++) {
		fprintf(fp, "%d %lf\n", i + 1, eval_teanl_node(curteanl, i, t));
	    }
	    break;
	case 1:		/* magnitudes of velocity */
	    fprintf(fp, "@title \"TEANL magnitudes of velocity at nodes\"\n");
	    fprintf(fp, "@subtitle \"Step %d, Time = %lf\"\n", c + 1, t);
	    fprintf(fp, "@xaxis label \"Node number\"\n");
	    fprintf(fp, "@yaxis label \"Magnitude\"\n");
	    for (i = 0; i < grid[flowf[curteanl].grid].nmnp; i++) {
		eval_teanl_flownode(curteanl, i, t, &u, &v);
		fprintf(fp, "%d %lf\n", i + 1, hypot(u, v));
	    }
	    break;
	}
	fclose(fp);
	sprintf(combuf, " /usr/local/ace/bin/xmgr5 -remove %s &", fname);
	system(combuf);
    } else {
	errwin("Unable to open file");
    }
}

void create_teanlprops_frame(void)
{
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

void inundate_accept_proc()
{
    g[cg].flowf[curteanl].wet.color = GetChoice(wet_color_item);
    g[cg].flowf[curteanl].dry.color = GetChoice(dry_color_item);
    g[cg].flowf[curteanl].wetdry.color = GetChoice(wetdry_color_item);
    g[cg].flowf[curteanl].display_inun = XmToggleButtonGetState(inundate_toggle_item) ? ON : OFF;
    g[cg].flowf[curteanl].display_irestrict = XmToggleButtonGetState(restrict_toggle_item) ? ON : OFF;
}

void update_inundate_flow()
{
    if (inundate_frame) {
	SetChoice(wet_color_item, g[cg].flowf[curteanl].wet.color);
	SetChoice(dry_color_item, g[cg].flowf[curteanl].dry.color);
	SetChoice(wetdry_color_item, g[cg].flowf[curteanl].wetdry.color);
	XmToggleButtonSetState(inundate_toggle_item, g[cg].flowf[curteanl].display_inun == ON, False);
	XmToggleButtonSetState(restrict_toggle_item, g[cg].flowf[curteanl].display_irestrict == ON, False);
    }
}

void create_inundate_frame(void)
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
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) inundate_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, inundate_frame);
	XtManageChild(rc);
	XtManageChild(inundate_panel);
    }
    update_inundate_flow();
    XtRaise(inundate_frame);
}
