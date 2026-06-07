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
 * status popup
 *
 */

#ifndef lint
static char RCSid[] = "$Id: statuswin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

#define SPAGESIZE 10
#define NPAGES (MAXPLOT/SPAGESIZE)

#define MAXITEMS 20

static int npages = NPAGES;

#define getdx(gno, setn)    getcol(gno, setn, 2)
#define getdy(gno, setn)    getcol(gno, setn, 3)

Widget status_frame = (Widget) 0;
Widget status_panel;
Widget header_item1;
Widget *select_status_item;
static int curpage = 0;
static int show_dxdy = 0;
extern XmStringCharSet charset;

static Widget header_w;
static Widget lab[MAXITEMS];

static char header[256];

static XFontStruct *f;
XmFontList xmf;

void status_item_proc(void);
static void update_all_status(void);

void update_set_status(void)
{
}

static void update_teanl_status(int flowno);
static void update_adcirc_status(int flowno);
static void update_adcirc3d_status(int flowno);
static void update_ela_status(int concno);
static void update_grid_status(int gridno);

static int cur_statusitem = TEANL;

static void status_done_proc(void)
{
    XtUnmanageChild(status_frame);
}

static void set_status_label(Widget w, char *buf)
{
    Arg al;
    XmString ls;
    ls = XmStringCreateLtoR(buf, charset);
    XtSetArg(al, XmNlabelString, ls);
    XtSetValues(w, &al, 1);
    XmStringFree(ls);
}

void update_status(int gno, int itemtype, int itemno)
{
    int i;
    void update_graph_status(int gno);

    if (status_frame) {
	set_status_label(header_w, header);
	switch (itemtype) {
	case ALL:
	    update_all_status();
	    break;
	case SETS:
	    if (itemno < 0) {
		for (i = 0; i < MAXPLOT; i++) {
		}
	    }
	    break;
	case GRID:
	    if (itemno < 0) {
		for (i = 0; i < MAXGRIDS; i++) {
		    update_grid_status(i);
		}
	    }
	    break;
	case TEANL:
	    if (itemno < 0) {
		for (i = 0; i < MAXTEANL; i++) {
		    update_teanl_status(i);
		}
	    }
	    break;
	case ADCIRC:
	    if (itemno < 0) {
		for (i = 0; i < MAXADCIRC; i++) {
		    update_adcirc_status(i);
		}
	    }
	    break;
	case ELA:
	    if (itemno < 0) {
		for (i = 0; i < MAXELA; i++) {
		    update_ela_status(i);
		}
	    }
	    break;
	case GRAPHS:
	    if (itemno < 0) {
		for (i = 0; i < MAXGRAPH; i++) {
		    update_graph_status(i);
		}
	    }
	    break;
	}
    }
}

void update_graph_status(int gno)
{
    char *graph_types(int it, int which);

    if (gno >= 0 && gno < MAXGRAPH) {
	if (status_frame && cur_statusitem == GRAPHS) {
	    if (gno == cg) {
		sprintf(buf, "  %2d    %3s    %3s    %6s    %d [Current graph]", gno, on_or_off(g[gno].active), yes_or_no((!g[gno].hidden)), graph_types(g[gno].type, 0), g[gno].maxplot);
	    } else {
		sprintf(buf, "  %2d    %3s    %3s    %6s    %d", gno, on_or_off(g[gno].active), yes_or_no((!g[gno].hidden)), graph_types(g[gno].type, 0), g[gno].maxplot);
	    }
	    set_status_label(lab[gno], buf);
	}
    }
}

static void update_teanl_status(int flowno)
{
    char tmpbuf[20];
    if (flowno >= 0 && flowno < MAXTEANL) {
	if (status_frame && cur_statusitem == TEANL) {
	    if (g[cg].flowf[flowno].display == OFF) {
		strcpy(tmpbuf, "OFF");
	    } else if (g[cg].flowf[flowno].display == CENTER) {
		strcpy(tmpbuf, "CENTER");
	    } else if (g[cg].flowf[flowno].display == NODES) {
		strcpy(tmpbuf, "NODES");
	    } else {
		strcpy(tmpbuf, "UNKNOWN");
	    }
	    sprintf(buf, "  %2d    %3s   %6s   %d   %d   %d   %lf   %3s",
		    flowno + 1, on_or_off(flowf[flowno].active), tmpbuf, flowf[flowno].nread, flowf[flowno].nfreq, flowf[flowno].grid, flowf[flowno].dcor, flowf[curteanl].sign_of_phase == 1 ? "ADD" : "SUB");
	    set_status_label(lab[flowno], buf);
	}
    }
}

static void update_adcirc_status(int flowno)
{
    char tmpbuf[20];
    if (flowno >= 0 && flowno < MAXADCIRC) {
	if (status_frame && cur_statusitem == ADCIRC) {
	    if (g[cg].flowt[flowno].display == OFF) {
		strcpy(tmpbuf, "OFF");
	    } else if (g[cg].flowt[flowno].display == CENTER) {
		strcpy(tmpbuf, "CENTER");
	    } else if (g[cg].flowt[flowno].display == NODES) {
		strcpy(tmpbuf, "NODES");
	    }
	    sprintf(buf, "  %2d    %3s   %6s  %.1lf  %.1lf", flowno + 1, on_or_off(flowt[flowno].active), tmpbuf, flowt[flowno].nsteps, flowt[flowno].start, flowt[flowno].stop);
	    set_status_label(lab[flowno], buf);
	}
    }
}

static void update_adcirc3d_status(int flowno)
{
    char tmpbuf[20];
    if (flowno >= 0 && flowno < MAXADCIRC3D) {
	if (status_frame && cur_statusitem == ADCIRC3DFLOW) {
	    if (g[cg].flow3d[flowno].display == OFF) {
		strcpy(tmpbuf, "OFF");
	    } else if (g[cg].flow3d[flowno].display == CENTER) {
		strcpy(tmpbuf, "CENTER");
	    } else if (g[cg].flow3d[flowno].display == NODES) {
		strcpy(tmpbuf, "NODES");
	    }
	    if (adc3d[flowno].nsteps > 0) {
		sprintf(buf, "  %2d    %3s   %6s   %d  %.1lf  %.1lf", flowno + 1, on_or_off(adc3d[flowno].active), tmpbuf, adc3d[flowno].nsteps, adc3d[flowno].time[0], adc3d[flowno].time[adc3d[flowno].nsteps - 1]);
		set_status_label(lab[flowno], buf);
	    }
	}
    }
}

static void update_ela_status(int concno)
{
    char tmpbuf[20];
    if (concno >= 0 && concno < MAXELA) {
	if (status_frame && cur_statusitem == ELA) {
	    sprintf(buf, "  %2d    %3s   %6s    %d", concno + 1, on_or_off(elaconc[concno].active), on_or_off(g[cg].elaconc[concno].display), elaconc[concno].nsteps);
	    set_status_label(lab[concno], buf);
	}
    }
}

static void update_grid_status(int gridno)
{
    char tmpbuf[20];
    if (gridno >= 0 && gridno < MAXGRIDS) {
	if (status_frame && cur_statusitem == GRID) {
	    sprintf(buf, "  %2d    %3s   %6s  %d   %d", gridno + 1, on_or_off(grid[gridno].active), on_or_off(g[cg].grid[gridno].display), grid[gridno].nmel, grid[gridno].nmnp);
	    set_status_label(lab[gridno], buf);
	}
    }
}

static void update_drogues_status(int drno)
{
}

void update_model_status(int gno, int itemtype, int itemno)
{
    int i;

    if (status_frame) {
	switch (itemtype) {
	case ALL:
	    if (itemno < 0) {
		update_all_status();
	    }
	    break;
	case GRAPHS:
	    if (itemno < 0) {
		for (i = 0; i < MAXGRAPH; i++) {
		    update_graph_status(i);
		}
	    }
	    break;
	case GRID:
	    if (itemno < 0) {
		for (i = 0; i < MAXGRIDS; i++) {
		    update_grid_status(i);
		}
	    }
	    break;
	case TEANL:
	    if (itemno < 0) {
		for (i = 0; i < MAXTEANL; i++) {
		    update_teanl_status(i);
		}
	    }
	    break;
	case ADCIRC:
	    if (itemno < 0) {
		for (i = 0; i < MAXADCIRC; i++) {
		    update_adcirc_status(i);
		}
	    }
	    break;
	case ADCIRC3DFLOW:
	    if (itemno < 0) {
		for (i = 0; i < MAXADCIRC3D; i++) {
		    update_adcirc3d_status(i);
		}
	    }
	    break;
	case ELA:
	    if (itemno < 0) {
		for (i = 0; i < MAXELA; i++) {
		    update_ela_status(i);
		}
	    }
	    break;
	case DROGUES:
	    if (itemno < 0) {
		for (i = 0; i < MAXPATHLINES; i++) {
		    update_drogues_status(i);
		}
	    }
	    break;
	}
    }
}

static void update_all_status(void)
{
    int i, cnt = 0, didone = 0;
    char buf1[128], buf2[128], tmpbuf[128];
/* display status of grids */
    for (i = 0; i < MAXGRIDS; i++) {
	if (object_isactive(GRID, i)) {
	    sprintf(buf, "GRID #  %2d | active = %3s | display = %6s | elements, nodes = %d, %d", i + 1, on_or_off(grid[i].active), on_or_off(g[cg].grid[i].display), grid[i].nmel, grid[i].nmnp);
	    set_status_label(lab[cnt++], buf);
	}
    }
/* TEANL */
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i)) {
	    didone = 1;
	    sprintf(buf, "TEANL #%1d ACTIVE, attached to grid %d", i + 1, flowf[i].grid);
	    set_status_label(lab[cnt++], buf);
	}
    }
    if (!didone) {
	set_status_label(lab[cnt++], "No TEANL flows active");
    }
/* ADCIRC */
    didone = 0;
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, ELEV, i)) {
	    didone = 1;
	    sprintf(buf, "ADCIRC #%1d ELEVATION ACTIVE", i + 1);
	    set_status_label(lab[cnt++], buf);
	}
	if (object_isactive(ADCIRC, FLOW, i)) {
	    didone = 1;
	    sprintf(buf, "ADCIRC #%1d FLOW ACTIVE", i + 1);
	    set_status_label(lab[cnt++], buf);
	}
    }
    if (!didone) {
	set_status_label(lab[cnt++], "No ADCIRC flows active");
    }
/* ADCIRC 3D */
    didone = 0;
    for (i = 0; i < MAXADCIRC3D; i++) {
	if (object_isactive(ADCIRC3DFLOW, i)) {
	    didone = 1;
	    sprintf(buf, "ADCIRC 3D #%1d ACTIVE", i + 1);
	    set_status_label(lab[cnt++], buf);
	}
    }
    if (!didone) {
	set_status_label(lab[cnt++], "No ADCIRC flows active");
    }
/* ELA */
    didone = 0;
    for (i = 0; i < MAXELA; i++) {
	if (object_isactive(ELA, i)) {
	    sprintf(buf, "ELA #%1d ACTIVE, attached to grid %d", i + 1, elaconc[i].grid);
	    set_status_label(lab[cnt++], buf);
	}
    }
    if (!didone) {
	set_status_label(lab[cnt++], "No ELA concentrations active");
    }
/* DROGUES */
    for (i = 0; i < MAXPATHLINES; i++) {
    }
/* Time histories of elevation */
    for (i = 0; i < MAXHISTMARKERS; i++) {
    }
/* Time histories of velocities */
    for (i = 0; i < MAXVELHIST; i++) {
    }
}

void clear_status(void)
{
    int i;
    for (i = 0; i < MAXITEMS; i++) {
	set_status_label(lab[i], " ");
    }
}

void update_status_popup(void)
{
    if (status_frame) {
    }
    status_item_proc();
}

static void page_status_proc(void)
{
    curpage = (curpage + 1) % npages;
    update_status(cg, cur_statusitem, -1);
}

static void home_status_proc(void)
{
    curpage = 0;
    update_status(cg, cur_statusitem, -1);
}

static void end_status_proc(void)
{
    curpage = npages - 1;
    update_status(cg, cur_statusitem, -1);
}

void status_item_proc(void)
{
    int cd, i;
    if (status_frame) {
	cd = GetChoice(select_status_item);
	clear_status();

	switch (cd) {
	case 0:
	    cur_statusitem = ALL;
	    sprintf(header, " Item # Show  Type ");
	    break;
	case 1:
	    cur_statusitem = GRID;
	    sprintf(header, " GRID # Active  Show  Type  Elements Nodes");
	    break;
	case 2:
	    cur_statusitem = TEANL;
	    sprintf(header, " TEANL # Active  Show  Type Steps   Start   Stop   Step");
	    break;
	case 3:
	    cur_statusitem = ADCIRC;
	    sprintf(header, " ADCIRC # Active  Show  Type Steps   Start   Stop   Step");
	    break;
	case 4:
	    cur_statusitem = ADCIRC3DFLOW;
	    sprintf(header, " ADCIRC # Active  Show  Type Steps   Start   Stop   Step");
	    break;
	case 5:
	    cur_statusitem = ELA;
	    sprintf(header, " ELA # Active  Show  Type Steps   Start   Stop   Step");
	    break;
	case 6:
	    cur_statusitem = GRAPHS;
	    sprintf(header, " Graph # Active  Show  Type");
	    break;
	case 7:
	    cur_statusitem = REGIONS;
	    sprintf(header, " Region # Active  Type");
	    break;
	}
	set_status_label(header_w, header);
	update_status(cg, cur_statusitem, -1);
    }
}

/*
 * write the status to the results file
 */
static void update_stuff_status(void)
{
}

void select_set(Widget w, int cd)
{
    printf("Got it %d\n", cd);
}

void define_status_popup(void)
{
    extern Widget app_shell;
    extern Display *disp;
    int i;
    Widget wbut, rc, rc2, sw, fr, form;

    if (!status_frame) {
	status_frame = XmCreateDialogShell(app_shell, "Status", NULL, 0);

	f = (XFontStruct *) XLoadQueryFont(disp, "fixed");
	xmf = XmFontListCreate(f, XmSTRING_DEFAULT_CHARSET);

	status_panel = XmCreateForm(status_frame, "form", NULL, 0);

	sw = XtVaCreateManagedWidget("sw", xmScrolledWindowWidgetClass, status_panel, XmNscrollingPolicy, XmAUTOMATIC, NULL);
	rc2 = XmCreateRowColumn(sw, "rc2", NULL, 0);
	header_w = XtVaCreateManagedWidget("header", xmLabelWidgetClass, rc2, XmNalignment, XmALIGNMENT_BEGINNING, XmNfontList, xmf, XmNrecomputeSize, True, NULL);
	for (i = 0; i < MAXITEMS; i++) {
	    lab[i] = XtVaCreateManagedWidget("X", xmLabelWidgetClass, rc2, XmNalignment, XmALIGNMENT_BEGINNING, XmNfontList, xmf, XmNrecomputeSize, True, NULL);
	    XtAddEventHandler(lab[i], ButtonPressMask, False, (XtEventHandler) select_set, (XtPointer) i);
	}
	XtManageChild(rc2);
	XtVaSetValues(sw, XmNworkWindow, rc2, NULL);

	rc = XmCreateRowColumn(status_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) status_done_proc, 0);

	wbut = XtVaCreateManagedWidget("Update", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) update_status_popup, 0);

	wbut = XtVaCreateManagedWidget("Write", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) update_stuff_status, 0);

	wbut = XtVaCreateManagedWidget("Page", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) page_status_proc, 0);

	wbut = XtVaCreateManagedWidget("Home", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) home_status_proc, 0);

	wbut = XtVaCreateManagedWidget("End", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) end_status_proc, 0);

	select_status_item = CreatePanelChoice1(rc, "Display: ", 9, "ALL", "GRIDS", "TEANL", "ADCIRC", "ADCIRC3D", "ELA", "Graphs", "Regions", NULL, 0);
	for (i = 0; i < 8; i++) {
	    XtAddCallback(select_status_item[2 + i], XmNactivateCallback, (XtCallbackProc) status_item_proc, (XtPointer) i);
	}
	XtManageChild(rc);


	XtVaSetValues(sw, XmNtopAttachment, XmATTACH_FORM, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_WIDGET, XmNbottomWidget, rc, NULL);
	XtVaSetValues(rc,
/*
		      XmNtopAttachment, XmATTACH_WIDGET,
		      XmNtopWidget, sw,
*/
		      XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_FORM, NULL);
	XtManageChild(status_panel);

    }
    XtRaise(status_frame);
    update_status_popup();
}
