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
 * status popup
 *
 */

#ifndef lint
static char RCSid[] = "$Id: statuswin.c,v 1.4 2003/09/05 18:49:45 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

void do_status_proc(void);

static void create_wstatus_frame(void);
static void clear_status(void);

static Widget status_frame, status_panel, text_w;

/*
 * Create the status Panel
 */
void create_status_popup(void)
{
    int x, y;
    Widget wbut, rc, fr;

    if (status_frame) {
	XtRaise(status_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    status_frame = XmCreateDialogShell(app_shell, "Status", NULL, 0);
    handle_close(status_frame);
    XtVaSetValues(status_frame, XmNx, x, XmNy, y, NULL);
    status_panel = XmCreateForm(status_frame, "status_form", NULL, 0);
    fr = XmCreateFrame(status_panel, "fr", NULL, 0);
    text_w = XmCreateScrolledText(fr, "text_w", NULL, 0);
    XtVaSetValues(text_w, XmNrows, 10, XmNcolumns, 60, XmNeditMode, XmMULTI_LINE_EDIT, XmNwordWrap, True, NULL);
    XtManageChild(text_w);
    XtManageChild(fr);

    rc = XmCreateRowColumn(status_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Save...", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_wstatus_frame, 0);
    wbut = XtVaCreateManagedWidget("Clear", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) clear_status, 0);
    wbut = XtVaCreateManagedWidget("Update", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_status_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, status_frame);
    XtManageChild(rc);

    XtVaSetValues(fr, XmNtopAttachment, XmATTACH_FORM, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_WIDGET, XmNbottomWidget, rc, NULL);
    XtVaSetValues(rc, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_FORM, NULL);

    XtManageChild(status_panel);
    XtManageChild(status_frame);
    do_status_proc();
}

static void clear_status(void)
{
    if (inwin) {
	create_status_popup();
	XmTextSetString(text_w, " ");
    }
}

void stuffstatus(char *s, int sp)
{
    extern int inwin;
    static XmTextPosition pos = 0;
    static XmTextPosition savepos = 0;

    if (inwin) {
	create_status_popup();
	if (sp == 1) {
	    pos = XmTextGetLastPosition(text_w);
	    savepos = pos;
	    XmTextSetTopCharacter(text_w, savepos);
	}
	XmTextInsert(text_w, pos, s);
	pos += strlen(s);
	if (sp == 2) {
	    XmTextSetTopCharacter(text_w, savepos);
	    savepos = pos;
	}
    } else {
	printf(s);
    }
}

static void wstatus_apply_notify_proc(void);
static Widget wstatus_text_item;

static Widget wstatus_frame;

/*
 * Create the wparam Frame and the wparam Panel
 */
static void create_wstatus_frame(void)
{
    int x, y;
    extern Widget app_shell;
    Widget wbut, rc, wstatus_panel;
    int n = 0;
    Widget CreateTextItem();

    if (wstatus_frame) {
	XtRaise(wstatus_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    wstatus_frame = XmCreateDialogShell(app_shell, "Save", NULL, 0);
    handle_close(wstatus_frame);
    XtVaSetValues(wstatus_frame, XmNx, x, XmNy, y, NULL);
    wstatus_panel = XmCreateRowColumn(wstatus_frame, "wstatus_rc", NULL, 0);

    wstatus_text_item = CreateTextItem2(wstatus_panel, 20, "Save to file: ");
    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, wstatus_panel, NULL);
    rc = XmCreateRowColumn(wstatus_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wstatus_apply_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Cancel", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, wstatus_frame);
    XtManageChild(rc);

    XtManageChild(wstatus_panel);
    XtManageChild(wstatus_frame);
}

static void wstatus_apply_notify_proc(void)
{
    int i, len;
    char s[256];
    char *text, *fname;
    strcpy(s, (char *) xv_getstr(wstatus_text_item));
    if (!fexists(s)) {
	FILE *pp = fopen(s, "w");
	if (pp != NULL) {
	    text = XmTextGetString(text_w);
	    len = XmTextGetLastPosition(text_w);
	    fwrite(text, sizeof(char), len, pp);
	    fclose(pp);
	    XtFree(text);
	} else {
	    errwin("Unable to open file");
	}
    }
}

void do_status_proc(void)
{
    char s[280];
    double tmp, tarea, amin, amax, area();
    double element_volume(), v;
    int i, j, itmp, minel, maxel, bno;

    if (status_frame) {
	clear_status();
	sprintf(s, "\nNumber of build points %d\n\n", build[curbuild].nbuild);
	stuffstatus(s, 1);
	for (i = 0; i <= MAXGRIDS; i++) {
	    if (grid[i].gactive) {
		if (i == MAXGRIDS) {
		    sprintf(s, "\nBackground grid\n");
		    stuffstatus(s, 0);
		} else {
		    sprintf(s, "\nGrid #%d\n", i + 1);
		    stuffstatus(s, 0);
		}
		sprintf(s, "Alpha Id: %s\n", grid[i].alphid);
		stuffstatus(s, 0);
		sprintf(s, "Number of elements, nodes = %d, %d\n", grid[i].nmel, grid[i].nmnp);
		stuffstatus(s, 0);
		sprintf(s, "Matrix half bandwidth = %d\n", halfbandwidth(i));
		stuffstatus(s, 0);
		amax = area(i, 0);
		amin = area(i, 0);
		tarea = amax;
		minel = maxel = 0;
		for (j = 1; j < grid[i].nmel; j++) {
		    tmp = area(i, j);
		    tarea += tmp;
		    if (amin > tmp) {
			amin = tmp;
			minel = j + 1;
		    }
		    if (amax < tmp) {
			amax = tmp;
			maxel = j + 1;
		    }
		}
		v = 0.0;
		for (j = 0; j < grid[i].nmel; j++) {
		    v += element_volume(i, j);
		}
		sprintf(s, "Minimum element area %15lg at element #%d\n", amin, minel);
		stuffstatus(s, 0);
		sprintf(s, "Maximum element area %15lg at element #%d\n", amax, maxel);
		stuffstatus(s, 0);
		sprintf(s, "Total area = %15lg\n", tarea);
		stuffstatus(s, 0);
		sprintf(s, "Total volume = %15lg\n", v);
		stuffstatus(s, 0);
		if (region_flag) {
		    double cx, cy;
		    v = 0.0;
		    tarea = 0.0;
		    for (j = 0; j < grid[i].nmel; j++) {
			get_center(i, j, &cx, &cy);
			if (inregion(regionx, regiony, nregion, cx, cy)) {
			    v += element_volume(i, j);
			    tarea += area(i, j);
			}
		    }
		    sprintf(s, "Total area in region = %.5lf\n", tarea);
		    stuffstatus(s, 0);
		    sprintf(s, "Total volume in region = %.5lf\n", v);
		    stuffstatus(s, 0);
		}
		bno = grid[i].boundaries[0];
		if (bno >= 0 && grid[i].nbounds > 0) {
		    j = boundary[bno].nbpts;
		    sprintf(s, "Number of points in boundary 1 (external) = %d\n", j);
		    stuffstatus(s, 0);
		} else {
		    sprintf(s, "No external boundary defined\n");
		    stuffstatus(s, 0);
		}
		if (grid[i].nbounds <= 1) {
		    sprintf(s, "No internal boundaries defined\n");
		    stuffstatus(s, 0);
		} else {
		    for (j = 1; j < grid[i].nbounds; j++) {
			bno = grid[i].boundaries[j];
			if (bno >= 0 && grid[i].nbounds > 0) {
			    itmp = boundary[bno].nbpts;
			    sprintf(s, "Number of points in boundary %d (internal) = %d\n", j + 1, itmp);
			    stuffstatus(s, 0);
			}
		    }
		}
		sprintf(s, "\n");
		stuffstatus(s, 0);
	    } else {
/*
		sprintf(s, "\nGrid %d not active\n", i + 1);
		stuffstatus(s, 0);
*/
	    }
	}
	stuffstatus(" ", 2);
    }
}

int halfbandwidth(int gridno)
{
    int i, j, k, maxbw = 0, elbw;
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    for (k = 0; k < grid[gridno].icon[i].nn; k++) {
		if (j != k) {
		    elbw = (grid[gridno].icon[i].nl[j] - grid[gridno].icon[i].nl[k]);
		    elbw *= (elbw < 0) ? -1 : 1;
		    if (elbw > maxbw) {
			maxbw = elbw;
		    }
		}
	    }
	}
    }
    return maxbw;
}

static Widget info_frame, info_panel;

void info_done_proc(void)
{
    XtUnmanageChild(info_frame);
}

void create_info_popup(void)
{
    Widget bt, rc, fr;
    char buf[512];
    extern char xmgredit_version[];
    if (info_frame) {
	XtRaise(info_frame);
	return;
    }
    info_frame = XmCreateDialogShell(app_shell, "Info", NULL, 0);
    info_panel = XmCreateRowColumn(info_frame, "rc", NULL, 0);
    sprintf(buf, "%s", xmgredit_version);
    bt = XtVaCreateManagedWidget(buf, xmLabelWidgetClass, info_panel, NULL);
    sprintf(buf, "Copyright (c) Oregon Health and Science University");
    bt = XtVaCreateManagedWidget(buf, xmLabelWidgetClass, info_panel, NULL);
    sprintf(buf, "All Rights Reserved");
    bt = XtVaCreateManagedWidget(buf, xmLabelWidgetClass, info_panel, NULL);
    sprintf(buf, "See http://www.ccalmr.ogi.edu/CORIE/software for more information.");
    bt = XtVaCreateManagedWidget(buf, xmLabelWidgetClass, info_panel, NULL);

    rc = XmCreateRowColumn(info_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) info_done_proc, NULL);
    XtManageChild(rc);

    XtManageChild(info_panel);
    XtManageChild(info_frame);
}
