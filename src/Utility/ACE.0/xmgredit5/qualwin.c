/*
 * ACE/gredit - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/*
 * Quality checks
 */

#ifndef lint
static char RCSid[] = "$Id: qualwin.c,v 1.4 2006/07/31 23:33:13 pturner Exp $";
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

typedef struct {
    Widget top;
    Widget angletoggle;
    Widget angle;
    Widget nnelstoggle;
    Widget *nnelsop;
    Widget nnels;
} QualUI;

static QualUI qui;

void create_qual_frame(Widget w, XtPointer clientd, XtPointer calld);
void accept_quals(Widget w, XtPointer clientd, XtPointer calld);

int display_maxangle = 0;
double maxAngle = 90.0;

int display_nnels = 0;
int nnelsop, nnels;

void drawMaxAngle(int gridno)
{
    double x[3], y[3], x1, y1, x2, y2, x3, y3;
    double a1, b1, a2, b2, a3, b3;
    double det, t1, t2, t3;
    int i, n1, n2, n3;
    if (display_maxangle) {
	setcolor(3);
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    x1 = grid[gridno].xord[n1];
	    y1 = grid[gridno].yord[n1];
	    n2 = grid[gridno].icon[i].nl[1];
	    x2 = grid[gridno].xord[n2];
	    y2 = grid[gridno].yord[n2];
	    n3 = grid[gridno].icon[i].nl[2];
	    x3 = grid[gridno].xord[n3];
	    y3 = grid[gridno].yord[n3];
	    a1 = x3 - x1;
	    b1 = y3 - y1;
	    a2 = x2 - x1;
	    b2 = y2 - y1;
	    det = sqrt((a1 * a1 + b1 * b1) * (a2 * a2 + b2 * b2));
	    t1 = acos((a1 * a2 + b1 * b2) / det) * 180.0 / M_PI;

	    a1 = x1 - x2;
	    b1 = y1 - y2;
	    a2 = x3 - x2;
	    b2 = y3 - y2;
	    det = sqrt((a1 * a1 + b1 * b1) * (a2 * a2 + b2 * b2));
	    t2 = acos((a1 * a2 + b1 * b2) / det) * 180.0 / M_PI;

	    a1 = x2 - x3;
	    b1 = y2 - y3;
	    a2 = x1 - x3;
	    b2 = y1 - y3;
	    det = sqrt((a1 * a1 + b1 * b1) * (a2 * a2 + b2 * b2));
	    t3 = acos((a1 * a2 + b1 * b2) / det) * 180.0 / M_PI;

	    if (t1 > maxAngle || t2 > maxAngle || t3 > maxAngle) {
		x[0] = x1;
		y[0] = y1;
		x[1] = x2;
		y[1] = y2;
		x[2] = x3;
		y[2] = y3;
		fillcolor(3, x, y);
	    }

	}
	setcolor(1);
    }
}


/*
 * Draw selected nodes with color
 */
void drawSelectedNodes(int gridno, int op, int nels, int color)
{
    int i, j, drawit;
    int bno;
    setcolor(4);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    grid[gridno].nlist[grid[gridno].icon[i].nl[j]]++;
	}
    }

    for (i = 0; i < grid[gridno].nmnp; i++) {
	drawit = 0;
	switch (op) {
	case 0:
	    if (grid[gridno].nlist[i] < nels) {
		grid[gridno].nlist[i] = 1;
	    } else {
		grid[gridno].nlist[i] = 0;
	    }
	    break;
	case 1:
	    if (grid[gridno].nlist[i] > nels) {
		grid[gridno].nlist[i] = 1;
	    } else {
		grid[gridno].nlist[i] = 0;
	    }
	    break;
	case 2:
	    if (grid[gridno].nlist[i] == nels) {
		grid[gridno].nlist[i] = 1;
	    } else {
		grid[gridno].nlist[i] = 0;
	    }
	    break;
	default:
	    grid[gridno].nlist[i] = 0;
	    drawit = 0;
	    break;
	}
	if (drawit && symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
	    my_circlefilled(grid[gridno].xord[i], grid[gridno].yord[i]);
	}
    }

    for (i = 0; i < grid[gridno].nbounds; i++) {
	bno = grid[gridno].boundaries[i];
	if (bno >= 0 && boundary[bno].bactive) {
	    for (j = 0; j < boundary[bno].nbpts; j++) {
		grid[gridno].nlist[boundary[bno].nodes[j]] = 0;
	    }
	}
    }
    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].nlist[i] && symok(grid[gridno].xord[i], grid[gridno].yord[i])) {
	    my_circlefilled(grid[gridno].xord[i], grid[gridno].yord[i]);
	}
    }
    setcolor(1);
}

void accept_quals(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256];
    QualUI *q = (QualUI *) clientd;
    display_maxangle = XmToggleButtonGetState(q->angletoggle);
    if (display_maxangle) {
	strcpy(buf, (char *) xv_getstr(q->angle));
	maxAngle = atof(buf);
    }
    display_nnels = XmToggleButtonGetState(q->nnelstoggle);
    if (display_nnels) {
	nnelsop = (int) GetChoice(q->nnelsop);
	strcpy(buf, (char *) xv_getstr(q->nnels));
	nnels = atoi(buf);
    }
}

void create_qual_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int x, y;
    Widget dialog;
    Widget wbut, lab, rc, rcl, rc1, rc2, form;

    set_wait_cursor();
    if (qui.top == NULL) {
	char *label1[2];
	Widget but1[2];
	label1[0] = "Accept";
	label1[1] = "Close";

	XmGetPos(app_shell, 0, &x, &y);
	qui.top = XmCreateDialogShell(app_shell, "Quality checks", NULL, 0);
	handle_close(qui.top);
	XtVaSetValues(qui.top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(qui.top, "rc", NULL, 0);

	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	qui.angletoggle = XmCreateToggleButton(rc, "Display elements", NULL, 0);
	XtManageChild(qui.angletoggle);
	qui.angle = CreateTextItem2(rc, 10, "having a maximum angle >= ");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	qui.nnelstoggle = XmCreateToggleButton(rc, "Display nodes attached to ", NULL, 0);
	XtManageChild(qui.nnelstoggle);
	qui.nnelsop = CreatePanelChoice1(rc, "", 4, " < ", " > ", " == ", NULL, NULL);

	qui.nnels = CreateTextItem2(rc, 10, "");
	XtVaCreateManagedWidget(" elements", xmLabelGadgetClass, rc, NULL);
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) accept_quals, (XtPointer) & qui);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) qui.top);
	XtManageChild(dialog);
    }
    XtRaise(qui.top);
    unset_wait_cursor();
}

typedef struct {
    Widget top;
    Widget angletoggle;
    Widget angle;
    Widget nnelstoggle;
    Widget *nnelsop;
    Widget nnels;
} VolUI;

static VolUI vui;

void create_vol_frame(Widget w, XtPointer clientd, XtPointer calld);
void accept_vol(Widget w, XtPointer clientd, XtPointer calld);

void accept_vol(Widget w, XtPointer clientd, XtPointer calld)
{
    char buf[256], s[1024];
    int i,j;
    double vol = 0.0;
    VolUI *v = (VolUI *) clientd;
    for (j = 0; j < grid[curgrid].nmel; j++) {
      if (region_flag) {
        double cx, cy;
             get_center(i, j, &cx, &cy);
             if (inregion(regionx, regiony, nregion, cx, cy)) {
                 vol += element_volume(i, j);
             }
        sprintf(s, "Total volume in region = %.5lf\n", vol);
      } else {
        for (j = 0; j < grid[curgrid].nmel; j++) {
            vol += element_volume(i, j);
        }
        sprintf(s, "Total volume = %15lg\n", vol);
    }
    }
}

void create_vol_frame(Widget w, XtPointer clientd, XtPointer calld)
{
    int x, y;
    Widget dialog;
    Widget wbut, lab, rc, rcl, rc1, rc2, form;

    set_wait_cursor();
    if (vui.top == NULL) {
        char *label1[2];
        Widget but1[2];
        label1[0] = "Accept";
        label1[1] = "Close";

        XmGetPos(app_shell, 0, &x, &y);
        vui.top = XmCreateDialogShell(app_shell, "Volume/depth checks", NULL, 0);
        handle_close(vui.top);
        XtVaSetValues(vui.top, XmNx, x, XmNy, y, NULL);
        dialog = XmCreateRowColumn(vui.top, "rc", NULL, 0);

        rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
        XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
        vui.angletoggle = XmCreateToggleButton(rc, "Display elements", NULL, 0);
        XtManageChild(vui.angletoggle);
        vui.angle = CreateTextItem2(rc, 10, "having a maximum angle >= ");
        XtManageChild(rc);

        XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

        rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
        XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
        vui.nnelstoggle = XmCreateToggleButton(rc, "Display nodes attached to ", NULL, 0);
        XtManageChild(vui.nnelstoggle);
        vui.nnelsop = CreatePanelChoice1(rc, "", 4, " < ", " > ", " == ", NULL, NULL);
        vui.nnels = CreateTextItem2(rc, 10, "");
        XtVaCreateManagedWidget(" elements", xmLabelGadgetClass, rc, NULL);
        XtManageChild(rc);

        XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

        CreateCommandButtons(dialog, 2, but1, label1);
        XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) accept_vol, (XtPointer) & vui);
        XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) vui.top);
        XtManageChild(dialog);
    }
    XtRaise(vui.top);
    unset_wait_cursor();
}

