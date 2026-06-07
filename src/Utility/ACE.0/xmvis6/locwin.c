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
 * Locate nodes, elements, goto xy, node, element
 *
 */

#ifndef lint
static char RCSid[] = "$Id: locwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
extern Widget app_shell;

void create_points_frame(void);
void create_goto_popup(void);
void create_gotoxy_popup(void);

Widget locate_point_message;

static Widget goto_pointx_item;
static Widget goto_pointy_item;

extern Widget locate_grid_item;

static XmString label_string;

void SetLabel(Widget w, char *buf)
{
    Arg al;
    XmStringFree(label_string);
    label_string = XmStringCreateLtoR(buf, charset);
    XtSetArg(al, XmNlabelString, label_string);
    XtSetValues(w, &al, 1);
}

static void close_points_frame(Widget w, Widget cd)
{
    XtUnmanageChild(cd);
    set_action(0);
}

void create_points_frame(void)
{
    int x, y;
    static Widget top;
    Widget dialog, rc, wbut;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Locate", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    label_string = XmStringCreateLtoR("Set, location, (X, Y): ", charset);
    locate_point_message = XtVaCreateManagedWidget("pointslabel", xmLabelGadgetClass, dialog, XmNlabelString, label_string, NULL);

    locate_grid_item = XtVaCreateManagedWidget("locator", xmTextWidgetClass, dialog, XmNtraversalOn, True, XmNcolumns, 60, NULL);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) close_points_frame, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtRaise(top);
}

void open_locate1_popup(void)
{
    set_action(0);
    set_action(GET_NEAREST_NODE);
    create_points_frame();
    SetLabel(locate_point_message, "Nearest node:");
    write_mode_str("Use the left mouse button to display the node nearest the pointer");
}

void open_locate2_popup(void)
{
    set_action(0);
    set_action(GET_ELEMENT);
    create_points_frame();
    SetLabel(locate_point_message, "Element containing the pointer:");
    write_mode_str("Use the left mouse button to display the element containing the pointer");
}

void open_locate3_popup(void)
{
    set_action(0);
    set_action(GET_NEAREST_ELEMENT);
    create_points_frame();
    SetLabel(locate_point_message, "Nearest element:");
    write_mode_str("Use the left mouse button to display the element nearest the pointer");
}

void open_locate4_popup(void)
{
    set_action(0);
    set_action(GET_GRID_DEPTH);
    create_points_frame();
    SetLabel(locate_point_message, "Depth in the grid (interpolated):");
    write_mode_str("Press the left mouse button to display the depth of the edit grid at the pointer");
}

void open_locate6_popup(void)
{
    write_mode_str("Goto node/element");
    create_goto_popup();
}

void open_locate7_popup(void)
{
    write_mode_str("Goto X, Y");
    create_gotoxy_popup();
}

static Widget goto_frame;
static Widget goto_panel;
static Widget goto_node_item;
static Widget goto_elem_item;

static void goto_node_proc(void)
{
    int nodetofind;
    int sx, sy;

    nodetofind = atoi(xv_getstr(goto_node_item)) - 1;
    if (nodetofind < 0 || nodetofind > grid[g[cg].curgrid].nmnp) {
	errwin("Node number out of range");
	return;
    }
    world2deviceabs(grid[g[cg].curgrid].xord[nodetofind], grid[g[cg].curgrid].yord[nodetofind], &sx, &sy);
    setpointer(sx, sy);
}

static void goto_elem_proc(void)
{
    int elemtofind;
    int sx, sy;
    double centx, centy;

    elemtofind = atoi((char *) xv_getstr(goto_elem_item)) - 1;
    if (elemtofind < 0 || elemtofind > grid[g[cg].curgrid].nmel) {
	errwin("Element number out of range");
	return;
    }
    get_center(g[cg].curgrid, elemtofind, &centx, &centy);
    world2deviceabs(centx, centy, &sx, &sy);
    setpointer(sx, sy);
}

void display_goto_proc(void)
{
    create_goto_popup();
}

void create_goto_popup(void)
{
    int x, y;
    static Widget top;
    Widget dialog, rc, bt;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Goto node/element", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    goto_node_item = CreateTextItem2(dialog, 10, "Node:");
    bt = XtVaCreateManagedWidget("Goto node", xmPushButtonWidgetClass, dialog, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) goto_node_proc, NULL);

    goto_elem_item = CreateTextItem2(dialog, 10, "Element:");
    bt = XtVaCreateManagedWidget("Goto element", xmPushButtonWidgetClass, dialog, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) goto_elem_proc, NULL);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, dialog, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtRaise(top);
}

static void do_goto_proc(int item, int event)
{
    double wx, wy;
    int sx, sy;

    wx = atof((char *) xv_getstr(goto_pointx_item));
    wy = atof((char *) xv_getstr(goto_pointy_item));
    world2deviceabs(wx, wy, &sx, &sy);
    setpointer(sx, sy);
}

void create_gotoxy_popup(void)
{
    int i, x, y;
    static Widget top;
    Widget dialog, rc, wbut;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Goto X, Y", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    goto_pointx_item = CreateTextItem2(dialog, 10, "X: ");
    goto_pointy_item = CreateTextItem2(dialog, 10, "Y: ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_goto_proc, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtRaise(top);
}
