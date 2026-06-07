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
 * modify nodes/elements popup
 *
 */

#ifndef lint
static char RCSid[] = "$Id: modwin.c,v 1.4 2003/08/06 20:51:50 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;

static Widget mod_frame;
static Widget mod_panel;

static void do_move_node(void);
static void do_add_element(void);
static void do_addq_element(void);
static void do_add_node(void);
static void do_delete_element(void);
static void do_delete_nearest_element(void);
static void set_swap_line(void);
static void do_snap_line(void);
static void do_split3(void);
static void do_split4(void);
static void do_makequad(void);
static void do_splitquad(void);

void set_modwin(void);
void mod_done_proc(void);

static Widget goto_elem_item;
static Widget goto_node_item;

static void goto_node_proc(void)
{
    int nodetofind;
    int sx, sy;

    nodetofind = atoi((char *) panel_getstr_value(goto_node_item)) - 1;
    if (nodetofind < 0 || nodetofind > grid[curgrid].nmnp) {
	errwin("Node number out of range");
	return;
    }
    get_device(grid[curgrid].xord[nodetofind], grid[curgrid].yord[nodetofind], &sx, &sy);
    setpointer(sx, sy);
}

static void goto_elem_proc(void)
{
    int elemtofind;
    int sx, sy;
    double centx, centy;

    elemtofind = atoi((char *) panel_getstr_value(goto_elem_item)) - 1;
    if (elemtofind < 0 || elemtofind > grid[curgrid].nmel) {
	errwin("Element number out of range");
	return;
    }
    get_center(curgrid, elemtofind, &centx, &centy);
    get_device(centx, centy, &sx, &sy);
    setpointer(sx, sy);
}

void create_mod_popup(void)
{
    Widget bt, fr, rc;
    Arg a[2];
    int ac;
    if (mod_frame) {
	XtRaise(mod_frame);
	return;
    }
    mod_frame = XmCreateDialogShell(app_shell, "Edit nodes/elements", NULL, 0);
    mod_panel = XmCreateRowColumn(mod_frame, "rc", NULL, 0);

    bt = XtVaCreateManagedWidget("Move node", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_move_node, NULL);

    bt = XtVaCreateManagedWidget("Add node", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_add_node, NULL);

    bt = XtVaCreateManagedWidget("Add triangle", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_add_element, NULL);

    bt = XtVaCreateManagedWidget("Add quad", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_addq_element, NULL);

    bt = XtVaCreateManagedWidget("Delete element", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_delete_element, NULL);

    bt = XtVaCreateManagedWidget("Delete nearest element", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_delete_nearest_element, NULL);

    bt = XtVaCreateManagedWidget("Snap nodes to line", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_snap_line, NULL);

    bt = XtVaCreateManagedWidget("Swap lines", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) set_swap_line, NULL);

    bt = XtVaCreateManagedWidget("Split 3", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_split3, NULL);

    bt = XtVaCreateManagedWidget("Split 4", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_split4, NULL);

    bt = XtVaCreateManagedWidget("Make Quadrangle", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_makequad, NULL);

    bt = XtVaCreateManagedWidget("Split Quadrangle", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_splitquad, NULL);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, mod_panel,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) mod_done_proc, NULL);
    XtManageChild(mod_panel);
    XtManageChild(mod_frame);
}

void set_modwin(void)
{
}

void mod_done_proc(void)
{
    XtUnmanageChild(mod_frame);
}

static void do_move_node(void)
{
    set_action(0);
    set_action(MOVE_NODE1ST);
}

static void do_add_node(void)
{
    set_action(0);
    set_action(ADD_NODE);
}

static void do_add_element(void)
{
    set_action(0);
    set_action(ADD_ELEMENT1);
}

static void do_addq_element(void)
{
    set_action(0);
    set_action(ADD_QUAD1);
}

static void do_delete_element(void)
{
    set_action(0);
    set_action(DELETE_ELEMENTS);
}

static void do_delete_nearest_element(void)
{
    set_action(0);
    set_action(DELETE_ELEMENT);
}

static void do_snap_line(void)
{
    set_action(0);
    set_action(SNAP_LINE1ST);
}

static void set_swap_line(void)
{
    set_action(0);
    set_action(SWAPLINE_1ST);
}

static void do_split3(void)
{
    set_action(0);
    set_action(SPLIT_ELEMENT3);
}

static void do_split4(void)
{
    set_action(0);
    set_action(SPLIT_ELEMENT4);
}

static void do_makequad(void)
{
    set_action(0);
    set_action(CONV_3TO4_1ST);
}

static void do_splitquad(void)
{
    set_action(0);
    set_action(CONV_4TO3);
}
