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
 * graph Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: graphwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
extern Widget app_shell;

static Widget graphs_frame;
static Widget graph_show_frame;
static Widget graph_focus_frame;
static Widget graph_type_frame;

static XmString gstring;

/*
 * Panel item declarations
 */
static Widget graph_curgraph_message_item;
static Widget *graph_setcur_choice_item;
static Widget *graph_rendsets_choice_item;
static Widget *graph_focus_choice_item;
static Widget graph_drawfocus_choice_item;
static Widget *graph_activate_choice_item;
static Widget *graph_kill_choice_item;
static Widget *graph_copyfrom_choice_item;
static Widget *graph_copyto_choice_item;
static Widget *graph_swapfrom_choice_item;
static Widget *graph_swapto_choice_item;
static Widget graph_show_choice_item[MAXGRAPH];

/*
 * Event and Notify proc declarations
 */
static void graph_setcur_notify_proc();
static void graph_rendsets_notify_proc(Widget w);
static void graph_focus_notify_proc(Widget w);
static void graph_activate_notify_proc(Widget w);
static void graph_kill_notify_proc(Widget w);
static void graph_copy_notify_proc(Widget w);
static void graph_show_notify_proc(Widget w);

static void update_focus_items(int gno);
static void update_show_items(void);
static void update_type_items(int gno);

void create_gactive_frame(void);
void create_gcopy_frame(void);
void create_gswap_frame(void);
void create_gkill_frame(void);
void create_gfocus_frame(void);
void create_gshow_frame(void);
void create_gtype_frame(void);
void create_arrange_frame(void);

static int gtypes[] = { XYFIXED, XY, LOGX, LOGY, LOGXY };

update_graph_items(void)
{
    update_focus_items(cg);
    update_show_items();
    update_type_items(cg);
}

void create_graph_frame(void)
{
    Widget wbut, graphs_panel;
    Arg a[2];
    int x, y, ac = 0;

    if (graphs_frame) {
	update_graph_items();
	XtRaise(graphs_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    graphs_frame = XmCreateDialogShell(app_shell, "Graph ops", NULL, 0);
    handle_close(graphs_frame);
    XtVaSetValues(graphs_frame, XmNx, x, XmNy, y, NULL);
    graphs_panel = XmCreateRowColumn(graphs_frame, "graphs_rc", NULL, 0);

    wbut = XtVaCreateManagedWidget("Activate...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gactive_frame, 0);

    wbut = XtVaCreateManagedWidget("Copy...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gcopy_frame, 0);

    wbut = XtVaCreateManagedWidget("Swap...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gswap_frame, 0);

    wbut = XtVaCreateManagedWidget("Kill...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gkill_frame, 0);

    wbut = XtVaCreateManagedWidget("Focus...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gfocus_frame, 0);

    wbut = XtVaCreateManagedWidget("Show...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gshow_frame, 0);

    wbut = XtVaCreateManagedWidget("Set graph type...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_gtype_frame, 0);

    wbut = XtVaCreateManagedWidget("Arrange...", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_arrange_frame, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, graphs_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, graphs_frame);

    update_graph_items();
    XtManageChild(graphs_panel);
    XtManageChild(graphs_frame);
}


/*
 * Notify and event procs
 */
static void graph_rendsets_notify_proc(Widget w)
{
    int graphtype = (int) GetChoice(graph_rendsets_choice_item);

    if (g[cg].type != gtypes[graphtype]) {
	g[cg].type = gtypes[graphtype];
	autoscale_graph(cg, -3);
	update_all(cg);
	drawgraph();
    }
}

static void graph_focus_notify_proc(Widget w)
{
    int newcg;

    switch ((int) GetChoice(graph_focus_choice_item)) {
    case 0:
	focus_policy = CLICK;
	break;
    case 1:
	focus_policy = SET;
	break;
    case 2:
	focus_policy = FOLLOWS;
	break;
    }
    draw_focus_flag = (int) XmToggleButtonGadgetGetState(graph_drawfocus_choice_item) ? ON : OFF;
    newcg = (int) GetChoice(graph_setcur_choice_item);
    if (newcg != cg) {
	switch_current_graph(cg, newcg);
    }
}

static void graph_activate_notify_proc(Widget w)
{
    int gno = (int) GetChoice(graph_activate_choice_item);

    set_graph_active(gno);
    update_all(cg);
    drawgraph();
}

static void graph_kill_notify_proc(Widget w)
{
    kill_graph((int) GetChoice(graph_kill_choice_item));
    update_all(cg);
    drawgraph();
}

static void graph_copy_notify_proc(Widget w)
{
    int from = (int) GetChoice(graph_copyfrom_choice_item);
    int to = (int) GetChoice(graph_copyto_choice_item);

    if (from == to) {
	errwin("Graph from and graph to are the same");
	return;
    }
    if (!isactive_graph(from)) {
	errwin("Graph from isn't active");
	return;
    }
    if (isactive_graph(to)) {
	if (!yesno("Graph to copy to is active, kill it?", "YES to kill or NO to abort", "YES", "NO")) {
	    return;
	}
    }
    copy_graph(from, to);
    update_all(cg);
    drawgraph();
}

static void graph_show_notify_proc(Widget w)
{
    int i;

    for (i = 0; i < MAXGRAPH; i++) {
	if (XmToggleButtonGadgetGetState(graph_show_choice_item[i]) == True) {
	    g[i].hidden = FALSE;
	} else {
	    g[i].hidden = TRUE;
	}
    }
    update_all(cg);
    drawgraph();
}

static void graph_swap_notify_proc(Widget w)
{
    int from = (int) GetChoice(graph_swapfrom_choice_item);
    int to = (int) GetChoice(graph_swapto_choice_item);

    if (from == to) {
	errwin("Graph from and graph to are the same");
	return;
    }
    swap_graph(from, to);
    update_all(cg);
    drawgraph();
}

/*
*/
static void update_type_items(int gno)
{
    int i;

    if (graph_type_frame) {
	i = 0;
	while (g[gno].type != gtypes[i])
	    i++;
	if (i > 8) {
	    errwin("Graph type not found");
	} else {
	    SetChoice(graph_rendsets_choice_item, i);
	}
    }
}

static void update_focus_items(int gno)
{
    int itest = 0;

    if (graph_focus_frame) {
	SetChoice(graph_setcur_choice_item, gno);
	if (focus_policy == SET) {
	    itest = 1;
	} else if (focus_policy == CLICK) {
	    itest = 0;
	} else if (focus_policy == FOLLOWS) {
	    itest = 2;
	}
	SetChoice(graph_focus_choice_item, itest);
	XmToggleButtonGadgetSetState(graph_drawfocus_choice_item, draw_focus_flag == ON ? True : False, False);
    }
}

static void update_show_items(void)
{
    int i, itest = 0;

    if (graph_show_frame) {
	for (i = 0; i < MAXGRAPH; i++) {
	    if (g[i].hidden) {
		XmToggleButtonGadgetSetState(graph_show_choice_item[i], False, False);
	    } else {
		XmToggleButtonGadgetSetState(graph_show_choice_item[i], True, False);
	    }
	}
    }
}

void create_gactive_frame(void)
{
    int i, x, y, ac = 0;
    Arg wargs[5];
    static Widget top, dialog;
    Widget lab, wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Activate graphs", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    graph_activate_choice_item = (Widget *) CreatePanelChoice1(dialog, "Activate graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_activate_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_gcopy_frame(void)
{
    int i, x, y, ac = 0;
    Arg wargs[5];
    static Widget top, dialog;
    Widget lab, wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Copy graphs", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    graph_copyfrom_choice_item = (Widget *) CreatePanelChoice1(dialog, "Copy graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    graph_copyto_choice_item = (Widget *) CreatePanelChoice1(dialog, "To graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_copy_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_gswap_frame(void)
{
    int i, x, y;
    static Widget top, dialog;
    Widget lab, wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Swap graphs", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    graph_swapfrom_choice_item = (Widget *) CreatePanelChoice1(dialog, "Swap graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    graph_swapto_choice_item = (Widget *) CreatePanelChoice1(dialog, "With graph: ", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_swap_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_gkill_frame(void)
{
    int i, x, y;
    static Widget top, dialog;
    Widget lab, wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Kill graphs", NULL, 0);
    handle_close(top);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    graph_kill_choice_item = (Widget *) CreatePanelChoice1(dialog, "Kill graph: ", 12, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "All", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_kill_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_gfocus_frame(void)
{
    int i, x, y;
    static Widget dialog;
    Widget lab, wbut, rc;

    if (graph_focus_frame) {
	update_focus_items(cg);
	XtRaise(graph_focus_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    graph_focus_frame = XmCreateDialogShell(app_shell, "Set focus", NULL, 0);
    handle_close(graph_focus_frame);
    XtVaSetValues(graph_focus_frame, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(graph_focus_frame, "rc", NULL, 0);

    graph_setcur_choice_item = (Widget *) CreatePanelChoice1(dialog, "Set current graph to", 11, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);

    graph_focus_choice_item = (Widget *) CreatePanelChoice1(dialog, "Graph focus", 4, "Button press", "As set", "Follows mouse", NULL, NULL);
    graph_drawfocus_choice_item = XtVaCreateManagedWidget("Display focus markers", xmToggleButtonGadgetClass, dialog, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_focus_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, graph_focus_frame);
    XtManageChild(rc);

    update_focus_items(cg);
    XtManageChild(dialog);
    XtManageChild(graph_focus_frame);
}

void create_gshow_frame(void)
{
    int i, x, y;
    static Widget dialog;
    Widget lab, wbut, rc;
    extern Widget app_shell;

    if (graph_show_frame) {
	update_show_items();
	XtRaise(graph_show_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    graph_show_frame = XmCreateDialogShell(app_shell, "Show graphs", NULL, 0);
    handle_close(graph_show_frame);
    XtVaSetValues(graph_show_frame, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(graph_show_frame, "rc", NULL, 0);

    lab = XtVaCreateManagedWidget("Select graphs for display (a graph must also be active to be seen):", xmLabelGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    for (i = 0; i < MAXGRAPH; i++) {
	sprintf(buf, "%d", i);
	graph_show_choice_item[i] = XtVaCreateManagedWidget(buf, xmToggleButtonGadgetClass, rc, NULL);
    }
    XtManageChild(rc);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_show_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, graph_show_frame);
    XtManageChild(rc);

    update_show_items();
    XtManageChild(dialog);
    XtManageChild(graph_show_frame);
}

void create_gtype_frame(void)
{
    int i, x, y;
    static Widget dialog;
    Widget lab, wbut, rc;
    extern Widget app_shell;

    if (graph_type_frame) {
	update_type_items(cg);
	XtRaise(graph_type_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    graph_type_frame = XmCreateDialogShell(app_shell, "Set graph type", NULL, 0);
    handle_close(graph_type_frame);
    XtVaSetValues(graph_type_frame, XmNx, x, XmNy, y, NULL);
    dialog = XmCreateRowColumn(graph_type_frame, "rc", NULL, 0);

    graph_rendsets_choice_item = (Widget *) CreatePanelChoice1(dialog, "Set current graph type to", 6, "Fixed scale XY", "XY graph", "Log-linear", "Linear-log", "Log-log", NULL, NULL);


    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) graph_rendsets_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, graph_type_frame);
    XtManageChild(rc);

    update_type_items(cg);
    XtManageChild(dialog);
    XtManageChild(graph_type_frame);
}
