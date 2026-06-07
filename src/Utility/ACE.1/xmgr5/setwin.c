/* $Id: setwin.c,v 1.4 2006/07/20 16:42:28 pturner Exp $
 *
 * setops - operations on sets
 *
 */

#include <stdio.h>
#include <math.h>

#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/FileSB.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/ToggleB.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/List.h>
#include <Xm/Separator.h>
#include <Xm/Protocols.h>

#include "globals.h"
#include "motifinc.h"

static Widget but1[2];
static Widget but2[3];

char format[128] = "%16lg %16lg";
char sformat[128] = "%16lg %16lg";

static Widget setops_frame;
static Widget setops_panel;
static Widget pickops_frame;
static Widget pickops_panel;
static Widget *activate_set_item;
static Widget *set_settype_item;
static Widget alength_to_set_item;
static Widget *deactivate_set_item;
static Widget *reactivate_set_item;
static Widget *set_length_item;
static Widget length_to_set_item;
static Widget *change_set_item;
static Widget *change_settype_item;
static Widget *copy_from_item;
static Widget *copy_graph_item;
static Widget *copy_to_item;
static Widget *move_from_item;
static Widget *move_graph_item;
static Widget *move_to_item;
static Widget *drop_points_item;
static Widget drop_start_item;
static Widget drop_end_item;
static Widget *join_from_item;
static Widget *join_to_item;
static Widget *part_set_item;
static Widget part_len_item;
static Widget *kill_set_item;
static Widget kill_soft_toggle;
static Widget *sort_set_item;
static Widget *sort_xy_item;
static Widget *sort_up_down_item;
static Widget *write_sets_item;
static Widget write_imbed_toggle;
static Widget *write_binary_toggle;
static Widget write_sets_format_item;
static Widget write_sets_file_item;
static Widget saveall_sets_format_item;
static Widget saveall_sets_file_item;
static Widget *reverse_sets_item;
static Widget *coalesce_sets_item;
static Widget *swap_from_item;
static Widget *swap_fgraph_item;
static Widget *swap_to_item;
static Widget *swap_tgraph_item;
static Widget *write_setsgraph_item;

static void do_activate_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_deactivate_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_reactivate_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_setlength_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_changetype_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_copy_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_setmove_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_drop_points_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_join_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_split_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_kill_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_sort_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_write_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_saveall_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_reverse_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_coalesce_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_swap_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_break_sets_proc(Widget w, XtPointer client_data, XtPointer call_data);

/* for pick ops */
static void define_pickops_popup(Widget w, XtPointer client_data, XtPointer call_data);

static void create_activate_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_deactivate_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_reactivate_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_change_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_copy_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_setlength_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_move_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_drop_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_join_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_split_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_break_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_kill_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_sort_popup(Widget w, XtPointer client_data, XtPointer call_data);
void create_write_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_reverse_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_coalesce_popup(Widget w, XtPointer client_data, XtPointer call_data);
static void create_swap_popup(Widget w, XtPointer client_data, XtPointer call_data);

static char errbuf[256];

extern int index_set_types[];	/* declared in setutils.c */
extern int index_set_ncols[];

void define_setops_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget wbut, rc;
    int x, y;
    set_wait_cursor();
    if (setops_frame == NULL) {
	XmGetPos(app_shell, 0, &x, &y);
	setops_frame = XmCreateDialogShell(app_shell, "Set ops", NULL, 0);
	handle_close(setops_frame);
	XtVaSetValues(setops_frame, XmNx, x, XmNy, y, NULL);
	setops_panel = XmCreateRowColumn(setops_frame, "setops_rc", NULL, 0);
	XtVaSetValues(setops_panel, XmNorientation, XmHORIZONTAL, NULL);

	rc = XmCreateRowColumn(setops_panel, "rc", NULL, 0);
	wbut = XtVaCreateManagedWidget("Pick ops...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_pickops_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Activate...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_activate_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("De-activate...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_deactivate_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Re-activate...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_reactivate_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Set length...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_setlength_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Change set type...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_change_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Copy...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_copy_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Move...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_move_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Drop points...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_drop_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Join sets...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_join_popup, (XtPointer) NULL);
	XtManageChild(rc);

	rc = XmCreateRowColumn(setops_panel, "rc", NULL, 0);
	wbut = XtVaCreateManagedWidget("Split...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_split_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Kill...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_kill_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Kill all", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_flush, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Sort...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_sort_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Write sets...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_write_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Reverse sets...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_reverse_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Coalesce sets...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_coalesce_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Swap sets...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_swap_popup, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Pack sets", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_packsets, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) setops_frame);
	XtManageChild(rc);

	XtManageChild(setops_panel);
    }
    XtRaise(setops_frame);
    unset_wait_cursor();
}

static void define_pickops_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget wbut;
    int x, y;

    set_wait_cursor();
    if (pickops_frame == NULL) {
	XmGetPos(app_shell, 0, &x, &y);
	pickops_frame = XmCreateDialogShell(app_shell, "Pick ops", NULL, 0);
	handle_close(pickops_frame);
	XtVaSetValues(pickops_frame, XmNx, x, XmNy, y, NULL);
	pickops_panel = XmCreateRowColumn(pickops_frame, "pickops_rc", NULL, 0);

	wbut = XtVaCreateManagedWidget("Kill nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) KILL_NEAREST);

	wbut = XtVaCreateManagedWidget("Copy nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) COPY_NEAREST1ST);

	wbut = XtVaCreateManagedWidget("Move nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) MOVE_NEAREST1ST);

	wbut = XtVaCreateManagedWidget("Reverse nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) REVERSE_NEAREST);

	wbut = XtVaCreateManagedWidget("De-activate nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) DEACTIVATE_NEAREST);

	wbut = XtVaCreateManagedWidget("Join nearest sets",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) JOIN_NEAREST1ST);

	wbut = XtVaCreateManagedWidget("Delete range in nearest set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) DELETE_NEAREST1ST);

	wbut = XtVaCreateManagedWidget("Break set",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) PICK_BREAK);

	wbut = XtVaCreateManagedWidget("Cancel operation",
				     xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) 0);

	wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, pickops_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog,
		      (XtPointer) pickops_frame);

	XtManageChild(pickops_panel);
    }
    XtRaise(pickops_frame);
    unset_wait_cursor();
}

static void create_activate_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;
    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Activate set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	activate_set_item = CreateSetChoice(dialog, "Activate set :", maxplot, 0);
	alength_to_set_item = CreateTextItem2(dialog, 10, "Length to set:");
	set_settype_item = CreatePanelChoice(dialog,
					     "Set type:",
					     11,
					     "XY",
					     "XY DX",
					     "XY DY",
					     "XY DX1 DX2",
					     "XY DY1 DY2",
					     "XY DX DY",
					     "XY Z",
					     "XY HILO",
					     "XY R",
					     "XY U V",
					     NULL, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_activate_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_deactivate_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "De-activate set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	deactivate_set_item = CreateSetChoice(dialog, "De-activate set :", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_deactivate_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_reactivate_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Re-activate set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	reactivate_set_item = CreateSetChoice(dialog, "Re-activate set :", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_reactivate_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_change_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Change set type", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	change_set_item = CreateSetChoice(dialog, "Change the type of set :", maxplot, 0);

	change_settype_item = CreatePanelChoice(dialog,
						"To type:",
						11,
						"XY",
						"XY DX",
						"XY DY",
						"XY DX1 DX2",
						"XY DY1 DY2",
						"XY DX DY",
						"XY Z",
						"XY HILO",
						"XY R",
						"XY U V",
						NULL, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_changetype_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
	XtManageChild(top);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_setlength_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Set length", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	set_length_item = CreateSetChoice(dialog, "Set length of set: ", maxplot, 0);
	length_to_set_item = CreateTextItem2(dialog, 10, "Length:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_setlength_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_copy_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Copy set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	copy_from_item = CreateSetChoice(dialog, "Copy set:", maxplot, 0);
	copy_to_item = CreateSetChoice(dialog, "To set:", maxplot, 4);

	copy_graph_item = CreateGraphChoice(dialog, "In graph:", maxgraph, 1);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_copy_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_move_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Move set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	move_from_item = CreateSetChoice(dialog, "Move set:", maxplot, 0);
	move_to_item = CreateSetChoice(dialog, "To set:", maxplot, 4);
	move_graph_item = CreateGraphChoice(dialog, "In graph:", maxgraph, 1);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_setmove_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_drop_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Drop points", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	drop_points_item = CreateSetChoice(dialog, "Drop points from set:", maxplot, 0);
	drop_start_item = CreateTextItem2(dialog, 10, "Start drop:");
	drop_end_item = CreateTextItem2(dialog, 10, "End drop:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_drop_points_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_join_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Join sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	join_from_item = CreateSetChoice(dialog, "Join set:", maxplot, 0);
	join_to_item = CreateSetChoice(dialog, "To set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_join_sets_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_split_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Split sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	part_set_item = CreateSetChoice(dialog, "Split set:", maxplot, 0);
	part_len_item = CreateTextItem2(dialog, 10, "Length:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_split_sets_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_kill_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Kill set", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	kill_set_item = CreateSetChoice(dialog, "Kill set:", maxplot, 1);
	kill_soft_toggle = XtVaCreateManagedWidget("Preserve parameters",
					  xmToggleButtonWidgetClass, dialog,
						   NULL);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_kill_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_sort_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Sort sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	sort_set_item = CreateSetChoice(dialog, "Sort set:", maxplot, 0);
	sort_xy_item = CreatePanelChoice(dialog,
					 "Sort on:",
					 9,
					 "X",
					 "Y", 
					 "Y1", 
					 "Y2", 
					 "Y3", 
					 "Y4", 
					 "Y5", 
					 "Y6", 
					0, 0);
	sort_up_down_item = CreatePanelChoice(dialog,
					      "Order:",
					      3,
					      "Ascending",
					      "Descending", 0,
					      0);
	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_sort_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

void create_write_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top;
    Widget dialog;
    Widget wbut, rc, fr;

    set_wait_cursor();
    if (top == NULL) {
	top = XmCreateFileSelectionDialog(app_shell, "write_sets", NULL, 0);
	XtVaSetValues(XtParent(top), XmNtitle, "Write sets", NULL);

	XtAddCallback(top, XmNokCallback, (XtCallbackProc) do_write_sets_proc, (XtPointer) top);
	XtAddCallback(top, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	fr = XmCreateFrame(top, "fr", NULL, 0);
	dialog = XmCreateRowColumn(fr, "dialog_rc", NULL, 0);

	write_sets_item = CreateSetChoice(dialog, "Write set:", maxplot, 1);
	write_setsgraph_item = CreateGraphChoice(dialog, "From graph:", maxgraph, 2);
	write_imbed_toggle = XtVaCreateManagedWidget("Imbed parameters",
					  xmToggleButtonWidgetClass, dialog,
						     NULL);
	write_binary_toggle = CreatePanelChoice(dialog,
					 "Output format:", 4,
					 "ASCII",
					 "Binary", 
					 "Netcdf", 0, 0);
	write_sets_format_item = CreateTextItem2(dialog, 15, "Format: ");

	XtManageChild(dialog);
	XtManageChild(fr);
	xv_setstr(write_sets_format_item, format);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_reverse_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Reverse sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	reverse_sets_item = CreateSetChoice(dialog, "Reverse set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_reverse_sets_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_coalesce_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Coalesce sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	coalesce_sets_item = CreateSetChoice(dialog,
			    "Coalesce all active sets to set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_coalesce_sets_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

static void create_swap_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Swap sets", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	swap_from_item = CreateSetChoice(dialog, "Swap set:", maxplot, 0);
	swap_fgraph_item = CreateGraphChoice(dialog, "In graph:", maxgraph, 1);

	swap_to_item = CreateSetChoice(dialog, "With set:", maxplot, 0);
	swap_tgraph_item = CreateGraphChoice(dialog, "In graph:", maxgraph, 1);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_swap_proc, (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/*
 * setops - combine, copy sets - callbacks
*/

/*
 * activate a set and set its length
 */
static void do_activate_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, len, type;
    setno = (int) GetChoice(activate_set_item);
    type = GetChoice(set_settype_item);
    len = atoi((char *) xv_getstr(alength_to_set_item));
    set_wait_cursor();
    do_activate(setno, type, len);
    unset_wait_cursor();
}

/*
 * de-activate a set
 */
static void do_deactivate_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno;
    setno = (int) GetChoice(deactivate_set_item);
    set_wait_cursor();
    do_deactivate(cg, setno);
    unset_wait_cursor();
}

/*
 * re-activate a set
 */
static void do_reactivate_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno;
    setno = (int) GetChoice(reactivate_set_item);
    set_wait_cursor();
    do_reactivate(cg, setno);
    unset_wait_cursor();
}

/*
 * change the type of a set
 */
static void do_changetype_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, type;
    setno = (int) GetChoice(change_set_item);
    type = GetChoice(change_settype_item);
    set_wait_cursor();
    do_changetype(setno, type);
    unset_wait_cursor();
}

/*
 * set the length of an active set - contents are destroyed
 */
static void do_setlength_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, len;
    setno = (int) GetChoice(set_length_item);
    len = atoi((char *) xv_getstr(length_to_set_item));
    set_wait_cursor();
    do_setlength(setno, len);
    unset_wait_cursor();
}

/*
 * copy a set to another set, if the to set doesn't exist
 * get a new one, if it does, ask if it is okay to overwrite
 */
static void do_copy_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int j1, j2, gto;
    j1 = (int) GetChoice(copy_from_item);
    j2 = (int) GetChoice(copy_to_item);
    gto = (int) GetChoice(copy_graph_item);
    set_wait_cursor();
    do_copy(j1, cg, j2, gto);
    unset_wait_cursor();
}

/*
 * move a set to another set, if the to set doesn't exist
 * get a new one, if it does, ask if it is okay to overwrite
 */
static void do_setmove_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int j1, j2, gto;

    j1 = (int) GetChoice(move_from_item);
    gto = (int) GetChoice(move_graph_item);
    j2 = (int) GetChoice(move_to_item) - 1;
    set_wait_cursor();
    do_move(j1, cg, j2, gto);
    unset_wait_cursor();
}

/*
 * swap a set with another set
 */
static void do_swap_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int j1, j2, gto, gfrom;

    j1 = (int) GetChoice(swap_from_item);
    gfrom = (int) GetChoice(swap_fgraph_item);
    j2 = (int) GetChoice(swap_to_item);
    gto = (int) GetChoice(swap_tgraph_item);
    set_wait_cursor();
    do_swap(j1, gfrom, j2, gto);
    unset_wait_cursor();
}

/*
 * drop points from an active set
 */
static void do_drop_points_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int startno, endno, setno;
    setno = (int) GetChoice(drop_points_item);
    startno = atoi((char *) xv_getstr(drop_start_item)) - 1;
    endno = atoi((char *) xv_getstr(drop_end_item)) - 1;
    set_wait_cursor();
    do_drop_points(setno, startno, endno);
    unset_wait_cursor();
}

/*
 * append one set to another
 */
static void do_join_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int j1, j2;
    j1 = (int) GetChoice(join_from_item);
    j2 = (int) GetChoice(join_to_item);
    set_wait_cursor();
    do_join_sets(cg, j1, cg, j2);
    unset_wait_cursor();
}

/*
 * reverse the order of a set
 */
static void do_reverse_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno;
    setno = (int) GetChoice(reverse_sets_item);
    set_wait_cursor();
    do_reverse_sets(setno);
    unset_wait_cursor();
}

/*
 * coalesce sets
 */
static void do_coalesce_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno;
    setno = (int) GetChoice(coalesce_sets_item);
    set_wait_cursor();
    do_coalesce_sets(setno);
    unset_wait_cursor();
}

/*
 * kill a set
 */
static void do_kill_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, soft;

    setno = (int) GetChoice(kill_set_item);
    if (setno == maxplot) {
	setno = -1;
    }
    soft = (int) XmToggleButtonGetState(kill_soft_toggle);
    set_wait_cursor();
    do_kill(cg, setno, soft);
    unset_wait_cursor();
}

/*
 sort sets
*/
static void do_sort_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, sorton, stype;
    int son[] = {X, Y, Y1, Y2, Y3, Y4, Y5};
    setno = (int) GetChoice(sort_set_item);
    if (setno == maxplot) {
	setno = -1;
    }
    sorton = son[(int) GetChoice(sort_xy_item)];
    stype = (int) GetChoice(sort_up_down_item);
    set_wait_cursor();
    do_sort(setno, sorton, stype);
    unset_wait_cursor();
}

/*
 *  write a set or sets to a file
 */
static void do_write_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int which_graph;
    int setno;
    int imbed, bin;
    char fn[512], *s;

    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR (cbs->value, charset, &s)) {
        errwin("Error converting XmString to char string in do_write_sets_proc()");
        return;
    }

    strcpy(fn, s);
    XtFree(s);

    imbed = (int) XmToggleButtonGetState(write_imbed_toggle);
    bin = (int) GetChoice(write_binary_toggle);
    setno = (int) GetChoice(write_sets_item);
    if (setno == maxplot) {
	setno = -1;
    }
    which_graph = (int) GetChoice(write_setsgraph_item) - 1;
    strcpy(format, (char *) xv_getstr(write_sets_format_item));
    set_wait_cursor();
    switch (bin) {
	case 0:
	do_writesets(which_graph, setno, imbed, fn, format);
	break;
	case 1:
	do_writesets_binary(cg, setno, fn);
	break;
#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)
	case 2:
	writesets_netcdf(fn);
	break;
#endif
    }
    unset_wait_cursor();
}

/*
 * split sets split by itmp, remainder in last set.
 */
static void do_split_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, lpart;

    setno = (int) GetChoice(part_set_item);
    lpart = atoi((char *) xv_getstr(part_len_item));
    set_wait_cursor();
    do_splitsets(cg, setno, lpart);
    unset_wait_cursor();
}

void create_saveall_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top;
    Widget dialog;
    Widget wbut, rc, fr;

    set_wait_cursor();
    if (top == NULL) {
	top = XmCreateFileSelectionDialog(app_shell, "save_all_sets", NULL, 0);
	XtVaSetValues(XtParent(top), XmNtitle, "Save all sets", NULL);

	XtAddCallback(top, XmNokCallback, (XtCallbackProc) do_saveall_sets_proc, (XtPointer) top);
	XtAddCallback(top, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	fr = XmCreateFrame(top, "fr", NULL, 0);
	dialog = XmCreateRowColumn(fr, "dialog_rc", NULL, 0);

	saveall_sets_format_item = CreateTextItem2(dialog, 15, "Format: ");

	XtManageChild(dialog);
	XtManageChild(fr);

	xv_setstr(saveall_sets_format_item, sformat);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/*
 *  write a set or sets to a file
 */
static void do_saveall_sets_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char fn[512], *s;

    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR (cbs->value, charset, &s)) {
        errwin("Error converting XmString to char string in do_write_sets_proc()");
        return;
    }
    strcpy(fn, s);
    XtFree(s);

    strcpy(sformat, (char *) xv_getstr(saveall_sets_format_item));
    set_wait_cursor();
    do_writesets(maxgraph, -1, 1, fn, sformat);
    unset_wait_cursor();
}
