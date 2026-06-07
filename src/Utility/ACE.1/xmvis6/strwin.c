/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/* $Id: strwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $
 *
 * strings, lines, and boxes
 *
 */

#include <stdio.h>
#include <math.h>
#include <sys/types.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

void define_objects_popup(void);
void define_strings_popup(void);
void define_lines_popup(void);
void define_boxes_popup(void);

static Widget objects_frame;
static Widget objects_panel;

static Widget strings_frame;
static Widget strings_panel;

static Widget lines_frame;
static Widget lines_panel;

static Widget boxes_frame;
static Widget boxes_panel;

static Widget string_item;
static Widget strings_item;
static Widget *strings_font_item;
static Widget strings_rot_item;
static Widget strings_size_item;
static Widget *strings_loc_item;
static Widget *strings_pen_item;
static Widget *strings_just_item;
Widget strings_x_item;
Widget strings_y_item;

static Widget *lines_arrow_item;
static Widget lines_asize_item;
static Widget *lines_atype_item;
static Widget *lines_pen_item;
static Widget *lines_style_item;
static Widget *lines_width_item;
static Widget *lines_loc_item;
static Widget *boxes_pen_item;
static Widget *boxes_lines_item;
static Widget *boxes_linew_item;
static Widget *boxes_fill_item;
static Widget *boxes_fillpat_item;
static Widget *boxes_fillcol_item;
static Widget *boxes_loc_item;

/* the following are defined in objutils.c */
void do_boxes_proc(void);
void do_lines_proc(void);
void do_move_proc(void);
void do_copy_object_proc(void);
void do_delete_object_proc(void);
void strings_loc_proc(void);
void strings_ang_proc(void);
void strings_edit_proc(void);
void do_clear_lines(void);
void do_clear_boxes(void);
void do_clear_text(void);

void boxes_def_proc(void)
{
    sysbox.color = (int) GetChoice(boxes_pen_item);
    sysbox.loctype = (int) GetChoice(boxes_loc_item) ? VIEW : WORLD;
    sysbox.lines = (int) GetChoice(boxes_lines_item) + 1;
    sysbox.linew = (int) GetChoice(boxes_linew_item) + 1;
    switch (GetChoice(boxes_fill_item)) {
    case 0:
	sysbox.fill = NONE;
	break;
    case 1:
	sysbox.fill = COLOR;
	break;
    case 2:
	sysbox.fill = PATTERN;
	break;
    }
    sysbox.fillcolor = (int) GetChoice(boxes_fillcol_item);
    sysbox.fillpattern = (int) GetChoice(boxes_fillpat_item);
}

void lines_def_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(lines_asize_item, &a, 1);
    sysline.asize = value / 50.0;
    sysline.color = (int) GetChoice(lines_pen_item);
    sysline.arrow = (int) GetChoice(lines_arrow_item);
    sysline.atype = (int) GetChoice(lines_atype_item);
    sysline.lines = (int) GetChoice(lines_style_item) + 1;
    sysline.linew = (int) GetChoice(lines_width_item) + 1;
    sysline.loctype = (int) GetChoice(lines_loc_item) ? VIEW : WORLD;
}

void updatestrings(void)
{
    Arg a;
    int iv;
    if (strings_frame) {
	SetChoice(strings_font_item, sysstr.font);
	SetChoice(strings_pen_item, sysstr.color);
	iv = (int) (100 * sysstr.charsize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(strings_size_item, &a, 1);
	XtSetArg(a, XmNvalue, sysstr.rot);
	XtSetValues(strings_rot_item, &a, 1);
	SetChoice(strings_loc_item, sysstr.loctype == VIEW ? 1 : 0);
	SetChoice(strings_just_item, sysstr.just);
    }
}

void update_lines(void)
{
    Arg a;
    int iv;
    if (lines_frame) {
	SetChoice(lines_pen_item, sysline.color);
	SetChoice(lines_style_item, sysline.lines - 1);
	SetChoice(lines_width_item, sysline.linew - 1);
	SetChoice(lines_arrow_item, sysline.arrow);
	SetChoice(lines_atype_item, sysline.atype);
	iv = (int) (50 * sysline.asize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(lines_asize_item, &a, 1);
	SetChoice(lines_loc_item, sysline.loctype == VIEW ? 1 : 0);
    }
}

void update_boxes(void)
{
    Arg a;
    int iv;
    if (boxes_frame) {
	SetChoice(boxes_pen_item, sysbox.color);
	SetChoice(boxes_lines_item, sysbox.lines - 1);
	SetChoice(boxes_linew_item, sysbox.linew - 1);
	switch (sysbox.fill) {
	case NONE:
	    SetChoice(boxes_fill_item, 0);
	    break;
	case COLOR:
	    SetChoice(boxes_fill_item, 1);
	    break;
	case PATTERN:
	    SetChoice(boxes_fill_item, 2);
	    break;
	}
	SetChoice(boxes_fillpat_item, sysbox.fillpattern);
	SetChoice(boxes_fillcol_item, sysbox.fillcolor);
	SetChoice(boxes_loc_item, sysbox.loctype == VIEW ? 1 : 0);
    }
}

void define_string_defaults(void)
{
    Arg a;
    int value;

    if (strings_frame) {
	sysstr.font = (int) GetChoice(strings_font_item);
	sysstr.color = (int) GetChoice(strings_pen_item);
	XtSetArg(a, XmNvalue, &value);
	XtGetValues(strings_size_item, &a, 1);
	sysstr.charsize = value / 100.0;
	XtSetArg(a, XmNvalue, &value);
	XtGetValues(strings_rot_item, &a, 1);
	sysstr.rot = value;
	sysstr.loctype = (int) GetChoice(strings_loc_item) ? VIEW : WORLD;
	sysstr.just = (int) GetChoice(strings_just_item);
    }
}

void define_objects_popup(void)
{
    extern Widget app_shell;
    Widget wbut;
    Widget wlabel;
    int x, y;

    if (objects_frame) {
	XtRaise(objects_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    objects_frame = XmCreateDialogShell(app_shell, "Objects", NULL, 0);
    handle_close(objects_frame);
    XtVaSetValues(objects_frame, XmNx, x, XmNy, y, NULL);
    objects_panel = XmCreateRowColumn(objects_frame, "ticks_rc", NULL, 0);

    wbut = XtVaCreateManagedWidget("Text", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_loc_proc, 0);

    wbut = XtVaCreateManagedWidget("Text at angle", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_ang_proc, 0);

    wbut = XtVaCreateManagedWidget("Edit Text", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_edit_proc, 0);

    wbut = XtVaCreateManagedWidget("Text props...", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_strings_popup, 0);

    wbut = XtVaCreateManagedWidget("Line", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_lines_proc, 0);

    wbut = XtVaCreateManagedWidget("Line props...", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_lines_popup, 0);

    wbut = XtVaCreateManagedWidget("Box", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_boxes_proc, 0);

    wbut = XtVaCreateManagedWidget("Box props...", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_boxes_popup, 0);

    wbut = XtVaCreateManagedWidget("Move object", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_move_proc, 0);

    wbut = XtVaCreateManagedWidget("Copy object", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_copy_object_proc, 0);

    wbut = XtVaCreateManagedWidget("Delete object", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_delete_object_proc, 0);

    wbut = XtVaCreateManagedWidget("Clear all text", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_text, 0);

    wbut = XtVaCreateManagedWidget("Clear all lines", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_lines, 0);

    wbut = XtVaCreateManagedWidget("Clear all boxes", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_boxes, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, objects_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, objects_frame);
    XtManageChild(objects_panel);
    XtManageChild(objects_frame);
}

void define_strings_popup(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    Widget wlabel;
    int x, y;

    if (strings_frame) {
	updatestrings();
	XtRaise(strings_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    strings_frame = XmCreateDialogShell(app_shell, "Strings", NULL, 0);
    handle_close(strings_frame);
    XtVaSetValues(strings_frame, XmNx, x, XmNy, y, NULL);
    strings_panel = XmCreateRowColumn(strings_frame, "strings_rc", NULL, 0);

    strings_font_item = CreatePanelChoice1(strings_panel, "Font:", 11, "Times-Roman", "Times-Bold", "Times-Italic", "Times-BoldItalic", "Helvetica", "Helvetica-Bold", "Helvetica-Oblique", "Helvetica-BoldOblique", "Greek", "Symbol", 0, 0);
    strings_pen_item = CreateColorChoice(strings_panel, "Color:", 1);
    strings_just_item = CreatePanelChoice1(strings_panel, "Justification:", 4, "Left", "Right", "Centered", 0, 0);
    strings_loc_item = CreatePanelChoice1(strings_panel, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);

    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Rotation:", xmLabelGadgetClass, rc, NULL);
    strings_rot_item = XtVaCreateManagedWidget("rotation", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 360, XmNvalue, 0, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Size:", xmLabelGadgetClass, rc, NULL);
    strings_size_item = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, strings_panel, NULL);
    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_string_defaults, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, strings_frame);
    XtManageChild(rc);

    updatestrings();
    XtManageChild(strings_panel);
    XtManageChild(strings_frame);
}

void define_lines_popup(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    Widget wlabel;
    int x, y;

    if (lines_frame) {
	update_lines();
	XtRaise(lines_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    lines_frame = XmCreateDialogShell(app_shell, "Lines", NULL, 0);
    handle_close(lines_frame);
    XtVaSetValues(lines_frame, XmNx, x, XmNy, y, NULL);
    lines_panel = XmCreateRowColumn(lines_frame, "lines_rc", NULL, 0);

    lines_pen_item = CreateColorChoice(lines_panel, "Color:", 1);
    lines_width_item = (Widget *) CreatePanelChoice1(lines_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", NULL, NULL);
    lines_style_item = (Widget *) CreatePanelChoice1(lines_panel, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);

    lines_arrow_item = CreatePanelChoice1(lines_panel, "Arrow:", 5, "None", "At start", "At end", "Both ends", 0, 0);
    lines_atype_item = CreatePanelChoice1(lines_panel, "Arrow head type:", 4, "Line", "Filled", "Hollow", 0, 0);

    rc = XmCreateRowColumn(lines_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Arrow head size:", xmLabelGadgetClass, rc, NULL);
    lines_asize_item = XtVaCreateManagedWidget("arrowsize", xmScaleWidgetClass, rc, XmNminimum, 0, XmNmaximum, 400, XmNvalue, 100, XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    lines_loc_item = CreatePanelChoice1(lines_panel, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, lines_panel, NULL);
    rc = XmCreateRowColumn(lines_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) lines_def_proc, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, lines_frame);
    XtManageChild(rc);

    update_lines();
    XtManageChild(lines_panel);
    XtManageChild(lines_frame);
}

void define_boxes_popup(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    Widget wlabel;
    int x, y;

    if (boxes_frame) {
	update_boxes();
	XtRaise(boxes_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    boxes_frame = XmCreateDialogShell(app_shell, "Boxes", NULL, 0);
    handle_close(boxes_frame);
    XtVaSetValues(boxes_frame, XmNx, x, XmNy, y, NULL);
    boxes_panel = XmCreateRowColumn(boxes_frame, "boxes_rc", NULL, 0);

    boxes_pen_item = CreateColorChoice(boxes_panel, "Color:", 1);
    boxes_linew_item = CreatePanelChoice1(boxes_panel, "Line width:", 10, "1", "2", "3", "4", "5", "6", "7", "8", "9", 0, NULL, 0);
    boxes_lines_item = (Widget *) CreatePanelChoice1(boxes_panel, "Line style:", 6, "Solid line", "Dotted line", "Dashed line", "Long Dashed", "Dot-dashed", NULL, NULL);
    boxes_fill_item = CreatePanelChoice1(boxes_panel, "Fill:", 4, "None", "Color", "Pattern", NULL, 0);
    boxes_fillpat_item = CreatePanelChoice2(boxes_panel, "Pattern:", 4, 17, "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 0, 0);
    boxes_fillcol_item = CreateColorChoice(boxes_panel, "Color:", 1);
    boxes_loc_item = CreatePanelChoice1(boxes_panel, "Position in:", 3, "World coordinates", "Viewport coordinates", 0, 0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, boxes_panel, NULL);
    rc = XmCreateRowColumn(boxes_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) boxes_def_proc, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, boxes_frame);
    XtManageChild(rc);

    update_boxes();
    XtManageChild(boxes_panel);
    XtManageChild(boxes_frame);
}
