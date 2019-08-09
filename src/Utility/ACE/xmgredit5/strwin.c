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
 * strings, lines, and boxes
 *
 */

#ifndef lint
static char RCSid[] = "$Id: strwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

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
static Widget strings_depth_item;
static Widget *strings_sym_item;
static Widget strings_symsize_item;
static Widget *strings_format_item;
static Widget strings_frame_item;
static Widget *strings_fill_item;
static Widget *strings_fillcol_item;
static Widget *strings_prec_item;
static Widget *strings_grid_item;
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
void do_boxes_proc();
void do_lines_proc();
void do_move_proc();
void do_copy_object_proc(void)
{
}
void do_delete_object_proc();
void strings_loc_proc();
void strings_depth_proc();
void strings_ang_proc();
void strings_edit_proc();
void do_clear_lines();
void do_clear_boxes();
void do_clear_text();
void do_place_label();

void boxes_def_proc(void)
{
    box_color = (int) GetChoice(boxes_pen_item);
    box_loctype = (int) GetChoice(boxes_loc_item) ? VIEW : WORLD;
    box_lines = (int) GetChoice(boxes_lines_item) + 1;
    box_linew = (int) GetChoice(boxes_linew_item) + 1;
    switch (GetChoice(boxes_fill_item)) {
    case 0:
	box_fill = NONE;
	break;
    case 1:
	box_fill = COLOR;
	break;
    case 2:
	box_fill = PATTERN;
	break;
    }
    box_fillcolor = (int) GetChoice(boxes_fillcol_item);
    box_fillpat = (int) GetChoice(boxes_fillpat_item);
}

void lines_def_proc(void)
{
    Arg a;
    int value;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(lines_asize_item, &a, 1);
    line_asize = value / 50.0;
    line_color = (int) GetChoice(lines_pen_item);
    line_arrow = (int) GetChoice(lines_arrow_item);
    line_atype = (int) GetChoice(lines_atype_item);
    line_lines = (int) GetChoice(lines_style_item) + 1;
    line_linew = (int) GetChoice(lines_width_item) + 1;
    line_loctype = (int) GetChoice(lines_loc_item) ? VIEW : WORLD;
}

void updatestrings(void)
{
    Arg a;
    int iv;
    if (strings_frame) {
	SetChoice(strings_font_item, string_font);
	SetChoice(strings_pen_item, string_color);
	iv = (int) (100 * string_size);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(strings_size_item, &a, 1);
	XtSetArg(a, XmNvalue, string_rot);
	XtSetValues(strings_rot_item, &a, 1);
	SetChoice(strings_loc_item, string_loctype == VIEW ? 1 : 0);
	SetChoice(strings_just_item, string_just);

	XmToggleButtonSetState(strings_depth_item, string_type, False);
	SetChoice(strings_grid_item, string_grid == 0 ? 0 : 1);
	SetChoice(strings_prec_item, string_prec);
	SetChoice(strings_sym_item, string_sym);
	iv = (int) (100 * string_symsize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(strings_symsize_item, &a, 1);
    }
}

void update_lines(void)
{
    Arg a;
    int iv;
    if (lines_frame) {
	SetChoice(lines_pen_item, line_color);
	SetChoice(lines_style_item, line_lines - 1);
	SetChoice(lines_width_item, line_linew - 1);
	SetChoice(lines_arrow_item, line_arrow);
	SetChoice(lines_atype_item, line_atype);
	iv = (int) (50 * line_asize);
	XtSetArg(a, XmNvalue, iv);
	XtSetValues(lines_asize_item, &a, 1);
	SetChoice(lines_loc_item, line_loctype == VIEW ? 1 : 0);
    }
}

void update_boxes(void)
{
    Arg a;
    int iv;
    if (boxes_frame) {
	SetChoice(boxes_pen_item, box_color);
	SetChoice(boxes_lines_item, box_lines - 1);
	SetChoice(boxes_linew_item, box_linew - 1);
	switch (box_fill) {
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
	SetChoice(boxes_fillpat_item, box_fillpat);
	SetChoice(boxes_fillcol_item, box_fillcolor);
	SetChoice(boxes_loc_item, box_loctype == VIEW ? 1 : 0);
    }
}

void define_string_defaults(void)
{
    Arg a;
    int value;

    if (strings_frame) {
	string_font = (int) GetChoice(strings_font_item);
	string_color = (int) GetChoice(strings_pen_item);
	XtSetArg(a, XmNvalue, &value);
	XtGetValues(strings_size_item, &a, 1);
	string_size = value / 100.0;
	XtSetArg(a, XmNvalue, &value);
	XtGetValues(strings_rot_item, &a, 1);
	string_rot = value;
	string_loctype = (int) GetChoice(strings_loc_item) ? VIEW : WORLD;
	string_just = (int) GetChoice(strings_just_item);

	string_type = XmToggleButtonGetState(strings_depth_item);
	if (string_type) {
	    string_grid = GetChoice(strings_grid_item);
	    if (string_grid) {
		string_grid = MAXGRIDS;
	    }
	    string_prec = GetChoice(strings_prec_item);
	    string_sym = GetChoice(strings_sym_item);
	    XtSetArg(a, XmNvalue, &value);
	    XtGetValues(strings_symsize_item, &a, 1);
	    string_symsize = value / 100.0;
	}
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

    wbut = XtVaCreateManagedWidget("Text", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_loc_proc, 0);

    wbut = XtVaCreateManagedWidget("Depth text", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_depth_proc, 0);

    wbut = XtVaCreateManagedWidget("Text at angle", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_ang_proc, 0);

    wbut = XtVaCreateManagedWidget("Edit Text", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) strings_edit_proc, 0);

    wbut = XtVaCreateManagedWidget("Text props...", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_strings_popup, 0);

    wbut = XtVaCreateManagedWidget("Line", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_lines_proc, 0);

    wbut = XtVaCreateManagedWidget("Line props...", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_lines_popup, 0);

    wbut = XtVaCreateManagedWidget("Box", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_boxes_proc, 0);

    wbut = XtVaCreateManagedWidget("Box props...", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_boxes_popup, 0);

    wbut = XtVaCreateManagedWidget("Move object", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_move_proc, 0);

    wbut = XtVaCreateManagedWidget("Copy object", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_copy_object_proc, 0);

    wbut = XtVaCreateManagedWidget("Delete object", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_delete_object_proc, 0);

    wbut = XtVaCreateManagedWidget("Clear all text", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_text, 0);

    wbut = XtVaCreateManagedWidget("Clear all lines", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_lines, 0);

    wbut = XtVaCreateManagedWidget("Clear all boxes", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_clear_boxes, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, objects_panel,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, objects_frame);
    XtManageChild(objects_panel);
    XtManageChild(objects_frame);
}

void define_strings_popup(void)
{
    extern Widget app_shell;
    Widget wbut, rc, fr, rc2;
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

    strings_font_item = CreatePanelChoice1(strings_panel, "Font:",
					   11,
				"Times-Roman", "Times-Bold", "Times-Italic",
					   "Times-BoldItalic", "Helvetica",
				      "Helvetica-Bold", "Helvetica-Oblique",
				 "Helvetica-BoldOblique", "Greek", "Symbol",
					   0,
					   0);
    strings_pen_item = CreateColorChoice(strings_panel, "Color:", 1);
    strings_just_item = CreatePanelChoice1(strings_panel, "Justification:",
					   4,
					   "Left",
					   "Right",
					   "Centered",
					   0,
					   0);
    strings_loc_item = CreatePanelChoice1(strings_panel, "Position in:",
					  3,
					  "World coordinates",
					  "Viewport coordinates",
					  0,
					  0);

    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Rotation:", xmLabelGadgetClass, rc,
				     NULL);
    strings_rot_item = XtVaCreateManagedWidget("rotation", xmScaleWidgetClass, rc,
					       XmNminimum, 0,
					       XmNmaximum, 360,
					       XmNvalue, 0,
					       XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
					       XmNorientation, XmHORIZONTAL,
					       NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Size:", xmLabelGadgetClass, rc,
				     NULL);
    strings_size_item = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, rc,
						XmNminimum, 0,
						XmNmaximum, 400,
						XmNvalue, 100,
						XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
						XmNorientation, XmHORIZONTAL,
						NULL);
    XtManageChild(rc);

    fr = XmCreateFrame(strings_panel, "fr", NULL, 0);
    rc2 = XmCreateRowColumn(fr, "rc2", NULL, 0);
    strings_depth_item = XmCreateToggleButton(rc2, "Depth string", NULL, 0);
    XtManageChild(strings_depth_item);
    strings_grid_item = CreatePanelChoice1(rc2, "Attach to grid:",
					   3,
					   "Edit",
					   "Background",
					   0,
					   0);
    strings_sym_item = CreatePanelChoice1(rc2, "Symbol:",
					  4,
					  "None",
					  "Circle",
					  "Square",
					  0,
					  0);
    strings_prec_item = CreatePanelChoice1(rc2, "Precision:",
					   10,
				"0", "1", "2", "3", "4", "5", "6", "7", "8",
					   0,
					   0);

    rc = XmCreateRowColumn(rc2, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Symbol size:", xmLabelGadgetClass, rc,
				     NULL);
    strings_symsize_item = XtVaCreateManagedWidget("stringsize", xmScaleWidgetClass, rc,
						   XmNminimum, 0,
						   XmNmaximum, 400,
						   XmNvalue, 100,
						   XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
					       XmNorientation, XmHORIZONTAL,
						   NULL);
    XtManageChild(rc);
    XtManageChild(rc2);
    XtManageChild(fr);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, strings_panel, NULL);
    rc = XmCreateRowColumn(strings_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_string_defaults, 0);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc,
				   NULL);
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
    lines_width_item = (Widget *) CreatePanelChoice1(lines_panel, "Line width:",
						     10,
				"1", "2", "3", "4", "5", "6", "7", "8", "9",
						     NULL,
						     NULL);
    lines_style_item = (Widget *) CreatePanelChoice1(lines_panel, "Line style:",
						     6,
						     "Solid line",
						     "Dotted line",
						     "Dashed line",
						     "Long Dashed",
						     "Dot-dashed",
						     NULL,
						     NULL);

    lines_arrow_item = CreatePanelChoice1(lines_panel, "Arrow:",
					  5,
					  "None",
					  "At start",
					  "At end",
					  "Both ends",
					  0,
					  0);
    lines_atype_item = CreatePanelChoice1(lines_panel, "Arrow head type:",
					  4,
					  "Line",
					  "Filled",
					  "Hollow",
					  0,
					  0);

    rc = XmCreateRowColumn(lines_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wlabel = XtVaCreateManagedWidget("Arrow head size:", xmLabelGadgetClass, rc,
				     NULL);
    lines_asize_item = XtVaCreateManagedWidget("arrowsize", xmScaleWidgetClass, rc,
					       XmNminimum, 0,
					       XmNmaximum, 400,
					       XmNvalue, 100,
					       XmNshowValue, True,
				     XmNprocessingDirection, XmMAX_ON_RIGHT,
					       XmNorientation, XmHORIZONTAL,
					       NULL);
    XtManageChild(rc);

    lines_loc_item = CreatePanelChoice1(lines_panel, "Position in:",
					3,
					"World coordinates",
					"Viewport coordinates",
					0,
					0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, lines_panel, NULL);
    rc = XmCreateRowColumn(lines_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) lines_def_proc, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc,
				   NULL);
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
    boxes_linew_item = CreatePanelChoice1(boxes_panel, "Line width:",
					  10,
			     "1", "2", "3", "4", "5", "6", "7", "8", "9", 0,
					  NULL,
					  0);
    boxes_lines_item = (Widget *) CreatePanelChoice1(boxes_panel, "Line style:",
						     6,
						     "Solid line",
						     "Dotted line",
						     "Dashed line",
						     "Long Dashed",
						     "Dot-dashed",
						     NULL,
						     NULL);
    boxes_fill_item = CreatePanelChoice1(boxes_panel, "Fill:",
					 4,
					 "None",
					 "Color",
					 "Pattern",
					 NULL,
					 0);
    boxes_fillpat_item = CreatePanelChoice2(boxes_panel,
					    "Pattern:", 4, 17,
			   "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
				  "10", "11", "12", "13", "14", "15", 0, 0);
    boxes_fillcol_item = CreateColorChoice(boxes_panel, "Color:", 1);
    boxes_loc_item = CreatePanelChoice1(boxes_panel, "Position in:",
					3,
					"World coordinates",
					"Viewport coordinates",
					0,
					0);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, boxes_panel, NULL);
    rc = XmCreateRowColumn(boxes_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) boxes_def_proc, 0);

    wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, boxes_frame);
    XtManageChild(rc);

    update_boxes();
    XtManageChild(boxes_panel);
    XtManageChild(boxes_frame);
}
