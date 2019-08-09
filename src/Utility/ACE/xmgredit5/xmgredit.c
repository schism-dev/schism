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
 * xmgredit.c - entry point for Motif
 *
 */
#ifndef lint
static char RCSid[] = "$Id: xmgredit.c,v 1.6 2008/01/17 02:56:08 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "gredit_icon.h"

#define MENU_HELP		200
#define MENU_EXIT		201

extern int readimage;
extern char image_filename[];

extern int debuglevel;

extern Colormap cmap;

Pixmap zoompm, shrinkpm, expandpm, autopm, elevpm;
Pixmap concpm, flowpm, pathlpm, pointspm, elevmpm, gridpm;
Pixmap uppm, leftpm, downpm, rightpm, isolpm, bgridpm, bisolpm;
Pixmap nodenpm, elemnpm, depthnpm, node2pm;
Pixmap eboundpm, bboundpm, gridfpm, buildppm;

Widget app_shell;
Widget main_frame;
Widget menu_pane;
Widget menu_bar;
Widget canvas;
Widget loclab;
Widget arealab;
Widget perimlab;
Widget command;
Widget redraw_button;
XmString string;		/* string for current location */
XmString astring;		/* string for area */
XmString pstring;		/* string for perimeter */

Widget main_display_item[30];

XmString mstring;		/* string for mode indicator */
Widget write_mode_item;

XmString statstring;		/* string for mode indicator */
Widget statlab;

XmString lstring;		/* string for grid locator */
extern Widget locate_grid_item;

void my_proc();
void reset_world();

void do_grad();
void do_test_tri();
void do_newboundcb();
void do_select_area();
void do_select_region();
void do_clear_region();
void do_showelems_region();
void do_shownodes_region();
void do_restrict_region();
void do_select_peri();
void place_build();
void delete_build();
void delete_build_inside_region();
void delete_build_outside_region();
void move_build();
void loadgrid_build();
void mergegrid_build();
void mergebound_build();
void clear_build();
void toggle_build();
void define_extbound();
void move_extbound();
void insert_boundpt();
void define_intbound();
void move_intbound();
void del_intbound();
void compute_bound();
void clear_extbound();
void clear_intbound();
void line_bound();
void points_bound();
void numbers_bound();
void define_openclosed_proc(void);
void create_gwind_popup();
void create_grad_popup();
void create_display_popup();
void create_auto_popup();
void create_interp_frame(void);
void create_status_popup();
void create_isolines_popup();
void create_slice_popup();
void create_files_popup();
void create_mod_popup();
void create_edit_popup();
void create_build_popup();
void create_genfd_popup();
void create_prop_popup();
void create_delh_popup();
void create_pom_popup();
void create_nodelist_popup();
void create_buildprops_frame();
void gwindleft_proc();
void gwindright_proc();
void gwinddown_proc();
void gwindup_proc();
void gwindshrink_proc();
void gwindexpand_proc();
void do_select_point(Widget w, int cd);
void create_file_popup();
void create_wdata_frame();
void create_autob_popup();
void create_vel_frame();
void create_info_popup();
void create_mindist_popup();
void create_triang_popup();
void create_spreadr_popup();
void create_printer_setup();
void create_dimw_popup();
void create_locator_frame();
void create_gotoxy_frame();
void create_worlddb_popup();
void create_graph2d_popup();
void create_cedit_popup();
void define_objects_popup();
void create_convgrid_frame();
void create_inform_popup(void);
void do_load_centers();
void do_check_gridboundary_intersect(void);
void do_drawgrid();
void do_random_placement();
void write_mode_str(char *buf);
void invert_pixmap(Pixmap p, int flag);
void create_image_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_world_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_view_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_adjd_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_qual_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_vol_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_adcircbound_frame(void);
void do_read_region(void);
void do_save_region();
void save_region_proc(Widget w, XtPointer clientd, XtPointer calld);
void init_pm(void);

/*
 * xlib objects for drawing
 */
extern Display *disp;
extern int use_colors;

int go_locateflag = 0;		/* locator */

int cancel_flag = 0;

/*
 * for locator
 */
int deltaflag = 0;
extern int pointset;
extern double dsx, dsy;

/* used to set up XmStrings */
XmStringCharSet charset = (XmStringCharSet) XmSTRING_DEFAULT_CHARSET;

XtAppContext app_con;

Dimension xmgredit_w, xmgredit_h;

typedef struct {
    XmFontList f;
} ApplicationData, *ApplicationDataPtr;

static XtResource res[] = {
    {XmNfontList, XmCFontList, XmRFontList, sizeof(XmFontList), XtOffset(ApplicationDataPtr, f), XmRString, "9x15"}
};

/*
 * open display and get defaults
 */
void initialize_screen(int *argc, char **argv)
{
    int screen_width, screen_height;
    ApplicationData data;
    app_shell = XtAppInitialize(&app_con, "XMgredit", NULL, 0, argc, argv, NULL, NULL, 0);
    disp = XtDisplay(app_shell);
    if (!disp) {
	XtWarning("XMgredit: can't open display, exiting...");
	exit(1);
    }
/*
    XtGetApplicationResources(app_shell, &data, res, XtNumber(res), NULL, 0);
*/
    use_colors = DisplayPlanes(disp, DefaultScreen(disp));
    initialize_cms_data();
    screen_width = DisplayWidth(disp, DefaultScreen(disp));
    screen_height = DisplayHeight(disp, DefaultScreen(disp));
    if (xmgredit_w == 0) {
	xmgredit_w = 3 * screen_width / 4;
    }
    if (xmgredit_h == 0) {
	xmgredit_h = 3 * screen_height / 4;
    }
}

void do_isolines_popup(Widget w, XtPointer cd, XtPointer cld)
{
    switch ((long) cd) {
    case 0:
	create_isolines_popup("Edit grid isolines", &grid[curgrid].ip);
	break;
    case 1:
	create_isolines_popup("Background grid isolines", &grid[backgrid].ip);
	break;
    }
}

/*
 * called from draw routine to set Window
*/
void set_canvas(void)
{
    set_window(canvas);
}

/*
 * cancel button on main window
 */
Widget abort_button;
Window abort_win;

void setup_abort(void)
{
    Arg al;
    XtSetArg(al, XmNsensitive, True);
    XtSetValues(abort_button, &al, 1);
    XmUpdateDisplay(abort_button);
    abort_win = XtWindow(abort_button);
}

void unsetup_abort(void)
{
    Arg al;
    XtSetArg(al, XmNsensitive, False);
    XtSetValues(abort_button, &al, 1);
    XmUpdateDisplay(abort_button);
}

void do_abort(void)
{
    cancel_flag = 1;
    unsetup_abort();
}

void set_abort_button(void)
{
    Arg al;
    XtSetArg(al, XmNlabelString, XmStringCreateLtoR("Cancel", charset));
    XtSetValues(redraw_button, &al, 1);
    XmUpdateDisplay(redraw_button);
}

/*
 * reset label on redraw button
 */
void set_redraw_button(void)
{
    Arg al;
    XtSetArg(al, XmNlabelString, XmStringCreateLtoR("Redraw", charset));
    XtSetValues(redraw_button, &al, 1);
    XmUpdateDisplay(redraw_button);
}


/*
 * restrict drawing to inside the current region
 */
int restrict_to_region;

void do_restrict_region(void)
{
    restrict_to_region = restrict_to_region ? 0 : 1;
}

/*
 * callback for exit and help buttons
 */
void MenuCB(Widget w, caddr_t client_data, caddr_t call_data)
{
    switch ((long) client_data) {
	case MENU_HELP:
	break;
    case MENU_EXIT:
	if (yesno("Exit xmgredit? Are you sure?", "Press Yes or No", "Yes", "No")) {
	    exit(0);
	}
	break;
    default:
	break;
    }

}

/*
 * set sensitive for menu items
 */
Widget bound_b[11];
Widget build_b[8], advf;

void update_fuzz_items(void)
{
    if (!bound_b[0]) {
	return;
    }
    if (ebound_defined) {
	XtVaSetValues(bound_b[0], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[2], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[9], XmNsensitive, True, NULL);
    } else {
	XtVaSetValues(bound_b[0], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[2], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[9], XmNsensitive, False, NULL);
    }
    if (ibound_defined) {
	XtVaSetValues(bound_b[1], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[3], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[4], XmNsensitive, True, NULL);
    } else {
	XtVaSetValues(bound_b[1], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[3], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[4], XmNsensitive, False, NULL);
    }
    if (!ebound_defined && !ibound_defined) {
	XtVaSetValues(bound_b[5], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[6], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[7], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[8], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[9], XmNsensitive, False, NULL);
	XtVaSetValues(bound_b[10], XmNsensitive, False, NULL);
    } else {
	XtVaSetValues(bound_b[5], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[6], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[7], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[8], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[9], XmNsensitive, True, NULL);
	XtVaSetValues(bound_b[10], XmNsensitive, True, NULL);
    }
    if (build_defined) {
	XtVaSetValues(build_b[0], XmNsensitive, True, NULL);
	XtVaSetValues(build_b[1], XmNsensitive, True, NULL);
	XtVaSetValues(build_b[3], XmNsensitive, True, NULL);
	XtVaSetValues(build_b[5], XmNsensitive, True, NULL);
/*
	XtVaSetValues(advf, XmNsensitive, True, NULL);
*/
    } else {
	XtVaSetValues(build_b[0], XmNsensitive, False, NULL);
	XtVaSetValues(build_b[1], XmNsensitive, False, NULL);
	XtVaSetValues(build_b[3], XmNsensitive, False, NULL);
	XtVaSetValues(build_b[5], XmNsensitive, False, NULL);
/*
	XtVaSetValues(advf, XmNsensitive, False, NULL);
*/
    }
    if (ebound_defined || ibound_defined) {
	XtVaSetValues(build_b[2], XmNsensitive, True, NULL);
    } else {
	XtVaSetValues(build_b[2], XmNsensitive, False, NULL);
    }
}

/*
 * callback for springs selection
 */
void do_springs(void)
{
	/*
    if (!grid_is_triangles(curgrid)) {
	errwin("Grid must be compose entirely of linear elements");
	return;
    }
    */
    compute_bound();
    write_mode_str("Adjusting nodes...");
    adjustnode(curgrid);
}

/*
 * callback for springs inside a region selection
 */
void do_region_springs(void)
{
	/*
    if (!grid_is_triangles(curgrid)) {
	errwin("Grid must be compose entirely of linear elements");
	return;
    }
    */
    compute_bound();
    write_mode_str("Adjusting nodes...");
    adjustnode_region(curgrid);
}

/*
    write_mode_str("Press the left mouse button to display the depth of the background grid at the pointer");
    write_mode_str("Press the left mouse button to display the depth of the edit grid at the pointer");
    write_mode_str("Use the left mouse button to display the element nearest the pointer");
    write_mode_str("Use the left mouse button to display the node nearest the pointer");
    write_mode_str("Use the left mouse button to display the element containing the pointer");
*/

void open_locate1_popup(void)
{
    set_action(0);
    set_action(GET_NEAREST_NODE);
}

void do_find_build(void)
{
    set_action(0);
    set_action(GET_NEAREST_BUILDPT);
}

void do_find_alldepth(void)
{
    set_action(0);
    set_action(GET_DEPTH_ALL);
}

void open_locate2_popup(void)
{
    set_action(0);
    set_action(GET_ELEMENT);
}

void open_locate3_popup(void)
{
    set_action(0);
    set_action(GET_NEAREST_ELEMENT);
}

void open_locate4_popup(void)
{
    set_action(0);
    curgrid = 0;
    set_action(GET_GRID_DEPTH);
}

void open_locate5_popup(void)
{
    set_action(0);
    set_action(GET_BACK_DEPTH);
}

void open_locate6_popup(void)
{
    write_mode_str("Goto node/element");
    create_goto_popup();
}

/*
 * zoom callback
 */
void do_zoom(void)
{
    set_action(0);
    set_action(ZOOM_1ST);
}

/*
 * bandwidth reduction callback
 */
void do_reducebw(void)
{
    int newb, oldb;
    char buf[256];
    set_wait_cursor(NULL);
    renum(curgrid, &oldb, &newb);
    unset_wait_cursor(NULL);
    if (oldb != -1) {
	sprintf(buf, "Old bandwidth = %d, New bandwidth = %d", oldb, newb);
	write_mode_str(buf);
    } else {
	write_mode_str("Grid bandwidth already optimized");
    }
}

void define_openclosed_proc(void)
{
}

/*
 * create the main menubar
 */
static Widget CreateMenuBar(Widget parent)
{
    Widget menu_bar;
    Widget cascade;
    Widget button;
    Arg al[10];
    int ac;

    XtSetArg(al[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);
    menu_bar = XmCreateMenuBar(parent, "menu_bar", NULL, 0);
    menu_pane = XmCreatePulldownMenu(menu_bar, "Files menu", al, 1);

    button = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_file_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_wdata_frame, 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Print...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_printer_setup, 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Exit", xmPushButtonWidgetClass, menu_pane,
		      XmNacceleratorText, XmStringCreateLtoR("F3", charset),
				     XmNaccelerator, "<Key>F3:",
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_EXIT);

    cascade = XtVaCreateManagedWidget("File", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'F',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Features", al, 1);

    button = XtVaCreateManagedWidget("Settings...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gwind_popup, 0);

/*
    button = XtVaCreateManagedWidget("Features...", xmPushButtonWidgetClass, menu_pane,
			  NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_display_popup, 0);
*/

    button = XtVaCreateManagedWidget("Autoscaling...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_auto_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Map scale...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_vel_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Isolines of bathymetry (Edit grid)...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_isolines_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Isolines of bathymetry (Background grid)...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_isolines_popup, (XtPointer) 1);
    button = XtVaCreateManagedWidget("World...", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_world_frame, (XtPointer) 0);
    button = XtVaCreateManagedWidget("Viewport...", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_view_frame, (XtPointer) 0);
    button = XtVaCreateManagedWidget("Read image...", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_image_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Colors...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_cedit_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Text...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_objects_popup, (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("Display", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'D',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Locate menu pane", al, 1);

    button = XtVaCreateManagedWidget("Find nearest node", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate1_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find element", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate2_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find nearest element", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate3_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find depth in edit grid", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate4_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find depth in background", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate5_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find nearest build point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_find_build, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Find depth in all", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_find_alldepth, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Goto node/element...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate6_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Goto X, Y...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gotoxy_frame, (XtPointer) NULL);

/*
    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
			  NULL);

    button = XtVaCreateManagedWidget("Distance between 2 points", xmPushButtonWidgetClass, menu_pane,
			  NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate6_popup, (XtPointer) NULL);
*/

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Set fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_point, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Clear fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_point, (XtPointer) 1);

    button = XtVaCreateManagedWidget("Locator props...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_locator_frame, (XtPointer) NULL);

    cascade = XtVaCreateManagedWidget("Locate", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'L',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Edit menu pane", al, 1);

    button = XtVaCreateManagedWidget("Edit triangles...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_mod_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Edit over grid/regions...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_edit_popup, (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("Edit", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'E',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Regions", al, 1);

    button = XtVaCreateManagedWidget("Create region", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Clear region", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_clear_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Read region...", xmPushButtonWidgetClass, menu_pane,
                                     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_read_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Save region...", xmPushButtonWidgetClass, menu_pane,
                                     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_save_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Show selected elements", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_showelems_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Show selected nodes", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_shownodes_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Restrict drawing to region", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_restrict_region, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Area", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_area, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Perimeter", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_peri, (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("Regions", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'R',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Boundary menu pane", al, 1);

    button = XtVaCreateManagedWidget("Create external boundary", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_extbound, (XtPointer) NULL);

    bound_b[0] = XtVaCreateManagedWidget("Move external boundary point", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[0], XmNactivateCallback, (XtCallbackProc) move_extbound, (XtPointer) NULL);

    bound_b[8] = XtVaCreateManagedWidget("Insert boundary points", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[8], XmNactivateCallback, (XtCallbackProc) insert_boundpt, (XtPointer) NULL);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    button = XtVaCreateManagedWidget("Create internal boundary", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_intbound, (XtPointer) NULL);

    bound_b[1] = XtVaCreateManagedWidget("Move internal boundary point", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[1], XmNactivateCallback, (XtCallbackProc) move_intbound, (XtPointer) NULL);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    button = XtVaCreateManagedWidget("Compute boundary", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) compute_bound, (XtPointer) NULL);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    bound_b[2] = XtVaCreateManagedWidget("Clear exterior boundary", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[2], XmNactivateCallback, (XtCallbackProc) clear_extbound, (XtPointer) NULL);

    bound_b[3] = XtVaCreateManagedWidget("Delete internal boundary", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[3], XmNactivateCallback, (XtCallbackProc) del_intbound, (XtPointer) NULL);

    bound_b[4] = XtVaCreateManagedWidget("Clear all interior boundaries", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[4], XmNactivateCallback, (XtCallbackProc) clear_intbound, (XtPointer) NULL);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    bound_b[5] = XtVaCreateManagedWidget("Boundary as line", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[5], XmNactivateCallback, (XtCallbackProc) line_bound, (XtPointer) NULL);

    bound_b[6] = XtVaCreateManagedWidget("Boundary as points", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[6], XmNactivateCallback, (XtCallbackProc) points_bound, (XtPointer) NULL);

    bound_b[7] = XtVaCreateManagedWidget("Boundary as numbers", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[7], XmNactivateCallback, (XtCallbackProc) numbers_bound, (XtPointer) NULL);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);
    bound_b[10] = XtVaCreateManagedWidget("Check boundary",
					  xmPushButtonWidgetClass, menu_pane,
					  XmNsensitive, False,
					  NULL);
    XtAddCallback(bound_b[10], XmNactivateCallback, (XtCallbackProc) do_check_gridboundary_intersect, (XtPointer) NULL);
    bound_b[9] = XtVaCreateManagedWidget("*Define open/closed boundaries",
					 xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(bound_b[9], XmNactivateCallback, (XtCallbackProc) define_openclosed_proc, (XtPointer) NULL);

    cascade = XtVaCreateManagedWidget("Boundaries", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'B',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Build menu pane", al, 1);

    button = XtVaCreateManagedWidget("Place build points",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) place_build,
		  (XtPointer) 0);

    build_b[0] = XtVaCreateManagedWidget("Delete build points inside region",
					 xmPushButtonWidgetClass, menu_pane,
					 NULL);
    XtAddCallback(build_b[0], XmNactivateCallback, (XtCallbackProc) delete_build_inside_region, (XtPointer) 0);

    build_b[0] = XtVaCreateManagedWidget("Delete build points outside region",
					 xmPushButtonWidgetClass, menu_pane,
					 NULL);
    XtAddCallback(build_b[0], XmNactivateCallback, (XtCallbackProc) delete_build_outside_region, (XtPointer) 0);

    build_b[0] = XtVaCreateManagedWidget("Delete build points",
					 xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(build_b[0], XmNactivateCallback, (XtCallbackProc) delete_build, (XtPointer) 0);

    build_b[1] = XtVaCreateManagedWidget("Move build points",
					 xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(build_b[1], XmNactivateCallback, (XtCallbackProc) move_build, (XtPointer) 0);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    button = XtVaCreateManagedWidget("Circular spread...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_build_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Rectangular spread...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_spreadr_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Finite difference grids...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_genfd_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Automatic placement...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_autob_popup, (XtPointer) 0);

/*    button = XtVaCreateManagedWidget("*Random placement",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_random_placement, (XtPointer) 0);
*/

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    button = XtVaCreateManagedWidget("Load grid to build points",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) loadgrid_build, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Merge grid to build points",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) mergegrid_build, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Load centers", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_load_centers, (XtPointer) 0);

    build_b[2] = XtVaCreateManagedWidget("Merge boundary to build points", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(build_b[2], XmNactivateCallback, (XtCallbackProc) mergebound_build, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    build_b[3] = XtVaCreateManagedWidget("Triangulate build points...", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(build_b[3], XmNactivateCallback, (XtCallbackProc) create_triang_popup, (XtPointer) 0);


    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    button = XtVaCreateManagedWidget("Set minimum distance...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_mindist_popup, (XtPointer) 0);

    build_b[6] = XtVaCreateManagedWidget("Springs", xmPushButtonWidgetClass, menu_pane,
					 NULL);
    XtAddCallback(build_b[6], XmNactivateCallback, (XtCallbackProc) do_springs, (XtPointer) 0);

    build_b[7] = XtVaCreateManagedWidget("Springs in region", xmPushButtonWidgetClass, menu_pane,
					 NULL);
    XtAddCallback(build_b[7], XmNactivateCallback, (XtCallbackProc) do_region_springs, (XtPointer) 0);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);

    build_b[5] = XtVaCreateManagedWidget("Clear current set of build points", xmPushButtonWidgetClass, menu_pane,
					 XmNsensitive, False,
					 NULL);
    XtAddCallback(build_b[5], XmNactivateCallback, (XtCallbackProc) clear_build, (XtPointer) 0);

    (void) XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				   NULL);
    button = XtVaCreateManagedWidget("Display properties...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_buildprops_frame,
		  (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("Build", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'u',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Modeling menu pane", al, 1);

    button = XtVaCreateManagedWidget("Load bathymetry", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_interp_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Create open/land boundaries...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_adcircbound_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane, NULL);

    button = XtVaCreateManagedWidget("Gradients & Slopes...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_grad_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane, NULL);

    button = XtVaCreateManagedWidget("Profiles...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_slice_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Reduce bandwidth (RENU2)", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_reducebw, (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("gridDEM", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'g',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "special_pane", al, 1);

    button = XtVaCreateManagedWidget("Properties...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_prop_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Select nodes...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_nodelist_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Dimensionless numbers...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_dimw_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Quality checks...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_qual_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Volume/depth check...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_vol_frame, (XtPointer) 0);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    cascade = XtVaCreateManagedWidget("Special", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'S',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "datab_pane", al, 1);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Status menu pane", NULL, 0);

    button = XtVaCreateManagedWidget("Status...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_status_popup, (XtPointer) 0);

    button = XtVaCreateManagedWidget("Info...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_info_popup, (XtPointer) 0);

    cascade = XtVaCreateManagedWidget("Status", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 't',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Help menu pane", NULL, 0);
    cascade = XtVaCreateManagedWidget("Help", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'H',
				      NULL);
    XtAddCallback(cascade, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_HELP);

    button = XtVaCreateManagedWidget("Help...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_info_popup, (XtPointer) 0);
    button = XtVaCreateManagedWidget("About...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_info_popup, (XtPointer) 0);

    ac = 0;
    XtSetArg(al[ac], XmNmenuHelpWidget, cascade);
    ac++;
    XtSetValues(menu_bar, al, ac);
    return (menu_bar);
}

/*
 * refresh callback
 */
void refresh(Widget w, caddr_t clid, caddr_t calld)
{
    Arg args;
    Dimension wh, ww;
    static int inc = 0;
    set_window(w);
    if (debuglevel == 1) {
	printf("In refresh, inwin = %d\n", inwin);
    }
    if (!inc) {
	initgraphics(0);
	set_defaults();
	set_up_world();
	do_drawgrid();
	inc++;
	inwin = 1;
    } else {
	refresh_from_backpix();
    }
}

/*
 * resize callback
 */
void resize(Widget w, caddr_t clid, caddr_t calld)
{
    Arg args;
    Dimension wh, ww;
    if (debuglevel == 1) {
	printf("In resize, inwin = %d\n", inwin);
    }
    if (inwin) {
	set_window(w);
	resize_backpix();
	initgraphics(0);
	set_up_world();
	do_drawgrid();
    }
}

/*
 * locator on main_panel
 */
void getpoints(int x, int y)
{
    double wx, wy, xtmp, ytmp, xconv(), yconv();
    char buf[256];
    Arg al;

    get_world(x, y, &wx, &wy);
    if (!go_locateflag) {
	return;
    }
    if (deltaflag) {
	switch (deltaflag) {
	case 1:
	    sprintf(buf, "DX, DY = [%.8g, %.8g]", wx - dsx, wy - dsy);
	    break;
	case 2:
	    sprintf(buf, "Distance = [%.8g]", hypot(dsx - wx, dsy - wy));
	    break;
	case 3:
	    if (dsx - wx != 0.0 && dsy - wy != 0.0) {
		sprintf(buf, "R, theta = [%.8g, %.8g]", hypot(dsx - wx, dsy - wy), 180.0 + 180.0 / M_PI * atan2(dsy - wy, dsx - wx));
	    }
	    break;
	case 4:
	    xtmp = xconv(wx);
	    ytmp = yconv(wy);
	    sprintf(buf, "VX, VY = [%.2lf, %.2lf]", xtmp, ytmp);
	    break;
	case 5:
	    sprintf(buf, "SX, SY = [%d, %d]\n", x, y);
	    break;

	}
    } else {
	sprintf(buf, "X, Y = [%.8g, %.8g]", wx, wy);
    }
    XmStringFree(string);
    string = XmStringCreateLtoR(buf, charset);
    XtSetArg(al, XmNlabelString, string);
    XtSetValues(loclab, &al, 1);
}

void panel_setmsgstr_value(Widget w, char *buf)
{
    Arg al;
    XmStringFree(string);
    string = XmStringCreateLtoR(buf, charset);
    XtSetArg(al, XmNlabelString, string);
    XtSetValues(loclab, &al, 1);
}

/*
 * informative message at the lower left of the canvas
 */
void write_mode_str(char *buf)
{
    Arg al;
    char tmpb[256];
    if (buf == NULL) {
	strcpy(tmpb, "Idle...");
    } else {
	strcpy(tmpb, buf);
    }
    XmStringFree(mstring);
    mstring = XmStringCreateLtoR(tmpb, charset);
    XtSetArg(al, XmNlabelString, mstring);
    XtSetValues(write_mode_item, &al, 1);
    XmUpdateDisplay(write_mode_item);
}

void do_select_point(Widget w, int cd)
{
    if (cd) {
	pointset = FALSE;
	dsx = dsy = 0.0;
	deltaflag = 0;
    } else {
	set_action(0);
	set_action(SEL_POINT);
	pointset = TRUE;
    }
}

void setbgcolor(int c)
{
    Arg al;
    XtSetArg(al, XmNbackground, c);
    XtSetValues(canvas, &al, 1);
}

/*
#define EDIT_GRID 0
#define BACKGROUND_GRID 1
#define BACKGROUND_BOUNDARY 2
#define EDIT_BOUNDARY 3
#define EDIT_GRID_NODE_NUMBERS 4
#define EDIT_GRID_ELEMENT_NUMBERS 5
#define EDIT_GRID_FILLED 6
#define EDIT_GRID_ISOLINES 7
#define BACKGROUND_GRID_ISOLINES 8
#define BUILD_POINTS 9
#define EDIT_GRID_DEPTHS 10
#define EDIT_GRID_NODES 11

display_flags[EDIT_GRID]
display_flags[EDIT_GRID_ISOLINES]
display_flags[EDIT_BOUNDARY]
display_flags[EDIT_GRID_FILLED]
display_flags[EDIT_GRID_NODE_NUMBERS]
display_flags[EDIT_GRID_ELEMENT_NUMBERS]
display_flags[EDIT_GRID_DEPTHS]
display_flags[BACKGROUND_GRID]
display_flags[BACKGROUND_GRID_ISOLINES]
display_flags[BACKGROUND_BOUNDARY]
display_flags[BUILD_POINTS]
display_flags[EDIT_GRID_NODES]
*/

Pixmap pm[12], savepm[12];

void togglegrid_proc(Widget w, int cd, XmAnyCallbackStruct * cbs)
{
    display_flags[cd] = XmToggleButtonGetState(main_display_item[cd]);
    if (display_flags[cd]) {
	XtVaSetValues(main_display_item[cd], XmNlabelPixmap, savepm[cd], NULL);
    } else {
	XtVaSetValues(main_display_item[cd], XmNlabelPixmap, pm[cd], NULL);
    }
    update_display_items();
}

/*
do_help()
{
printf("Got some help!\n");
}
*/

/*
 * create main window, menubar, initialize canvas, etc.
 * then call XtAppMainLoop();
 */
void do_main_loop(void)
{
    Widget bt, bt2, fr, fr2, fr3, rc, rc2, rc3, bb, form, rctmp;
    Arg al[10];
    int ac;
    XSetWindowAttributes sw;
    Pixmap icon;
    XEvent event;


    ac = 0;
    XtSetArg(al[ac], XmNshadowThickness, 0);
    ac++;
    XtSetArg(al[ac], XmNwidth, xmgredit_w);
    ac++;
    XtSetArg(al[ac], XmNheight, xmgredit_h + 100);
    ac++;
    main_frame = XmCreateMainWindow(app_shell, "main", al, ac);

    menu_bar = CreateMenuBar(main_frame);
    XtManageChild(menu_bar);

    form = XmCreateForm(main_frame, "form", NULL, 0);

    fr = XtVaCreateManagedWidget("fr", xmFrameWidgetClass, form,
				 NULL);
    rc = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, fr,
				 XmNorientation, XmVERTICAL,
				 XmNpacking, XmPACK_TIGHT,
				 XmNspacing, 0,
				 XmNentryBorder, 0,
				 XmNmarginWidth, 0,
				 XmNmarginHeight, 0,
				 NULL);

    canvas = XtVaCreateManagedWidget("canvas", xmDrawingAreaWidgetClass, form,
				     XmNwidth, xmgredit_w,
				     XmNheight, xmgredit_h,
				   XmNnavigationType, XmEXCLUSIVE_TAB_GROUP,
				     XmNresizeWidth, False,
				     XmNresizeHeight, False,

			    XmNbackground, WhitePixel(XtDisplay(main_frame),
				      DefaultScreen(XtDisplay(main_frame))),
				     NULL);

    XtAddCallback(canvas, XmNexposeCallback, (XtCallbackProc) refresh, (XtPointer) NULL);
    XtAddCallback(canvas, XmNresizeCallback, (XtCallbackProc) resize, (XtPointer) NULL);
    XtAddEventHandler(canvas,
		      EnterWindowMask |
		      LeaveWindowMask |
		      ButtonPressMask |
		      ButtonReleaseMask |
		      PointerMotionMask |
		      Button1MotionMask |
		      KeyPressMask,
		      FALSE, (XtEventHandler) my_proc, NULL);

/*
    statstring = XmStringCreateLtoR("Idle...                                                                           ", charset);
    statlab = XtVaCreateManagedWidget("statlab", xmLabelGadgetClass, form,
                          XmNlabelString, statstring,
                          XmNalignment, XmALIGNMENT_BEGINNING,
                          NULL);
*/

    fr2 = XtVaCreateManagedWidget("fr2", xmFrameWidgetClass, form, NULL);
    rc2 = XmCreateForm(fr2, "form", NULL, 0);
    write_mode_item = XtVaCreateManagedWidget("modelabel", xmLabelGadgetClass, rc2,
					      XmNlabelString, mstring = XmStringCreateLtoR("Idle...                                       ", charset),
					XmNalignment, XmALIGNMENT_BEGINNING,
					      XmNrecomputeSize, True,
					      NULL);

    string = XmStringCreateLtoR("X, Y = [,]", charset);
    loclab = XtVaCreateManagedWidget("label Locate", xmLabelGadgetClass, rc2,
				     XmNlabelString, string,
				     XmNalignment, XmALIGNMENT_END,
				     XmNrecomputeSize, True,
				     NULL);
    XtVaSetValues(write_mode_item,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(loclab,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);
    XtManageChild(rc2);

    XtVaSetValues(fr,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);

    XtVaSetValues(canvas,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNleftAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_WIDGET,
		  XmNleftWidget, fr,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);
    init_pm();
/*
    make_buttons(rc);
*/
    rctmp = rc;

    XtVaSetValues(fr2,
		  XmNbottomAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtManageChild(form);

    redraw_button = XtVaCreateManagedWidget("Draw", xmPushButtonWidgetClass, rctmp,
					    NULL);
    XtAddCallback(redraw_button, XmNactivateCallback, (XtCallbackProc) do_drawgrid, (XtPointer) 0);
    /* XtAddCallback(bt, XmNhelpCallback, do_help, 0); */

/* zoom and autoscale */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    bt = XtVaCreateManagedWidget("Zoom", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, zoompm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_zoom, (XtPointer) NULL);
    bt = XtVaCreateManagedWidget("AS", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, autopm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) reset_world, (XtPointer) NULL);

/* expand/shrink */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    bt = XtVaCreateManagedWidget("Z", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, expandpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindshrink_proc, (XtPointer) NULL);
    bt = XtVaCreateManagedWidget("z", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, shrinkpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindexpand_proc, (XtPointer) NULL);

/*
 * scrolling buttons
 */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    bt = XtVaCreateManagedWidget("Left", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, leftpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindright_proc, (XtPointer) NULL);
    bt = XtVaCreateManagedWidget("Right", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, rightpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindleft_proc, (XtPointer) NULL);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);

    bt = XtVaCreateManagedWidget("Down", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, downpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindup_proc, (XtPointer) NULL);
    bt2 = XtVaCreateManagedWidget("Up", xmPushButtonWidgetClass, rc3,
				  NULL);
    XtVaSetValues(bt2,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, uppm,
		  NULL);
    XtAddCallback(bt2, XmNactivateCallback, (XtCallbackProc) gwinddown_proc, (XtPointer) NULL);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);

/* display features */
    XtVaCreateManagedWidget("EdGr", xmLabelGadgetClass, rctmp, NULL);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);

    main_display_item[EDIT_GRID] = XtVaCreateManagedWidget("grid", xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						     XmNlabelPixmap, gridpm,
						     XmNselectPixmap, gridpm,
						     XmNfillOnSelect, False,
							   XmNspacing, 0,
						      XmNshadowThickness, 3,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							   XmNentryBorder, 0,
							   XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							   XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							   XmNmarginRight, 0,
							   NULL);
    XtAddCallback(main_display_item[EDIT_GRID],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID);
    display_flags[EDIT_GRID] = 1;
    XtVaSetValues(main_display_item[EDIT_GRID], XmNlabelPixmap, savepm[EDIT_GRID], NULL);

    main_display_item[EDIT_GRID_ISOLINES] = XtVaCreateManagedWidget("isol",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						     XmNlabelPixmap, isolpm,
						     XmNselectPixmap, isolpm,
						     XmNfillOnSelect, False,
							      XmNspacing, 0,
						      XmNshadowThickness, 3,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								    NULL);
    XtAddCallback(main_display_item[EDIT_GRID_ISOLINES],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_ISOLINES);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginTop, 0,
				  XmNmarginBottom, 0,
				  XmNmarginHeight, 0,
				  XmNmarginWidth, 0,
				  XmNmarginLeft, 0,
				  XmNmarginRight, 0,
				  NULL);
    main_display_item[EDIT_BOUNDARY] = XtVaCreateManagedWidget("editbound",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						   XmNlabelPixmap, eboundpm,
						     XmNselectPixmap, eboundpm,
						     XmNfillOnSelect, False,
							       XmNspacing, 0,
						      XmNshadowThickness, 3,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
							       NULL);
    XtAddCallback(main_display_item[EDIT_BOUNDARY], XmNvalueChangedCallback,
	       (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_BOUNDARY);
    display_flags[EDIT_BOUNDARY] = 1;
    XtVaSetValues(main_display_item[EDIT_BOUNDARY], XmNlabelPixmap, savepm[EDIT_BOUNDARY], NULL);

    main_display_item[EDIT_GRID_FILLED] = XtVaCreateManagedWidget("gridfilled",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, gridfpm,
						     XmNselectPixmap, gridfpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								  NULL);
    XtAddCallback(main_display_item[EDIT_GRID_FILLED],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_FILLED);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNpacking, XmPACK_COLUMN,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    main_display_item[EDIT_GRID_NODE_NUMBERS] = XtVaCreateManagedWidget("noden",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, nodenpm,
						     XmNselectPixmap, nodenpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								      NULL);
    XtAddCallback(main_display_item[EDIT_GRID_NODE_NUMBERS],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_NODE_NUMBERS);
    main_display_item[EDIT_GRID_ELEMENT_NUMBERS] = XtVaCreateManagedWidget("elevn",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, elemnpm,
						     XmNselectPixmap, elemnpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								      NULL);
    XtAddCallback(main_display_item[EDIT_GRID_ELEMENT_NUMBERS],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_ELEMENT_NUMBERS);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    main_display_item[EDIT_GRID_DEPTHS] = XtVaCreateManagedWidget("depthn",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						   XmNlabelPixmap, depthnpm,
						     XmNselectPixmap, depthnpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								  NULL);
    XtAddCallback(main_display_item[EDIT_GRID_DEPTHS],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_DEPTHS);
    main_display_item[EDIT_GRID_NODES] = XtVaCreateManagedWidget("noden2",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, node2pm,
						     XmNselectPixmap, node2pm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								 NULL);
    XtAddCallback(main_display_item[EDIT_GRID_NODES],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) EDIT_GRID_NODES);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);

    XtVaCreateManagedWidget("BkGr", xmLabelGadgetClass, rctmp, NULL);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    main_display_item[BACKGROUND_GRID] = XtVaCreateManagedWidget("backgrid",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, bgridpm,
						     XmNselectPixmap, bgridpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								 NULL);
    XtAddCallback(main_display_item[BACKGROUND_GRID],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) BACKGROUND_GRID);
    main_display_item[BACKGROUND_GRID_ISOLINES] = XtVaCreateManagedWidget("isol",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, bisolpm,
						     XmNselectPixmap, bisolpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								      NULL);
    XtAddCallback(main_display_item[BACKGROUND_GRID_ISOLINES],
		  XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) BACKGROUND_GRID_ISOLINES);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    main_display_item[BACKGROUND_BOUNDARY] = XtVaCreateManagedWidget("backbound",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						   XmNlabelPixmap, bboundpm,
						     XmNselectPixmap, bboundpm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
								     NULL);
    XtAddCallback(main_display_item[BACKGROUND_BOUNDARY], XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) BACKGROUND_BOUNDARY);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);
    XtVaCreateManagedWidget("BldP", xmLabelGadgetClass, rctmp, NULL);
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    main_display_item[BUILD_POINTS] = XtVaCreateManagedWidget("buildpoints",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						   XmNlabelPixmap, buildppm,
						     XmNselectPixmap, buildppm,
						     XmNfillOnSelect, False,
						      XmNshadowThickness, 3,
							      XmNspacing, 0,
						      XmNindicatorOn, False,
							XmNindicatorSize, 0,
							  XmNentryBorder, 0,
							    XmNmarginTop, 0,
							 XmNmarginBottom, 0,
							 XmNmarginHeight, 0,
							  XmNmarginWidth, 0,
							   XmNmarginLeft, 0,
							  XmNmarginRight, 0,
							      NULL);
    XtAddCallback(main_display_item[BUILD_POINTS], XmNvalueChangedCallback, (XtCallbackProc) togglegrid_proc, (XtPointer) BUILD_POINTS);
    display_flags[BUILD_POINTS] = 1;
    invert_pixmap(pm[BUILD_POINTS], 1);
    XtVaSetValues(main_display_item[BUILD_POINTS], XmNlabelPixmap, savepm[BUILD_POINTS], NULL);
    update_display_items();

    abort_button = XtVaCreateManagedWidget("Stop", xmPushButtonWidgetClass, rctmp,
					   XmNsensitive, False,
					   NULL);
    XtAddCallback(abort_button, XmNactivateCallback, (XtCallbackProc) do_abort, (XtPointer) NULL);
    abort_win = XtWindow(abort_button);
    bt = XtVaCreateManagedWidget("Exit", xmPushButtonWidgetClass, rctmp,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_EXIT);

/*
    locate_grid_item = XtVaCreateManagedWidget("Locate grid", xmLabelGadgetClass, bb,
				    XmNlabelString, lstring = XmStringCreateLtoR("Locate:                                       ", charset),
				    XmNy, 45,
				    NULL);
*/
    XtManageChild(main_frame);
    XmMainWindowSetAreas(main_frame, menu_bar, NULL, NULL, NULL, form);
    XtRealizeWidget(app_shell);

    icon = XCreateBitmapFromData(XtDisplay(app_shell),
				 DefaultRootWindow(XtDisplay(app_shell)),
				 gredit_icon_bits, gredit_icon_width,
				 gredit_icon_height);
    XtSetArg(al[0], XtNiconPixmap, icon);
    XtSetArg(al[1], XtNiconMask, icon);
    XtSetValues(app_shell, al, 2);

    update_fuzz_items();
    inwin = 1;
    xlibinitgc(canvas);
    set_window(canvas);
    xlibinitcmap();
    xlibinit_pm();
    init_cursors();
    /*
     * if an image was placed on the command line, read it in
     */
    if (readimage) {
        read_image(image_filename);
    }
    XtAppMainLoop(app_con);
}

#include "bitmaps.h"

void init_pm(void)
{
    Display *disp = XtDisplay(app_shell);
    Window cwin = RootWindowOfScreen(XtScreen(app_shell));
    GC gc = DefaultGC(disp, DefaultScreen(disp));
    Pixmap ptmp;
    XGCValues gc_val;
    GC my_gc;
    int i;
    gc_val.function = GXcopy;
/*
	gc_val.background = BlackPixel(disp, DefaultScreen(disp));
*/
    my_gc = XCreateGC(disp, cwin, GCFunction, &gc_val);
    zoompm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, zoom_bits, 16, 16);
    XCopyPlane(disp, ptmp, zoompm, gc, 0, 0, 16, 16, 0, 0, 1);
    autopm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, auto_bits, 16, 16);
    XCopyPlane(disp, ptmp, autopm, gc, 0, 0, 16, 16, 0, 0, 1);
    shrinkpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, shrink_bits, 16, 16);
    XCopyPlane(disp, ptmp, shrinkpm, gc, 0, 0, 16, 16, 0, 0, 1);
    expandpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, expand_bits, 16, 16);
    XCopyPlane(disp, ptmp, expandpm, gc, 0, 0, 16, 16, 0, 0, 1);
    rightpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, right_bits, 16, 16);
    XCopyPlane(disp, ptmp, rightpm, gc, 0, 0, 16, 16, 0, 0, 1);
    leftpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, left_bits, 16, 16);
    XCopyPlane(disp, ptmp, leftpm, gc, 0, 0, 16, 16, 0, 0, 1);
    uppm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, up_bits, 16, 16);
    XCopyPlane(disp, ptmp, uppm, gc, 0, 0, 16, 16, 0, 0, 1);
    downpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, down_bits, 16, 16);
    XCopyPlane(disp, ptmp, downpm, gc, 0, 0, 16, 16, 0, 0, 1);
    gridpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, grid_bits, 16, 16);
    XCopyPlane(disp, ptmp, gridpm, gc, 0, 0, 16, 16, 0, 0, 1);
    isolpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, isol_bits, 16, 16);
    XCopyPlane(disp, ptmp, isolpm, gc, 0, 0, 16, 16, 0, 0, 1);
    bisolpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, isol_bits, 16, 16);
    XCopyPlane(disp, ptmp, bisolpm, gc, 0, 0, 16, 16, 0, 0, 1);
    bgridpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, grid_bits, 16, 16);
    XCopyPlane(disp, ptmp, bgridpm, gc, 0, 0, 16, 16, 0, 0, 1);
    nodenpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, noden_bits, 16, 16);
    XCopyPlane(disp, ptmp, nodenpm, gc, 0, 0, 16, 16, 0, 0, 1);
    node2pm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, node2_bits, 16, 16);
    XCopyPlane(disp, ptmp, node2pm, gc, 0, 0, 16, 16, 0, 0, 1);
    elemnpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, elemn_bits, 16, 16);
    XCopyPlane(disp, ptmp, elemnpm, gc, 0, 0, 16, 16, 0, 0, 1);
    depthnpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, depthn_bits, 16, 16);
    XCopyPlane(disp, ptmp, depthnpm, gc, 0, 0, 16, 16, 0, 0, 1);
    gridfpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, gridf_bits, 16, 16);
    XCopyPlane(disp, ptmp, gridfpm, gc, 0, 0, 16, 16, 0, 0, 1);
    eboundpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, ebound_bits, 16, 16);
    XCopyPlane(disp, ptmp, eboundpm, gc, 0, 0, 16, 16, 0, 0, 1);
    bboundpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, ebound_bits, 16, 16);
    XCopyPlane(disp, ptmp, bboundpm, gc, 0, 0, 16, 16, 0, 0, 1);
    buildppm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, buildp_bits, 16, 16);
    XCopyPlane(disp, ptmp, buildppm, gc, 0, 0, 16, 16, 0, 0, 1);
    pm[EDIT_GRID] = gridpm;
    pm[EDIT_GRID_ISOLINES] = isolpm;
    pm[EDIT_BOUNDARY] = eboundpm;
    pm[EDIT_GRID_FILLED] = gridfpm;
    pm[EDIT_GRID_NODE_NUMBERS] = nodenpm;
    pm[EDIT_GRID_NODES] = node2pm;
    pm[EDIT_GRID_ELEMENT_NUMBERS] = elemnpm;
    pm[EDIT_GRID_DEPTHS] = depthnpm;
    pm[BACKGROUND_GRID] = gridpm;
    pm[BACKGROUND_GRID_ISOLINES] = bisolpm;
    pm[BACKGROUND_BOUNDARY] = bboundpm;
    pm[BUILD_POINTS] = buildppm;
    for (i = 0; i < 12; i++) {
	savepm[i] = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
	XCopyArea(disp, pm[i], savepm[i], my_gc, 0, 0, 16, 16, 0, 0);
	{
	    XImage *ximage;
	    int h, w, pix;
	    ximage = XGetImage(disp, savepm[i], 0, 0, 16, 16, AllPlanes, ZPixmap);
	    for (h = 0; h < ximage->height; h++) {
		for (w = 0; w < ximage->width; w++) {
		    pix = XGetPixel(ximage, w, h);
		    if (pix == 1L) {
			XPutPixel(ximage, w, h, 0L);
		    } else {
			XPutPixel(ximage, w, h, 1L);
		    }
		}
	    }
	    XPutImage(disp, savepm[i], my_gc, ximage, 0, 0, 0, 0, 16, 16);
	}
    }
}

void invert_pixmap(Pixmap p, int flag)
{
    Display *disp = XtDisplay(app_shell);
    Window cwin = RootWindowOfScreen(XtScreen(app_shell));
    XGCValues gc_val;
    static GC my_gc;
    if (!my_gc) {
	gc_val.function = GXcopy;
	my_gc = XCreateGC(disp, cwin, GCFunction, &gc_val);
    }
/*
    if (flag) {
	XSetForeground(disp, my_gc, WhitePixel(disp, DefaultScreen(disp)));
	XSetBackground(disp, my_gc, BlackPixel(disp, DefaultScreen(disp)));
    }
    else {
	XSetForeground(disp, my_gc, BlackPixel(disp, DefaultScreen(disp)));
	XSetBackground(disp, my_gc, WhitePixel(disp, DefaultScreen(disp)));
    }
    XCopyArea(disp, p, p, my_gc, 0, 0, 15, 15, 0, 0);
    if (flag) {
    XSetWindowBackground(disp, p, BlackPixel(disp, DefaultScreen(disp)));
    }
    else {
    XSetWindowBackground(disp, p, WhitePixel(disp, DefaultScreen(disp)));
    }
*/
}
