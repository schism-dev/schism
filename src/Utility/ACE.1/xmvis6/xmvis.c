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
 * Xmvis.c - initialize UI
 *
 */

#ifndef lint
static char RCSid[] = "$Id: xmvis.c,v 1.9 2005/04/07 17:30:16 pturner Exp $";
#endif

#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <X11/cursorfont.h>

#include <Xm/Xm.h>
#include <Xm/CascadeB.h>
#include <Xm/DrawingA.h>
#include <Xm/Frame.h>
#include <Xm/DialogS.h>
#include <Xm/BulletinB.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/ArrowB.h>
#include <Xm/ArrowBG.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/ToggleB.h>
#include <Xm/MainW.h>
#include <Xm/MessageB.h>
#include <Xm/Label.h>
#include <Xm/LabelG.h>
#include <Xm/RowColumn.h>
#include <Xm/SelectioB.h>
#include <Xm/SeparatoG.h>
#include <Xm/Scale.h>

Cursor wait_cursor;

#include "bitmaps.h"

Pixmap zoompm, shrinkpm, expandpm, autopm, elevpm;
Pixmap concpm, flowpm, pathlpm, pointspm, elevmpm, gridpm;
Pixmap uppm, leftpm, downpm, rightpm, isolpm, bgridpm;
Pixmap nodenpm, elemnpm, depthnpm;
Pixmap eboundpm, bboundpm, gridfpm, buildppm;

Widget main_display_item[30];

static Pixel fg, bg;

void init_pm(void);

#include "defines.h"
#include "globals.h"
#include "vis_icon.h"

#define MENU_HELP		200
#define MENU_EXIT		201

/* defines a temporary file for file transfers */

Widget app_shell;		/* ApplicationShell 	 */
Widget main_frame;		/* MainWindow */
Widget menu_pane;		/* RowColumn	 		 */
Widget menu_bar;		/* RowColumn */
Widget canvas;			/* drawing area */
Widget loclab;			/* label for locator */
XmString string;		/* string for current location */
static Widget steps_item;

XmString mstring;		/* string for mode indicator */
Widget write_mode_item;

XmString statstring;		/* string for mode indicator */
Widget statlab;

XmString lstring;		/* string for grid locator */
extern Widget locate_grid_item;

XmString sdstring;		/* string for stack depth */
Widget stack_depth_item;

XmString cystring;		/* string for stack cycle */
Widget curw_item;

void page_left(void);
void page_right(void);
void page_up(void);
void page_down(void);
void do_shrink(void);
void do_expand(void);
void my_proc(Widget w, XtPointer data, XEvent * event);
void refresh(Widget w, XtPointer clid, XtPointer calld);
void resize(Widget w, XtPointer clid, XtPointer calld);
void setirun(void);
void create_area_frame(void);
void create_teanl_frame(void);
void create_sadcirc_frame(void);
void create_steanl_frame(void);
void create_shist_frame(void);
void create_adcirc_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_elcirc_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_sampflow_frame(Widget w, XtPointer clid, XtPointer calld);
void create_setadc3d_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_settrans_frame(Widget w, XtPointer clid, XtPointer calld);
void create_surf_frame(Widget w, XtPointer clientd, XtPointer calld);
void create_track_frame(Widget w, XtPointer clientd, XtPointer calld);
void NewTransFrame(Widget w, XtPointer clid, XtPointer calld);
void NewTransDisplayFrame(Widget w, XtPointer clid, XtPointer calld);
void NewSigmaTransFrame(Widget w, XtPointer clid, XtPointer calld);
void NewSigmaTransDisplayFrame(Widget w, XtPointer clid, XtPointer calld);
void SigmaZTransFrame(Widget w, XtPointer clid, XtPointer calld);
void SigmaZTransDisplayFrame(Widget w, XtPointer clid, XtPointer calld);
void create_hist_frame(void);
void create_svhist_frame(void);
void create_ela_frame(void);
void create_sela_frame(void);
void create_pathl_frame(void);
void create_pathls_frame(void);
void create_tide_frame(void);
void create_stide_frame(void);
void create_drogues_frame(void);
void create_slice_frame(void);
void create_flux_frame(void);
void create_zoom_frame(void);
void create_config_frame();
void create_vel_frame(void);
void create_map_frame(void);
void create_fluxscale_frame(void);
void create_bisolines_popup();
void create_place_popup(void);
void create_panel_items(void);
void create_display_popup(void);
void create_cedit_popup(void);
void create_timeline_frame(void);
void create_tidalclock_frame(void);
void create_timeinfo_frame(void);
void create_world_frame(void);
void create_view_frame(void);
void create_frame_frame(void);
void create_printer_setup(void);
void create_ticks_frame(void);
void create_graph_frame(void);
void create_time_frame(void);
void create_bound_frame(void);
void create_boundio_frame(void);
void create_grid_frame(void);
void create_gridio_frame(void);
void create_gridt_frame(void);
void create_gridtio_frame(void);
void create_analyze_frame(void);
void create_graph2d_popup(int sno);
void create_new_frame(void);
void create_gredit_frame(void);
void create_record_frame(void);
void create_ship_frame(Widget w, XtPointer clientd, XtPointer calld);
void define_status_popup(void);
void open_command(void);
void setrewind(void);
void setreverse(void);
void setfaster(void);
void setslower(void);
void setzoom(void);
void setistep(Widget w, int cd);
void setistop(void);
void setfastforward(void);
void setfastreverse(void);
void setforward(void);
void setreset_world(void);
void setredraw_world(void);
void do_run(void);
void open_command(void);
void do_setreverse(void);
void do_setfastforward(void);
void do_setfastreverse(void);
void do_setforward(void);
void push_world(void);
void pop_world(void);
void push_and_zoom(void);
void cycle_world_stack(void);
void autoticks_proc(void);
void define_objects_popup(void);
void create_rparams_popup(void);
void create_wparam_frame(void);
void create_sslice_frame(void);
void create_sflux_frame(void);
void create_region_frame(void);
void do_select_single_region(void);
void do_clear_region2(void);
void do_select_area(void);
void do_select_peri(void);
/*
 * locator callbacks
 */
void open_locate1_popup(void);
void open_locate2_popup(void);
void open_locate3_popup(void);
void open_locate4_popup(void);
void open_locate6_popup(void);
void open_locate7_popup(void);
void do_select_point(void);
void do_clear_point(void);
void create_locator_frame(void);

extern int readimage;
extern char image_filename[];
void create_image_frame(Widget w, XtPointer client_data, XtPointer call_data);

/*
 * xlib objects for drawing
 */
extern Display *disp;

/* used to set up XmStrings */
XmStringCharSet charset = (XmStringCharSet) XmSTRING_DEFAULT_CHARSET;

XtAppContext app_con;
XVisualInfo vinfo;
Visual *visual;

void initialize_screen(int *argc, char **argv)
{
    app_shell = XtAppInitialize(&app_con, "XMvis", NULL, 0, argc, argv, NULL, NULL, 0);
    disp = XtDisplay(app_shell);
    if (!disp) {
	XtWarning("Can't open display, exiting...");
	exit(1);
    }
/*
    if (!XMatchVisualInfo(disp, DefaultScreen(disp), 8, PseudoColor, &vinfo)) {
        fprintf(stderr, "Need an 8 bit pseudocolor visual to work, server doesn't support it.\n");
        exit(1);
    }
    visual = vinfo.visual;
*/
    initialize_cms_data();
}

/*
 * cancel button on main window
 */
Widget abort_button;
Window abort_win;
int cancel_flag = 0;

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

void set_steps_label(int step)
{
    Arg al;
    XmString ls;
    char buf[256];
    double get_current_time(void);
    if (timeclock.nsteps == 0) {
	sprintf(buf, "Step 0 of 0: Current time 0.0s");
    } else {
	sprintf(buf, "Step %d of %d: Current time %.1lfs",
	      get_current_step() + 1, timeclock.nsteps, get_current_time());
    }
    ls = XmStringCreateLtoR(buf, charset);
    XtSetArg(al, XmNlabelString, ls);
    XtSetValues(steps_item, &al, 1);
    XmStringFree(ls);
    XmUpdateDisplay(steps_item);
}

Pixmap pm[11], savepm[11];
void set_display_pixmaps(void);

static Widget current_grid_item;

void set_current_grid(void)
{
    g[cg].curgrid = (g[cg].curgrid + 1) % MAXGRIDS;
    set_display_pixmaps();
}

void togglegrid_proc(Widget w, int cd, XmAnyCallbackStruct * cbs)
{
    void invert_pixmap(Pixmap p, int flag);
    int c = g[cg].curgrid;
    g[cg].grid[c].display_flags[cd] = XmToggleButtonGetState(main_display_item[cd]);
    if (g[cg].grid[c].display_flags[cd]) {
	XtVaSetValues(main_display_item[cd], XmNlabelPixmap, savepm[cd], NULL);
    } else {
	XtVaSetValues(main_display_item[cd], XmNlabelPixmap, pm[cd], NULL);
    }
    set_display_pixmaps();
}

void set_display_pixmaps(void)
{
    int i;
    XmString s;
    int c = g[cg].curgrid;
    sprintf(buf, "Grid %1d", g[cg].curgrid + 1);
    s = XmStringCreateLtoR(buf, charset);
    XtVaSetValues(current_grid_item, XmNlabelString, s,
		  NULL);
    /*XtFree(s);*/
    if (main_display_item[0]) {
	for (i = 0; i < ndisplay; i++) {
	    if (g[cg].grid[c].display_flags[i]) {
		XtVaSetValues(main_display_item[i], XmNlabelPixmap, savepm[i], NULL);
	    } else {
		XtVaSetValues(main_display_item[i], XmNlabelPixmap, pm[i], NULL);
	    }
	}
    }
    update_display_items();
}

void MenuCB(Widget w, XtPointer client_data, XmAnyCallbackStruct * cbs)
 /* widget id		 */

 /* data from application   */


{
    register int ac;		/* arg count		    */
    Arg al[10];			/* arg list		    */
    char *command;		/* command used in printing */
    Widget helpw;
    Cursor c;
    char pfile[1024];

    switch ((int) client_data) {
    case MENU_EXIT:
	if (yesno("Exit ACE/vis? Are you sure?", "Press Yes or No", "Yes", "No")) {
	    exit(0);
	}
	break;
    case MENU_HELP:
/*
	system("xmosaic help/xmvis.test &");
    c = XCreateFontCursor(XtDisplay(w), XC_hand2);
    if (helpw = XmTrackingLocate(w, c, False)) {
	cbs->reason = XmCR_HELP;
	XtCallCallbacks(helpw, XmNhelpCallback, &cbs);
    }
    XFreeCursor(XtDisplay(w), c);
	
*/
/*
	create_help_frame(w, client_data);
	errwin("Sorry, no help available in this version of ACE/vis");
*/
	break;
    default:
	fprintf(stderr, "Warning: in menu callback\n");
	break;
    }
}

static Widget CreateMenuBar(Widget parent)
{
    Widget menu_bar;
    Widget cascade;
    Widget button;

    Arg al[10];
    register int ac;

    XtSetArg(al[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);

    ac = 0;
    menu_bar = XmCreateMenuBar(parent, "menu_bar", al, ac);
    menu_pane = XmCreatePulldownMenu(menu_bar, "Files menu", al, 1);

    button = XmCreatePushButton(menu_pane, "New ACE/vis...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_new_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Read Parameters...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_rparams_popup, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Save Parameters...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_wparam_frame, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "TEA-NL...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_teanl_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "ADCIRC...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_adcirc_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Slabs...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_elcirc_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Samples...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_setadc3d_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Surface/Bottom...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_surf_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Old SELFE Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_settrans_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Z Grid Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) NewTransFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Sigma Grid Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) NewSigmaTransFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Sigma+Z Grid Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) SigmaZTransFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Drogues...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_pathl_frame, 0);
    XtManageChild(button);
/*  Not implemented or broken options.
    button = XmCreatePushButton(menu_pane, "Track...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_track_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Tide stations...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_tide_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Time dependent grids...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gridtio_frame, 0);
    XtManageChild(button);
*/

    button = XmCreatePushButton(menu_pane, "Time histories...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_hist_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Grids...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gridio_frame, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("Read image...", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_image_frame, (XtPointer) 0);

    button = XmCreatePushButton(menu_pane, "Boundaries...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_bound_frame, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Status...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_status_popup, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Print...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_printer_setup, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    ac = 0;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Exit", charset));
    ac++;
    XtSetArg(al[ac], XmNacceleratorText, XmStringCreateLtoR("F3", charset));
    ac++;
    XtSetArg(al[ac], XmNaccelerator, "<Key>F3:");
    ac++;
    button = XmCreatePushButton(menu_pane, "Exit", al, ac);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_EXIT);
    XtManageChild(button);

    ac = 0;
    XtSetArg(al[ac], XmNsubMenuId, menu_pane);
    ac++;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Files", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'F');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Files", al, ac);
    XtManageChild(cascade);

/*
 * Create Models
 */
    ac = 0;
    XtSetArg(al[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);
    menu_pane = XmCreatePulldownMenu(menu_bar, "Models menu", al, 1);

    button = XmCreatePushButton(menu_pane, "TEANL setup...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_steanl_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE setup...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_sadcirc_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Sample settings...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_setadc3d_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Sample flow...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_sampflow_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "SELFE Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) NewTransDisplayFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Sigma Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) NewSigmaTransDisplayFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Sigma+Z Transects...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) SigmaZTransDisplayFrame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Drogues setup...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_pathls_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Track setup...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_track_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Tide stations...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_stide_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Time histories of flows...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_svhist_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Time histories...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_shist_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Time dependent grids...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gridt_frame, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Slice...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_sslice_frame, 0);
    XtManageChild(button);
/*
    button = XmCreatePushButton(menu_pane, "*Particle tracking...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_drogues_frame, 0);
    XtManageChild(button);
*/
    button = XmCreatePushButton(menu_pane, "Region", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_single_region, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Clear region", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_clear_region2, 0);
    XtManageChild(button);
    button = XtVaCreateManagedWidget("Area/perimeter...", xmPushButtonGadgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_area_frame, 0);

    ac = 0;
    XtSetArg(al[ac], XmNsubMenuId, menu_pane);
    ac++;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Models", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'M');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Models", al, ac);
    XtManageChild(cascade);

/*
 * Create view
 */
    XtSetArg(al[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);
    menu_pane = XmCreatePulldownMenu(menu_bar, "View menu", al, 1);

    button = XmCreatePushButton(menu_pane, "Toggle display...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_display_popup, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Time...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_time_frame, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Time line...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_timeline_frame, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Tidal clock...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_tidalclock_frame, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Time info string...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_timeinfo_frame, 0);
    XtManageChild(button);
    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Controls...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_panel_items, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Map Scale...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_map_frame, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Velocity Scale...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_vel_frame, 0);
    XtManageChild(button);
    button = XmCreatePushButton(menu_pane, "Flux Scale...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_fluxscale_frame, 0);
    XtManageChild(button);
    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Grids...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_grid_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Boundaries...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_bound_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Zoom box...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_zoom_frame, 0);
    XtManageChild(button);

#ifndef CORPS
    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);
    button = XmCreatePushButton(menu_pane, "Graphs...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_graph_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "World...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_world_frame, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Viewport...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_view_frame, 0);
    XtManageChild(button);

    button = XtVaCreateManagedWidget("Tick labels/tick marks...", xmPushButtonGadgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_ticks_frame, 0);

    button = XmCreatePushButton(menu_pane, "Frame...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_frame_frame, 0);
    XtManageChild(button);
#endif

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XmCreatePushButton(menu_pane, "Text/lines/boxes...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_objects_popup, 0);
    XtManageChild(button);

    button = XmCreatePushButton(menu_pane, "Edit colors...", NULL, 0);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_cedit_popup, 0);
    XtManageChild(button);

    ac = 0;
    XtSetArg(al[ac], XmNsubMenuId, menu_pane);
    ac++;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("View", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'V');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "View", al, ac);
    XtManageChild(cascade);

    XtSetArg(al[0], XmNtearOffModel, XmTEAR_OFF_ENABLED);
    menu_pane = XmCreatePulldownMenu(menu_bar, "Locate menu", al, 1);
    button = XtVaCreateManagedWidget("Find nearest node", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate1_popup, NULL);

    button = XtVaCreateManagedWidget("Find element", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate2_popup, NULL);

    button = XtVaCreateManagedWidget("Find nearest element", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate3_popup, NULL);

    button = XtVaCreateManagedWidget("Find depth in grid", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate4_popup, NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Goto node/element...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate6_popup, NULL);

    button = XtVaCreateManagedWidget("Goto X, Y...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_locate7_popup, NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Set fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_select_point, 0);

    button = XtVaCreateManagedWidget("Clear fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_clear_point, 0);

    button = XtVaCreateManagedWidget("Locator props...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_locator_frame, NULL);

    ac = 0;
    XtSetArg(al[ac], XmNsubMenuId, menu_pane);
    ac++;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Locate", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'L');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Locate", al, ac);
    XtManageChild(cascade);

/*
    ac = 0;
    menu_pane = XmCreatePulldownMenu(menu_bar, "edit menu_pane", al, ac);
    button = XtVaCreateManagedWidget("Grids (ACE/gredit)...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_gredit_frame, NULL);
*/

    ac = 0;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Commands", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'C');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Commands", al, ac);
    XtAddCallback(cascade, XmNactivateCallback, (XtCallbackProc) open_command, 0);
    XtManageChild(cascade);

    ac = 0;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Info", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'I');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Info", al, ac);
/* TODO
    XtAddCallback(cascade, XmNactivateCallback, (XtCallbackProc) do_credits, 0);
*/
    XtManageChild(cascade);

    /*
     * Create "Help" button.
     */
    menu_pane = XmCreatePulldownMenu(menu_bar, "Help menu pane", NULL, 0);

    button = XtVaCreateManagedWidget("On context", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_HELP);

    ac = 0;
    XtSetArg(al[ac], XmNsubMenuId, menu_pane);
    ac++;
    XtSetArg(al[ac], XmNlabelString, XmStringCreateLtoR("Help", charset));
    ac++;
    XtSetArg(al[ac], XmNmnemonic, 'H');
    ac++;
    cascade = XmCreateCascadeButton(menu_bar, "Help", al, ac);
    XtManageChild(cascade);

    ac = 0;
    XtSetArg(al[ac], XmNmenuHelpWidget, cascade);
    ac++;
    XtSetValues(menu_bar, al, ac);

    return (menu_bar);
}

void set_canvas_size(int w, int h)
{

    Dimension dax = w, day = h;
    XtVaSetValues(main_frame,
		  XmNwidth, dax,
		  XmNheight, day,
		  NULL);
    resize(canvas, 0, 0);
}

int resize_policy = 1;

static void CreateDrawingArea(Widget parent)
{
    Dimension dax, day;
    XFontStruct *newfont;

    dax = winsetwidth;
    day = winsetheight;

    canvas = XtVaCreateManagedWidget("canvas", xmDrawingAreaWidgetClass, parent,
				     XmNwidth, dax,
				     XmNheight, day,
				   XmNnavigationType, XmEXCLUSIVE_TAB_GROUP,
	      XmNresizePolicy, resize_policy ? XmRESIZE_ANY : XmRESIZE_NONE,
				     XmNresizeWidth, False,
				     XmNresizeHeight, False,
			    XmNbackground, WhitePixel(XtDisplay(main_frame),
				      DefaultScreen(XtDisplay(main_frame))),
				     NULL);

    XtAddCallback(canvas, XmNexposeCallback, (XtCallbackProc) refresh, NULL);
    XtAddCallback(canvas, XmNresizeCallback, (XtCallbackProc) resize, NULL);
    XtAddEventHandler(canvas,
		      EnterWindowMask |
		      LeaveWindowMask |
		      ButtonPressMask |
		      PointerMotionMask |
		      Button1MotionMask |
		      KeyPressMask,
		      FALSE, (XtEventHandler) my_proc, NULL);
    return;
}

void do_main_loop(int argc, char **argv)
{
    Widget frtop, bt, bt2, fr, fr2, fr3, rctop, rc, rc2, rc3, bb, form, form2, rctmp;
    Arg al[10];
    int ac;
    XSetWindowAttributes sw;
    Pixmap icon;
    XEvent event;
    Dimension dax, day;
    int i, wx1, wy1;

    wx1 = DisplayWidth(XtDisplay(app_shell), DefaultScreen(XtDisplay(app_shell)));
    wy1 = DisplayHeight(XtDisplay(app_shell), DefaultScreen(XtDisplay(app_shell)));

    if (winsetwidth == 0) {	/* initialize the size if not given on the
				 * command line */
	winsetwidth = wx1 * 0.65;
	winsetheight = wy1 * 0.65;
    } else {
	resize_policy = 0;
    }
    dax = winsetwidth;
    day = winsetheight;

    ac = 0;
    XtSetArg(al[ac], XmNwidth, dax);
    ac++;
    XtSetArg(al[ac], XmNheight, day);
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
    frtop = XtVaCreateManagedWidget("frtop", xmFrameWidgetClass, form,
				    NULL);
    form2 = XmCreateForm(frtop, "form", NULL, 0);
    rctop = XtVaCreateManagedWidget("rctop", xmRowColumnWidgetClass, form2,
				    XmNorientation, XmHORIZONTAL,
				    XmNpacking, XmPACK_TIGHT,
				    XmNspacing, 0,
				    XmNentryBorder, 0,
				    XmNmarginWidth, 0,
				    XmNmarginHeight, 0,
				    NULL);

    CreateDrawingArea(form);
    XtManageChild(canvas);

    fr2 = XtVaCreateManagedWidget("fr2", xmFrameWidgetClass, form, NULL);
    XtManageChild(fr2);
    ac = 0;
/*
    XtSetArg(al[ac], XmNrecomputeSize, True);
    ac++;
*/
    rc2 = XmCreateForm(fr2, "form", al, ac);
    write_mode_item = XtVaCreateManagedWidget("modelabel", xmLabelGadgetClass, rc2,
					      XmNlabelString,
			   mstring = XmStringCreateLtoR("Idle...", charset),
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
		  XmNtopAttachment, XmATTACH_WIDGET,
		  XmNtopWidget, frtop,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(frtop,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(canvas,
		  XmNtopAttachment, XmATTACH_WIDGET,
		  XmNtopWidget, frtop,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, fr2,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_WIDGET,
		  XmNleftWidget, fr,
		  NULL);
    XtVaSetValues(fr2,
		  XmNbottomAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);

    init_pm();
    rctmp = rc;
    bt = XtVaCreateManagedWidget("Draw", xmPushButtonWidgetClass, rctmp,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setredraw_world, 0);

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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setzoom, 0);

    bt = XtVaCreateManagedWidget("AS", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, autopm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setreset_world, 0);

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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_expand, NULL);

    bt = XtVaCreateManagedWidget("z", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, shrinkpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_shrink, NULL);

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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) page_left, NULL);
    bt = XtVaCreateManagedWidget("Right", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, rightpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) page_right, NULL);

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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) page_down, NULL);
    bt = XtVaCreateManagedWidget("Up", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, uppm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) page_up, NULL);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);

	    bt = XtVaCreateManagedWidget("AT", xmPushButtonWidgetClass, rctmp,
					 NULL);
	    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoticks_proc, NULL);

	    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
					  XmNorientation, XmHORIZONTAL,
					  XmNpacking, XmPACK_TIGHT,
					  XmNspacing, 0,
					  XmNentryBorder, 0,
					  XmNmarginWidth, 0,
					  XmNmarginHeight, 0,
					  NULL);
	    bt = XtVaCreateManagedWidget("PZ", xmPushButtonWidgetClass, rc3,
					 NULL);
	    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) push_and_zoom, NULL);

	    bt = XtVaCreateManagedWidget("Pu", xmPushButtonWidgetClass, rc3,
					 NULL);
	    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) push_world, NULL);

	    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rctmp,
					  XmNorientation, XmHORIZONTAL,
					  XmNpacking, XmPACK_TIGHT,
					  XmNspacing, 0,
					  XmNentryBorder, 0,
					  XmNmarginWidth, 0,
					  XmNmarginHeight, 0,
					  NULL);
	    bt = XtVaCreateManagedWidget("Po", xmPushButtonWidgetClass, rc3,
					 NULL);
	    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) pop_world, NULL);

	    bt = XtVaCreateManagedWidget("Cy", xmPushButtonWidgetClass, rc3,
					 NULL);
	    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) cycle_world_stack, NULL);

	    sdstring = XmStringCreateLtoR("SD:0 ", charset);
	    stack_depth_item = XtVaCreateManagedWidget("stackdepth", xmLabelGadgetClass, rctmp,
						   XmNlabelString, sdstring,
						       NULL);

	    cystring = XmStringCreateLtoR("CW:0 ", charset);
	    curw_item = XtVaCreateManagedWidget("curworld", xmLabelGadgetClass, rctmp,
						XmNlabelString, cystring,
						NULL);
/* display features */

	ndisplay = 7;

	XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
				NULL);

	current_grid_item = XtVaCreateManagedWidget("Grid", xmLabelWidgetClass, rctmp, NULL);
	XtAddEventHandler(current_grid_item, ButtonPressMask, False,
			  (XtEventHandler) set_current_grid, 0);

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

	main_display_item[EDIT_GRID_ISOLINES] = XtVaCreateManagedWidget("isol",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						     XmNlabelPixmap, isolpm,
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

	main_display_item[EDIT_GRID_FILLED] = XtVaCreateManagedWidget("gridfilled",
					     xmToggleButtonWidgetClass, rc3,
						     XmNlabelType, XmPIXMAP,
						    XmNlabelPixmap, gridfpm,
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
    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);

    abort_button = XtVaCreateManagedWidget("Abort", xmPushButtonWidgetClass, rctmp,
					   XmNsensitive, False,
					   NULL);
    XtAddCallback(abort_button, XmNactivateCallback, (XtCallbackProc) do_abort, NULL);
    abort_win = XtWindow(abort_button);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, rctmp,
			    NULL);

    bt = XtVaCreateManagedWidget("Exit", xmPushButtonWidgetClass, rctmp,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_EXIT);

    bt = XtVaCreateManagedWidget("|<-", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setrewind, NULL);

    bt = XtVaCreateManagedWidget(" <- ", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_setreverse, NULL);

    bt = XtVaCreateManagedWidget(" 1<- ", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setistep, (XtPointer) - 1);

    bt = XtVaCreateManagedWidget("Stop", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setistop, NULL);

    bt = XtVaCreateManagedWidget("->1 ", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setistep, (XtPointer) 1);

    bt = XtVaCreateManagedWidget(" -> ", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_setforward, NULL);

    bt = XtVaCreateManagedWidget("->|", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_setfastforward, NULL);

    bt = XtVaCreateManagedWidget("Record", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_record_frame, NULL);

    steps_item = XtVaCreateManagedWidget("steps",
					 xmLabelWidgetClass, form2,
					 NULL);
    XtVaSetValues(rctop,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(steps_item,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);
    XtManageChild(form2);
    set_steps_label(0);

/*
    bt = XtVaCreateManagedWidget("Step", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setistep, NULL);

    bt = XtVaCreateManagedWidget("Rewind", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setrewind, NULL);
*/

/*
    steps_item = XtVaCreateManagedWidget("time", xmScaleWidgetClass, rctop,
                                    XmNwidth, 180,
                                    XmNminimum, 0,
                                    XmNmaximum, 1,
                                    XmNvalue, 0,
                                    XmNshowValue, False,
                                    XmNprocessingDirection, XmMAX_ON_RIGHT,
                                    XmNorientation, XmHORIZONTAL,
                                    NULL);

    maxsteps_item = XtVaCreateManagedWidget("0", xmLabelGadgetClass, rctop,
                          NULL);
*/
    bt = XtVaCreateManagedWidget("+", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setfaster, NULL);

    bt = XtVaCreateManagedWidget("-", xmPushButtonWidgetClass, rctop, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) setslower, NULL);

    XtManageChild(rctop);
    XtManageChild(frtop);

    XtManageChild(rctmp);
    XtManageChild(fr);

    XtManageChild(form);
    XtManageChild(main_frame);

    XmMainWindowSetAreas(main_frame, menu_bar, NULL, NULL, NULL, form);

    XtRealizeWidget(app_shell);

    icon = XCreateBitmapFromData(XtDisplay(app_shell),
				 DefaultRootWindow(XtDisplay(app_shell)),
				 vis_icon_bits, vis_icon_width,
				 vis_icon_height);
    XtSetArg(al[0], XtNiconPixmap, icon);
    XtSetArg(al[1], XtNiconMask, icon);
    XtSetValues(app_shell, al, 2);

    use_colors = DisplayPlanes(disp, DefaultScreen(disp));
    if (!g[cg].parmsread) {
	set_defaults(cg);
    } else {
	set_display_items();
    }
    set_display_pixmaps();

    inwin = 1;
    xlibinitgc(canvas);
    set_window(canvas);
    xlibinitcmap();
    xlibinit_pm();
    init_cursors();
    inwin = 0;
    XtAppMainLoop(app_con);
}

void do_set_window(void)
{
    set_window(canvas);
}

void init_pm(void)
{
    Display *disp = XtDisplay(app_shell);
    Window cwin = RootWindowOfScreen(XtScreen(app_shell));
    GC gc = DefaultGC(disp, DefaultScreen(disp));
    Pixmap ptmp;
    XGCValues gc_val;
    GC my_gc;
    int i;
    XtVaGetValues(menu_bar,
		  XmNforeground, &fg,
		  XmNbackground, &bg,
		  NULL);
    gc_val.function = GXcopy;
    my_gc = XCreateGC(disp, cwin, GCFunction, &gc_val);

    zoompm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, zoom_bits, 16, 16);
    XCopyPlane(disp, ptmp, zoompm, gc, 0, 0, 16, 16, 0, 0, 1);

    autopm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, auto_bits, 16, 16);
    XCopyPlane(disp, ptmp, autopm, gc, 0, 0, 16, 16, 0, 0, 1);

    shrinkpm = XCreatePixmap(disp, cwin, 16, 16,
			     DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, shrink_bits, 16, 16);
    XCopyPlane(disp, ptmp, shrinkpm, gc, 0, 0, 16, 16, 0, 0, 1);

    expandpm = XCreatePixmap(disp, cwin, 16, 16,
			     DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, expand_bits, 16, 16);
    XCopyPlane(disp, ptmp, expandpm, gc, 0, 0, 16, 16, 0, 0, 1);

    rightpm = XCreatePixmap(disp, cwin, 16, 16,
			    DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, right_bits, 16, 16);
    XCopyPlane(disp, ptmp, rightpm, gc, 0, 0, 16, 16, 0, 0, 1);

    leftpm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, left_bits, 16, 16);
    XCopyPlane(disp, ptmp, leftpm, gc, 0, 0, 16, 16, 0, 0, 1);

    uppm = XCreatePixmap(disp, cwin, 16, 16,
			 DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, up_bits, 16, 16);
    XCopyPlane(disp, ptmp, uppm, gc, 0, 0, 16, 16, 0, 0, 1);

    downpm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, down_bits, 16, 16);
    XCopyPlane(disp, ptmp, downpm, gc, 0, 0, 16, 16, 0, 0, 1);

    gridpm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, grid_bits, 16, 16);
    XCopyPlane(disp, ptmp, gridpm, gc, 0, 0, 16, 16, 0, 0, 1);

    isolpm = XCreatePixmap(disp, cwin, 16, 16,
			   DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, isol_bits, 16, 16);
    XCopyPlane(disp, ptmp, isolpm, gc, 0, 0, 16, 16, 0, 0, 1);

    nodenpm = XCreatePixmap(disp, cwin, 16, 16,
			    DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, noden_bits, 16, 16);
    XCopyPlane(disp, ptmp, nodenpm, gc, 0, 0, 16, 16, 0, 0, 1);

    elemnpm = XCreatePixmap(disp, cwin, 16, 16,
			    DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, elemn_bits, 16, 16);
    XCopyPlane(disp, ptmp, elemnpm, gc, 0, 0, 16, 16, 0, 0, 1);

    depthnpm = XCreatePixmap(disp, cwin, 16, 16,
			     DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, depthn_bits, 16, 16);
    XCopyPlane(disp, ptmp, depthnpm, gc, 0, 0, 16, 16, 0, 0, 1);

    gridfpm = XCreatePixmap(disp, cwin, 16, 16,
			    DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, gridf_bits, 16, 16);
    XCopyPlane(disp, ptmp, gridfpm, gc, 0, 0, 16, 16, 0, 0, 1);

    eboundpm = XCreatePixmap(disp, cwin, 16, 16,
			     DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, ebound_bits, 16, 16);
    XCopyPlane(disp, ptmp, eboundpm, gc, 0, 0, 16, 16, 0, 0, 1);

    pm[EDIT_GRID] = gridpm;
    pm[EDIT_GRID_ISOLINES] = isolpm;
    pm[EDIT_BOUNDARY] = eboundpm;
    pm[EDIT_GRID_FILLED] = gridfpm;
    pm[EDIT_GRID_NODE_NUMBERS] = nodenpm;
    pm[EDIT_GRID_ELEMENT_NUMBERS] = elemnpm;
    pm[EDIT_GRID_DEPTHS] = depthnpm;

    for (i = 0; i < 7; i++) {
	savepm[i] = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
	XCopyArea(disp, pm[i], savepm[i], my_gc, 0, 0, 16, 16, 0, 0);
	XFlush(disp);
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

Pixmap make_button_image(char *bits, int wx, int wy, int fg, int bg)
{
    Pixmap ptmp;
    Pixmap retval = (Pixmap) NULL;
    return retval;
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
}
