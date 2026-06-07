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
 * Toggle streamlines and wrap
 *
 */

#ifndef lint
static char RCSid[] = "$Id: controlwin.c,v 1.3 2003/08/05 05:47:02 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;

static Widget frame;
static Widget wrap_button;
static Widget quit_button;
static Widget overlay_button;
static Widget dim_item;
static Widget save_button;
static Widget save_fname;
static Widget save_type;
static Widget save_start_count;
static Widget gl_tog;
static Widget shrink_item;
static Widget gwindscroll_item;
static Widget startstop_item;

/* used to set up XmStrings */
extern XmStringCharSet charset;

void page_right(void)
{
    page(1.0, 1);
}

void page_left(void)
{
    page(1.0, 0);
}

void page_up(void)
{
    page(1.0, 3);
}

void page_down(void)
{
    page(1.0, 2);
}

void do_shrink(void)
{
    page(1.0, 5);
}

void do_expand(void)
{
    page(1.0, 4);
}

void update_controls(void)
{
    char buf[128];
    extern int win_w, win_h;
    int value;
    if (frame) {
	XmToggleButtonSetState(overlay_button, overlay, False);
	XmToggleButtonSetState(wrap_button, get_wrap(), False);
	XmToggleButtonSetState(save_button, save_images, False);
	xv_setstr(save_fname, save_images_fname);
	sprintf(buf, "%d %d\n", win_w, win_h);
	value = shexper * 100.0;
	XtVaSetValues(shrink_item, XmNvalue, value, NULL);
	value = scrollper * 100.0;
	XtVaSetValues(gwindscroll_item, XmNvalue, value, NULL);
    }
}

void control_accept_proc(void)
{
    Arg a;
    extern int tdevice, overlay;
    int w, h;
    char buf[128];
    int value;
    overlay = XmToggleButtonGetState(overlay_button);
    set_wrap(XmToggleButtonGetState(wrap_button));
    tdevice = XmToggleButtonGetState(gl_tog) ? 7 : 0;
    save_images = XmToggleButtonGetState(save_button) ? 1 : 0;
    if (save_images) {
	strcpy(save_images_fname, (char *) xv_getstr(save_fname));
	strcpy(buf, (char *) xv_getstr(dim_item));
	sscanf(buf, "%d %d\n", &w, &h);
	if (w != 0 && h != 0) {
	    set_canvas_size(w, h);
	}
    }
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(shrink_item, &a, 1);
    shexper = value / 100.0;
    XtSetArg(a, XmNvalue, &value);
    XtGetValues(gwindscroll_item, &a, 1);
    scrollper = value / 100.0;
}

void set_save_images(int d)
{
    save_images = d;
    update_controls();
}

void create_panel_items(void)
{
    int n, rowx = 10, rowy = 20;
    Widget bt, fr, rc, panel;
    if (frame) {
	update_controls();
	XtRaise(frame);
	return;
    }
    frame = XmCreateDialogShell(app_shell, "Controls", NULL, 0);
    handle_close(frame);
    panel = XmCreateRowColumn(frame, "controls_rc", NULL, 0);
    wrap_button = XtVaCreateManagedWidget("Wrap", xmToggleButtonWidgetClass, panel, NULL);
    overlay_button = XtVaCreateManagedWidget("Overlay frames", xmToggleButtonWidgetClass, panel, NULL);
    gl_tog = XtVaCreateManagedWidget("*GL mode", xmToggleButtonWidgetClass, panel, NULL);
    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    gwindscroll_item = XtVaCreateManagedWidget("Scroll %", xmScaleWidgetClass, rc,
					       XmNwidth, 150,
					       XmNminimum, 0,
					       XmNmaximum, 100, XmNvalue, (int) (scrollper * 100), XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNtitleString, XmStringCreateLtoR("Scroll %", charset), XmNorientation, XmHORIZONTAL, NULL);

    shrink_item = XtVaCreateManagedWidget("shrink", xmScaleWidgetClass, rc,
					  XmNwidth, 150,
					  XmNminimum, 0,
					  XmNmaximum, 100, XmNvalue, (int) (100 * shexper), XmNshowValue, True, XmNprocessingDirection, XmMAX_ON_RIGHT, XmNtitleString, XmStringCreateLtoR("Expand/shrink %", charset), XmNorientation, XmHORIZONTAL, NULL);
    XtManageChild(rc);

    fr = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    save_button = XtVaCreateManagedWidget("Save images", xmToggleButtonWidgetClass, rc, NULL);
    save_fname = CreateTextItem2(rc, 15, "Save images basename: ");
    dim_item = CreateTextItem2(rc, 15, "*Canvas dimensions (w h): ");
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) control_accept_proc, 0);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, frame);
    XtManageChild(rc);

    XtManageChild(panel);
    XtRaise(frame);
    update_controls();
}

static Widget flcfname_item;
static Widget flcstart_item;
static Widget flcstop_item;
static Widget *flcres_item;
static Widget *flcborder_color_item;

static void update_genflc(void)
{
}

static void genflc_accept_proc(void)
{
    char buf[255];
    int start, stop, res;
    strcpy(buf, (char *) xv_getstr(flcstart_item));
    sscanf(buf, "%d\n", &start);
    strcpy(buf, (char *) xv_getstr(flcstop_item));
    sscanf(buf, "%d\n", &stop);
    strcpy(buf, (char *) xv_getstr(flcfname_item));
    res = GetChoice(flcres_item);
    if (!fexists(buf)) {
	create_fli(buf, start, stop, res);
    }
}

static void do_edit_image(void)
{
    extern Display *disp;
    extern int win_w, win_h;
    extern Pixmap displaybuff;
    extern Window cwin;
    save_image_on_disk(disp, cwin, displaybuff, 0, 0, win_w, win_h, "tmp.xwd", NULL);
    system("xpaint tmp.xwd");
}

void create_record_frame(void)
{
    int n;
    Widget bt, fr, rc, panel;
    static Widget frame;
    if (frame) {
	update_genflc();
	XtRaise(frame);
	return;
    }
    frame = XmCreateDialogShell(app_shell, "Controls", NULL, 0);
    handle_close(frame);
    panel = XmCreateRowColumn(frame, "controls_rc", NULL, 0);
    fr = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    flcfname_item = CreateTextItem2(rc, 15, "FLC filename: ");
    flcstart_item = CreateTextItem2(rc, 15, "Start at frame: ");
    flcstop_item = CreateTextItem2(rc, 15, "Stop at frame: ");
    flcres_item = CreatePanelChoice1(rc, "Resolution: ", 7, "320x200", "640x400", "640x480", "800x600", "1024x768", "1280x1024", 0, 0);

    flcborder_color_item = CreateColorChoice(rc, "Border color: ", 1);
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Record/close", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Record/pause", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Include image", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Edit image", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_edit_image, 0);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, frame);
    XtManageChild(rc);

    XtManageChild(panel);
    XtRaise(frame);
    SetChoice(flcres_item, 2);
    update_genflc();
}
