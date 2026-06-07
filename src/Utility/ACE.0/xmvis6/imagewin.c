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
 * Display an image along with the graph
 *
 */
#ifndef lint
static char RCSid[] = "$Id: imagewin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern Widget canvas;
extern XmStringCharSet charset;

extern int readimage;
extern char image_filename[];
extern double imagex;
extern double imagey;

static Widget image_frame;
static Widget image_name_item;
static Widget image_x_item;
static Widget image_y_item;
static Widget image_display_item;
XImage *img = NULL;
int imagew, imageh;
int drawimage_flag = 1;
static int open_image_dialog = 0;
void do_rimage_proc(Widget w, XtPointer client_data, XtPointer call_data);

void create_rimage_popup(Widget w, XtPointer client_data, XtPointer call_data);

XImage *read_gif(Display * disp, Window win_id, char *name_of_file, int *width, int *height, int *depth, int pflag);

XImage *read_image_from_disk(Display * disp, Window win_id, char *name_of_file, int *width, int *height, int *depth, int pflag);

void drawimage()
{
    int sx, sy;
    Window xwin;
    Display *disp;
    GC gc;
    if (img != NULL && drawimage_flag) {
	xwin = XtWindow(canvas);
	disp = XtDisplay(canvas);
	gc = DefaultGC(disp, DefaultScreen(disp));
	world2deviceabs(imagex, imagey, &sx, &sy);
/*
printf("imagex %lf, imagey %lf, sx %d, sy %d, wx %d, wy %d\n", imagex, imagey, sx, sy, imagew, imageh); 
*/
	XPutImage(disp, xwin, gc, img, 0, 0, sx, sy, imagew, imageh);
    }
}

static void update_image()
{
    char buf[256];
    XmToggleButtonSetState(image_display_item, drawimage_flag, False);
    sprintf(buf, "%lf", imagex);
    xv_setstr(image_x_item, buf);
    sprintf(buf, "%lf", imagey);
    xv_setstr(image_y_item, buf);
}

void do_accept_image_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    drawimage_flag = XmToggleButtonGetState(image_display_item);
    imagex = atof((char *) xv_getstr(image_x_item));
    imagey = atof((char *) xv_getstr(image_y_item));
    do_drawgrid();
}

void create_image_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget dialog;
    Widget wbut, rc;
    Widget but3[3];

    if (img == NULL) {
	open_image_dialog = 1;
	create_rimage_popup(w, client_data, call_data);
	return;
    }
    set_wait_cursor();
    if (image_frame == NULL) {
	char *label3[3];
	label3[0] = "Accept";
	label3[1] = "Read image...";
	label3[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	image_frame = XmCreateDialogShell(app_shell, "Image", NULL, 0);
	handle_close(image_frame);
	XtVaSetValues(image_frame, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(image_frame, "dialog_rc", NULL, 0);

/*    image_name_item = CreateTextItem2(dialog, 15, "Image file name: ");*/
	image_x_item = CreateTextItem2(dialog, 20, "Anchor top left to X (world coords) = : ");
	image_y_item = CreateTextItem2(dialog, 20, "Anchor top left to Y (world coords) = : ");
	image_display_item = XtVaCreateManagedWidget("Display image", xmToggleButtonWidgetClass, dialog, NULL);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 3, but3, label3);
	XtAddCallback(but3[0], XmNactivateCallback, (XtCallbackProc) do_accept_image_proc, (XtPointer) NULL);
	XtAddCallback(but3[1], XmNactivateCallback, (XtCallbackProc) create_rimage_popup, (XtPointer) NULL);
	XtAddCallback(but3[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) image_frame);

	XtManageChild(dialog);
    }
    XtRaise(image_frame);
    update_image();
    unset_wait_cursor();
}

static Widget rimage_dialog;

void close_rimage_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    XtUnmanageChild(rimage_dialog);
}

void read_image(char *fname)
{
    int width, height, depth;
    Window xwin;
    Display *disp;
    GC gc;
    if (img != NULL) {
	XDestroyImage(img);
    }
    xwin = XtWindow(canvas);
    disp = XtDisplay(canvas);
    gc = DefaultGC(disp, DefaultScreen(disp));
    set_wait_cursor();
    if ((img = read_gif(disp, xwin, fname, &width, &height, &depth, 0)) == NULL) {
	img = read_image_from_disk(disp, xwin, fname, &width, &height, &depth, 0);
    }

    unset_wait_cursor();
    if (img != NULL) {
	imagew = width;
	imageh = height;
    } else {
	errwin("Unable to load image");
    }
}

void do_rimage_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int width, height, depth;
    Display *disp;
    Window xwin;
    GC gc;
    char *s;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in do_rimage_proc()");
	return;
    }
    strcpy(image_filename, s);
    XtFree(s);
    read_image(image_filename);
    if (img != NULL) {
	do_drawgrid();
    }
    unset_wait_cursor();

    XtUnmanageChild(rimage_dialog);
    if (open_image_dialog) {
	open_image_dialog = 0;
	create_image_frame((Widget) NULL, (XtPointer) NULL, (XtPointer) NULL);
    }
}

void create_rimage_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (rimage_dialog == NULL) {
	rimage_dialog = XmCreateFileSelectionDialog(app_shell, "rimage_dialog", NULL, 0);
	XtVaSetValues(XtParent(rimage_dialog), XmNtitle, "Read image (.xwd format)", NULL);
	XtVaSetValues(rimage_dialog, XmNdirMask, XmStringCreate("*.xwd", charset), NULL);
	XtAddCallback(rimage_dialog, XmNcancelCallback, (XtCallbackProc) close_rimage_popup, (XtPointer) NULL);
	XtAddCallback(rimage_dialog, XmNokCallback, (XtCallbackProc) do_rimage_proc, (XtPointer) NULL);
    }
    XtRaise(rimage_dialog);
}
