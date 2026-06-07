/* $Id: viswin.c,v 1.1.1.1 2003/07/21 16:18:42 pturner Exp $
 * 
 * popup for 1d visualization
 */
#include <stdio.h>
#include <math.h>

#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Frame.h>
#include <Xm/LabelG.h>
#include <Xm/PushBG.h>
#include <Xm/PushB.h>
#include <Xm/ToggleBG.h>
#include <Xm/Scale.h>
#include <Xm/SeparatoG.h>
#include <Xm/RowColumn.h>

#include "defines.h"
#include "globals.h"
#include "vis.h"

Widget *CreatePanelChoice();
Widget CreateTextItem2();
Widget *CreateColorChoice();

/* viswin.c */
void updatevis(int vno);
void set_vis_proc(Widget w, int vno);
void run_vis_proc(void);
void pause_vis_proc(void);
void stop_vis_proc(void);
void create_vis_frame(Widget w, Widget cd);
int vis_error(FILE *fp, int src, char *msg);
void create_fli(char *fname, int start, int stop, int res);
void create_record_frame(void);

extern vis_struct vis[];

#define RUN 0
#define VIS_STOP 1
#define PAUSE 2
#define MAXVIS 20

Widget vis_frame = (Widget) NULL;
Widget vis_panel;

/*
 * Panel item declarations
 */
Widget vis_command;
Widget vis_file_item;
Widget vis_time_item;
Widget *vis_select_item;
Widget vis_skip_item;
Widget *vis_graph_item;
Widget *vis_set_item;
/*
 *  Needed definition for starting work processes under later releases of X11.
 */
extern XtAppContext app_con;
extern Widget app_shell;
/*
 * Event and Notify proc declarations
 */
static void vis_Done_notify_proc(void);
void run_vis_proc(void);

XtWorkProcId wd = NULL;

static int busy, mode;
XtWorkProc read_vis();
static int first_pic = 1;
static int cur_vis = 0;
int vis_time = 0;

void updatevis(int vno)
{
    char s[25];
    XmToggleButtonGadgetSetState(vis_file_item, vis[cur_vis].src == DISK, False);
    XmToggleButtonGadgetSetState(vis_time_item, vis_time, False);
    xv_setstr(vis_command, vis[cur_vis].fname);
    sprintf(s, "%d", vis[cur_vis].skip - 1);
    xv_setstr(vis_skip_item, s);
    if (vis[cur_vis].active == ON) {
	SetChoice(vis_set_item, vis[cur_vis].set1 + 1);
	SetChoice(vis_graph_item, vis[cur_vis].gno + 1);
    } else {
	SetChoice(vis_set_item, 0);
	SetChoice(vis_graph_item, 0);
    }
}

void set_vis_proc(Widget w, int vno)
{
    cur_vis = vno;
    updatevis(vno);
}

/*
 * accept this file/command setting
 */
static void accept_vis_proc(void)
{
    int v = GetChoice(vis_select_item);
    int s = GetChoice(vis_set_item) - 1;
    int g = GetChoice(vis_graph_item) - 1;
    if (g < 0) {
	g = cg;
    }
    if (s < 0) {
	if ((s = nextset(g)) == -1) {
	    return;
	}
    }
    vis[cur_vis].set1 = s;
    SetChoice(vis_set_item, vis[cur_vis].set1 + 1);
    vis[cur_vis].active = ON;
    vis[cur_vis].gno = g;
    SetChoice(vis_graph_item, vis[cur_vis].gno + 1);
    vis[cur_vis].src = XmToggleButtonGadgetGetState(vis_file_item) ? DISK : PIPE;
    vis_time = XmToggleButtonGadgetGetState(vis_time_item);
    strcpy(vis[cur_vis].fname, (char *) xv_getstr(vis_command));
    vis[cur_vis].skip = atoi((char *) xv_getstr(vis_skip_item)) + 1;
    if (vis[cur_vis].skip <= 0) {
	vis[cur_vis].skip = 1;
    }
    readbin_vis(cur_vis, RUN);
    readbin_vis(cur_vis, VIS_STOP);
}

/*
 * run vis
 */
void run_vis_proc(void)
{
    if (!wd) {
        set_vismode(RUN);
	wd = XtAppAddWorkProc(app_con, read_vis, RUN);
	mode = RUN;
	first_pic = 1;
    } else {
	errwin("Simulation is already running");
    }
}

void pause_vis_proc(void)
{
    int i;
    set_vismode(PAUSE);
    for (i = 0; i < MAXVIS; i++) {
	if (vis[i].active == ON) {
	    readbin_vis(i, PAUSE);
	}
    }
    if (wd) {
	XtRemoveWorkProc(wd);
    }
    wd = NULL;
}

void stop_vis_proc(void)
{
    int i;
    set_vismode(VIS_STOP);
    for (i = 0; i < MAXVIS; i++) {
	if (vis[i].active == ON) {
	    readbin_vis(i, VIS_STOP);
	}
    }
    if (wd) {
	XtRemoveWorkProc(wd);
    }
    wd = NULL;
}

/*
 * Create the draw Frame and the draw Panel
 */
void create_vis_frame(Widget w, Widget cd)
{
    Widget wbut, rc;
    Widget wlabel;
    int i;

    if (vis_frame) {
	updatevis(cur_vis);
	XtManageChild(vis_panel);
	XtManageChild(vis_frame);
	return;
    }
    vis_frame = XmCreateDialogShell(app_shell, "Vis", NULL, 0);
    vis_panel = XmCreateRowColumn(vis_frame, "vis panel", NULL, 0);

    vis_select_item = CreatePanelChoice(vis_panel, "Use visualization set: ",
					21,
					"0", "1", "2", "3", "4",
					"5", "6", "7", "8", "9",
					"10", "11", "12", "13", "14",
					"15", "16", "17", "18", "19",
					NULL,
					NULL);
    for (i = 0; i < MAXVIS; i++) {
	XtAddCallback(vis_select_item[2 + i], XmNactivateCallback, (XtCallbackProc) set_vis_proc, (XtPointer) i);
    }
    vis_command = CreateTextItem2(vis_panel, 40, "File/command line:  ");

    rc = XmCreateRowColumn(vis_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    vis_skip_item = CreateTextItem2(rc, 5, "Time step skip:");
    vis_file_item = XtVaCreateManagedWidget("Data is from a file", xmToggleButtonGadgetClass, rc,
				 NULL);
    vis_time_item = XtVaCreateManagedWidget("Display time", xmToggleButtonGadgetClass, rc, NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(vis_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    vis_set_item = CreatePanelChoice(rc, "Use gr set: ",
				     12,
				     "Next", "0", "1", "2", "3", "4",
				     "5", "6", "7", "8", "9",
				     NULL,
				     NULL);

    vis_graph_item = CreatePanelChoice(rc, "in graph: ",
				       12,
				       "Current", "0", "1", "2", "3", "4",
				       "5", "6", "7", "8", "9",
				       NULL,
				       NULL);
    XtManageChild(rc);

    rc = XmCreateRowColumn(vis_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) accept_vis_proc, 0);
    XtManageChild(rc);

    XtVaCreateManagedWidget("sep1", xmSeparatorGadgetClass, vis_panel, NULL);
    rc = XmCreateRowColumn(vis_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Run", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) run_vis_proc, 0);

    wbut = XtVaCreateManagedWidget("Pause", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) pause_vis_proc, 0);

    wbut = XtVaCreateManagedWidget("Stop", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) stop_vis_proc, 0);

    wbut = XtVaCreateManagedWidget("*Record...", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_record_frame, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc,
			NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) vis_Done_notify_proc, 0);
    XtManageChild(rc);

    updatevis(cur_vis);
    XtManageChild(vis_panel);
    XtManageChild(vis_frame);
}				/* end create_vis_panel */

static void vis_Done_notify_proc(void)
{
     XtUnmanageChild(vis_frame);
}

vis_error(FILE *fp, int src, char *msg)
{
    if (src == PIPE) {
        if (fp != NULL) {
	    pclose(fp);
        }
    } else if (src == DISK) {
        if (fp != NULL) {
	    fclose(fp);
        }
    }
    if (msg != NULL) {
	errwin(msg);
    }
}

int save_fli;

void create_fli(char *fname, int start, int stop, int res)
{
    int i, nf = 0;
    for (i = start - 1; i < stop; i++) {
	nf++;
    }
    save_fli = 1;
    initialize_fli(fname, nf, res);
    initmake_fli();
    set_vismode(RUN);
    for (i = start - 1; i < stop; i++) {
	read_vis();
	add_image(i - start + 2); /* start with image counting from 1 */
    }
    save_fli = 0;
    close_fli();
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
    extern Window xwin;
    save_image_on_disk(disp, xwin, displaybuff, 0, 0, win_w, win_h, "tmp.xwd", NULL);
    system("xpaint tmp.xwd");
}

void destroy_dialog();

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
    flcres_item = CreatePanelChoice(rc, "Resolution: ",
                                                7,
                                                "320x200", 
						"640x400",
						"640x480",
						"800x600",
						"1024x768",
						"1280x1024",
                                                0, 0);

    flcborder_color_item = CreateColorChoice(rc, "Border color: ", 1);
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Record/close", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Record/pause", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Include image", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) genflc_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Edit image", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_edit_image, 0);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, frame);
    XtManageChild(rc);

    XtManageChild(panel);
    XtRaise(frame);
    SetChoice(flcres_item, 2);
    update_genflc();
}
