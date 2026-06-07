/* $Id: xmgr.c,v 1.3 2004/06/16 18:34:11 pturner Exp $

 * main loop
 *
 * Has Motif and X specific variable declarations
 *
 */

#include <stdio.h>
#include <stdlib.h>
#ifndef WIN32
#include <unistd.h>
#endif
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <signal.h>
#include <errno.h>

#include <X11/X.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <X11/keysym.h>
#include <X11/StringDefs.h>

#include <Xm/Xm.h>
#include <Xm/ArrowB.h>
#include <Xm/CascadeB.h>
#include <Xm/DialogS.h>
#include <Xm/DrawingA.h>
#include <Xm/BulletinB.h>
#include <Xm/FileSB.h>
#include <Xm/Frame.h>
#include <Xm/Form.h>
#include <Xm/MainW.h>
#include <Xm/MessageB.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/RowColumn.h>
#include <Xm/SelectioB.h>
#include <Xm/ToggleB.h>
#include <Xm/Separator.h>
#include <Xm/ScrolledW.h>

#include "xmgr_icon.h"
#include "globals.h"
#include "motifinc.h"

void create_workingdir_popup(Widget w, XtPointer client_data, XtPointer call_data);
void createADPdb(Widget w, XtPointer client_data, XtPointer call_data);
void create_vis_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_velocityp_frame(Widget w, XtPointer client_data, XtPointer call_data);
void CreateElcircTimeHistoriesFrame(Widget w, XtPointer client_data, XtPointer call_data);
void CreateElcircProfilesFrame(Widget w, XtPointer client_data, XtPointer call_data);

/* for Monitos in monitors.c */
void settimer(void);
void stoptimer(void);

/* used to set up XmStrings */
#if defined (XmFONTLIST_DEFAULT_TAG)
XmStringCharSet charset = (XmStringCharSet) XmFONTLIST_DEFAULT_TAG;
#else
XmStringCharSet charset = (XmStringCharSet) XmSTRING_DEFAULT_CHARSET;
#endif

XtAppContext app_con;

/* used globally */
Widget app_shell;
Widget canvas;
Widget pagew[4];		/* toggle buttons for page layouts */

static Widget scrollw;		/* container for drawing area */

Widget loclab;			/* locator label */
Widget statlab;			/* status line at the bottom */
Widget stack_depth_item;	/* stack depth item on the main panel */
Widget curw_item;		/* current world stack item on the main panel */
XmString string;		/* string for current location */
XmString sdstring;		/* string for stack depth */
XmString cystring;		/* string for stack cycle */
XmString statstring;		/* string for pointer status */

extern Colormap mycmap;		/* colormap for canvas */
extern Display *disp;
extern GC gc;
extern GC gcxor;
extern GC gcclr;
extern Window xwin;
extern unsigned long colors[];

extern XGCValues gc_val;

/* used locally */
static Widget main_frame;
static Widget menu_pane;
static Widget menu_bar;
static Widget frleft, frtop, frbot;	/* dialogs along canvas edge */
static Widget form;		/* form for mainwindow */

static void MenuCB(Widget w, XtPointer client_data, XtPointer call_data);
static Widget CreateMenuBar(Widget parent);
static void init_pm(void);

/*
 * for buttons on front panel
 */
#include "bitmaps.h"

static Pixmap zoompm, shrinkpm, expandpm, autopm;
static Pixmap uppm, leftpm, downpm, rightpm;

/*
 * establish resource stuff
 */
typedef struct {
    Boolean invert;
    Boolean revflag;
    Boolean backingstore;
    Boolean allow_dc;
    Boolean autoscale_onread;
    Boolean verify_action;
    int maxplot;
    int maxgraph;
    int maxcolors;
    Boolean noask;
    Boolean logwindow;
} ApplicationData, *ApplicationDataPtr;

static XtResource resources[] =
{
    {"invertDraw", "InvertDraw", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, invert), XtRImmediate,
     (XtPointer) FALSE},
    {"reverseVideo", "ReverseVideo", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, revflag), XtRImmediate,
     (XtPointer) FALSE},
    {"backingstore", "Backingstore", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, backingstore), XtRImmediate,
     (XtPointer) FALSE},
    {"allowDoubleClick", "AllowDoubleClick", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, allow_dc), XtRImmediate,
     (XtPointer) TRUE},
    {"autoscaleOnRead", "AutoscaleOnRead", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, autoscale_onread), XtRImmediate,
     (XtPointer) FALSE},
    {"verifyAction", "VerifyAction", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, verify_action), XtRImmediate,
     (XtPointer) FALSE},
    {"maxSets", "MaxSets", XtRInt, sizeof(int),
     XtOffset(ApplicationDataPtr, maxplot), XtRImmediate,
     (XtPointer) MAXPLOT},
    {"maxGraphs", "MaxGraphs", XtRInt, sizeof(int),
     XtOffset(ApplicationDataPtr, maxgraph), XtRImmediate,
     (XtPointer) MAXGRAPH},
    {"maxColors", "MaxColors", XtRInt, sizeof(int),
     XtOffset(ApplicationDataPtr, maxcolors), XtRImmediate,
     (XtPointer) MAXCOLORS},
    {"noAsk", "NoAsk", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, noask), XtRImmediate,
     (XtPointer) FALSE},
    {"logWindow", "LogWindow", XtRBoolean, sizeof(Boolean),
     XtOffset(ApplicationDataPtr, logwindow), XtRImmediate,
     (XtPointer) FALSE},
};

/*
 * put the current working directory in the title bar
 */
void set_title(char *ts)
{
    if (ts == NULL) {
	if (getcwd(buf, 1023) != NULL) {
	    strcpy(workingdir, buf);
	    XtVaSetValues(app_shell, XtNtitle, workingdir, NULL);
	}
    } else {
	XtVaSetValues(app_shell, XtNtitle, ts, NULL);
    }
}

/*
 * initialize the X-Toolkit
 */
void initialize_screen(int *argc, char **argv)
{
    ApplicationData rd;
    int i;

#ifdef WIN32
    HCLXmInit();
#endif

    app_shell = XtAppInitialize(&app_con, "XMgr", NULL, 0, argc, argv, NULL, NULL, 0);
    savewidget(app_shell);
    disp = XtDisplay(app_shell);
    if (!disp) {
	sprintf(buf, "%s: can't open display, exiting...", argv[0]);
	XtWarning(buf);
	exit(0);
    }
    use_colors = DisplayPlanes(disp, DefaultScreen(disp));
    if (use_colors < 8) {
	use_colors = 1;
    }
    XtGetApplicationResources(app_shell, &rd, resources,
			      XtNumber(resources), NULL, 0);
    invert = rd.invert;
    revflag = rd.revflag;
    backingstore = rd.backingstore;
    allow_dc = rd.allow_dc;
    autoscale_onread = rd.autoscale_onread;
    verify_action = rd.verify_action;
    maxplot = rd.maxplot;
    maxgraph = rd.maxgraph;
    maxcolors = rd.maxcolors;
    logwindow = rd.logwindow;
}

/*
 * stuff results, etc. into a text window
 */
void log_results(char *buf)
{
    char tmpbuf[512];
    if (logwindow) {
	strcpy(tmpbuf, buf);
	if (tmpbuf[strlen(tmpbuf) - 1] != '\n') {
	    strcat(tmpbuf, "\n");
	}
	stufftext(tmpbuf, 1);
    }
}

/*
 * exit ACE/gr
 */
void bailout(void)
{
    if (resfp) {
	fclose(resfp);
    }
    exit(0);
}

/*
 * main menubar
 */
#define MENU_HELP	200
#define MENU_EXIT	201
#define MENU_CLEAR	202
#define MENU_NEW	203
#define MENU_OPEN	204
#define MENU_SAVE	205
#define MENU_SAVEAS	206
#define MENU_PRINT	207

/*
   void HelpCB(Widget w, XtPointer client_data, XtPointer call_data);
 */

static void MenuCB(Widget w, XtPointer client_data, XtPointer call_data)
{
    extern int use_help;
    switch ((int) client_data) {
    case MENU_HELP:
	switch (use_help) {
	case 1:		/* dumb text window */
	    create_help_frame(NULL, NULL, NULL);
	    break;
	case 2:		/* xmosaic or other HTML viewer */
	    sprintf(buf, "cd %s/doc ; %s %s &", acegrdir, help_viewer, help_file);
	    set_wait_cursor();
	    system(buf);
	    unset_wait_cursor();
	    break;
	}
	break;
    case MENU_EXIT:
	if (yesno("Exit xmgr? Are you sure?", NULL, NULL, NULL)) {
	    bailout();
	}
	break;
    case MENU_CLEAR:
	wipeout(1);
	break;
    case MENU_NEW:
	break;
    case MENU_OPEN:
	break;
    case MENU_SAVE:
	break;
    case MENU_SAVEAS:
	break;
    case MENU_PRINT:
	set_wait_cursor();
	do_hardcopy();
	unset_wait_cursor();
	break;
    default:
	break;
    }
}

/*
 * service the autoscale buttons on the main panel
 */
void autoscale_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (activeset(cg)) {
	switch ((int) client_data) {
	case 0:
	    autoscale_graph(cg, -3);
	    break;
	case 1:
	    autoscale_graph(cg, -2);
	    break;
	case 2:
	    autoscale_graph(cg, -1);
	    break;
	}
	drawgraph();
    } else {
	errwin("No active sets!");
    }
}

/*
 * service the auticks button on the main panel
 */
void autoticks_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    default_axis(cg, g[cg].auto_type, X_AXIS);
    default_axis(cg, g[cg].auto_type, ZX_AXIS);
    default_axis(cg, g[cg].auto_type, Y_AXIS);
    default_axis(cg, g[cg].auto_type, ZY_AXIS);
    update_all(cg);
    drawgraph();
}

/*
 * set the message in the left footer
 */
void set_left_footer(char *s)
{
    Arg al;

    XmStringFree(statstring);
    statstring = XmStringCreateLtoR(s, charset);
    XtSetArg(al, XmNlabelString, statstring);
    XtSetValues(statlab, &al, 1);
    XmUpdateDisplay(statlab);
}

/*
 * clear the locator reference point
 */
void do_clear_point(Widget w, XtPointer client_data, XtPointer call_data)
{
    g[cg].pointset = FALSE;
    g[cg].pt_type = 0;
    g[cg].dsx = g[cg].dsy = 0.0;
}

/*
 * set visibility of the toolbars
 */
int toolbar_visible = 1;
int statusbar_visible = 1;
int locbar_visible = 1;

static Widget windowbarw[3];

static void set_view_items()
{
    if (statusbar_visible) {
	XmToggleButtonSetState(windowbarw[1], True, False);
	XtManageChild(frbot);
	XtVaSetValues(scrollw,
		      XmNbottomAttachment, XmATTACH_WIDGET,
		      XmNbottomWidget, frbot,
		      NULL);
	if (toolbar_visible) {
	    XtVaSetValues(frleft,
			  XmNbottomAttachment, XmATTACH_WIDGET,
			  XmNbottomWidget, frbot,
			  NULL);
	}
    } else {
	XmToggleButtonSetState(windowbarw[1], False, False);
	XtVaSetValues(scrollw,
		      XmNbottomAttachment, XmATTACH_FORM,
		      NULL);
	XtUnmanageChild(frbot);
	if (toolbar_visible) {
	    XtVaSetValues(frleft,
			  XmNbottomAttachment, XmATTACH_FORM,
			  NULL);
	}
    }
    if (toolbar_visible) {
	XmToggleButtonSetState(windowbarw[2], True, False);
	XtManageChild(frleft);
	if (statusbar_visible) {
	    XtVaSetValues(frleft,
			  XmNbottomAttachment, XmATTACH_WIDGET,
			  XmNbottomWidget, frbot,
			  NULL);
	}
	if (locbar_visible) {
	    XtVaSetValues(frleft,
			  XmNtopAttachment, XmATTACH_WIDGET,
			  XmNtopWidget, frtop,
			  NULL);
	}
	XtVaSetValues(scrollw,
		      XmNleftAttachment, XmATTACH_WIDGET,
		      XmNleftWidget, frleft,
		      NULL);
    } else {
	XmToggleButtonSetState(windowbarw[2], False, False);
	XtUnmanageChild(frleft);
	XtVaSetValues(scrollw,
		      XmNleftAttachment, XmATTACH_FORM,
		      NULL);
    }
    if (locbar_visible) {
	XmToggleButtonSetState(windowbarw[0], True, False);
	XtManageChild(frtop);
	XtVaSetValues(scrollw,
		      XmNtopAttachment, XmATTACH_WIDGET,
		      XmNtopWidget, frtop,
		      NULL);
	if (toolbar_visible) {
	    XtVaSetValues(frleft,
			  XmNtopAttachment, XmATTACH_WIDGET,
			  XmNtopWidget, frtop,
			  NULL);
	}
    } else {
	XmToggleButtonSetState(windowbarw[0], False, False);
	XtUnmanageChild(frtop);
	XtVaSetValues(scrollw,
		      XmNtopAttachment, XmATTACH_FORM,
		      NULL);
	if (toolbar_visible) {
	    XtVaSetValues(frleft,
			  XmNtopAttachment, XmATTACH_FORM,
			  NULL);
	}
    }
}

/*
 * called from the parser
 */
void set_toolbars(int bar, int onoff)
{
    switch (bar) {
    case TOOLBAR:
	toolbar_visible = onoff;
	break;
    case STATUSBAR:
	statusbar_visible = onoff;
	break;
    case LOCATORBAR:
	locbar_visible = onoff;
	break;
    }
    if (inwin) {
	set_view_items();
    }
}

/*
 * service routines for the View pulldown
 */
void set_statusbar(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (XmToggleButtonGetState(w)) {
	statusbar_visible = 1;
    } else {
	statusbar_visible = 0;
    }
    set_view_items();
}

void set_toolbar(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (XmToggleButtonGetState(w)) {
	toolbar_visible = 1;
    } else {
	toolbar_visible = 0;
    }
    set_view_items();
}

void set_locbar(Widget w, XtPointer client_data, XtPointer call_data)
{
    if (XmToggleButtonGetState(w)) {
	locbar_visible = 1;
    } else {
	locbar_visible = 0;
    }
    set_view_items();
}

/*
 * set the canvas size
 */
void set_canvas_size(int w, int h, int o)
{
    Dimension px, py;
    px = w;
    py = h;
    XtVaSetValues(canvas,
		  XmNwidth, px,
		  XmNheight, py,
		  NULL);
}

void get_default_canvas_size(int *w, int *h)
{
    Dimension ww, wh;
    Arg args;
    XtSetArg(args, XmNwidth, &ww);
    XtGetValues(scrollw, &args, 1);
    XtSetArg(args, XmNheight, &wh);
    XtGetValues(scrollw, &args, 1);
    *w = ww - 5;
    *h = wh - 5;
}

/*
 * service the Page pulldown item
 */
void set_page(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i;
    double wx1, wx2, wy1, wy2;
    Dimension px, py;
    int pageorient = (int) client_data;
    wx1 = DisplayWidth(disp, DefaultScreen(disp));
    wx2 = DisplayWidthMM(disp, DefaultScreen(disp));
    wy1 = DisplayHeight(disp, DefaultScreen(disp));
    wy2 = DisplayHeightMM(disp, DefaultScreen(disp));
    px = (Dimension) (wx1 / wx2 * (8.5 * 25.4));
    py = (Dimension) (wy1 / wy2 * (11.5 * 25.4));

    switch (pageorient) {
    case LANDSCAPE:
	page_layout = LANDSCAPE;
	XtVaSetValues(canvas,
		      XmNwidth, py,
		      XmNheight, px,
		      NULL);
	break;
    case PORTRAIT:
	page_layout = PORTRAIT;
	XtVaSetValues(canvas,
		      XmNwidth, px,
		      XmNheight, py,
		      NULL);
	break;
    case FIXED:
	page_layout = FIXED;
	XtVaSetValues(canvas,
		      XmNwidth, canvasw,
		      XmNheight, canvash,
		      NULL);
	break;
    case FREE:			/* falls through */
    default:
	page_layout = FREE;
	{
	    int w, h;
	    get_default_canvas_size(&w, &h);
	    px = w;
	    py = h;
	    XtVaSetValues(canvas,
			  XmNwidth, px,
			  XmNheight, py,
			  NULL);
	}
	break;
    }
    for (i = 0; i < 4; i++) {
	XmToggleButtonSetState(pagew[i], False, False);
    }
    XmToggleButtonSetState(pagew[get_pagelayout(pageorient)], True, False);
}

/*
 * get/set page layouts
 */
int get_pagelayout()
{
    switch (page_layout) {
    case FREE:
	return 0;
    case LANDSCAPE:
	return 1;
    case PORTRAIT:
	return 2;
    case FIXED:
	return 3;
    }
}

int set_pagelayout(int layout)
{
    switch (layout) {
    case LANDSCAPE:
	page_layout = LANDSCAPE;
	break;
    case PORTRAIT:
	page_layout = PORTRAIT;
	break;
    case FIXED:
	page_layout = FIXED;
	break;
    case FREE:			/* falls through */
    default:
	page_layout = FREE;
	break;
    }
    if (inwin) {
	set_page(NULL, (XtPointer) page_layout, NULL);
    }
    return page_layout;
}

/*
 */
void create_isolines_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    get_xybox_minmax(cg);
    create_isolines_popup("Isolines", &g[cg].isol, 0);
}

/*
 * create the main menubar
 */
static Widget CreateMenuBar(Widget parent)
{
    Widget menu_bar;
    Widget cascade;
    Widget button;

    menu_bar = XmCreateMenuBar(parent, "menu_bar", NULL, 0);
    menu_pane = XmCreatePulldownMenu(menu_bar, "Files menu", NULL, 0);

/* TODO as soon as a decent binary file format is found */
/*
   button = XtVaCreateManagedWidget("New", xmPushButtonWidgetClass, menu_pane,
   NULL);
   button = XtVaCreateManagedWidget("Open...", xmPushButtonWidgetClass, menu_pane,
   NULL);
   button = XtVaCreateManagedWidget("Save", xmPushButtonWidgetClass, menu_pane,
   NULL);
   button = XtVaCreateManagedWidget("Save as...", xmPushButtonWidgetClass, menu_pane,
   NULL);

   button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
   NULL);
 */

    button = XtVaCreateManagedWidget("Read sets...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_file_popup, NULL);

#ifdef HAVE_MFHDF
    button = XtVaCreateManagedWidget("Read netCDF/HDF...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_netcdfs_popup, NULL);
#else

#ifdef HAVE_NETCDF
    button = XtVaCreateManagedWidget("Read netCDF...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_netcdfs_popup, NULL);
#endif

#endif

    button = XtVaCreateManagedWidget("Read parameters...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_rparams_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Read block data...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_block_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Write sets...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_write_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Write parameters...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_wparam_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
				     NULL);
    button = XtVaCreateManagedWidget("Working directory...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_workingdir_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("1D visualization...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_vis_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Save all...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_saveall_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Clear all...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_CLEAR);

    button = XtVaCreateManagedWidget("sep3", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Print", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_PRINT);

    button = XtVaCreateManagedWidget("Printer setup...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_printer_setup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Commands...", xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("F4", charset),
				     XmNaccelerator, "<Key>F4:",
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) open_command, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep3", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Status...", xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("F5", charset),
				     XmNaccelerator, "<Key>F5:",
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_status_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Results...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_monitor_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep4", xmSeparatorWidgetClass, menu_pane,
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

/*
   menu_pane = XmCreatePulldownMenu(menu_bar, "Editmenu_pane", NULL, 0);
   cascade = XtVaCreateManagedWidget("Edit", xmCascadeButtonWidgetClass, menu_bar,
   XmNsubMenuId, menu_pane,
   XmNmnemonic, 'E',
   NULL);
 */

    menu_pane = XmCreatePulldownMenu(menu_bar, "Datamenu_pane", NULL, 0);

    button = XtVaCreateManagedWidget("Status...", xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("F5", charset),
				     XmNaccelerator, "<Key>F5:",
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_status_popup, (XtPointer) NULL);
    button = XtVaCreateManagedWidget("Results...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_monitor_frame, (XtPointer) NULL);

    XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane, NULL);

    button = XtVaCreateManagedWidget("Transformations...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_comp_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Set operations...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_setops_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Edit/create set...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_editp_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Region operations...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_region_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Point operations...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_points_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Block data...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_eblock_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Image data...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_image_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Hot links...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_hotlinks_popup, (XtPointer) NULL);

    cascade = XtVaCreateManagedWidget("Data", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'D',
				      NULL);

/* CCALMR specific options */
    menu_pane = XmCreatePulldownMenu(menu_bar, "ccalmr_pane", NULL, 0);

/*
    button = XtVaCreateManagedWidget("ADP/DB...", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) createADPdb, (XtPointer) NULL);
*/
    button = XtVaCreateManagedWidget("Start monitors", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) settimer, (XtPointer) NULL);
    button = XtVaCreateManagedWidget("Stop monitors", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) stoptimer, (XtPointer) NULL);
    button = XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane, NULL);
    button = XtVaCreateManagedWidget("Visualization", xmPushButtonWidgetClass, menu_pane, NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_vis_frame, (XtPointer) NULL);

    cascade = XtVaCreateManagedWidget("CCALMR", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'C',
				      NULL);

/* Graph menu */
    menu_pane = XmCreatePulldownMenu(menu_bar, "Graphmenu_pane", NULL, 0);
    button = XtVaCreateManagedWidget("Graph operations...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_graph_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("World scaling...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_world_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Viewport...", xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("^V", charset),
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_view_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Autoscale...",
				     xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("^A", charset),
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_autos_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Draw options...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_draw_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane, NULL);

    button = XtVaCreateManagedWidget("Flip XY", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_flipxy, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Invert X", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_invertx, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Invert Y", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_inverty, (XtPointer) NULL);
/*
   XtAddCallback(button, XmNvalueChangedCallback, do_inverty, (XtPointer) NULL);
   XtVaSetValues(button, XmNvisibleWhenOff, True, NULL);
 */

    button = XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Titles...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_label_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Tick labels/tick marks...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_ticks_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Frame...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_frame_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Symbols...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_symbols_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Legends...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_legend_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_isolines_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Velocities...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_velocityp_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Strings & things...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) define_objects_popup, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Time stamp...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_misc_frame, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Misc...",
				     xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_props_frame, (XtPointer) NULL);
    cascade = XtVaCreateManagedWidget("Graph", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'G',
				      NULL);

/* page layout */
    menu_pane = XmCreatePulldownMenu(menu_bar, "Pagemenu_pane", NULL, 0);
    XtVaSetValues(menu_pane, XmNradioBehavior, True, NULL);
    XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane, NULL);
    {
	pagew[0] = XmCreateToggleButton(menu_pane, "Free", NULL, 0);
	pagew[1] = XmCreateToggleButton(menu_pane, "Landscape", NULL, 0);
	pagew[2] = XmCreateToggleButton(menu_pane, "Portrait", NULL, 0);
	pagew[3] = XmCreateToggleButton(menu_pane, "Fixed", NULL, 0);
	XtAddCallback(pagew[0], XmNvalueChangedCallback, set_page, (XtPointer) FREE);
	XtVaSetValues(pagew[0], XmNvisibleWhenOff, True, NULL);
	XtAddCallback(pagew[1], XmNvalueChangedCallback, set_page, (XtPointer) LANDSCAPE);
	XtVaSetValues(pagew[1], XmNvisibleWhenOff, True, NULL);
	XtAddCallback(pagew[2], XmNvalueChangedCallback, set_page, (XtPointer) PORTRAIT);
	XtVaSetValues(pagew[2], XmNvisibleWhenOff, True, NULL);
	XtAddCallback(pagew[3], XmNvalueChangedCallback, set_page, (XtPointer) FIXED);
	XtVaSetValues(pagew[3], XmNvisibleWhenOff, True, NULL);
	XtManageChildren(pagew, 4);

	XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, menu_pane, NULL);

	button = XtVaCreateManagedWidget("Size...", xmPushButtonWidgetClass, menu_pane,
					 NULL);
	XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_page_frame,
		      (XtPointer) NULL);
    }
    cascade = XtVaCreateManagedWidget("Page", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'P',
				      NULL);

/* view */
    menu_pane = XmCreatePulldownMenu(menu_bar, "Viewmenu_pane", NULL, 0);
    windowbarw[0] = XtVaCreateManagedWidget("Locator bar", xmToggleButtonWidgetClass, menu_pane,
					    XmNindicatorOn, True,
					    XmNvisibleWhenOff, False,
					    NULL);
    XtAddCallback(windowbarw[0], XmNvalueChangedCallback,
		  (XtCallbackProc) set_locbar, (XtPointer) & frtop);
    windowbarw[1] = XtVaCreateManagedWidget("Status bar", xmToggleButtonWidgetClass, menu_pane,
					    XmNindicatorOn, True,
					    XmNvisibleWhenOff, False,
					    NULL);
    XtAddCallback(windowbarw[1], XmNvalueChangedCallback,
		  (XtCallbackProc) set_statusbar, (XtPointer) & frbot);
    windowbarw[2] = XtVaCreateManagedWidget("Tool bar", xmToggleButtonWidgetClass, menu_pane,
					    XmNindicatorOn, True,
					    XmNvisibleWhenOff, False,
					    NULL);
    XtAddCallback(windowbarw[2], XmNvalueChangedCallback,
		  (XtCallbackProc) set_toolbar, (XtPointer) & frleft);

/*
   button = XtVaCreateManagedWidget("Full screen (N/A)", xmToggleButtonWidgetClass, menu_pane,
   XmNindicatorOn, True,
   XmNvisibleWhenOff, False,
   NULL);
 */

    button = XtVaCreateManagedWidget("sep4", xmSeparatorWidgetClass, menu_pane,
				     NULL);

    button = XtVaCreateManagedWidget("Set locator fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) SEL_POINT);

    button = XtVaCreateManagedWidget("Clear locator fixed point", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) do_clear_point, (XtPointer) NULL);

    button = XtVaCreateManagedWidget("Locator props...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_locator_frame, (XtPointer) NULL);
    cascade = XtVaCreateManagedWidget("View", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'V',
				      NULL);

    menu_pane = XmCreatePulldownMenu(menu_bar, "Help menu pane", NULL, 0);
    button = XtVaCreateManagedWidget("Help...", xmPushButtonWidgetClass, menu_pane,
		   XmNacceleratorText, XmStringCreateLtoR("F1", charset),
				     XmNaccelerator, "<Key>F1:",
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_HELP);
    button = XtVaCreateManagedWidget("About...", xmPushButtonWidgetClass, menu_pane,
				     NULL);
    XtAddCallback(button, XmNactivateCallback, (XtCallbackProc) create_about_grtool, (XtPointer) NULL);
    cascade = XtVaCreateManagedWidget("Help", xmCascadeButtonWidgetClass, menu_bar,
				      XmNsubMenuId, menu_pane,
				      XmNmnemonic, 'H',
				      NULL);

    XtVaSetValues(menu_bar,
		  XmNmenuHelpWidget, cascade,
		  NULL);

    return (menu_bar);
}

/*
 * build the UI here
 */
void do_main_loop(void)
{
    Widget bt, rc3, rcleft, rctop, formbot;
    Pixmap icon;
    Arg al[10];
    int ac;
    int i;
    XSetWindowAttributes sw;

    main_frame = XtVaCreateManagedWidget("main", xmMainWindowWidgetClass, app_shell,
					 XmNshadowThickness, 0,
					 XmNwidth, 800,
					 XmNheight, 700,
/*
   XmNscrollBarDisplayPolicy, XmAS_NEEDED,
   XmNscrollingPolicy, XmAUTOMATIC,
 */
					 NULL);

    menu_bar = CreateMenuBar(main_frame);
    XtManageChild(menu_bar);

    form = XmCreateForm(main_frame, "form", NULL, 0);

    frleft = XtVaCreateManagedWidget("fr", xmFrameWidgetClass, form,
				     NULL);
    rcleft = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, frleft,
				     XmNorientation, XmVERTICAL,
				     XmNpacking, XmPACK_TIGHT,
				     XmNspacing, 0,
				     XmNentryBorder, 0,
				     XmNmarginWidth, 0,
				     XmNmarginHeight, 0,
				     NULL);
    frtop = XtVaCreateManagedWidget("frtop", xmFrameWidgetClass, form,
				    NULL);
    rctop = XtVaCreateManagedWidget("rctop", xmRowColumnWidgetClass, frtop,
				    XmNorientation, XmHORIZONTAL,
				    XmNpacking, XmPACK_TIGHT,
				    XmNspacing, 0,
				    XmNentryBorder, 0,
				    XmNmarginWidth, 0,
				    XmNmarginHeight, 0,
				    NULL);

    frbot = XtVaCreateManagedWidget("frbot", xmFrameWidgetClass, form, NULL);
    XtManageChild(frbot);
    /* formbot = XmCreateForm(frbot, "form", NULL, 0); */
    formbot = XmCreateRowColumn(frbot, "rc", NULL, 0);
    set_default_message(buf);
    statstring = XmStringCreateLtoR(buf, charset);
    statlab = XtVaCreateManagedWidget("statlab", xmLabelWidgetClass, formbot,
				      XmNlabelString, statstring,
				      XmNalignment, XmALIGNMENT_BEGINNING,
				      XmNrecomputeSize, True,
				      NULL);

    string = XmStringCreateLtoR("G0:[X, Y] =                                           ",
				charset);
    loclab = XtVaCreateManagedWidget("label Locate", xmLabelWidgetClass, rctop,
				     XmNlabelString, string,
				     XmNalignment, XmALIGNMENT_END,
				     XmNrecomputeSize, True,
				     NULL);
    XtManageChild(formbot);

    scrollw = XtVaCreateManagedWidget("scrollw",
				      xmScrolledWindowWidgetClass, form,
				XmNnavigationType, XmEXCLUSIVE_TAB_GROUP,
				      XmNscrollingPolicy, XmAUTOMATIC,
				      XmNvisualPolicy, XmVARIABLE,
				      NULL);

    canvas = XtVaCreateManagedWidget("canvas", xmDrawingAreaWidgetClass, scrollw,
				     XmNwidth, (Dimension) canvasw,
				     XmNheight, (Dimension) canvash,
				     XmNbackground,
				     WhitePixel(XtDisplay(main_frame),
				   DefaultScreen(XtDisplay(main_frame))),
				     NULL);
    XtAddCallback(canvas, XmNexposeCallback, (XtCallbackProc) refresh, NULL);

    XtAddEventHandler(canvas, EnterWindowMask
		      | LeaveWindowMask
		      | ButtonPressMask
		      | ButtonReleaseMask
		      | PointerMotionMask
		      | KeyPressMask
		      | ColormapChangeMask,
		      FALSE,
		      (XtEventHandler) my_proc, NULL);

    XtVaSetValues(frleft,
		  XmNtopAttachment, XmATTACH_WIDGET,
		  XmNtopWidget, frtop,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, frbot,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(frtop,
		  XmNtopAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  NULL);
    XtVaSetValues(scrollw,
		  XmNtopAttachment, XmATTACH_WIDGET,
		  XmNtopWidget, frtop,
		  XmNbottomAttachment, XmATTACH_WIDGET,
		  XmNbottomWidget, frbot,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_WIDGET,
		  XmNleftWidget, frleft,
		  NULL);
    XtVaSetValues(frbot,
		  XmNbottomAttachment, XmATTACH_FORM,
		  XmNrightAttachment, XmATTACH_FORM,
		  XmNleftAttachment, XmATTACH_FORM,
		  NULL);

    XtManageChild(form);

/*
   fr = XtVaCreateManagedWidget("fr", xmFrameWidgetClass, main_frame,
   XmNheight, 80,
   NULL);

   rc = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, fr,
   XmNorientation, XmHORIZONTAL,
   NULL);
 */

    init_pm();
    bt = XtVaCreateManagedWidget("Draw", xmPushButtonWidgetClass, rcleft,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) doforce_redraw, (XtPointer) NULL);

/* zoom and autoscale */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) ZOOM_1ST);

    bt = XtVaCreateManagedWidget("AS", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, autopm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_proc, (XtPointer) 0);

/* expand/shrink */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindshrink_proc, NULL);

    bt = XtVaCreateManagedWidget("z", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, shrinkpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindexpand_proc, NULL);

/*
 * scrolling buttons
 */
    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindleft_proc, NULL);
    bt = XtVaCreateManagedWidget("Right", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, rightpm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindright_proc, NULL);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwinddown_proc, NULL);
    bt = XtVaCreateManagedWidget("Up", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtVaSetValues(bt,
		  XmNlabelType, XmPIXMAP,
		  XmNlabelPixmap, uppm,
		  NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) gwindup_proc, NULL);

    XtVaCreateManagedWidget("sep1", xmSeparatorWidgetClass, rcleft,
			    NULL);

/*
   abort_button = XtVaCreateManagedWidget("Abort", xmPushButtonWidgetClass, rcleft,
   XmNsensitive, False,
   NULL);
   XtAddCallback(abort_button, XmNactivateCallback, (XtCallbackProc) do_abort, NULL);
   abort_win = XtWindow(abort_button);
 */

    bt = XtVaCreateManagedWidget("AutoT", xmPushButtonWidgetClass, rcleft,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoticks_proc, NULL);
    bt = XtVaCreateManagedWidget("AutoO", xmPushButtonWidgetClass, rcleft,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoon_proc, NULL);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    bt = XtVaCreateManagedWidget("ZX", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) ZOOMX_1ST);

    bt = XtVaCreateManagedWidget("ZY", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) set_actioncb, (XtPointer) ZOOMY_1ST);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
				  XmNorientation, XmHORIZONTAL,
				  XmNpacking, XmPACK_TIGHT,
				  XmNspacing, 0,
				  XmNentryBorder, 0,
				  XmNmarginWidth, 0,
				  XmNmarginHeight, 0,
				  NULL);
    bt = XtVaCreateManagedWidget("AX", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_proc, (XtPointer) 1);

    bt = XtVaCreateManagedWidget("AY", xmPushButtonWidgetClass, rc3,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_proc, (XtPointer) 2);

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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

    rc3 = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rcleft,
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
    stack_depth_item = XtVaCreateManagedWidget("stackdepth", xmLabelWidgetClass, rcleft,
					       XmNlabelString, sdstring,
					       NULL);

    cystring = XmStringCreateLtoR("CW:0 ", charset);
    curw_item = XtVaCreateManagedWidget("curworld", xmLabelWidgetClass, rcleft,
					XmNlabelString, cystring,
					NULL);

    bt = XtVaCreateManagedWidget("Exit", xmPushButtonWidgetClass, rcleft,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) MenuCB, (XtPointer) MENU_EXIT);

/*
 * initialize the tool bars
 */
    set_view_items();

    XmMainWindowSetAreas(main_frame, menu_bar, NULL, NULL, NULL, form);
    XtRealizeWidget(app_shell);

    xwin = XtWindow(canvas);
    disp = XtDisplay(canvas);

    set_page(NULL, (XtPointer) page_layout, NULL);

    sw.backing_store = Always;
    XChangeWindowAttributes(disp, xwin, CWBackingStore, &sw);

    XtAddCallback(canvas, XmNresizeCallback, (XtCallbackProc) refresh, (XtPointer) 1);

    if (named_pipe) {
/*
	set_pipetimer();
*/
    }
    xlibinitcmap();
    if (use_colors > 2) {
	ac = 0;
	XtSetArg(al[ac], XmNcolormap, mycmap);
	ac++;
	XtSetValues(canvas, al, ac);
	XSetWindowColormap(disp, xwin, mycmap);

/* moved to initialize_screen() above
   for (i = 0; i < maxgraph; i++) {
   if (g[i].parmsread != TRUE) {
   setdefaultcolors(i);
   }
   }
 */
    }
    gc = DefaultGC(disp, DefaultScreen(disp));
    gc_val.foreground = WhitePixel(disp, DefaultScreen(disp));
    gc_val.foreground = BlackPixel(disp, DefaultScreen(disp)) ^ WhitePixel(disp, DefaultScreen(disp));
    if (invert) {
	gc_val.function = GXinvert;
    } else {
	gc_val.function = GXxor;
    }
    gcxor = XCreateGC(disp, xwin, GCFunction | GCForeground, &gc_val);
    gc_val.foreground = WhitePixel(disp, DefaultScreen(disp));
    gc_val.function = GXcopy;
    gcclr = XCreateGC(disp, xwin, GCFunction | GCForeground, &gc_val);

    icon = XCreateBitmapFromData(XtDisplay(app_shell),
				 DefaultRootWindow(XtDisplay(app_shell)),
				 (char *) xmgr_icon_bits, xmgr_icon_width,
				 xmgr_icon_height);
    XtVaSetValues(app_shell,
		  XtNiconPixmap, icon,
		  XtNiconMask, icon,
		  NULL);

    /*
     * initialize cursors
     */
    init_cursors();

    /*
     * if an image was placed on the command line, read it in
     */
    if (readimage) {
	read_image(image_filename);
    }
    inwin = 1;
    log_results("Startup");
    inwin = 0;

    /*
     * set the title to the working directory
     */
    set_title(NULL);

    /*
     * Process events.
     */
    XtAppMainLoop(app_con);
}

/*
 * initialize pixmaps for buttons on front
 */
static void init_pm(void)
{
    Display *disp = XtDisplay(app_shell);
    Window cwin = RootWindowOfScreen(XtScreen(app_shell));
/*
   GC gc = DefaultGC(disp, DefaultScreen(disp));
 */
    GC gc;
    Pixmap ptmp;
    Pixel fg, bg;

    XtVaGetValues(menu_bar,
		  XmNforeground, &fg,
		  XmNbackground, &bg,
		  NULL);

    gc = XCreateGC(disp, cwin, 0, NULL);
    XSetForeground(disp, gc, fg);
    XSetBackground(disp, gc, bg);

    zoompm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) zoom_bits, 16, 16);
    XCopyPlane(disp, ptmp, zoompm, gc, 0, 0, 16, 16, 0, 0, 1);
    autopm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) auto_bits, 16, 16);
    XCopyPlane(disp, ptmp, autopm, gc, 0, 0, 16, 16, 0, 0, 1);
    shrinkpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) shrink_bits, 16, 16);
    XCopyPlane(disp, ptmp, shrinkpm, gc, 0, 0, 16, 16, 0, 0, 1);
    expandpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) expand_bits, 16, 16);
    XCopyPlane(disp, ptmp, expandpm, gc, 0, 0, 16, 16, 0, 0, 1);
    rightpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) right_bits, 16, 16);
    XCopyPlane(disp, ptmp, rightpm, gc, 0, 0, 16, 16, 0, 0, 1);
    leftpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) left_bits, 16, 16);
    XCopyPlane(disp, ptmp, leftpm, gc, 0, 0, 16, 16, 0, 0, 1);
    uppm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) up_bits, 16, 16);
    XCopyPlane(disp, ptmp, uppm, gc, 0, 0, 16, 16, 0, 0, 1);
    downpm = XCreatePixmap(disp, cwin, 16, 16, DisplayPlanes(disp, DefaultScreen(disp)));
    ptmp = XCreateBitmapFromData(disp, cwin, (char *) down_bits, 16, 16);
    XCopyPlane(disp, ptmp, downpm, gc, 0, 0, 16, 16, 0, 0, 1);
}
