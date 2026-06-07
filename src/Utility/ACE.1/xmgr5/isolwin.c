/*
 *
 * set the parameters for isolines
 *
 */

#ifndef lint
static char RCSid[] = "$Id: isolwin.c,v 1.2 2004/07/07 04:00:27 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/ToggleB.h>
 
#include "globals.h"
#include "motifinc.h"

extern Widget app_shell;

void create_isolines_popup(char *title, Isolparms * cd, int cm);
static void autoscale_isolines_proc(void);

static int isolinespopup = 0;

char isoltitle[256] = "Isolines";

typedef struct _Isol_ui {
    int type;
    int nitem;
    int stype;
    char title[256];
    Widget top;
    Widget *n_item;
    Widget *precision_item;
    Widget start_item;
    Widget step_item;
    Widget minmax_item;
    Widget *line_item;
    Widget *layout_item;
    Widget *val_item;
    Widget *linepat_item;
    Widget *type_item;
    int ndefs;
    Widget *defs_item;
    Widget *defslabels_item;
    Widget *defs_frame;
    Widget *ssdefs_frame;
} Isol_ui;

static Widget isolines_frame;
static Widget *isol_n_item;
static Widget *isol_precision_item;
static Widget isol_start_item;
static Widget isol_step_item;
static Widget isol_minmax_item;
static Widget *isol_line_item;
static Widget *isol_layout_item;
static Widget *isol_val_item;
static Widget *isol_linepat_item;
static Widget *isol_type_item;
static Widget isol_defs_item[20];
static Widget isol_defslabels_item[20];
static Widget isol_defs_frame;
static Widget isol_ssdefs_frame;

static Widget legend_toggle;

char *xv_getstr(Widget w);
void update_isolines(void);
void create_cedit_popup(void);
Widget CreateTextItem3();

void linep_accept_proc(void);
void update_linep(void);
void create_linep_frame(void);

static Isolparms *iptr;

void isolleg_place_proc(void)
{
    set_action(0);
    set_action(PLACE_ISOLINES_LEGEND);
}

static void define_isolines_proc(void)
{
    char buf[80];
    int n, i;

    curip.nisol = GetChoice(isol_n_item) + 1;
    if (curip.nisol <= 0) {
	curip.nisol = 1;
    }
    if (curip.nisol > 20) {
	curip.nisol = 20;
    }
    curip.isoltype = (int) GetChoice(isol_type_item);
    if (curip.isoltype == 0) {
	curip.cis[0] = atof(xv_getstr(isol_start_item));
	curip.cint = atof(xv_getstr(isol_step_item));
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = curip.cis[0] + i * curip.cint;
	}
    } else {
	for (i = 0; i < curip.nisol; i++) {
	    curip.cis[i] = atof(xv_getstr(isol_defs_item[i]));
	}
    }
    curip.type = (int) GetChoice(isol_line_item);
    curip.layout = (int) GetChoice(isol_layout_item) ? HORIZONTAL : VERTICAL;
/*
    curip.valflag = (int) GetChoice(isol_val_item) ? UP : DOWN;
*/
    curip.p.prec = (int) GetChoice(isol_precision_item);
    curip.lactive = XmToggleButtonGetState(legend_toggle) ? ON : OFF;
    *iptr = curip;
    update_isolines();
    update_linep();
    /* XtUnmanageChild(isolines_frame); */
}

static void set_isolines_type_proc(void)
{
    int i;
    char s[80];
    Arg a;

    curip.type = (int) GetChoice(isol_type_item);
    curip.p.prec = (int) GetChoice(isol_precision_item);
    if (curip.isoltype == 0) {
	sprintf(s, "%.*lf", curip.p.prec, curip.cmin);
	xv_setstr(isol_start_item, s);
	curip.cis[0] = curip.cmin;
	curip.cint = (curip.cmax - curip.cmin) / curip.nisol;
	sprintf(s, "%.*lf", curip.p.prec, curip.cint);
	xv_setstr(isol_step_item, s);
    } else {
	for (i = 0; i < 20; i++) {
	    sprintf(s, "%.*lf", curip.p.prec, curip.cis[i]);
	    xv_setstr(isol_defs_item[i], s);
	}
    }
}

void update_isolines_items(int gno)
{
    get_xybox_minmax(gno);
    if (isolinespopup) {
        create_isolines_popup("Isolines", &g[gno].isol, 0);
    }
}

void update_isolines(void)
{
    int i;
    char buf[256];

    if (isolines_frame) {
	sprintf(buf, "%lf %lf", curip.cmin, curip.cmax);
	xv_setstr(isol_minmax_item, buf);
	if (curip.nisol <= 0) {
	    curip.nisol = 1;
	}
	if (curip.nisol > 20) {
	    curip.nisol = 20;
	}
	SetChoice(isol_n_item, curip.nisol - 1);

	SetChoice(isol_type_item, curip.isoltype);
	if (curip.isoltype == 0) {
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cis[0]);
	    xv_setstr(isol_start_item, buf);
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cint);
	    xv_setstr(isol_step_item, buf);
	    for (i = 0; i < curip.nisol; i++) {
		sprintf(buf, "%.*lf", curip.p.prec, curip.cis[0] + i * curip.cint);
		xv_setstr(isol_defs_item[i], buf);
	    }
	    XtSetSensitive(isol_defs_frame, False);
	    XtSetSensitive(isol_ssdefs_frame, True);
	} else {
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cis[0]);
	    xv_setstr(isol_start_item, buf);
	    sprintf(buf, "%.*lf", curip.p.prec, curip.cint);
	    xv_setstr(isol_step_item, buf);
	    for (i = 0; i < curip.nisol; i++) {
		sprintf(buf, "%.*lf", curip.p.prec, curip.cis[i]);
		xv_setstr(isol_defs_item[i], buf);
	    }
	    XtSetSensitive(isol_defs_frame, True);
	    XtSetSensitive(isol_ssdefs_frame, False);
	}
	SetChoice(isol_line_item, curip.type);
	SetChoice(isol_precision_item, curip.p.prec);
	XmToggleButtonSetState(legend_toggle, curip.lactive == ON, False);
	SetChoice(isol_layout_item, curip.layout == HORIZONTAL ? 1 : 0);
/*
	SetChoice(isol_val_item, curip.valflag == UP ? 1 : 0);
*/
	update_linep();
    }
}

static void set_curtype_proc(Widget w, int cd)
{
    curip.isoltype = cd;
    update_isolines();
}

static void isolines_done_proc(void)
{
    isolinespopup = 0;
    XtUnmanageChild(isolines_frame);
}

static void autoscale_isolines_proc(void)
{
    int i;
    if (isolines_frame) {
	curip.p.prec = (int) GetChoice(isol_precision_item);
    }
    default_isolines(cg, &curip);
    update_isolines();
}

void create_isolines_popup(char *title, Isolparms * cd, int cm)
{
    int i;
    char s[80];
    Widget panel, bt, rc, rc1, fr, lab;
    XmString xms;
    extern long colors[];

    curip = *cd;
    iptr = cd;
    if (title != NULL) {
	strcpy(isoltitle, title);
    } else {
	strcpy(isoltitle, "Isolines");
    }
    if (isolines_frame) {
	update_isolines();
	XtVaSetValues(isolines_frame,
		      XmNtitle, isoltitle,
		      NULL);
	for (i = 0; i < 20; i++) {

	    if (cm == 0) {	/* bathymetry */
		XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 7], NULL);
		sprintf(s, "#%2d (C%2d): ", i + 1, i + 7);
	    } else {
		XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
		sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	    }
	    xms = XmStringCreateLtoR(s, charset);
	    XtVaSetValues(isol_defslabels_item[i], XmNlabelString, xms, NULL);
	    XmStringFree(xms);
	}
	XtRaise(isolines_frame);
	isolinespopup = 1;
	return;
    }
    isolines_frame = XmCreateDialogShell(app_shell, "Isolines", NULL, 0);
    handle_close(isolines_frame);
    panel = XmCreateRowColumn(isolines_frame, "rc", NULL, 0);
    isol_minmax_item = CreateTextItem2(panel,
				       20, "Limits (status only) ");
    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_line_item = CreatePanelChoice(rc, "Draw isolines:",
				     4, "As lines", "Filled", "Both", 0, 0);
    isol_precision_item = CreatePanelChoice0(rc, "Precision:",
					     2, 11,
			"0", "1", "2", "3", "4", "5", "6", "7", "8", "9", 0,
					     0);
    XtManageChild(rc);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_n_item = CreatePanelChoice0(rc, "# of isolines:",
				     4, 21,
				"1", "2", "3", "4", "5", "6", "7", "8", "9",
				   "10", "11", "12", "13", "14", "15", "16",
				   "17", "18", "19", "20",
				     0, 0);

    isol_type_item = CreatePanelChoice(rc, "Spacing by:",
					3,
					"Start, step",
					"Specified",
					0, 0);
    XtAddCallback(isol_type_item[2], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, 0);
    XtAddCallback(isol_type_item[3], XmNactivateCallback, (XtCallbackProc) set_curtype_proc, (XtPointer) 1);
    XtManageChild(rc);

    isol_ssdefs_frame = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_ssdefs_frame, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    isol_start_item = CreateTextItem2(rc, 8, "Start:");
    isol_step_item = CreateTextItem2(rc, 8, "Step:");
    XtManageChild(rc);
    XtManageChild(isol_ssdefs_frame);

    isol_defs_frame = XmCreateFrame(panel, "fr", NULL, 0);
    rc = XmCreateRowColumn(isol_defs_frame, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = 0; i < 10; i++) {
	if (cm == 0) {		/* bathymetry */
	    sprintf(s, "#%2d (C%2d): ", i + 1, i + 7);
	    isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 7], NULL);
	} else {
	    sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	    isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
	}
    }
    XtManageChild(rc1);
    rc1 = XmCreateRowColumn(rc, "rc1", NULL, 0);
    for (i = 10; i < 20; i++) {
	sprintf(s, "#%2d (C%2d): ", i + 1, i + 32);
	isol_defs_item[i] = CreateTextItem3(rc1, 10, s, &isol_defslabels_item[i]);
	if (cm == 0) {		/* bathymetry */
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 7], NULL);
	} else {
	    XtVaSetValues(isol_defslabels_item[i], XtNbackground, colors[i + 32], NULL);
	}
    }
    XtManageChild(rc1);
    XtManageChild(rc);
    XtManageChild(isol_defs_frame);

    legend_toggle = XtVaCreateManagedWidget("Display legend", xmToggleButtonWidgetClass,
					    panel,
					    NULL);
    isol_layout_item = CreatePanelChoice(panel, "Legend layout:",
					  3,
					  "Vertical",
					  "Horizontal",
					  0, 0);
    isol_val_item = CreatePanelChoice(panel, "Legend values:",
				       3,
				       "Ascending",
				       "Descending",
				       0, 0);

    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) define_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Autoscale", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) autoscale_isolines_proc, NULL);

    bt = XtVaCreateManagedWidget("Place legend", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolleg_place_proc, NULL);

/*
    bt = XtVaCreateManagedWidget("Colors...", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_cedit_popup, NULL);
*/

    bt = XtVaCreateManagedWidget("Line props...", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_linep_frame, NULL);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) isolines_done_proc, NULL);
    XtManageChild(rc);
    XtManageChild(panel);
    XtRaise(isolines_frame);
    XtVaSetValues(isolines_frame,
		  XmNtitle, isoltitle,
		  NULL);
    update_isolines();
    update_linep();
    if (curip.nisol == 0) {
	autoscale_isolines_proc();
    }
    isolinespopup = 1;
}

static Widget *linep_color_item[20];
static Widget *linep_linew_item[20];
static Widget *linep_lines_item[20];
static Widget linep_frame;

void create_linep_frame(void)
{
    Widget panel, wbut, fr, bb, sep, rc, rc2;
    int i;
    if (!linep_frame) {
	linep_frame = XmCreateDialogShell(app_shell, "Line properties", NULL, 0);
	panel = XmCreateRowColumn(linep_frame, "lineprc", NULL, 0);
	for (i = 0; i < 20; i++) {
	    rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	    linep_color_item[i] = CreateColorChoice(rc, "Line color:", 1);
	    linep_linew_item[i] = CreatePanelChoice0(rc,
						     "Width: ",
						     3, 10,
				"1", "2", "3", "4", "5", "6", "7", "8", "9",
						     0, 0);
	    linep_lines_item[i] = CreatePanelChoice(rc, "Style:",
						     6,
						     "Solid line",
						     "Dotted line",
						     "Dashed line",
						     "Long Dashed",
						     "Dot-dashed",
						     NULL,
						     NULL);
	    XtManageChild(rc);
	}

	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) linep_accept_proc, NULL);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, linep_frame);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    update_linep();
    XtRaise(linep_frame);
}

void update_linep(void)
{
    char buf[256];
    int i;
    if (linep_frame) {
	for (i = 0; i < 20; i++) {
	    SetChoice(linep_color_item[i], curip.color[i]);
	    SetChoice(linep_linew_item[i], curip.linew[i] - 1);
	    SetChoice(linep_lines_item[i], curip.lines[i] - 1);
	}
    }
}

void linep_accept_proc(void)
{
    char buf[256];
    int i;
    for (i = 0; i < 20; i++) {
	curip.color[i] = GetChoice(linep_color_item[i]);
	curip.linew[i] = GetChoice(linep_linew_item[i]) + 1;
	curip.lines[i] = GetChoice(linep_lines_item[i]) + 1;
    }
    *iptr = curip;
    XtUnmanageChild(linep_frame);
}

static Widget vel_frame;
static Widget vscale_item;
static Widget *vscale_units_item;
static Widget vscale_val_item;
static Widget vscale_length_item;
static Widget *vscale_arrow_item;
static Widget vscale_string_item;
static Widget vscale_llength_item;
static Widget *vscale_color_item;
static Widget *vscale_font_item;
static Widget *vscale_loctype_item;
static void define_vplace_proc(void);
void update_scaling(void);
static void vel_accept_proc(void);

/*
 * define the velocity scale for this graph and properties of the
 * velocity scale legend
 */
void create_velocityp_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    Arg wargs[8];
    int i;
    Widget vel_panel, wbut, vfr2, fr, rc, rb, bb, sep;

    if (!vel_frame) {
	vel_frame = XmCreateDialogShell(app_shell, "Scaling", NULL, 0);
	handle_close(vel_frame);
	vel_panel = XmCreateRowColumn(vel_frame, "rc", NULL, 0);

	vfr2 = XmCreateFrame(vel_panel, "frame", NULL, 0);
	bb = XmCreateRowColumn(vfr2, "rc", NULL, 0);
	vscale_length_item = CreateTextItem2(bb, 15, "Scale factor = ");
	vscale_llength_item = CreateTextItem2(bb, 15, "Velocity scale legend length: ");
	vscale_string_item = CreateTextItem2(bb, 15, "Velocity scale label: ");
	vscale_color_item = CreateColorChoice(bb, "Velocity scale legend color: ", 1);
	vscale_arrow_item = CreatePanelChoice(bb, "Arrow type:",
					      4,
						"->", "-|>", "-|> filled",
					      0, 0);
	vscale_font_item = CreatePanelChoice(bb, "Font:",
					      11,
				"Times-Roman", "Times-Bold", "Times-Italic",
					    "Times-BoldItalic", "Helvetica",
				      "Helvetica-Bold", "Helvetica-Oblique",
				 "Helvetica-BoldOblique", "Greek", "Symbol",
					      0, 0);
	vscale_item = XtVaCreateManagedWidget("Display velocity scale legend",
					      xmToggleButtonWidgetClass, bb,
					      NULL);
	XtManageChild(bb);
	XtManageChild(vfr2);

	rc = XmCreateRowColumn(vel_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) vel_accept_proc, NULL);

	wbut = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) define_vplace_proc, NULL);

	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, vel_frame);
	XtManageChild(rc);
	XtManageChild(vel_panel);
    }
    XtRaise(vel_frame);
    update_scaling();
}

/*
 * place the velocity scale
 */
static void define_vplace_proc(void)
{
    set_action(0);
    set_action(PLACE_VSCALE);
}

/*
 * update the velocity scale popup
 */
void update_scaling(void)
{
    char buf[256];
    if (vel_frame) {
	sprintf(buf, "%lf\n", g[cg].vp.vscale);
	xv_setstr(vscale_length_item, buf);
	sprintf(buf, "%.2lf\n", g[cg].vp.userlength);
	xv_setstr(vscale_llength_item, buf);
	xv_setstr(vscale_string_item, g[cg].vp.vstr.s);
	XmToggleButtonSetState(vscale_item, g[cg].vp.active == ON, False);
	SetChoice(vscale_color_item, g[cg].vp.color);
	SetChoice(vscale_arrow_item, g[cg].vp.arrowtype);
    }
}

/*
 * accept the current state of the velocity popup
 */
static void vel_accept_proc(void)
{
    char buf[256];
    g[cg].vp.vscale = atof(xv_getstr(vscale_length_item));
    g[cg].vp.userlength = atof(xv_getstr(vscale_llength_item));
    g[cg].vp.active = XmToggleButtonGetState(vscale_item) == True ? ON : OFF;
    g[cg].vp.color = GetChoice(vscale_color_item);
    g[cg].vp.arrowtype = GetChoice(vscale_arrow_item);
    set_plotstr_string(&g[cg].vp.vstr, (char *) xv_getstr(vscale_string_item));
    drawgraph();
}
