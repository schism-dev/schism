/* $Id: compwin.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * transformations, curve fitting, etc.
 *
 * formerly, this was all one big popup, now it is several.
 * All are created as needed
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/LabelG.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/ToggleB.h>
#include <Xm/Text.h>

#include "globals.h"
#include "motifinc.h"

static Widget but1[2];
static Widget but2[3];

void create_ntiles_frame(Widget w, XtPointer client_data, XtPointer call_data);
void create_geom_frame(Widget w, XtPointer client_data, XtPointer call_data);

static void do_compute_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_load_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_compute_proc2(Widget w, XtPointer client_data, XtPointer call_data);
static void do_digfilter_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_linearc_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_xcor_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_spline_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_int_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_differ_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_seasonal_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_interp_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_regress_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_runavg_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_fourier_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_fft_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_window_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_histo_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_sample_proc(Widget w, XtPointer client_data, XtPointer call_data);

/*
 * Widget declarations
 */
static Widget comp_frame;
static Widget comp_panel;

static Widget compute_formula_item;
static Widget *compute_set_item;
static Widget *compute_load_item;
static Widget *compute_loadgraph_item;
static Widget *compute_region_item;
static Widget compute_rinvert_item;
static Widget *load_set_item;
static Widget load_start_item;
static Widget load_step_item;
static Widget *load_to_item;
static Widget histo_binw_item;
static Widget histo_hxmin_item;
static Widget histo_hxmax_item;
static Widget *histo_type_item;
static Widget *histo_set_item;
static Widget *histo_to_item;
static Widget *histo_graph_item;
static Widget *toggle_set_fourier_item;
static Widget *toggle_load_fourier_item;
static Widget *toggle_window_fourier_item;
static Widget *toggle_loadx_fourier_item;
static Widget *toggle_inv_fourier_item;
static Widget *toggle_type_fourier_item;
static Widget toggle_ravglen_item;
static Widget *toggle_run_item;
static Widget *toggle_set_runavg_item;
static Widget *runavg_region_item;
static Widget runavg_rinvert_item;
static Widget *toggle_degree_item;
static Widget toggle_zero_item;
static Widget *toggle_resid_item;
static Widget *toggle_set_regress_item;
static Widget *regress_region_item;
static Widget regress_rinvert_item;
static Widget regress_start_item;
static Widget regress_stop_item;
static Widget regress_step_item;
static Widget *toggle_set_differ_item;
static Widget *toggle_differ_type_item;
static Widget *toggle_set_int_item;
static Widget *toggle_int_type_item;
static Widget int_sum_item;
static Widget *toggle_set_seasonal_item;
static Widget seasonal_period_item;
static Widget *interp_from_item;
static Widget *interp_to_item;
static Widget interp_type_item;
static Widget *spline_item;
static Widget spline_start_item;
static Widget spline_stop_item;
static Widget spline_step_item;
static Widget *toggle_set_xcor1_item;
static Widget *toggle_set_xcor2_item;
static Widget *toggle_xcor_type_item;
static Widget xcor_lag_item;
static Widget *toggle_set_sample_item;
static Widget sample_start_item;
static Widget sample_step_item;
static Widget *sample_type_item;
static Widget sample_expr_item;
static Widget *toggle_set_digf1_item;
static Widget *toggle_set_digf2_item;
static Widget *toggle_set_linc1_item;
static Widget *toggle_set_linc2_item;
static Widget compute_formulax_item;
static Widget compute_formulay_item;
static Widget load_start2_item;
static Widget load_stop2_item;
static Widget load_npts_item;
static Widget *load_to2_item;

/*
 * Create the comp Widget
 */
void create_comp_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget wbut;
    int x, y;

    set_wait_cursor();
    if (comp_frame == NULL) {
	XmGetPos(app_shell, 0, &x, &y);
	comp_frame = XmCreateDialogShell(app_shell, "Transformations", NULL, 0);
	handle_close(comp_frame);
	XtVaSetValues(comp_frame, XmNx, x, XmNy, y, NULL);
	comp_panel = XmCreateRowColumn(comp_frame, "comp_rc", NULL, 0);

	wbut = XtVaCreateManagedWidget("Evaluate expression...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_eval_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Load values...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_load_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Load & evaluate...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_leval_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Histograms...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_histo_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Fourier transforms...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_fourier_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Running averages...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_run_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Regression...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_reg_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Differences...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_diff_frame, (XtPointer) NULL);
	wbut = XtVaCreateManagedWidget("Seasonal differences...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_seasonal_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Integration...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_int_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Cross/auto correlation...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_xcor_frame, (XtPointer) NULL);

/*
	wbut = XtVaCreateManagedWidget("Interpolation...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_interp_frame, (XtPointer) NULL);
*/

	wbut = XtVaCreateManagedWidget("Splines...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_spline_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Sample points...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_samp_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Digital filter...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_digf_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Linear convolution...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_lconv_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Geometric transformations...", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_geom_frame, (XtPointer) NULL);

	wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, comp_panel,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) comp_frame);
	XtManageChild(comp_panel);
    }
    XtRaise(comp_frame);
    unset_wait_cursor();
}

void do_pick_compose(Widget w, XtPointer client_data, XtPointer call_data)
{
    set_action(0);
    set_action((int) client_data);
}

void create_eval_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;
    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Pick";
	label2[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Evaluate expression", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	compute_set_item = CreateSetChoice(dialog, "Apply to set:", maxplot, 1);

	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	compute_load_item = CreatePanelChoice(rc,
					      "Result to:", 3,
					 "Same set", "New set", NULL, 0);
	compute_loadgraph_item = CreateGraphChoice(rc, "In graph: ", maxgraph, 1);
	XtManageChild(rc);

/*
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    compute_region_item = CreatePanelChoice(rc,
					    "Restrictions:",
					    9,
					    "None",
					    "Region 0",
					    "Region 1",
					    "Region 2",
					    "Region 3",
					    "Region 4",
					    "Inside graph",
					    "Outside graph",
					    0,
					    0);
    compute_rinvert_item = XmCreateToggleButton(rc, "Invert region", NULL, 0);
    XtManageChild(compute_rinvert_item);
    XtManageChild(rc);
*/

	compute_formula_item = CreateTextItem2(dialog, 30, "Formula:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 3, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_compute_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) do_pick_compose, (XtPointer) PICK_EXPR);
	XtAddCallback(but2[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* load a set */

void create_load_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "Load values", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 4,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set: ", xmLabelWidgetClass, rc, NULL);
	load_set_item = CreateSetChoice(rc, " ", maxplot, 1);

	XtVaCreateManagedWidget("Load to: ", xmLabelWidgetClass, rc, NULL);
	load_to_item = CreatePanelChoice(rc,
					 " ",
					 7,
					 "Set X",
					 "Set Y",
					 "Scratch A",
					 "Scratch B",
					 "Scratch C",
					 "Scratch D", 0,
					 0);
	load_start_item = CreateTextItem4(rc, 10, "Start:");
	load_step_item = CreateTextItem4(rc, 10, "Step:");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_load_proc, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* histograms */

void create_histo_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc, rc1, form;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Histograms", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 4,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set: ", xmLabelWidgetClass, rc, NULL);
	histo_set_item = CreateSetChoice(rc, " ", maxplot, 0);

	XtVaCreateManagedWidget("Bin width: ", xmLabelWidgetClass, rc, NULL);
	histo_binw_item = XtVaCreateManagedWidget("binwidth", xmTextWidgetClass, rc, NULL);
	XtVaSetValues(histo_binw_item, XmNcolumns, 10, NULL);
	XtVaCreateManagedWidget("Start value: ", xmLabelWidgetClass, rc, NULL);
	histo_hxmin_item = XtVaCreateManagedWidget("xmin", xmTextWidgetClass, rc, NULL);
	XtVaSetValues(histo_hxmin_item, XmNcolumns, 10, NULL);
	XtVaCreateManagedWidget("Ending value: ", xmLabelWidgetClass, rc, NULL);
	histo_hxmax_item = XtVaCreateManagedWidget("xmax", xmTextWidgetClass, rc, NULL);
	XtVaSetValues(histo_hxmax_item, XmNcolumns, 10, NULL);
	XtManageChild(rc);

	histo_type_item = CreatePanelChoice(dialog, "Load: ",
					    3,
					    "Histogram",
					    "Cumulative histogram",
					    0,
					    0);

	rc1 = XmCreateRowColumn(dialog, "rc1", NULL, 0);
	XtVaSetValues(rc1, XmNorientation, XmHORIZONTAL, NULL);
	XtVaCreateManagedWidget("Result to set:", xmLabelWidgetClass, rc1, NULL);
	histo_to_item = CreateSetChoice(rc1, "", maxplot, 4);
	histo_graph_item = CreateGraphChoice(rc1, "In graph:", maxgraph, 1);
	XtManageChild(rc1);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_histo_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* DFTs */

void create_fourier_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;
    Widget buts[5];

    set_wait_cursor();
    if (top == NULL) {
	char *l[5];
	l[0] = "DFT";
	l[1] = "FFT";
	l[2] = "Window only";
	l[3] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Fourier transforms", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 6,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set: ", xmLabelWidgetClass, rc, NULL);
	toggle_set_fourier_item = CreateSetChoice(rc, " ", maxplot, 0);

	XtVaCreateManagedWidget("Data window: ", xmLabelWidgetClass, rc, NULL);
	toggle_window_fourier_item = CreatePanelChoice(rc,
						       " ",
						       8,
						    "None (Rectangular)",
						       "Triangular",
						       "Hanning",
						       "Welch",
						       "Hamming",
						       "Blackman",
						       "Parzen",
						       NULL,
						       NULL);

	XtVaCreateManagedWidget("Load result as: ", xmLabelWidgetClass, rc, NULL);

	toggle_load_fourier_item = CreatePanelChoice(rc,
						     " ",
						     4,
						     "Magnitude",
						     "Phase",
						     "Coefficients",
						     0,
						     0);

	XtVaCreateManagedWidget("Let result X = ", xmLabelWidgetClass, rc, NULL);
	toggle_loadx_fourier_item = CreatePanelChoice(rc,
						      " ",
						      4,
						      "Index",
						      "Frequency",
						      "Period",
						      0,
						      0);

	XtVaCreateManagedWidget("Perform: ", xmLabelWidgetClass, rc, NULL);
	toggle_inv_fourier_item = CreatePanelChoice(rc,
						    " ",
						    3,
						    "Transform",
						    "Inverse transform",
						    0,
						    0);

	XtVaCreateManagedWidget("Data is: ", xmLabelWidgetClass, rc, NULL);
	toggle_type_fourier_item = CreatePanelChoice(rc,
						     " ",
						     3,
						     "Real",
						     "Complex",
						     0,
						     0);
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);
	CreateCommandButtons(dialog, 4, buts, l);
	XtAddCallback(buts[0], XmNactivateCallback, (XtCallbackProc) do_fourier_proc, (XtPointer) top);
	XtAddCallback(buts[1], XmNactivateCallback, (XtCallbackProc) do_fft_proc, (XtPointer) top);
	XtAddCallback(buts[2], XmNactivateCallback, (XtCallbackProc) do_window_proc, (XtPointer) top);
	XtAddCallback(buts[3], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* running averages */

void create_run_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Pick";
	label2[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Running averages", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 5,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set:", xmLabelWidgetClass, rc, NULL);
	toggle_set_runavg_item = CreateSetChoice(rc, " ", maxplot, 0);

	XtVaCreateManagedWidget("Running:", xmLabelWidgetClass, rc, NULL);
	toggle_run_item = CreatePanelChoice(rc,
					    " ",
					    6,
					    "Average",
					    "Median",
					    "Minimum",
					    "Maximum",
					    "Std. dev.", 0,
					    0);
	toggle_ravglen_item = CreateTextItem4(rc, 10, "Length of average:");

	XtVaCreateManagedWidget("Restrictions:", xmLabelWidgetClass, rc, NULL);
	runavg_region_item = CreatePanelChoice(rc,
					       " ",
					       9,
					       "None",
					       "Region 0",
					       "Region 1",
					       "Region 2",
					       "Region 3",
					       "Region 4",
					       "Inside graph",
					       "Outside graph",
					       0,
					       0);

	XtVaCreateManagedWidget("Invert region:", xmLabelWidgetClass, rc, NULL);
	runavg_rinvert_item = XmCreateToggleButton(rc, " ", NULL, 0);
	XtManageChild(runavg_rinvert_item);

	XtManageChild(rc);


	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 3, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_runavg_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) do_pick_compose, (XtPointer) PICK_RUNAVG);
	XtAddCallback(but2[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* TODO finish this */
void do_eval_regress()
{
}

void create_reg_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i, x, y;
    static Widget top, dialog;
    Widget wbut, rc;
    Widget buts[4];

    set_wait_cursor();
    if (top == NULL) {
	char *label1[4];
	label1[0] = "Accept";
	label1[1] = "Pick";
	label1[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Regression", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 5,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set:", xmLabelWidgetClass, rc, NULL);
	toggle_set_regress_item = CreateSetChoice(rc, " ", maxplot, 1);

	XtVaCreateManagedWidget("Type of fit:", xmLabelWidgetClass, rc, NULL);
	toggle_degree_item = CreatePanelChoice(rc,
					       " ",
					       16,
					       "Linear",
					       "Quadratic",
					       "Cubic",
					       "4th degree",
					       "5th degree",
					       "6th degree",
					       "7th degree",
					       "8th degree",
					       "9th degree",
					       "10th degree",
					       "1-10",
					       "Power y=A*x^B",
					       "Exponential y=A*exp(B*x)",
					       "Logarithmic y=A+B*ln(x)",
					       "Inverse y=1/(A+Bx)",
					       0,
					       0);

/*
    toggle_zero_item = XmCreateToggleButton(rc, "Force fit through X = 0", NULL, 0);
    XtManageChild(toggle_zero_item);
					      "Evaluate fit",
*/

	XtVaCreateManagedWidget("Load:", xmLabelWidgetClass, rc, NULL);
	toggle_resid_item = CreatePanelChoice(rc,
					      " ",
					      3,
					      "Fitted values",
					      "Residuals",
					      0,
					      0);

/*
	regress_start_item = CreateTextItem4(rc, 10, "Start:");
	regress_stop_item = CreateTextItem4(rc, 10, "Stop:");
	regress_step_item = CreateTextItem4(rc, 10, "Number of points:");
*/

	XtVaCreateManagedWidget("Restrictions:", xmLabelWidgetClass, rc, NULL);
	regress_region_item = CreatePanelChoice(rc,
						" ",
						9,
						"None",
						"Region 0",
						"Region 1",
						"Region 2",
						"Region 3",
						"Region 4",
						"Inside graph",
						"Outside graph",
						0,
						0);
	XtVaCreateManagedWidget("Invert region:", xmLabelWidgetClass, rc, NULL);
	regress_rinvert_item = XmCreateToggleButton(rc, " ", NULL, 0);
	XtManageChild(regress_rinvert_item);
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 3, buts, label1);
	XtAddCallback(buts[0], XmNactivateCallback, (XtCallbackProc) do_regress_proc, (XtPointer) top);
	XtAddCallback(buts[1], XmNactivateCallback, (XtCallbackProc) do_pick_compose, (XtPointer) PICK_REG);
	XtAddCallback(buts[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* finite differencing */

void create_diff_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Differences", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_differ_item = CreateSetChoice(dialog, "Apply to set:", maxplot, 0);

	toggle_differ_type_item = CreatePanelChoice(dialog,
						    "Method:",
						    4,
						    "Forward difference",
						    "Backward difference",
						    "Centered difference",
						    0,
						    0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_differ_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* numerical integration */

void create_int_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Integration", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_int_item = CreateSetChoice(dialog, "Aplly to set:", maxplot, 0);

	toggle_int_type_item = CreatePanelChoice(dialog,
						 "Load:",
						 3,
						 "Cumulative sum",
						 "Sum only",
						 0,
						 0);
	int_sum_item = CreateTextItem2(dialog, 10, "Sum:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_int_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* seasonal differencing */

void create_seasonal_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Seasonal differences", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_seasonal_item = CreateSetChoice(dialog, "Seasonal difference set:", maxplot, 0);
	seasonal_period_item = CreateTextItem2(dialog, 10, "Period:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_seasonal_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* interpolation */

void create_interp_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Pick";
	label2[2] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Interpolation", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	interp_from_item = CreateSetChoice(dialog, "Interpolate from set:", maxplot, 0);
	interp_to_item = CreateSetChoice(dialog, "To set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 3, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_interp_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) do_pick_compose, (XtPointer) PICK_INTERP);
	XtAddCallback(but2[2], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* cross correlation */

void create_xcor_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "X-correlation", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_xcor1_item = CreateSetChoice(dialog, "Select set:", maxplot, 0);
	toggle_set_xcor2_item = CreateSetChoice(dialog, "Select set:", maxplot, 0);

	toggle_xcor_type_item = CreatePanelChoice(dialog,
						  "Load:",
						  3,
						  "Biased estimate",
						  "Unbiased estimate",
						  0,
						  0);
	xcor_lag_item = CreateTextItem2(dialog, 10, "Lag:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_xcor_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* splines */

void create_spline_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Splines", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 4,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set:", xmLabelWidgetClass, rc, NULL);
	spline_item = CreateSetChoice(rc, " ", maxplot, 0);

	spline_start_item = CreateTextItem4(rc, 10, "Start:");
	spline_stop_item = CreateTextItem4(rc, 10, "Stop:");
	spline_step_item = CreateTextItem4(rc, 10, "Number of points:");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_spline_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
	XtManageChild(top);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* sample a set */

void create_samp_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, rc;

    set_wait_cursor();
    if (top == NULL) {
	char *label2[3];
	label2[0] = "Accept";
	label2[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	top = XmCreateDialogShell(app_shell, "Sample points", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 5,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set:", xmLabelWidgetClass, rc, NULL);
	toggle_set_sample_item = CreateSetChoice(rc, " ", maxplot, 0);

	XtVaCreateManagedWidget("Sample type:", xmLabelWidgetClass, rc, NULL);
	sample_type_item = CreatePanelChoice(rc,
					     " ",
					     3,
					     "Start/step",
					     "Expression",
					     0,
					     0);
	sample_start_item = CreateTextItem4(rc, 10, "Start:");
	sample_step_item = CreateTextItem4(rc, 10, "Step:");
	sample_expr_item = CreateTextItem4(rc, 10, "Logical expression:");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but2, label2);
	XtAddCallback(but2[0], XmNactivateCallback, (XtCallbackProc) do_sample_proc, (XtPointer) top);
	XtAddCallback(but2[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* apply a digital filter in set 2 to set 1 */

void create_digf_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "Digital filter", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_digf1_item = CreateSetChoice(dialog, "Filter set:", maxplot, 0);
	toggle_set_digf2_item = CreateSetChoice(dialog, "With weights from set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_digfilter_proc, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* linear convolution */

void create_lconv_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "Linear convolution", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	toggle_set_linc1_item = CreateSetChoice(dialog, "Convolve set:", maxplot, 0);
	toggle_set_linc2_item = CreateSetChoice(dialog, "With set:", maxplot, 0);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_linearc_proc, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/* evaluate a formula - load the next set */

void create_leval_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "Load & evaluate", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 6,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	compute_formulax_item = CreateTextItem4(rc, 10, "X = ");
	compute_formulay_item = CreateTextItem4(rc, 10, "Y = ");

	XtVaCreateManagedWidget("Load:", xmLabelWidgetClass, rc, NULL);
	load_to2_item = CreatePanelChoice(rc,
					  " ",
					  7,
					  "Set X",
					  "Set Y",
					  "Scratch A",
					  "Scratch B",
					  "Scratch C",
					  "Scratch D", 0,
					  0);
	load_start2_item = CreateTextItem4(rc, 10, "Start load at:");
	load_stop2_item = CreateTextItem4(rc, 10, "Stop load at:");
	load_npts_item = CreateTextItem4(rc, 10, "# of points:");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_compute_proc2, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/*
 * Compute n-tiles
 */

static Widget *ntiles_set_item;
static Widget *ntiles_nt_item;
static Widget ntiles_ntval_item;
static void do_ntiles_proc(Widget w, XtPointer client_data, XtPointer call_data);

void create_ntiles_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "N-tiles", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	ntiles_set_item = CreateSetChoice(dialog, "Apply to set:", maxplot, 0);

	ntiles_nt_item = CreatePanelChoice(dialog,
					   "Compute:",
					   5,
					   "Quartiles",
					   "Deciles",
					   "Percentiles",
					   "N-tiles:",
					   0,
					   0);

	ntiles_ntval_item = CreateTextItem2(dialog, 10, "N tiles:");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_ntiles_proc, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
    }
    XtRaise(top);
    unset_wait_cursor();
}

/*
 * Rotate, scale, translate
 */

static Widget *geom_set_item;
static Widget *geom_order_item;
static Widget geom_degrees_item;
static Widget geom_rotx_item;
static Widget geom_roty_item;
static Widget geom_scalex_item;
static Widget geom_scaley_item;
static Widget geom_transx_item;
static Widget geom_transy_item;

static void do_geom_proc(Widget w, XtPointer client_data, XtPointer call_data);

void create_geom_frame(Widget w, XtPointer client_data, XtPointer call_data)
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
	top = XmCreateDialogShell(app_shell, "Geometric transformations", NULL, 0);
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	rc = XtVaCreateWidget("rc", xmRowColumnWidgetClass, dialog,
			      XmNpacking, XmPACK_COLUMN,
			      XmNnumColumns, 9,
			      XmNorientation, XmHORIZONTAL,
			      XmNisAligned, True,
			      XmNadjustLast, False,
			      XmNentryAlignment, XmALIGNMENT_END,
			      NULL);

	XtVaCreateManagedWidget("Apply to set: ", xmLabelWidgetClass, rc, NULL);
	geom_set_item = CreateSetChoice(rc, " ", maxplot, 1);
	geom_order_item = CreatePanelChoice(dialog,
					    "Apply in order:",
					    7,
					    "Rotate, translate, scale",
					    "Rotate, scale, translate",
					    "Translate, scale, rotate",
					    "Translate, rotate, scale",
					    "Scale, translate, rotate",
					    "Scale, rotate, translate",
					    0,
					    0);

	geom_degrees_item = CreateTextItem4(rc, 10, "Rotation (degrees):");
	geom_rotx_item = CreateTextItem4(rc, 10, "Rotate about X = :");
	geom_roty_item = CreateTextItem4(rc, 10, "Rotate about Y = :");
	geom_scalex_item = CreateTextItem4(rc, 10, "Scale X:");
	geom_scaley_item = CreateTextItem4(rc, 10, "Scale Y:");
	geom_transx_item = CreateTextItem4(rc, 10, "Translate X:");
	geom_transy_item = CreateTextItem4(rc, 10, "Translate Y:");
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 2, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_geom_proc, (XtPointer) top);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);

	XtManageChild(dialog);
	xv_setstr(geom_degrees_item, "0.0");
	xv_setstr(geom_rotx_item, "0.0");
	xv_setstr(geom_roty_item, "0.0");
	xv_setstr(geom_scalex_item, "1.0");
	xv_setstr(geom_scaley_item, "1.0");
	xv_setstr(geom_transx_item, "0.0");
	xv_setstr(geom_transy_item, "0.0");
    }
    XtRaise(top);
    unset_wait_cursor();
}

/*
 * compute geom
 */
static void do_geom_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i, j, k, order[3], setno, nt, loadto, graphto, ord;
    int minset, maxset;
    double degrees, sx, sy, rotx, roty, tx, ty, xtmp, ytmp, *x, *y;
    char buf[256];
    setno = (int) GetChoice(geom_set_item);
    if (setno < maxplot) {
	if (!isactive(cg, setno)) {
	    errwin("Set not active");
	    return;
	}
	minset = maxset = setno;
    } else if (setno == maxplot) {
	minset = 0;
	maxset = g[cg].maxplot - 1;
    }
    ord = (int) GetChoice(geom_order_item);
    switch (ord) {
    case 0:
	order[0] = 0;		/* rotate */
	order[1] = 1;		/* translate */
	order[2] = 2;		/* scale */
	break;
    case 1:
	order[0] = 0;
	order[1] = 2;
	order[2] = 1;
    case 2:
	order[0] = 1;
	order[1] = 2;
	order[2] = 0;
	break;
    case 3:
	order[0] = 1;
	order[1] = 0;
	order[2] = 2;
	break;
    case 4:
	order[0] = 2;
	order[1] = 1;
	order[2] = 0;
	break;
    case 5:
	order[0] = 2;
	order[1] = 0;
	order[2] = 1;
	break;
    }
    set_wait_cursor();
    for (k = minset; k <= maxset; k++) {
	if (isactive(cg, k)) {
	    x = getx(cg, k);
	    y = gety(cg, k);
	    for (j = 0; j < 3; j++) {
		switch (order[j]) {
		case 0:
		    strcpy(buf, (char *) xv_getstr(geom_degrees_item));
		    degrees = atof(buf);
		    if (degrees == 0.0) {
			break;
		    }
		    degrees = M_PI / 180.0 * degrees;
		    strcpy(buf, (char *) xv_getstr(geom_rotx_item));
		    rotx = atof(buf);
		    strcpy(buf, (char *) xv_getstr(geom_roty_item));
		    roty = atof(buf);
		    for (i = 0; i < getsetlength(cg, k); i++) {
			xtmp = x[i] - rotx;
			ytmp = y[i] - roty;
			x[i] = rotx + cos(degrees) * xtmp - sin(degrees) * ytmp;
			y[i] = roty + sin(degrees) * xtmp + cos(degrees) * ytmp;
		    }
		    break;
		case 1:
		    strcpy(buf, (char *) xv_getstr(geom_transx_item));
		    tx = atof(buf);
		    strcpy(buf, (char *) xv_getstr(geom_transy_item));
		    ty = atof(buf);
		    for (i = 0; i < getsetlength(cg, k); i++) {
			x[i] -= tx;
			y[i] -= ty;
		    }
		    break;
		case 2:
		    strcpy(buf, (char *) xv_getstr(geom_scalex_item));
		    sx = atof(buf);
		    strcpy(buf, (char *) xv_getstr(geom_scaley_item));
		    sy = atof(buf);
		    for (i = 0; i < getsetlength(cg, k); i++) {
			x[i] *= sx;
			y[i] *= sy;
		    }
		    break;
		}		/* end case */
	    }			/* end for j */
	    updatesetminmax(cg, k);
	    update_set_status(cg, k);
	}			/* end if */
    }				/* end for k */
    unset_wait_cursor();
    drawgraph();
}

/*
 * compute ntiles
 */
static void do_ntiles_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, nt, loadto, graphto;
    char fstr[256];
    setno = (int) GetChoice(ntiles_set_item);
    nt = (int) GetChoice(ntiles_nt_item);
    strcpy(fstr, (char *) xv_getstr(ntiles_ntval_item));
    switch (nt) {
    case 0:
	nt = 4;
	break;
    case 1:
	nt = 10;
	break;
    case 2:
	nt = 100;
	break;
    case 3:
	nt = atoi(fstr);
	break;
    }
    set_wait_cursor();
    do_ntiles(cg, setno, nt);
    unset_wait_cursor();
}

/*
 * evaluate a formula
 */
static void do_compute_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, loadto, graphto;
    char fstr[256];
    setno = (int) GetChoice(compute_set_item);
    loadto = (int) GetChoice(compute_load_item);
    graphto = (int) GetChoice(compute_loadgraph_item) - 1;
    strcpy(fstr, (char *) xv_getstr(compute_formula_item));
    set_wait_cursor();
    do_compute(setno, loadto, graphto, fstr);
    unset_wait_cursor();
}

/*
 * load a set
 */
static void do_load_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, toval;
    char startstr[256], stepstr[256];
    setno = (int) GetChoice(load_set_item);
    toval = (int) GetChoice(load_to_item) + 1;
    strcpy(stepstr, (char *) xv_getstr(load_step_item));
    strcpy(startstr, (char *) xv_getstr(load_start_item));
    set_wait_cursor();
    do_load(setno, toval, startstr, stepstr);
    unset_wait_cursor();
}

/*
 * evaluate a formula loading the next set
 */
static void do_compute_proc2(Widget w, XtPointer client_data, XtPointer call_data)
{
    int npts, toval;
    char fstrx[256], fstry[256];
    char startstr[256], stopstr[256];
    npts = atoi((char *) xv_getstr(load_npts_item));
    strcpy(fstrx, (char *) xv_getstr(compute_formulax_item));
    strcpy(fstry, (char *) xv_getstr(compute_formulay_item));
    strcpy(startstr, (char *) xv_getstr(load_start2_item));
    strcpy(stopstr, (char *) xv_getstr(load_stop2_item));
    toval = (int) GetChoice(load_to2_item) + 1;
    set_wait_cursor();
    do_compute2(fstrx, fstry, startstr, stopstr, npts, toval);
    unset_wait_cursor();
}

/*
 * apply a digital filter
 */
static void do_digfilter_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int set1, set2, digfiltset;
    set1 = (int) GetChoice(toggle_set_digf1_item);
    set2 = (int) GetChoice(toggle_set_digf2_item);
    set_wait_cursor();
    do_digfilter(set1, set2);
    unset_wait_cursor();
}

/*
 * linear convolution
 */
static void do_linearc_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int set1, set2, linearcset, i, itmp;
    double *xtmp;
    set1 = (int) GetChoice(toggle_set_linc1_item);
    set2 = (int) GetChoice(toggle_set_linc2_item);
    set_wait_cursor();
    do_linearc(set1, set2);
    unset_wait_cursor();
}

/*
 * cross correlation
 */
static void do_xcor_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int set1, set2, itype, lag;
    set1 = (int) GetChoice(toggle_set_xcor1_item);
    set2 = (int) GetChoice(toggle_set_xcor2_item);
    itype = (int) GetChoice(toggle_xcor_type_item);
    lag = atoi((char *) xv_getstr(xcor_lag_item));
    set_wait_cursor();
    do_xcor(set1, set2, itype, lag);
    unset_wait_cursor();
}

/*
 * splines
 */
static void do_spline_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int set, n;
    double start, stop;
    set = GetChoice(spline_item);
    start = atof((char *) xv_getstr(spline_start_item));
    stop = atof((char *) xv_getstr(spline_stop_item));
    n = atoi((char *) xv_getstr(spline_step_item));
    set_wait_cursor();
    do_spline(set, start, stop, n);
    unset_wait_cursor();
}

/*
 * numerical integration
 */
static void do_int_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, itype;
    double sum, do_int(int setno, int itype);
    setno = (int) GetChoice(toggle_set_int_item);
    itype = (int) GetChoice(toggle_int_type_item);
    set_wait_cursor();
    sum = do_int(setno, itype);
    sprintf(buf, "%lf", sum);
    xv_setstr(int_sum_item, buf);
    unset_wait_cursor();
}

/*
 * finite differences
 */
static void do_differ_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, itype;
    setno = (int) GetChoice(toggle_set_differ_item);
    itype = (int) GetChoice(toggle_differ_type_item);
    set_wait_cursor();
    do_differ(setno, itype);
    unset_wait_cursor();
}

/*
 * seasonal differences
 */
static void do_seasonal_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, period;
    setno = (int) GetChoice(toggle_set_seasonal_item);
    period = atoi(xv_getstr(seasonal_period_item));
    set_wait_cursor();
    do_seasonal_diff(setno, period);
    unset_wait_cursor();
}

/*
 * interpolation
 */
static void do_interp_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno;
    set_wait_cursor();
    unset_wait_cursor();
}

/*
 * regression
 */
static void do_regress_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, ideg, iresid, i, j, start, stop, *doit = (int *) NULL;
    int rno = GetChoice(regress_region_item) - 1;
    int invr = XmToggleButtonGetState(regress_rinvert_item);

    setno = (int) GetChoice(toggle_set_regress_item);
    if (setno == maxplot) {
	start = 0;
	stop = g[cg].maxplot - 1;
	doit = (int *) calloc(g[cg].maxplot, sizeof(int));
	for (i = 0; i < g[cg].maxplot; i++) {
	    if (isactive(cg, i)) {
		doit[i] = 1;
	    } else {
		doit[i] = 0;
	    }
	}
    } else {
	start = stop = setno;
    }
    ideg = (int) GetChoice(toggle_degree_item) + 1;
    iresid = (int) GetChoice(toggle_resid_item);
    set_wait_cursor();
    for (j = start; j <= stop; j++) {
	if (start == stop || doit[j]) {
	    if (ideg == 11) {
		for (i = 1; i <= ideg - 1; i++) {
		    do_regress(j, i, iresid, rno, invr);
		}
	    } else {
		do_regress(j, ideg, iresid, rno, invr);
	    }
	}
    }
    unset_wait_cursor();
    if (doit) {
	free(doit);
    }
    drawgraph();
}

/*
 * running averages, medians, min, max, std. deviation
 */
static void do_runavg_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int runset, runlen, runtype, setno, rno, invr;
    setno = GetChoice(toggle_set_runavg_item);
    runlen = atoi((char *) xv_getstr(toggle_ravglen_item));
    runtype = GetChoice(toggle_run_item);
    rno = GetChoice(runavg_region_item) - 1;
    invr = XmToggleButtonGetState(runavg_rinvert_item);
    set_wait_cursor();
    do_runavg(setno, runlen, runtype, rno, invr);
    unset_wait_cursor();
    drawgraph();
}

/*
 * DFT
 */
static void do_fourier_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, load, loadx, invflag, type, wind;
    setno = (int) GetChoice(toggle_set_fourier_item);
    wind = (int) GetChoice(toggle_window_fourier_item);
    load = (int) GetChoice(toggle_load_fourier_item);
    loadx = (int) GetChoice(toggle_loadx_fourier_item);
    invflag = (int) GetChoice(toggle_inv_fourier_item);
    type = (int) GetChoice(toggle_type_fourier_item);
    set_wait_cursor();
    do_fourier(0, setno, load, loadx, invflag, type, wind);
    unset_wait_cursor();
    drawgraph();
}

/*
 * DFT by FFT
 */
static void do_fft_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int load, loadx, invflag;
    int setno, type, wind;
    setno = (int) GetChoice(toggle_set_fourier_item);
    wind = (int) GetChoice(toggle_window_fourier_item);
    load = (int) GetChoice(toggle_load_fourier_item);
    loadx = (int) GetChoice(toggle_loadx_fourier_item);
    invflag = (int) GetChoice(toggle_inv_fourier_item);
    type = (int) GetChoice(toggle_type_fourier_item);
    set_wait_cursor();
    do_fourier(1, setno, load, loadx, invflag, type, wind);
    unset_wait_cursor();
    drawgraph();
}

/*
 * Apply data window only
 */
static void do_window_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, type, wind;
    setno = (int) GetChoice(toggle_set_fourier_item);
    wind = (int) GetChoice(toggle_window_fourier_item);
    type = (int) GetChoice(toggle_type_fourier_item);
    set_wait_cursor();
    do_window(setno, type, wind);
    unset_wait_cursor();
    drawgraph();
}

/*
 * histograms
 */
static void do_histo_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int fromset, toset, tograph, hist_type;
    double binw, xmin, xmax;
    fromset = (int) GetChoice(histo_set_item);
    toset = (int) GetChoice(histo_to_item) - 1;
    tograph = (int) GetChoice(histo_graph_item) - 1;
    binw = atof((char *) xv_getstr(histo_binw_item));
    xmin = atof((char *) xv_getstr(histo_hxmin_item));
    xmax = atof((char *) xv_getstr(histo_hxmax_item));
    hist_type = (int) GetChoice(histo_type_item);
    set_wait_cursor();
    do_histo(fromset, toset, tograph, binw, xmin, xmax, hist_type);
    unset_wait_cursor();
}

/*
 * sample a set, by start/step or logical expression
 */
static void do_sample_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, typeno;
    char exprstr[256];
    int startno, stepno;
    setno = (int) GetChoice(toggle_set_sample_item);
    typeno = (int) GetChoice(sample_type_item);
    strcpy(exprstr, (char *) xv_getstr(sample_expr_item));
    startno = atoi((char *) xv_getstr(sample_start_item));
    stepno = atoi((char *) xv_getstr(sample_step_item));
    set_wait_cursor();
    do_sample(setno, typeno, exprstr, startno, stepno);
    unset_wait_cursor();
}

void execute_pick_compute(int gno, int setno, int function)
{
    int loadto, graphto;
    char fstr[256];
    switch (function) {
    case PICK_EXPR:
	SetChoice(compute_set_item, setno);
	do_compute_proc((Widget) NULL, (XtPointer) 0, (XtPointer) 0);
	break;
    case PICK_RUNAVG:
	SetChoice(toggle_set_runavg_item, setno);
	do_runavg_proc((Widget) NULL, (XtPointer) 0, (XtPointer) 0);
	break;
    case PICK_REG:
	SetChoice(toggle_set_regress_item, setno);
	do_regress_proc((Widget) NULL, (XtPointer) 0, (XtPointer) 0);
	break;
    case PICK_SEASONAL:
	SetChoice(toggle_set_seasonal_item, setno);
	do_seasonal_proc((Widget) NULL, (XtPointer) 0, (XtPointer) 0);
	break;
    }
}
