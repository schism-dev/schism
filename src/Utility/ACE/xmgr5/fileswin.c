/* $Id: fileswin.c,v 1.2 2006/07/20 16:42:28 pturner Exp $

 * read/write data/parameter files
 *
 */

#include <stdio.h>
#ifndef WIN32
#include <sys/param.h>
#endif

#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/BulletinB.h>
#include <Xm/FileSB.h>
#include <Xm/Frame.h>
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/List.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/SelectioB.h>
#include <Xm/Separator.h>
#include <Xm/ToggleB.h>

#include "globals.h"
#include "motifinc.h"
#include "noxprotos.h"

static Widget rdata_dialog;	/* read data popup */
static Widget *read_graph_item;	/* graph choice item */
static Widget *read_ftype_item;	/* set type choice item */
static Widget read_auto_item;	/* autoscale on read button */
static Widget wparam_frame;	/* write params popup */
static Widget wparam_panel;
static Widget wparam_text_item;
static Widget *wparam_choice_item;
static void set_type_proc(int data);
static void set_src_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void rdata_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void do_rparams_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void wparam_apply_notify_proc(Widget w, XtPointer client_data, XtPointer call_data);

static Widget rparams_dialog;	/* read params popup */

static void set_type_proc(int data)
{
    switch (data) {
    case 0:
	curtype = XY;
	break;
    case 1:
	curtype = NXY;
	break;
    case 2:
	curtype = IHL;
	break;
    case 3:
	curtype = BIN;
	break;
    case 4:
	curtype = XYDX;
	break;
    case 5:
	curtype = XYDY;
	break;
    case 6:
	curtype = XYDXDX;
	break;
    case 7:
	curtype = XYDYDY;
	break;
    case 8:
	curtype = XYDXDY;
	break;
    case 9:
	curtype = XYZ;
	break;
    case 10:
	curtype = XYHILO;
	break;
    case 11:
	curtype = XYRT;
	break;
    case 12:
	curtype = XYBOX;
	break;
    case 13:
	curtype = XYBOXPLOT;
	break;
    case 14:
	curtype = CTD;
	break;
    case 15:
	curtype = XYUV;
	break;
    }
}

static void set_src_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int data = (int) client_data;

    switch (data) {
    case 0:
	cursource = DISK;
	break;
    case 1:
	cursource = PIPE;
	break;
    }
}

static void rdata_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int graphno, autoflag;
    char *s;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in rdata_proc()");
	return;
    }
    graphno = GetChoice(read_graph_item) - 1;
    autoflag = XmToggleButtonGetState(read_auto_item);
    if (graphno == -1) {
	graphno = cg;
    }
    if (g[graphno].active == OFF) {
	set_graph_active(graphno);
    }
    set_type_proc(GetChoice(read_ftype_item));
    set_wait_cursor();
    if (getdata(graphno, s, cursource, curtype)) {
	if (autoscale_onread || autoflag) {
	    autoscale_proc((Widget) NULL, (XtPointer) 0, (XtPointer) NULL);
	} else {
	    drawgraph();
	}
    }
    XtFree(s);
    unset_wait_cursor();
}

void create_file_popup(Widget wid, XtPointer client_data, XtPointer call_data)
{
    int i;
    Widget lab, rc, rc2, fr, rb, w[3];

    set_wait_cursor();

    if (rdata_dialog == NULL) {
	rdata_dialog = XmCreateFileSelectionDialog(app_shell, "rdata_dialog", NULL, 0);
	XtVaSetValues(XtParent(rdata_dialog), XmNtitle, "Read sets", NULL);
	XtAddCallback(rdata_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, rdata_dialog);
	XtAddCallback(rdata_dialog, XmNokCallback, (XtCallbackProc) rdata_proc, 0);

	curtype = XY;

	rc = XmCreateRowColumn(rdata_dialog, "Read data main RC", NULL, 0);

	fr = XmCreateFrame(rc, "frame_1", NULL, 0);
	rc2 = XmCreateRowColumn(fr, "Read data main RC", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	read_ftype_item = CreatePanelChoice(rc2, "File format: ", 17,
					    "X Y",
					    "X Y1 Y2 ... ",
					    "IHL",
					    "Binary",
					    "X Y DX",
					    "X Y DY",
					    "X Y DX1 DX2",
					    "X Y DY1 DY2",
					    "X Y DX DY",
					    "X Y Z",
					    "X HI LO OPEN CLOSE",
					    "X Y RADIUS",
					    "X Y BOX",
					    "X Y BOXPLOT",
					    "CTD",
					    "XY U V",
					    NULL, NULL);

	XtManageChild(rc2);
	XtManageChild(fr);

	fr = XmCreateFrame(rc, "frame_2", NULL, 0);
	rc2 = XmCreateRowColumn(fr, "Read data main RC", NULL, 0);
	XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
	lab = XmCreateLabel(rc2, "File Source:", NULL, 0);
	rb = XmCreateRadioBox(rc2, "radio_box_2", NULL, 0);
	XtVaSetValues(rb, XmNorientation, XmHORIZONTAL, NULL);
	w[0] = XmCreateToggleButton(rb, "Disk", NULL, 0);
	w[1] = XmCreateToggleButton(rb, "Pipe", NULL, 0);
	for (i = 0; i < 2; i++) {
	    XtAddCallback(w[i], XmNvalueChangedCallback, set_src_proc, (XtPointer) i);
	}
	XtManageChild(lab);
	XtManageChild(rb);
	XtManageChildren(w, 2);
	XtManageChild(rc2);
	XtManageChild(fr);
	XmToggleButtonSetState(w[0], True, False);

	fr = XmCreateFrame(rc, "frame_3", NULL, 0);
	rc2 = XmCreateRowColumn(fr, "Read data main RC", NULL, 0);
	read_graph_item = CreateGraphChoice(rc2, "Read to graph: ", maxgraph, 1);
	read_auto_item = XmCreateToggleButton(rc2, "Autoscale on read", NULL, 0);
	XtManageChild(read_auto_item);
	XtManageChild(rc2);
	XtManageChild(fr);
	XtManageChild(rc);

	XtManageChild(rc);
    }
    XtRaise(rdata_dialog);
    unset_wait_cursor();
}

static void do_rparams_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char *s;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in do_rparams_proc()");
	return;
    }
    set_wait_cursor();
    getparms(cg, s);
    unset_wait_cursor();
    XtFree(s);
}

void create_rparams_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget tw;
    set_wait_cursor();
    if (rparams_dialog == NULL) {
	rparams_dialog = XmCreateFileSelectionDialog(app_shell, "rparams_dialog", NULL, 0);
	XtVaSetValues(XtParent(rparams_dialog), XmNtitle, "Read parameters", NULL);
	XtAddCallback(rparams_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, rparams_dialog);
	XtAddCallback(rparams_dialog, XmNokCallback, (XtCallbackProc) do_rparams_proc, 0);
	if (plfile[0]) {
	    tw = XmFileSelectionBoxGetChild(rparams_dialog, XmDIALOG_TEXT);
	    xv_setstr(tw, plfile);
	}
    }
    XtRaise(rparams_dialog);
    unset_wait_cursor();
}

/*
 * Create the wparam Frame and the wparam Panel
 */
void create_wparam_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget wbut, fr, tw;

    set_wait_cursor();
    if (wparam_frame == NULL) {
	wparam_frame = XmCreateFileSelectionDialog(app_shell, "wparam_frame", NULL, 0);
	XtVaSetValues(XtParent(wparam_frame), XmNtitle, "Write plot parameters", NULL);
	XtAddCallback(wparam_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, wparam_frame);
	XtAddCallback(wparam_frame, XmNokCallback, (XtCallbackProc) wparam_apply_notify_proc, 0);

/* may not be needed
   handle_close(wparam_frame);
 */

	fr = XmCreateFrame(wparam_frame, "fr", NULL, 0);
	wparam_panel = XmCreateRowColumn(fr, "wparam_rc", NULL, 0);
	wparam_choice_item = CreateGraphChoice(wparam_panel, "Write parameters from graph: ", maxgraph, 2);

	if (plfile[0]) {
	    tw = XmFileSelectionBoxGetChild(wparam_frame, XmDIALOG_TEXT);
	    xv_setstr(tw, plfile);
	}
	XtManageChild(fr);
	XtManageChild(wparam_panel);
    }
    XtRaise(wparam_frame);
    unset_wait_cursor();
}

static void wparam_apply_notify_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i;
    char fname[512], *s;
    int wparamno = (int) GetChoice(wparam_choice_item);
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in wparam_apply_notify_proc()");
	return;
    }
    wparamno--;

    strcpy(fname, s);

    if (!fexists(fname)) {
	FILE *pp = fopen(fname, "w");

	if (pp != NULL) {
	    set_wait_cursor();
	    if (wparamno == -1) {
		wparamno = cg;
		putparms(wparamno, pp, 0);
		fclose(pp);
	    } else if (wparamno == maxgraph) {
		putparms(-1, pp, 0);
		fclose(pp);
	    } else {
		putparms(wparamno, pp, 0);
		fclose(pp);
	    }
	    unset_wait_cursor();
	} else {
	    errwin("Unable to open file");
	}
    }
    XtFree(s);
}

static Widget workingd_dialog;

static Widget dir_item;

static void workingdir_apply_notify_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char buf[MAXPATHLEN];
    char *s;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in workingdir_apply_notify_proc()");
	return;
    }
    strcpy(buf, s);
    XtFree(s);

    if (buf[0] == '~') {
	expand_tilde(buf);
    }
    if (chdir(buf) >= 0) {
	strcpy(workingdir, buf);
	set_title(workingdir);
	XmFileSelectionDoSearch(workingd_dialog, NULL);
    } else {
	errwin("Can't change to directory");
    }
    XtUnmanageChild(workingd_dialog);
}

static void select_dir(Widget w, XtPointer cd, XmListCallbackStruct * cbs)
{
    char buf[MAXPATHLEN], *str;

    XmStringGetLtoR(cbs->item, charset, &str);
    strcpy(buf, str);
    XtFree(str);

    xv_setstr(dir_item, buf);
    XmFileSelectionDoSearch(workingd_dialog, NULL);
}

void create_workingdir_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget filelist;
    XmString str;

    set_wait_cursor();
    if (workingd_dialog == NULL) {
	workingd_dialog = XmCreateFileSelectionDialog(app_shell, "workingd_dialog", NULL, 0);
	XtVaSetValues(XtParent(workingd_dialog), XmNtitle, "Set working directory", NULL);
	XtAddCallback(workingd_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) workingd_dialog);
	XtAddCallback(workingd_dialog, XmNokCallback, (XtCallbackProc) workingdir_apply_notify_proc, (XtPointer) 0);

/* unmanage unneeded items */
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_LIST);
	XtUnmanageChild(XtParent(w));
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_LIST_LABEL);
	XtUnmanageChild(w);
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_FILTER_LABEL);
	XtUnmanageChild(w);
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_FILTER_TEXT);
	XtUnmanageChild(w);
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_APPLY_BUTTON);
	XtUnmanageChild(w);

/* save the name of the text item used for definition */
	dir_item = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_TEXT);

/* Add a callback to the dir list */
	w = XmFileSelectionBoxGetChild(workingd_dialog, XmDIALOG_DIR_LIST);
	XtAddCallback(w, XmNsingleSelectionCallback, (XtCallbackProc) select_dir, (XtPointer) 0);
	XtVaSetValues(w, XmNselectionPolicy, XmSINGLE_SELECT, NULL);
    }
    xv_setstr(dir_item, workingdir);
    XtVaSetValues(workingd_dialog, XmNdirectory,
		  str = XmStringCreateLtoR(workingdir, charset), NULL);
    XmFileSelectionDoSearch(workingd_dialog, NULL);
    XmStringFree(str);
    XtRaise(workingd_dialog);
    unset_wait_cursor();
}

#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)

#include "netcdf.h"

/*

 * netcdf reader
 *
 */

extern int readcdf;		/* declared in main.c */

extern char netcdf_name[], xvar_name[], yvar_name[];

static Widget netcdf_frame = (Widget) NULL;

static Widget *netcdf_graph_item;
static Widget *netcdf_set_item;
static Widget netcdf_listx_item;
static Widget netcdf_listy_item;
static Widget netcdf_file_item;
static Widget netcdf_start_item;
static Widget netcdf_stop_item;
static Widget netcdf_stride_item;
static Widget netcdf_2dindex_item;
static Widget netcdf_basex_item;	/* base for X if index */
static Widget netcdf_incrx_item;	/* increment for X if index */
static Widget netcdf_auto_item;

void create_netcdffiles_popup(Widget w, XtPointer client_data, XtPointer call_data);

static void do_netcdfquery_proc(Widget w, XtPointer client_data, XtPointer call_data);

void update_netcdfs(void);

int getnetcdfvars(void);

static void do_netcdf_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int gno, setno, src;
    char fname[256];
    char buf[256], xvar[256], yvar[256];
    XmString xms;
    XmString *s, cs;
    int *pos_list;
    int j, pos_cnt, cnt, autoflag, retval;
    int stride, index2d, start, stop;
    char *cstr;

    set_wait_cursor();
    autoflag = XmToggleButtonGetState(netcdf_auto_item);

/*
 * setno == -1, then next set
 */
    gno = GetChoice(netcdf_graph_item) - 1;
    if (gno < 0) {
	gno = cg;
    }
    setno = GetChoice(netcdf_set_item) - 1;
    strcpy(fname, xv_getstr(netcdf_file_item));
    start = atoi(xv_getstr(netcdf_start_item));
    stop = atoi(xv_getstr(netcdf_stop_item));
    stride = atoi(xv_getstr(netcdf_stride_item));
    index2d = atoi(xv_getstr(netcdf_2dindex_item));
/*
   printf("start, stop, stride, index2d = %d %d %d %d\n", start, stop, stride, index2d);
 */
    if (XmListGetSelectedPos(netcdf_listx_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(netcdf_listx_item,
		      XmNselectedItemCount, &cnt,
		      XmNselectedItems, &s,
		      NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(xvar, cstr);
	    XtFree(cstr);
	}
	XmStringFree(cs);
    } else {
	errwin("Need to select X, either variable name or INDEX");
	unset_wait_cursor();
	return;
    }
    if (XmListGetSelectedPos(netcdf_listy_item, &pos_list, &pos_cnt)) {
	j = pos_list[0];
	XtVaGetValues(netcdf_listy_item,
		      XmNselectedItemCount, &cnt,
		      XmNselectedItems, &s,
		      NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(yvar, cstr);
	    XtFree(cstr);
	}
	XmStringFree(cs);
    } else {
	errwin("Need to select Y");
	unset_wait_cursor();
	return;
    }
    if (strcmp(xvar, "INDEX") == 0) {
	retval = readnetcdf(gno, setno, fname, NULL, yvar, start, stop, stride, index2d);
    } else {
	retval = readnetcdf(gno, setno, fname, xvar, yvar, start, stop, stride, index2d);
    }
    if (retval) {
	if (autoflag) {
	    autoscale_graph(gno, -3);
	} else {
	    drawgraph();
	}
    } else {			/* error from readnetcdf() */
    }

    unset_wait_cursor();
}

void update_netcdfs(void)
{
    int i, j;
    char buf[256], fname[512];
    XmString xms;
    int cdfid, ncstatus;
    int ndims, nvars, ngatts, recdim;
    int var_id;
    size_t start[2];
    size_t count[2];
    char varname[256];
    nc_type datatype;
    int dim[100], natts;
    size_t dimlen[100];
    size_t len;

    if (netcdf_frame != NULL) {
	strcpy(fname, xv_getstr(netcdf_file_item));
	set_wait_cursor();
	XmListDeleteAllItems(netcdf_listx_item);
	XmListDeleteAllItems(netcdf_listy_item);
	xms = XmStringCreateLtoR("INDEX", charset);
	XmListAddItemUnselected(netcdf_listx_item, xms, 0);
	XmStringFree(xms);

	if (strlen(fname) < 2) {
	    unset_wait_cursor();
	    return;
	}
	if ((ncstatus = nc_open(fname, NC_NOWRITE, &cdfid)) != NC_NOERR) {
	    errwin("Can't open file.");
	    unset_wait_cursor();
	    return;
	}
	nc_inq(cdfid, &ndims, &nvars, &ngatts, &recdim);
/*
   printf("%d %d %d %d\n", ndims, nvars, ngatts, recdim);
 */
	for (i = 0; i < ndims; i++) {
	    nc_inq_dim(cdfid, i, NULL, &dimlen[i]);
	}
	for (i = 0; i < nvars; i++) {
	    nc_inq_var(cdfid, i, varname, &datatype, &ndims, dim, &natts);
	    if ((ncstatus = nc_inq_varid(cdfid, varname, &var_id)) != NC_NOERR) {
		char ebuf[256];
		sprintf(ebuf, "update_netcdfs(): No such variable %s", varname);
		errwin(ebuf);
		continue;
	    }
	    if (ndims > 2) {	/* only 2d or less variables */
		continue;
	    }
	    nc_inq_dim(cdfid, dim[0], (char *) NULL, &len);
	    sprintf(buf, "%s", varname);
	    xms = XmStringCreateLtoR(buf, charset);
	    XmListAddItemUnselected(netcdf_listx_item, xms, 0);
	    XmListAddItemUnselected(netcdf_listy_item, xms, 0);
	    XmStringFree(xms);
	}
	nc_close(cdfid);
	unset_wait_cursor();
    }
}

static void do_netcdfupdate_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i;
    char buf[256];

    set_wait_cursor();
    update_netcdfs();
    unset_wait_cursor();
}

void create_netcdfs_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top, dialog;
    Widget wbut, lab, rc, rcl, rc1, rc2, form;
    Arg args[3];

    set_wait_cursor();
    if (top == NULL) {
	char *label1[5];
	Widget but1[5];

	label1[0] = "Accept";
	label1[1] = "Files...";
	label1[2] = "Update";
	label1[3] = "Query";
	label1[4] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
#ifdef HAVE_MFHDF
	top = XmCreateDialogShell(app_shell, "netCDF/HDF", NULL, 0);
#else
#ifdef HAVE_NETCDF
	top = XmCreateDialogShell(app_shell, "netCDF", NULL, 0);
#endif

#endif
	handle_close(top);
	XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
	dialog = XmCreateRowColumn(top, "dialog_rc", NULL, 0);

	netcdf_file_item = CreateTextItem2(dialog, 30, "netCDF file:");
	lab = XmCreateLabel(dialog, "Select set X:", NULL, 0);
	XtManageChild(lab);
	netcdf_listx_item = XmCreateScrolledList(dialog, "list", NULL, 0);
	XtVaSetValues(netcdf_listx_item,
		      XmNvisibleItemCount, 5,
		      NULL);
	XtManageChild(netcdf_listx_item);

	lab = XmCreateLabel(dialog, "Select set Y:", NULL, 0);
	XtManageChild(lab);
	netcdf_listy_item = XmCreateScrolledList(dialog, "list", NULL, 0);
	XtVaSetValues(netcdf_listy_item,
		      XmNvisibleItemCount, 5,
		      NULL);
	XtManageChild(netcdf_listy_item);

	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	netcdf_graph_item = CreateGraphChoice(rc, "Read to graph: ", maxgraph, 1);
	netcdf_set_item = CreateSetChoice(rc, "set:", maxplot, 4);
	XtManageChild(rc);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

 	XtVaCreateManagedWidget("All indices run from 0 to N-1", xmLabelWidgetClass, dialog, NULL);
	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	netcdf_start_item = CreateTextItem2(rc, 9, "Start index:");
	netcdf_stop_item = CreateTextItem2(rc, 9, "Stop index:");
	XtManageChild(rc);
	rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	netcdf_2dindex_item = CreateTextItem2(rc, 9, "2d index:");
	netcdf_stride_item = CreateTextItem2(rc, 9, "*Stride:");
	XtManageChild(rc);

	netcdf_auto_item = XmCreateToggleButton(dialog, "Autoscale on read", NULL, 0);
	XtManageChild(netcdf_auto_item);

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, dialog, NULL);

	CreateCommandButtons(dialog, 5, but1, label1);
	XtAddCallback(but1[0], XmNactivateCallback, (XtCallbackProc) do_netcdf_proc,
		      (XtPointer) NULL);
	XtAddCallback(but1[1], XmNactivateCallback, (XtCallbackProc) create_netcdffiles_popup,
		      (XtPointer) NULL);
	XtAddCallback(but1[2], XmNactivateCallback, (XtCallbackProc) do_netcdfupdate_proc,
		      (XtPointer) NULL);
	XtAddCallback(but1[3], XmNactivateCallback, (XtCallbackProc) do_netcdfquery_proc,
		      (XtPointer) NULL);
	XtAddCallback(but1[4], XmNactivateCallback, (XtCallbackProc) destroy_dialog,
		      (XtPointer) top);

	XtManageChild(dialog);
	netcdf_frame = top;
	if (strlen(netcdf_name)) {
	    xv_setstr(netcdf_file_item, netcdf_name);
	}
    }
    update_netcdfs();
    XtRaise(top);
    unset_wait_cursor();
}

static void do_netcdffile_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    Widget dialog = (Widget) client_data;
    char *s;
    char fname[256];
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in do_netcdffile_proc()");
	return;
    }
    set_wait_cursor();
    xv_setstr(netcdf_file_item, s);
    XtFree(s);
    unset_wait_cursor();
    XtUnmanageChild(dialog);
    update_netcdfs();
}

void create_netcdffiles_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    static Widget top;
    Widget dialog;
    Widget wbut, rc, fr;
    Arg args[2];

    set_wait_cursor();
    if (top == NULL) {
	top = XmCreateFileSelectionDialog(app_shell, "netcdfs", NULL, 0);
	XtVaSetValues(XtParent(top), XmNtitle, "Select netCDF file", NULL);

	XtAddCallback(top, XmNokCallback, (XtCallbackProc) do_netcdffile_proc, (XtPointer) top);
	XtAddCallback(top, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) top);
    }
    XtRaise(top);
    unset_wait_cursor();
}

char *getcdf_type(nc_type datatype)
{
    switch (datatype) {
    case NC_CHAR:
	return "NC_CHAR";
    case NC_BYTE:
	return "NC_BYTE";
    case NC_SHORT:
	return "NC_SHORT";
    case NC_INT:
	return "NC_INT";
    case NC_FLOAT:
	return "NC_FLOAT";
    case NC_DOUBLE:
	return "NC_DOUBLE";
    default:
	return "UNKNOWN (can't read this)";
    }
}

/*
 * TODO, lots of declared, but unused variables here
 */
static void do_netcdfquery_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int setno, src;
    char xvar[256], yvar[256];
    char buf[256], fname[512];
    XmString xms;
    XmString *s, cs;
    int *pos_list;
    int i, j, pos_cnt, cnt;
    char *cstr;

    int cdfid, ncstatus;
    int ndims, nvars, ngatts, recdim;
    int var_id;
    size_t start[2];
    size_t count[2];
    char varname[256];
    nc_type datatype;
    int dim[100], natts;
    size_t dimlen[100];
    size_t len;

    int x_id, y_id;
    nc_type xdatatype;
    nc_type ydatatype;
    int xndims, xdim[10], xnatts;
    int yndims, ydim[10], ynatts;
    size_t nx, ny;

    size_t atlen;
    char attname[256];
    char atcharval[256];

    int atint;
    double atdouble;

    set_wait_cursor();

    strcpy(fname, xv_getstr(netcdf_file_item));

    if ((ncstatus = nc_open(fname, NC_NOWRITE, &cdfid)) != NC_NOERR) {
	errwin("Can't open file.");
	goto out2;
    }
    if (XmListGetSelectedPos(netcdf_listx_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(netcdf_listx_item,
		      XmNselectedItemCount, &cnt,
		      XmNselectedItems, &s,
		      NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(xvar, cstr);
	    XtFree(cstr);
	}
	XmStringFree(cs);
    } else {
	errwin("Need to select X, either variable name or INDEX");
	goto out1;
    }
    if (XmListGetSelectedPos(netcdf_listy_item, &pos_list, &pos_cnt)) {
	XtVaGetValues(netcdf_listy_item,
		      XmNselectedItemCount, &cnt,
		      XmNselectedItems, &s,
		      NULL);
	cs = XmStringCopy(*s);
	if (XmStringGetLtoR(cs, charset, &cstr)) {
	    strcpy(yvar, cstr);
	    XtFree(cstr);
	}
	XmStringFree(cs);
    } else {
	errwin("Need to select Y");
	goto out1;
    }
/*
    nc_inq(cdfid, &ndims, &nvars, &ngatts, &recdim);
    if (ngatts <= 0) {
	stufftext("No global attributes.\n", STUFF_START);
    } else {
	stufftext("File global attributes:\n", STUFF_START);
    }
    for (i = 0; i < ngatts; i++) {
	atcharval[0] = 0;
	nc_inq_attname(cdfid, NC_GLOBAL, i, attname);
	nc_inq_att(cdfid, NC_GLOBAL, attname, &datatype, &atlen);
	switch (datatype) {
	case NC_CHAR:
	    nc_get_att_text(cdfid, NC_GLOBAL, attname, atcharval);
	    atcharval[atlen] = 0;
	    sprintf(buf, "\t%s: %s\n", attname, atcharval);
	    break;
	case NC_INT:
	    nc_get_att_int(cdfid, NC_GLOBAL, attname, &atint);
	    sprintf(buf, "\t%s: %d\n", attname, atint);
	    break;
	case NC_DOUBLE:
	    nc_get_att_double(cdfid, NC_GLOBAL, attname, &atdouble);
	    sprintf(buf, "\t%s: %lf\n", attname, atdouble);
	    break;
	}
	stufftext(buf, STUFF_TEXT);
    }
*/
/*
   printf("%d %d %d %d\n", ndims, nvars, ngatts, recdim);
 */
    if (strcmp(xvar, "INDEX") == 0) {
	stufftext("X is the index of the Y variable\n", STUFF_START);
    } else {
	if ((ncstatus = nc_inq_varid(cdfid, xvar, &x_id)) != NC_NOERR) {
	    char ebuf[256];
	    sprintf(ebuf, "do_query(): No such variable %s for X", xvar);
	    errwin(ebuf);
	    goto out1;
	}
	nc_inq_var(cdfid, x_id, NULL, &xdatatype, &xndims, xdim, &xnatts);
	nc_inq_dim(cdfid, xdim[0], NULL, &nx);
	sprintf(buf, "X is %s, data type %s \t length [%d]\n", xvar, getcdf_type(xdatatype), nx);
	stufftext(buf, STUFF_START);
	sprintf(buf, "\t%d Attributes:\n", xnatts);
	stufftext(buf, STUFF_TEXT);
	for (i = 0; i < xnatts; i++) {
	    atcharval[0] = 0;
	    nc_inq_attname(cdfid, x_id, i, attname);
	    nc_inq_att(cdfid, x_id, attname, &datatype, &atlen);
	    switch (datatype) {
	    case NC_CHAR:
		nc_get_att_text(cdfid, x_id, attname, atcharval);
		atcharval[atlen] = 0;
		break;
	    }
	    sprintf(buf, "\t\t%s: %s\n", attname, atcharval);
	    stufftext(buf, STUFF_TEXT);
	}
    }
    if ((ncstatus = nc_inq_varid(cdfid, yvar, &y_id)) != NC_NOERR) {
	char ebuf[256];
	sprintf(ebuf, "do_query(): No such variable %s for Y", yvar);
	errwin(ebuf);
	goto out1;
    }
    nc_inq_var(cdfid, y_id, NULL, &ydatatype, &yndims, ydim, &ynatts);
    nc_inq_dim(cdfid, ydim[0], NULL, &ny);
    sprintf(buf, "Y is %s, data type %s \t length [%d]\n", yvar, getcdf_type(ydatatype), ny);
    stufftext(buf, STUFF_TEXT);
    sprintf(buf, "\t%d Attributes:\n", ynatts);
    stufftext(buf, STUFF_TEXT);
    for (i = 0; i < ynatts; i++) {
	atcharval[0] = 0;
	nc_inq_attname(cdfid, y_id, i, attname);
	nc_inq_att(cdfid, y_id, attname, &datatype, &atlen);
	switch (datatype) {
	case NC_CHAR:
	    nc_get_att_text(cdfid, y_id, attname, atcharval);
	    atcharval[atlen] = 0;
	    break;
	}
	sprintf(buf, "\t\t%s: %s\n", attname, atcharval);
	stufftext(buf, STUFF_TEXT);
    }

  out1:;
    nc_close(cdfid);

  out2:;
    stufftext("\n", STUFF_STOP);
    unset_wait_cursor();
}

#endif
