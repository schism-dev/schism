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

 * read/write objects
 *
 */
#ifndef lint
static char RCSid[] = "$Id: fileswin.c,v 1.4 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

/* fileswin.c */
static void set_format_proc(Widget w, int data, caddr_t call_data);
static void rdata_proc(void);
void create_file_popup(void);
void create_wdata_frame(void);
static void wdata_apply_notify_proc(void);
static void backgrid_writefile_proc(char *fname, int form);
static void grid_writefile_proc(char *fname, int form);
static void grid_readfile_proc(char *fname, int form);
static void writebfile_proc(char *fname, int form, int val);
static void write_backbound_proc(char *fname, int form, int val);
static void readbfile_proc(char *fname, int form, int val);
static void writebufile_proc(char *fname, int form);
static void readbufile_proc(char *fname, int form);
static void read_backgrid_proc(char *fname, int form);
static void read_backbound_proc(char *fname, int form);
void read_region_proc(Widget w, XtPointer clientd, XtPointer calld);
void do_read_region(void);
void save_region_proc(Widget w, XtPointer clientd, XtPointer calld);
void do_save_region(void);

static int curtype = 0;
static int curformat = 0;

static Widget rdata_dialog;
static Widget *rdata_choice_item;
static Widget rdata_swap_item;

static void set_format_proc(Widget w, int data, caddr_t call_data)
{
    curformat = data;
}

static void rdata_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[1024];
    FILE *fp;
    int c;
    struct stat statb;
    swapBytes = XmToggleButtonGetState(rdata_swap_item);

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rdata_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    /* check to make sure this is a file and not a dir */
    if (!isfile(s)) {
	sprintf(buf, "File %s is not a regular file", s);
	errwin(buf);
	return;
    }
    set_wait_cursor(rdata_dialog);
    curtype = GetChoice(rdata_choice_item);
    switch (curtype) {
    case 0:
	grid_readfile_proc(s, curformat);
	break;
    case 1:
	readbfile_proc(s, curformat, 0);
	break;
    case 2:
	readbfile_proc(s, curformat, 1);
	break;
    case 3:
	read_backgrid_proc(s, curformat);
	break;
    case 4:
	read_backbound_proc(s, curformat);
	break;
    case 5:
	readbufile_proc(s, curformat);
	break;
    case 6:
	getparms_binary(0, s);
	break;
    case 7:
	mergegrid(0, s);
	break;
    }
    update_fuzz_items();
    unset_wait_cursor(rdata_dialog);
}

void create_file_popup(void)
{
    int i, nitems = 2;
    Widget lab, rc, fr, rb, w[10];
    Arg a;

    if (rdata_dialog) {
	XmToggleButtonSetState(rdata_swap_item, swapBytes, False);
	XtRaise(rdata_dialog);
	return;
    }
    rdata_dialog = XmCreateFileSelectionDialog(app_shell, "Read data", NULL, 0);
    XtSetArg(a, XmNdialogTitle, XmStringCreateLtoR("Read data", charset));
    XtSetValues(rdata_dialog, &a, 1);
    XtAddCallback(rdata_dialog, XmNcancelCallback,
	      (XtCallbackProc) destroy_dialog, (XtPointer) rdata_dialog);
    XtAddCallback(rdata_dialog, XmNokCallback, (XtCallbackProc) rdata_proc, 0);

    rc = XmCreateRowColumn(rdata_dialog, "rc", NULL, 0);
    rdata_choice_item = CreatePanelChoice1(rc, "Read object: ",
					   9,
					   "Edit grid",
					   "Edit boundary (X, Y)",
					   "Edit boundary (nodes)",
					   "Background grid",
					   "Coastalboundary",
					   "Build points",
					   "Parameters",
					   "Merge grid",
					   NULL,
					   NULL);
    lab = XmCreateLabelGadget(rc, "File format:", NULL, 0);
    fr = XmCreateFrame(rc, "frame_2", NULL, 0);
    rb = XmCreateRadioBox(fr, "radio_box_2", NULL, 0);
    w[0] = XmCreateToggleButton(rb, "ASCII", NULL, 0);
    w[1] = XmCreateToggleButton(rb, "Binary", NULL, 0);
#ifdef HAVE_NETCDF
    w[2] = XmCreateToggleButton(rb, "Netcdf", NULL, 0);
    nitems = 3;
#endif
    for (i = 0; i < nitems; i++) {
	XtAddCallback(w[i], XmNvalueChangedCallback, (XtCallbackProc) set_format_proc, (XtPointer) ((long) i));
    }
    XtManageChild(lab);
    XtManageChild(fr);
    XtManageChild(rb);
    XtManageChildren(w, nitems);
    XmToggleButtonSetState(w[0], True, False);
    rdata_swap_item = XtVaCreateManagedWidget("Swap bytes",
				    xmToggleButtonWidgetClass, rc, NULL);

    XtManageChild(rc);
    XtManageChild(rdata_dialog);
}

static Widget wdata_frame;
static Widget wdata_panel;

/*
 * Panel item declarations
 */
static Widget wdata_text_item;
static Widget *wdata_choice_item;
static Widget *wdata_format_item;

/*
 * Event and Notify proc declarations
 */
static void wdata_apply_notify_proc(void);

/*
 * Create the wdata Frame and the wdata Panel
 */
void create_wdata_frame(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    Widget CreateTextItem2();
    Widget *CreatePanelChoice1();

    if (wdata_frame) {
	XtRaise(wdata_frame);
	return;
    }
    wdata_frame = XmCreateDialogShell(app_shell, "Write data", NULL, 0);
    wdata_panel = XmCreateRowColumn(wdata_frame, "write_rc", NULL, 0);

    wdata_choice_item = CreatePanelChoice1(wdata_panel, "Write object: ",
					   8,
					   "Edit grid",
					   "Edit boundary (X, Y)",
					   "Edit boundary (nodes)",
					   "Background grid",
					   "Background boundary",
					   "Build points",
					   "Parameters",
					   NULL,
					   NULL);

#ifdef HAVE_NETCDF
    wdata_format_item = CreatePanelChoice1(wdata_panel, "Using format: ",
					   4,
					   "ASCII",
					   "Binary",
					   "Netcdf",
					   NULL,
					   NULL);
#else
    wdata_format_item = CreatePanelChoice1(wdata_panel, "Using format: ",
					   3,
					   "ASCII",
					   "Binary",
					   NULL,
					   NULL);
#endif

    wdata_text_item = CreateTextItem2(wdata_panel, 20, "To file: ");

    rc = XmCreateRowColumn(wdata_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wdata_apply_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) wdata_frame);
    XtManageChild(rc);

    XtManageChild(wdata_panel);
    XtManageChild(wdata_frame);
}

static void wdata_apply_notify_proc(void)
{
    char s[256];
    int wdataobj = (int) GetChoice(wdata_choice_item);
    int wdataformat = (int) GetChoice(wdata_format_item);

    strcpy(s, (char *) xv_getstr(wdata_text_item));
    switch (wdataobj) {
    case 0:
	grid_writefile_proc(s, wdataformat);
	break;
    case 1:
	writebfile_proc(s, wdataformat, 1);
	break;
    case 2:
	writebfile_proc(s, wdataformat, 0);
	break;
    case 3:
	backgrid_writefile_proc(s, wdataformat);
	break;
    case 4:
	write_backbound_proc(s, wdataformat, 1);
	break;
    case 5:
	writebufile_proc(s, wdataformat);
	break;
    case 6:
	if (s[0]) {
	    if (!fexists(s)) {
		putparms_binary(0, s);
	    }
	} else {
	    errwin("Define file name first");
	}
	break;
    }

}

static void backgrid_writefile_proc(char *fname, int form)
{
    if (fname[0]) {
	if (!fexists(fname)) {
	    switch (form) {
	    case 0:
		writegrid(backgrid, fname);
		break;
	    case 1:
		writegridbin(backgrid, fname);
		break;
#ifdef HAVE_NETCDF
	    case 2:
		CreateGridnetcdf(fname, backgrid);
		WriteGridnetcdf(fname, backgrid);
		break;
#endif
	    }
	}
    } else {
	errwin("Define file name first");
    }
}

static void grid_writefile_proc(char *fname, int form)
{
    if (fname[0]) {
	if (!fexists(fname)) {
	    switch (form) {
	    case 0:
		writegrid(curgrid, fname);
		break;
	    case 1:
		writegridbin(curgrid, fname);
		break;
#ifdef HAVE_NETCDF
	    case 2:
		CreateGridnetcdf(fname, curgrid);
		WriteGridnetcdf(fname, curgrid);
		break;
#endif
	    }
	}
    } else {
	errwin("Define file name first");
    }
}

static void grid_readfile_proc(char *fname, int form)
{
    if (fname[0]) {
	switch (form) {
	case 0:
	    if (readgrid(curgrid, fname)) {
	    }
	    break;
	case 1:
	    if (readgridbin(curgrid, fname)) {
	    }
	    break;
#ifdef HAVE_NETCDF
	case 2:
	    ReadGridnetcdf(fname, curgrid);
	    break;
#endif
	}
    } else {
	errwin("Define file name first");
    }
}

static void writebfile_proc(char *fname, int form, int val)
{
    if (fname[0]) {
	if (!fexists(fname)) {
	    switch (form) {
	    case 0:
		writeboundary(curgrid, fname, val);
		break;
	    case 1:
		errwin("No binary format for boundaries as yet");
		break;
#ifdef HAVE_NETCDF
	    case 2:
		CreateBoundarynetcdf(fname, curgrid);
		WriteBoundarynetcdf(fname, curgrid);
		break;
#endif
	    }
	}
    } else {
	errwin("Define file name first");
    }
}

static void write_backbound_proc(char *fname, int form, int val)
{
    if (fname[0]) {
	if (!fexists(fname)) {
	    switch (form) {
	    case 0:
		writeboundary(MAXGRIDS, fname, val);
		break;
	    case 1:
		errwin("No binary format for boundaries as yet");
		break;
#ifdef HAVE_NETCDF
	    case 2:
		CreateBoundarynetcdf(fname, MAXGRIDS);
		WriteBoundarynetcdf(fname, MAXGRIDS);
		break;
#endif
	    }
	}
    } else {
	errwin("Define file name first");
    }
}

static void readbfile_proc(char *fname, int form, int val)
{
    if (fname[0]) {
	switch (form) {
	case 0:
	    readboundary(curgrid, fname, val);
	    break;
	case 1:
	    readbinboundary(curgrid, fname, val);
	    break;
#ifdef HAVE_NETCDF
	case 2:
	    ReadBoundarynetcdf(fname, curgrid);
	    break;
#endif
	}
    } else {
	errwin("Define file name first");
    }
}

static void writebufile_proc(char *fname, int form)
{
    if (fname[0]) {
	if (!fexists(fname)) {
	    switch (form) {
	    case 0:
		writebuild(curbuild, fname);
		break;
	    case 1:
		writebuildbinary(curbuild, fname);
		break;
#ifdef HAVE_NETCDF
	    case 2:
		CreateBuildnetcdf(fname, curbuild);
		WriteBuildnetcdf(fname, curbuild);
		break;
#endif
	    }
	}
    } else {
	errwin("Define file name first");
    }
}

static void readbufile_proc(char *fname, int form)
{

    if (fname[0]) {
	switch (form) {
	case 0:
	    readbuild(curbuild, fname);
	    break;
	case 1:
	    readbuildbinary(curbuild, fname);
	    break;
#ifdef HAVE_NETCDF
	case 2:
	    ReadBuildnetcdf(fname, curbuild);
	    break;
#endif
	}
    } else {

	errwin("Define file name first");
    }
}

static void read_backgrid_proc(char *fname, int form)
{
    if (fname[0]) {
	switch (form) {
	case 0:
	    if (readgrid(backgrid, fname)) {
	    }
	    break;
	case 1:
	    if (readgridbin(backgrid, fname)) {
	    }
	    break;
#ifdef HAVE_NETCDF
	case 2:
	    ReadGridnetcdf(fname, backgrid);
	    break;
#endif
	}
    } else {
	errwin("Define file name first");
    }
}

static void read_backbound_proc(char *fname, int form)
{
    if (fname[0]) {
	switch (form) {
	case 0:
	    readboundary(backgrid, fname, 0);
	    break;
	case 1:
	    readbinboundary(backgrid, fname, 0);
	    break;
#ifdef HAVE_NETCDF
	case 2:
	    ReadBoundarynetcdf(fname, backgrid);
	    break;
#endif
	}
    } else {

	errwin("Define file name first");
    }
}

/*
 * read a region
 */
void read_region_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    char *s, buf[1024];
    int i;
    FILE *fp;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) calld;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in read_region_proc()");
	return;
    }
    if ((fp = fopen(s, "r")) != NULL) {
	fgets(buf, 1023, fp);
	fgets(buf, 1023, fp);
	fgets(buf, 1023, fp);
	sscanf(buf, "%d", &nregion);
	for (i = 0; i < nregion; i++) {
	    fgets(buf, 1023, fp);
	    sscanf(buf, "%lf %lf", &regionx[i], &regiony[i]);
	}
	fclose(fp);
	region_flag = 1;
    } else {
	errwin("Unable to open file");
	return;

    }
    XtFree(s);
}

/*
 * Read a region popup
 */
void do_read_region(void)
{
    static Widget top;
    int i;
    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateFileSelectionDialog(app_shell, "Read region", NULL, 0);
    XtVaSetValues(XtParent(top), XmNtitle, "Read region", NULL);
    XtAddCallback(top, XmNcancelCallback,
		  (XtCallbackProc) destroy_dialog, (XtPointer) top);
    XtAddCallback(top, XmNokCallback, (XtCallbackProc) read_region_proc,
		  (XtPointer) top);
    XtManageChild(top);
}

/*
 * Save a region
 */
void save_region_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    char *s, buf[1024];
    int i;
    FILE *fp;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) calld;
    if (!XmStringGetLtoR(cbs->value, charset, &s)) {
	errwin("Error converting XmString to char string in read_region_proc()");
	return;
    }
    if ((fp = fopen(s, "w")) != NULL) {
	fprintf(fp, "Region written by ACE/gredit\n");
	fprintf(fp, "1\n");
	fprintf(fp, "%d 1\n", nregion);
	for (i = 0; i < nregion; i++) {
	    fprintf(fp, "%lf %lf\n", regionx[i], regiony[i]);
	}
	fclose(fp);
    } else {
	errwin("Unable to open file");
	return;

    }
    XtFree(s);
}

void do_save_region(void)
{
    static Widget top;
    int i;
    if (region_flag == 0) {
	errwin("No region to save, operation cancelled");
	return;
    }
    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateFileSelectionDialog(app_shell, "Save region", NULL, 0);
    XtVaSetValues(XtParent(top), XmNtitle, "Save region", NULL);
    XtAddCallback(top, XmNcancelCallback,
		  (XtCallbackProc) destroy_dialog, (XtPointer) top);
    XtAddCallback(top, XmNokCallback, (XtCallbackProc) save_region_proc,
		  (XtPointer) top);
    XtManageChild(top);
}
