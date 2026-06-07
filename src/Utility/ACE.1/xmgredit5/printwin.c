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
 * Printer initialization
 */

#ifndef lint
static char RCSid[] = "$Id: printwin.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

void cancel_ok();		/* defined in statuswin.c */
void do_hardcopy();

static Widget *printto_item;	/* for printer select popup */
static Widget printstring_item;
static Widget printscaletoggle_item;
static Widget printscale_item;
static void update_printer_setup(void);
static void set_printer_proc(void);
static void do_print(void);

Widget psetup_frame;
Widget psetup_panel;
Widget *devices_item;

int ptofile = 0;		/* flag to indicate destination of hardcopy
				 * output, ptofile = 0 means print to printer
				 * non-zero print to file */
char printstr[128] = "pout.dat";/* hardcopy to this file */
int device = 0, tdevice = 0, hdevice = GR_PS_L;
char noprint[] = "No printer installed";
int hardcopyflag = FALSE;	/* TRUE if printing out a hardcopy */
char *curprint = ps_prstr;

extern double mapscale;

static char buf[256];

#define FILEP 100

/*
 * set the current print options
 */
void set_printer(int device, char *prstr)
{
    if (device == FILEP) {
	if (prstr != NULL) {
	    strcpy(printstr, prstr);
	}
	ptofile = TRUE;
    } else {
	switch (device) {
	case GR_PS_L:
	case GR_PS_P:
	    if (prstr != NULL) {
		strcpy(ps_prstr, prstr);
	    }
	    curprint = ps_prstr;
	    break;
	default:
	    sprintf(buf, "Unknown printer device %d, printer unchanged", device);
	    errwin(buf);
	    return;
	}
	hdevice = device;
	ptofile = FALSE;
    }
    if (psetup_frame) {
	update_printer_setup();
    }
}

/*
 * set the print options
 */
void do_prstr_toggle(Widget item, int value)
{
    set_printer(value + 1, NULL);
    if ((int) GetChoice(printto_item) == 0) {
	xv_setstr(printstring_item, curprint);
    }
}

void do_pr_toggle(Widget item, int value)
{
    if (value) {
	xv_setstr(printstring_item, printstr);
    } else {
	xv_setstr(printstring_item, curprint);
    }
}

void create_printer_setup(void)
{
    extern Widget app_shell;
    Widget wbut, rc;
    int i, x, y;

    if (psetup_frame) {
	update_printer_setup();
	XtRaise(psetup_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    psetup_frame = XmCreateDialogShell(app_shell, "Printer setup", NULL, 0);
    handle_close(psetup_frame);
    XtVaSetValues(psetup_frame, XmNx, x, XmNy, y, NULL);
    psetup_panel = XmCreateRowColumn(psetup_frame, "psetup_rc", NULL, 0);

    devices_item = CreatePanelChoice1(psetup_panel, "Device:",
				      3,
				      "PostScript landscape",
				      "PostScript portrait",
				      0, 0);
    for (i = 0; i < 2; i++) {
	XtAddCallback(devices_item[2 + i], XmNactivateCallback, (XtCallbackProc) do_prstr_toggle, (XtPointer) ((long)i));
    }

    printto_item = CreatePanelChoice1(psetup_panel, "Print to:",
				      3,
				      "Printer",
				      "File", 0, 0);
    for (i = 0; i < 2; i++) {
	XtAddCallback(printto_item[2 + i], XmNactivateCallback, (XtCallbackProc) do_pr_toggle, (XtPointer) ((long)i));
    }

    printstring_item = CreateTextItem2(psetup_panel, 20, "Print control string:");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, psetup_panel, NULL);
    printscaletoggle_item = XmCreateToggleButton(psetup_panel, "Print at a given scale", NULL, 0);
    XtManageChild(printscaletoggle_item);
    printscale_item = CreateTextItem2(psetup_panel, 20, "Map scale 1 to :");

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, psetup_panel, NULL);
    rc = XmCreateRowColumn(psetup_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) set_printer_proc, 0);
    wbut = XtVaCreateManagedWidget("Print", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_print, 0);
    wbut = XtVaCreateManagedWidget("Cancel", xmPushButtonGadgetClass, rc,
				   NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, psetup_frame);
    XtManageChild(rc);

    update_printer_setup();
    XtManageChild(psetup_panel);
    XtManageChild(psetup_frame);
}

static void update_printer_setup(void)
{
    SetChoice(devices_item, hdevice - 1);
    SetChoice(printto_item, ptofile);
    if (ptofile) {
	xv_setstr(printstring_item, printstr);
    } else {
	xv_setstr(printstring_item, curprint);
    }
}

static void set_printer_proc(void)
{
    char tmpstr[128];
    hdevice = (int) GetChoice(devices_item) + 1;
    ptofile = (int) GetChoice(printto_item);
    strcpy(tmpstr, (char *) xv_getstr(printstring_item));
    if (ptofile) {
	strcpy(printstr, tmpstr);
    } else {
	strcpy(curprint, tmpstr);
    }
}

/*
 * Print button
 */
static void do_print(void)
{
    int doscale = XmToggleButtonGetState(printscaletoggle_item);
    if (doscale) {
	mapscale = atof((char *) xv_getstr(printscale_item));
    } else {
	mapscale = 1.0;
    }
    set_printer_proc();
    do_hardcopy();
    XtUnmanageChild(psetup_frame);
}
