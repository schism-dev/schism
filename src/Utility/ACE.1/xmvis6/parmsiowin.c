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
 * read/write data/parameter files
 *
 */

#ifndef lint
static char RCSid[] = "$Id: parmsiowin.c,v 1.3 2004/02/26 19:32:48 pturner Exp $";
#endif

#include <stdio.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
extern Widget app_shell;

Widget wparam_frame = (Widget) 0;
Widget wparam_panel;
Widget wparam_text_item;
Widget *wparam_choice_item;
static void wparam_Done_notify_proc(void);
static void wparam_apply_notify_proc(void);

Widget rparams_dialog = (Widget) NULL;	/* read params popup */

void close_rparams_popup(void)
{
    XtUnmanageChild(rparams_dialog);
}

void do_rparams_proc(void)
{
    Arg args;
    XmString list_item;
    char *s;
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rparams_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    getparms(cg, s);
    XtFree(s);
}

void create_rparams_popup(void)
{
    if (rparams_dialog) {
	XtRaise(rparams_dialog);
	return;
    }
    rparams_dialog = XmCreateFileSelectionDialog(app_shell, "Read Parameters", NULL, 0);
    XtVaSetValues(rparams_dialog, 
       XmNdialogTitle, XmStringCreate("Read Parameters", charset), 
       XmNtitle, "Read Parameters", 
       XmNcancelLabelString, XmStringCreate("Done", charset), 
       XmNdirMask, XmStringCreate("*.par", charset), 
       NULL);

    XtAddCallback(rparams_dialog, XmNcancelCallback, (XtCallbackProc) close_rparams_popup, 0);
    XtAddCallback(rparams_dialog, XmNokCallback, (XtCallbackProc) do_rparams_proc, 0);
    XtRaise(rparams_dialog);
}

/*
 * Create the wparam Frame and the wparam Panel
 */
void create_wparam_frame(void)
{
    int x, y;
    Widget wbut, rc;
    int n = 0;
    if (wparam_frame) {
	XtRaise(wparam_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    wparam_frame = XmCreateDialogShell(app_shell, "Save Parameters", NULL, 0);
    XtVaSetValues(wparam_frame, XmNx, x, XmNy, y, NULL);
    wparam_panel = XmCreateRowColumn(wparam_frame, "wparam_rc", NULL, 0);

    wparam_text_item = CreateTextItem2(wparam_panel, 20, "Save parameters to: ");

    wparam_choice_item = CreatePanelChoice1(wparam_panel, "From graph: ", 13, "All active", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "Current", NULL, NULL);

    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, wparam_panel, NULL);
    rc = XmCreateRowColumn(wparam_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wparam_apply_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wparam_Done_notify_proc, 0);
    XtManageChild(rc);

    XtManageChild(wparam_panel);
    XtRaise(wparam_frame);
}

static void wparam_Done_notify_proc(void)
{
    XtUnmanageChild(wparam_frame);
}

static void wparam_apply_notify_proc(void)
{
    int i;
    int wparamno = (int) GetChoice(wparam_choice_item);
    char s[256];
    wparamno--;
    if (wparamno == -1) {
	wparamno = -1;
    } else if (wparamno == maxgraph) {
	wparamno = cg;
    }
    strcpy(s, (char *) xv_getstr(wparam_text_item));
    if (strlen(s) > 0) {
	if (!fexists(s)) {
	    putparms(wparamno, s);
	}
    } else {
	errwin("No output file defined, enter or select filename");
    }
}
