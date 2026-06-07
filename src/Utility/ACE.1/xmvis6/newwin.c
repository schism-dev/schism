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
 * Open another copy of xmvis
 */

#ifndef lint
static char RCSid[] = "$Id: newwin.c,v 1.4 2004/02/26 19:32:48 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;

Widget new_frame;
Widget new_panel;

/*
 * Event and Notify proc declarations
 */
void new_done_proc(void);
void new_accept_proc(void);

extern Widget app_shell;

void create_new_frame(void)
{
    Widget wbut, rc;
    int i;
    setistop();
    if (!new_frame) {
	new_frame = XmCreateFileSelectionDialog(app_shell, "New", NULL, 0);
	XtVaSetValues(new_frame, XmNdialogTitle, XmStringCreate("New", charset), XmNtitle, "New", XmNcancelLabelString, XmStringCreate("Done", charset), XmNdirMask, XmStringCreate("*.gr*", charset), NULL);

	XtAddCallback(new_frame, XmNcancelCallback, (XtCallbackProc) new_done_proc, NULL);
	XtAddCallback(new_frame, XmNokCallback, (XtCallbackProc) new_accept_proc, NULL);
    }
    XtRaise(new_frame);
}

void new_done_proc(void)
{
    XtUnmanageChild(new_frame);
}

void new_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(new_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    set_wait_cursor(new_frame);
    sprintf(buf, "xmvis6 %s", s);
    system(buf);
    unset_wait_cursor(new_frame);
    XtUnmanageChild(new_frame);
}

/*
 * Open a copy of xmgredit
 */

Widget gredit_frame;
Widget gredit_panel;

/*
 * Event and Notify proc declarations
 */
void gredit_done_proc(void);
void gredit_accept_proc(void);

extern Widget app_shell;

void create_gredit_frame(void)
{
    Widget wbut, rc;
    int i;
    setistop();
    if (!gredit_frame) {
	gredit_frame = XmCreateFileSelectionDialog(app_shell, "gredit_frame", NULL, 0);
	XtVaSetValues(gredit_frame, XmNtitle, "New", XmNcancelLabelString, XmStringCreate("Done", charset), XmNdirMask, XmStringCreate("*.gr*", charset), NULL);

	XtAddCallback(gredit_frame, XmNcancelCallback, (XtCallbackProc) gredit_done_proc, NULL);
	XtAddCallback(gredit_frame, XmNokCallback, (XtCallbackProc) gredit_accept_proc, NULL);
    }
    XtRaise(gredit_frame);
}

void gredit_done_proc(void)
{
    XtUnmanageChild(gredit_frame);
}

void gredit_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(gredit_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);

    set_wait_cursor(gredit_frame);
    sprintf(buf, "xmgredit %s", s);
    system(buf);
    unset_wait_cursor(gredit_frame);
    XtUnmanageChild(gredit_frame);
}
