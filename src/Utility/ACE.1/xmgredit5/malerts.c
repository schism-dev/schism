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
 * alerts for Motif
 */

#ifndef lint
static char RCSid[] = "$Id: malerts.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/MessageB.h>

static Widget error_popup;
static Widget yesno_popup;

extern Widget app_shell;
extern XmStringCharSet charset;
extern int inwin;

extern XtAppContext app_con;

static int yesno_retval = 0;
static Boolean keep_grab = True;

void yesnoCB(Widget w, Boolean * keep_grab, XmAnyCallbackStruct * reason)
{
    int why = reason->reason;
    *keep_grab = False;
    XtRemoveGrab(XtParent(w));
    XtUnmanageChild(w);
    switch (why) {
    case XmCR_OK:
	yesno_retval = 1;
	/* process ok action */
	break;
    case XmCR_CANCEL:
	yesno_retval = 0;
	/* process cancel action */
	break;
    case XmCR_HELP:
	yesno_retval = 0;
	/* process help action */
	break;
    }
}

int yesno(char *msg1, char *msg2, char *s1, char *s2)
{
    Arg al[5];
    int ac;
    char buf[256];
    static XmString str, str1, str2;
    XEvent event;
    keep_grab = True;
    if (!inwin) {
	fprintf(stderr, "%s\n", msg1);
	fprintf(stderr, "%s\n", "abort? (y/n)");
	fgets(buf, 255, stdin);
	if (buf[0] == 'y') {
	    return 1;
	} else {
	    return 0;
	}
    }
    if (yesno_popup) {
	XmStringFree(str);
	XmStringFree(str1);
	XmStringFree(str2);
	str = XmStringCreateLtoR(msg1, charset);
	XtVaSetValues(yesno_popup,
		      XmNmessageString, str,
		      NULL);
	if (s1 != NULL) {
	    XtVaSetValues(yesno_popup, XmNokLabelString, str1 = XmStringCreateLtoR(s1, charset),
			  NULL);
	} else {
	    XtVaSetValues(yesno_popup, XmNokLabelString, str1 = XmStringCreateLtoR("OK", charset),
			  NULL);
	}
	if (s2 != NULL) {
	    XtVaSetValues(yesno_popup, XmNcancelLabelString, str2 = XmStringCreateLtoR(s2, charset),
			  NULL);
	} else {
	    XtVaSetValues(yesno_popup, XmNcancelLabelString, str2 = XmStringCreateLtoR("Cancel", charset),
			  NULL);
	}
    } else {
	str = XmStringCreateLtoR(msg1, charset);
	yesno_popup = XmCreateErrorDialog(app_shell, "Yes/No", NULL, 0);
	XtVaSetValues(yesno_popup,
		      XmNmessageString, str,
		      XmNdialogTitle, XmStringCreateLtoR("Warning", charset),
		      NULL);
	if (s1 != NULL) {
	    XtVaSetValues(yesno_popup, XmNokLabelString, str1 = XmStringCreateLtoR(s1, charset),
			  NULL);
	} else {
	}
	if (s2 != NULL) {
	    XtVaSetValues(yesno_popup, XmNcancelLabelString, str2 = XmStringCreateLtoR(s2, charset),
			  NULL);
	} else {
	}
	XtAddCallback(yesno_popup, XmNokCallback, (XtCallbackProc) yesnoCB, &keep_grab);
	XtAddCallback(yesno_popup, XmNcancelCallback, (XtCallbackProc) yesnoCB, &keep_grab);
	XtAddCallback(yesno_popup, XmNhelpCallback, (XtCallbackProc) yesnoCB, &keep_grab);
    }
    XtManageChild(yesno_popup);
    XtAddGrab(XtParent(yesno_popup), True, False);
    while (keep_grab || XtAppPending(app_con)) {
	XtAppNextEvent(app_con, &event);
	XtDispatchEvent(&event);
    }
    return yesno_retval;
}

void errwin(char *s)
{
    Arg al[3];
    int ac;
    static XmString str;
    if (!inwin) {
	fprintf(stderr, "%s\n", s);
	return;
    }
    if (error_popup) {
	XmStringFree(str);
	str = XmStringCreateLtoR(s, charset);
	ac = 0;
	XtSetArg(al[ac], XmNmessageString, str);
	ac++;
	XtSetValues(error_popup, al, ac);
	XtManageChild(error_popup);
	return;
    }
    ac = 0;
    XtSetArg(al[ac], XmNmessageString, str = XmStringCreateLtoR(s, charset));
    ac++;
    XtSetArg(al[ac], XmNdialogStyle, XmDIALOG_APPLICATION_MODAL);
    ac++;
    XtSetArg(al[ac], XmNdialogTitle, XmStringCreate("Error", charset));
    ac++;
    error_popup = XmCreateErrorDialog(app_shell, "error", al, ac);
    XtUnmanageChild(XmMessageBoxGetChild(error_popup, XmDIALOG_CANCEL_BUTTON));
    XtManageChild(error_popup);
    XmUpdateDisplay(error_popup);
}
