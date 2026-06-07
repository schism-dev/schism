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
 * monitor Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: monwin.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include "motifinc.h"

extern XmStringCharSet charset;

static XFontStruct *f;
static XmFontList xmf;

void create_wmon_frame(void);
static void clear_results(Widget w, XtPointer client_data, XtPointer call_data);

static Widget mon_frame, mon_panel, text_w;
static Widget wmon_text_item;
static Widget wmon_frame;

static void mon_Done_notify_proc(void);

/*
 * Create the mon Panel
 */
void create_monitor_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc, fr;
    set_wait_cursor();
    if (mon_frame == NULL) {
	Display *disp = XtDisplay(app_shell);
	f = (XFontStruct *) XLoadQueryFont(disp, "fixed");
	xmf = XmFontListCreate(f, charset);
	XmGetPos(app_shell, 0, &x, &y);
	mon_frame = XmCreateDialogShell(app_shell, "Results", NULL, 0);
	handle_close(mon_frame);
	XtVaSetValues(mon_frame, XmNx, x, XmNy, y, NULL);
	mon_panel = XmCreateForm(mon_frame, "mon_form", NULL, 0);
	fr = XmCreateFrame(mon_panel, "fr", NULL, 0);
	text_w = XmCreateScrolledText(fr, "text_w", NULL, 0);
	XtVaSetValues(XtParent(text_w), XmNscrollVertical, True, NULL);
	XtVaSetValues(text_w, XmNrows, 10, XmNcolumns, 80, XmNeditMode, XmMULTI_LINE_EDIT, XmNfontList, xmf, XmNwordWrap, True, NULL);
	XtManageChild(text_w);
	XtManageChild(fr);

	rc = XmCreateRowColumn(mon_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_wmon_frame, (XtPointer) NULL);
	wbut = XtVaCreateManagedWidget("Clear", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) clear_results, (XtPointer) NULL);
	wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) mon_frame);
	XtManageChild(rc);

	XtVaSetValues(fr, XmNtopAttachment, XmATTACH_FORM, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_WIDGET, XmNbottomWidget, rc, NULL);
	XtVaSetValues(rc, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_FORM, NULL);

	XtManageChild(mon_panel);
    }
    XtRaise(mon_frame);
    unset_wait_cursor();
}

/*
 * Notify and event procs
 */
static void mon_Done_notify_proc(void)
{
    XtUnmanageChild(mon_frame);
}

static void clear_results(Widget w, XtPointer client_data, XtPointer call_data)
{
    XmTextSetString(text_w, "");
}

/*

#define STUFF_TEXT 0
#define STUFF_START 1
#define STUFF_STOP 2

 * sp = 0, just put text
 * sp = 1, place text at end initialize savepos (for sp = 2)
 * sp = 2, place text at end and go to the beginning of the sequence
 */
void stufftext(char *s, int sp)
{
    extern int inwin;
    static XmTextPosition pos = 0;
    static XmTextPosition savepos = 0;

    if (inwin) {
	create_monitor_frame(NULL, NULL, NULL);
	if (sp == 1) {
	    pos = XmTextGetLastPosition(text_w);
	    savepos = pos;
	    XmTextSetTopCharacter(text_w, savepos);
	}
	XmTextInsert(text_w, pos, s);
	pos += strlen(s);
	if (sp == 2) {
	    XmTextSetTopCharacter(text_w, savepos);
	    savepos = pos;
	}
    } else {
	printf(s);
    }
}

static void wmon_apply_notify_proc(void);
static void wmon_Done_notify_proc(void);

/*
 * Create the wparam Frame and the wparam Panel
 */
void create_wmon_frame(void)
{
    int x, y;
    extern Widget app_shell;
    Widget wbut, rc, wmon_panel;
    int n = 0;
    Widget CreateTextItem(Widget parent, int x, int y, int len, char *s);

    if (wmon_frame) {
	XtRaise(wmon_frame);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    wmon_frame = XmCreateDialogShell(app_shell, "Save", NULL, 0);
    XtVaSetValues(wmon_frame, XmNx, x, XmNy, y, NULL);
    wmon_panel = XmCreateRowColumn(wmon_frame, "wmon_rc", NULL, 0);

    wmon_text_item = CreateTextItem2(wmon_panel, 20, "Save to file: ");
    XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, wmon_panel, NULL);
    rc = XmCreateRowColumn(wmon_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wmon_apply_notify_proc, 0);
    wbut = XtVaCreateManagedWidget("Cancel", xmPushButtonGadgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) wmon_Done_notify_proc, 0);
    XtManageChild(rc);

    XtManageChild(wmon_panel);
    XtRaise(wmon_frame);
}

static void wmon_Done_notify_proc(void)
{
    XtUnmanageChild(wmon_frame);
}

static void wmon_apply_notify_proc(void)
{
    int i, len;
    char s[256];
    char *text, *fname;
    strcpy(s, (char *) xv_getstr(wmon_text_item));
    if (!fexists(s)) {
	FILE *pp = fopen(s, "w");
	if (pp != NULL) {
	    text = XmTextGetString(text_w);
	    len = XmTextGetLastPosition(text_w);
	    fwrite(text, sizeof(char), len, pp);
	    fclose(pp);
	    XtFree(text);
	} else {
	    errwin("Unable to open file");
	}
    }
}
