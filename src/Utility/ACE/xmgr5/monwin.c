/* $Id: monwin.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * monitor Panel
 *
 */
#include <stdio.h>
#include <Xm/Xm.h>
#include <Xm/Text.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>

#include "globals.h"
#include "motifinc.h"

extern Display *disp;

static XFontStruct *f;
static XmFontList xmf;

static Widget mon_frame, mon_panel, text_w;
static Widget wmon_text_item;
static Widget wmon_frame;

static void clear_results(Widget w, XtPointer client_data, XtPointer call_data);
static void create_wmon_frame(Widget w, XtPointer client_data, XtPointer call_data);
static void wmon_apply_notify_proc(Widget w, XtPointer client_data, XtPointer call_data);

/*
 * Create the mon Panel
 */
void create_monitor_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc, fr;
    set_wait_cursor();
    if (mon_frame == NULL) {
	f = (XFontStruct *) XLoadQueryFont(disp, "fixed");
	xmf = XmFontListCreate(f, charset);
	XmGetPos(app_shell, 0, &x, &y);
	mon_frame = XmCreateDialogShell(app_shell, "Results", NULL, 0);
	handle_close(mon_frame);
	XtVaSetValues(mon_frame, XmNx, x, XmNy, y, NULL);
	mon_panel = XmCreateForm(mon_frame, "mon_form", NULL, 0);
	fr = XmCreateFrame(mon_panel, "fr", NULL, 0);
	text_w = XmCreateScrolledText(fr, "text_w", NULL, 0);
	XtVaSetValues(text_w,
		      XmNrows, 10,
		      XmNcolumns, 80,
		      XmNeditMode, XmMULTI_LINE_EDIT,
		      XmNfontList, xmf,
		      XmNwordWrap, True,
		      NULL);
	XtManageChild(text_w);
	XtManageChild(fr);

	rc = XmCreateRowColumn(mon_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_wmon_frame, (XtPointer) NULL);
	wbut = XtVaCreateManagedWidget("Clear", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) clear_results, (XtPointer) NULL);
	wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) mon_frame);
	XtManageChild(rc);

	XtVaSetValues(fr,
		      XmNtopAttachment, XmATTACH_FORM,
		      XmNleftAttachment, XmATTACH_FORM,
		      XmNrightAttachment, XmATTACH_FORM,
		      XmNbottomAttachment, XmATTACH_WIDGET,
		      XmNbottomWidget, rc,
		      NULL);
	XtVaSetValues(rc,
		      XmNleftAttachment, XmATTACH_FORM,
		      XmNrightAttachment, XmATTACH_FORM,
		      XmNbottomAttachment, XmATTACH_FORM,
		      NULL);

	XtManageChild(mon_panel);
    }
    XtRaise(mon_frame);
    unset_wait_cursor();
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
    static XmTextPosition pos = 0; /* total character count */
    static XmTextPosition savepos = 0; /* for this sequence */

    if (inwin) {
	create_monitor_frame(NULL, NULL, NULL);
	if (sp == STUFF_START) {
	    pos = XmTextGetLastPosition(text_w);
	    savepos = pos;
	    XmTextSetTopCharacter(text_w, savepos);
	}
	XmTextInsert(text_w, pos, s);
	pos += strlen(s);
	if (sp == STUFF_STOP) {
	    XmTextSetTopCharacter(text_w, savepos);
	    savepos = pos;
	}
    } else {
	printf(s);
    }
    if (resfp != NULL) {	/* results file opened in main.c */
	fprintf(resfp, s);
    }
}

/*
 * Create the wparam Frame and the wparam Panel
 */
static void create_wmon_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc, wmon_panel;
    Widget buts[2];

    int n = 0;
    unset_wait_cursor();
    if (wmon_frame == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	wmon_frame = XmCreateDialogShell(app_shell, "Save", NULL, 0);
	handle_close(wmon_frame);
	XtVaSetValues(wmon_frame, XmNx, x, XmNy, y, NULL);
	wmon_panel = XmCreateRowColumn(wmon_frame, "wmon_rc", NULL, 0);

	wmon_text_item = CreateTextItem2(wmon_panel, 20, "Save to file: ");
	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, wmon_panel, NULL);

	CreateCommandButtons(wmon_panel, 2, buts, label1);
	XtAddCallback(buts[0], XmNactivateCallback,
		    (XtCallbackProc) wmon_apply_notify_proc, (XtPointer) 0);
	XtAddCallback(buts[1], XmNactivateCallback,
		   (XtCallbackProc) destroy_dialog, (XtPointer) wmon_frame);

	XtManageChild(wmon_panel);
    }
    XtRaise(wmon_frame);
    unset_wait_cursor();
}

static void wmon_apply_notify_proc(Widget w, XtPointer client_data, XtPointer call_data)
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
