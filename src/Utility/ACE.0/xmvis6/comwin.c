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
 * Command Panel
 *
 */

#ifndef lint
static char RCSid[] = "$Id: comwin.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

/* all declared in pars.yacc */
extern int gotbatch;
extern int gotparams;
extern int gotread;
extern int readsrc, readtype;
extern char batchfile[];
extern char paramfile[];
extern char readfile[];
extern double result;

/*
 * Widget item declarations
 */

static Widget command;
static Widget comshell;
void create_rhist_popup(void);
void create_whist_frame(void);
extern XmStringCharSet charset;

void comcall(Widget w, XtPointer cd, XmCommandCallbackStruct * s)
{
    static int errpos;
    static char val[256];
    double x, y, a, b, c, d;
    int i = 0, setno = 0, lcnt = 1;
    extern int gotbatch;
    extern char batfile[];
    char *ts;
    extern double result;
    errpos = 0;

    XmStringGetLtoR(s->value, charset, &ts);
    strcpy(val, ts);
    XtFree(ts);

    lowtoupper(val);
    strcat(val, "\n");
    scanner(val, &x, &y, 1, &a, &b, &c, &d, 1, i, setno, &errpos);
    if (gotbatch && batfile[0]) {
	runbatch(batfile);
	gotbatch = 0;
    }
}

void clear_history(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i, n;
    int ac = 0, hc;
    Arg al[5];
    Widget h = XmCommandGetChild(command, XmDIALOG_HISTORY_LIST);
    ac = 0;
    XtSetArg(al[ac], XmNhistoryItemCount, &hc);
    ac++;
    XtGetValues(command, al, ac);
    for (i = 0; i < hc; i++) {
	XmListDeletePos(h, 0);
    }
}

void replay_history(Widget w, XtPointer client_data, XtPointer call_data)
{
    static int errpos, errcount;
    static char val[128];
    extern int gotbatch;
    extern char batfile[];
    extern double result;
    char buf[256], s[256], *ts;
    int n;
    int ac = 0, hc;
    double x, y, a, b, c, d;
    int i, setno = 0, lcnt = 1;
    XmStringTable xmstrs;
    Arg al[5];
    Widget h = XmCommandGetChild(command, XmDIALOG_HISTORY_LIST);
    ac = 0;
    XtSetArg(al[ac], XmNhistoryItems, &xmstrs);
    ac++;
    XtSetArg(al[ac], XmNhistoryItemCount, &hc);
    ac++;
    XtGetValues(command, al, ac);
    if (hc > 100) {
	hc = 100;
    }
    errcount = 0;
    for (i = 0; i < hc; i++) {
	errpos = 0;
	XmStringGetLtoR(xmstrs[i], charset, &ts);
	strcpy(buf, ts);
	XtFree(ts);
	lowtoupper(buf);
	strcat(buf, "\n");
	scanner(buf, &x, &y, 1, &a, &b, &c, &d, 1, i, setno, &errpos);
	if (gotbatch && batfile[0]) {
	    runbatch(batfile);
	    gotbatch = 0;
	}
    }
}

void execute_command(char *s)
{
    static int errpos, errcount;
    static char val[128];
    char buf[256], *ts;
    extern int gotbatch;
    extern char batfile[];
    extern double result;
    int n;
    int ac = 0, hc;
    double x, y, a, b, c, d;
    int i = 0, setno = 0, lcnt = 1;
    strcpy(buf, s);
    errpos = 0;
    lowtoupper(buf);
    strcat(buf, "\n");
    scanner(buf, &x, &y, 1, &a, &b, &c, &d, 1, i, setno, &errpos);
    if (gotbatch && batfile[0]) {
	runbatch(batfile);
	gotbatch = 0;
    }
}

void open_command(void)
{
    Widget bt, fr, rc;
    extern Widget app_shell, menu_pane;
    Arg al[10];
    int ac;
    if (command) {
	XtRaise(comshell);
	return;
    }
    comshell = XmCreateDialogShell(app_shell, "Commands", NULL, 0);
    ac = 0;
    XtSetArg(al[ac], XmNpromptString, XmStringCreateLtoR("Command", (XmStringCharSet) XmSTRING_DEFAULT_CHARSET));
    ac++;
    XtSetArg(al[ac], XmNwidth, 400);
    ac++;
    XtSetArg(al[ac], XmNheight, 250);
    ac++;
    command = XmCreateCommand(comshell, "command", al, ac);

    fr = XmCreateFrame(command, "commandframe", NULL, 0);
    ac = 0;
    XtSetArg(al[ac], XmNorientation, XmHORIZONTAL);
    ac++;
    rc = XmCreateRowColumn(fr, "commandrc", al, ac);
    XtManageChild(fr);
    XtManageChild(rc);

    bt = XtVaCreateManagedWidget("Save...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_whist_frame, 0);

    bt = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_rhist_popup, 0);

    bt = XtVaCreateManagedWidget("Clear", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) clear_history, 0);

    bt = XtVaCreateManagedWidget("Replay", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) replay_history, 0);

    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) comshell);

    XtAddCallback(command, XmNcommandEnteredCallback, (XtCallbackProc) comcall, NULL);
    XtManageChild(command);
    XtRaise(comshell);
}

Widget rhist_dialog = (Widget) NULL;	/* read params popup */

void do_rhist_proc(void)
{
    Arg args, al[2];
    XmString list_item;
    char *s, buf[256];
    FILE *fp;
    int sl;
    Widget h = XmCommandGetChild(command, XmDIALOG_HISTORY_LIST);
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rhist_dialog, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    strcpy(buf, s);
    XtFree(s);
    if ((fp = fopen(buf, "r")) != NULL) {
	while (fgets(buf, 255, fp) != NULL) {
	    sl = strlen(buf);
	    buf[sl - 1] = 0;
	    list_item = XmStringCreateLtoR(buf, charset);
	    XmListAddItemUnselected(h, list_item, 0);
	    XtFree((XtPointer) list_item);
	}
	fclose(fp);
    } else {
	errwin("Unable to open file");
    }
}

void create_rhist_popup(void)
{
    Arg al[10];
    int ac;
    extern Widget app_shell;

    if (rhist_dialog) {
	XtRaise(rhist_dialog);
	return;
    }
    rhist_dialog = XmCreateFileSelectionDialog(app_shell, "Read history", NULL, 0);
    XtAddCallback(rhist_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) rhist_dialog);
    XtAddCallback(rhist_dialog, XmNokCallback, (XtCallbackProc) do_rhist_proc, 0);
    XtRaise(rhist_dialog);
}

Widget whist_frame = (Widget) 0;
Widget whist_panel;

/*
 * Panel item declarations
 */
static Widget whist_text_item;
static Widget *whist_choice_item;

/*
 * Event and Notify proc declarations
 */
static whist_Done_notify_proc(void);
static whist_apply_notify_proc(void);

/*
 * Create the whist Frame and the whist Panel
 */
void create_whist_frame(void)
{
    extern Widget app_shell;
    Widget wbut;
    Widget wlabel, rc;
    Arg wargs[50];
    XmString string;
    int n = 0, i;
    int itmp;

    if (whist_frame) {
	XtRaise(whist_frame);
	return;
    }
    whist_frame = XmCreateDialogShell(app_shell, "Write history", NULL, 0);
    whist_panel = XtVaCreateManagedWidget("history panel", xmRowColumnWidgetClass, whist_frame, NULL);

    whist_text_item = CreateTextItem2(whist_panel, 30, "Write history to:");

    rc = XmCreateRowColumn(whist_panel, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

    wbut = XtVaCreateManagedWidget("Cancel", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) whist_frame);

    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) whist_apply_notify_proc, 0);
    XtManageChild(rc);

    XtManageChild(whist_panel);
    XtRaise(whist_frame);
}

static int whist_apply_notify_proc(void)
{
    int i, ac = 0, hc;
    char s[256], *ts;
    XmStringTable xmstrs;
    Arg al[5];
    Widget h = XmCommandGetChild(command, XmDIALOG_HISTORY_LIST);
    strcpy(s, (char *) xv_getstr(whist_text_item));
    if (!fexists(s)) {
	FILE *pp = fopen(s, "w");
	if (pp != NULL) {
	    ac = 0;
	    XtSetArg(al[ac], XmNhistoryItems, &xmstrs);
	    ac++;
	    XtSetArg(al[ac], XmNhistoryItemCount, &hc);
	    ac++;
	    XtGetValues(command, al, ac);
	    for (i = 0; i < hc; i++) {
		XmStringGetLtoR(xmstrs[i], charset, &ts);
		fprintf(pp, "%s\n", ts);
		XtFree(ts);
	    }
	    fclose(pp);
	} else {
	    errwin("Unable to open file");
	}
    }
}
