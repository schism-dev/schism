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
 *
 * utilities for Motif
 *
 */

#ifndef lint
static char RCSid[] = "$Id: motifutils.c,v 1.4 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/LabelG.h>
#include <Xm/PushBG.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/Text.h>
#include <Xm/Protocols.h>

#define MAXITEMS 100
#define MAXARGS 100

extern Widget app_shell;

static XFontStruct *f;

/*
int tolower(c)
int c;
{
   if (c >= 'A' && c <= 'Z') {
        return(c - ' ');
    }
    else {
        return c;
    }
}
*/

void SetChoice(Widget * w, int value)
{
    Arg a;
    XtSetArg(a, XmNmenuHistory, w[value + 2]);
    XtSetValues(w[0], &a, 1);
}

int GetChoice(Widget * w)
{
    Arg a;
    Widget warg;
    int i;
    XtSetArg(a, XmNmenuHistory, &warg);
    XtGetValues(w[0], &a, 1);
    i = 0;
    while ((w[i + 2] != warg) && i < MAXITEMS)
	i++;
    return i;
}

Widget CreateTextItem(Widget parent, int x, int y, int len, char *s)
{
    static Widget w;
    XmString str;
    char buf[256];
    Dimension ww, wh;

    strcpy(buf, "text");
    str = XmStringCreateLtoR(s, XmSTRING_DEFAULT_CHARSET);
    w = XtVaCreateManagedWidget(buf, xmLabelGadgetClass, parent,
				XmNlabelString, str,
				NULL);
    w = XtVaCreateManagedWidget(buf, xmTextWidgetClass, parent,
				XmNtraversalOn, True,
				XmNcolumns, len,
				NULL);
    return w;
}

Widget CreateTextItem2(Widget parent, int len, char *s)
{
    static Widget w;
    Widget rc;
    XmString str;
    char buf[256];
    Arg a;

    strcpy(buf, "text");
    XtSetArg(a, XmNorientation, XmHORIZONTAL);
    rc = XmCreateRowColumn(parent, buf, &a, 1);
    strcpy(buf, "text");
    str = XmStringCreateLtoR(s, XmSTRING_DEFAULT_CHARSET),
    w = XtVaCreateManagedWidget(buf, xmLabelGadgetClass, rc,
				XmNlabelString, str,
				NULL);
    strcpy(buf, "text");
    w = XtVaCreateManagedWidget(buf, xmTextWidgetClass, rc,
				XmNtraversalOn, True,
				XmNcolumns, len,
				NULL);
    XtManageChild(rc);
    return w;
}

Widget CreateTextItem3(Widget parent, int len, char *s, Widget * lab)
{
    static Widget w;
    Widget rc;
    XmString str;
    char buf[256];
    Arg a;

    strcpy(buf, "text");
    XtSetArg(a, XmNorientation, XmHORIZONTAL);
    rc = XmCreateRowColumn(parent, buf, &a, 1);
    strcpy(buf, "text");
    str = XmStringCreateLtoR(s, XmSTRING_DEFAULT_CHARSET),
    *lab = XtVaCreateManagedWidget(buf, xmLabelGadgetClass, rc,
				   XmNlabelString, str,
				   NULL);
    strcpy(buf, "text");
    w = XtVaCreateManagedWidget(buf, xmTextWidgetClass, rc,
				XmNtraversalOn, True,
				XmNcolumns, len,
				NULL);
    XtManageChild(rc);
    return w;
}

char *xv_getstr(Widget w)
{
    char *s;
    static char buf[256];

    strcpy(buf, s = XmTextGetString(w));
    XtFree(s);
    return buf;
}

void xv_setstr(Widget w, char *s)
{
    XmTextSetString(w, s);
}

/*
 * set position
 */
void XmGetPos(Widget w, int type, int *x, int *y)
{
    Arg a[4];
    int ac = 0;
    static int curpos = 0;
    Position xt, yt;
    Dimension wt, ht;
/*
 * get the location of the whole app
 */
    ac = 0;
    XtSetArg(a[ac], XmNx, &xt);
    ac++;
    XtSetArg(a[ac], XmNy, &yt);
    ac++;
    XtSetArg(a[ac], XmNheight, &ht);
    ac++;
    XtSetArg(a[ac], XmNwidth, &wt);
    ac++;
    XtGetValues(w, a, ac);

    switch (type) {
    case 0:			/* position at upper left */
	*x = xt - 50;
	*y = yt - 50 + (curpos % 5) * 50;
	curpos++;
	if (*x < 0) {
	    *x = 0;
	}
	if (*y < 0) {
	    *y = 0;
	}
	break;
    case 1:			/* position at upper right */
	*x = xt + 50 + wt;
	*y = yt + (curpos % 5 + 1) * 50;
	break;
    case 2:			/* center on w */
	*x = xt + wt / 2;
	*y = yt + ht / 2;
	break;
    }
}

void GetWH(Widget w, Dimension * ww, Dimension * wh)
{
    Arg args;
    XtSetArg(args, XmNwidth, ww);
    XtGetValues(w, &args, 1);
    XtSetArg(args, XmNheight, wh);
    XtGetValues(w, &args, 1);
}

void GetXY(Widget w, Position * x, Position * y)
{
    Arg args;
    XtSetArg(args, XmNx, x);
    XtGetValues(w, &args, 1);
    XtSetArg(args, XmNy, y);
    XtGetValues(w, &args, 1);
}

/**********************************************************************
 * XtFlush - Flushes all Xt events.
 **********************************************************************/
void XtFlush(void)
{
    extern Widget app_shell;
    XtAppContext app;
    app = XtWidgetToApplicationContext(app_shell);
    while (XtAppPending(app) & XtIMXEvent) {
	XtAppProcessEvent(app, XtIMXEvent);
	XFlush(XtDisplay(app_shell));
    }
}

static void set_parent_color(Widget w, int cd)
{
    extern long colors[];
/*
    XtVaSetValues(XtParent(w),
		  XmNbackground, colors[cd],
		  NULL);
*/

}

Widget *CreateColorChoice(Widget parent, char *s, int map)
{
    va_list var;
    Arg args[MAXARGS];
    int nchoices, ncols, err = 0, nargs, i = 0, managed = 1;
    XmString str;
    String argstr;
    XtArgVal argval;
    char *name, labelstr[32], buf[32];
    WidgetClass class;
    Widget lab, w[MAXITEMS], *retval;
    extern long colors[];

    if (s == NULL) {
	strcpy(labelstr, "");
    } else {
	strcpy(labelstr, s);
    }

    ncols = 4;
    nchoices = 16;
    retval = (Widget *) XtMalloc((16 + 2) * sizeof(Widget));

    strcpy(buf, "text");
    name = buf;
    nargs = 0;
    XtSetArg(args[nargs], XmNorientation, XmVERTICAL);
    nargs++;
    XtSetArg(args[nargs], XmNpacking, XmPACK_COLUMN);
    nargs++;
    XtSetArg(args[nargs], XmNnumColumns, ncols);
    nargs++;
    w[1] = XmCreatePulldownMenu(parent, name, args, nargs);
    for (i = 0; i < 16; i++) {
	sprintf(buf, "%2d", i);
	w[i + 2] = XmCreatePushButton(w[1], buf, NULL, 0);
	XtVaSetValues(w[i + 2],
		      XmNbackground, colors[i + (map - 1) * 16],
		      NULL);
    }
    XtManageChildren(w + 2, nchoices);

    strcpy(buf, "text");
    name = buf;

    str = XmStringCreate(labelstr, XmSTRING_DEFAULT_CHARSET);

    nargs = 0;
    XtSetArg(args[nargs], XmNlabelString, str);
    nargs++;
    XtSetArg(args[nargs], XmNsubMenuId, w[1]);
    nargs++;
    XtSetArg(args[nargs], XmNentryBorder, 2);
    nargs++;
    XtSetArg(args[nargs], XmNwhichButton, 1);
    nargs++;
    w[0] = XmCreateOptionMenu(parent, name, args, nargs);
    XtManageChild(w[0]);
    for (i = 0; i < nchoices + 2; i++) {
	retval[i] = w[i];
    }
    for (i = 0; i < nchoices; i++) {
	XtAddCallback(w[i + 2], XmNactivateCallback, (XtCallbackProc) set_parent_color, (XtPointer) ((long) i));
    }
    return retval;
}

Widget *CreatePanelChoice1(Widget parent, char *labelstr, int nchoices,...)
{
    va_list var;
    Arg args[MAXARGS];
    int err = 0, nargs, i = 0, managed = 1;
    String argstr;
    XmString str;
    XtArgVal argval;
    char *name, *s, buf[32];
    WidgetClass class;
    Widget lab, w[MAXITEMS], *retval;

    nchoices--;			/* always one too many */
    retval = (Widget *) XtMalloc((nchoices + 2) * sizeof(Widget));

    strcpy(buf, "choice");
    name = buf;
    w[1] = XmCreatePulldownMenu(parent, name, NULL, 0);

    va_start(var, nchoices);
    i = 0;
    while ((s = va_arg(var, char *)) != NULL && i < nchoices) {
	w[i + 2] = XmCreatePushButtonGadget(w[1], s, NULL, 0);
	i++;
    }
    XtManageChildren(w + 2, nchoices);
    va_end(var);

    str = XmStringCreate(labelstr, XmSTRING_DEFAULT_CHARSET);

    nargs = 0;
    XtSetArg(args[nargs], XmNlabelString, str);
    nargs++;
    XtSetArg(args[nargs], XmNsubMenuId, w[1]);
    nargs++;
    XtSetArg(args[nargs], XmNentryBorder, 2);
    nargs++;
    XtSetArg(args[nargs], XmNwhichButton, 1);
    nargs++;
    w[0] = XmCreateOptionMenu(parent, name, args, nargs);
    XtManageChild(w[0]);
    for (i = 0; i < nchoices + 2; i++) {
	retval[i] = w[i];
    }
    return retval;
}

Widget *CreatePanelChoice2(Widget parent, char *labelstr, int ncols, int nchoices,...)
{
    va_list var;
    Arg args[MAXARGS];
    int err = 0, nargs, i = 0, managed = 1;
    XmString str;
    String argstr;
    XtArgVal argval;
    char *name, *s, buf[32];
    WidgetClass class;
    Widget lab, w[MAXITEMS], *retval;

    nchoices--;			/* always one too many */
    retval = (Widget *) XtMalloc((nchoices + 2) * sizeof(Widget));

    strcpy(buf, "choice");
    name = buf;
    nargs = 0;
    XtSetArg(args[nargs], XmNorientation, XmVERTICAL);
    nargs++;
    XtSetArg(args[nargs], XmNpacking, XmPACK_COLUMN);
    nargs++;
    XtSetArg(args[nargs], XmNnumColumns, ncols);
    nargs++;
    w[1] = XmCreatePulldownMenu(parent, name, args, nargs);

    i = 0;
    va_start(var, nchoices);
    while ((s = va_arg(var, char *)) != NULL && i < nchoices) {
	w[i + 2] = XmCreatePushButtonGadget(w[1], s, NULL, 0);
	i++;
    }
    XtManageChildren(w + 2, nchoices);
    va_end(var);

    str = XmStringCreate(labelstr, XmSTRING_DEFAULT_CHARSET);

    nargs = 0;
    XtSetArg(args[nargs], XmNlabelString, str);
    nargs++;
    XtSetArg(args[nargs], XmNsubMenuId, w[1]);
    nargs++;
    XtSetArg(args[nargs], XmNentryBorder, 2);
    nargs++;
    XtSetArg(args[nargs], XmNwhichButton, 1);
    nargs++;
    w[0] = XmCreateOptionMenu(parent, name, args, nargs);
    XtManageChild(w[0]);
    for (i = 0; i < nchoices + 2; i++) {
	retval[i] = w[i];
    }
    return retval;
}

/*
 * generic unmanage popup routine, used elswhere
 */
void destroy_dialog(Widget w, Widget p)
{
    XtUnmanageChild(p);
}

/*
 * handle the close item on the WM menu
 */
void handle_close(Widget w)
{
    Atom WM_DELETE_WINDOW;
    XtVaSetValues(w,
		  XmNdeleteResponse, XmDO_NOTHING,
		  NULL);
    WM_DELETE_WINDOW = XmInternAtom(XtDisplay(app_shell), "WM_DELETE_WINDOW", False);
    XmAddWMProtocolCallback(w, WM_DELETE_WINDOW, (XtCallbackProc) destroy_dialog, (XtPointer) w);
}

Widget CreateTextItem0(Widget parent, int len, char *s, Widget * label)
{
    static Widget w;
    XmString str;
    char buf[256];
    Dimension ww, wh;

    strcpy(buf, "text");
    str = XmStringCreateLtoR(s, XmSTRING_DEFAULT_CHARSET);
    *label = XtVaCreateManagedWidget(buf, xmLabelGadgetClass, parent,
				     XmNlabelString, str,
				     NULL);
    w = XtVaCreateManagedWidget(buf, xmTextWidgetClass, parent,
				XmNtraversalOn, True,
				XmNcolumns, len,
				NULL);
    return w;
}

void panel_setstr_value(Widget w, char *s)
{
    XmTextSetString(w, s);
}

char *panel_getstr_value(Widget w)
{
    char *s;
    static char buf[256];
    strcpy(buf, s = XmTextGetString(w));
    XtFree(s);
    return buf;
}

int check_action(void)
{
    extern XtAppContext app_con;
    return XtAppPending(app_con);
}

void cancel_action(void)
{
    extern XtAppContext app_con;
    XEvent event;
    while (XtAppPending(app_con)) {
	XtAppNextEvent(app_con, &event);
	if (event.type == ButtonPress && event.xbutton.button == Button1) {
	    extern Window abort_win;
	    if (event.xbutton.window == abort_win) {
		XtAppProcessEvent(app_con, XtIMXEvent);
	    }
	}
    }
}

/*
 * Manage and raise
 */
void XtRaise(Widget w)
{
    XtManageChild(w);
    XMapRaised(XtDisplay(w), XtWindow(w));
}

void CreateCommandButtons(Widget parent, int n, Widget * buts, char **labels)
{
    int i;
    Widget form;
    form = XtVaCreateWidget("form", xmFormWidgetClass, parent,
                            XmNfractionBase, n,
                            NULL);
 
    for (i = 0; i < n; i++) {
        buts[i] = XtVaCreateManagedWidget(labels[i],
                                          xmPushButtonWidgetClass, form,
                                          XmNtopAttachment, XmATTACH_FORM,
                                          XmNbottomAttachment, XmATTACH_FORM,
                                       XmNleftAttachment, XmATTACH_POSITION,
                                          XmNleftPosition, i,
                                      XmNrightAttachment, XmATTACH_POSITION,
                                          XmNrightPosition, i + 1,
                                    XmNshowAsDefault, i == 0 ? True : False,
                                          XmNdefaultButtonShadowThickness, 1,
                                          NULL);
    }
    XtManageChild(form);
    {
        Dimension h;
        XtVaGetValues(buts[0], XmNheight, &h, NULL);
        XtVaSetValues(form, XmNpaneMaximum, h, XmNpaneMinimum, h, NULL);
    }
}
