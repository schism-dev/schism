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
 * Help
 *
 */

#ifndef lint
static char RCSid[] = "$Id: helpwin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "motifinc.h"

extern XmStringCharSet charset;

static void clear_results(void);

static Widget help_frame, help_panel, text_w;

/*
 * Create the help Panel
 */
void create_help_frame(Widget w, int cd)
{
    int x, y;
    Widget wbut, rc, fr;
    extern Widget app_shell;
    struct stat sb;
    FILE *fp;
    char *t, hname[256];

    if (!help_frame) {
	XmGetPos(app_shell, 0, &x, &y);
	help_frame = XmCreateDialogShell(app_shell, "Help", NULL, 0);
	handle_close(help_frame);
	XtVaSetValues(help_frame, XmNx, x, XmNy, y, NULL);
	help_panel = XmCreateForm(help_frame, "help_form", NULL, 0);
	fr = XmCreateFrame(help_panel, "fr", NULL, 0);
	text_w = XmCreateScrolledText(fr, "text_w", NULL, 0);
	XtVaSetValues(text_w,
		      XmNrows, 10,
		      XmNcolumns, 40,
		      XmNeditable, False,
		      XmNeditMode, XmMULTI_LINE_EDIT,
		      XmNwordWrap, True,
		      NULL);
	XtManageChild(text_w);
	XtManageChild(fr);

	rc = XmCreateRowColumn(help_panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Close", xmPushButtonGadgetClass, rc,
				       NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, help_frame);
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

	XtManageChild(help_panel);
    }
    XtRaise(help_frame);
    switch (cd) {
    case 10:
	strcpy(hname, "adcirc.help");
	break;
    case 11:
	strcpy(hname, "rma2.help");
	break;
    }
    /*if (stat(hname, &sb) == -1 || !S_ISREG(sb.st_mode) || TODO */
    if (stat(hname, &sb) == -1 || !(fp = fopen(hname, "r"))) {
    } else {
	if (!(t = XtMalloc((unsigned) (sb.st_size + 1)))) {
	    return;
	}
	if (!fread(t, sizeof(char), sb.st_size + 1, fp)) {
	}
	t[sb.st_size] = 0;
	XmTextSetString(text_w, t);
	XtFree(t);
	fclose(fp);
    }
}

static void clear_results(void)
{
    XmTextSetString(text_w, "");
}

void stuffhelp(char *s, int sp)
{
    extern int inwin;
    static XmTextPosition pos = 0;
    static XmTextPosition savepos = 0;

    if (inwin) {
	create_help_frame(NULL, 0);
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
