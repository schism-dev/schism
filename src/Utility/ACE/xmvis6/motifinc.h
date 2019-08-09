/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/* $Id: motifinc.h,v 1.2 2003/07/24 15:44:06 pturner Exp $
 *
 * includes for motif
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <sys/types.h>
#ifndef WIN32
#include <sys/time.h>
#include <sys/signal.h>
#else
#include <time.h>
#include <signal.h>
#endif

#include <X11/X.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <X11/Shell.h>
#include <Xm/Xm.h>
#include <Xm/ArrowB.h>
#include <Xm/BulletinB.h>
#include <Xm/CascadeB.h>
#include <Xm/CascadeBG.h>
#include <Xm/Command.h>
#include <Xm/DrawingA.h>
#include <Xm/DialogS.h>
#include <Xm/FileSB.h>
#include <Xm/Form.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/LabelG.h>
#include <Xm/List.h>
#include <Xm/MainW.h>
#include <Xm/MenuShell.h>
#include <Xm/MessageB.h>
#include <Xm/PanedW.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/RowColumn.h>
#include <Xm/Scale.h>
#include <Xm/ScrollBar.h>
#include <Xm/ScrolledW.h>
#include <Xm/SelectioB.h>
#include <Xm/Separator.h>
#include <Xm/SeparatoG.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/Text.h>
#include <Xm/TextF.h>

Widget *CreatePanelChoice1(Widget w, char *lab, int count, ...);
Widget *CreatePanelChoice2(Widget w, char *lab, int ncols, int count, ...);
Widget *CreateColorChoice(Widget parent, char *s, int map);
Widget CreateTextItem(Widget parent, int x, int y, int len, char *s);
Widget CreateTextItem2(Widget parent, int len, char *s);
Widget CreateTextItem3(Widget parent, int len, char *s, Widget * lab);
void destroy_dialog(Widget w, XtPointer clientd, XtPointer calld);
void do_browser(Widget w, XtPointer client_data, XtPointer call_data);
char *xv_getstr(Widget w);

extern XmStringCharSet charset;
extern Widget app_shell;

/*
 * a file popup and a text item
 */
typedef struct {
    Widget file;
    Widget browser;
} Browser;
