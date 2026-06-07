/* $Id: pagewin.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * Set page size and orientation
 */
#include <stdio.h>
#include <math.h>

#include <Xm/Xm.h>
#include <Xm/BulletinB.h>
#include <Xm/DialogS.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/ToggleB.h>
#include <Xm/Scale.h>
#include <Xm/Separator.h>

#include "globals.h"
#include "motifinc.h"

extern int canvasw, canvash;

static Widget page_frame;
static Widget page_panel;

/*
 * Panel item declarations
 */
static Widget *page_layout_item;
static Widget page_width_item;
static Widget page_height_item;

/*
 * Event and Notify proc declarations
 */

void update_page(void)
{
    char buf[256];
    int ww, wh;
    if (page_frame) {
	SetChoice(page_layout_item, get_pagelayout());
	sprintf(buf, "%d", canvasw);
	xv_setstr(page_width_item, buf);
	sprintf(buf, "%d", canvash);
	xv_setstr(page_height_item, buf);
    }
}

/*
 * define the draw options
 */
static void define_page_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char buf[256];
    int ww, wh;
    set_page(NULL, (XtPointer) get_pagelayout(), NULL);
    switch (GetChoice(page_layout_item)) {
        case 1:
        page_layout = LANDSCAPE;
        break;
    case 2:
        page_layout = PORTRAIT;
        break;
    case 3:
        page_layout = FIXED;
        break;
    case 0:                  /* falls through */
    default:
        page_layout = FREE;
        break;
    }

    if (page_layout == FIXED) {
        strcpy(buf, (char *) xv_getstr(page_width_item));
        canvasw = atoi(buf);
        strcpy(buf, (char *) xv_getstr(page_height_item));
        canvash = atoi(buf);
    }
    set_page(NULL, (XtPointer) page_layout, NULL);
    XtUnmanageChild(page_frame);
    drawgraph();
}

/*
 * Create the draw Frame and the draw Panel
 */
void create_page_frame(Widget w, XtPointer client_data, XtPointer call_data)
{
    int x, y;
    Widget wbut, rc, buts[2];
    Widget wlabel;

    set_wait_cursor();
    if (page_frame == NULL) {
	char *label1[2];
	label1[0] = "Accept";
	label1[1] = "Close";
	XmGetPos(app_shell, 0, &x, &y);
	page_frame = XmCreateDialogShell(app_shell, "Page size", NULL, 0);
	handle_close(page_frame);
	XtVaSetValues(page_frame,
		      XmNx, x,
		      XmNy, y,
		      NULL);
	page_panel = XmCreateRowColumn(page_frame, "page_rc", NULL, 0);
	page_layout_item = CreatePanelChoice(page_panel, "Page layout",
					     5,
					     "Free",
					     "Landscape",
					     "Portrait",
					     "Fixed",
					     NULL, NULL);
	page_width_item = CreateTextItem2(page_panel, 10, "Page width (pixels)");
	page_height_item = CreateTextItem2(page_panel, 10, "Page height (pixels)");

	XtVaCreateManagedWidget("sep", xmSeparatorWidgetClass, page_panel, NULL);

	CreateCommandButtons(page_panel, 2, buts, label1);
	XtAddCallback(buts[0], XmNactivateCallback,
		      (XtCallbackProc) define_page_proc, (XtPointer) 0);
	XtAddCallback(buts[1], XmNactivateCallback,
		   (XtCallbackProc) destroy_dialog, (XtPointer) page_frame);

	XtManageChild(page_panel);
    }
    XtRaise(page_frame);
    update_page();
    unset_wait_cursor();
}
