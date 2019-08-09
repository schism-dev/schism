/* $Id: blockwin.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * read block data files
 *
 */

#include <stdio.h>

#include <Xm/Xm.h>
#include <Xm/DialogS.h>
#include <Xm/BulletinB.h>
#include <Xm/FileSB.h>
#include <Xm/Frame.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include <Xm/RowColumn.h>
#include <Xm/SelectioB.h>
#include <Xm/ToggleB.h>

#include "globals.h"
#include "motifinc.h"

static Widget block_dialog;	/* read data popup */

extern double *blockdata[];	/* TODO => globals.h */
extern int blocklen;
extern int blockncols;

static int blocksrc = DISK;

static void set_src_proc(Widget w, XtPointer client_data, XtPointer call_data);
static void block_proc(Widget w, XtPointer client_data, XtPointer call_data);

static void set_src_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    int data = (int) client_data;
    switch (data) {
    case 0:
	blocksrc = DISK;
	break;
    case 1:
	blocksrc = PIPE;
	break;
    }
}

static void block_proc(Widget w, XtPointer client_data, XtPointer call_data)
{
    char *s;
    XmFileSelectionBoxCallbackStruct *cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
    if (!XmStringGetLtoR (cbs->value, charset, &s)) {
        errwin("Error converting XmString to char string in block_proc()");
        return;
    }
    if (getdata(cg, s, blocksrc, BLOCK)) {
	if (blocklen == 0) {
	    errwin("Block data length = 0");
	} else if (blockncols == 0) {
	    errwin("Number of columns in block data = 0");
	} else {
	    XtUnmanageChild(block_dialog);
	    create_eblock_frame(NULL, NULL, NULL);
	}
    }
    XtFree(s);
}

void create_block_popup(Widget w, XtPointer client_data, XtPointer call_data)
{
    int i;
    Widget lab, rc, fr, rb, rw[5];

    set_wait_cursor();
    if (block_dialog == NULL) {
	block_dialog = XmCreateFileSelectionDialog(app_shell, "read_block_data", NULL, 0);
	XtVaSetValues(XtParent(block_dialog), XmNtitle, "Read block data", NULL);

	XtAddCallback(block_dialog, XmNcancelCallback, (XtCallbackProc) destroy_dialog, block_dialog);
	XtAddCallback(block_dialog, XmNokCallback, (XtCallbackProc) block_proc, 0);

	fr = XmCreateFrame(block_dialog, "frame", NULL, 0);

	rc = XmCreateRowColumn(fr, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	lab = XmCreateLabel(rc, "Data source:", NULL, 0);
	XtManageChild(lab);

	rb = XmCreateRadioBox(rc, "rb", NULL, 0);
	XtVaSetValues(rb, XmNorientation, XmHORIZONTAL, NULL);

	rw[0] = XmCreateToggleButton(rb, "Disk", NULL, 0);
	rw[1] = XmCreateToggleButton(rb, "Pipe", NULL, 0);
	for (i = 0; i < 2; i++) {
	    XtAddCallback(rw[i], XmNvalueChangedCallback, (XtCallbackProc) set_src_proc, (XtPointer) i);
	}

	XtManageChildren(rw, 2);
	XtManageChild(rb);
	XtManageChild(rc);
	XtManageChild(fr);
	XmToggleButtonSetState(rw[0], True, False);
    }
    XtRaise(block_dialog);
    unset_wait_cursor();
}
