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
 * Read Sample TEA-NL flows
 */

#ifndef lint
static char RCSid[] = "$Id: sampflowwin.c,v 1.11 2007/03/12 21:25:06 pturner Exp $";
#endif

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern XmStringCharSet charset;
Widget CreateTextItem(Widget parent, int x, int y, int len, char *s);

Widget sampflow_frame;
Widget sampflow_panel;

void create_sampflow_frame(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc sampflow_accept_proc(Widget w, XtPointer clid, XtPointer calld);
XtCallbackProc sampflow_pickxy_proc(Widget w, XtPointer clid, XtPointer calld);
void ReadSampleFlow(DisplayFlow * f, const char *s);
void SetMinSampleFlow(int flowno, double mindist);
void ReadSampleFlowXY(DisplayFlow * f, const char *s);
void WriteSampleFlow(DisplayFlow * f, const char *s);
void AddSampleFlowNode(DisplayFlow * f, int node);
void AddSampleFlowElem(DisplayFlow * f, int elem);
void AddSampleFlowXY(DisplayFlow * f, double x, double y);
void DeleteSampleFlowNode(DisplayFlow * f, int node);
void ClearSampleFlowNode(DisplayFlow * f);
void DeleteSampleFlowElem(DisplayFlow * f, int elem);
void DeleteSampleFlowXY(DisplayFlow * f, double x, double y);
void ShowSampleLocations(DisplayFlow * f, Grid * g);

/* used for communication with event handler to define nodes */
static DisplayFlow *f;

/* Read a list of nodes to sample */
typedef struct {
    Widget top;
    Widget file;
    Widget *flow;
} SampReadUI;

XtCallbackProc sampflow_read(Widget w, XtPointer clid, XtPointer calld);

XtCallbackProc sampflow_read_proc(Widget w, XtPointer clid, XtPointer calld)
{
    static SampReadUI ui;
    Widget panel;
    if (ui.top != (Widget) NULL) {
	XtRaise(ui.top);
	return;
    }
    ui.top = XmCreateFileSelectionDialog(app_shell, "Read Samples", NULL, 0);
    panel = XmCreateRowColumn(ui.top, "panel", NULL, 0);
    ui.flow = CreatePanelChoice1(panel, "For ADCIRC flow: ", 6, "1", "2", "3", "4", "5", 0, 0);
    XtManageChild(panel);
    XtAddCallback(ui.top, XmNcancelCallback, (XtCallbackProc) destroy_dialog, (XtPointer) ui.top);
    XtAddCallback(ui.top, XmNokCallback, (XtCallbackProc) sampflow_read, (XtPointer) & ui);
    XtRaise(ui.top);
}

XtCallbackProc sampflow_read(Widget w, XtPointer clid, XtPointer calld)
{
    Arg args;
    XmString list_item;
    char *s;
    int flowno;
    DisplayFlow *f;
    SampReadUI *ui = (SampReadUI *) clid;
    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(ui->top, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    flowno = GetChoice(ui->flow);
    f = &g[cg].flowt[flowno];
    ReadSampleFlowXY(f, s);
    XtFree(s);
}

/* Write a list of nodes to sample */
typedef struct {
    Widget top;
    Widget file;
    Widget *flow;
} SampWriteUI;

static void sampflow_write(Widget w, XtPointer clid, XtPointer calld);

XtCallbackProc sampflow_write_proc(Widget w, XtPointer clid, XtPointer calld)
{
    static SampWriteUI ui;
    int x, y;
    Widget panel, wbut, rc;

    int n = 0;
    if (ui.top == (Widget) NULL) {
	ui.top = XmCreateDialogShell(app_shell, "Write samples", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "panel", NULL, 0);
	ui.flow = CreatePanelChoice1(panel, "From ADCIRC flow: ", 6, "1", "2", "3", "4", "5", 0, 0);
	ui.file = CreateTextItem2(panel, 20, "Write parameters to: ");
	XtVaCreateManagedWidget("sep", xmSeparatorGadgetClass, panel, NULL);
	rc = XmCreateRowColumn(panel, "rc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);

	wbut = XtVaCreateManagedWidget("Accept", xmPushButtonGadgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_write, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonGadgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) ui.top);
	XtManageChild(rc);

	XtManageChild(panel);
    }
    XtRaise(ui.top);
}

static void sampflow_write(Widget w, XtPointer clid, XtPointer calld)
{
    int i, flowno;
    char s[1024];
    DisplayFlow *f;
    SampWriteUI *ui = (SampWriteUI *) clid;
    strcpy(s, (char *) xv_getstr(ui->file));
    if (strlen(s) > 0) {
	if (!fexists(s)) {
	    flowno = GetChoice(ui->flow);
	    f = &g[cg].flowt[flowno];
	    WriteSampleFlow(f, s);
	}
    } else {
	errwin("No output file defined, enter or select filename");
    }
}

typedef struct {
    Widget top;
    Widget *flow;
    Widget display;
    Widget displayLocations;
} SampUI;

XtCallbackProc sampflow_accept_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SampUI *ui = (SampUI *) clid;
    int flowno = GetChoice(ui->flow);
    g[cg].flowt[flowno].sample = XmToggleButtonGetState(ui->display) ? ON : OFF;
}

XtCallbackProc sampflow_clear_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SampUI *ui = (SampUI *) clid;
    int flowno = GetChoice(ui->flow);
    DisplayFlow *f = &g[cg].flowt[flowno];
    ClearSampleFlowNode(f);
}

void add_sample_node(int node)
{
    AddSampleFlowNode(f, node);
}

void del_sample_node(int node)
{
}

void pick_samp_node(void)
{
    set_restart();
    set_action(0);
    set_action(SAMPLE_ADCIRC);
}

void del_samp_node(void)
{
    set_restart();
    set_action(0);
    set_action(DEL_SAMPLE_ADCIRC);
}

XtCallbackProc sampflow_pick_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SampUI *ui = (SampUI *) clid;
    int flowno = GetChoice(ui->flow);
    f = &g[cg].flowt[flowno];
    pick_samp_node();
}

XtCallbackProc sampflow_delete_proc(Widget w, XtPointer clid, XtPointer calld)
{
    SampUI *ui = (SampUI *) clid;
}

void update_sampflow(SampUI * ui)
{
    int flowno = GetChoice(ui->flow);
    XmToggleButtonSetState(ui->display, g[cg].flowt[flowno].sample == ON, False);
}

void create_sampflow_frame(Widget w, XtPointer clid, XtPointer calld)
{
    static SampUI ui;
    Widget frame, panel, wbut, fr, bb, sep, rc;
    setistop();
    if (!ui.top) {
	ui.top = XmCreateDialogShell(app_shell, "Sample Flow setup", NULL, 0);
	panel = XmCreateRowColumn(ui.top, "sampflowrc", NULL, 0);
	XtVaSetValues(panel, XmNorientation, XmVERTICAL, NULL);

	ui.display = XtVaCreateManagedWidget("Use samples", xmToggleButtonWidgetClass, panel, NULL);
	ui.displayLocations = XtVaCreateManagedWidget("Display locations", xmToggleButtonWidgetClass, panel, NULL);
	ui.flow = CreatePanelChoice1(panel, "Apply to ADCIRC flow: ", 6, "1", "2", "3", "4", "5", 0, 0);

	rc = XmCreateRowColumn(panel, "sampflowrc", NULL, 0);
	XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
	wbut = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_read_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Write...", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_write_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_pick_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("*Delete", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_delete_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Clear", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_clear_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Apply", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) sampflow_accept_proc, (XtPointer) & ui);
	wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
	XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, ui.top);
	XtManageChild(rc);
	XtManageChild(panel);
    }
    XtRaise(ui.top);
    update_sampflow(&ui);
}

/*
 * Sample velocities to cut down the clutter
 */
void ReadSampleFlow(DisplayFlow * f, const char *s)
{
    FILE *fp;
    char buf[256];
    int node;
    if (s != NULL) {
	if ((fp = fopen(s, "rb")) != NULL) {
	    while (fgets(buf, 255, fp) != NULL) {
		sscanf(buf, "%d", &node);
		AddSampleFlowNode(f, node);
	    }
	    fclose(fp);
	    f->samptype = NODE;
	} else {
	    /* couldn't open file */
	}
    } else {
	/* Filename error */
    }
}

/*
 * Sample velocities to cut down the clutter
 */
void SetMinSampleFlow(int flowno, double mindis)
{
    int i, j, cnt, *imap, *icol;
    double dis;
    DisplayFlow * f = &g[cg].flowt[flowno];
    int np = flowt[flowno].g.nmnp;
    double *x = flowt[flowno].g.xord, *y = flowt[flowno].g.yord;
    cnt = 1;
    imap = (int *) malloc(np * sizeof(int));
    icol = (int *) malloc(np * sizeof(int));
    imap[0] = 0;
    for (i=0;i<np;i++) {
        icol[i] = 0;
    }
    for (i=0;i<np -1;i++) {
        if (icol[i] == 0) {
            icol[i] = 2;
            for (j=i + 1;j<np;j++) {
                 if (icol[j] == 0 ) {
                     dis=hypot(x[i]-x[j], y[i]-y[j]);
                     if (dis < mindis) {
                       icol[j]=1;
                     }
                 }
             }
        }
    }
    for (i=0;i<np;i++) {
        if (icol[i] == 2) {
            AddSampleFlowNode(f, i);
        }
    }
    free(imap);
    free(icol);
}

/*
 * Sample velocities to cut down the clutter
 */
void ReadSampleFlowXY(DisplayFlow * f, const char *s)
{
    FILE *fp;
    char buf[256];
    int i, n;
    double x, y, flag;
    if (s != NULL) {
	if ((fp = fopen(s, "rb")) != NULL) {
	    if (fgets(buf, 255, fp) == NULL) {
		return;
	    }
	    if (fgets(buf, 255, fp) == NULL) {
		return;
	    }
	    sscanf(buf, "%d", &n);
	    for (i = 0; i < n; i++) {
		if (fgets(buf, 255, fp) != NULL) {
		    sscanf(buf, "%*d %lf %lf %lf", &x, &y, &flag);
		    if ((int) flag != 0) {
			AddSampleFlowXY(f, x, y);
		    }
		} else {
		    break;
		}
	    }
	    fclose(fp);
	    f->samptype = XY;
	} else {
	    /* couldn't open file */
	}
    } else {
	/* Filename error */
    }
}

void WriteSampleFlow(DisplayFlow * f, const char *s)
{
    FILE *fp;
    char buf[256];
    int i;
    if (s != NULL) {
	if ((fp = fopen(s, "wb")) != NULL) {
	    for (i = 0; i < f->nsamples; i++) {
		fprintf(fp, "%d\n", f->samples[i]);
	    }
	    fclose(fp);
	} else {
	    /* couldn't open file */
	}
    } else {
	/* Filename error */
    }
}

void AddSampleFlowNode(DisplayFlow * f, int node)
{
    int i;
    if (node < 0) {
	errwin("Sample node < 0");
    }
    if (f->nsamples == 0) {
	f->samples = (int *) malloc(sizeof(int));
    } else {
	for (i = 0; i < f->nsamples; i++) {
	    if (f->samples[i] == node) {	/* aleady in list */
		return;
	    }
	}
	f->samples = (int *) realloc(f->samples, (f->nsamples + 1) * sizeof(int));
    }
    if (f->samples != NULL) {
	f->samples[f->nsamples] = node;
	f->nsamples++;
    } else {
	errwin("Samples array is NULL");
    }
}

void AddSampleFlowElem(DisplayFlow * f, int elem)
{
}

void AddSampleFlowXY(DisplayFlow * f, double x, double y)
{
    int i;
    if (f->nsamples == 0) {
	f->sampx = (double *) malloc(sizeof(double));
	f->sampy = (double *) malloc(sizeof(double));
    } else {
	f->sampx = (double *) realloc(f->sampx, (f->nsamples + 1) * sizeof(double));
	f->sampy = (double *) realloc(f->sampy, (f->nsamples + 1) * sizeof(double));
    }
    if (f->sampx != NULL || f->sampy != NULL) {
	f->sampx[f->nsamples] = x;
	f->sampy[f->nsamples] = y;
	f->nsamples++;
    } else {
	errwin("Samples array is NULL");
    }
}

void DeleteSampleFlowNode(DisplayFlow * f, int node)
{
    int *tmp, i, cnt = 0;
    if (f->nsamples == 0) {
    } else {
	tmp = (int *) malloc(f->nsamples * sizeof(int));
	for (i = 0; i < f->nsamples; i++) {
	    if (f->samples[i] != node) {
		tmp[cnt] = f->samples[i];
		cnt++;
	    }
	}
    }
    free(f->samples);
    f->samples = tmp;
    f->nsamples = cnt;
}

void ClearSampleFlowNode(DisplayFlow * f)
{
    if (f->nsamples == 0) {
	return;
    } else {
	f->nsamples = 0;
	free(f->samples);
    }
}

void ClearSampleFlowXY(DisplayFlow * f)
{
    if (f->nsamples == 0) {
	return;
    } else {
	f->nsamples = 0;
	free(f->sampx);
	free(f->sampy);
    }
}

void DeleteSampleFlowElem(DisplayFlow * f, int elem)
{
}

void DeleteSampleFlowXY(DisplayFlow * f, double x, double y)
{
}

void ShowSampleLocations(DisplayFlow * f, Grid * g)
{
    int i, node;
    double x, y;
    for (i = 0; i < f->nsamples; i++) {
	node = f->samples[i];
	if (node >= 0 && node < g->nmnp) {
	    x = g->xord[f->samples[i]];
	    y = g->yord[f->samples[i]];
	    my_circle(x, y);
	}
    }
}

void LoadSampleLocations(DisplayFlow * f, Grid * g)
{
    int i, node;
    if (f->samptype == XY) {
	f->samples = (int *) malloc(f->nsamples * sizeof(int));
	if (f->nsamples > 0 && f->sampx != NULL && f->sampy != NULL) {
	    for (i = 0; i < f->nsamples; i++) {
		FindNearestNode(g, f->sampx[i], f->sampy[i], &node);
		f->samples[i] = node;
	    }
	}
    }
}

int DisplaySample(DisplayFlow * f, int node)
{
    int i;
    if (f->sample == OFF) {
	return 1;		/* if off, then always display */
    }
    for (i = 0; i < f->nsamples; i++) {
	if (f->samples[i] == node) {
	    return 1;
	}
    }
    return 0;
}

int DisplaySampleXY(DisplayFlow * f, double x, double y)
{
    int i;
    double tmp;
    if (f->sample == OFF) {
	return 1;		/* if off, then always display */
    }
    for (i = 0; i < f->nsamples; i++) {
    }
    return 0;
}
