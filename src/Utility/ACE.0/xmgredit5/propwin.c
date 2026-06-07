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
 * Set properties of nodes and elements
 *
 */

#ifndef lint
static char RCSid[] = "$Id: propwin.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget prop_frame;
static Widget prop_panel;
static Widget prop_write_item;
static Widget *prop_type_item;	/* type of prop */
static Widget *prop_propn_item;	/* which prop struct to use */
static Widget *prop_color_item;	/* which color to use */
static Widget prop_default_item;
static Widget prop_current_item;
static Widget display_prop;
static Widget writep_frame;
static Widget writep_panel;
void create_writep_popup(void);
void writep_accept_proc(void);
void drawgrid_prop_filled(int gridno);

static FILE *propfp;
static int propcnt;
static double defval;
static double curval;

typedef struct _Prop_struct {
    int type;			/* nodes, elements, element type, node type,
				 * or just XY */
    int stype;			/* storage type, real or integer */
    int active;
    Isolparms propip;
    int *nodes;
    int *els;
    int *ncolor;
    int *ecolor;
    double *dprop;
    int n;
    double *x, *y;
} Prop;

Prop p[10];
int curprop;
int curcolor;

extern int draw_prop;

void create_prop_popup(void);

void propread_accept_proc(void);

static Widget propread_frame;
static Widget *rprop_choice_item;

void create_readprop_frame(void)
{
    Widget wbut, rc;
    int i;
    Arg a;

    if (!propread_frame) {
	propread_frame = XmCreateFileSelectionDialog(app_shell, "propread_frame", NULL, 0);
	XtVaSetValues(propread_frame,
		      XmNtitle, XmStringCreate("Read property", charset),
		      XmNdirMask, XmStringCreate("*.prop", charset),
		      NULL);

	rc = XmCreateRowColumn(propread_frame, "rc", NULL, 0);

        rprop_choice_item = CreatePanelChoice1(rc, "Read property: ",
					   3,
	   				   "Property at nodes",
					   "Propery at elements",
					   NULL,
					   NULL);

	XtManageChild(rc);
	XtAddCallback(propread_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, propread_frame);
	XtAddCallback(propread_frame, XmNokCallback, (XtCallbackProc) propread_accept_proc, NULL);

    }
    XtRaise(propread_frame);
}

void prop_isolines_proc(void)
{
/*
PBP
*/
    int i,propn;
    curprop = propn = GetChoice(prop_propn_item);
    p[propn].propip.active=1;p[propn].propip.loctype = WORLD;p[propn].propip.p.prec=2;
    p[propn].propip.p.format = DECIMAL;p[propn].propip.p.color = 1;
    create_isolines_popup("Element properties", &(p->propip));
}

void propread_accept_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];
    int flowno, gridno,ptype,propn,i;
    float jf;
    double dp,min,max;
    Widget textw;
    FILE *pf;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(propread_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    set_wait_cursor(propread_frame);

    curprop = propn = GetChoice(prop_propn_item);
    ptype = GetChoice(rprop_choice_item);

    switch (ptype) {
    case 0:
	p[propn].type = ptype;
	p[propn].active = 1;
	p[propn].nodes = (int *) malloc(grid[curgrid].nmnp * sizeof(int));
	p[propn].ncolor = (int *) malloc(grid[curgrid].nmnp * sizeof(int));
	p[propn].dprop = (double *) malloc(grid[curgrid].nmnp * sizeof(double));
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    p[propn].dprop[i] = dp;
	    p[propn].ncolor[i] = curcolor;
	}
	break;
     case 1:
	p[propn].type = ptype;
	p[propn].active = 1;
	p[propn].stype = 1;
	p[propn].els = (int *) malloc(grid[curgrid].nmel * sizeof(int));
	p[propn].ecolor = (int *) malloc(grid[curgrid].nmel * sizeof(int));
	p[propn].dprop = (double *) malloc(grid[curgrid].nmel * sizeof(double));
        strcpy(buf,s);
        if ((pf = fopen(buf,"r")) == NULL) {
		errwin("Unable to open file");
		return;
	}
	for (i = 0; i < grid[curgrid].nmel; i++) {
            if (fgets(buf,250,pf) == NULL) {
		errwin("Not enough lines in depth file, operation cancelled.");
		fclose(pf);
		return;
	    }
            sscanf(buf,"%*d %le", &dp); 
	    if (i == 0) {
		min = max = dp;
	    }
            if (dp>max) {
		max=dp; 
	    }
	    if (dp<min) {
		min=dp;
	    }
	    p[propn].dprop[i] = dp;
	}
        fclose(pf);
        p[propn].propip.cmin=min;
	p[propn].propip.cmax=max;
        p[propn].propip.cint=1.0;
	p[propn].propip.nisol=16;
        break;
        }

    unset_wait_cursor(propread_frame);
}

void writep_done_proc(void)
{
    XtUnmanageChild(writep_frame);
}

void create_writep_popup(void)
{
    extern Widget app_shell;
    Widget bt, fr, rc, rc1;

    if (writep_frame) {
	XtRaise(writep_frame);
	return;
    }
    writep_frame = XmCreateDialogShell(app_shell, "Write properties", NULL, 0);
    writep_panel = XmCreateRowColumn(writep_frame, "rc1", NULL, 0);
    rc1 = XmCreateRowColumn(writep_panel, "rc1", NULL, 0);
    prop_write_item = CreateTextItem2(rc1, 15, "File:");

    rc = XtVaCreateManagedWidget("rc", xmRowColumnWidgetClass, rc1,
				 XmNorientation, XmHORIZONTAL,
				 NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) writep_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc,
				 NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) writep_done_proc, 0);
    XtManageChild(rc1);
    XtManageChild(writep_panel);
    XtManageChild(writep_frame);
}

void writep_accept_proc(void)
{
    int i, stype;
    char fname[256];
    strcpy(fname, (char *) xv_getstr(prop_write_item));
    if (!fexists(fname)) {
	propfp = fopen((char *) xv_getstr(prop_write_item), "w");
	if (propfp != NULL) {
	    stype = p[curprop].stype;
	    if (p[curprop].type) {
		/*fprintf(propfp, "Element property\n");*/
		/*fprintf(propfp, "%d %d %d %d\n", curprop, 1, grid[curgrid].nmel, stype);*/
		if (stype) {
		    for (i = 0; i < grid[curgrid].nmel; i++) {
			fprintf(propfp, "%d %lf\n", i + 1, p[curprop].dprop[i]);
		    }
		} else {
		    for (i = 0; i < grid[curgrid].nmel; i++) {
			fprintf(propfp, "%d %d\n", i + 1, (int) p[curprop].dprop[i]);
		    }
		}
	    } else {
		fprintf(propfp, "Nodal property\n");
		fprintf(propfp, "%d %d %d %d\n", curprop, 0, grid[curgrid].nmnp, stype);
		if (stype) {
		    for (i = 0; i < grid[curgrid].nmnp; i++) {
			fprintf(propfp, "%d %lf\n", i + 1, p[curprop].dprop[i]);
		    }
		} else {
		    for (i = 0; i < grid[curgrid].nmnp; i++) {
			fprintf(propfp, "%d %d\n", i + 1, (int) p[curprop].dprop[i]);
		    }
		}
	    }
	    fclose(propfp);
	}
    }
    XtUnmanageChild(writep_frame);
}

double prop_get(int gno, int i)
{
    double dp = 0.0, x, y;
    if (p[curprop].type) { /* element */
	dp = p[curprop].dprop[i];
    } else {
    }
    return dp;
}

void prop_set(int gno, int i)
{
    double dp, x, y;
    curcolor = GetChoice(prop_color_item);
    curval = atof((char *) xv_getstr(prop_current_item));
    if (p[curprop].type) {
	p[curprop].dprop[i] = curval;
	p[curprop].ecolor[i] = curcolor;
	setcolor(p[curprop].ecolor[i]);
	draw_element_filled(gno, i, 1.0);
	setcolor(1);
    } else {
	p[curprop].dprop[i] = curval;
	p[curprop].ncolor[i] = curcolor;
	setcolor(p[curprop].ncolor[i]);
	x = grid[gno].xord[i];
	y = grid[gno].yord[i];
	my_circle(x, y);
	setcolor(1);
    }
}

void prop_select_proc(void)
{
    do_select_region();
}

void prop_display_proc(void)
{
    draw_prop = XmToggleButtonGetState(display_prop);
}

void prop_pick_proc(void)
{
    if (p[curprop].type) {
	set_action(0);
	set_action(PICK_PROP_ELEM);
    } else {
	set_action(0);
	set_action(PICK_PROP_NODE);
    }
}

void prop_query_proc(void)
{
    if (p[curprop].type) {
	set_action(0);
	set_action(QUERY_PROP_ELEM);
    } else {
	set_action(0);
	set_action(QUERY_PROP_NODE);
    }
}

void prop_accept_proc(void)
{
    int i;
    double xg, yg;
    char buf[128];
    extern double regionx[], regiony[];
    extern int nregion;
    curval = atof((char *) xv_getstr(prop_current_item));
    curcolor = GetChoice(prop_color_item);
    if (p[curprop].type) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    get_center(0, i, &xg, &yg);
	    if (inregion(regionx, regiony, nregion, xg, yg)) {
		p[curprop].ecolor[i] = curcolor;
		p[curprop].dprop[i] = curval;
		my_move2(xg, yg);
		writetext(buf);
	    }
	}
    } else {
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    xg = grid[curgrid].xord[i];
	    yg = grid[curgrid].yord[i];
	    if (inregion(regionx, regiony, nregion, xg, yg)) {
		p[curprop].ncolor[i] = curcolor;
		p[curprop].dprop[i] = curval;
		my_move2(xg, yg);
		writetext(buf);
	    }
	}
    }
    do_clear_region();
}

void prop_initialize_proc(void)
{
    char buf[1024];
    int i, propn, ptype;
    float jf;
    double dp;
    FILE *tip;
    curprop = propn = GetChoice(prop_propn_item);
    if (p[propn].active) {
	if (yesno("Property active, clear?", "Press Yes or No", "Yes", "No")) {
	} else {
	    return;
	}
    }
    ptype = GetChoice(prop_type_item);
    dp = atof((char *) xv_getstr(prop_default_item));
    curcolor = GetChoice(prop_color_item);

    if (p[propn].els) {
	free(p[propn].els);
    }
    if (p[propn].nodes) {
	free(p[propn].nodes);
    }
    if (p[propn].ncolor) {
	free(p[propn].ncolor);
    }
    if (p[propn].dprop) {
	free(p[propn].dprop);
    }
    switch (ptype) {
    case 0:
	p[propn].type = ptype;
	p[propn].active = 1;
	p[propn].nodes = (int *) calloc(grid[curgrid].nmnp, sizeof(int));
	p[propn].ncolor = (int *) calloc(grid[curgrid].nmnp, sizeof(int));
	p[propn].dprop = (double *) calloc(grid[curgrid].nmnp, sizeof(double));
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    p[propn].dprop[i] = dp;
	    p[propn].ncolor[i] = curcolor;
	}
	break;
    case 1:
	p[propn].type = ptype;
	p[propn].active = 1;
	p[propn].els = (int *) calloc(grid[curgrid].nmel, sizeof(int));
	p[propn].ecolor = (int *) calloc(grid[curgrid].nmel, sizeof(int));
	p[propn].dprop = (double *) calloc(grid[curgrid].nmel, sizeof(double));
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    p[propn].dprop[i] = dp;
	    p[propn].ecolor[i] = curcolor;
	}
	break;
    }
}

void prop_kill_proc(void)
{
    int i, propn, ptype;
    double dp;
    curprop = propn = GetChoice(prop_propn_item);
    if (p[propn].active) {
	if (yesno("Property active, clear?", "Press Yes or No", "Yes", "No")) {
	} else {
	    return;
	}
    } else {
	errwin("Property already inactive");
	return;
    }
    p[propn].active = 0;
    if (p[propn].els) {
	free(p[propn].els);
    }
    if (p[propn].nodes) {
	free(p[propn].nodes);
    }
    if (p[propn].ncolor) {
	free(p[propn].ncolor);
    }
    if (p[propn].dprop) {
	free(p[propn].dprop);
    }
}

void update_prop(void)
{
    SetChoice(prop_color_item, curcolor);
    SetChoice(prop_propn_item, curprop);
}

void prop_done_proc(void)
{
    XtUnmanageChild(prop_frame);
}

/*
 * Main properties popup
 */
void create_prop_popup(void)
{
    extern Widget app_shell;
    Widget bt, rc;

    if (prop_frame) {
	update_prop();
	XtRaise(prop_frame);
	return;
    }
    prop_frame = XmCreateDialogShell(app_shell, "Properties", NULL, 0);
    prop_panel = XmCreateRowColumn(prop_frame, "rc1", NULL, 0);

    prop_propn_item = CreatePanelChoice2(prop_panel,
					 "Property #: ",
					 2, 11,
					 "1", "2", "3",
					 "4", "5", "6",
					 "7", "8", "9", "10",
					 NULL, 0);

    prop_type_item = CreatePanelChoice1(prop_panel,
					"Define: ",
					4,
					"Property at nodes",
					"Property at elements",
					"Element type (IMAT)",
					NULL, 0);

    prop_default_item = CreateTextItem2(prop_panel, 15, "Default property value: ");
    prop_current_item = CreateTextItem2(prop_panel, 15, "Current property value: ");
    prop_color_item = CreateColorChoice(prop_panel, "Color: ", 1);

    display_prop = XtVaCreateManagedWidget("Display current",
				   xmToggleButtonWidgetClass, prop_panel,
					   NULL);
    XtAddCallback(display_prop, XmNvalueChangedCallback, (XtCallbackProc) prop_display_proc, NULL);

    rc = XmCreateRowColumn(prop_panel, "rc1", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Initialize", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_initialize_proc, NULL);
    bt = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_pick_proc, NULL);
    bt = XtVaCreateManagedWidget("Region", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_select_proc, NULL);
    bt = XtVaCreateManagedWidget("Accept in region", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_accept_proc, NULL);
    bt = XtVaCreateManagedWidget("Query", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_query_proc, NULL);
    XtManageChild(rc);
    rc = XmCreateRowColumn(prop_panel, "rc1", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_isolines_proc, NULL);
    bt = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_readprop_frame, NULL);
    bt = XtVaCreateManagedWidget("Write...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) create_writep_popup, NULL);
    bt = XtVaCreateManagedWidget("Kill", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_kill_proc, NULL);
    bt = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) prop_done_proc, 0);
    XtManageChild(rc);
    update_prop();
    XtManageChild(prop_panel);
    XtManageChild(prop_frame);
}

/*
 * included for prop
 */
void drawgrid_prop_filled(int gridno)
{
    extern int mapisolbath[];
    char buf[1024];
    int i, j, k, n0, n1, n2,curcolor;
    double xg, yg, x[4], y[4],dp;
    for (j = 0; j < 10; j++) {
	if (p[j].active) {
/*PBP*/
	    if (p[j].type) {
                if (p[j].propip.active == 1){
                    for (i=0; i < grid[gridno].nmel; i++){
                        dp = p[j].dprop[i];
                        for (k=0;k<p[j].propip.nisol-1;k++){
                        if (dp>=p[j].propip.cis[k] && dp<p[j].propip.cis[k+1]){
                            p[j].ecolor[i] = k+16;k = p[j].propip.nisol+1;}
                        }
		        setcolor(p[j].ecolor[i]);
      			draw_element_filled(gridno, i, 1.0);
                    }
                    if (p[j].propip.lactive == ON) {
                        dolegend(p[j].propip, mapisolbath, p[j].propip.nisol);
		    }
                } else {
		    for (i = 0; i < grid[gridno].nmel; i++) {
		        setcolor(p[j].ecolor[i]);
    			draw_element_filled(gridno, i, 1.0);
	    	    }
                }
	    } else {
		for (i = 0; i < grid[gridno].nmnp; i++) {
		    setcolor(p[j].ncolor[i]);
		    x[0] = grid[gridno].xord[i];
		    y[0] = grid[gridno].yord[i];
		    my_circle(x[0], y[0]);
		}
	    }
	}
    }
    setcolor(1);
}

Widget cour_frame;
Widget cour_panel;
Widget cour_dt_item;
Widget cour1_item;
Widget cour2_item;
Widget cour3_item;
Widget display_cour_item;
Widget display_courval_item;

double cour_dt = 300.0, cour1_level = 0.8, cour2_level = 1.0;
int display_cour = 0, display_courn = 0;

Widget dimw_frame;
Widget dimw_panel;
Widget dimw_dt_item;
Widget dimw1_item;
Widget dimw2_item;
Widget dimw3_item;
Widget display_dimw_item;
Widget display_dimwval_item;

double dimw_dt = 4.1 * 3600.0, dimw1_level = 40.0, dimw2_level = 20.0;
int display_dimw = 0, display_dimwn = 0;

void dimw_accept_proc(void)
{
    dimw_dt = atof((char *) xv_getstr(dimw_dt_item)) * 3600.0;
    dimw1_level = atof((char *) xv_getstr(dimw1_item));
    dimw2_level = atof((char *) xv_getstr(dimw2_item));
    display_dimw = XmToggleButtonGetState(display_dimw_item);
    display_dimwn = XmToggleButtonGetState(display_dimwval_item);
    cour_dt = atof((char *) xv_getstr(cour_dt_item));
    cour1_level = atof((char *) xv_getstr(cour1_item));
    cour2_level = atof((char *) xv_getstr(cour2_item));
    display_cour = XmToggleButtonGetState(display_cour_item);
    display_courn = XmToggleButtonGetState(display_courval_item);
    do_drawgrid();
}

void update_dimw(void)
{
    char buf[128];
    if (dimw_frame) {
	sprintf(buf, "%.3lf", dimw_dt / 3600.0);
	panel_setstr_value(dimw_dt_item, buf);
	sprintf(buf, "%.3lf", dimw1_level);
	panel_setstr_value(dimw1_item, buf);
	sprintf(buf, "%.3lf", dimw2_level);
	panel_setstr_value(dimw2_item, buf);
	XmToggleButtonSetState(display_dimw_item, display_dimw, False);
	XmToggleButtonSetState(display_dimwval_item, display_dimwn, False);
	sprintf(buf, "%.3lf", cour_dt);
	panel_setstr_value(cour_dt_item, buf);
	sprintf(buf, "%.3lf", cour1_level);
	panel_setstr_value(cour1_item, buf);
	sprintf(buf, "%.3lf", cour2_level);
	panel_setstr_value(cour2_item, buf);
	XmToggleButtonSetState(display_cour_item, display_cour, False);
	XmToggleButtonSetState(display_courval_item, display_courn, False);
    }
/*
   SetChoice(dimw_color_item, curcolor);
   SetChoice(dimw_dimwn_item, curdimw);
 */
}

void dimw_done_proc(void)
{
    XtUnmanageChild(dimw_frame);
}

void create_dimw_popup(void)
{
    extern Widget app_shell;
    Widget bt, fr, rc, rc2, lab;

    if (dimw_frame) {
	update_dimw();
	XtRaise(dimw_frame);
	return;
    }
    dimw_frame = XmCreateDialogShell(app_shell, "Dimensionless numbers", NULL, 0);
    dimw_panel = XmCreateRowColumn(dimw_frame, "rc1", NULL, 0);
    XtVaSetValues(dimw_panel, XmNorientation, XmHORIZONTAL, NULL);

    fr = XmCreateFrame(dimw_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);
    lab = XmCreateLabel(rc, "Dimensionless wavelength", NULL, 0);
    XtManageChild(lab);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "T (in hours) = ", NULL, 0);
    dimw_dt_item = XmCreateText(rc2, "time", NULL, 0);
    XtVaSetValues(dimw_dt_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(dimw_dt_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Warning value", NULL, 0);
    dimw1_item = XmCreateText(rc2, "warn", NULL, 0);
    XtVaSetValues(dimw1_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(dimw1_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Unacceptable value", NULL, 0);
    dimw2_item = XmCreateText(rc2, "unaccept", NULL, 0);
    XtVaSetValues(dimw2_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(dimw2_item);
    XtManageChild(rc2);

    display_dimw_item = XtVaCreateManagedWidget("Display filled", xmToggleButtonWidgetClass, rc,
						NULL);
    display_dimwval_item = XtVaCreateManagedWidget("Display dimensionless wavelength", xmToggleButtonWidgetClass, rc,
						   NULL);


    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(dimw_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);
    lab = XmCreateLabel(rc, "Courant number", NULL, 0);
    XtManageChild(lab);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "dt (seconds) = ", NULL, 0);
    cour_dt_item = XmCreateText(rc2, "time", NULL, 0);
    XtVaSetValues(cour_dt_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(cour_dt_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Warning value", NULL, 0);
    cour1_item = XmCreateText(rc2, "warn", NULL, 0);
    XtVaSetValues(cour1_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(cour1_item);
    XtManageChild(rc2);

    rc2 = XmCreateRowColumn(rc, "rc2", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    lab = XmCreateLabel(rc2, "Unacceptable value", NULL, 0);
    cour2_item = XmCreateText(rc2, "unaccept", NULL, 0);
    XtVaSetValues(cour2_item, XmNcolumns, 10, NULL);
    XtManageChild(lab);
    XtManageChild(cour2_item);
    XtManageChild(rc2);

    display_cour_item = XtVaCreateManagedWidget("Display filled", xmToggleButtonWidgetClass, rc,
						NULL);
    display_courval_item = XtVaCreateManagedWidget("Display Courant number", xmToggleButtonWidgetClass, rc,
						   NULL);

    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(rc, "frame1", NULL, 0);
    rc2 = XmCreateRowColumn(fr, "rc2", NULL, 0);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dimw_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dimw_done_proc, 0);
    XtManageChild(rc2);
    XtManageChild(fr);

    update_dimw();
    XtManageChild(dimw_panel);
    XtManageChild(dimw_frame);
}

Widget grad_frame;
Widget grad_panel;
Widget grad_length_item;
Widget grad_length_label;
Widget grad_llength_item;
Widget grad_llength_label;
Widget *grad_units_item;
Widget grad_legend_item;
Widget display_grad_item;
Widget display_slope_item;
Widget slope_cutoff_item;

#define M 362
#define CM 317
#define MM 368
#define KM 350

int display_grad = 0, display_slope = 0, display_units = MM;
int display_grad_legend;
double slope_cutoff = 0.0;
double grad_legx, grad_legy;
double grad_units = 1000.0, grad_scale = 1.0;
double grad_legend_length = 2.0;

void grad_accept_proc(void)
{
/*
   int u = GetChoice(grad_units_item);
 */
    int u = 0;
    grad_scale = 1.0 / atof((char *) xv_getstr(grad_length_item));
    grad_legend_length = atof((char *) xv_getstr(grad_llength_item));
    slope_cutoff = atof((char *) xv_getstr(slope_cutoff_item));
    display_grad = XmToggleButtonGetState(display_grad_item);
    display_slope = XmToggleButtonGetState(display_slope_item);
    display_grad_legend = XmToggleButtonGetState(grad_legend_item);
    switch (u) {
    case 0:
	display_units = MM;
	grad_units = 1000.0;
	break;
    case 1:
	display_units = CM;
	grad_units = 10.0;
	break;
    case 2:
	display_units = M;
	grad_units = 1.0;
	break;
    case 3:
	display_units = KM;
	grad_units = 0.001;
	break;
    }
    do_drawgrid();
}

void update_grad(void)
{
    char buf[128];
    if (grad_frame) {
	sprintf(buf, "%.3lf", 1.0 / grad_scale * 10.0);
	panel_setstr_value(grad_length_item, buf);
	sprintf(buf, "%.3lf", grad_legend_length);
	panel_setstr_value(grad_llength_item, buf);
	sprintf(buf, "%.3lf", slope_cutoff);
	panel_setstr_value(slope_cutoff_item, buf);
	XmToggleButtonSetState(display_grad_item, display_grad, False);
	XmToggleButtonSetState(display_slope_item, display_slope, False);
/*
   switch (display_units) {
   case MM:
   SetChoice(grad_units_item, 0);
   break;
   case CM:
   SetChoice(grad_units_item, 1);
   break;
   case M:
   SetChoice(grad_units_item, 2);
   break;
   case KM:
   SetChoice(grad_units_item, 3);
   break;
   }
 */
    }
}

void grad_done_proc(void)
{
    XtUnmanageChild(grad_frame);
}

static void define_vplace_proc(void)
{
    set_action(0);
    set_action(PLACE_VSCALE);
}

void create_grad_popup(void)
{
    Widget bt, fr, rc, rc2;

    if (grad_frame) {
	update_grad();
	XtRaise(grad_frame);
	return;
    }
    grad_frame = XmCreateDialogShell(app_shell, "Gradients & Slopes", NULL, 0);
    handle_close(grad_frame);
    grad_panel = XmCreateRowColumn(grad_frame, "panel", NULL, 0);

    fr = XmCreateFrame(grad_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);

    display_grad_item = XtVaCreateManagedWidget("Display gradient",
					   xmToggleButtonWidgetClass, rc,
						NULL);
    display_slope_item = XtVaCreateManagedWidget("*Display slopes",
					   xmToggleButtonWidgetClass, rc,
						 NULL);
    grad_legend_item = XtVaCreateManagedWidget("Display gradient scale legend",
					   xmToggleButtonWidgetClass, rc,
					       NULL);

    grad_length_item = CreateTextItem2(rc, 10, "1 cm on screen =");
    grad_llength_item = CreateTextItem2(rc, 10, "Legend length (mm/m):");
    slope_cutoff_item = CreateTextItem2(rc, 10, "Slope cutoff =");

    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(grad_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) grad_accept_proc, 0);
    bt = XtVaCreateManagedWidget("Place", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) define_vplace_proc, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, grad_frame);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(grad_panel);
    XtManageChild(grad_frame);
    update_grad();
}

void drawgradlegend(int gno)
{
    double x, y;
    int units;
    double unitfac;
    char buf[50];
    units = display_units;
    unitfac = grad_units;
    switch (display_units) {
    case MM:
	sprintf(buf, "%.1lf mm/m", grad_legend_length);
	break;
    case CM:
	sprintf(buf, "%.1lf m/km", grad_legend_length);
	break;
    case M:
	sprintf(buf, "%.1lf m/s", grad_legend_length);
	break;
    case KM:
	sprintf(buf, "%.1lf km/s", grad_legend_length);
	break;
    }
    x = grad_legx;
    y = grad_legy;
    setcolor(1);
    drawvellegend(x, y, grad_legend_length, 0.0, grad_units * grad_scale, buf);
}

typedef struct _PNodeList {
    int type;
    int active;
    int n;
    int *nodes;
    int *ival;
    double *dval;
    int *color;
} PNodeList;

#define MAXNODELISTS 150

int curnl = 0;
int curnlcolor = 1;
int curnlival = 1;
double curnldval = 0.0;
int displaynl = 0;

PNodeList nl[MAXNODELISTS];

int DestroyPNodeList(int nlno);
int IsActivePNodeList(int nlno);
int CreatePNodeList(int nlno);
int DestroyPNodeList(int nlno);
void DestroyAllPNodeLists();
int GetNewPNodeList();

int IsActivePNodeList(int nlno)
{
    return nl[nlno].active;
}

int CreatePNodeList(int nlno)
{
    DestroyPNodeList(nlno);
    nl[nlno].active = 1;
    return 0;
}

int DestroyPNodeList(int nlno)
{
    if (nl[nlno].active && nl[nlno].n) {
	free(nl[nlno].nodes);
	free(nl[nlno].color);
	free(nl[nlno].ival);
    }
    nl[nlno].active = 0;
    nl[nlno].n = 0;
    nl[nlno].nodes = NULL;
    nl[nlno].color = NULL;
    nl[nlno].ival = NULL;
    return 0;
}

void DestroyAllPNodeLists()
{
    int i;
    for (i = 0; i < MAXNODELISTS; i++) {
	if (IsActivePNodeList(i)) {
	    DestroyPNodeList(i);
	}
    }
}

int GetNewPNodeList()
{
    int i;
    for (i = 0; i < MAXNODELISTS; i++) {
	if (!IsActivePNodeList(i)) {
	    return i;
	}
    }
    return -1;
}

int GetNumNodes(int nlno)
{
    return nl[nlno].n;
}

int *GetPNodeList(int nlno)
{
    return nl[nlno].nodes;
}

int GetPNodeListColor(int nlno, int pos)
{
    return nl[nlno].nodes[pos];
}

void AddToPNodeList(int nlno, int node, int ival, double dval, int color)
{
    if (GetNumNodes(nlno) == 0) {
	nl[nlno].active = 1;
	nl[nlno].nodes = (int *) malloc(sizeof(int));
	nl[nlno].color = (int *) malloc(sizeof(int));
	nl[nlno].ival = (int *) malloc(sizeof(int));
	nl[nlno].dval = (double *) malloc(sizeof(double));
    } else {
	nl[nlno].nodes = (int *) realloc(nl[nlno].nodes, (nl[nlno].n + 1) * sizeof(int));
	nl[nlno].color = (int *) realloc(nl[nlno].color, (nl[nlno].n + 1) * sizeof(int));
	nl[nlno].ival = (int *) realloc(nl[nlno].ival, (nl[nlno].n + 1) * sizeof(int));
	nl[nlno].dval = (double *) realloc(nl[nlno].dval, (nl[nlno].n + 1) * sizeof(double));
    }
    nl[nlno].nodes[nl[nlno].n] = node;
    nl[nlno].color[nl[nlno].n] = color;
    nl[nlno].ival[nl[nlno].n] = ival;
    nl[nlno].dval[nl[nlno].n] = dval;
    nl[nlno].n++;
}

void PickPNodeList()
{
    set_action(0);
    set_action(PICK_NODE_LIST);
}

int WritePNodeList(char *fname)
{
    int i, j, cnt = 0;
    FILE *fp;
    if (!fexists(fname)) {
	fp = fopen(fname, "w");
	if (fp != NULL) {
	    for (i = 0; i < MAXNODELISTS; i++) {
		if (nl[i].active) {
		    cnt++;
		}
	    }
	    fprintf(fp, "Node list file\n");
	    fprintf(fp, "%d Number of lists in this file\n", cnt);
	    for (i = 0; i < MAXNODELISTS; i++) {
		if (nl[i].active) {
		    fprintf(fp, "%d Number of nodes in this list\n", nl[i].n);
		    for (j = 0; j < nl[i].n; j++) {
			fprintf(fp, "%d %d %lf %d\n",
				nl[i].nodes[j] + 1,
				nl[i].ival[j],
				nl[i].dval[j],
				nl[i].color[j]);
		    }
		}
	    }
	    fclose(fp);
	} else {
	    return 0;
	}
    } else {
	return 0;
    }
    return 1;
}

int ReadPNodeList(char *fname)
{
    int i, j, cnt, nlno;
    int n, node, ival, color;
    double dval;
    char buf[256];
    FILE *fp = fopen(fname, "r");
    if (fp == NULL) {
	errwin("Unable to open file");
	return 1;
    }
    DestroyAllPNodeLists();
/*
   GetNewPNodeList();
   if (yesno("Replace current node lists?", "Press Yes or No", "Yes", "No")) {
   } else {
   }
 */
    fgets(buf, 255, fp);
    fgets(buf, 255, fp);
    sscanf(buf, "%d\n", &cnt);
    for (i = 0; i < cnt; i++) {
	fgets(buf, 255, fp);
	sscanf(buf, "%d", &n);
	for (j = 0; j < n; j++) {
	    fgets(buf, 255, fp);
	    sscanf(buf, "%d %d %lf %d\n", &node,
		   &ival,
		   &dval,
		   &color);
	    node--;
	    AddToPNodeList(0, node, ival, dval, color);
	}
    }
    return 0;
}

#include "symdefs.h"

void DrawPNodeLists(int gridno)
{
    int i, j, nn;
    double x, y;
    for (j = 0; j < MAXNODELISTS; j++) {
	if (nl[j].active) {
	    for (i = 0; i < nl[j].n; i++) {
		if (nl[j].nodes[i] >= 0 && nl[j].nodes[i] < grid[gridno].nmnp) {
		    nn = nl[j].nodes[i];
		    setcolor(nl[j].color[i]);
		    x = grid[gridno].xord[nn];
		    y = grid[gridno].yord[nn];
		    drawpolysym(&x, &y, 1, SYM_CIRCLE, 0, 0, 0.5 + nl[j].color[i] / 15.0);
		}
	    }
	}
    }
    setcolor(1);
}

void create_nodelist_frame(void);
void create_readnodelist_frame(void);
void create_writenodelist_frame(void);

void nodelist_read_proc(void);
void nodelist_write_proc(void);
static Widget rnodelist_frame;
static Widget wnodelist_frame;

void create_readnodelist_frame(void)
{
    Widget wbut, rc;
    int i;
    if (!rnodelist_frame) {
	rnodelist_frame = XmCreateFileSelectionDialog(app_shell, "rnodelist_frame", NULL, 0);
	XtVaSetValues(rnodelist_frame,
		      XmNtitle, XmStringCreate("Read node list", charset),
		      XmNdirMask, XmStringCreate("*.nlist", charset),
		      NULL);

	rc = XmCreateRowColumn(rnodelist_frame, "rc", NULL, 0);
	XtManageChild(rc);
	XtAddCallback(rnodelist_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, rnodelist_frame);
	XtAddCallback(rnodelist_frame, XmNokCallback, (XtCallbackProc) nodelist_read_proc, NULL);
    }
    XtRaise(rnodelist_frame);
}

void nodelist_read_proc(void)
{
    Arg args;
    XmString list_item;
    char *s, buf[256];

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(rnodelist_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    set_wait_cursor(rnodelist_frame);
    if (ReadPNodeList(s)) {
	sprintf(buf, "Error reading file %s", s);
	errwin(buf);
    } else {
	XtUnmanageChild(rnodelist_frame);
    }
    unset_wait_cursor(rnodelist_frame);
}

void create_writenodelist_frame(void)
{
    extern Widget app_shell;
    Widget bt, fr, rc, rc1;
    Widget wbut;
    int i;
    if (!wnodelist_frame) {
	wnodelist_frame = XmCreateFileSelectionDialog(app_shell, "wnodelist_frame", NULL, 0);
	XtVaSetValues(wnodelist_frame,
		    XmNtitle, XmStringCreate("Write node list", charset),
		      XmNdirMask, XmStringCreate("*.nlist", charset),
		      NULL);

	rc = XmCreateRowColumn(wnodelist_frame, "rc", NULL, 0);
	XtManageChild(rc);
	XtAddCallback(wnodelist_frame, XmNcancelCallback, (XtCallbackProc) destroy_dialog, wnodelist_frame);
	XtAddCallback(wnodelist_frame, XmNokCallback, (XtCallbackProc) nodelist_write_proc, NULL);
    }
    XtRaise(wnodelist_frame);
}

void nodelist_write_proc(void)
{
    Arg args;
    XmString list_item;
    char *s;

    XtSetArg(args, XmNtextString, &list_item);
    XtGetValues(wnodelist_frame, &args, 1);
    XmStringGetLtoR(list_item, charset, &s);
    set_wait_cursor(wnodelist_frame);
    if (WritePNodeList(s)) {
    } else {
	XtUnmanageChild(wnodelist_frame);
    }
    unset_wait_cursor(wnodelist_frame);
}

static Widget nodelist_frame;
static Widget nodelist_panel;
static Widget display_nodelist_item;
static Widget nodelist_ival_item;
static Widget nodelist_dval_item;
static Widget *nodelist_color_item;
static Widget nlbut[10];

void nodelist_accept_proc(void)
{
    int i;
    curnlival = atoi((char *) xv_getstr(nodelist_ival_item));
    curnldval = atof((char *) xv_getstr(nodelist_dval_item));
    curnlcolor = GetChoice(nodelist_color_item);
    displaynl = XmToggleButtonGetState(display_nodelist_item);
/*
   for (i = 0; i < 4; i++) {
   XtSetSensitive(nlbut[i], True);
   }
 */
}

void update_nodelist(void)
{
    char buf[128];
    if (nodelist_frame) {
	sprintf(buf, "%d", curnlival);
	xv_setstr(nodelist_ival_item, buf);
	sprintf(buf, "%lf", curnldval);
	xv_setstr(nodelist_dval_item, buf);
	SetChoice(nodelist_color_item, curnlcolor);
	XmToggleButtonSetState(display_nodelist_item, displaynl, False);
    }
}

void nodelist_done_proc(void)
{
    XtUnmanageChild(nodelist_frame);
}

void create_nodelist_popup(void)
{
    Widget bt, fr, rc, rc2;
    int i;

    if (nodelist_frame) {
	update_nodelist();
	XtRaise(nodelist_frame);
	return;
    }
    CreatePNodeList(0);
    nodelist_frame = XmCreateDialogShell(app_shell, "Select nodes", NULL, 0);
    handle_close(nodelist_frame);
    nodelist_panel = XmCreateRowColumn(nodelist_frame, "panel", NULL, 0);

    fr = XmCreateFrame(nodelist_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);

    nodelist_ival_item = CreateTextItem2(rc, 15, "Current node ival: ");
    nodelist_dval_item = CreateTextItem2(rc, 15, "Current node dval: ");
    nodelist_color_item = CreateColorChoice(rc, "Color: ", 1);
    display_nodelist_item = XtVaCreateManagedWidget("Display selected nodes",
					   xmToggleButtonWidgetClass, rc,
						    NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(nodelist_panel, "frame1", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc1", NULL, 0);
    XtVaSetValues(rc,
		  XmNorientation, XmHORIZONTAL,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) nodelist_accept_proc, NULL);
    nlbut[0] = XtVaCreateManagedWidget("Pick", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(nlbut[0], XmNactivateCallback, (XtCallbackProc) PickPNodeList, NULL);
    nlbut[1] = XtVaCreateManagedWidget("Clear All", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(nlbut[1], XmNactivateCallback, (XtCallbackProc) DestroyAllPNodeLists, NULL);
    nlbut[2] = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(nlbut[2], XmNactivateCallback, (XtCallbackProc) create_readnodelist_frame, NULL);
    nlbut[3] = XtVaCreateManagedWidget("Write...", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(nlbut[3], XmNactivateCallback, (XtCallbackProc) create_writenodelist_frame, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, nodelist_frame);
    XtManageChild(rc);
    XtManageChild(fr);
    XtManageChild(nodelist_panel);
    XtManageChild(nodelist_frame);
/*
   for (i = 0; i < 4; i++) {
   XtSetSensitive(nlbut[i], False);
   }
 */
    update_nodelist();
}

typedef struct {
    Widget top;
    Widget display;
    Widget *item;
    Isolparms delhip;
    Isolparms delhhip;
    Isolparms hip;
} DelHUI;

static DelHUI dui;
int display_delh = 0;

void delh_accept_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    DelHUI *dui = (DelHUI *) clientd;
    display_delh = XmToggleButtonGetState(dui->display);
}

void drawdelh(int gridno)
{
    extern int mapisolconc[];
    int cd = GetChoice(dui.item);
    int i, n1, n2, n3;
    double hmin, hmax, delh, h;
    switch ((int) cd) {
    case 0:
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    hmax = hmin = grid[gridno].depth[n1];
	    if (hmin > grid[gridno].depth[n2]) {
		hmin = grid[gridno].depth[n2];
	    }
	    if (hmin > grid[gridno].depth[n3]) {
		hmin = grid[gridno].depth[n3];
	    }
	    if (hmax < grid[gridno].depth[n2]) {
		hmax = grid[gridno].depth[n2];
	    }
	    if (hmax < grid[gridno].depth[n3]) {
		hmax = grid[gridno].depth[n3];
	    }
	    delh = hmax - hmin;
/*PBPFOUND*/
	    fillel(i, delh, dui.delhip, mapisolconc, dui.delhip.nisol);
	}
	if (dui.delhip.lactive == ON) {
	    dolegend(dui.delhip, mapisolconc, dui.delhip.nisol);
	}
	break;
    case 1:
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    hmax = hmin = grid[gridno].depth[n1];
	    if (hmin > grid[gridno].depth[n2]) {
		hmin = grid[gridno].depth[n2];
	    }
	    if (hmin > grid[gridno].depth[n3]) {
		hmin = grid[gridno].depth[n3];
	    }
	    if (hmax < grid[gridno].depth[n2]) {
		hmax = grid[gridno].depth[n2];
	    }
	    if (hmax < grid[gridno].depth[n3]) {
		hmax = grid[gridno].depth[n3];
	    }
	    delh = hmax - hmin;
	    h = (grid[gridno].depth[n1] + grid[gridno].depth[n2] + grid[gridno].depth[n3]);
	    h *= 0.333333333333;
	    if (h < 1.0) {
		h = 1.0;
	    }
	    fillel(i, delh / h, dui.delhhip, mapisolconc, dui.delhhip.nisol);
	}
	if (dui.delhhip.lactive == ON) {
	    dolegend(dui.delhhip, mapisolconc, dui.delhhip.nisol);
	}
	break;
    case 2:
	for (i = 0; i < grid[gridno].nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    h = (grid[gridno].depth[n1] + grid[gridno].depth[n2] + grid[gridno].depth[n3]);
	    h *= 0.333333333333;
	    fillel(i, h, dui.hip, mapisolconc, dui.hip.nisol);
	}
	if (dui.hip.lactive == ON) {
	    dolegend(dui.hip, mapisolconc, dui.hip.nisol);
	}
	break;
    }
}

void set_delh_limits(int gridno)
{
    int i, n1, n2, n3;
    double h, delh, delhh;
    double tmin, tmax;
    double hmin, hmax;
    double delhmin, delhmax;
    double delhhmin, delhhmax;
    delhmin = delhhmin = 1e38;
    delhmax = delhhmax = -1e38;
    for (i = 0; i < grid[gridno].nmel; i++) {
	n1 = grid[gridno].icon[i].nl[0];
	n2 = grid[gridno].icon[i].nl[1];
	n3 = grid[gridno].icon[i].nl[2];
	tmax = tmin = grid[gridno].depth[n1];
	if (tmin > grid[gridno].depth[n2]) {
	    tmin = grid[gridno].depth[n2];
	}
	if (tmin > grid[gridno].depth[n3]) {
	    tmin = grid[gridno].depth[n3];
	}
	if (tmax < grid[gridno].depth[n2]) {
	    tmax = grid[gridno].depth[n2];
	}
	if (tmax < grid[gridno].depth[n3]) {
	    tmax = grid[gridno].depth[n3];
	}
	delh = tmax - tmin;
	if (delhmin > delh) {
	    delhmin = delh;
	}
	if (delhmax < delh) {
	    delhmax = delh;
	}
	h = (grid[gridno].depth[n1] + grid[gridno].depth[n2] + grid[gridno].depth[n3]);
	h *= 0.333333333333;
	if (h < 1.0) {
	    h = 1.0;
	}
	if (hmin > h) {
	    hmin = h;
	}
	if (hmax < h) {
	    hmax = h;
	}
	delhh = delh / h;
	if (delhhmin > delhh) {
	    delhhmin = delhh;
	}
	if (delhhmax < delhh) {
	    delhhmax = delhh;
	}
    }
    dui.delhip.cmin = delhmin;
    dui.delhip.cmax = delhmax;
    dui.delhhip.cmin = delhhmin;
    dui.delhhip.cmax = delhhmax;
    dui.hip.cmin = hmin;
    dui.hip.cmax = hmax;
}

void delh_isol_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    DelHUI *dui = (DelHUI *) clientd;
    int cd = GetChoice(dui->item);
    switch ((int) cd) {
    case 0:
	create_isolines_popup("Delta H", &(dui->delhip));
	break;
    case 1:
	create_isolines_popup("Delta H / H", &(dui->delhhip));
	break;
    case 2:
	create_isolines_popup("H", &(dui->hip));
	break;
    }
}

void update_delh(DelHUI dui)
{
    if (dui.top) {
    }
}

void create_delh_popup(void)
{
    Widget panel, bt, fr, rc, rc2, lab;
    if (dui.top) {
	set_delh_limits(curgrid);
	update_delh(dui);
	XtRaise(dui.top);
	return;
    }
    dui.top = XmCreateDialogShell(app_shell, "DelH/H", NULL, 0);
    panel = XmCreateRowColumn(dui.top, "panel", NULL, 0);
    dui.display = XtVaCreateManagedWidget("Display selected quantity",
				 xmToggleButtonWidgetClass, panel, NULL);
    dui.item = CreatePanelChoice1(panel, "Show: ", 4,
				  "Delta H",
				  "Delta H / H",
				  "H", NULL, 0);
    rc2 = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) delh_accept_proc, (XtPointer) & dui);
    bt = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) delh_isol_proc, (XtPointer) & dui);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) dui.top);
    XtManageChild(rc2);

    XtManageChild(panel);
    XtManageChild(dui.top);
    set_isolines_defaults(&dui.delhip);
    set_isolines_defaults(&dui.delhhip);
    set_isolines_defaults(&dui.hip);
    set_delh_limits(curgrid);
    update_delh(dui);
}

typedef struct {
    Widget top;
    Widget display;
    Widget *item;
    Isolparms ip;
} ElDepthsUI;

void update_dpe(ElDepthsUI dui)
{
    if (dui.top) {
    }
}

void dpe_isol_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    ElDepthsUI *dui = (ElDepthsUI *) clientd;
    set_isolines_defaults(dui->ip);
    create_isolines_popup("Element depths", &(dui->ip));
}

void dpe_read_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    ElDepthsUI *dui = (ElDepthsUI *) clientd;
}

void dpe_write_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    ElDepthsUI *dui = (ElDepthsUI *) clientd;
}

void dpe_accept_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    ElDepthsUI *dui = (ElDepthsUI *) clientd;
}

void dpe_query_proc(Widget w, XtPointer clientd, XtPointer calld)
{
    ElDepthsUI *dui = (ElDepthsUI *) clientd;
}

/*
 * Depths at elements
 */
void create_dpe_frame(void)
{
    static ElDepthsUI ui;
    Widget panel, bt, fr, rc, rc2, lab;
    if (ui.top) {
	update_dpe(ui);
	XtRaise(ui.top);
	return;
    }
    ui.top = XmCreateDialogShell(app_shell, "Depths at Elements", NULL, 0);
    panel = XmCreateRowColumn(ui.top, "panel", NULL, 0);
    ui.display = XtVaCreateManagedWidget("Display selected quantity",
				 xmToggleButtonWidgetClass, panel, NULL);
    ui.item = CreatePanelChoice1(panel, "Show: ", 3,
				  "Isolines",
				  "Values",
				  NULL, 0);
    rc2 = XmCreateRowColumn(panel, "rc", NULL, 0);
    XtVaSetValues(rc2,
		  XmNorientation, XmHORIZONTAL,
		  XmNpacking, XmPACK_TIGHT,
		  NULL);
    bt = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dpe_accept_proc, (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Read...", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dpe_read_proc, (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Write...", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dpe_write_proc, (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Query...", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dpe_write_proc, (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Isolines...", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) dpe_isol_proc, (XtPointer) & ui);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, (XtPointer) ui.top);
    XtManageChild(rc2);

    XtManageChild(panel);
    XtManageChild(ui.top);
    set_isolines_defaults(&ui.ip);
    update_dpe(ui);
}
