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
 * edit popup
 *
 */

#ifndef lint
static char RCSid[] = "$Id: editwin.c,v 1.11 2008/01/17 16:27:48 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Xm/Text.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"

extern Widget app_shell;
extern XmStringCharSet charset;

static Widget edit_frame;
static Widget edit_panel;
static Widget equiv_radius_item;
static Widget *op_select_item;
static Widget *grid_select_item;
static Widget *move_select1_item;
static Widget *move_select2_item;
static Widget *join_select1_item;
static Widget *join_select2_item;
static Widget *delete_select_item;
static Widget *current_select_item;
static Widget *extract_select1_item;
static Widget *extract_select2_item;
static Widget *extract_method_item;
static Widget rotate_degrees_item;
static Widget rotate_aboutx_item;
static Widget rotate_abouty_item;
static Widget translate_x_item;
static Widget translate_y_item;
static Widget evaluate_item;
static Widget *eval_region_item;
static Widget *split_item;
static Widget goto_node_item;
static Widget goto_elem_item;
static Widget goto_pointx_item;
static Widget goto_pointy_item;
static Widget copy_from_item;
static Widget copy_to_item;
static Widget *delete_toggle_item;

void create_shapiro_frame(void);
void create_ftri_frame(void);
void create_eval_frame(void);
void create_rot_frame(void);
void create_trans_frame(void);
void create_equiv_frame(void);
void create_delr_frame(void);
void create_splitr_frame(void);
void create_goto_popup(void);
void create_copy_popup(void);
void create_extract_frame(void);
void create_cpp_frame(void);
void create_volume_frame();

void create_tri2quad_frame(void);
void create_quad2tri_frame(void);
void create_qqual_frame(void);
void create_highlight_frame(void);

void do_select_region(void);
void do_clear_region(void);
void do_evaluate_function(void);
void do_rotate_grid(void);
void do_translate_grid(void);

void do_copy_grid(void);
void do_move_grid(void);
void do_delete_grid(void);
void do_join_grid(void);
void do_acceptjoin_grid(void);
void do_extract_grid(void);

void do_make_current(void);

void do_mark_proc(void);
void do_marked_delete_proc(void);

void do_region_delete_proc(void);
void do_region_split3_proc(void);
void do_region_split4_proc(void);
void do_split3_proc(void);
void do_split4_proc(void);
void do_marked_split3_proc(void);
void do_marked_split4_proc(void);

void do_fix_areas_proc(void);

void set_editwin(void);

/*
 * Event and Notify proc declarations
 */

void not_active(void)
{
    errwin("Operation not implemented");
}

/*
 * Create the region Frame and Panel
 */
void create_edit_popup(void)
{
    Widget wbut;
    Widget fr, rc;
    Arg a[2];
    int ac = 0;

    if (edit_frame) {
	XtRaise(edit_frame);
	return;
    }
    edit_frame = XmCreateDialogShell(app_shell, "Edit region", NULL, 0);

    edit_panel = XmCreateRowColumn(edit_frame, "rc", NULL, 0);

    wbut = XtVaCreateManagedWidget("Evaluate...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_eval_frame, 0);

    wbut = XtVaCreateManagedWidget("Translate and rotate...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_trans_frame, 0);

    wbut = XtVaCreateManagedWidget("Acceptable skewness...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_equiv_frame, 0);

    wbut = XtVaCreateManagedWidget("Mark Shapiro filter violations...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_shapiro_frame, 0);

    wbut = XtVaCreateManagedWidget("Delete elements in region...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_delr_frame, 0);
    wbut = XtVaCreateManagedWidget("Filter triangles...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_ftri_frame, 0);

    wbut = XtVaCreateManagedWidget("Split elements...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_splitr_frame, 0);

    wbut = XtVaCreateManagedWidget("Copy grids...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_copy_popup, 0);

    wbut = XtVaCreateManagedWidget("Quality check for quadrangles...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_qqual_frame, 0);

    wbut = XtVaCreateManagedWidget("Convert Tris to Quads...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_tri2quad_frame, 0);
    wbut = XtVaCreateManagedWidget("Convert Quads to Tris...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_quad2tri_frame, 0);
    wbut = XtVaCreateManagedWidget("Highlight Quads/Tris...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_highlight_frame, 0);

    wbut = XtVaCreateManagedWidget("Extract sub-grid...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_extract_frame, 0);
    wbut = XtVaCreateManagedWidget("Fix negative areas", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_fix_areas_proc, NULL);
    wbut = XtVaCreateManagedWidget("Volumes...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_volume_frame, NULL);
    wbut = XtVaCreateManagedWidget("CPP projection...", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) create_cpp_frame, NULL);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, edit_panel, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, edit_frame);
    XtManageChild(edit_panel);
    XtManageChild(edit_frame);
}

void do_mark_proc(void)
{
    double dist, xg, yg, cutoff, er, tmp = 1e-300, a, area();
    int i, n1, n2, n3;
    char s[128];

    strcpy(s, (char *) panel_getstr_value(equiv_radius_item));
    cutoff = atof(s);
    if (cutoff <= 0.0) {
	errwin("Cutoff <= 0.0, operation cancelled");
	return;
    }
    for (i = 0; i < grid[curgrid].nmel; i++) {
	if (grid[curgrid].icon[i].type == 3) {
	    a = area(curgrid, i);
	    er = sqrt(a / M_PI);
	    n1 = grid[curgrid].icon[i].nl[0];
	    n2 = grid[curgrid].icon[i].nl[1];
	    n3 = grid[curgrid].icon[i].nl[2];
	    tmp = hypot(grid[curgrid].xord[n1] - grid[curgrid].xord[n2], grid[curgrid].yord[n1] - grid[curgrid].yord[n2]);
	    dist = hypot(grid[curgrid].xord[n2] - grid[curgrid].xord[n3], grid[curgrid].yord[n2] - grid[curgrid].yord[n3]);
	    if (dist > tmp) {
		tmp = dist;
	    }
	    dist = hypot(grid[curgrid].xord[n1] - grid[curgrid].xord[n3], grid[curgrid].yord[n1] - grid[curgrid].yord[n3]);
	    if (dist > tmp) {
		tmp = dist;
	    }
	    if (er > 0.0) {
		if (tmp / er > cutoff) {
		    get_center(curgrid, i, &xg, &yg);
		    solidbox(xg, yg);
		    grid[curgrid].ellist[i] = 1;
		} else {
		    grid[curgrid].ellist[i] = 0;
		}
	    }
	}
    }
}

void do_marked_delete_proc(void)
{
    compact_grid(curgrid);
}

void do_fix_areas_proc(void)
{
    fix_areas(curgrid);
}

void set_editwin(void)
{
}

void do_select_region(void)
{
    if (region_flag) {
	if (yesno("Region defined, clear?", "Press Yes or No", "Yes", "No")) {
	    do_clear_region();
	} else {
	    errwin("Region not cleared");
	    return;
	}
    }
    set_action(0);
    set_action(DEFINE_REGION);
}

void do_clear_region(void)
{
    setcolor(0);
    draw_region(regionx, regiony, nregion);
    setcolor(1);
    nregion = 0;
    region_flag = 0;
}

/*
 * write elements inside (inout==1) a region, outside if inout==0
 * inside region means center of element is in the region.
 * TODO allow for varying definitions of inside.
 *       1. center
 *       2. all nodes
 *       3. one node
 */
int do_writeelems_region(char *fname, int inout)
{
    int i;
    double xg, yg;
    FILE *fp;
    if ((fp = fopen(fname, "wb")) == NULL) {
	errwin("Unable to open file.");
	return 1;
    }
    for (i = 0; i < grid[0].nmel; i++) {
        get_center(0, i, &xg, &yg);
	if (inout) {
            if (region_flag && inregion(regionx, regiony, nregion, xg, yg)) {
		fprintf(fp, "%d\n", i + 1);
	    }
	} else {
            if (region_flag && !inregion(regionx, regiony, nregion, xg, yg)) {
		fprintf(fp, "%d\n", i + 1);
            }
	}
    }
    fclose(fp);
    return 0;
}

/*
 * write nodes inside (inout==1) a region, outside if inout==0
 */
int do_writenodes_region(char *fname, int inout)
{
    int i;
    char buf[256];
    double x, y;
    FILE *fp;
    if ((fp = fopen(fname, "wb")) == NULL) {
	errwin("Unable to open file.");
	return 1;
    }
    for (i = 0; i < grid[0].nmnp; i++) {
        x = grid[0].xord[i];
        y = grid[0].yord[i];
	if (inout) {
            if (region_flag && inregion(regionx, regiony, nregion, x, y)) {
		fprintf(fp, "%d\n", i + 1);
            }
	} else {
            if (region_flag && !inregion(regionx, regiony, nregion, x, y)) {
		fprintf(fp, "%d\n", i + 1);
            }
	}
    }
    fclose(fp);
    return 0;
}

int inregion(double *boundx, double *boundy, int nbpts, double x, double y)
{
    return (inbound(x, y, boundx, boundy, nbpts));
}

void do_evaluate_function(void)
{
    int i, j, itmp, inode, ib, ind, gridno = curgrid, errpos;
    double xtmp, ytmp, dtmp, b = 0.0, c = 0.0, d = 0.0;
    char s[128];

    strcpy(s, (char *) panel_getstr_value(evaluate_item));
    itmp = (int) GetChoice(eval_region_item);
    fixupstr(s);
    switch (itmp) {
    case 0:
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    xtmp = grid[gridno].xord[i];
	    ytmp = grid[gridno].yord[i];
	    dtmp = grid[gridno].depth[i];
	    scanner(s, &xtmp, &ytmp, &dtmp, &b, &c, &d, i + 1, gridno, &errpos);
	    grid[gridno].xord[i] = xtmp;
	    grid[gridno].yord[i] = ytmp;
	    grid[gridno].depth[i] = dtmp;
	}
	break;
    case 1:
	if (region_flag) {
	    for (i = 0; i < grid[curgrid].nmnp; i++) {
		xtmp = grid[gridno].xord[i];
		ytmp = grid[gridno].yord[i];
		dtmp = grid[gridno].depth[i];
		if (inregion(regionx, regiony, nregion, xtmp, ytmp)) {
		    scanner(s, &xtmp, &ytmp, &dtmp, &b, &c, &d, i + 1, gridno, &errpos);
		    grid[gridno].xord[i] = xtmp;
		    grid[gridno].yord[i] = ytmp;
		    grid[gridno].depth[i] = dtmp;
		}
	    }
	} else {
	    errwin("Region not defined");
	}
	break;
    case 2:
	for (i = 0; i < grid[curgrid].nbounds; i++) {
	    ib = grid[curgrid].boundaries[i];
	    for (j = 0; j < boundary[ib].nbpts; j++) {
		inode = boundary[ib].nodes[j];
		xtmp = grid[gridno].xord[inode];
		ytmp = grid[gridno].yord[inode];
		dtmp = grid[gridno].depth[inode];
		scanner(s, &xtmp, &ytmp, &dtmp, &b, &c, &d, i + 1, gridno, &errpos);
		grid[gridno].xord[inode] = xtmp;
		grid[gridno].yord[inode] = ytmp;
		grid[gridno].depth[inode] = dtmp;
	    }
	}
	break;
    }
    set_grid_limits(gridno);
    do_drawgrid();
}

void do_rotate_grid(void)
{
    int i;
    double rotx, roty, degrees, xtmp, ytmp;
    char buf[128];
    strcpy(buf, (char *) panel_getstr_value(rotate_degrees_item));
    degrees = atof(buf);
    if (degrees == 0.0) {
	return;
    }
    degrees = M_PI / 180.0 * degrees;
    strcpy(buf, (char *) panel_getstr_value(rotate_aboutx_item));
    rotx = atof(buf);
    strcpy(buf, (char *) panel_getstr_value(rotate_abouty_item));
    roty = atof(buf);
    for (i = 0; i < grid[curgrid].nmnp; i++) {
	xtmp = grid[curgrid].xord[i] - rotx;
	ytmp = grid[curgrid].yord[i] - roty;
	grid[curgrid].xord[i] = rotx + cos(degrees) * xtmp - sin(degrees) * ytmp;
	grid[curgrid].yord[i] = roty + sin(degrees) * xtmp + cos(degrees) * ytmp;
    }
    set_grid_limits(curgrid);
}

void do_translate_grid(void)
{
    int i;
    double tx, ty;
    char buf[128];
    strcpy(buf, (char *) panel_getstr_value(translate_x_item));
    tx = atof(buf);
    strcpy(buf, (char *) panel_getstr_value(translate_y_item));
    ty = atof(buf);
    for (i = 0; i < grid[curgrid].nmnp; i++) {
	grid[curgrid].xord[i] -= tx;
	grid[curgrid].yord[i] -= ty;
    }
    do_rotate_grid();
    set_grid_limits(curgrid);
    do_drawgrid();
}

static int copyfrom, copyto;

void do_copy_grid(void)
{
    int from, to;
    copy_grid(copyfrom, copyto);
}

void set_copyfrom_proc(int w, int cd)
{
    copyfrom = cd;
}

void set_copyto_proc(int w, int cd)
{
    copyto = cd;
}

void create_copy_popup(void)
{
    static Widget top, dialog;
    Widget wbut, rc;
    Widget bt, fr2, rc2;
    int x, y;
    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Copy grid", NULL, 0);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    handle_close(top);

    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    copy_from_item = XmVaCreateSimpleOptionMenu(dialog, "option_menu",
						XmStringCreateSimple("Copy grid: "), 'a', 0, (XtCallbackProc) set_copyfrom_proc,
						XmVaPUSHBUTTON, XmStringCreateSimple("0"), 0, NULL, NULL,
						XmVaPUSHBUTTON, XmStringCreateSimple("1"), 0, NULL, NULL,
						XmVaPUSHBUTTON, XmStringCreateSimple("2"), 0, NULL, NULL,
						XmVaPUSHBUTTON, XmStringCreateSimple("3"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("4"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("Background"), 0, NULL, NULL, NULL);
    XtManageChild(copy_from_item);
    copy_to_item = XmVaCreateSimpleOptionMenu(dialog, "option_menu",
					      XmStringCreateSimple("To grid: "), 'a', 0, (XtCallbackProc) set_copyto_proc,
					      XmVaPUSHBUTTON, XmStringCreateSimple("0"), 0, NULL, NULL,
					      XmVaPUSHBUTTON, XmStringCreateSimple("1"), 0, NULL, NULL,
					      XmVaPUSHBUTTON, XmStringCreateSimple("2"), 0, NULL, NULL,
					      XmVaPUSHBUTTON, XmStringCreateSimple("3"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("4"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("Background"), 0, NULL, NULL, NULL);
    XtManageChild(copy_to_item);

    fr2 = XmCreateFrame(dialog, "fr", NULL, 0);
    rc2 = XmCreateRowColumn(fr2, "rc", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Apply", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_copy_grid, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc2);
    XtManageChild(fr2);

    XtManageChild(dialog);
    XtManageChild(top);
}

void do_move_grid(void)
{
}

void do_delete_grid(void)
{
}

void do_join_grid(void)
{
}

void do_acceptjoin_grid(void)
{
}

void do_make_current(void)
{
}

void do_region_delete_proc(void)
{
    double xg, yg;
    int i;
    int which = GetChoice(delete_toggle_item);
    if (region_flag) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    get_center(curgrid, i, &xg, &yg);
	    if (!which) {	/* delete in region */
		if (inregion(regionx, regiony, nregion, xg, yg)) {
		    solidbox(xg, yg);
		    grid[curgrid].ellist[i] = 1;
		} else {
		    grid[curgrid].ellist[i] = 0;
		}
	    } else {		/* delete outside region */
		if (!inregion(regionx, regiony, nregion, xg, yg)) {
		    solidbox(xg, yg);
		    grid[curgrid].ellist[i] = 1;
		} else {
		    grid[curgrid].ellist[i] = 0;
		}
	    }
	}
    } else {
	errwin("Region not defined");
    }
}

static int doall = 0;

void do_marked_split3_proc(void)
{
    int i, itmp = grid[curgrid].nmel;
    if (region_flag || doall) {
	for (i = 0; i < itmp; i++) {
	    if (grid[curgrid].ellist[i] && (grid[curgrid].icon[i].type == 3)) {
		split_elem3(curgrid, i);
	    }
	}
    } else {
	errwin("Region not defined");
    }
    doall = 0;
}

void do_marked_split4_proc(void)
{
    int i, itmp = grid[curgrid].nmel;
    if (region_flag || doall) {
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    grid[curgrid].nlist[i] = 1;
	}
	split_allelem4(curgrid);
    } else {
	errwin("Region not defined");
    }
    doall = 0;
}

void do_split3_proc(void)
{
    int i;
    for (i = 0; i < grid[curgrid].nmel; i++) {
	if (grid[curgrid].icon[i].type == 3) {
	    grid[curgrid].ellist[i] = 1;
	} else {
	    grid[curgrid].ellist[i] = 0;
	}
    }
    doall = 1;
}

void do_split4_proc(void)
{
    int i;
    for (i = 0; i < grid[curgrid].nmel; i++) {
	if (grid[curgrid].icon[i].type == 3) {
	    grid[curgrid].ellist[i] = 1;
	} else {
	    grid[curgrid].ellist[i] = 0;
	}
    }
    doall = 1;
}

void do_region_split3_proc(void)
{
    double xg, yg;
    int i;
    if (region_flag) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    grid[curgrid].ellist[i] = 0;
	    if (grid[curgrid].icon[i].type == 3) {
		get_center(curgrid, i, &xg, &yg);
		if (inregion(regionx, regiony, nregion, xg, yg)) {
		    writestr(xg, yg, 0, 0, "3");
		    grid[curgrid].ellist[i] = 1;
		}
	    }
	}
    } else {
	errwin("Region not defined");
    }
}

void do_region_split4_proc(void)
{
    double xg, yg;
    int i;
    if (region_flag) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    grid[curgrid].ellist[i] = 0;
	    get_center(curgrid, i, &xg, &yg);
	    if (grid[curgrid].icon[i].type == 3) {
		if (inregion(regionx, regiony, nregion, xg, yg)) {
		    writestr(xg, yg, 0, 0, "4");
		    grid[curgrid].ellist[i] = 1;
		} else {
		    grid[curgrid].ellist[i] = 0;
		}
	    }
	}
    } else {
	errwin("Region not defined");
    }
}

void do_marked_split_proc(void)
{
    int type = GetChoice(split_item);
    switch (type) {
    case 0:
	doall = 0;
	do_region_split3_proc();
	do_marked_split3_proc();
	break;
    case 1:
	doall = 1;
	do_split3_proc();
	do_marked_split3_proc();
	break;
    case 2:
	doall = 0;
	do_region_split4_proc();
	do_marked_split4_proc();
	break;
    case 3:
	doall = 1;
	do_split4_proc();
	do_marked_split4_proc();
	break;
    }
}

void do_preview_split_proc(void)
{
    int type = GetChoice(split_item);
    switch (type) {
    case 0:
	doall = 0;
	do_region_split3_proc();
	break;
    case 1:
	doall = 1;
	do_split3_proc();
	break;
    case 2:
	doall = 0;
	do_region_split4_proc();
	break;
    case 3:
	doall = 1;
	do_split4_proc();
	break;
    }
}

static void goto_node_proc(void)
{
    int nodetofind;
    int sx, sy;

    nodetofind = atoi(panel_getstr_value(goto_node_item)) - 1;
    if (nodetofind < 0 || nodetofind > grid[curgrid].nmnp) {
	errwin("Node number out of range");
	return;
    }
    get_device(grid[curgrid].xord[nodetofind], grid[curgrid].yord[nodetofind], &sx, &sy);
    setpointer(sx, sy);

}

static void goto_elem_proc(void)
{
    int elemtofind;
    int sx, sy;
    double centx, centy;

    elemtofind = atoi((char *) panel_getstr_value(goto_elem_item)) - 1;
    if (elemtofind < 0 || elemtofind > grid[curgrid].nmel) {
	errwin("Element number out of range");
	return;
    }
    get_center(curgrid, elemtofind, &centx, &centy);
    get_device(centx, centy, &sx, &sy);
    setpointer(sx, sy);
}

static void do_goto_proc(int item, int event)
{
    double wx, wy;
    int sx, sy;

    wx = atof((char *) xv_getstr(goto_pointx_item));
    wy = atof((char *) xv_getstr(goto_pointy_item));
    world2deviceabs(wx, wy, &sx, &sy);
    setpointer(sx, sy);
    getpoints(sx, sy);
}

void display_goto_proc(void)
{
    create_goto_popup();
}

void create_goto_popup(void)
{
    static Widget top, dialog;
    Widget wbut, rc, fr;
    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Goto node/element", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    fr = XmCreateFrame(dialog, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    goto_node_item = CreateTextItem2(rc, 10, "Node:");
    wbut = XtVaCreateManagedWidget("Goto node", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) goto_node_proc, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    fr = XmCreateFrame(dialog, "fr", NULL, 0);
    rc = XmCreateRowColumn(fr, "rc", NULL, 0);
    goto_elem_item = CreateTextItem2(rc, 10, "Element:");
    wbut = XtVaCreateManagedWidget("Goto element", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) goto_elem_proc, NULL);
    XtManageChild(rc);
    XtManageChild(fr);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Close", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_gotoxy_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtManageChild(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Goto", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    goto_pointx_item = CreateTextItem2(dialog, 10, "X: ");
    goto_pointy_item = CreateTextItem2(dialog, 10, "Y: ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_goto_proc, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_eval_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Evaluate", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    eval_region_item = CreatePanelChoice1(dialog, "Operate on:", 4, "Entire grid", "Defined region", "Boundary", NULL, 0);

    evaluate_item = CreateTextItem2(dialog, 20, "Function: ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_evaluate_function, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_rot_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Rotate", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    rotate_degrees_item = CreateTextItem2(dialog, 10, "Rotate by degrees = ");
    rotate_aboutx_item = CreateTextItem2(dialog, 10, "About x = ");
    rotate_abouty_item = CreateTextItem2(dialog, 10, "About y = ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_rotate_grid, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_trans_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Translate and rotate", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    translate_x_item = CreateTextItem2(dialog, 10, "Translate by X = ");
    translate_y_item = CreateTextItem2(dialog, 10, "Y = ");
    rotate_degrees_item = CreateTextItem2(dialog, 10, "Rotate by degrees = ");
    rotate_aboutx_item = CreateTextItem2(dialog, 10, "About x = ");
    rotate_abouty_item = CreateTextItem2(dialog, 10, "About y = ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_translate_grid, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

Widget *t2q_region_item;
Widget cutoff_item;

void do_triquad(void)
{
    find_quadrangles(curgrid);
}

void create_tri2quad_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Convert Tris to Quads", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    cutoff_item = CreateTextItem2(dialog, 10, "Cutoff angle (degrees) = ");
    t2q_region_item = CreatePanelChoice1(dialog, "Operate on:", 3, "Entire grid", "Defined region", NULL, 0);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_triquad, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

Widget *q2t_region_item;

void do_quadtri(void)
{
    int i;
    int a1, a2, a3, a4;
    int n = 0;
    double xg, yg;
    int doregion = GetChoice(q2t_region_item);
    int *els = (int *) malloc(grid[curgrid].nmel * sizeof(int));
    if (els == NULL) {
	errwin("Unable to allocate memory, operation cancelled");
	return;
    }
    if (doregion && !region_flag) {
	errwin("Region selected but no region defined, operation cancelled.");
	return;
    }
    if (doregion && region_flag) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    if (grid[curgrid].icon[i].nn == 4) {
		get_center(curgrid, i, &xg, &yg);
		if (inregion(regionx, regiony, nregion, xg, yg)) {
		    a1 = grid[curgrid].icon[i].nl[0];
		    a2 = grid[curgrid].icon[i].nl[1];
		    a3 = grid[curgrid].icon[i].nl[2];
		    a4 = grid[curgrid].icon[i].nl[3];
		    add_element(curgrid, a1, a2, a3);
		    add_element(curgrid, a1, a3, a4);
		    els[n] = i;
		    n++;
		}
	    }
	}
    } else {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    if (grid[curgrid].icon[i].nn == 4) {
		a1 = grid[curgrid].icon[i].nl[0];
		a2 = grid[curgrid].icon[i].nl[1];
		a3 = grid[curgrid].icon[i].nl[2];
		a4 = grid[curgrid].icon[i].nl[3];
		add_element(curgrid, a1, a2, a3);
		add_element(curgrid, a1, a3, a4);
		els[n] = i;
		n++;
	    }
	}
    }
    if (n) {
	delete_elements(curgrid, els, n);
	do_drawgrid();
    } else {
	errwin("No quads found to convert, operation cancelled");
    }
    free(els);
}

void create_quad2tri_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Convert Quads to Tris", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    q2t_region_item = CreatePanelChoice1(dialog, "Operate on:", 3, "Entire grid", "Defined region", NULL, 0);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_quadtri, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

Widget quadtoggle_item;
Widget tritoggle_item;
Widget *quadcolor_item;
Widget *tricolor_item;

int drawquadflag = 0;
int drawtriflag = 0;
int quadcolor = 2;
int tricolor = 1;

void do_highlight(void)
{
    drawquadflag = XmToggleButtonGetState(quadtoggle_item);
    drawtriflag = XmToggleButtonGetState(tritoggle_item);
    quadcolor = GetChoice(quadcolor_item);
    tricolor = GetChoice(tricolor_item);
    do_drawgrid();
}

void create_highlight_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Hightlight Quads/Tris", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    quadtoggle_item = XtVaCreateManagedWidget("Display quads", xmToggleButtonWidgetClass, dialog, NULL);
    quadcolor_item = CreateColorChoice(dialog, "Quad color:", 1);
    tritoggle_item = XtVaCreateManagedWidget("Display tris", xmToggleButtonWidgetClass, dialog, NULL);
    tricolor_item = CreateColorChoice(dialog, "Tri color:", 1);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_highlight, 0);

    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_equiv_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Acceptable skewness", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    equiv_radius_item = CreateTextItem2(dialog, 10, "Acceptable skewness: ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_marked_delete_proc, 0);
    wbut = XtVaCreateManagedWidget("Preview", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_mark_proc, NULL);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

double quadqual(int gridno, int elem)
{
    char buf[1024];
    int i, j, k, n1, n2, n3, n4, nn[4];
    double cutoff, test, x1, y1, x2, y2, x3, y3, ang, q[4], minq, f1, g1, f2, g2, xg, yg;
    if (grid[gridno].icon[elem].nn == 4) {
	nn[0] = grid[gridno].icon[elem].nl[0];
	nn[1] = grid[gridno].icon[elem].nl[1];
	nn[2] = grid[gridno].icon[elem].nl[2];
	nn[3] = grid[gridno].icon[elem].nl[3];
	q[k] = 0;
	for (k = 0; k < 4; k++) {
	    x1 = grid[gridno].xord[nn[k]];
	    x2 = grid[gridno].xord[nn[(k + 1) % 4]];
	    x3 = grid[gridno].xord[nn[(k + 2) % 4]];
	    y1 = grid[gridno].yord[nn[k]];
	    y2 = grid[gridno].yord[nn[(k + 1) % 4]];
	    y3 = grid[gridno].yord[nn[(k + 2) % 4]];
	    f1 = (x2 - x1);
	    f2 = (x3 - x2);
	    g1 = (y2 - y1);
	    g2 = (y3 - y2);
	    if (fabs(f1 * f2 + g1 * g2 - (hypot(f1, g1) * hypot(f2, g2))) < 0.0000005) {
		ang = M_PI;
	    } else {
		ang = acos((f1 * f2 + g1 * g2) / (hypot(f1, g1) * hypot(f2, g2)));
	    }
	    test = 1.0 - fabs(1.0 - (2.0 * ang / M_PI));
	    q[k] = test > 0.0 ? test : 0.0;
	}			/* for k */
	minq = q[0];
	minq = (minq > q[1]) ? q[1] : minq;
	minq = (minq > q[2]) ? q[2] : minq;
	minq = (minq > q[3]) ? q[3] : minq;
    }				/* if ngeom == 4 */
    return minq;
}

static Widget qual_item;

int WhichSide(double x1, double y1, double x2, double y2, double x3, double y3)
	/* p q r */
{
    double result;
    result = (x1 - x2) * (y2 - y3) - (y1 - y2) * (x2 - x3);
    if (result < 0)
	return -1;		/* q lies to the left  (qr turns CW).   */
    if (result > 0)
	return 1;		/* q lies to the right (qr turns CCW).  */
    return 0;			/* q lies on the line from p to r.      */
}

void qstufftext(char *s, int sp);
static void qclear_results(void);
static int qcount;

void do_qual_proc(void)
{
    char s[31];
    char buf[1024];
    int i, j, k, n1, n2, n3, n4, nn[4];
    double cutoff, test, x1, y1, x2, y2, x3, y3, ang, q[4], minq, f1, g1, f2, g2, xg, yg;
    strncpy(s, (char *) panel_getstr_value(qual_item), 30);
    cutoff = atof(s);
    if (cutoff <= 0.0 || cutoff > 1.0) {
	errwin("Cutoff must be between 0.0 and 1.0, operation cancelled");
	return;
    }
    qcount = 0;
    qclear_results();
    sprintf(buf, "\nChecking quad elements for quality at cutoff %.3lf\n", cutoff);
    qstufftext(buf, 1);

    for (i = 0; i < grid[curgrid].nmel; i++) {
	if (grid[curgrid].icon[i].nn == 4) {
	    nn[0] = grid[curgrid].icon[i].nl[0];
	    nn[1] = grid[curgrid].icon[i].nl[1];
	    nn[2] = grid[curgrid].icon[i].nl[2];
	    nn[3] = grid[curgrid].icon[i].nl[3];
	    for (k = 0; k < 4; k++) {
		x1 = grid[curgrid].xord[nn[k]];
		x2 = grid[curgrid].xord[nn[(k + 1) % 4]];
		x3 = grid[curgrid].xord[nn[(k + 2) % 4]];
		y1 = grid[curgrid].yord[nn[k]];
		y2 = grid[curgrid].yord[nn[(k + 1) % 4]];
		y3 = grid[curgrid].yord[nn[(k + 2) % 4]];
		if (WhichSide(x1, y1, x2, y2, x3, y3) < 1) {
		    sprintf(buf, "Concave or degenerate at element %d\n", i + 1);
		    qstufftext(buf, 0);
		    q[k] = 0;
		    break;
		}

		f1 = (x2 - x1);
		f2 = (x3 - x2);
		g1 = (y2 - y1);
		g2 = (y3 - y2);
		if (fabs(f1 * f2 + g1 * g2 - (hypot(f1, g1) * hypot(f2, g2))) < 0.0000005) {
		    ang = M_PI;
		} else {
		    ang = acos((f1 * f2 + g1 * g2) / (hypot(f1, g1) * hypot(f2, g2)));
		}
		test = 1.0 - fabs(1.0 - (2.0 * ang / M_PI));
		q[k] = test > 0.0 ? test : 0.0;
		/*
		   printf("elem, k, nn, q, ang, side = %d %d %d %lf %lf %d\n", i + 1, k, nn[k], q[k], ang * 180.0 / M_PI,
		   WhichSide(x1, y1, x2, y2, x3, y3));
		 */
	    }			/* for k */
	    minq = q[0];
	    minq = (minq > q[1]) ? q[1] : minq;
	    minq = (minq > q[2]) ? q[2] : minq;
	    minq = (minq > q[3]) ? q[3] : minq;
	    if (minq < cutoff) {
		sprintf(buf, "Element %d fails, minq = %.3lf\n", i + 1, minq);
		qstufftext(buf, 0);
		get_center(curgrid, i, &xg, &yg);
		solidbox(xg, yg);
	    }
	}			/* if ngeom == 4 */
    }				/* for i */
    sprintf(buf, "Done check.\n");
    qstufftext(buf, 2);
}

Widget qtext;

static void qclear_results(void)
{
    extern int inwin;
    if (inwin) {
	XmTextSetString(qtext, "");
    }
}

void qstufftext(char *s, int sp)
{
    extern int inwin;
    static XmTextPosition pos = 0;
    static XmTextPosition savepos = 0;
    qcount++;
    if (qcount > 100) {
	    return;
    }

    if (inwin) {
	if (sp == 1) {
	    pos = XmTextGetLastPosition(qtext);
	    savepos = pos;
	    XmTextSetTopCharacter(qtext, savepos);
	}
	XmTextInsert(qtext, pos, s);
	pos += strlen(s);
	if (sp == 2) {
	    XmTextSetTopCharacter(qtext, savepos);
	    savepos = pos;
	}
    } else {
	printf(s);
    }
}

void create_qqual_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc, fr;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Quality check for quadrangles", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    qual_item = CreateTextItem2(dialog, 10, "Cutoff (0 < c <= 1): ");

    fr = XmCreateFrame(dialog, "fr", NULL, 0);
    qtext = XmCreateScrolledText(fr, "mon", NULL, 0);
    XtVaSetValues(qtext, XmNrows, 10, XmNcolumns, 60, XmNeditMode, XmMULTI_LINE_EDIT, XmNwordWrap, True, NULL);
    XtManageChild(qtext);

    XtVaSetValues(fr, XmNtopAttachment, XmATTACH_FORM, XmNleftAttachment, XmATTACH_FORM, XmNrightAttachment, XmATTACH_FORM, XmNbottomAttachment, XmATTACH_FORM, NULL);
    XtManageChild(fr);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_qual_proc, 0);
    wbut = XtVaCreateManagedWidget("Clear", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) qclear_results, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}


void create_splitr_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Split in region", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    split_item = CreatePanelChoice1(dialog, "Split elements:", 5, "into 3 in region", "into 3 in entire grid", "into 4 in region", "into 4 in entire grid", NULL, 0);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_marked_split_proc, NULL);
    wbut = XtVaCreateManagedWidget("Preview", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_preview_split_proc, NULL);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

void create_delr_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Delete elements in region", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    delete_toggle_item = CreatePanelChoice1(dialog, "Delete:", 3, "All elements inside region:", "All elements outside region:", NULL, 0);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_marked_delete_proc, NULL);

    wbut = XtVaCreateManagedWidget("Mark elements to delete", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_region_delete_proc, NULL);

    wbut = XtVaCreateManagedWidget("Cancel", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

static Widget *vol_region_item;
static Widget vol_val_item;
static Widget *vol_grid_item;

double areapts(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double a[2], b[2], ar;

    a[0] = x3 - x2;
    a[1] = x1 - x3;
    b[0] = y2 - y3;
    b[1] = y3 - y1;
    return 0.5 * (a[1] * b[0] - a[0] * b[1]);
}

double volumepts(double x1, double y1, double d1, double x2, double y2, double d2, double x3, double y3, double d3)
{
    double ar = areapts(x1, y1, x2, y2, x3, y3);
    return ar * (d1 + d2 + d3) * 0.3333333333333;
}

void do_volume(void)
{
    int i, j, itmp, inode, ib, ind, gridno = curgrid, errpos;
    int n1, n2, n3;
    double v, x[3], y[3];
    double element_volume(int gridno, int elno);
    char buf[128];

    itmp = (int) GetChoice(vol_region_item);
    gridno = (int) GetChoice(vol_grid_item);
    if (gridno == 1) {
	gridno = MAXGRIDS;
    } else {
	gridno = curgrid;
    }
    switch (itmp) {
    case 0:
	v = 0.0;
	for (i = 0; i < grid[gridno].nmel; i++) {
	    v += element_volume(gridno, i);
	}
	sprintf(buf, "%lf", v);
	xv_setstr(vol_val_item, buf);
	break;
    case 1:
	if (region_flag) {
	    int t1, t2, t3;
	    v = 0.0;
	    for (i = 0; i < grid[gridno].nmel; i++) {
		n1 = grid[gridno].icon[i].nl[0];
		n2 = grid[gridno].icon[i].nl[1];
		n3 = grid[gridno].icon[i].nl[2];
		x[0] = grid[gridno].xord[n1];
		y[0] = grid[gridno].yord[n1];
		x[1] = grid[gridno].xord[n2];
		y[1] = grid[gridno].yord[n2];
		x[2] = grid[gridno].xord[n3];
		y[2] = grid[gridno].yord[n3];
		t1 = inregion(regionx, regiony, nregion, x[0], y[0]);
		t2 = inregion(regionx, regiony, nregion, x[1], y[1]);
		t3 = inregion(regionx, regiony, nregion, x[2], y[2]);
		if (t1 && t2 && t3) {
		    v += element_volume(gridno, i);
		    fillelement(gridno, i, 2);
		} else {
		}
	    }
	    sprintf(buf, "%lf", v);
	    xv_setstr(vol_val_item, buf);
	}
	break;
    }
}

void integrate(int gridno, int el1, double x1, double y1, double d1, double x2, double y2, double d2, double x3, double y3, double d3, int limit, double *sum)
{
    if (limit > 0) {
	if (gridbelel(gridno, el1, x1, y1) && gridbelel(gridno, el1, x2, y2) && gridbelel(gridno, el1, x3, y3)) {
	    *sum += volumepts(x1, y1, d1, x2, y2, d2, x3, y3, d3);
	} else {
/* p1 p(1+2) p(1+3) */
	    integrate(gridno, el1, x1, y1, d1, 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (d1 + d2), 0.5 * (x1 + x3), 0.5 * (y1 + y3), 0.5 * (d1 + d3), limit - 1, sum);
/* p2 p(3+2) p(1+2) */
	    integrate(gridno, el1, x2, y2, d2, 0.5 * (x3 + x2), 0.5 * (y3 + y1), 0.5 * (d3 + d2), 0.5 * (x2 + x1), 0.5 * (y2 + y1), 0.5 * (d2 + d1), limit - 1, sum);
/* p3 p(3+1) p(3+2) */
	    integrate(gridno, el1, x3, y3, d3, 0.5 * (x1 + x3), 0.5 * (y1 + y3), 0.5 * (d1 + d3), 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (d1 + d2), limit - 1, sum);
/* p(3+1) p(2+1) p(3+2) */
	    integrate(gridno, el1, 0.5 * (x1 + x3), 0.5 * (y1 + y3), 0.5 * (d1 + d3), 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (d1 + d2), 0.5 * (x3 + x2), 0.5 * (y3 + y1), 0.5 * (d3 + d2), limit - 1, sum);
/*
   if (gridbelel(gridno, el1, x1, y1)) {
   integrate(gridno, el1, x1, y1, d1, 0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (d1 + d2),
   0.5 * (x1 + x3), 0.5 * (y1 + y3), 0.5 * (d1 + d3), limit - 1, sum);
   }
   if (gridbelel(gridno, el1, x2, y2)) {
   integrate(gridno, el1, x2, y2, d2, 0.5 * (x3 + x2), 0.5 * (y3 + y1), 0.5 * (d3 + d2),
   0.5 * (x2 + x1), 0.5 * (y2 + y1), 0.5 * (d2 + d1), limit - 1, sum);
   }
   if (gridbelel(gridno, el1, x3, y3)) {
   integrate(gridno, el1, x3, y3, d3, 0.5 * (x1 + x3), 0.5 * (y1 + y3), 0.5 * (d1 + d3),
   0.5 * (x1 + x2), 0.5 * (y1 + y2), 0.5 * (d1 + d2), limit - 1, sum);
   }
 */
	}
    }
}

void create_volume_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	return;
    }
    top = XmCreateDialogShell(app_shell, "Volumes", NULL, 0);
    handle_close(top);
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    vol_grid_item = CreatePanelChoice1(dialog, "Grid:", 3, "Edit", "Background grid", NULL, 0);
    vol_region_item = CreatePanelChoice1(dialog, "Operate on:", 3, "Entire grid", "Defined region", NULL, 0);

    vol_val_item = CreateTextItem2(dialog, 20, "Volume: ");

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_volume, 0);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    XtManageChild(dialog);
    XtManageChild(top);
}

/*
 * Added CPP project routines 10-95 pturner
 * The proj library doesn't contain this one.
 */
static Widget cpp_frame;

static Widget *cppfrom_item;
static Widget *cppto_item;
static Widget *cppapply_item;
static Widget cpp_text;
static Widget cppinv_item;
static Widget cpplam0_text;
static Widget cppphi1_text;
static int cppinv_flag = 0;
static double cppradius = 6378206.4;
static double cpplam0, cppphi1;

void cpp_forward(double lon, double lat, double *x, double *y)
{
    *x = cppradius * (lon - cpplam0) * cos(cppphi1);
    *y = lat * cppradius;
}

void cpp_inverse(double x, double y, double *lon, double *lat)
{
    *lat = y / cppradius;
    *lon = cpplam0 + x / (cppradius * cos(cppphi1));
}

void do_cpp(double *x, double *y, double *u, double *v, int n)
{
    int i, inverse = cppinv_flag;
    double utmp, vtmp, xres, yres;
    for (i = 0; i < n; i++) {
	if (inverse) {
	    utmp = x[i];
	    vtmp = y[i];
	} else {
	    utmp = x[i] * M_PI / 180.0;
	    vtmp = y[i] * M_PI / 180.0;
	}
	if (inverse) {
	    cpp_inverse(utmp, vtmp, &xres, &yres);
	    u[i] = xres * 180.0 / M_PI;
	    v[i] = yres * 180.0 / M_PI;
	} else {
	    cpp_forward(utmp, vtmp, &xres, &yres);
	    u[i] = xres;
	    v[i] = yres;
	}
    }
}

void do_cpp_proc(void)
{
    char buf[256];
    int i, j, bno, gno;
    int applyto = GetChoice(cppapply_item);
    cpplam0 = atof(xv_getstr(cpplam0_text));
    cppphi1 = atof(xv_getstr(cppphi1_text));
    cppinv_flag = XmToggleButtonGetState(cppinv_item);
    switch (applyto) {
    case 0:
	for (i = 0; i <= MAXGRIDS; i++) {
	    if (grid[i].gactive) {
		do_cpp(grid[i].xord, grid[i].yord, grid[i].xord, grid[i].yord, grid[i].nmnp);
		for (j = 0; j < grid[i].nbounds; j++) {
		    bno = grid[i].boundaries[j];
		    if (bno >= 0 && boundary[bno].bactive) {
			do_cpp(boundary[bno].boundx, boundary[bno].boundy, boundary[bno].boundx, boundary[bno].boundy, boundary[bno].nbpts);
		    }
		}
	    }
	}
	break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
	gno = applyto - 1;
	if (grid[gno].gactive) {
	    do_cpp(grid[gno].xord, grid[gno].yord, grid[gno].xord, grid[gno].yord, grid[gno].nmnp);
	    for (j = 0; j < grid[gno].nbounds; j++) {
		bno = grid[gno].boundaries[j];
		if (bno >= 0 && boundary[bno].bactive) {
		    do_cpp(boundary[bno].boundx, boundary[bno].boundy, boundary[bno].boundx, boundary[bno].boundy, boundary[bno].nbpts);
		}
	    }
	} else {
	    errwin("Grid not defined");
	}
	break;
    case 7:
	if (build[0].nbuild > 0) {
	    do_cpp(build[0].bx, build[0].by, build[0].bx, build[0].by, build[0].nbuild);
	}
	break;
    case 8:
	for (j = 0; j < grid[MAXGRIDS].nbounds; j++) {
	    bno = grid[MAXGRIDS].boundaries[j];
	    if (bno >= 0 && boundary[bno].bactive) {
		do_cpp(boundary[bno].boundx, boundary[bno].boundy, boundary[bno].boundx, boundary[bno].boundy, boundary[bno].nbpts);
	    }
	}
	break;
    }
}

void update_cpp(void)
{
    if (cpp_frame) {
    }
}

void create_cpp_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	update_cpp();
	return;
    }
    top = XmCreateDialogShell(app_shell, "CPP map projection", NULL, 0);
    handle_close(top);
    cpp_frame = top;
    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    cppapply_item = CreatePanelChoice1(dialog, "Apply to:", 10, "All", "Grid 1 (Edit)", "Grid 2", "Grid 3", "Grid 4", "Grid 5", "Background grid", "Build points", "Coastal boundaries", NULL, 0);
    cpplam0_text = CreateTextItem2(dialog, 9, "Central meridian: ");
    cppphi1_text = CreateTextItem2(dialog, 9, "Central latitude: ");
    cppinv_item = XtVaCreateManagedWidget("Invert projection (X, Y -> Lon, Lat)", xmToggleButtonWidgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_cpp_proc, NULL);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    update_cpp();
    XtManageChild(dialog);
    XtManageChild(top);
}

static int extractto;
static Widget extract_to_item;

void do_extract_grid(void)
{
    if (extractto == curgrid) {
	errwin("Grid to extract from and extract to are the same, operation cancelled!");
	set_action(0);
    } else {
	set_action(0);
	set_action(EXTRACT_GRID);
    }
}

void set_extractto_proc(int w, int cd)
{
    extractto = cd;
}

void extract_grid(double *polyx, double *polyy, int n)
{
    int i, j, elcnt = 0, nodcnt = 0, nod, nn;
    double xg, yg;
    if (n > 2) {
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    grid[curgrid].ellist[i] = 0;
	    get_center(curgrid, i, &xg, &yg);
	    if (inregion(polyx, polyy, n, xg, yg)) {
		writestr(xg, yg, 0, 0, "1");
		grid[curgrid].ellist[i] = elcnt;
		elcnt++;
	    } else {
		grid[curgrid].ellist[i] = -1;
	    }
	}
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    grid[curgrid].nlist[i] = 0;
	}
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    if (grid[curgrid].ellist[i] != -1) {
		for (j = 0; j < numnodes(curgrid, i); j++) {
		    nod = grid[curgrid].icon[i].nl[j];
		    grid[curgrid].nlist[nod] = 1;
		}
	    }
	}
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    if (grid[curgrid].nlist[i]) {
		grid[curgrid].nlist[i] = nodcnt++;
	    } else {
		grid[curgrid].nlist[i] = -1;
	    }
	}
	Free_grid(extractto);
	Allocate_grid(extractto, elcnt, nodcnt);
	nodcnt = elcnt = 0;
	for (i = 0; i < grid[curgrid].nmnp; i++) {
	    if (grid[curgrid].nlist[i] != -1) {
		grid[extractto].xord[nodcnt] = grid[curgrid].xord[i];
		grid[extractto].yord[nodcnt] = grid[curgrid].yord[i];
		grid[extractto].depth[nodcnt] = grid[curgrid].depth[i];
		nodcnt++;
	    }
	}
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    if (grid[curgrid].ellist[i] != -1) {
		for (j = 0; j < numnodes(curgrid, i); j++) {
		    nn = grid[curgrid].icon[i].nl[j];
		    nod = grid[curgrid].nlist[nn];
		    grid[extractto].icon[elcnt].nl[j] = nod;
		}
		grid[extractto].icon[elcnt].type = grid[curgrid].icon[i].type;
		elcnt++;
	    }
	}
	set_grid_limits(extractto);
	default_isolines(&grid[extractto].ip);
    } else {
	errwin("Region not defined");
    }
}

void create_extract_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;
    Widget bt, fr2, rc2;
    int x, y;
    if (top) {
	XtRaise(top);
	return;
    }
    XmGetPos(app_shell, 0, &x, &y);
    top = XmCreateDialogShell(app_shell, "Extract grid", NULL, 0);
    XtVaSetValues(top, XmNx, x, XmNy, y, NULL);
    handle_close(top);

    dialog = XmCreateRowColumn(top, "rc", NULL, 0);

    extract_to_item = XmVaCreateSimpleOptionMenu(dialog, "option_menu",
						 XmStringCreateSimple("Extract sub-grid to grid: "), 'a', 0,
						 (XtCallbackProc) set_extractto_proc,
						 XmVaPUSHBUTTON, XmStringCreateSimple("0"), 0, NULL, NULL,
						 XmVaPUSHBUTTON, XmStringCreateSimple("1"), 0, NULL, NULL,
						 XmVaPUSHBUTTON, XmStringCreateSimple("2"), 0, NULL, NULL,
						 XmVaPUSHBUTTON, XmStringCreateSimple("3"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("4"), 0, NULL, NULL, XmVaPUSHBUTTON, XmStringCreateSimple("Background"), 0, NULL, NULL, NULL);
    XtManageChild(extract_to_item);

    fr2 = XmCreateFrame(dialog, "fr", NULL, 0);
    rc2 = XmCreateRowColumn(fr2, "rc", NULL, 0);
    XtVaSetValues(rc2, XmNorientation, XmHORIZONTAL, NULL);
    bt = XtVaCreateManagedWidget("Apply", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) do_extract_grid, NULL);
    bt = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc2, NULL);
    XtAddCallback(bt, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc2);
    XtManageChild(fr2);

    XtManageChild(dialog);
    XtManageChild(top);
}

static Widget ftri_frame;
static Widget ftri_len_item;
static Widget *ftri_region_item;

void do_ftri_proc(void)
{
    compact_grid(curgrid);
}

void do_ftrimark_proc(void)
{
    char buf[256];
    int i, j, n1, n2, n3, which;
    double len, tlen1, tlen2, x[3], y[3];
    double xg, yg;
    strcpy(buf, (char *) panel_getstr_value(ftri_len_item));
    sscanf(buf, "%lf", &len);
    which = GetChoice(ftri_region_item);
    switch (which) {
    case 0:
	for (i = 0; i < grid[curgrid].nmel; i++) {
	    n1 = grid[curgrid].icon[i].nl[0];
	    n2 = grid[curgrid].icon[i].nl[1];
	    n3 = grid[curgrid].icon[i].nl[2];
	    x[0] = grid[curgrid].xord[n1];
	    y[0] = grid[curgrid].yord[n1];
	    x[1] = grid[curgrid].xord[n2];
	    y[1] = grid[curgrid].yord[n2];
	    x[2] = grid[curgrid].xord[n3];
	    y[2] = grid[curgrid].yord[n3];
	    tlen1 = hypot((x[1] - x[0]), (y[1] - y[0]));
	    tlen2 = hypot((x[2] - x[1]), (y[2] - y[1]));
	    if (tlen1 < tlen2)
		tlen1 = tlen2;
	    tlen2 = hypot((x[0] - x[2]), (y[0] - y[2]));
	    if (tlen1 < tlen2)
		tlen1 = tlen2;
	    if (tlen1 > len) {	/* delete element */
		grid[curgrid].ellist[i] = 1;
	    } else {
		grid[curgrid].ellist[i] = 0;
	    }
	}
	break;
    case 1:
	if (region_flag) {
	    for (i = 0; i < grid[curgrid].nmel; i++) {
		get_center(curgrid, i, &xg, &yg);
		if (inregion(regionx, regiony, nregion, xg, yg)) {
		    n1 = grid[curgrid].icon[i].nl[0];
		    n2 = grid[curgrid].icon[i].nl[1];
		    n3 = grid[curgrid].icon[i].nl[2];
		    x[0] = grid[curgrid].xord[n1];
		    y[0] = grid[curgrid].yord[n1];
		    x[1] = grid[curgrid].xord[n2];
		    y[1] = grid[curgrid].yord[n2];
		    x[2] = grid[curgrid].xord[n3];
		    y[2] = grid[curgrid].yord[n3];
		    tlen1 = hypot((x[1] - x[0]), (y[1] - y[0]));
		    tlen2 = hypot((x[2] - x[1]), (y[2] - y[1]));
		    if (tlen1 < tlen2)
			tlen1 = tlen2;
		    tlen2 = hypot((x[0] - x[2]), (y[0] - y[2]));
		    if (tlen1 < tlen2)
			tlen1 = tlen2;
		    if (tlen1 > len) {	/* delete element */
			solidbox(xg, yg);
			grid[curgrid].ellist[i] = 1;
		    } else {
			grid[curgrid].ellist[i] = 0;
		    }
		}
	    }
	} else {
	    errwin("Region not defined");
	}
	break;
    case 2:
	if (region_flag) {
	    for (i = 0; i < grid[curgrid].nmel; i++) {
		get_center(curgrid, i, &xg, &yg);
		if (!inregion(regionx, regiony, nregion, xg, yg)) {
		    n1 = grid[curgrid].icon[i].nl[0];
		    n2 = grid[curgrid].icon[i].nl[1];
		    n3 = grid[curgrid].icon[i].nl[2];
		    x[0] = grid[curgrid].xord[n1];
		    y[0] = grid[curgrid].yord[n1];
		    x[1] = grid[curgrid].xord[n2];
		    y[1] = grid[curgrid].yord[n2];
		    x[2] = grid[curgrid].xord[n3];
		    y[2] = grid[curgrid].yord[n3];
		    tlen1 = hypot((x[1] - x[0]), (y[1] - y[0]));
		    tlen2 = hypot((x[2] - x[1]), (y[2] - y[1]));
		    if (tlen1 < tlen2)
			tlen1 = tlen2;
		    tlen2 = hypot((x[0] - x[2]), (y[0] - y[2]));
		    if (tlen1 < tlen2)
			tlen1 = tlen2;
		    if (tlen1 > len) {	/* delete element */
			solidbox(xg, yg);
			grid[curgrid].ellist[i] = 1;
		    } else {
			grid[curgrid].ellist[i] = 0;
		    }
		}
	    }
	} else {
	    errwin("Region not defined");
	}
	break;
    }
}

void update_ftri(void)
{
    if (ftri_frame) {
    }
}

void create_ftri_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
	XtRaise(top);
	update_ftri();
	return;
    }
    top = XmCreateDialogShell(app_shell, "Filter triangles", NULL, 0);
    handle_close(top);
    ftri_frame = top;

    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    ftri_len_item = CreateTextItem2(dialog, 30, "Minimum side length: ");
    ftri_region_item = CreatePanelChoice1(dialog, "Apply to:", 4, "entire grid", "inside region", "outside region", NULL, 0);

    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_ftri_proc, NULL);
    wbut = XtVaCreateManagedWidget("Mark", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_ftrimark_proc, NULL);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    update_ftri();
    XtManageChild(dialog);
    XtManageChild(top);
}

static Widget shapiro_frame;
static Widget toggle_shapiro_item;

int display_shapiro = 0;

void do_shapiro_proc(void)
{
//     do_drawgrid();
//     find_sh(curgrid);
    display_shapiro = XmToggleButtonGetState(toggle_shapiro_item);
    do_drawgrid();
}

void update_shapiro(void)
{
    if (shapiro_frame) {
	XmToggleButtonSetState(toggle_shapiro_item, display_shapiro, False);
    }
}

void create_shapiro_frame(void)
{
    static Widget top, dialog;
    Widget wbut, rc;

    if (top) {
        XtRaise(top);
        update_shapiro();
        return;
    }
    top = XmCreateDialogShell(app_shell, "Show edge midpoints", NULL, 0);
    handle_close(top);
    shapiro_frame = top;

    dialog = XmCreateRowColumn(top, "rc", NULL, 0);
    toggle_shapiro_item = XtVaCreateManagedWidget("Display Shapiro filter violations", xmToggleButtonWidgetClass, dialog, NULL);
    rc = XmCreateRowColumn(dialog, "rc", NULL, 0);
    XtVaSetValues(rc, XmNorientation, XmHORIZONTAL, NULL);
    //wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    //XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_shapiro_proc, NULL);
    wbut = XtVaCreateManagedWidget("Accept", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) do_shapiro_proc, NULL);
    wbut = XtVaCreateManagedWidget("Done", xmPushButtonWidgetClass, rc, NULL);
    XtAddCallback(wbut, XmNactivateCallback, (XtCallbackProc) destroy_dialog, top);
    XtManageChild(rc);

    update_shapiro();
    XtManageChild(dialog);
    XtManageChild(top);
}
