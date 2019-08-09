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
 * canvas event proc and set_action()
 *
 */

#ifndef lint
static char RCSid[] = "$Id: events.c,v 1.6 2011/09/14 17:44:21 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "graphics.h"

#define MOVE_WORLD1ST 997
#define MOVE_WORLD2ND 998
#define CENTER_WORLD 999

Widget mode_item;
Widget locate_grid_item;
Widget calc_item;
Widget locate_item;		/* locator on main_panel */
Widget comparea_item;
Widget compperi_item;

extern int go_locateflag;

extern Display *disp;
extern Window cwin;
extern GC gc;
extern GC gcxor;
extern GC gcclr;

extern int spreadflag;
extern int cercnum;

void DrawPNodeLists(int gridno);
void AddToPNodeList(int nlno, int node, int ival, double dval, int color);
extern int curnl;
extern int curnlcolor;
extern int curnlival;
extern double curnldval;

double area(), area_from_nodes(), get_depth_element();
void set_action(int act);

static char buf[256];

static int action_flag = 0;
static int call_flag = 0;

static int rectflag = 0;
static int rotrectflag = 0;
static int rubber_flag = 0;
int mbox_flag = 0;		/* moving box attached to cursor */
int mline_flag = 0;		/* moving line attached to cursor */

static double wx1, wx2, wy1, wy2, wx3, wy3, wx, wy;
static int sx, sy, old_x, old_y;	/* for rubber banding  and rects */
static int sx1, sy1;		/* for rotated rect rubberbanding */
static int sx2, sy2;		/* for rotated rect rubberbanding */
static int sx3, sy3;		/* for rotated rect rubberbanding */
static int tmpind;
static int bnd1, bnd2, ind1, ind2, ind3, ind4, ind5;	/* for adding elements */
static double sm, sb;		/* for snap lines */
static int bnof, indf;
static int n0, n1, n2, n3, n4, n5, n6, n7;

extern int deltaflag;
extern int pointset;
double dsx, dsy;

/* for 1-d slices */
extern double slicex[], slicey[];
extern int slice_npts;
extern int win_h, win_w;

/*
 * variables for the text handling routine
 */
static int strx = 0, stry = 0;
static int drawx = 0, drawy = 0;
static char tmpstr[256];
static int justflag = 0;
static double si = 0.0;
static double co = 1.0;
static int xs, ys;
static int cg = 0;

void doSync(void)
{
    XSync(disp, False);
}

void setpointer(int x, int y)
{
    XWarpPointer(disp, None, cwin, 0, None, (unsigned int) win_w, (unsigned int) win_h, x, y);
    XFlush(disp);
}

/*
 * rubber band line
 */
void select_line(int x1, int y1, int x2, int y2)
{
    XDrawLine(disp, cwin, gcxor, x1, y1, x2, y2);
}

/*
 * draw rotated box
 */
void select_rotated_rect(int x1, int y1, int x2, int y2, int x3, int y3)
{
    double rot;
    XDrawLine(disp, cwin, gcxor, x1, y1, x2, y2);
    XDrawLine(disp, cwin, gcxor, x2, y2, x3, y3);
    XDrawLine(disp, cwin, gcxor, x3, y3, x1, y1);
}

/*
 * draw a box on the display
 */
void draw_rectangle(int x1, int y1, int x2, int y2)
{
    XDrawRectangle(disp, cwin, gc, x1, y1, x2, y2);
}

/*
 * draw an xor'ed box
 */
void select_region(int x1, int y1, int x2, int y2)
{
    int dx = x2 - x1;
    int dy = y2 - y1;

    if (dx < 0) {
	iswap(&x1, &x2);
	dx = -dx;
    }
    if (dy < 0) {
	iswap(&y1, &y2);
	dy = -dy;
    }
    XDrawRectangle(disp, cwin, gcxor, x1, y1, dx, dy);
}

/*
 * switch on the area calculator
 */
void do_select_area(void)
{
    set_action(0);
    set_action(COMP_AREA);
}

/*
 * switch on the perimeter calculator
 */
void do_select_peri(void)
{
    set_action(0);
    set_action(COMP_PERIMETER);
}

/*
 * evaluate an expression for the calculator
 */
void do_calc_proc(void)
{
    static double a = 0.0;
    static double b = 0.0;
    static double c = 0.0;
    static double d = 0.0;
    static double x = 0.0;
    static double y = 0.0;
    static int errpos;
    static char val[128];
    extern double result;

    errpos = 0;
    strcpy(val, (char *) panel_getstr_value(calc_item));
    fixupstr(val);
    scanner(val, &x, &y, &a, &b, &c, &d, 0, 0, &errpos);
    if (errpos) {
    } else {
	sprintf(val, "%.5g", result);
    }
}

/*
 * select the type of display for locator
 */
void do_setdeltype(int value)
{
    deltaflag = value;
}

void reset_action(void)
{
    set_action(action_flag);
}

/*
 * set the action_flag to the desired action (actions are
 * defined in defines.h), if 0 then cleanup the results
 * from previous actions.
 */
void set_action(int act)
{
    int i;
    extern Widget canvas;

    set_window(canvas);

    if ((action_flag = act) == 0) {	/* clean up */
	nibuf = 0;
	npoly = 0;
	nipoly = 0;
	write_mode_str("Idle ...");
	if (rectflag) {
	    select_region(sx, sy, old_x, old_y);
	    rectflag = 0;
	}
	if (rotrectflag) {
	    select_rotated_rect(sx1, sy1, sx2, sy2, sx3, sy3);
	    rotrectflag = 0;
	}
	if (rubber_flag) {
	    select_line(sx, sy, old_x, old_y);
	    rubber_flag = 0;
	}
	if (mbox_flag) {
	    select_region(sx, sy, xs, ys);
	    mbox_flag = 0;
	}
	if (mline_flag) {
	    select_line(sx, sy, xs, ys);
	    mline_flag = 0;
	}
	set_cursor(-1);
    } else {
	switch (act) {
	case DEFINE_REGION:
	    set_cursor(0);
	    break;
	case SELECT_SNAP:
	    set_cursor(0);
	    write_mode_str("Click near nodes to snap");
	    break;
	case SNAP_LINE1ST:
	    set_cursor(0);
	    write_mode_str("Click at the beginning of the line to snap to");
	    break;
	case SNAP_LINE2ND:
	    set_cursor(0);
	    write_mode_str("Click at the end of the line to snap to");
	    break;
	case DEF_EXT_BOUND:
	    write_mode_str("Click at a point to include into the exterior boundary");
	    set_cursor(0);
	    break;
	case DEF_INT_BOUND:
	    write_mode_str("Click at a point to include into an interior boundary");
	    set_cursor(0);
	    break;
	case MOVE_EXT_BOUND1ST:
	    write_mode_str("Click near the external boundary point to move");
	    set_cursor(0);
	    break;
	case MOVE_EXT_BOUND2ND:
	    write_mode_str("Click at the new location of the external boundary point");
	    set_cursor(0);
	    break;
	case MOVE_INT_BOUND1ST:
	    write_mode_str("Click near the internal boundary point to move");
	    set_cursor(0);
	    break;
	case MOVE_INT_BOUND2ND:
	    write_mode_str("Click at the new location of the internal boundary point");
	    set_cursor(0);
	    break;
	case DEL_EXT_BOUNDPT:
	    set_cursor(3);
	    break;
	case DEL_INT_BOUNDPT:
	    set_cursor(3);
	    break;
	case DEL_INT_BOUND:
	    set_cursor(3);
	    break;
	case DEL_BOUND_PT:
	    set_cursor(3);
	    break;
	case ADD_BOUND1:
	    set_cursor(0);
	    write_mode_str("Click at first boundary node (CCW)");
	    break;
	case ADD_BOUND2:
	    set_cursor(0);
	    write_mode_str("Click at second boundary node (CCW)");
	    break;
	case ADD_BOUND3:
	    set_cursor(0);
	    write_mode_str("Click at location for new boundary node");
	    break;
	case ADD_BUILD1:
	    set_cursor(0);
	    break;
	case ADD_BUILD2:
	    set_cursor(0);
	    break;
	case ADD_BUILD3:
	    set_cursor(0);
	    break;
	case GENFD_1ST:
	    set_cursor(0);
	    write_mode_str("Click at the first point of the base line for the rectangle");
	    break;
	case GENFD_2ND:
	    set_cursor(0);
	    write_mode_str("Click at the end of the base line for the rectangle");
	    break;
	case GENFD_3RD:
	    set_cursor(0);
	    write_mode_str("Click at the point off the base line to define the rectangle");
	    break;
	case ADD_NODE_BOUND1:
	    set_cursor(0);
	    break;
	case ADD_NODE_BOUND2:
	    set_cursor(0);
	    break;
	case ADD_NODE_BOUND3:
	    set_cursor(0);
	    break;
	case PLACE_BUILD:
	    set_cursor(0);
	    write_mode_str("Click at a point for the new build point location");
	    break;
	case PLACE_BUILD_ARC:
	    break;
	case DELETE_BUILD:
	    set_cursor(3);
	    write_mode_str("Click near a build point to delete (accept with the middle button, cancel with the right)");
	    break;
	case MOVE_BUILD1ST:
	    set_cursor(0);
	    write_mode_str("Click on a build point to move");
	    break;
	case MOVE_BUILD2ND:
	    write_mode_str("Click at the new location of the build point");
	    break;
	case ADD_NODE:
	    set_cursor(0);
	    write_mode_str("Click at a point for the new node location");
	    break;
	case GET_NEAREST_NODE:
	    write_mode_str("Click near a node");
	    set_cursor(1);
	    break;
	case GET_NEAREST_ELEMENT:
	    write_mode_str("Click near an element");
	    set_cursor(1);
	    break;
	case GET_ELEMENT:
	    write_mode_str("Click in the interior of an element");
	    set_cursor(1);
	    break;
	case GET_GRID_DEPTH:
	    write_mode_str("Click at a point in the interior of the edit grid");
	    set_cursor(0);
	    break;
	case GET_BACK_DEPTH:
	    write_mode_str("Click at a point in the interior of the background grid");
	    set_cursor(0);
	    break;
	case GET_NEAREST_BUILDPT:
	    write_mode_str("Click near a build point");
	    set_cursor(1);
	    break;
	case GET_DEPTH_ALL:
	    write_mode_str("Click near a point");
	    set_cursor(1);
	    break;
	case GET_ALL:
	    write_mode_str("Click at a point in the interior of the edit grid");
	    set_cursor(1);
	    break;
	case MOVE_NODE1ST:
	    set_cursor(0);
	    write_mode_str("Click on a node to move");
	    break;
	case MOVE_NODE2ND:
	    set_cursor(0);
	    write_mode_str("Click at a point for the new position of the node");
	    break;
	case SWAPLINE_1ST:
	    set_cursor(0);
	    write_mode_str("Click on one element with shared line to swap");
	    break;
	case SWAPLINE_2ND:
	    set_cursor(0);
	    write_mode_str("Click on adjacent element with shared line to swap");
	    break;
	case CUT_GRID1ST:
	    set_cursor(0);
	    break;
	case CUT_GRID2ND:
	    set_cursor(0);
	    break;
	case ADD_ELEMENT1:
	    set_cursor(0);
	    write_mode_str("Click on the first node (add nodes in a counterclockwise fashion)");
	    break;
	case ADD_ELEMENT2:
	    set_cursor(0);
	    write_mode_str("Click on the second node");
	    break;
	case ADD_ELEMENT3:
	    set_cursor(0);
	    write_mode_str("Click on the third node");
	    break;
	case DELETE_ELEMENT:
	    set_cursor(3);
	    write_mode_str("Click near the element to delete (accept with middle button, cancel with right)");
	    break;
	case DELETE_ELEMENTS:
	    set_cursor(3);
	    write_mode_str("Click in the interior of the element to delete (accept with middle button, cancel with right)");
	    break;
	case SPLIT_ELEMENT3:
	    set_cursor(0);
	    write_mode_str("Click on elements to split in 3");
	    break;
	case SPLIT_ELEMENT4:
	    set_cursor(0);
	    write_mode_str("Click on elements to split in 4");
	    break;
	case SPLIT_ELEMENTS3:
	    set_cursor(0);
	    break;
	case SPLIT_ELEMENTS4:
	    set_cursor(0);
	    break;
	case ZOOM_1ST:
	    set_cursor(0);
	    write_mode_str("Click on one corner of a rectangle to enlarge");
	    break;
	case ZOOM_2ND:
	    set_cursor(0);
	    write_mode_str("Click on the opposite corner of rectangle to enlarge");
	    break;
	case COMP_AREA:
	    set_cursor(0);
	    break;
	case COMP_PERIMETER:
	    set_cursor(0);
	    break;
	case SEL_POINT:
	    set_cursor(0);
	    break;
	case PICK_NODE_LIST:
	    write_mode_str("Click on a node to add to the current list of nodes");
	    set_cursor(0);
	    break;
/*
   case PICK_ELEMENT_LIST:
   set_cursor(0);
   break;
 */
	case PICK_ISTOK:
	    set_cursor(0);
	    break;
	case PICK_PROP_NODE:
	    write_mode_str("Click on a node to set the nodal property");
	    set_cursor(0);
	    break;
	case PICK_PROP_ELEM:
	    write_mode_str("Click on an element to set the elemental property");
	    set_cursor(0);
	    break;
	case QUERY_PROP_ELEM:
	    write_mode_str("Click on an element to get the elemental property");
	    set_cursor(0);
	    break;
	case PLACE_GRIDBATH_LEGEND:
	    write_mode_str("Click at the location of the legend");
	    set_cursor(0);
	    break;
	case PLACE_BACKBATH_LEGEND:
	    write_mode_str("Click at the location of the legend");
	    set_cursor(0);
	    break;
	case PLACE_CONC_LEGEND:
	    write_mode_str("Click at the location of the legend");
	    set_cursor(0);
	    break;
	case DEL_OBJECT:
	    set_cursor(3);
	    write_mode_str("Delete object");
	    break;
	case MOVE_OBJECT_1ST:
	    set_cursor(4);
	    write_mode_str("Pick object");
	    break;
	case MOVE_OBJECT_2ND:
	    write_mode_str("Place object");
	    break;
	case MAKE_BOX_1ST:
	    set_cursor(0);
	    write_mode_str("First corner of box");
	    break;
	case MAKE_BOX_2ND:
	    write_mode_str("Second corner of box");
	    break;
	case STR_LOC1ST:
	    set_cursor(0);
	    write_mode_str("Pick start of text line");
	    break;
	case STR_LOC2ND:
	    write_mode_str("Pick end of text line");
	    break;
	case MAKE_LINE_1ST:
	    set_cursor(0);
	    write_mode_str("Pick beginning of line");
	    break;
	case MAKE_LINE_2ND:
	    write_mode_str("Pick end of line");
	    break;
	case STR_EDIT:
	    set_cursor(2);
	    write_mode_str("Edit string");
	    break;
	case STR_LOC:
	    set_cursor(2);
	    write_mode_str("Pick beginning of text");
	    break;
	case EXTRACT_GRID:
	    set_cursor(0);
	    write_mode_str("Define a region of the edit grid to extract");
	    break;
	case PLACE_VSCALE:
	    set_cursor(0);
	    write_mode_str("Click at a point for the gradient scale");
	    break;
	case ADD_QUAD1:
	    set_cursor(0);
	    write_mode_str("Click near the 1st node of the quadrangular element");
	    break;
	case ADD_QUAD2:
	    set_cursor(0);
	    write_mode_str("Click near the 2nd node of the quadrangular element");
	    break;
	case ADD_QUAD3:
	    set_cursor(0);
	    write_mode_str("Click near the 3rd node of the quadrangular element");
	    break;
	case ADD_QUAD4:
	    set_cursor(0);
	    write_mode_str("Click near the 4th node of the quadrangular element");
	    break;
	case CONV_4TO3:
	    set_cursor(0);
	    write_mode_str("Click inside the quadrangle to convert to 2 triangles");
	    break;
	case CONV_3TO4_1ST:
	    set_cursor(0);
	    write_mode_str("Click inside the first triangle");
	    break;
	case CONV_3TO4_2ND:
	    set_cursor(0);
	    write_mode_str("Click inside the second triangle (must share an edge with the first)");
	    break;
	case EDIT_NODE:
	    set_cursor(0);
	    write_mode_str("Click near the node to edit");
	    break;
	case CONV_NODE2COL:
	    set_cursor(0);
	    write_mode_str("Click near the node to convert to 3 colocated nodes");
	    break;
	case ADD_NODELIST:
	case ADD_BCNODELIST:
	    set_cursor(0);
	    write_mode_str("Click near node to add to the current list of nodes");
	    break;
	case ADD_BCNODELIST1ST:
	case ADD_NODELIST1ST:
	    set_cursor(0);
	    write_mode_str("Click near the starting node");
	    break;
	case ADD_BCNODELIST2ND:
	case ADD_NODELIST2ND:
	    set_cursor(0);
	    write_mode_str("Click near the ending node");
	    break;
	case DELETE_NODE:
	    set_cursor(3);
	    write_mode_str("Click near the node to delete");
	    break;
	case PLACE_ISOLINE_LABEL:
	    set_cursor(0);
	    write_mode_str("Click at a point in the grid to place a label");
	    break;
	case PLACE_ISOLINES_LEGEND:
	    set_cursor(0);
	    write_mode_str("Click on the location of the top of the legend");
	    break;
	case VIEW_1ST:
	    set_cursor(0);
	    write_mode_str("Pick first corner of viewport");
	    break;
	case VIEW_2ND:
	    write_mode_str("Pick second corner of viewport");
	    break;
	case OPENB1ST:
	    set_cursor(0);
	    write_mode_str("Select starting open boundary node");
	    break;
	case OPENB2ND:
	    set_cursor(0);
	    write_mode_str("Select ending open boundary node");
	    break;
	case LANDB1ST:
	    set_cursor(0);
	    write_mode_str("Select starting land boundary node");
	    break;
	case LANDB2ND:
	    set_cursor(0);
	    write_mode_str("Select ending land boundary node");
	    break;
	default:
	    break;
	}
    }
}

/*
 * draw a cursor for text writing
 * TODO: fix the rotation problems (cursor doesn't track)
 */
void update_text_cursor(char *s, int x, int y)
{
    int hgt, tx, xtx, ytx, xhgt, yhgt;

    hgt = stringextenty(charsize * xlibcharsize, "N") / 2;
    tx = stringextentx(charsize * xlibcharsize, s);
    xtx = (int) tx *co;
    ytx = (int) tx *si;

    xhgt = (int) -hgt * si;
    yhgt = (int) hgt *co;

/*    select_line(x + tx, win_h - y + hgt, x + tx, win_h - y - hgt); */
    select_line(x + xtx + xhgt, win_h - (y + ytx + yhgt), x + xtx - xhgt, win_h - (y + ytx - yhgt));
}

/*
 * update string drawn on the canvas
 */
void do_text_string(int op, int c)
{
    char stmp[2];

    drawx = strx;
    drawy = stry;

    update_text_cursor(tmpstr, drawx, drawy);
    set_write_mode(0);
    dispstrxlib(drawx, drawy, string_rot, tmpstr, justflag, 0);
    switch (op) {
    case 0:
	if (strlen(tmpstr) > 0) {
	    tmpstr[strlen(tmpstr) - 1] = 0;
	}
	break;
    case 1:
	sprintf(stmp, "%c", c);
	strcat(tmpstr, stmp);
	break;
    case 2:
	break;
    }
    set_write_mode(1);
    dispstrxlib(drawx, drawy, string_rot, tmpstr, justflag, 0);
    update_text_cursor(tmpstr, drawx, drawy);
}



/*
 * canvas event proc
 */
void my_proc(Widget w, caddr_t data, XEvent * event)
{
    static int x, y, boxno, lineno, ox, oy;
    static double wx1, wx2, wy1, wy2;
    static double wx, wy, dx, dy;
    static int ty, no, c;
    static KeySym keys;
    static XComposeStatus compose;
    extern Widget canvas;
    double xconv(), yconv();

/*
 * hot keys
 */
    x = event->xmotion.x;
    y = event->xmotion.y;
    switch (event->type) {
    case KeyPress:
	buf[0] = 0;
	XLookupString((XKeyEvent *) event, buf, 1, &keys, &compose);
	switch (c = buf[0]) {
	case 1:		/* ^A */
	    break;
	case 2:		/* ^B */
	    break;
	case 3:		/* ^C */
	    set_action(0);
	    set_action(CENTER_WORLD);
	    break;
	case 4:		/* ^D */
	    break;
	case 5:		/* ^E */
	    break;
	case 6:		/* ^F */
	    break;
	case 7:		/* ^G */
	    break;
	    /* stay off 8 (^H) - needed by text routines */
	case 12:		/* ^L */
	    break;
	case 14:		/* ^N */
	    break;
	case 16:		/* ^P */
	    break;
	case 18:		/* ^R */
	    break;
	case 19:		/* ^S */
	    break;
	case 20:		/* ^T */
	    break;
	case 22:		/* ^V */
	    break;
	case 23:		/* ^W */
	    set_action(0);
	    set_action(MOVE_WORLD1ST);
	    break;
	case 24:		/* ^X */
	    break;
	case 26:		/* ^Z */
	    set_action(0);
	    set_action(ZOOM_1ST);
	    break;
	case 27:		/* ESC cancel any operation */
	    set_action(0);
	    break;
	case 8:
	case 127:
	    if (action_flag == STR_LOC) {
		do_text_string(0, 0);
	    }
	    break;
	case '\r':
	case '\n':
	    if (action_flag == STR_LOC) {
		int itmp;

		update_text_cursor(tmpstr, drawx, drawy);
		if (tmpstr[0]) {
		    device2world(drawx, win_h - drawy, &wx, &wy);
		    define_string(tmpstr, wx, wy);
		}
		itmp = (int) (1.25 * stringextenty(charsize * xlibcharsize, "Ny"));
		strx = strx + si * itmp;
		stry = stry - co * itmp;
		tmpstr[0] = 0;
		update_text_cursor(tmpstr, strx, stry);
	    }
	    break;
	default:
	    if (action_flag == STR_LOC) {
		if (c >= 32 && c < 128) {
		    do_text_string(1, c);
		}
	    }
	    break;
	}
	break;
    case EnterNotify:
	set_window(w);
	defineworld(xg1, yg1, xg2, yg2, 0, 0);
	viewport(0.0, 0.0, 1.0, 1.0);
	break;
    case LeaveNotify:
	break;
    case ButtonRelease:
	switch (event->xbutton.button) {
	case Button1:
	    break;
	case Button2:
	    break;
	case Button3:
	    break;
	}
	break;
    case ButtonPress:
	switch (event->xbutton.button) {

	case Button3:
	    switch (action_flag) {
	    case COMP_AREA:
	    case COMP_PERIMETER:
		if (nipoly >= 3) {
		    int i;
		    for (i = 0; i < nipoly; i++) {
		    }
		}
		break;
	    case EXTRACT_GRID:
	    case DEF_EXT_BOUND:
	    case DEF_INT_BOUND:
		npoly = 0;
		nipoly = 0;
		break;
	    case SELECT_SNAP:
	    case DELETE_BUILD:
	    case DELETE_ELEMENTS:
	    case DELETE_ELEMENT:
	    case DELETE_NODE:
		nibuf = 0;
		break;
	    case DEFINE_REGION:
		nregion = 0;
		region_flag = 0;
		break;
	    default:
		break;
	    }
	    npoly = 0;
	    nipoly = 0;
	    nibuf = 0;
	    slice_npts = 0;
	    set_action(0);
	    break;
	case Button1:
	    c = go_locateflag;
	    go_locateflag = TRUE;
	    getpoints(x, y);
	    go_locateflag = c;
	    switch (action_flag) {
	    case DEL_OBJECT:	/* delete a box or a line */
		set_action(0);
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty >= 0) {
		    switch (ty) {
		    case BOX:
			set_write_mode(0);
			draw_box(-2, no);
			set_write_mode(1);
			kill_box(no);
			break;
		    case LINE:
			set_write_mode(0);
			draw_line(-2, no);
			set_write_mode(1);
			kill_line(no);
			break;
		    case STRING:
			set_write_mode(0);
			draw_string(-2, no);
			set_write_mode(1);
			kill_string(no);
			break;
		    }
		}
		break;

/*
 * select a box or a line to move
 */
	    case MOVE_OBJECT_1ST:
		set_action(MOVE_OBJECT_2ND);
		device2world(x, y, &wx, &wy);
/* ZZZZ */
		find_item(cg, wx, wy, &ty, &no);
		if (ty < 0) {
		    set_action(0);
		} else {
		    switch (ty) {
		    case BOX:
			if (boxes[no].loctype == VIEW) {
			    sx = (int) (win_w * boxes[no].x1);
			    sy = (int) (win_h - win_h * boxes[no].y1);
			    xs = (int) (win_w * boxes[no].x2);
			    ys = (int) (win_h - win_h * boxes[no].y2);
			} else {
			    world2deviceabs(boxes[no].x1, boxes[no].y1, &sx, &sy);
			    world2deviceabs(boxes[no].x2, boxes[no].y2, &xs, &ys);
			}
			select_region(sx, sy, xs, ys);
			mbox_flag = 1;
			break;
		    case LINE:
			if (lines[no].loctype == VIEW) {
			    sx = (int) (win_w * lines[no].x1);
			    sy = (int) (win_h - win_h * lines[no].y1);
			    xs = (int) (win_w * lines[no].x2);
			    ys = (int) (win_h - win_h * lines[no].y2);
			} else {
			    world2deviceabs(lines[no].x1, lines[no].y1, &sx, &sy);
			    world2deviceabs(lines[no].x2, lines[no].y2, &xs, &ys);
			}
			select_line(sx, sy, xs, ys);
			mline_flag = 1;
			break;
		    case STRING:
			xs = stringextentx(charsize * xlibcharsize, pstr[no].s);
			ys = stringextenty(charsize * xlibcharsize, pstr[no].s);
			if (pstr[no].loctype == VIEW) {
			    sx = (int) (win_w * pstr[no].x);
			    sy = (int) (win_h - win_h * pstr[no].y);
			} else {
			    world2device(pstr[no].x, pstr[no].y, &sx, &sy);
			}
			xs = sx + xs;
			ys = sy + ys;
			mbox_flag = 1;
			select_region(sx, sy, xs, ys);
			break;
		    }
		}

		break;
/*
 * box has been selected and new position found
 */
	    case MOVE_OBJECT_2ND:
		dx = sx - x;
		dy = sy - y;

		set_action(0);
		sx = x;
		sy = y;
		xs = xs - dx;
		ys = ys - dy;
		device2world(sx, sy, &wx1, &wy1);
		device2world(xs, ys, &wx2, &wy2);
		switch (ty) {
		case BOX:
		    set_write_mode(0);
		    draw_box(-2, no);
		    if (boxes[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			boxes[no].gno = cg;
		    }
		    boxes[no].x1 = wx1;
		    boxes[no].x2 = wx2;
		    boxes[no].y1 = wy1;
		    boxes[no].y2 = wy2;
		    set_write_mode(1);
		    draw_box(-2, no);
		    break;
		case LINE:
		    set_write_mode(0);
		    draw_line(-2, no);
		    if (lines[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			lines[no].gno = cg;
		    }
		    lines[no].x1 = wx1;
		    lines[no].x2 = wx2;
		    lines[no].y1 = wy1;
		    lines[no].y2 = wy2;
		    set_write_mode(1);
		    draw_line(-2, no);
		    break;
		case STRING:
		    set_write_mode(0);
		    draw_string(-2, no);
		    if (pstr[no].loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
		    } else {
			pstr[no].gno = cg;
		    }
		    pstr[no].x = wx1;
		    pstr[no].y = wy1;
		    set_write_mode(1);
		    draw_string(-2, no);
		    break;
		}
		set_action(MOVE_OBJECT_1ST);
		break;
/*
 * select a box or a line to move
 */
	    case EDIT_OBJECT:
		set_action(EDIT_OBJECT);
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if (ty < 0) {
		    set_action(0);
		} else {
		    switch (ty) {
		    case BOX:
			break;
		    case LINE:
			break;
		    case STRING:
			break;
		    }
		}

		break;
/*
 * make a new box, select first corner
 */
	    case MAKE_BOX_1ST:
		set_action(MAKE_BOX_2ND);
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(sx, sy, x, y);
		break;
/*
 * make a new box, select opposite corner
 */
	    case MAKE_BOX_2ND:
		set_action(0);
		if ((boxno = next_box()) >= 0) {
		    device2world(sx, sy, &wx1, &wy1);
		    device2world(x, y, &wx2, &wy2);
		    if (box_loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			boxes[no].gno = cg;
		    }
		    boxes[boxno].loctype = box_loctype;
		    boxes[boxno].x1 = wx1;
		    boxes[boxno].x2 = wx2;
		    boxes[boxno].y1 = wy1;
		    boxes[boxno].y2 = wy2;
		    boxes[boxno].color = box_color;
		    boxes[boxno].linew = box_linew;
		    boxes[boxno].lines = box_lines;
		    boxes[boxno].fill = box_fill;
		    boxes[boxno].fillcolor = box_fillcolor;
		    boxes[boxno].fillpattern = box_fillpat;
		    draw_box(-2, boxno);
		}
		break;
/*
 * locate angled string
 */
	    case STR_LOC1ST:
		set_action(STR_LOC2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		break;
	    case STR_LOC2ND:
		device2world(sx, sy, &wx1, &wy1);
		device2world(x, y, &wx2, &wy2);
		wx1 = xconv(wx1);
		wy1 = yconv(wy1);
		wx2 = xconv(wx2);
		wy2 = yconv(wy2);
		string_rot = (int) ((atan2((wy2 - wy1) * win_h, (wx2 - wx1) * win_w) * 180.0 / M_PI) + 360.0) % 360;
		updatestrings();
		strings_loc_proc();
		break;
/*
 * make a new line, select start point
 */
	    case MAKE_LINE_1ST:
		sx = x;
		sy = y;
		set_action(MAKE_LINE_2ND);
		rubber_flag = 1;
		select_line(sx, sy, x, y);
		break;
/*
 * make a new line, select end point
 */
	    case MAKE_LINE_2ND:
		set_action(0);
		if ((lineno = next_line()) >= 0) {
		    device2world(sx, sy, &wx1, &wy1);
		    device2world(x, y, &wx2, &wy2);
		    if (line_loctype == VIEW) {
			wx1 = xconv(wx1);
			wy1 = yconv(wy1);
			wx2 = xconv(wx2);
			wy2 = yconv(wy2);
		    } else {
			lines[no].gno = cg;
		    }
		    lines[lineno].loctype = line_loctype;
		    lines[lineno].x1 = wx1;
		    lines[lineno].x2 = wx2;
		    lines[lineno].y1 = wy1;
		    lines[lineno].y2 = wy2;
		    lines[lineno].color = line_color;
		    lines[lineno].lines = line_lines;
		    lines[lineno].linew = line_linew;
		    lines[lineno].arrow = line_arrow;
		    lines[lineno].asize = line_asize;
		    lines[lineno].atype = line_atype;
		    draw_line(-2, lineno);
		}
		break;
/*
 * Edit an existing string
 */
	    case STR_EDIT:
		device2world(x, y, &wx, &wy);
		find_item(cg, wx, wy, &ty, &no);
		if ((ty >= 0) && (ty == STRING)) {
		    int ilenx, ileny;

		    wx1 = pstr[no].x;
		    wy1 = pstr[no].y;
		    if (pstr[no].loctype == VIEW) {	/* in viewport coords */
			view2world(wx1, wy1, &wx2, &wy2);
			wx1 = wx2;
			wy1 = wy2;
		    }
		    world2device(wx1, wy1, &strx, &stry);
		    drawx = strx;
		    drawy = stry;
		    strcpy(tmpstr, pstr[no].s);
		    setcharsize(pstr[no].charsize);
		    setfont(pstr[no].font);
		    setcolor(pstr[no].color);
		    string_just = pstr[no].just;
		    justflag = string_just;
		    string_size = pstr[no].charsize;
		    string_font = pstr[no].font;
		    string_color = pstr[no].color;
		    string_linew = pstr[no].linew;
		    string_rot = pstr[no].rot;
		    string_loctype = pstr[no].loctype;
		    updatestrings();
		    kill_string(no);
		    si = sin(M_PI / 180.0 * string_rot) * ((double) win_w) / ((double) win_h);
		    co = cos(M_PI / 180.0 * string_rot);

		    ilenx = stringextentx(charsize * xlibcharsize, tmpstr);
		    ileny = stringextenty(charsize * xlibcharsize, tmpstr);

		    switch (justflag) {
		    case 1:
			strx = drawx + co * ilenx - si * ileny;
			stry = drawy + si * ilenx + co * ileny;
			break;
		    case 2:
			strx = drawx + (co * ilenx - si * ileny) / 2;
			stry = drawy + (si * ilenx + co * ileny) / 2;
			break;
		    }
		    update_text_cursor(tmpstr, drawx, drawy);
		    do_text_string(2, 0);
		    action_flag = STR_LOC;
		} else {
		    set_action(0);
		}
		break;
/*
 * locate a string on the canvas
 */
	    case STR_LOC:
		if (tmpstr[0]) {
		    device2world(strx, win_h - stry, &wx, &wy);
		    define_string(tmpstr, wx, wy);
		}
		strx = x;
		stry = win_h - y;
		drawx = strx;
		drawy = stry;
		tmpstr[0] = 0;
		define_string_defaults();
		justflag = string_just;
		setcharsize(string_size);
		xlibsetfont(string_font);
		xlibsetcolor(string_color);
		xlibsetlinewidth(string_linew);
		si = sin(M_PI / 180.0 * string_rot) * ((double) win_w) / ((double) win_h);
		co = cos(M_PI / 180.0 * string_rot);
		update_text_cursor(tmpstr, strx, stry);
		break;
/*
 * Compute the area of a polygon
 */
	    case COMP_AREA:
		{
		    double area, comp_area();

		    get_world(x, y, &polyx[npoly], &polyy[npoly]);
		    ipolyx[npoly] = x;
		    ipolyy[npoly] = y;
		    npoly++;
		    nipoly++;
		    if (npoly <= 2) {
			area = 0.0;
		    } else {
			area = comp_area(npoly, polyx, polyy);
		    }
		    sprintf(buf, "Area [%lf]", fabs(area));
		    panel_setmsgstr_value(locate_grid_item, buf);
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(COMP_AREA);
		}
		break;
	    case COMP_PERIMETER:
		{
		    double area, comp_perimeter();

		    get_world(x, y, &polyx[npoly], &polyy[npoly]);
		    ipolyx[npoly] = x;
		    ipolyy[npoly] = y;
		    npoly++;
		    nipoly++;
		    if (npoly <= 1) {
			area = 0.0;
		    } else {
			area = comp_perimeter(npoly, polyx, polyy);
		    }
		    sprintf(buf, "Perimeter [%lf]", fabs(area));
		    panel_setmsgstr_value(locate_grid_item, buf);
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(COMP_PERIMETER);
		}
		break;
/*
 * place an isoline label
 */
	    case PLACE_ISOLINE_LABEL:
		{
		    get_world(x, y, &wx, &wy);
		    find_element(string_grid, wx, wy, &tmpind);
		    if (tmpind >= 0) {
			define_depth_label(string_grid, tmpind, wx, wy);
		    } else {
			errwin("No element found\n");
		    }
		    set_action(PLACE_ISOLINE_LABEL);
		}
		break;
/*
 * Pick prop
 */
	    case PICK_PROP_NODE:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    diamond(grid[curgrid].xord[tmpind], grid[curgrid].yord[tmpind]);
		    flush_pending();
		    prop_set(0, tmpind);
		}
		set_action(PICK_PROP_NODE);
		break;
	    case PICK_PROP_ELEM:
		get_world(x, y, &wx, &wy);
		find_nearest_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    get_center(curgrid, tmpind, &wx, &wy);
		    writestr(wx, wy, 0, 0, "s");
		    flush_pending();
		    prop_set(0, tmpind);
		}
		set_action(PICK_PROP_ELEM);
		break;
	    case QUERY_PROP_ELEM:
		get_world(x, y, &wx, &wy);
		find_nearest_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    double d, prop_get(int, int);
		    d = prop_get(0, tmpind);
		    sprintf(buf, "Element %d: Prop = %.4lf", tmpind + 1, d);
		    panel_setmsgstr_value(locate_grid_item, buf);
		}
		set_action(QUERY_PROP_ELEM);
		break;
	    case PICK_ISTOK:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    get_center(curgrid, tmpind, &wx, &wy);
		    writestr(wx, wy, 0, 0, "w");
		}
		set_action(PICK_ISTOK);
		break;
	    case PICK_NODE_LIST:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    AddToPNodeList(curnl, tmpind, curnlival, curnldval, curnlcolor);
		    DrawPNodeLists(curgrid);
		}
		set_action(PICK_NODE_LIST);
		break;
/*
 * select a reference point for the locator in main_panel
 */
	    case SEL_POINT:
		get_world(x, y, &wx, &wy);
		dsx = wx;
		dsy = wy;
		my_circle(dsx, dsy);
		set_action(0);
		break;
/*
 * Center the world
 */
	    case CENTER_WORLD:
		get_world(x, y, &wx1, &wy1);
		write_mode_str("Click on the center of the new view");
		dx = (xg1 + xg2) * 0.5 - wx1;
		dy = (yg1 + yg2) * 0.5 - wy1;
		xg1 = xg1 - dx;
		xg2 = xg2 - dx;
		yg1 = yg1 - dy;
		yg2 = yg2 - dy;
		set_up_world();
		do_drawgrid();
		set_action(CENTER_WORLD);
		break;
/*
 * Move the world
 */
	    case MOVE_WORLD1ST:
		get_world(x, y, &wx1, &wy1);
		ox = x;
		oy = y;
		set_action(MOVE_WORLD2ND);
		write_mode_str("Move to the new view");
		world2deviceabs(xg1, yg1, &sx, &sy);
		world2deviceabs(xg2, yg2, &xs, &ys);
		select_region(sx, sy, xs, ys);
		break;
	    case MOVE_WORLD2ND:
		set_action(0);
		get_world(x, y, &wx2, &wy2);
		dx = wx2 - wx1;
		dy = wy2 - wy1;
		xg1 = xg1 - dx;
		xg2 = xg2 - dx;
		yg1 = yg1 - dy;
		yg2 = yg2 - dy;
		set_up_world();
		do_drawgrid();
		break;
/*
 * set one corner of zoom
 */
	    case ZOOM_1ST:
		set_action(ZOOM_2ND);
		write_mode_str("Click on the opposite corner of rectangle to enlarge");
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(x, y, x, y);
		break;
/*
 * set opposing corner of zoom
 */
	    case ZOOM_2ND:
		set_action(0);
		select_region(sx, sy, old_x, old_y);
		get_world(sx, sy, &wx1, &wy1);
		get_world(old_x, old_y, &wx2, &wy2);
		if (wx1 > wx2)
		    fswap(&wx1, &wx2);
		if (wy1 > wy2)
		    fswap(&wy1, &wy2);
		xg1 = wx1;
		xg2 = wx2;
		yg1 = wy1;
		yg2 = wy2;
		set_up_world();
		do_drawgrid();
		break;
/*
 * set one corner of viewport
 */
	    case VIEW_1ST:
		set_action(VIEW_2ND);
		rectflag = 1;
		sx = x;
		sy = y;
		select_region(x, y, x, y);
		break;
/*
 * set opposing corner of viewport
 */
	    case VIEW_2ND:
		{
		    double vx1, vx2, vy1, vy2;
		    extern double xv1, xv2, yv1, yv2;

		    set_action(0);
		    select_region(sx, sy, old_x, old_y);
		    if (sx == old_x || sy == old_y) {
			errwin("Viewport size incorrect, not changed");
		    } else {
			device2world(sx, sy, &wx1, &wy1);
			device2world(old_x, old_y, &wx2, &wy2);
			world2view(wx1, wy1, &vx1, &vy1);
			world2view(wx2, wy2, &vx2, &vy2);
			if (vx1 > vx2) {
			    fswap(&vx1, &vx2);
			}
			if (vy1 > vy2) {
			    fswap(&vy1, &vy2);
			}
			xv1 = vx1;
			yv1 = vy1;
			xv2 = vx2;
			yv2 = vy2;
			set_up_world();
			update_view();
			do_drawgrid();
		    }
		}
		break;
/*
 * The following are specific for gredit
 */
/*
 * Define an exterior boundary
 */
	    case EXTRACT_GRID:
		get_world(x, y, &polyx[npoly], &polyy[npoly]);
		ipolyx[npoly] = x;
		ipolyy[npoly] = y;
		npoly++;
		nipoly++;
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(EXTRACT_GRID);
		break;
	    case DEF_EXT_BOUND:
		get_world(x, y, &polyx[npoly], &polyy[npoly]);
		ipolyx[npoly] = x;
		ipolyy[npoly] = y;
		npoly++;
		nipoly++;
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(DEF_EXT_BOUND);
		break;
/*
 * Define an interior boundary
 */
	    case DEF_INT_BOUND:
		get_world(x, y, &polyx[npoly], &polyy[npoly]);
		ipolyx[npoly] = x;
		ipolyy[npoly] = y;
		npoly++;
		nipoly++;
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(DEF_INT_BOUND);
		break;
/*
 * Move an exterior boundary point
 */
	    case MOVE_EXT_BOUND1ST:
		get_world(x, y, &wx, &wy);
		find_external_boundary_point(curgrid, wx, wy, &bnof, &indf);
		if (indf < 0) {
		    set_action(0);
		    return;
		}
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(MOVE_EXT_BOUND2ND);
		break;
	    case MOVE_EXT_BOUND2ND:
		set_action(0);
		setcolor(0);
		diamond(wx, wy);
		setcolor(1);
		get_world(x, y, &wx, &wy);
		move_boundary_node(curgrid, bnof, indf, wx, wy);
		diamond(wx, wy);
		set_action(MOVE_EXT_BOUND1ST);
		break;
/*
 * Move an interior boundary point
 */
	    case MOVE_INT_BOUND1ST:
		get_world(x, y, &wx, &wy);
		find_internal_boundary_point(curgrid, wx, wy, &bnof, &indf);
		if (indf < 0) {
		    set_action(0);
		    return;
		}
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(MOVE_INT_BOUND2ND);
		break;
	    case MOVE_INT_BOUND2ND:
		set_action(0);
		setcolor(0);
		diamond(wx, wy);
		setcolor(1);
		get_world(x, y, &wx, &wy);
		move_boundary_node(curgrid, bnof, indf, wx, wy);
		diamond(wx, wy);
		set_action(MOVE_INT_BOUND1ST);
		break;
/*
 * Delete an exterior boundary point
 */
	    case DEL_EXT_BOUNDPT:
		get_world(x, y, &wx, &wy);
		break;
/*
 * Delete an interior boundary point
 */
	    case DEL_INT_BOUNDPT:
		get_world(x, y, &wx, &wy);
		break;
/*
 * Delete an interior boundary
 */
	    case DEL_INT_BOUND:
		get_world(x, y, &wx, &wy);
		find_internal_boundary_point(curgrid, wx, wy, &ind1, &ind2);
		if (ind1 >= 0) {
		    if (ind2 >= 0) {
			disassociate_grid(curgrid, ind1);
		    }
		}
		break;
/*
 * Add a point to a boundary
 */
	    case ADD_BOUND1:
		/* select 1st boundary node */
		get_world(x, y, &wx, &wy);
		find_boundary_point(curgrid, wx, wy, &bnd1, &ind1);
		if (bnd1 < 0 || ind1 < 0) {
		    set_action(0);
		} else {
		    my_circlefilled(wx, wy);
		    flush_pending();
		    set_action(ADD_BOUND2);
		}
		break;
	    case ADD_BOUND2:
		/* select 2nd boundary node */
		get_world(x, y, &wx, &wy);
		find_boundary_point(curgrid, wx, wy, &bnd2, &ind2);
		if (bnd1 != bnd2 || bnd2 < 0 || ind2 < 0) {
		    set_action(0);
		} else {
		    my_circlefilled(wx, wy);
		    flush_pending();
		    set_action(ADD_BOUND3);
		}
		break;
	    case ADD_BOUND3:
		/* select point to insert */
		get_world(x, y, &wx, &wy);
		my_circlefilled(wx, wy);
		flush_pending();
		add_boundary_node(curgrid, bnd1, ind1, ind2, wx, wy);
		set_action(ADD_BOUND1);
		break;
/*
 * Spread points
 */
	    case ADD_BUILD1:
		get_world(x, y, &wx1, &wy1);
		switch (spreadflag) {
		case 0:
		    rectflag = 1;
		    sx = x;
		    sy = y;
		    select_region(x, y, x, y);
		    break;
		case 1:
		    rectflag = 1;
		    sx = x;
		    sy = y;
		    select_region(x, y, x, y);
		    break;
		case 2:
		    rectflag = 1;
		    sx = x;
		    sy = y;
		    select_region(x, y, x, y);
		    break;
		case 3:
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    break;
		case 4:
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    break;
		}
		set_action(ADD_BUILD2);
		break;
	    case ADD_BUILD2:
		set_action(0);
		switch (spreadflag) {
		case 0:
		    select_region(sx, sy, old_x, old_y);
		    get_world(sx, sy, &wx1, &wy1);
		    get_world(old_x, old_y, &wx2, &wy2);
		    if (wx1 > wx2)
			fswap(&wx1, &wx2);
		    if (wy1 > wy2)
			fswap(&wy1, &wy2);
		    spread_rectangular(curbuild, wx1, wy1, wx2, wy2);
		    break;
		case 1:
		    select_region(sx, sy, old_x, old_y);
		    get_world(sx, sy, &wx1, &wy1);
		    get_world(old_x, old_y, &wx2, &wy2);
		    if (wx1 > wx2)
			fswap(&wx1, &wx2);
		    if (wy1 > wy2)
			fswap(&wy1, &wy2);
		    spread_rectangular_offset(curbuild, wx1, wy1, wx2, wy2);
		    break;
		case 2:
		    select_region(sx, sy, old_x, old_y);
		    get_world(sx, sy, &wx1, &wy1);
		    get_world(old_x, old_y, &wx2, &wy2);
		    if (wx1 > wx2)
			fswap(&wx1, &wx2);
		    if (wy1 > wy2)
			fswap(&wy1, &wy2);
		    spread_random(curbuild, wx1, wy1, wx2, wy2);
		    break;
		case 3:
		    get_world(sx, sy, &wx1, &wy1);
		    get_world(old_x, old_y, &wx2, &wy2);
		    setlinewidth(3);
		    my_move2(wx1, wy1);
		    my_draw2(wx2, wy2);
		    setlinewidth(1);
		    set_action(ADD_BUILD3);
		    break;
		case 4:
		    get_world(sx, sy, &wx1, &wy1);
		    get_world(old_x, old_y, &wx2, &wy2);
		    setlinewidth(3);
		    my_move2(wx1, wy1);
		    my_draw2(wx2, wy2);
		    setlinewidth(1);
		    set_action(ADD_BUILD3);
		    break;
		}
		break;
	    case ADD_BUILD3:
		get_world(x, y, &wx3, &wy3);
		switch (spreadflag) {
		case 3:
		    spread_rotated_rect(curbuild, wx1, wy1, wx2, wy2, wx3, wy3);
		    break;
		case 4:
		    spread_rotated_rect_offset(curbuild, wx1, wy1, wx2, wy2, wx3, wy3);
		    break;
		}
		set_action(0);
		break;
	    case GENFD_1ST:
		sx = x;
		sy = y;
		rubber_flag = 1;
		select_line(x, y, x, y);
		set_action(GENFD_2ND);
		break;
	    case GENFD_2ND:
		set_action(0);
		get_world(sx, sy, &wx1, &wy1);
		get_world(old_x, old_y, &wx2, &wy2);
		setlinewidth(3);
		my_move2(wx1, wy1);
		my_draw2(wx2, wy2);
		setlinewidth(1);
		set_action(GENFD_3RD);
		sx = old_x;
		sy = old_y;
		rubber_flag = 1;
		select_line(old_x, old_y, old_x, old_y);
		break;
	    case GENFD_3RD:
		get_world(old_x, old_y, &wx3, &wy3);
		set_action(0);
		setlinewidth(3);
		my_move2(wx2, wy2);
		my_draw2(wx3, wy3);
		setlinewidth(1);
		genfd_grid(wx1, wy1, wx2, wy2, wx3, wy3);
		break;
/*
 * Add a nodal point (in the current grid) to a boundary
 */
	    case ADD_NODE_BOUND1:
		/* select 1st boundary node */
		get_world(x, y, &wx, &wy);
		break;
	    case ADD_NODE_BOUND2:
		/* select 2nd boundary node */
		get_world(x, y, &wx, &wy);
		break;
	    case ADD_NODE_BOUND3:
		/* select nodal point */
		get_world(x, y, &wx, &wy);
		break;
/*
 * Place a single build point
 */
	    case PLACE_BUILD:
		get_world(x, y, &wx1, &wy1);
		if (goodpoint(curgrid, wx1, wy1)) {
		    add_build(curbuild, wx1, wy1, 1.0);
		    setcolor(2);
		    box(wx1, wy1);
		    flush_pending();
		    setcolor(1);
		}
		set_action(PLACE_BUILD);
		break;
/*
 * Delete a single build point
 */
	    case DELETE_BUILD:
		get_world(x, y, &wx, &wy);
		accumulate_nearest_buildpts(curbuild, wx, wy);
		set_action(DELETE_BUILD);
		break;
/*
 * Move a single build point
 */
	    case MOVE_BUILD1ST:
		get_world(x, y, &wx, &wy);
		find_nearest_buildpt(curbuild, wx, wy, &indf);
		if (indf < 0) {
		    set_action(0);
		    return;
		}
		wx = build[curbuild].bx[indf];
		wy = build[curbuild].by[indf];
		get_device(wx, wy, &x, &y);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(MOVE_BUILD2ND);
		break;
	    case MOVE_BUILD2ND:
		setcolor(0);
		box(wx, wy);
		setcolor(1);
		get_world(x, y, &wx, &wy);
		move_buildpt(curbuild, indf, wx, wy);
		set_action(0);
		box(wx, wy);
		flush_pending();
		set_action(MOVE_BUILD1ST);
		break;
/*
 * Place a build point arc
 */
	    case PLACE_BUILD_ARC:
		get_world(x, y, &wx, &wy);
		break;
/*
 * Report on the nearest node to the mouse location
 */
	    case GET_NEAREST_NODE:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind1);
		wx = grid[curgrid].xord[ind1];
		wy = grid[curgrid].yord[ind1];
		if (ind1 >= 0) {
		    sprintf(buf, "Node %d: [%.2lf,%.2lf], depth = %.3lf", ind1 + 1, wx, wy, grid[curgrid].depth[ind1]);
		}
		panel_setmsgstr_value(locate_grid_item, buf);
		set_action(GET_NEAREST_NODE);
		break;
	    case GET_NEAREST_BUILDPT:
		get_world(x, y, &wx, &wy);
		find_nearest_buildpt(0, wx, wy, &ind1);
		if (ind1 >= 0) {
		    wx = build[0].bx[ind1];
		    wy = build[0].by[ind1];
		    sprintf(buf, "Point %d: [%.2lf,%.2lf], depth = %.3lf", ind1 + 1, wx, wy, build[0].db[ind1]);
		    panel_setmsgstr_value(locate_grid_item, buf);
		}
		set_action(GET_NEAREST_BUILDPT);
		break;
	    case GET_DEPTH_ALL:
		{
		    double backd = 0.0, editd = 0.0, buildd = 0.0;
		    int backi, editi, buildi;
		    char tbuf[256];
		    get_world(x, y, &wx, &wy);
		    find_nearest_buildpt(0, wx, wy, &buildi);
		    if (buildi >= 0) {
			wx = build[0].bx[buildi];
			wy = build[0].by[buildi];
			buildd = build[0].db[buildi];
		    }
		    find_element(curgrid, wx, wy, &editi);
		    if (editi >= 0) {
			editd = get_depth_element(curgrid, editi, wx, wy);
		    }
		    find_element(backgrid, wx, wy, &backi);
		    if (backi >= 0) {
			backd = get_depth_element(backgrid, backi, wx, wy);
		    }
		    sprintf(buf, "Point %d, Edit el %d, Back el %d: [%.2lf,%.2lf], %.3lf, %.3lf, %.3lf", buildi + 1, editi + 1, backi + 1, wx, wy, buildd, editd, backd);
		    panel_setmsgstr_value(locate_grid_item, buf);
		    set_action(GET_DEPTH_ALL);
		}
		break;
/*
 * Report on the nearest element to the mouse location
 */
	    case GET_NEAREST_ELEMENT:
		get_world(x, y, &wx, &wy);
		find_nearest_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    double quadqual(int gno, int e);
		    switch (grid[curgrid].icon[ind1].type) {
		    case 3:
		    case 6:
			n0 = grid[curgrid].icon[ind1].nl[0] + 1;
			n1 = grid[curgrid].icon[ind1].nl[1] + 1;
			n2 = grid[curgrid].icon[ind1].nl[2] + 1;
			wx = area(curgrid, ind1);
			if (wx < 0.0) {
			    wy = sqrt(-wx / M_PI);
			} else {
			    wy = sqrt(wx / M_PI);
			}
			sprintf(buf, "Element %d: [%d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, wx, wy);
			break;
		    case 4:
		    case 8:
			n0 = grid[curgrid].icon[ind1].nl[0] + 1;
			n1 = grid[curgrid].icon[ind1].nl[1] + 1;
			n2 = grid[curgrid].icon[ind1].nl[2] + 1;
			n3 = grid[curgrid].icon[ind1].nl[3] + 1;
			sprintf(buf, "Element %d: [%d %d %d %d] Q = %.3lf", ind1 + 1, n0, n1, n2, n3, quadqual(curgrid, ind1));
			break;
		    }
		} else {
		    sprintf(buf, "No element found");
		}
		get_center(curgrid, ind1, &wx, &wy);
		my_circlefilled(wx, wy);
		panel_setmsgstr_value(locate_grid_item, buf);
		set_action(GET_NEAREST_ELEMENT);
		break;
/*
 * Report on the element containing the current mouse location
 */
	    case GET_ELEMENT:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    double quadqual(int gno, int e);
		    switch (grid[curgrid].icon[ind1].type) {
		    case 3:
		    case 6:
			n0 = grid[curgrid].icon[ind1].nl[0] + 1;
			n1 = grid[curgrid].icon[ind1].nl[1] + 1;
			n2 = grid[curgrid].icon[ind1].nl[2] + 1;
			wx = area(curgrid, ind1);
			if (wx < 0.0) {
			    wy = sqrt(-wx / M_PI);
			} else {
			    wy = sqrt(wx / M_PI);
			}
			sprintf(buf, "Element %d: [%d %d %d], A = %.2lf, ER = %.2lf", ind1 + 1, n0, n1, n2, wx, wy);
			break;
		    case 4:
		    case 8:
			n0 = grid[curgrid].icon[ind1].nl[0] + 1;
			n1 = grid[curgrid].icon[ind1].nl[1] + 1;
			n2 = grid[curgrid].icon[ind1].nl[2] + 1;
			n3 = grid[curgrid].icon[ind1].nl[3] + 1;
			sprintf(buf, "Element %d: [%d %d %d %d] Q = %.3lf", ind1 + 1, n0, n1, n2, n3, quadqual(curgrid, ind1));
			break;
		    }
		} else {
		    sprintf(buf, "No element found");
		}
		panel_setmsgstr_value(locate_grid_item, buf);
		set_action(GET_ELEMENT);
		break;
/*
 * Verbose report on the node and element
 */
	    case GET_ALL:
		get_world(x, y, &wx, &wy);
		break;
/*
 * Get the depth at the mouse location
 */
	    case GET_GRID_DEPTH:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {

		    sprintf(buf, "Element %d: at (%lf, %lf) Depth = %lf", ind1 + 1, wx, wy, get_depth_element(curgrid, ind1, wx, wy));
		    panel_setmsgstr_value(locate_grid_item, buf);
		    set_action(GET_GRID_DEPTH);
		}
		break;
	    case GET_BACK_DEPTH:
		get_world(x, y, &wx, &wy);
		find_element(backgrid, wx, wy, &ind1);
		if (ind1 >= 0) {

		    sprintf(buf, "Element %d: at (%lf, %lf) Depth = %lf", ind1 + 1, wx, wy, get_depth_element(backgrid, ind1, wx, wy));
		    panel_setmsgstr_value(locate_grid_item, buf);
		    set_action(GET_BACK_DEPTH);
		}
		get_world(x, y, &wx, &wy);
		break;
/*
 * Move the node nearest the mouse location
 */
	    case MOVE_NODE1ST:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &tmpind);
		wx = grid[curgrid].xord[tmpind];
		wy = grid[curgrid].yord[tmpind];
		get_device(wx, wy, &sx, &sy);
		rubber_flag = 1;
		select_line(sx, sy, sx, sy);
		set_action(MOVE_NODE2ND);
		break;
	    case MOVE_NODE2ND:
		set_action(0);
		get_world(x, y, &wx, &wy);
		setcolor(0);
		find_assoc_elements(curgrid, tmpind, &nibuf, ibuf);
		draw_elements(curgrid, nibuf, ibuf, redfact);
		setcolor(1);
		grid[curgrid].xord[tmpind] = wx;
		grid[curgrid].yord[tmpind] = wy;
		draw_elements(curgrid, nibuf, ibuf, redfact);
		flush_pending();
		set_action(MOVE_NODE1ST);
		break;
/*
 * Delete the elements containing the mouse pointer
 */
	    case DELETE_ELEMENTS:
		get_world(x, y, &wx, &wy);
		accumulate_elements(curgrid, wx, wy);
		find_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    get_center(curgrid, ind1, &wx, &wy);
		    solidbox(wx, wy);
		    flush_pending();
		}
		set_action(DELETE_ELEMENTS);
		break;
/*
 * Delete the element with center of mass nearest mouse pointer
 */
	    case DELETE_ELEMENT:
		get_world(x, y, &wx, &wy);
		accumulate_nearest_elements(curgrid, wx, wy);
		set_action(DELETE_ELEMENT);
		break;
/*
 * Add a node to the table of nodes
 */
	    case ADD_NODE:
		get_world(x, y, &wx, &wy);
		my_circlefilled(wx, wy);
		my_move2(wx, wy);
		add_node(curgrid, wx, wy, 1.0);
		set_action(ADD_NODE);
		break;
/*
 * Delete the nearest node and all associated elements
 */
	    case DELETE_NODE:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    writestr(wx, wy, 0, 0, "x");
		    delete_node(curgrid, ind1);
		    set_action(DELETE_NODE);
		} else {
		    errwin("No node found!");
		    set_action(0);
		}
		break;
/*
 * Add an element to the table of elements
 */
	    case ADD_ELEMENT1:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind1);
		writestr(wx, wy, 0, 0, "1");
		set_action(ADD_ELEMENT2);
		break;
	    case ADD_ELEMENT2:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind2);
		writestr(wx, wy, 0, 0, "2");
		set_action(ADD_ELEMENT3);
		break;
	    case ADD_ELEMENT3:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind3);
		if (area_from_nodes(curgrid, ind1, ind2, ind3) <= 0.0) {
		    iswap(&ind1, &ind3);
		    writestr(wx, wy, 0, 0, "3 (swapped, area <= 0.0)");
		} else {
		    writestr(wx, wy, 0, 0, "3");
		}
		tmpind = grid[curgrid].nmel;
		add_element(curgrid, ind1, ind2, ind3);
		draw_element(curgrid, tmpind, redfact);
		my_move2(wx, wy);
		set_action(ADD_ELEMENT1);
		break;
/*
 * Add a quadratic element to the table of elements
 */
	    case ADD_QUAD1:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind1);
		writestr(wx, wy, 0, 0, "1");
		set_action(ADD_QUAD2);
		break;
	    case ADD_QUAD2:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind2);
		writestr(wx, wy, 0, 0, "2");
		set_action(ADD_QUAD3);
		break;
	    case ADD_QUAD3:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind3);
		writestr(wx, wy, 0, 0, "3");
		set_action(ADD_QUAD4);
		break;
	    case ADD_QUAD4:
		get_world(x, y, &wx, &wy);
		find_nearest_node(curgrid, wx, wy, &ind4);
/*
   if (area_from_nodes(curgrid, ind1, ind2, ind3) <= 0.0) {
   iswap(&ind1, &ind3);
   writestr(wx, wy, 0, 0, "3 (swapped, area <= 0.0)");
   } else {
   writestr(wx, wy, 0, 0, "3");
   }
 */
		writestr(wx, wy, 0, 0, "4");
		tmpind = grid[curgrid].nmel;
		add_element2(curgrid, 4, ind1, ind2, ind3, ind4);
		draw_element(curgrid, tmpind, redfact);
		my_move2(wx, wy);
		set_action(ADD_QUAD1);
		break;
	    case CONV_4TO3:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    get_center(curgrid, ind1, &wx, &wy);
		    writestr(wx, wy, 0, 0, "q");
		    flush_pending();
		    quad_to_tris(curgrid, ind1);
		    set_action(CONV_4TO3);
		} else {
		    set_action(0);
		}
		break;
	    case CONV_3TO4_1ST:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind1);
		if (ind1 >= 0) {
		    get_center(curgrid, ind1, &wx, &wy);
		    writestr(wx, wy, 0, 0, "1");
		    flush_pending();
		    set_action(CONV_3TO4_2ND);
		} else {
		    set_action(0);
		}
		break;
	    case CONV_3TO4_2ND:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind2);
		if (ind2 >= 0) {
		    get_center(curgrid, ind2, &wx, &wy);
		    writestr(wx, wy, 0, 0, "2");
		    flush_pending();
		    tris_to_quad(curgrid, ind1, ind2);
		    set_action(CONV_3TO4_1ST);
		} else {
		    set_action(0);
		}
		break;
	    case EDIT_NODE:
		break;
	    case CONV_NODE2COL:
		break;
/*
 * Swap the line between two elements
 */
	    case SWAPLINE_1ST:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind1);
		set_action(0);
		if (ind1 >= 0) {
		    get_center(curgrid, ind1, &wx, &wy);
		    writestr(wx, wy, 0, 0, "1");
		    set_action(SWAPLINE_2ND);
		}
		break;
	    case SWAPLINE_2ND:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &ind2);
		set_action(0);
		if (ind2 >= 0) {
		    get_center(curgrid, ind2, &wx, &wy);
		    writestr(wx, wy, 0, 0, "2");
		    setcolor(0);
		    draw_element(curgrid, ind1, redfact);
		    draw_element(curgrid, ind2, redfact);
		    setcolor(1);
		    swapline(curgrid, ind1, ind2);
		    draw_element(curgrid, ind1, redfact);
		    draw_element(curgrid, ind2, redfact);
		    my_move2(wx, wy);
		}
		set_action(SWAPLINE_1ST);
		break;
/*
 * Snap nodes to a line
 */
	    case SELECT_SNAP:
		get_world(x, y, &wx, &wy);
		accumulate_nearest_nodes(curgrid, wx, wy);
		set_action(SELECT_SNAP);
		break;
	    case SNAP_LINE1ST:
		get_world(x, y, &wx1, &wy1);
		set_action(0);
		set_action(SNAP_LINE2ND);
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		break;
	    case SNAP_LINE2ND:
		get_world(x, y, &wx2, &wy2);
		if (wx2 - wx1 == 0.0) {
		    sm = 1.0e37;
		    sb = wx1;
		} else {
		    sm = (wy2 - wy1) / (wx2 - wx1);
		    sb = wy2 - sm * wx2;
		}
		nibuf = 0;
		set_action(0);
		write_mode_str("Select nodes to snap with left button, middle button to register, right button to abort");
		set_action(SELECT_SNAP);
		select_line(sx, sy, x, y);
		break;
/*
 * Split an element into 3 elements
 */
	    case SPLIT_ELEMENT3:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
/*
   get_center(curgrid, tmpind, &wx, &wy);
   writestr(wx, wy, 0, 0, "3");
 */
		    split_elem3(curgrid, tmpind);
		}
		break;
	    case SPLIT_ELEMENTS3:
		break;
/*
 * Split an element into 4 elements
 */
	    case SPLIT_ELEMENT4:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
/*
   get_center(curgrid, tmpind, &wx, &wy);
   writestr(wx, wy, 0, 0, "4");
 */
		    split_elem4(curgrid, tmpind);
		}
		break;
	    case SPLIT_ELEMENTS4:
		break;
/*
 * Cut a grid by a line
 */
	    case CUT_GRID1ST:
		rubber_flag = 1;
		sx = x;
		sy = y;
		select_line(sx, sy, x, y);
		set_action(CUT_GRID2ND);
		break;
	    case CUT_GRID2ND:
		set_action(0);
		get_world(sx, sy, &wx1, &wy1);
		get_world(x, y, &wx2, &wy2);
		cut_grid(curgrid, wx1, wy1, wx2, wy2, cutgrid_type);
		do_drawgrid();
		break;
/*
 * Define a region for operations
 */
	    case DEFINE_REGION:
		{
		    get_world(x, y, &regionx[nregion], &regiony[nregion]);
		    ipolyx[npoly] = x;
		    ipolyy[npoly] = y;
		    nipoly++;
		    nregion++;
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(DEFINE_REGION);
		}
		break;
/*
 * Slice a grid
 */
	    case SLICE_BATH1ST:
		get_world(x, y, &wx1, &wy1);
		find_element(curgrid, wx1, wy1, &tmpind);
		if (tmpind >= 0) {
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(SLICE_BATH2ND);
		} else {
		    set_action(0);
		}
		break;
	    case SLICE_BATH2ND:
		get_world(x, y, &wx2, &wy2);
		find_element(curgrid, wx2, wy2, &tmpind);
		if (tmpind >= 0) {
		    extern int slice_back;
		    get_slice(curgrid, wx1, wy1, wx2, wy2, 1);
		    if (slice_back) {
			get_slice(MAXGRIDS, wx1, wy1, wx2, wy2, 0);
		    }
		} else {
		}
		set_action(0);
		break;
	    case SLICE_GRIDBATH1ST:
		get_world(x, y, &wx1, &wy1);
		find_element(curgrid, wx1, wy1, &tmpind);
		if (tmpind >= 0) {
		    rubber_flag = 1;
		    sx = x;
		    sy = y;
		    select_line(sx, sy, x, y);
		    set_action(SLICE_GRIDBATH2ND);
		} else {
		    set_action(0);
		}
		break;
	    case SLICE_GRIDBATH2ND:
		get_world(x, y, &wx2, &wy2);
		find_element(curgrid, wx2, wy2, &tmpind);
		if (tmpind >= 0) {
		    set_action(0);
		    my_move2(wx1, wy1);
		    my_draw2(wx2, wy2);
		} else {
		    set_action(0);
		}
		set_action(SLICE_GRIDBATH3RD);
		break;
	    case SLICE_GRIDBATH3RD:
		get_world(x, y, &wx, &wy);
		find_element(curgrid, wx, wy, &tmpind);
		if (tmpind >= 0) {
		    box(wx, wy);
		    flush_pending();
		    slice_bygrid(curgrid, wx1, wy1, wx2, wy2, wx, wy);
		} else {
		}
		set_action(0);
		break;

	    case OPENB1ST:
		{
		    int bno = grid[curgrid].boundaries[0];
		    set_action(OPENB2ND);
		    get_world(x, y, &wx1, &wy1);
		    find_nearest_boundary_point2(bno, wx1, wy1, &ind1);
		    if (ind1 >= 0) {
			diamond(boundary[bno].boundx[ind1], boundary[bno].boundy[ind1]);
			flush_pending();
		    } else {
			errwin("No boundary point found");
			set_action(0);
		    }
		}
		break;
	    case OPENB2ND:
		{
		    int bno = grid[curgrid].boundaries[0];
		    set_action(OPENB1ST);
		    get_world(x, y, &wx2, &wy2);
		    find_nearest_boundary_point2(bno, wx2, wy2, &ind2);
		    diamond(boundary[bno].boundx[ind2], boundary[bno].boundy[ind2]);
		    define_adcirc_openb(0, ind1, ind2);
		    draw_adcirc_openb(0);
		    flush_pending();
		}
		break;
	    case LANDB1ST:
		{
		    int bno = grid[curgrid].boundaries[0];
		    set_action(LANDB2ND);
		    get_world(x, y, &wx1, &wy1);
		    find_nearest_boundary_point2(bno, wx1, wy1, &ind1);
		    if (ind1 >= 0) {
			diamond(boundary[bno].boundx[ind1], boundary[bno].boundy[ind1]);
			flush_pending();
		    } else {
			errwin("No boundary point found");
			set_action(0);
		    }
		}
		break;
	    case LANDB2ND:
		{
		    int bno = grid[curgrid].boundaries[0];
		    set_action(LANDB1ST);
		    get_world(x, y, &wx2, &wy2);
		    find_nearest_boundary_point2(bno, wx2, wy2, &ind2);
		    diamond(boundary[bno].boundx[ind2], boundary[bno].boundy[ind2]);
		    define_adcirc_landb(0, ind1, ind2);
		    draw_adcirc_landb(0);
		    flush_pending();
		}
		break;
	    case PICK_BATH:
		{
		    rubber_flag = 1;
		    sx = x;
		    sy = y;

		    get_world(x, y, &wx, &wy);
		    slicex[slice_npts] = wx;
		    slicey[slice_npts] = wy;
		    slice_npts++;
		    set_action(PICK_BATH);
		}
		break;
	    case PLACE_MAPSCALE:
		set_action(0);
		device2world(x, y, &wx, &wy);
		set_mapscale_loc(wx, wy);
		break;
	    case PLACE_ISOLINES_LEGEND:
		device2world(x, y, &wx, &wy);
		set_isolleg_loc(wx, wy);
		set_action(0);
		break;
	    case PLACE_VSCALE:
		{
		    extern double grad_legx, grad_legy;
		    get_world(x, y, &grad_legx, &grad_legy);
		    set_action(0);
		}
		break;
	    }
	    break;
	case Button2:
	    switch (action_flag) {
	    case EXTRACT_GRID:
		get_world(x, y, &polyx[npoly], &polyy[npoly]);
		npoly++;
		if (npoly > 2) {
		    extract_grid(polyx, polyy, npoly);
		}
		set_action(0);
		do_drawgrid();
		break;
	    case DEF_EXT_BOUND:
		{
		    int ib, i;

		    get_world(x, y, &polyx[npoly], &polyy[npoly]);
		    npoly++;
		    if (npoly > 2) {
			Free_boundary(grid[curgrid].boundaries[0]);
			if (grid[curgrid].nbounds > 0) {
			    grid[curgrid].nbounds--;
			}
			ib = nextboundary();
			grid[curgrid].boundaries[0] = ib;
			Allocate_boundary(ib, npoly, 0);
			grid[curgrid].nbounds++;
			for (i = 0; i < npoly; i++) {
			    boundary[ib].boundx[i] = polyx[i];
			    boundary[ib].boundy[i] = polyy[i];
			}
			ebound_defined = 1;
			update_fuzz_items();
		    }
		}
		set_action(0);
		do_drawgrid();
		break;
	    case DEF_INT_BOUND:
		{
		    int ib, i;

		    get_world(x, y, &polyx[npoly], &polyy[npoly]);
		    npoly++;
		    if (npoly > 2) {
			ib = nextboundary();
			Allocate_boundary(ib, npoly, 1);
			for (i = 0; i < npoly; i++) {
			    boundary[ib].boundx[i] = polyx[i];
			    boundary[ib].boundy[i] = polyy[i];
			}
			associate_grid(curgrid, ib);
			ibound_defined = 1;
			update_fuzz_items();
		    }
		}
		set_action(0);
		do_drawgrid();
		break;
	    case DEFINE_REGION:
		if (nregion >= 3) {
		    region_flag = 1;
		    draw_region(regionx, regiony, nregion);
		} else {
		    nregion = 0;
		    region_flag = 0;
		}
		break;
	    case PICK_BATH:
		{
		    extern int slice_back;
		    write_slice(curgrid, 1);
		    if (slice_back) {
			write_slice(backgrid, 0);
		    }
		}
		break;
	    case DELETE_BUILD:
		delete_build_points(curbuild, ibuf, nibuf);
		do_drawgrid();
		break;
	    case DELETE_ELEMENT:
	    case DELETE_ELEMENTS:
		delete_elements(curgrid, ibuf, nibuf);
		break;
	    case SPLIT_ELEMENTS3:
		break;
	    case SPLIT_ELEMENTS4:
		break;
	    case SELECT_SNAP:
		{
		    int i;
		    double x, y, r, at, xinter;

		    for (i = 0; i < nibuf; i++) {
			x = grid[curgrid].xord[ibuf[i]];
			y = grid[curgrid].yord[ibuf[i]];
			if (sm == 1e37) {
			    grid[curgrid].xord[ibuf[i]] = sb;
			    grid[curgrid].yord[ibuf[i]] = y;
			} else if (sm == 0.0) {
			    grid[curgrid].xord[ibuf[i]] = x;
			    grid[curgrid].yord[ibuf[i]] = sb;
			} else {
			    r = y + 1.0 / sm * x;
			    xinter = (r - sb) / (sm + 1.0 / sm);
			    grid[curgrid].xord[ibuf[i]] = xinter;
			    grid[curgrid].yord[ibuf[i]] = sm * xinter + sb;
			}
		    }
		    do_drawgrid();
		    set_action(0);
		}
		break;
	    default:
		break;
	    }
	    npoly = 0;
	    nipoly = 0;
	    nibuf = 0;
	    slice_npts = 0;
	    set_action(0);
	    break;
	}
	break;
    case MotionNotify:
	x = event->xmotion.x;
	y = event->xmotion.y;
	if (event->xmotion.state & Button1MotionMask) {
	    switch (action_flag) {
/*
 * Delete build points
 */
	    case DELETE_BUILD:
		get_world(x, y, &wx, &wy);
		accumulate_nearest_buildpts(curbuild, wx, wy);
		set_action(DELETE_BUILD);
		break;
/*
 * Delete the element with center of mass nearest mouse pointer
 */
	    case DELETE_ELEMENT:
		get_world(x, y, &wx, &wy);
		accumulate_nearest_elements(curgrid, wx, wy);
		set_action(DELETE_ELEMENT);
		break;
	    }
	}
	getpoints(x, y);
	break;
    default:
	break;
    }
/*
 * some mouse tracking stuff
 */
    switch (action_flag) {
    case MOVE_WORLD2ND:
	select_region(sx, sy, xs, ys);
	sx = sx - (ox - x);
	sy = sy - (oy - y);
	xs = xs - (ox - x);
	ys = ys - (oy - y);
	ox = x;
	oy = y;
	select_region(sx, sy, xs, ys);
	break;
    case MOVE_OBJECT_2ND:
	dx = sx - x;
	dy = sy - y;

	switch (ty) {
	case BOX:
	    select_region(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_region(sx, sy, xs, ys);
	    break;
	case LINE:
	    select_line(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_line(sx, sy, xs, ys);
	    break;
	case STRING:
	    select_region(sx, sy, xs, ys);
	    sx = x;
	    sy = y;
	    xs = xs - dx;
	    ys = ys - dy;
	    select_region(sx, sy, xs, ys);
	    break;
	}
	break;
    case STR_LOC:
	break;
    }

    if (rectflag) {
	select_region(sx, sy, old_x, old_y);
	select_region(sx, sy, x, y);
    }
    if (rotrectflag) {
	select_rotated_rect(sx1, sy1, sx2, sy2, sx3, sy3);
	select_rotated_rect(sx1, sy1, sx2, sy2, sx3, sy3);
    }
    if (rubber_flag) {
	select_line(sx, sy, old_x, old_y);
	select_line(sx, sy, x, y);
    }
    old_x = x;
    old_y = y;
}
