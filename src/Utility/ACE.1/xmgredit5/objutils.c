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
 * operations on objects (strings, lines, and boxes)
 *
 */

#ifndef lint
static char RCSid[] = "$Id: objutils.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "globals.h"
#include "graphics.h"

double xconv(), yconv();

static int cg = 0;

void get_graph_box(int i, boxtype * b)
{
    memcpy(b, &boxes[i], sizeof(boxtype));
}

void get_graph_line(int i, linetype * l)
{
    memcpy(l, &lines[i], sizeof(linetype));
}

void get_graph_string(int i, plotstr * s)
{
    memcpy(s, &pstr[i], sizeof(plotstr));
}

/*
 * find the nearest object to (x,y)
 */
void find_item(int gno, double x, double y, int *type, int *numb)
{
    int i;
    double tmp, xtmp1, ytmp1, xtmp2, ytmp2, m = 1e307;
    double dx, dy;
    boxtype box;
    linetype line;
    plotstr str;

    x = xconv(x);
    y = yconv(y);
    *type = -1;
    for (i = 0; i < MAXBOXES; i++) {
	get_graph_box(i, &box);
	if (isactive_box(i)) {
	    if (box.loctype == VIEW) {
		xtmp1 = box.x1;
		ytmp1 = box.y1;
		xtmp2 = box.x2;
		ytmp2 = box.y2;
	    } else {
		if (gno == box.gno) {
		    xtmp1 = xconv(box.x1);
		    ytmp1 = yconv(box.y1);
		    xtmp2 = xconv(box.x2);
		    ytmp2 = yconv(box.y2);
		} else {
		    continue;
		}
	    }
	    tmp = hypot((x - xtmp1), (y - ytmp1));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = hypot((x - xtmp1), (y - ytmp2));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = hypot((x - xtmp2), (y - ytmp1));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = hypot((x - xtmp2), (y - ytmp2));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	}
    }
    for (i = 0; i < MAXLINES; i++) {
	get_graph_line(i, &line);
	if (isactive_line(i)) {
	    if (line.loctype == VIEW) {
		xtmp1 = line.x1;
		ytmp1 = line.y1;
		xtmp2 = line.x2;
		ytmp2 = line.y2;
	    } else {
		if (gno == line.gno) {
		    xtmp1 = xconv(line.x1);
		    ytmp1 = yconv(line.y1);
		    xtmp2 = xconv(line.x2);
		    ytmp2 = yconv(line.y2);
		} else {
		    continue;
		}
	    }
	    tmp = hypot((x - xtmp1), (y - ytmp1));
	    if (m > tmp) {
		*type = LINE;
		*numb = i;
		m = tmp;
	    }
	    tmp = hypot((x - xtmp2), (y - ytmp2));
	    if (m > tmp) {
		*type = LINE;
		*numb = i;
		m = tmp;
	    }
	}
    }
    for (i = 0; i < MAXSTR; i++) {
	get_graph_string(i, &str);
	if (isactive_string(i)) {
	    if (str.loctype == VIEW) {
		xtmp1 = str.x;
		ytmp1 = str.y;
	    } else {
		if (gno == str.gno) {
		    xtmp1 = xconv(str.x);
		    ytmp1 = yconv(str.y);
		} else {
		    continue;
		}
	    }
	    tmp = hypot((x - xtmp1), (y - ytmp1));
	    if (m > tmp) {
		*type = STRING;
		*numb = i;
		m = tmp;
	    }
	}
    }
}

int isactive_line(int lineno)
{
    if (0 <= lineno && lineno < MAXLINES)
	return (lines[lineno].active == ON);
    return (0);
}

int isactive_box(int boxno)
{
    if (0 <= boxno && boxno < MAXBOXES)
	return (boxes[boxno].active == ON);
    return (0);
}

int isactive_string(int strno)
{
    if (0 <= strno && strno < MAXSTR)
	return ((pstr[strno].active == ON && pstr[strno].s[0])
		|| (pstr[strno].active == ON && pstr[strno].type == 1));
    return (0);
}

int next_line(void)
{
    int i;

    for (i = 0; i < MAXLINES; i++) {
	if (!isactive_line(i)) {
	    lines[i].active = ON;
	    return (i);
	}
    }
    errwin("Error - no lines available");
    return (-1);
}

int next_box(void)
{
    int i;

    for (i = 0; i < MAXBOXES; i++) {
	if (!isactive_box(i)) {
	    boxes[i].active = ON;
	    return (i);
	}
    }
    errwin("Error - no boxes available");
    return (-1);
}

int next_string(void)
{
    int i;

    for (i = 0; i < MAXSTR; i++) {
	if (!isactive_string(i)) {
	    return (i);
	}
    }
    errwin("Error - no strings available");
    return (-1);
}

void kill_box(int boxno)
{
    boxes[boxno].active = OFF;
}

void kill_line(int lineno)
{
    lines[lineno].active = OFF;
}

void kill_string(int stringno)
{
    pstr[stringno].active = OFF;
    pstr[stringno].type = 0;
    pstr[stringno].s[0] = 0;
}

void do_boxes_proc(void)
{
    set_action(0);
    set_action(MAKE_BOX_1ST);
}

void do_lines_proc(void)
{
    set_action(0);
    set_action(MAKE_LINE_1ST);
}

void do_move_proc(void)
{
    set_action(0);
    set_action(MOVE_OBJECT_1ST);
}

void do_delete_object_proc(void)
{
    set_action(0);
    set_action(DEL_OBJECT);
}

void edit_objects_proc(void)
{
    set_action(0);
    set_action(EDIT_OBJECT);
}

int define_string(char *s, double wx, double wy)
{
    int i;

    i = next_string();
    if (i >= 0) {
	if (s != NULL) {
	    strcpy(pstr[i].s, s);
	} else {
	    pstr[i].s[0] = 0;
	    pstr[i].active = OFF;
	}
	pstr[i].font = string_font;
	pstr[i].color = string_color;
	pstr[i].linew = string_linew;
	pstr[i].rot = string_rot;
	pstr[i].charsize = string_size;
	pstr[i].loctype = string_loctype;
	pstr[i].just = string_just;
	pstr[i].type = 0;
	pstr[i].sym = string_sym;
	pstr[i].active = ON;
	if (string_loctype == VIEW) {
	    pstr[i].x = xconv(wx);
	    pstr[i].y = yconv(wy);
	    pstr[i].gno = -1;
	} else {
	    pstr[i].x = wx;
	    pstr[i].y = wy;
	    pstr[i].gno = cg;
	}
	return i;
    }
    return -1;
}

int define_depth_label(int gridno, int el, double wx, double wy)
{
    int i;

    i = next_string();
    if (i >= 0) {
	pstr[i].font = string_font;
	pstr[i].color = string_color;
	pstr[i].linew = string_linew;
	pstr[i].rot = string_rot;
	pstr[i].charsize = string_size;
	pstr[i].loctype = string_loctype;
	pstr[i].just = string_just;
	pstr[i].active = ON;
	pstr[i].type = 1;
	pstr[i].sym = string_sym;
	pstr[i].symsize = string_symsize;
	pstr[i].x = wx;
	pstr[i].y = wy;
	pstr[i].gno = gridno;
	pstr[i].el = el;
	return i;
    }
    return -1;
}

void strings_loc_proc(void)
{
    set_action(0);
    set_action(STR_LOC);
}

void strings_depth_proc(void)
{
    set_action(0);
    set_action(PLACE_ISOLINE_LABEL);
}

void strings_ang_proc(void)
{
    set_action(0);
    set_action(STR_LOC1ST);
}

void strings_edit_proc(void)
{
    set_action(0);
    set_action(STR_EDIT);
}

void do_clear_lines(void)
{
    int i;

    for (i = 0; i < MAXLINES; i++) {
	kill_line(i);
    }
}

void do_clear_boxes(void)
{
    int i;

    for (i = 0; i < MAXBOXES; i++) {
	kill_box(i);
    }
}

void do_clear_text(void)
{
    int i;

    for (i = 0; i < MAXSTR; i++) {
	kill_string(i);
    }
}

/*
 * draw annotative text
 */
void draw_string(int gno, int i)
{
    double xtmp1, ytmp1;
    int f, c, w, tmpind;
    double s, d;
    plotstr pstr;
    double get_depth_element();

    get_graph_string(i, &pstr);

    if (pstr.type == 1) {
	if (gno == -1)
	    return;
	find_element(pstr.gno, pstr.x, pstr.y, &tmpind);
	if (tmpind >= 0) {
	    d = get_depth_element(pstr.gno, tmpind, pstr.x, pstr.y);
	} else {
	    return;
	}
	if (pstr.sym) {
	    sprintf(pstr.s, " %.*lf", string_prec, d);
	} else {
	    sprintf(pstr.s, "%.*lf", string_prec, d);
	}
    } else {
	if (gno != -2) {
	    if (pstr.loctype == WORLD && pstr.gno != gno) {
		return;
	    }
	    if (pstr.loctype == VIEW && gno != -1) {
		return;
	    }
	}
    }
    if (strlen(pstr.s) && (pstr.charsize > 0.0) && (pstr.active == ON)) {
	c = setcolor(pstr.color);
	w = setlinewidth(pstr.linew);
	(void) setlinestyle(1);
	s = setcharsize(pstr.charsize);
	f = setfont(pstr.font);
	if (pstr.loctype == WORLD) {
	    writestr(pstr.x, pstr.y, pstr.rot, pstr.just, pstr.s);
	} else {
	    view2world(pstr.x, pstr.y, &xtmp1, &ytmp1);
	    writestr(xtmp1, ytmp1, pstr.rot, pstr.just, pstr.s);
	}
	if (pstr.type == 1 && pstr.sym > 0) {
	    my_circlefilled(pstr.x, pstr.y);
	}
	(void) setcolor(c);
	(void) setlinewidth(w);
	(void) setlinestyle(1);
	(void) setcharsize(s);
	(void) setfont(f);
    }
}

/*
 * draw annotative boxes
 */
void draw_box(int gno, int i)
{
    double xtmp1, ytmp1;
    double xtmp2, ytmp2;
    int c, l, w;
    boxtype b;

    get_graph_box(i, &b);
    if (gno != -2) {
	if (b.loctype == WORLD && b.gno != gno) {
	    return;
	}
	if (b.loctype == VIEW && gno != -1) {
	    return;
	}
    }
    if (b.active == ON) {
	setclipping(0);

	if (b.fill == COLOR) {
	    c = setcolor(b.fillcolor);
	    if (b.loctype == WORLD) {
		fillrectcolor(b.x1, b.y1, b.x2, b.y2);
	    } else {
		view2world(b.x1, b.y1, &xtmp1, &ytmp1);
		view2world(b.x2, b.y2, &xtmp2, &ytmp2);
		fillrectcolor(xtmp1, ytmp1, xtmp2, ytmp2);
	    }
	    setcolor(c);
	} else if (b.fill == PATTERN) {
	    c = setpattern(b.fillpattern);
	    if (b.loctype == WORLD) {
		fillrectpat(b.x1, b.y1, b.x2, b.y2);
	    } else {
		view2world(b.x1, b.y1, &xtmp1, &ytmp1);
		view2world(b.x2, b.y2, &xtmp2, &ytmp2);
		fillrectpat(xtmp1, ytmp1, xtmp2, ytmp2);
	    }
	}
	c = setcolor(b.color);
	l = setlinestyle(b.lines);
	w = setlinewidth(b.linew);
	if (b.loctype == WORLD) {
	    rect(b.x1, b.y1, b.x2, b.y2);
	} else {
	    view2world(b.x1, b.y1, &xtmp1, &ytmp1);
	    view2world(b.x2, b.y2, &xtmp2, &ytmp2);
	    rect(xtmp1, ytmp1, xtmp2, ytmp2);
	}
	setclipping(1);
	setcolor(c);
	setlinewidth(w);
	setlinestyle(l);
    }
}

/*
 * draw annotative lines
 */
void draw_line(int gno, int i)
{
    double xtmp1, ytmp1;
    double xtmp2, ytmp2;
    int c, ll, w;
    linetype l;

    get_graph_line(i, &l);
    if (gno != -2) {
	if (l.loctype == WORLD && l.gno != gno) {
	    return;
	}
	if (l.loctype == VIEW && gno != -1) {
	    return;
	}
    }
    if (l.active == ON) {
	setclipping(0);
	c = setcolor(l.color);
	ll = setlinestyle(l.lines);
	w = setlinewidth(l.linew);
	if (l.loctype == WORLD) {
	    draw_arrow(l.x1, l.y1, l.x2, l.y2, l.arrow, l.asize, l.atype);
	} else {
	    view2world(l.x1, l.y1, &xtmp1, &ytmp1);
	    view2world(l.x2, l.y2, &xtmp2, &ytmp2);
	    draw_arrow(xtmp1, ytmp1, xtmp2, ytmp2, l.arrow, l.asize, l.atype);
	}
	setclipping(1);
	setcolor(c);
	setlinewidth(w);
	setlinestyle(ll);
    }
}

void draw_annotation(int gno)
{
    int i;

    setclipping(0);		/* shut down clipping for strings, boxes,
				 * lines, and legends */
    for (i = 0; i < MAXBOXES; i++) {
	if (isactive_box(i)) {
	    draw_box(gno, i);
	}
    }
    for (i = 0; i < MAXLINES; i++) {
	if (isactive_line(i)) {
	    draw_line(gno, i);
	}
    }
    for (i = 0; i < MAXSTR; i++) {
	if (isactive_string(i)) {
	    if (debuglevel == 5) {
		printf("String %d\n", i);
	    }
	    draw_string(gno, i);
	}
    }
    setclipping(1);
}
