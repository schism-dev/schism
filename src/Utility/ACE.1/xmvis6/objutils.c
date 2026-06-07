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
 *
 * operations on objects (strings, lines, and boxes)
 *
 */

#ifndef lint
static char RCSid[] = "$Id: objutils.c,v 1.2 2003/07/24 15:44:06 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

double xconv(), yconv();

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
	return (pstr[strno].s[0]);
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

void copy_object(int type, int from, int to)
{
    switch (type) {
    case BOX:
	boxes[to] = boxes[from];
	break;
    case LINE:
	lines[to] = lines[from];
	break;
    case STRING:
	pstr[to] = pstr[from];
	break;
    }

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

void do_copy_object_proc(void)
{
    set_action(0);
    set_action(COPY_OBJECT1ST);
}

void do_cut_object_proc(void)
{
    set_action(0);
    set_action(CUT_OBJECT);
}

void edit_objects_proc(void)
{
    set_action(0);
    set_action(EDIT_OBJECT);
}

/*
 * print the data part of a plot string
 * for debugging.
 */
void print_plotstr(plotstr str)
{
    printf("%s\n", str.s);
    printf("%s\n", str.active == ON ? "ON" : "OFF");
    printf("font = %d\n", str.font);
    printf("color = %d\n", str.color);
    printf("linew = %d\n", str.linew);
    printf("rot = %d\n", str.rot);
    printf("just = %d\n", str.just);
    printf("gno = %d\n", str.gno);
    printf("loctype = %s\n", str.loctype == VIEW ? "VIEW" : "WORLD");
    printf("loc = (%lf, %lf)\n", str.x, str.y);
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
	pstr[i].font = sysstr.font;
	pstr[i].color = sysstr.color;
	pstr[i].linew = sysstr.linew;
	pstr[i].rot = sysstr.rot;
	pstr[i].charsize = sysstr.charsize;
	pstr[i].loctype = sysstr.loctype;
	pstr[i].just = sysstr.just;
	pstr[i].active = ON;
	if (sysstr.loctype == VIEW) {
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

void strings_loc_proc(void)
{
    set_action(0);
    set_action(STR_LOC);
/*
    xv_set(xv_get(canvas, CANVAS_NTH_PAINT_WINDOW, 0), WIN_CURSOR, cursor_strloc, 0);
*/
}

void strings_ang_proc(void)
{
    set_action(0);
    set_action(STR_LOC1ST);
    /*
     * xv_set(xv_get(canvas, CANVAS_NTH_PAINT_WINDOW, 0),WIN_CURSOR,
     * cursor_strloc, 0);
     */
}

void strings_edit_proc(void)
{
    set_action(0);
    set_action(STR_EDIT);
/*
    xv_set(xv_get(canvas, CANVAS_NTH_PAINT_WINDOW, 0), WIN_CURSOR, cursor_strloc, 0);
*/
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

void set_plotstr_string(plotstr * pstr, char *buf)
{
    strcpy(pstr->s, buf);
}

/* TODO dynmically allocate strings.
void kill_string(int stringno)
{
    if (pstr[stringno].s != NULL) {
        free(pstr[stringno].s);
    }
    pstr[stringno].s = (char *) malloc(sizeof(char));
    pstr[stringno].s[0] = 0;
    pstr[stringno].active = OFF;
}

void set_plotstr_string(plotstr * pstr, char *buf)
{
    if (pstr->s != NULL) {
        free(pstr->s);
    }
    pstr->s = NULL;
    if (buf != NULL) {
        pstr->s = (char *) malloc(sizeof(char) * (strlen(buf) + 1));
        strcpy(pstr->s, buf);
    } else {
        pstr->s = (char *) malloc(sizeof(char));
        pstr->s[0] = 0;
    }
}

plotstr copy_plotstr(plotstr p)
{
    static plotstr pto;
    pto = p;
    if (p.s != NULL) {
        pto.s = (char *) malloc((strlen(p.s) + 1) * sizeof(char));
        if (pto.s != NULL) {
            strcpy(pto.s, p.s);
        }
    } else {
        pto.s = NULL;
    }
    return pto;
}
*/
