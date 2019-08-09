/* $Id: objutils.c,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 *
 * operations on objects (strings, lines, and boxes)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "globals.h"
#include "noxprotos.h"

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
	    tmp = my_hypot((x - xtmp1), (y - ytmp1));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = my_hypot((x - xtmp1), (y - ytmp2));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = my_hypot((x - xtmp2), (y - ytmp1));
	    if (m > tmp) {
		*type = BOX;
		*numb = i;
		m = tmp;
	    }
	    tmp = my_hypot((x - xtmp2), (y - ytmp2));
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
	    tmp = my_hypot((x - xtmp1), (y - ytmp1));
	    if (m > tmp) {
		*type = LINE;
		*numb = i;
		m = tmp;
	    }
	    tmp = my_hypot((x - xtmp2), (y - ytmp2));
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
	    tmp = my_hypot((x - xtmp1), (y - ytmp1));
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
    char *tmpbuf;
    switch (type) {
	case BOX:
	boxes[to] = boxes[from];
	break;
    case LINE:
	lines[to] = lines[from];
	break;
    case STRING:
	kill_string(to);
	free(pstr[to].s);
	tmpbuf = (char *) malloc((strlen(pstr[from].s) + 1) * sizeof(char));
	pstr[to] = pstr[from];
	pstr[to].s = tmpbuf;
	strcpy(pstr[to].s, pstr[from].s);
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

int define_string(char *s, double wx, double wy)
{
    int i;

    i = next_string();
    if (i >= 0) {
	if (s != NULL) {
	    free(pstr[i].s);
	}
	if (s != NULL) {
	    pstr[i].s = (char *) malloc(sizeof(char) * (strlen(s) + 1));
	    strcpy(pstr[i].s, s);
	} else {
	    pstr[i].s = (char *) malloc(sizeof(char));
	    pstr[i].s[0] = 0;
	}
	pstr[i].font = string_font;
	pstr[i].color = string_color;
	pstr[i].linew = string_linew;
	pstr[i].rot = string_rot;
	pstr[i].charsize = string_size;
	pstr[i].loctype = string_loctype;
	pstr[i].just = string_just;
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

void do_clear_lines(void)
{
    int i;

    for (i = 0; i < MAXLINES; i++) {
	kill_line(i);
    }
    if (inwin) {
	drawgraph();
    }
}

void do_clear_boxes(void)
{
    int i;

    for (i = 0; i < MAXBOXES; i++) {
	kill_box(i);
    }
    if (inwin) {
	drawgraph();
    }
}

void do_clear_text(void)
{
    int i;

    for (i = 0; i < MAXSTR; i++) {
	kill_string(i);
    }
    if (inwin) {
	drawgraph();
    }
}

void realloc_lines(int n)
{
    int i;
    if (n > maxlines) {
	lines = (linetype *) realloc(lines, n * sizeof(linetype));
	for (i = maxlines; i < n; i++) {
	    set_default_line(&lines[i]);
	}
	maxlines = n;
    }
}

void realloc_boxes(int n)
{
    int i;
    if (n > maxboxes) {
	boxes = (boxtype *) realloc(boxes, n * sizeof(boxtype));
	for (i = maxboxes; i < n; i++) {
	    set_default_box(&boxes[i]);
	}
	maxboxes = n;
    }
}

void realloc_strings(int n)
{
    int i;
    if (n > maxstr) {
	pstr = (plotstr *) realloc(pstr, n * sizeof(plotstr));
	for (i = maxstr; i < n; i++) {
	    set_default_string(&pstr[i]);
	}
	maxstr = n;
    }
}
