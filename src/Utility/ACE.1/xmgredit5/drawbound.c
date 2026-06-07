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
 * draw the clipping boundary
 *
 */

#ifndef lint
static char RCSid[] = "$Id: drawbound.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>

#include "defines.h"
#include "globals.h"

int dispbound_flag;
extern int bound_color;

/*
 * draw the boundaries associated with a grid
 */
void draw_boundary(int gridno)
{
    int i, j, bno;
    char s[128];
    for (i = 0; i < grid[gridno].nbounds; i++) {
	bno = grid[gridno].boundaries[i];
	if (bno >= 0 && boundary[bno].bactive) {
	    switch (dispbound_flag) {
	    case 0:
		if (boundary[bno].nbpts == 1) {
		    my_circlefilled(boundary[bno].boundx[0], boundary[bno].boundy[0]);
		} else {
		    my_move2(boundary[bno].boundx[0], boundary[bno].boundy[0]);
		    for (j = 1; j < boundary[bno].nbpts; j++) {
			my_draw2(boundary[bno].boundx[j], boundary[bno].boundy[j]);
		    }
		    if (boundary[bno].boundtype != 3) {
			my_draw2(boundary[bno].boundx[0], boundary[bno].boundy[0]);
		    } else {
			if (boundary[bno].nbpts == 1) {
			    my_circle(boundary[bno].boundx[0],
				      boundary[bno].boundy[0]);
			}
		    }
		}
		break;
	    case 1:
		for (j = 0; j < boundary[bno].nbpts; j++) {
		    my_circle(boundary[bno].boundx[j], boundary[bno].boundy[j]);
		}
		break;
	    case 2:
		setcolor(1);
		for (j = 0; j < boundary[bno].nbpts; j++) {
		    if (symok(boundary[bno].boundx[j], boundary[bno].boundy[j])) {
			sprintf(s, "%d:%d, %d", i + 1, j + 1, boundary[bno].nodes[j] + 1);
			writestr(boundary[bno].boundx[j],
				 boundary[bno].boundy[j],
				 0, 0, s);
		    }
		}
		break;
	    }
	}
    }
    setcolor(1);
}
