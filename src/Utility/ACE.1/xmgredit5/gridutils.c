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
 * Utilities for manipulating grids, constructing edge lists, splitting
 * elements, ...
 *
 */

#ifndef lint
static char RCSid[] = "$Id: gridutils.c,v 1.11 2008/01/17 16:27:48 pturner Exp $";

#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "globals.h"

static double sm, sb;

char item_buf[256];

double get_depth_element(int gno, int ind, double x, double y);
void split_elem4(int gridno, int i);
void split_elem3(int gridno, int i);
void get_center(int gridno, int i, double *cx, double *cy);
void interpdepths(int g1, int g2);
void end_find(void);
void prep_find(int gno);
void getadjnodes(int gridno);

int IsCCW(Grid * g, int elno);
double ElementArea(Grid * g, int elno);
double ElementVolume(Grid * g, int elno);
void FixAreas(Grid * g);
void FixArea(Grid * g, int elno);

void compute_circumcenter(double xk, double yk, double xl, double yl, double xm, double ym, double *xc, double *yc);

static int findfxmin(double s);

/*
 * Return the number of nodes for each element type
 */
int numnodes(int gridno, int el)
{
    switch (grid[gridno].icon[el].type) {
    case 0:
	return 3;
    case 1:
	return 3;
    case 2:
	return 2;
    case 3:
	return 3;
    case 4:
	return 4;
    case 5:
	return 5;
    case 6:
	return 6;
    case 7:
	return -1;
    case 8:
	return 8;
    }
    return -1;
}

/*
 * Return the number of nodes in el and a list of
 * nodes in nl
 */
int getallnodes(int gridno, int el, int *nl)
{
    int i;
    for (i = 0; i < 8; i++) {
	nl[i] = grid[gridno].icon[el].nl[i];
    }
    return numnodes(gridno, el);
}

/*
 * Return the number of nodes and the list of geometry nodes
 * on nl
 */
int getallgeomnodes(int gridno, int el, int *nl)
{
    int i;
    for (i = 0; i < 8; i++) {
	nl[i] = grid[gridno].icon[el].nl[i];
    }
    return numnodes(gridno, el);
}

/*
 * Return the type of element el
 */
int getelementtype(int gridno, int el)
{
    return grid[gridno].icon[el].type;
}

/*
 * Return the type of node n
 */
int getnodetype(int gridno, int n)
{
    return grid[gridno].ntype[n];
}

/*
 * test the grid for composition of purely triangles
 */
int grid_is_triangles(int gridno)
{
    int i, j;
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type != 3) {
	    return 0;
	}
    }
    return 1;
}

/* 
 * test the node for colocation
 */
int node_is_colocated(int gridno, int n)
{
    int i, j;
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    if (grid[gridno].icon[i].nl[j] == n) {
		if (grid[gridno].icon[i].type == 0) {
		    return i;
		}
	    }
	}
    }
    return -1;
}

/*

 */
int next_free_colocated_node(int gridno, int el)
{
    int i, j, k, found, n;
    for (k = 0; k < grid[gridno].icon[el].nn; k++) {
	found = 0;
	for (i = 0; i < grid[gridno].nmel; i++) {
	    if (i != el) {
		for (j = 0; j < grid[gridno].icon[i].nn; j++) {
		    if (grid[gridno].icon[i].nl[j] == n) {
			if (grid[gridno].icon[i].type == 0) {
			    found = 1;
			}
		    }
		}
	    } else {
		found = 1;	/* element i == element el, so skip */
	    }
	}
	if (!found) {
	    return k;
	}
    }
    /* shouldn't happen */
    return -1;
}

/*
 * test a grid for existence of 1-d elements
 */
int grid_has_1delements(int gridno)
{
    int i;
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type == 2 || grid[gridno].icon[i].type == 1 || grid[gridno].icon[i].type == 5 || grid[gridno].icon[i].type == 0) {
	    return 1;
	}
    }
    return 0;
}

/*
 * test an element for type
 */
int is_1delement(int gridno, int i)
{
    if (grid[gridno].icon[i].type == 2 || grid[gridno].icon[i].type == 1 || grid[gridno].icon[i].type == 5 || grid[gridno].icon[i].type == 0) {
	return 1;
    }
    return 0;
}

/*
 * Test a grid for being composed of quadrangles
 */
int grid_is_quadrangles(int gridno)
{
    int i, j;
    int etypes[8];
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type != 4) {
	    return 0;
	}
    }
    return 1;
}

/*
 * delete the elements given in the array elar
 */
void delete_elements(int gridno, int *elar, int npts)
{
    int i, j;

    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
    }
    for (i = 0; i < npts; i++) {
	grid[gridno].ellist[elar[i]] = 1;
    }
    compact_grid(gridno);
}

/*
 * Nodes have been renumbered, load them back to the grid
 */
void load_renumbered_nodes(int gridno)
{
    int i, j, n = 0;
    double *xtmp, *ytmp, *dtmp;
    int *ntmp;

    xtmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));
    ytmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));
    dtmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));
    ntmp = (int *) calloc(grid[gridno].nmnp, sizeof(int));

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].nlist[i] == -1) {
	    grid[gridno].nlist[i] = n;
	    n++;
	}
    }
    for (i = 0; i < grid[gridno].nmnp; i++) {
	xtmp[grid[gridno].nlist[i]] = grid[gridno].xord[i];
	ytmp[grid[gridno].nlist[i]] = grid[gridno].yord[i];
	dtmp[grid[gridno].nlist[i]] = grid[gridno].depth[i];
	ntmp[grid[gridno].nlist[i]] = grid[gridno].ntype[i];
    }
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].xord[i] = xtmp[i];
	grid[gridno].yord[i] = ytmp[i];
	grid[gridno].depth[i] = dtmp[i];
	grid[gridno].ntype[i] = ntmp[i];
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    n = grid[gridno].icon[i].nl[j];
	    grid[gridno].icon[i].nl[j] = grid[gridno].nlist[n];
	}
    }
    free(xtmp);
    free(ytmp);
    free(dtmp);
    free(ntmp);
}

/*
 * Delete nodes given in the array nodes
 */
void delete_nodes(int gridno, int *nodes, int ntodelete)
{
    int i, j, newnmnp = 0, newnmel = 0;
    int saveelem;

    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    for (i = 0; i < ntodelete; i++) {
	grid[gridno].nlist[nodes[i]] = 1;
    }
    compact_grid(gridno);
}

/*
 * cut a grid along a given line
 */
void cut_grid(int gridno, double x1, double y1, double x2, double y2, int cutgrid_type)
{
    int i, j;
    double m, b, *xord = grid[gridno].xord, *yord = grid[gridno].yord;
    int ic;

    if (x2 - x1 != 0.0) {
	m = (y2 - y1) / (x2 - x1);
    } else {
	return;
    }
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    b = y1 - m * x1;
    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (cutgrid_type == 0) {
	    ic = (yord[i] > (m * xord[i] + b));
	} else {
	    ic = (yord[i] < (m * xord[i] + b));
	}
	if (ic) {
	    grid[gridno].nlist[i] = 1;
	}
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    if (grid[gridno].nlist[grid[gridno].icon[i].nl[j]]) {
		grid[gridno].ellist[i] = 1;
	    }
	}
    }
    compact_grid(gridno);
}

/*
 * Loop through a grid and find elements that are not
 * oriented in a counter-clockwise fashion - mark them
 * on the drawing area and sawp nodes to fix them.
 */
void fix_areas(int gridno)
{
    int i, ne;
    double xg, yg;

    ne = grid[gridno].nmel;
    for (i = 0; i < ne; i++) {
	if (!comparea(gridno, grid[gridno].icon[i].nl[0], grid[gridno].icon[i].nl[1], grid[gridno].icon[i].nl[2])) {
	    get_center(gridno, i, &xg, &yg);
	    my_move2(xg, yg);
	    writetext("-");
	    iswap(&grid[gridno].icon[i].nl[0], &grid[gridno].icon[i].nl[2]);
	}
    }
}

/*
 * Split every element in a grid into 3 or 4 elements
 * depending on the value of it
 */
void subdividegrid(int gridno, int it)
{
    int i, ne;

    ne = grid[gridno].nmel;
    for (i = 0; i < grid[gridno].nmnp; i++) {
	grid[gridno].nlist[i] = 0;
    }
    for (i = 0; i < ne; i++) {
	if (it == 3) {
	    split_elem3(gridno, i);
	} else if (it == 4) {
	    split_elem4(gridno, i);
	}
    }
}

/*
 * queue up nodes and elements for bulk operations
 */
#define BUFSIZE 1000

/*
 * queue for nodes
 */
static double *qx, *qy, *qd;
static int qn;

/*
 * queue for elements
 */
static int *q1, *q2, *q3;
static int qe;

/*
 * initialize q for nodes
 */
void init_qadd_node(int gridno)
{
    if (qx) {
	free(qx);
	free(qy);
	free(qd);
    }
    qn = 0;
    qx = (double *) calloc(BUFSIZE, sizeof(double));
    qy = (double *) calloc(BUFSIZE, sizeof(double));
    qd = (double *) calloc(BUFSIZE, sizeof(double));
}

/*
 * initialize q for elements
 */
void init_qadd_element(int gridno)
{
    if (q1) {
	free(q1);
	free(q2);
	free(q3);
    }
    qe = 0;
    q1 = (int *) calloc(BUFSIZE, sizeof(int));
    q2 = (int *) calloc(BUFSIZE, sizeof(int));
    q3 = (int *) calloc(BUFSIZE, sizeof(int));
}

/*
 * Add nodes
 */
void flush_qadd_node(int gridno)
{
    add_nodes(gridno, qn, qx, qy, qd);
    if (qx) {
	free(qx);
	free(qy);
	free(qd);
	qx = NULL;
    }
}

/*
 * Add elements
 */
void flush_qadd_element(int gridno)
{
    add_elements(gridno, qe, q1, q2, q3);
    if (q1) {
	free(q1);
	free(q2);
	free(q3);
	q1 = NULL;
    }
}

/*
 * Add a single node to the q
 */
void qadd_node(int gridno, double x, double y, double d)
{
    qx[qn] = x;
    qy[qn] = y;
    qd[qn] = d;
    qn++;
    if (qn % BUFSIZE == 0) {
	qx = (double *) realloc(qx, (qn + BUFSIZE) * sizeof(double));
	qy = (double *) realloc(qy, (qn + BUFSIZE) * sizeof(double));
	qd = (double *) realloc(qd, (qn + BUFSIZE) * sizeof(double));
    }
}

/*
 * Add a single element to the q
 */
void qadd_element(int gridno, int n1, int n2, int n3)
{
    q1[qe] = n1;
    q2[qe] = n2;
    q3[qe] = n3;
    qe++;
    if (qe % BUFSIZE == 0) {
	q1 = (int *) realloc(q1, (qe + BUFSIZE) * sizeof(int));
	q2 = (int *) realloc(q2, (qe + BUFSIZE) * sizeof(int));
	q3 = (int *) realloc(q3, (qe + BUFSIZE) * sizeof(int));
    }
}

/*
 * Split all elements into 4, using an edgelist for performance
 */
void split_allelem4(int gridno)
{
    double x1, y1, x2, y2, x3, y3;
    double d1, d2, d3;
    double x4, y4, x5, y5, x6, y6;
    double d4, d5, d6;
    int nmnpm1, nmnporg, i, j, n1, n2, n3, n4, n5, n6, el[3], nmnp, nmel = grid[gridno].nmel;
    int cnt1, cnt2, addn1n2, addn2n3, addn1n3;
    struct edgelist {
	int *list;
    } *e;

    nmnporg = nmnp = grid[gridno].nmnp;
    init_qadd_node(gridno);
    init_qadd_element(gridno);
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) calloc(2, sizeof(int));
	e[i].list[0] = -1;
	e[i].list[1] = -1;
    }

    nmnpm1 = nmnp - 1;
    for (i = 0; i < nmel; i++) {
	if (grid[gridno].ellist[i]) {
	    n1 = grid[gridno].icon[i].nl[0];
	    x1 = grid[gridno].xord[n1];
	    y1 = grid[gridno].yord[n1];
	    d1 = grid[gridno].depth[n1];
	    n2 = grid[gridno].icon[i].nl[1];
	    x2 = grid[gridno].xord[n2];
	    y2 = grid[gridno].yord[n2];
	    d2 = grid[gridno].depth[n2];
	    n3 = grid[gridno].icon[i].nl[2];
	    x3 = grid[gridno].xord[n3];
	    y3 = grid[gridno].yord[n3];
	    d3 = grid[gridno].depth[n3];
	    x4 = 0.5 * (x1 + x2);
	    y4 = 0.5 * (y1 + y2);
	    d4 = 0.5 * (d1 + d2);
	    x6 = 0.5 * (x1 + x3);
	    y6 = 0.5 * (y1 + y3);
	    d6 = 0.5 * (d1 + d3);
	    x5 = 0.5 * (x2 + x3);
	    y5 = 0.5 * (y2 + y3);
	    d5 = 0.5 * (d2 + d3);

	    cnt1 = 0;
	    cnt2 = 0;
	    while (e[n1].list[2 * cnt1] != -1 && e[n1].list[2 * cnt1] != n2)
		cnt1++;
	    n4 = e[n1].list[2 * cnt1 + 1];
	    if (n4 == -1) {
		while (e[n2].list[2 * cnt2] != -1 && e[n2].list[2 * cnt2] != n1)
		    cnt2++;
		n4 = e[n2].list[2 * cnt2 + 1];
	    }
	    if (n4 == -1) {
		n4 = nmnp++;

		e[n1].list = (int *) realloc(e[n1].list, 2 * (cnt1 + 2) * sizeof(int));
		e[n1].list[2 * cnt1] = n2;
		e[n1].list[2 * cnt1 + 1] = n4;
		cnt1++;
		e[n1].list[2 * cnt1] = -1;
		e[n1].list[2 * cnt1 + 1] = -1;

		e[n2].list = (int *) realloc(e[n2].list, 2 * (cnt2 + 2) * sizeof(int));
		e[n2].list[2 * cnt2] = n1;
		e[n2].list[2 * cnt2 + 1] = n4;
		cnt2++;
		e[n2].list[2 * cnt2] = -1;
		e[n2].list[2 * cnt2 + 1] = -1;

		qadd_node(gridno, x4, y4, d4);
	    }
	    cnt1 = 0;
	    cnt2 = 0;
	    while (e[n2].list[2 * cnt1] != -1 && e[n2].list[2 * cnt1] != n3)
		cnt1++;
	    n5 = e[n2].list[2 * cnt1 + 1];
	    if (n5 == -1) {
		while (e[n3].list[2 * cnt2] != -1 && e[n3].list[2 * cnt2] != n2)
		    cnt2++;
		n5 = e[n3].list[2 * cnt2 + 1];
	    }
	    if (n5 == -1) {
		n5 = nmnp++;

		e[n2].list = (int *) realloc(e[n2].list, 2 * (cnt1 + 2) * sizeof(int));
		e[n2].list[2 * cnt1] = n3;
		e[n2].list[2 * cnt1 + 1] = n5;
		cnt1++;
		e[n2].list[2 * cnt1] = -1;
		e[n2].list[2 * cnt1 + 1] = -1;

		e[n3].list = (int *) realloc(e[n3].list, 2 * (cnt2 + 2) * sizeof(int));
		e[n3].list[2 * cnt2] = n2;
		e[n3].list[2 * cnt2 + 1] = n5;
		cnt2++;
		e[n3].list[2 * cnt2] = -1;
		e[n3].list[2 * cnt2 + 1] = -1;

		qadd_node(gridno, x5, y5, d5);
	    }
	    cnt1 = 0;
	    cnt2 = 0;
	    while (e[n3].list[2 * cnt1] != -1 && e[n3].list[2 * cnt1] != n1)
		cnt1++;
	    n6 = e[n3].list[2 * cnt1 + 1];
	    if (n6 == -1) {
		while (e[n1].list[2 * cnt2] != -1 && e[n1].list[2 * cnt2] != n3)
		    cnt2++;
		n6 = e[n1].list[2 * cnt2 + 1];
	    }
	    if (n6 == -1) {
		n6 = nmnp++;

		e[n3].list = (int *) realloc(e[n3].list, 2 * (cnt1 + 2) * sizeof(int));
		e[n3].list[2 * cnt1] = n1;
		e[n3].list[2 * cnt1 + 1] = n6;
		cnt1++;
		e[n3].list[2 * cnt1] = -1;
		e[n3].list[2 * cnt1 + 1] = -1;

		e[n1].list = (int *) realloc(e[n1].list, 2 * (cnt2 + 2) * sizeof(int));
		e[n1].list[2 * cnt2] = n3;
		e[n1].list[2 * cnt2 + 1] = n6;
		cnt2++;
		e[n1].list[2 * cnt2] = -1;
		e[n1].list[2 * cnt2 + 1] = -1;

		qadd_node(gridno, x6, y6, d6);
	    }
	    grid[gridno].icon[i].type = 3;
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	    grid[gridno].icon[i].nl[0] = n1;
	    grid[gridno].icon[i].nl[1] = n4;
	    grid[gridno].icon[i].nl[2] = n6;
	    qadd_element(gridno, n2, n5, n4);
	    qadd_element(gridno, n3, n6, n5);
	    qadd_element(gridno, n4, n5, n6);
	}
    }
    for (i = 0; i < nmnporg; i++) {
	if (e[i].list != NULL) {
	    free(e[i].list);
	}
    }
    free(e);
    flush_qadd_node(gridno);
    flush_qadd_element(gridno);
}

/*
 * determine if a node given by x, y is already in
 * the table of nodes, possible problems here with
 * the 0.00001 business.
 */
int inlist(int gridno, double x, double y)
{
    int i;
    double *xtmp = grid[gridno].xord, *ytmp = grid[gridno].yord;

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (hypot((xtmp[i] - x), (ytmp[i] - y)) < 0.00001) {
	    return i;
	}
    }
    return -1;
}

/*
 * split a single element into 4, being careful not
 * to include the added nodes more than once to the table of nodes.
 */
void split_elem4(int gridno, int i)
{
    double x1, y1, x2, y2, x3, y3;
    double d1, d2, d3;
    double x4, y4, x5, y5, x6, y6;
    double d4, d5, d6;
    int n1, n2, n3, n4, n5, n6, el[3], nmnp;

    nmnp = grid[gridno].nmnp;
    n1 = grid[gridno].icon[i].nl[0];
    x1 = grid[gridno].xord[n1];
    y1 = grid[gridno].yord[n1];
    d1 = grid[gridno].depth[n1];
    n2 = grid[gridno].icon[i].nl[1];
    x2 = grid[gridno].xord[n2];
    y2 = grid[gridno].yord[n2];
    d2 = grid[gridno].depth[n2];
    n3 = grid[gridno].icon[i].nl[2];
    x3 = grid[gridno].xord[n3];
    y3 = grid[gridno].yord[n3];
    d3 = grid[gridno].depth[n3];
    x4 = 0.5 * (x1 + x2);
    y4 = 0.5 * (y1 + y2);
    d4 = 0.5 * (d1 + d2);
    x6 = 0.5 * (x1 + x3);
    y6 = 0.5 * (y1 + y3);
    d6 = 0.5 * (d1 + d3);
    x5 = 0.5 * (x2 + x3);
    y5 = 0.5 * (y2 + y3);
    d5 = 0.5 * (d2 + d3);
    if ((n4 = inlist(gridno, x4, y4)) < 0) {
	add_node(gridno, x4, y4, d4);
	n4 = nmnp++;
    }
    if ((n5 = inlist(gridno, x5, y5)) < 0) {
	add_node(gridno, x5, y5, d5);
	n5 = nmnp++;
    }
    if ((n6 = inlist(gridno, x6, y6)) < 0) {
	add_node(gridno, x6, y6, d6);
	n6 = nmnp++;
    }
    grid[gridno].icon[i].type = 3;
    grid[gridno].icon[i].nn = 3;
    grid[gridno].icon[i].ngeom = 3;
    grid[gridno].icon[i].nl[0] = n1;
    grid[gridno].icon[i].nl[1] = n4;
    grid[gridno].icon[i].nl[2] = n6;
    add_element(gridno, n2, n5, n4);
    add_element(gridno, n3, n6, n5);
    add_element(gridno, n4, n5, n6);
    my_move2(x4, y4);
    my_draw2(x5, y5);
    my_draw2(x6, y6);
    my_draw2(x4, y4);
    my_move2(x5, y5);
}

/*
 * split an element into 3 elements by adding a middle node
 */
void split_elem3(int gridno, int i)
{
    double x1, y1, x2, y2, x3, y3;
    double d1, d2, d3;
    double x4, y4, d4;
    int n1, n2, n3, n4;

    n4 = grid[gridno].nmnp;
    n1 = grid[gridno].icon[i].nl[0];
    x1 = grid[gridno].xord[n1];
    y1 = grid[gridno].yord[n1];
    d1 = grid[gridno].depth[n1];
    n2 = grid[gridno].icon[i].nl[1];
    x2 = grid[gridno].xord[n2];
    y2 = grid[gridno].yord[n2];
    d2 = grid[gridno].depth[n2];
    n3 = grid[gridno].icon[i].nl[2];
    x3 = grid[gridno].xord[n3];
    y3 = grid[gridno].yord[n3];
    d3 = grid[gridno].depth[n3];
    x4 = 0.33333333333 * (x1 + x2 + x3);
    y4 = 0.33333333333 * (y1 + y2 + y3);
    d4 = 0.33333333333 * (d1 + d2 + d3);
    add_node(gridno, x4, y4, d4);
    grid[gridno].icon[i].type = 3;
    grid[gridno].icon[i].nn = 3;
    grid[gridno].icon[i].ngeom = 3;
    grid[gridno].icon[i].nl[0] = n1;
    grid[gridno].icon[i].nl[1] = n2;
    grid[gridno].icon[i].nl[2] = n4;
    if (!comparea(gridno, n1, n2, n4)) {
	fprintf(stderr, "Negative area in split3 n1, n2, n4\n");
    }
    add_element(gridno, n2, n3, n4);
    if (!comparea(gridno, n2, n3, n4)) {
	fprintf(stderr, "Negative area in split3 n2, n3, n4\n");
    }
    add_element(gridno, n1, n4, n3);
    if (!comparea(gridno, n1, n4, n3)) {
	fprintf(stderr, "Negative area in split3 n1, n4, n3\n");
    }
    my_move2(x4, y4);
    my_draw2(x1, y1);
    my_move2(x4, y4);
    my_draw2(x2, y2);
    my_move2(x4, y4);
    my_draw2(x3, y3);
    my_move2(x4, y4);
}

/*
 * swap a shared line between 2 adjacent elements
 */
void swapline(int gridno, int n1, int n2)
{
    int i, j, n[4], en1[3], en2[3], eq1[3], eq2[3], eqnum = 1;

    for (i = 0; i < 3; i++) {
	eq1[i] = 0;
	eq2[i] = 0;
	en1[i] = grid[gridno].icon[n1].nl[i];
	en2[i] = grid[gridno].icon[n2].nl[i];
    }
/*      find the 2 edges in common or return if none */
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    if (en1[i] == en2[j]) {
		eq1[i] = eqnum;
		eq2[j] = eqnum++;
		goto breakout;
	    }
	}
      breakout:;
    }
    if (eqnum != 3) {
	errwin("Elements do not share a common edge");
	return;
    }
    for (i = 0; i < 3; i++) {
	if (eq1[i] == 0) {
	    n[0] = en1[i];
	}
	if (eq2[i] == 0) {
	    n[1] = en2[i];
	}
	if (eq1[i] == 1) {
	    n[2] = en1[i];
	}
	if (eq1[i] == 2) {
	    n[3] = en1[i];
	}
    }
    iswap(&n[2], &n[3]);
    iswap(&n[0], &n[2]);
    iswap(&n[1], &n[3]);
    en1[0] = n[0];
    en1[1] = n[2];
    en1[2] = n[3];
    if (!comparea(gridno, en1[0], en1[1], en1[2])) {
	iswap(&en1[0], &en1[2]);
    }
    en2[0] = n[1];
    en2[1] = n[2];
    en2[2] = n[3];
    if (!comparea(gridno, en2[0], en2[1], en2[2])) {
	iswap(&en2[0], &en2[2]);
    }
    for (i = 0; i < 3; i++) {
	grid[gridno].icon[n1].nl[i] = en1[i];
	grid[gridno].icon[n2].nl[i] = en2[i];
    }
}

/*
 * convert a quadrangle to triangles
 */
void quad_to_tris(int gridno, int e)
{
    int a1, a2, a3, a4;
    a1 = grid[gridno].icon[e].nl[0];
    a2 = grid[gridno].icon[e].nl[1];
    a3 = grid[gridno].icon[e].nl[2];
    a4 = grid[gridno].icon[e].nl[3];
    add_element(gridno, a1, a2, a3);
    add_element(gridno, a1, a3, a4);
    delete_elements(gridno, &e, 1);
}

/*
 * loop through a grid and convert triangles to
 * quadrangles by looking at neighboring elements
 */
void tris_to_quad(int gridno, int e1, int e2)
{
    int i, j, w, n1[3], n2[3], p1, p2, q1, q2, addn, startn, edel[2];
    int a1, a2, a3, a4;

    for (i = 0; i < 3; i++) {
	n1[i] = grid[gridno].icon[e1].nl[i];
	n2[i] = grid[gridno].icon[e2].nl[i];
    }
/* find the shared edge */
    w = 1;
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    if (w == 1 && n1[i] == n2[j]) {
		w = 2;
		p1 = j;
		q1 = i;
	    } else if (w == 2 && n1[i] == n2[j]) {
		w = 3;
		p2 = j;
		q2 = i;
		goto breakout;
	    }
	}
    }
  breakout:;
    if (w != 3) {
	errwin("Elements do not share a common edge");
	return;
    }
    addn = -1;
    if (p1 == 0 && p2 == 1 || p1 == 1 && p2 == 0) {
	addn = 2;
    }
    if (p1 == 1 && p2 == 2 || p1 == 2 && p2 == 1) {
	addn = 0;
    }
    if (p1 == 0 && p2 == 2 || p1 == 2 && p2 == 0) {
	addn = 1;
    }
    if (addn < 0) {
	errwin("Major failure in tris_to_quad");
	return;
    }
    startn = -1;
    if (q1 == 0 && q2 == 1 || q1 == 1 && q2 == 0) {
	startn = 2;
    }
    if (q1 == 1 && q2 == 2 || q1 == 2 && q2 == 1) {
	startn = 0;
    }
    if (q1 == 0 && q2 == 2 || q1 == 2 && q2 == 0) {
	startn = 1;
    }
    if (startn < 0) {
	errwin("Major failure in tris_to_quad, can't find starting node");
	return;
    }
    addn = grid[gridno].icon[e2].nl[addn];
    edel[0] = e1;
    edel[1] = e2;
    switch (startn) {
    case 0:
	a1 = grid[gridno].icon[e1].nl[0];
	a2 = grid[gridno].icon[e1].nl[1];
	a3 = addn;
	a4 = grid[gridno].icon[e1].nl[2];
	break;
    case 1:
	a1 = grid[gridno].icon[e1].nl[1];
	a2 = grid[gridno].icon[e1].nl[2];
	a3 = addn;
	a4 = grid[gridno].icon[e1].nl[0];
	break;
    case 2:
	a1 = grid[gridno].icon[e1].nl[2];
	a2 = grid[gridno].icon[e1].nl[0];
	a3 = addn;
	a4 = grid[gridno].icon[e1].nl[1];
	break;
    }
/*
   printf("case %d: %d %d %d %d\n", startn, a1, a2, a3, a4);
 */
    add_quad_element(gridno, a1, a2, a3, a4);
    delete_elements(gridno, edel, 2);
}

/*
 * test for the sign of the area, test for orientation of
 * nodes
 */
int comparea(int gridno, int n1, int n2, int n3)
{
    double a[2], b[2], ar;

    a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
    a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
    b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
    b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    if (ar <= 0) {
	return 0;
    }
    return 1;
}

/*
 * compute area for 2d elements
 */
double area(int gridno, int elem)
{
    int n1, n2, n3;
    double a[2], b[2], ar;

    switch (grid[gridno].icon[elem].type) {
    case 3:
    case 6:
	n1 = grid[gridno].icon[elem].nl[0];
	n2 = grid[gridno].icon[elem].nl[1];
	n3 = grid[gridno].icon[elem].nl[2];
	a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
	a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
	b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
	b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
	ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
	break;
    case 4:
    case 8:
	n1 = grid[gridno].icon[elem].nl[0];
	n2 = grid[gridno].icon[elem].nl[1];
	n3 = grid[gridno].icon[elem].nl[2];
	a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
	a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
	b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
	b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
	ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
	n1 = grid[gridno].icon[elem].nl[0];
	n2 = grid[gridno].icon[elem].nl[2];
	n3 = grid[gridno].icon[elem].nl[3];
	a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
	a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
	b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
	b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
	ar = ar + 0.5 * (a[1] * b[0] - a[0] * b[1]);
	break;
    }
    return ar;
}

/*
 * compute the area given the nodes
 */
double area_from_nodes(int gridno, int n1, int n2, int n3)
{
    double a[2], b[2], ar;

    a[0] = grid[gridno].xord[n3] - grid[gridno].xord[n2];
    a[1] = grid[gridno].xord[n1] - grid[gridno].xord[n3];
    b[0] = grid[gridno].yord[n2] - grid[gridno].yord[n3];
    b[1] = grid[gridno].yord[n3] - grid[gridno].yord[n1];
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    return ar;
}

/*
 * compute the area given the vertices
 */
double area_from_vertices(double x1, double y1, double x2, double y2, double x3, double y3)
{
    double a[2], b[2], ar;

    a[0] = x3 - x2;
    a[1] = x1 - x3;
    b[0] = y2 - y3;
    b[1] = y3 - y1;
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    return ar;
}

/*
 * get the bounding box for a particular element
 */
void get_bounding_box(int gridno, int elem, double *bx1, double *bx2, double *by1, double *by2)
{
    int n1, n2, n3;
    double a[2], b[2], ar;

    n1 = grid[gridno].icon[elem].nl[0];
    n2 = grid[gridno].icon[elem].nl[1];
    n3 = grid[gridno].icon[elem].nl[2];
    *bx1 = grid[gridno].xord[n1];
    *bx2 = grid[gridno].xord[n1];
    *by1 = grid[gridno].yord[n1];
    *by2 = grid[gridno].yord[n1];
    if (*bx1 > grid[gridno].xord[n2]) {
	*bx1 = grid[gridno].xord[n2];
    }
    if (*bx1 > grid[gridno].xord[n3]) {
	*bx1 = grid[gridno].xord[n3];
    }
    if (*bx2 < grid[gridno].xord[n2]) {
	*bx2 = grid[gridno].xord[n2];
    }
    if (*bx2 < grid[gridno].xord[n3]) {
	*bx2 = grid[gridno].xord[n3];
    }
    if (*by1 > grid[gridno].yord[n2]) {
	*by1 = grid[gridno].yord[n2];
    }
    if (*by1 > grid[gridno].yord[n3]) {
	*by1 = grid[gridno].yord[n3];
    }
    if (*by2 < grid[gridno].yord[n2]) {
	*by2 = grid[gridno].yord[n2];
    }
    if (*by2 < grid[gridno].yord[n3]) {
	*by2 = grid[gridno].yord[n3];
    }
}

/*
 * move the node nearsest to oldx, oldy to newx, newy
 */
void do_move_node(int gridno, double oldx, double oldy, double newx, double newy)
{
    int ind;

    find_nearest_node(gridno, oldx, oldy, &ind);
    grid[gridno].xord[ind] = newx;
    grid[gridno].yord[ind] = newy;
    return;
}

/*
 * Find the nearset element and as a side effect, load item_buf 
 * with info about the element.
 */
void getelement(int gridno, double wx, double wy)
{
    int elem, n0, n1, n2;

    find_element(gridno, wx, wy, &elem);
    if (elem) {
	n0 = grid[gridno].icon[elem].nl[0];
	n1 = grid[gridno].icon[elem].nl[1];
	n2 = grid[gridno].icon[elem].nl[2];
	sprintf(item_buf, "(x, y, element, [n1, n2, n3]) = (%7.2lf, %7.2lf, %d, [%d %d %d])", wx, wy, elem, n0, n1, n2);
    } else {
	sprintf(item_buf, "(x, y, No element found)", wx, wy);
    }
}

/*
 * report on the nearest node to wx, wy
 * load item_buf with info about the node.
 */
void getnodes(int gridno, double wx, double wy)
{
    int ind;
    double px, py, d;

    find_nearest_node(gridno, wx, wy, &ind);
    if (ind >= 0) {
	px = grid[gridno].xord[ind];
	py = grid[gridno].yord[ind];
	d = grid[gridno].depth[ind];
	sprintf(item_buf, "(x,y), depth, node) = (%7.2lf,%7.2lf), %9.3lf, %d)", px, py, d, ind + 1);
    }
}

void getall(int gridno, double wx, double wy)
{
    int ind, elem;
    double px, py;

    find_nearest_node(gridno, wx, wy, &ind);
    find_element(gridno, wx, wy, &elem);
    sprintf(item_buf, "(node,element) = (%1d, %1d)", ind, elem);
}

void select_snap_line(double wx1, double wy1, double wx2, double wy2)
{
    sm = (wy2 - wy1) / (wx2 - wx1);
    sb = wy2 - sm * wx2;
    nibuf = 0;
    write_mode_str("Select nodes to snap with left button, middle button to register");
}

void snap_nodes(int gridno, int *ibuf, int nibuf)
{
    int i;
    double x, y, r, at, xinter;

    for (i = 0; i < nibuf; i++) {
	at = atan(sm);
	x = grid[gridno].xord[ibuf[i]];
	y = grid[gridno].yord[ibuf[i]];
	xinter = (y - sb) / sm;
	r = x - xinter;
	grid[gridno].xord[ibuf[i]] = xinter + r * cos(at);
	grid[gridno].yord[ibuf[i]] = y + r * sin(at);
    }
}

void do_snap(int gridno, double wx, double wy)
{
    int ind;

    solidbox(wx, wy);
    find_nearest_node(gridno, wx, wy, &ind);
    if (ind >= 0)
	ibuf[nibuf++] = ind;
}

void accumulate_nearest_elements(int gridno, double wx, double wy)
{
    int ind;
    double x, y;

    find_nearest_element(gridno, wx, wy, &ind);
    get_center(gridno, ind, &x, &y);
    solidbox(x, y);
    if (ind >= 0) {
	ibuf[nibuf++] = ind;
    }
}

void accumulate_nearest_nodes(int gridno, double wx, double wy)
{
    int ind;

    find_nearest_node(gridno, wx, wy, &ind);
    if (ind >= 0) {
	ibuf[nibuf++] = ind;
	wx = grid[curgrid].xord[ind];
	wy = grid[curgrid].yord[ind];
	solidbox(wx, wy);
    }
}

void accumulate_elements(int gridno, double wx, double wy)
{
    int ind;

    find_element(gridno, wx, wy, &ind);
    if (ind >= 0) {
	ibuf[nibuf++] = ind;
    }
}

void do_split_element3(int gridno, double wx, double wy)
{
    int ind, elem;

    find_element(gridno, wx, wy, &elem);
    if (elem >= 0) {
	split_elem3(gridno, elem);
    }
}

void do_split_element4(int gridno, double wx, double wy)
{
    int ind, elem;

    find_element(gridno, wx, wy, &elem);
    if (elem >= 0) {
	split_elem4(gridno, elem);
    }
}

void do_swap_line(int gridno, double wx1, double wy1, double wx2, double wy2)
{
    int i1, i2;

    find_element(gridno, wx1, wy1, &i1);
    find_element(gridno, wx2, wy2, &i2);
    swapline(gridno, i1, i2);
}

void get_center(int gridno, int elem, double *cx, double *cy)
{
    double x = 0.0, y = 0.0;
    int i, n;

    for (i = 0; i < grid[gridno].icon[elem].ngeom; i++) {
	n = grid[gridno].icon[elem].nl[i];
	x += grid[gridno].xord[n];
	y += grid[gridno].yord[n];
    }
    *cx = x / grid[gridno].icon[elem].ngeom;
    *cy = y / grid[gridno].icon[elem].ngeom;
}

/*
void get_center(int gridno, int i, double *cx, double *cy)
{
    double x1, y1, x2, y2, x3, y3, x4, y4;
    int n1, n2, n3, n4;

    switch (grid[gridno].icon[i].type) {
    case 3:
    case 6:
	n1 = grid[gridno].icon[i].nl[0];
	x1 = grid[gridno].xord[n1];
	y1 = grid[gridno].yord[n1];
	n2 = grid[gridno].icon[i].nl[1];
	x2 = grid[gridno].xord[n2];
	y2 = grid[gridno].yord[n2];
	n3 = grid[gridno].icon[i].nl[2];
	x3 = grid[gridno].xord[n3];
	y3 = grid[gridno].yord[n3];
	*cx = 0.33333333333333 * (x1 + x2 + x3);
	*cy = 0.33333333333333 * (y1 + y2 + y3);
	break;
    case 4:
    case 8:
	n1 = grid[gridno].icon[i].nl[0];
	x1 = grid[gridno].xord[n1];
	y1 = grid[gridno].yord[n1];
	n2 = grid[gridno].icon[i].nl[1];
	x2 = grid[gridno].xord[n2];
	y2 = grid[gridno].yord[n2];
	n3 = grid[gridno].icon[i].nl[2];
	x3 = grid[gridno].xord[n3];
	y3 = grid[gridno].yord[n3];
	n4 = grid[gridno].icon[i].nl[3];
	x4 = grid[gridno].xord[n4];
	y4 = grid[gridno].yord[n4];
	*cx = 0.25 * (x1 + x2 + x3 + x4);
	*cy = 0.25 * (y1 + y2 + y3 + y4);
	break;
    }
}
*/

void get_average_depth(int gridno, int elem, double *depth)
{
    double d = 0.0;
    int i, n;

    for (i = 0; i < grid[gridno].icon[elem].type; i++) {
	n = grid[gridno].icon[elem].nl[i];
	d += grid[gridno].depth[n];
    }
    *depth = d / grid[gridno].icon[elem].ngeom;
}

void copy_grid(int from, int to)
{
    int i, j, ib;

    Free_grid(to);
    Allocate_grid(to, grid[from].nmel, grid[from].nmnp);
    for (i = 0; i < grid[from].nmnp; i++) {
	grid[to].xord[i] = grid[from].xord[i];
	grid[to].yord[i] = grid[from].yord[i];
	grid[to].depth[i] = grid[from].depth[i];
    }
    grid[to].gridtype = grid[from].gridtype;
    for (i = 0; i < grid[from].nmel; i++) {
	grid[to].icon[i] = grid[from].icon[i];
    }
    grid[to].xmin = grid[from].xmin;
    grid[to].xmax = grid[from].xmax;
    grid[to].ymin = grid[from].ymin;
    grid[to].ymax = grid[from].ymax;
    grid[to].dmin = grid[from].dmin;
    grid[to].dmax = grid[from].dmax;
    if (grid[from].nbounds) {
	grid[to].nbounds = grid[from].nbounds;
	for (j = 0; j < grid[from].nbounds; j++) {
	    ib = grid[from].boundaries[j];
	    grid[to].boundaries[j] = copy_boundary(ib);
	}
    }
    default_isolines(&grid[to].ip);
}

void do_interp_depths(void)
{
    write_mode_str("Interpolating depths from background grid, please wait...");
    interpdepths(MAXGRIDS, curgrid);
    write_mode_str(NULL);
}

/*
 * interpolate depths defined in grid 1 (g1) to grid 2 (g2)
 *
 */
void interpdepths(int g1, int g2)
{
    char buf[256];
    int i, ind, err;
    double coefs[3], x, y;
    extern int interp_warn_level;
    extern int interp_in_region;
    extern int noskip;

    if (interp_in_region && !region_flag) {
	errwin("No region defined, operation cancelled");
	return;
    }
    prep_find(g1);
    err = 0;
    if (interp_warn_level) {
	sprintf(buf, "\nDepth interpolation: reporting on nodes outside background grid\n");
	stufftext(buf, 1);
    }
    for (i = 0; i < grid[g2].nmnp; i++) {
	x = grid[g2].xord[i];
	y = grid[g2].yord[i];
	if (interp_in_region && region_flag && (!inregion(regionx, regiony, nregion, x, y))) {
	    continue;
	}
	ind = find_element3(g1, x, y);
	if (ind < 0) {
	    if (noskip) {
		double close;
		double tmp;
		int j, closenode;
		close = hypot(grid[g1].xord[0] - x, grid[g1].yord[0] - y);
		for (j = 1; j < grid[g1].nmnp; j++) {
		    tmp = hypot(grid[g1].xord[j] - x, grid[g1].yord[j] - y);
		    if (close > tmp) {
			close = tmp;
			closenode = j;
		    }
		}
/*
   sprintf(buf, "   Node %d assigned a depth using node %d\n", i + 1, closenode);
   stufftext(buf, 0);
 */
		grid[g2].depth[i] = grid[g1].depth[closenode];
	    } else {
		if (interp_warn_level == 1) {
		    sprintf(buf, "   Node %d not interpolated\n", i + 1);
		    stufftext(buf, 0);
		    err++;
		} else if (interp_warn_level == 2) {
		    err++;
		}
	    }
	} else {
	    grid[g2].depth[i] = get_depth_element(g1, ind, x, y);
	}
    }

    if (interp_warn_level) {
	if (err == 0) {
	    sprintf(buf, "All nodes were interpolated\n");
	} else if (err == 1) {
	    sprintf(buf, "Depth interpolation: 1 node not interpolated\n");
	} else {
	    sprintf(buf, "Depth interpolation: %d nodes were not interpolated\n", err);
	}
	stufftext(buf, 2);
    }
    end_find();
}

double get_depth_element(int gno, int ind, double x, double y)
{
    double x1, x2, x3, y1, y2, y3, x4, y4;
    double h, r, s, area2, a[3], b[3], g[3];
    int n0, n1, n2, n3;

    if (grid[gno].depthflag) {
	if (ind >= 0 && ind < grid[gno].nmel && grid[gno].edepth != NULL) {
	    return grid[gno].edepth[ind];
	} else {
	}
    }
    n0 = grid[gno].icon[ind].nl[0];
    n1 = grid[gno].icon[ind].nl[1];
    n2 = grid[gno].icon[ind].nl[2];
    x1 = grid[gno].xord[n0];
    x2 = grid[gno].xord[n1];
    x3 = grid[gno].xord[n2];
    y1 = grid[gno].yord[n0];
    y2 = grid[gno].yord[n1];
    y3 = grid[gno].yord[n2];
    if (grid[gno].icon[ind].ngeom == 3) {
	a[0] = x2 * y3 - x3 * y2;
	a[1] = x3 * y1 - x1 * y3;
	a[2] = x1 * y2 - x2 * y1;
	b[0] = x2 - y3;
	b[1] = y3 - y1;
	b[2] = y1 - y2;
	g[0] = x3 - x2;
	g[1] = x1 - x3;
	g[2] = x2 - x1;
	area2 = a[0] + a[1] + a[2];
	s = (a[2] + b[2] * x + g[2] * y) / area2;
	r = (a[1] + b[1] * x + g[1] * y) / area2;
	h = 1.0 - r - s;
	return (grid[gno].depth[n0] * h + grid[gno].depth[n1] * r + grid[gno].depth[n2] * s);
    } else if (grid[gno].icon[ind].ngeom == 4) {
	double d = 0.0, xi, eta, w[4];
	int i;
	n3 = grid[gno].icon[ind].nl[3];
	x4 = grid[gno].xord[n3];
	y4 = grid[gno].yord[n3];
	ibilinear(x1, x2, x3, x4, y1, y2, y3, y4, x, y, &xi, &eta, w);
	/*
	   printf("       x1 = %lf\n       x2 = %lf\n       x3 = %lf\n       x4 = %lf\n       y1 = %lf\n       y2 = %lf\n        y3 = %lf\n       y4 = %lf\n       x = %lf\n       y = %lf\n", x1, x2, x3, x4, y1, y2, y3, y4, x, y);
	 */
	for (i = 0; i < 4; i++) {
	    d += w[i] * grid[gno].depth[grid[gno].icon[ind].nl[i]];
	    /*
	       printf("%lf %lf %lf\n", d, w[i], grid[gno].depth[grid[gno].icon[ind].nl[i]]);
	     */
	}
	return d;
    } else {
	return -999999.0;
    }
}

/*
double get_depth_element(int gno, int ind, double x, double y)
{
    double x1, x2, x3, x4, y1, y2, y3, y4, xi, eta;
    double h, r, s, area2, a[4], b[4], g[4], shapef[4];
    int n0, n1, n2, n3;

    if (grid[gno].depthflag) {
	if (ind >= 0 && ind < grid[gno].nmel && grid[gno].edepth != NULL) {
	    return grid[gno].edepth[ind];
	} else {
	}
    }
    if (grid[gno].icon[ind].nn = 3) {

	n0 = grid[gno].icon[ind].nl[0];
	n1 = grid[gno].icon[ind].nl[1];
	n2 = grid[gno].icon[ind].nl[2];
	x1 = grid[gno].xord[n0];
	x2 = grid[gno].xord[n1];
	x3 = grid[gno].xord[n2];
	y1 = grid[gno].yord[n0];
	y2 = grid[gno].yord[n1];
	y3 = grid[gno].yord[n2];
	a[0] = x2 * y3 - x3 * y2;
	a[1] = x3 * y1 - x1 * y3;
	a[2] = x1 * y2 - x2 * y1;
	b[0] = x2 - y3;
	b[1] = y3 - y1;
	b[2] = y1 - y2;
	g[0] = x3 - x2;
	g[1] = x1 - x3;
	g[2] = x2 - x1;
	area2 = a[0] + a[1] + a[2];
	area2 = 1.0 / area2;
	s = (a[2] + b[2] * x + g[2] * y) * area2;
	r = (a[1] + b[1] * x + g[1] * y) * area2;
	h = 1.0 - r - s;
	return (grid[gno].depth[n0] * h + grid[gno].depth[n1] * r + grid[gno].depth[n2] * s);
    } else if (grid[gno].icon[ind].nn = 4) {

	n0 = grid[gno].icon[ind].nl[0];
	n1 = grid[gno].icon[ind].nl[1];
	n2 = grid[gno].icon[ind].nl[2];
	n3 = grid[gno].icon[ind].nl[3];
	x1 = grid[gno].xord[n0];
	x2 = grid[gno].xord[n1];
	x3 = grid[gno].xord[n2];
	x4 = grid[gno].xord[n3];
	y1 = grid[gno].yord[n0];
	y2 = grid[gno].yord[n1];
	y3 = grid[gno].yord[n2];
	y4 = grid[gno].yord[n3];
	ibilinear(x1, x2, x3, x4, y1, y2, y3, y4, x, y, &xi, &eta, shapef);
	return (grid[gno].depth[n0] * shapef[0] + grid[gno].depth[n1] * shapef[1] + grid[gno].depth[n2] * shapef[2] + grid[gno].depth[n3] * shapef[3]);
    } else {
	fprintf(stderr, "Error, unknown element type %d in get_depth_element() setting depth to -1.0\n", ind);
	return -1.0;
    }
}
*/

void find_assoc_elements(int gridno, int node, int *nel, int *ellist)
{
    int i, j, itmp = 0;

    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < 3; j++) {
	    if (grid[gridno].icon[i].nl[j] == node) {
		ellist[itmp] = i;
		itmp++;
	    }
	}
    }
    *nel = itmp;
}

void compute_nodel_con(int gridno)
{
    int i, j, itmp, ntmp;

    for (i = 0; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].nodecon[i].el != NULL) {
	    free(grid[gridno].nodecon[i].el);
	}
	grid[gridno].nodecon[i].el = NULL;
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < 3; j++) {
	    ntmp = grid[gridno].icon[i].nl[j];
	    if (grid[gridno].nodecon[ntmp].el == NULL) {
		grid[gridno].nodecon[ntmp].el = (int *) calloc(1, sizeof(int));
		grid[gridno].nodecon[ntmp].cnt = 1;
		grid[gridno].nodecon[ntmp].el[0] = i;
	    } else {
		itmp = grid[gridno].nodecon[ntmp].cnt;
		grid[gridno].nodecon[ntmp].el = (int *) realloc(grid[gridno].nodecon[ntmp].el, (itmp + 1) * sizeof(int));
		grid[gridno].nodecon[ntmp].cnt = itmp + 1;
		grid[gridno].nodecon[ntmp].el[itmp] = i;
	    }
	}
    }
}

int gengrid(int gridno, int bno)
{
    char buf[256];
    int i, j, itmp;
    double xmin, xmax, ymin, ymax, dmin, dmax;

    Free_grid(gridno);
    strcpy(grid[gridno].alphid, "Generated grid");
    grid[gridno].nmel = 2;
    grid[gridno].nmnp = 4;
    if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	return 0;
    }
    build_minmax(bno, &xmin, &xmax, &ymin, &ymax, &dmin, &dmax);
    grid[gridno].xmin = xmin;
    grid[gridno].xmax = xmax;
    grid[gridno].ymin = ymin;
    grid[gridno].ymax = ymax;
    grid[gridno].dmin = dmin;
    grid[gridno].dmax = dmax;
    grid[gridno].xord[0] = xmin;
    grid[gridno].xord[1] = xmax;
    grid[gridno].xord[2] = xmax;
    grid[gridno].xord[3] = xmin;
    grid[gridno].yord[0] = ymin;
    grid[gridno].yord[1] = ymin;
    grid[gridno].yord[2] = ymax;
    grid[gridno].yord[3] = ymax;
    grid[gridno].depth[0] = dmin;
    grid[gridno].depth[1] = dmin;
    grid[gridno].depth[2] = dmax;
    grid[gridno].depth[3] = dmax;
    grid[gridno].icon[0].type = 3;
    grid[gridno].icon[0].nn = 3;
    grid[gridno].icon[0].ngeom = 3;
    grid[gridno].icon[0].nl[0] = 0;
    grid[gridno].icon[0].nl[1] = 1;
    grid[gridno].icon[0].nl[2] = 3;
    grid[gridno].icon[1].type = 3;
    grid[gridno].icon[1].nn = 3;
    grid[gridno].icon[1].ngeom = 3;
    grid[gridno].icon[1].nl[0] = 1;
    grid[gridno].icon[1].nl[1] = 2;
    grid[gridno].icon[1].nl[2] = 3;
    return 1;
}

static int *minx, *maxx;
static int *miny, *maxy;
static int *ix1;

static double *gx, *gy, *gd;
static Element *icon;

static int comparemin(const void *x1ptr, const void *x2ptr);
int my_nmel;

/*
 */

int find_element3(int gno, double x, double y)
{
    int mvforward, endel, i, nel, ac, nmel = grid[gno].nmel;
    double pmintest, pmaxtest, x1, y1, x2, y2, x3, y3;

    x1 = y1 = x2 = y2 = x3 = y3 = 0.0;

/* get the lub of the minimums wrt to the min array */
    my_nmel = nmel;
    endel = findfxmin(x);
    if (endel >= nmel) {
	endel = nmel - 1;
    }
    if (endel == -1) {
	fprintf(stderr, "Fatal error in find_element3(): endel == -1\n");
    }
    mvforward = endel;
    while (mvforward < nmel && gx[minx[ix1[mvforward]]] < x)
	mvforward++;
    if (endel != mvforward) {
	endel = mvforward;
    }
    if (endel >= nmel) {
	endel = nmel - 1;
    }
    for (i = endel; i >= 0; i--) {
/* check for containment in y using the max array */
	if (y >= gy[miny[ix1[i]]] && y <= gy[maxy[ix1[i]]]) {
	    nel = ix1[i];
	    x1 = gx[icon[nel].nl[0]];
	    x2 = gx[icon[nel].nl[1]];
	    x3 = gx[icon[nel].nl[2]];
	    y1 = gy[icon[nel].nl[0]];
	    y2 = gy[icon[nel].nl[1]];
	    y3 = gy[icon[nel].nl[2]];
	    if (belel(x, y, x1, y1, x2, y2, x3, y3)) {
		return nel;
	    }
	}
    }
    return -1;
}

static int comparemin(const void *x1ptr, const void *x2ptr)
{
    int *x1 = (int *) x1ptr;
    int *x2 = (int *) x2ptr;
    if (gx[minx[*x1]] < gx[minx[*x2]])
	return -1;
    if (gx[minx[*x1]] > gx[minx[*x2]])
	return 1;
    return 0;
}

static int findfxmin(double s)
{

    int low, high, mid;

    low = 0;
    high = my_nmel - 1;
    if (s < gx[minx[ix1[low]]]) {
	return low;
    } else if (s > gx[minx[ix1[high]]]) {
	return high;
    }
    while (low <= high) {
	mid = (low + high) / 2;
	if (s < gx[minx[ix1[mid]]] && s < gx[minx[ix1[mid - 1]]]) {
	    high = mid - 1;
	} else {
	    if (s > gx[minx[ix1[mid]]] && s > gx[minx[ix1[mid + 1]]]) {
		low = mid + 1;
	    } else {
		return mid;
	    }
	}
    }
    return (-1);
}

void end_find(void)
{
    free(minx);
    free(miny);
    free(maxy);
    free(ix1);
    minx = NULL;
    miny = NULL;
    maxy = NULL;
    ix1 = NULL;
}

void prep_find(int gno)
{
    int i, nmel = grid[gno].nmel;
    double pmintest, pmaxtest;

    gx = grid[gno].xord;
    gy = grid[gno].yord;
    gd = grid[gno].depth;
    icon = grid[gno].icon;
    my_nmel = grid[gno].nmel;

/* local storage allocation */
    minx = (int *) calloc(nmel, sizeof(int));
    miny = (int *) calloc(nmel, sizeof(int));
    maxy = (int *) calloc(nmel, sizeof(int));
    ix1 = (int *) calloc(nmel, sizeof(int));

/* find min and max node numbers for each element */
    for (i = 0; i < nmel; i++) {
	ix1[i] = i;
	pmintest = gx[icon[i].nl[0]];
	minx[i] = icon[i].nl[0];
	if (gx[icon[i].nl[1]] < pmintest) {
	    pmintest = gx[icon[i].nl[1]];
	    minx[i] = icon[i].nl[1];
	}
	if (gx[icon[i].nl[2]] < pmintest) {
	    minx[i] = icon[i].nl[2];
	}
	pmintest = pmaxtest = gy[icon[i].nl[0]];
	miny[i] = maxy[i] = icon[i].nl[0];
	if (gy[icon[i].nl[1]] < pmintest) {
	    pmintest = gy[icon[i].nl[1]];
	    miny[i] = icon[i].nl[1];
	}
	if (gy[icon[i].nl[1]] > pmaxtest) {
	    pmaxtest = gy[icon[i].nl[1]];
	    maxy[i] = icon[i].nl[1];
	}
	if (gy[icon[i].nl[2]] < pmintest) {
	    miny[i] = icon[i].nl[2];
	}
	if (gy[icon[i].nl[2]] > pmaxtest) {
	    maxy[i] = icon[i].nl[2];
	}
	if (miny[i] > grid[gno].nmnp) {
	    printf("ERROR %d %d\n", i, miny[i]);
	}
    }
    qsort(ix1, nmel, sizeof(int), comparemin);
}

/*
 * routines to reduce bandwidth
 *
 */
void renum(int gridno, int *idiff, int *ndiff)
{
    int *jnt, *jt, *memjt, *jmem;
    int i, j, ii, n1, n2, n3, nmnp = grid[gridno].nmnp, nmel = grid[gridno].nmel;
    double *xtmp, *ytmp, *dtmp;

    xtmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));
    ytmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));
    dtmp = (double *) calloc(grid[gridno].nmnp, sizeof(double));

    for (i = 0; i < nmnp; i++) {
	xtmp[i] = grid[gridno].xord[i];
	ytmp[i] = grid[gridno].yord[i];
	dtmp[i] = grid[gridno].depth[i];
    }
    jnt = (int *) calloc(4 * nmel, sizeof(int));
    jt = (int *) calloc(4 * nmel, sizeof(int));
    memjt = (int *) calloc(8 * nmel, sizeof(int));
    jmem = (int *) calloc(nmel, sizeof(int));
    n1 = nmel;
    n2 = nmel * 2;
    n3 = nmel * 3;
    for (i = 0; i < nmel; i++) {
	jt[i] = grid[gridno].icon[i].nl[0] + 1;
	jt[i + n1] = grid[gridno].icon[i].nl[1] + 1;
	jt[i + n2] = grid[gridno].icon[i].nl[2] + 1;
	jt[i + n3] = 0;
    }
    *idiff = *ndiff = 0;
    setup(nmnp, nmel, jt, memjt, jmem, jnt, idiff, ndiff);
    optnum(nmnp, nmel, jt, memjt, jmem, jnt, idiff, ndiff);
    if (*idiff == *ndiff) {
	errwin("Grid bandwidth already optimized");
    } else {
	for (ii = 0; ii < nmnp; ii++) {
	    for (i = 0; i < nmnp; i++) {
		j = jnt[i];
		if (ii != j - 1) {
		    continue;
		}
		grid[gridno].xord[ii] = xtmp[i];
		grid[gridno].yord[ii] = ytmp[i];
		grid[gridno].depth[ii] = dtmp[i];
	    }
	}
	for (i = 0; i < nmel; i++) {
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    grid[gridno].icon[i].type = 3;
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	    grid[gridno].icon[i].nl[0] = jnt[n1] - 1;
	    grid[gridno].icon[i].nl[1] = jnt[n2] - 1;
	    grid[gridno].icon[i].nl[2] = jnt[n3] - 1;
/*
   printf("%d %d %d - %d %d %d\n", n1, n2, n3, jnt[n1] - 1, jnt[n2] - 1, jnt[n3] - 1);
 */
	}
    }
    free(jnt);
    free(jt);
    free(memjt);
    free(jmem);
    free(xtmp);
    free(ytmp);
    free(dtmp);
}

/* greditrenu2.f -- translated by f2c (version of 10 December 1991  17:07:07).
   You must link the resulting object file with the libraries:
   -lF77 -lI77 -lm -lc   (in that order)
 */

#include "f2c.h"

 /* Subroutine */ int setup(integer nodes, integer lments, integer * jt,
			    integer * memjt, integer * jmem, integer * jnt, integer * idiff, integer * ndiff)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer jsub, jnti, i, j, ii, iii, jjt, mem1;

    /* Parameter adjustments */
    --jnt;
    --jmem;
    --memjt;
    --jt;

    *idiff = 0;
    i__1 = nodes;
    for (j = 1; j <= i__1; ++j) {
	jmem[j] = 0;
/* L10: */
    }
    i__1 = lments;
    for (j = 1; j <= i__1; ++j) {
	for (i = 1; i <= 4; ++i) {
	    jnti = jt[lments * (i - 1) + j];
	    if (jnti == 0) {
		goto L60;
	    }
	    jsub = (jnti - 1) << 3;
	    for (ii = 1; ii <= 4; ++ii) {
		if (ii == i) {
		    goto L40;
		}
		jjt = jt[lments * (ii - 1) + j];
		if (jjt == 0) {
		    goto L50;
		}
		mem1 = jmem[jnti];
		if (mem1 == 0) {
		    goto L30;
		}
		i__2 = mem1;
		for (iii = 1; iii <= i__2; ++iii) {
		    if (memjt[jsub + iii] == jjt) {
			goto L40;
		    }
/* L20: */
		}
	      L30:
		++jmem[jnti];
		memjt[jsub + jmem[jnti]] = jjt;
		if ((i__2 = jnti - jjt, abs(i__2)) > *idiff) {
		    *idiff = (i__3 = jnti - jjt, abs(i__3));
		}
	      L40:
		;
	    }
	  L50:
	    ;
	}
      L60:
	;
    }
    return 0;
}				/* setup_ */

 /* Subroutine */ int optnum(int nodes, int lments, int *jt, int *memjt,
			     int *jmem, int *jnt, int *idiff, int *ndiff)
{
    /* System generated locals */
    int i__1, i__2, i__3;

    /* Local variables */
    static int jsub, njts, i, j, k, *joint, *newjt, k4, k5, ik, jj, nmdiff, minmax, max_;

    joint = (int *) calloc(2 * nodes, sizeof(int));
    newjt = (int *) calloc(2 * nodes, sizeof(int));
/* cccc */
    /* Parameter adjustments */
    --jnt;
    --jmem;
    --memjt;
    --jt;

    /* Function Body */
    njts = nodes;
    minmax = *idiff;
    i__1 = nodes;
    for (ik = 1; ik <= i__1; ++ik) {
	i__2 = nodes;
	for (j = 1; j <= i__2; ++j) {
	    joint[j - 1] = 0;
/* L20: */
	    newjt[j - 1] = 0;
	}
	max_ = 0;
	i = 1;
	newjt[0] = ik;
	joint[ik - 1] = 1;
	k = 1;
      L30:
	k4 = jmem[newjt[i - 1]];
	if (k4 == 0) {
	    goto L45;
	}
	jsub = (newjt[i - 1] - 1) << 3;
	i__2 = k4;
	for (jj = 1; jj <= i__2; ++jj) {
	    k5 = memjt[jsub + jj];
	    if (joint[k5 - 1] > 0) {
		goto L40;
	    }
	    ++k;
	    newjt[k - 1] = k5;
	    joint[k5 - 1] = k;
	    *ndiff = (i__3 = i - k, abs(i__3));
	    if (*ndiff >= minmax) {
		goto L60;
	    }
	    if (*ndiff > max_) {
		max_ = *ndiff;
	    }
	  L40:
	    ;
	}
	if (k == njts) {
	    goto L50;
	}
      L45:
	++i;
	goto L30;
      L50:
	minmax = max_;
	i__2 = nodes;
	for (j = 1; j <= i__2; ++j) {
/* L55: */
	    jnt[j] = joint[j - 1];
	}
      L60:
	;
    }
    nmdiff = *ndiff + 1;
    free(joint);
    free(newjt);
    return 0;
}				/* optnum_ */

struct edgelist {
    int *list;
} *edgel, *connecl, *adjn, *adjel;

struct edge {
    int n1, n2;
    int el1, el2;
    int pos1, pos2;
    int visited;
};

void free_edgelist(int gridno, struct edgelist *e)
{
    int i, nmnp;

    nmnp = grid[gridno].nmnp;
    for (i = 0; i < nmnp; i++) {
	if (e[i].list != NULL) {
	    free(e[i].list);
	}
    }
    free(e);
}

struct edgelist2 {
    int n;
    struct edge *edges;
} elist;

struct elem {
    int n;
    int e[3];
};

void getedgelist2(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[3], v, nmnp, nmel, edgecount = 0;
    struct elem *gel;
    FILE *fp = fopen("grid.tmp", "w");
    nmnp = grid[gridno].nmnp;
    nmel = grid[gridno].nmel;
    gel = (struct elem *) malloc(nmel * sizeof(struct elem));
    elist.edges = (struct edge *) malloc(nmnp * (nmel - 1) * sizeof(struct edge));
    elist.n = 0;
    for (i = 0; i < nmel; i++) {
	for (j = 0; j < grid[gridno].icon[i].ngeom; j++) {
	    if (elist.n > nmnp * (nmel - 1)) {
		fprintf(stderr, "Too many edges\n");
		exit(1);
	    }
	    n1 = grid[gridno].icon[i].nl[j];
	    n2 = grid[gridno].icon[i].nl[(j + 1) % grid[gridno].icon[i].ngeom];
	    if (n1 > n2) {
		int tmp = n1;
		n1 = n2;
		n2 = tmp;
	    }
	    elist.edges[elist.n].n1 = n1;
	    elist.edges[elist.n].n2 = n2;
	    elist.edges[elist.n].el1 = i;
	    gel[i].e[gel[i].n] = elist.n;
	    gel[i].n++;
	    elist.n++;
	}
    }
    fprintf(fp, "%d %d %d\n", nmnp, elist.n, nmel);
    for (i = 0; i < nmnp; i++) {
	fprintf(fp, "%d %lf %lf %lf\n", i, grid[gridno].xord[i], grid[gridno].yord[i], grid[gridno].depth[i]);
    }
    for (i = 0; i < elist.n; i++) {
	fprintf(fp, "%d %d %d\n", i, elist.edges[i].n1, elist.edges[i].n2);
    }
    for (i = 0; i < nmel; i++) {
	fprintf(fp, "%d %d %d %d\n", i, gel[i].e[0], gel[i].e[1], gel[i].e[2]);

    }
    fclose(fp);
    free(elist.edges);
    free(gel);
}

/*
   void writeedgelist2(int gridno, char *fn)
   {
   getedgelist2(gridno);
   for (i=0;i<elist.n;i++) {
   fprintf(fp, "%d %d %d\n", i, elist.edges[i].n1, elist.edges[i].n2);
   }
   for (i=0;i<grid[gridno].nmel;i++) {
   fprintf(fp, "%d %d %d %d\n", i, elist.edges[i].n1, elist.edges[i].n2);
   }
   }
 */

void getedgelist(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[3], v, nmnp;
    struct edgelist *e;

    nmnp = grid[gridno].nmnp;
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) calloc(2, sizeof(int));
	e[i].list[0] = -1;
	e[i].list[1] = 0;
    }

    for (i = 0; i < grid[gridno].nmel; i++) {
	el[0] = grid[gridno].icon[i].nl[0];
	el[1] = grid[gridno].icon[i].nl[1];
	el[2] = grid[gridno].icon[i].nl[2];

	for (j = 0; j < 3; j++) {
/* search for edge */
	    cnt1 = 0;
	    cnt2 = 0;
	    f1 = 0;
	    f2 = 0;
	    n1 = el[j];
	    n2 = el[(j + 1) % 3];

	    while (e[n1].list[2 * cnt1] != -1 && e[n1].list[2 * cnt1] != n2) {
		cnt1++;
	    }
	    while (e[n2].list[2 * cnt2] != -1 && e[n2].list[2 * cnt2] != n1) {
		cnt2++;
	    }

	    if (e[n1].list[2 * cnt1] == -1) {
		f1 = 0;
	    } else {
		f1 = 1;
	    }
	    if (e[n2].list[2 * cnt2] == -1) {
		f2 = 0;
	    } else {
		f2 = 1;
	    }

/* check cnt1 and cnt2, should only appear once */
	    if (f1 && f2) {
		printf("Error: duplicate edge, send mail to ace_bugs\n");
		exit(1);
	    }
/* if cnt1 != 0 then increment element counter */
	    if (f1 && !f2) {
		e[n1].list[2 * cnt1 + 1]++;
	    }
/* if cnt2 != 0 then increment element counter */
	    else if (f2 && !f1) {
		e[n2].list[2 * cnt2 + 1]++;
	    }
/* if cnt1 == 0 && cnt2 == 0 then add edge n1-n2 and increment counter */
	    else if (!f1 && !f2) {
		e[n1].list = (int *) realloc(e[n1].list, 2 * (cnt1 + 2) * sizeof(int));
		e[n1].list[2 * cnt1] = n2;
		e[n1].list[2 * cnt1 + 1]++;
		cnt1++;
		e[n1].list[2 * cnt1] = -1;
		e[n1].list[2 * cnt1 + 1] = 0;
	    } else {
		printf("Error: duplicate edge (phase 2), send mail to ace_bugs\n");
		exit(1);
	    }
	}
    }
    edgel = e;
}

void getalledgelist(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[8], v, nmnp;
    struct edgelist *e;

    nmnp = grid[gridno].nmnp;
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
/*
 * initilize edgelist
 */
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) calloc(3, sizeof(int));
	e[i].list[0] = -1;
	e[i].list[1] = 0;
	e[i].list[2] = -1;
    }
/*
 * construct edgelist - this one will have duplicate edges
 */
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!(grid[gridno].icon[i].type == 3 || grid[gridno].icon[i].type == 4)) {
	    continue;
	}
	for (j = 0; j < grid[gridno].icon[i].ngeom; j++) {
	    el[j] = grid[gridno].icon[i].nl[j];
	}

	for (j = 0; j < grid[gridno].icon[i].ngeom; j++) {
	    cnt1 = 0;
/* loop through the edges of the i'th element */
	    n1 = el[j];
	    n2 = el[(j + 1) % grid[gridno].icon[i].ngeom];
/* find new spot */
	    while (e[n1].list[3 * cnt1] != -1) {
		cnt1++;
	    }
/* add this edge and save element number */
	    e[n1].list = (int *) realloc(e[n1].list, 3 * (cnt1 + 3) * sizeof(int));
/* other vertex of this edge */
	    e[n1].list[3 * cnt1] = n2;
/* element number */
	    e[n1].list[3 * cnt1 + 1] = i;
/* position of edge in element, 0->2 in the case of triangles */
	    e[n1].list[3 * cnt1 + 2] = j;
	    cnt1++;
/* terminate the edge list for this node */
	    e[n1].list[3 * cnt1] = -1;
	    e[n1].list[3 * cnt1 + 1] = 0;
	    e[n1].list[3 * cnt1 + 2] = -1;
	}
    }
    adjel = e;
}

static int *adjels[4];

void getassoc_elements(int gridno)
{
    int i, j, cnt1, cnt2, n1, n2, e1, e2, p1, p2;
    adjels[0] = (int *) calloc(grid[gridno].nmel, sizeof(int));
    adjels[1] = (int *) calloc(grid[gridno].nmel, sizeof(int));
    adjels[2] = (int *) calloc(grid[gridno].nmel, sizeof(int));
    adjels[3] = (int *) calloc(grid[gridno].nmel, sizeof(int));
    for (i = 0; i < grid[gridno].nmel; i++) {
	adjels[0][i] = -1;
	adjels[1][i] = -1;
	adjels[2][i] = -1;
	adjels[3][i] = -1;
    }
    getalledgelist(gridno);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	cnt1 = 0;
	while (adjel[i].list[3 * cnt1] != -1) {
/* first vertext */
	    n1 = i;
/* second vertext */
	    n2 = adjel[n1].list[3 * cnt1];
/* this element */
	    e1 = adjel[n1].list[3 * cnt1 + 1];
/* position of edge in element */
	    p1 = adjel[n1].list[3 * cnt1 + 2];
	    e2 = -1;
	    cnt2 = 0;
/* look in the list of vertext n2 to see if n1 appears, if it does,
   then this edge has a neighbor, if not, then it is a boundary edge */
	    while (adjel[n2].list[3 * cnt2] != -1 && adjel[n2].list[3 * cnt2] != n1) {
		cnt2++;
	    }
/* if cnt2 points to a legitimate node, then we have a neighbor */
	    if (adjel[n2].list[3 * cnt2] == n1) {
		e2 = adjel[n2].list[3 * cnt2 + 1];
		p2 = adjel[n2].list[3 * cnt2 + 2];
	    }
/* if e2 points to a legitimate element, then add it to the list of neighbors */
	    if (e2 >= 0) {
		adjels[p2][e2] = e1;
	    }
	    cnt1++;
	}
    }
}

int findsharednodes(int gridno, int e1, int e2, int *nn, int *or)
{
    int nn1[3], nn2[3];
    int n1 = -1, n2 = -1;
    int i, j;
    for (i = 0; i < 3; i++) {
	nn1[i] = grid[gridno].icon[e1].nl[i];
	nn2[i] = grid[gridno].icon[e2].nl[i];
    }
    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++) {
	    if (nn1[i] == nn2[j]) {
		if (n1 < 0) {
		    n1 = nn1[i];
		    or[0] = i;
		} else if (n2 < 0) {
		    n2 = nn2[j];
		    nn[0] = n1;
		    nn[1] = n2;
		    or[1] = j;
		    return 0;
		}
	    }
	}
    }
    return -1;
}

void find_sh(int gridno)
{
    int i, j, k, bel;
    int n1, n2, n3, n4, nn[2], or[2];
    int nn1[3], nn2[3];
    double mid1[3], mid2[3], midx, midy, qx[5], qy[5];
    getassoc_elements(gridno);
    for (i = 0; i < grid[gridno].nmel; i++) {
	bel = -1;
	for (j = 0; j < 3; j++) {
	    int el = adjels[j][i];
	    if (el < 0) {
		bel = j;
		break;
	    }
	}
// only if a boundary element
	if (bel >= 0) {
	    int ns[8];
	    for (j = 0; j < 3; j++) {
		int el = adjels[j][i];
                int ipos = -1;
		//if (el >= 0 && (el == 142 || el == 145)) {
		//if (el >= 0 && (i + 1 == 133 && el + 1 == 131)) {
		if (el >= 0) {
		    for (k=0;k<3;k++) {
			if (adjels[k][el] == i) {
			    ipos = k;
			    break;
			}
		    }
		    if (j == 0) {
			nn[0] = grid[gridno].icon[i].nl[0];
			nn[1] = grid[gridno].icon[i].nl[1];
			or[0] = 0;
			n1 = grid[gridno].icon[i].nl[1];
			n2 = grid[gridno].icon[i].nl[2];
			n3 = grid[gridno].icon[i].nl[2];
			n4 = grid[gridno].icon[i].nl[0];
		    } else if (j == 1) {
			nn[0] = grid[gridno].icon[i].nl[1];
			nn[1] = grid[gridno].icon[i].nl[2];
			or[0] = 1;
			n1 = grid[gridno].icon[i].nl[2];
			n2 = grid[gridno].icon[i].nl[0];
			n3 = grid[gridno].icon[i].nl[0];
			n4 = grid[gridno].icon[i].nl[1];
		    } else if (j == 2) {
			nn[0] = grid[gridno].icon[i].nl[2];
			nn[1] = grid[gridno].icon[i].nl[0];
			or[0] = 2;
			n1 = grid[gridno].icon[i].nl[0];
			n2 = grid[gridno].icon[i].nl[1];
			n3 = grid[gridno].icon[i].nl[1];
			n4 = grid[gridno].icon[i].nl[2];
		    }
		    
		    //findsharednodes(gridno, i, el, nn, or);
		    midx = (grid[gridno].xord[nn[0]] + grid[gridno].xord[nn[1]]) * 0.5;
		    midy = (grid[gridno].yord[nn[0]] + grid[gridno].yord[nn[1]]) * 0.5;
		    qx[0] = (grid[gridno].xord[n1] + grid[gridno].xord[n2]) * 0.5;
		    qy[0] = (grid[gridno].yord[n1] + grid[gridno].yord[n2]) * 0.5;
		    qx[1] = (grid[gridno].xord[n3] + grid[gridno].xord[n4]) * 0.5;
		    qy[1] = (grid[gridno].yord[n3] + grid[gridno].yord[n4]) * 0.5;
		    ns[0] = n1;
		    ns[1] = n2;
		    ns[2] = n3;
		    ns[3] = n4;
//printf("iel %d %d: %d %d %d %d (%d %d %d %d)\n", i+1, bel, n1+1, n2+1, n3+1, n4+1, nn[0] + 1, nn[1] + 1, or[0], or[1]);
		    if (ipos == 0) {
			or[1] = 0;
			n1 = grid[gridno].icon[el].nl[1];
			n2 = grid[gridno].icon[el].nl[2];
			n3 = grid[gridno].icon[el].nl[2];
			n4 = grid[gridno].icon[el].nl[0];
		    } else if (ipos == 1) {
			or[1] = 1;
			n1 = grid[gridno].icon[el].nl[2];
			n2 = grid[gridno].icon[el].nl[0];
			n3 = grid[gridno].icon[el].nl[0];
			n4 = grid[gridno].icon[el].nl[1];
		    } else if (ipos == 2) {
			or[1] = 2;
			n1 = grid[gridno].icon[el].nl[0];
			n2 = grid[gridno].icon[el].nl[1];
			n3 = grid[gridno].icon[el].nl[1];
			n4 = grid[gridno].icon[el].nl[2];
		    }
		    qx[2] = (grid[gridno].xord[n1] + grid[gridno].xord[n2]) * 0.5;
		    qy[2] = (grid[gridno].yord[n1] + grid[gridno].yord[n2]) * 0.5;
		    qx[3] = (grid[gridno].xord[n3] + grid[gridno].xord[n4]) * 0.5;
		    qy[3] = (grid[gridno].yord[n3] + grid[gridno].yord[n4]) * 0.5;
		    qx[4] = qx[0];
		    qy[4] = qy[0];
		    ns[4] = n1;
		    ns[5] = n2;
		    ns[6] = n3;
		    ns[7] = n4;

//printf("el %d %d: %d %d %d %d (%d %d %d %d)\n", el+1, ipos, n1+1, n2+1, n3+1, n4+1, nn[0] + 1, nn[1] + 1, or[0], or[1]);
//
/*
for (k=0;k<4;k++) {
    //char buf[100];
    //sprintf(buf, "%d: %d %d", k, ns[2 * k] + 1, ns[2 * k + 1] + 1);
    printf("(%d: %d %d)", k, ns[2 * k] + 1, ns[2 * k + 1] + 1);
    //writestr(qx[k], qy[k], 0, 0, buf);
}
printf("\n");
*/
		    if (!inregion(qx, qy, 5, midx, midy) ) {
		    //if (i == 32438 || i == 32577) {
		    setcolor(2);
		    my_circle(midx, midy);
		    //if (!inregion(qx, qy, 5, midx, midy) ) {
		        //printf("No good %d %d\n", i + 1, el + 1);
		    //} else {
		        //printf("El OK %d %d\n", i + 1, el + 1);
		    //}
setlinewidth(4);
		    setcolor(3);
		    drawpoly(qx, qy, 5);
setlinewidth(1);
		    }
		}
	    }
	}
    }
    for (j = 0; j < 3; j++) {
	free(adjels[j]);
    }
}

void find_quadrangles(int gridno)
{
    int i, j, k;
    int n1, n2, n3, n4, nn[4];
    int e1, e2, e3, eladj, opel, addn, form, oldnmel;
    double x1, y1, x2, y2, x3, y3, f1, f2, g1, g2, ang;
    extern double qcutoff;
    getassoc_elements(gridno);
/*
   for (i = 0; i < grid[gridno].nmel; i++) {
   grid[gridno].ellist[i] = 0;
   if (grid[gridno].icon[i].type != 3) {
   errwin("Grid must be composed entirely of 3 node triangles");
   free_edgelist(gridno, adjel);
   for (j = 0; j < 3; j++) {
   free(adjels[j]);
   }
   return;

   }
   }
 */
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
    }
    oldnmel = grid[gridno].nmel;
    for (i = 0; i < oldnmel; i++) {
	if (grid[gridno].icon[i].type != 3) {
	    grid[gridno].ellist[i] = 0;
	    continue;
	}
/* element already taken */
	if (grid[gridno].ellist[i]) {
	    continue;
	}
	for (j = 0; j < 3; j++) {
	    if ((eladj = adjels[j][i]) != -1 && grid[gridno].ellist[eladj] != 1) {
/* find the side the i element crosses wrt to eladj */
		for (k = 0; k < 3; k++) {
		    if (i == adjels[k][eladj]) {
			break;
		    }
		}
/*
   printf("at el = %d, adjel = %d, locs = %d %d\n", i, eladj, j, k);
 */
		switch (j) {
		case 0:
		    switch (k) {
		    case 0:
			addn = grid[gridno].icon[eladj].nl[2];
			break;
		    case 1:
			addn = grid[gridno].icon[eladj].nl[0];
			break;
		    case 2:
			addn = grid[gridno].icon[eladj].nl[1];
			break;
		    }
		    n1 = grid[gridno].icon[i].nl[0];
		    n2 = addn;
		    n3 = grid[gridno].icon[i].nl[1];
		    n4 = grid[gridno].icon[i].nl[2];
		    break;
		case 1:
		    switch (k) {
		    case 0:
			addn = grid[gridno].icon[eladj].nl[2];
			break;
		    case 1:
			addn = grid[gridno].icon[eladj].nl[0];
			break;
		    case 2:
			addn = grid[gridno].icon[eladj].nl[1];
			break;
		    }
		    n1 = grid[gridno].icon[i].nl[1];
		    n2 = addn;
		    n3 = grid[gridno].icon[i].nl[2];
		    n4 = grid[gridno].icon[i].nl[0];
		    break;
		case 2:
		    switch (k) {
		    case 0:
			addn = grid[gridno].icon[eladj].nl[2];
			break;
		    case 1:
			addn = grid[gridno].icon[eladj].nl[0];
			break;
		    case 2:
			addn = grid[gridno].icon[eladj].nl[1];
			break;
		    }
		    n1 = grid[gridno].icon[i].nl[2];
		    n2 = addn;
		    n3 = grid[gridno].icon[i].nl[0];
		    n4 = grid[gridno].icon[i].nl[1];
		    break;
		}		/* end switch (j) */
		setcolor(j + 2);
		nn[0] = n1;
		nn[1] = n2;
		nn[2] = n3;
		nn[3] = n4;
		form = 1;
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
		    /*
		       printf("Element %d: %lf %lf %lf\n", i+1,  f1 * f2 + g1 * g2, (hypot(f1, g1) * hypot(f2, g2)), acos((f1 * f2 + g1 * g2) / (hypot(f1, g1) * hypot(f2, g2))));
		     */
		    if (fabs(f1 * f2 + g1 * g2 - (hypot(f1, g1) * hypot(f2, g2))) < 0.0000005) {
			ang = M_PI;
		    } else {
			ang = acos((f1 * f2 + g1 * g2) / (hypot(f1, g1) * hypot(f2, g2)));
		    }
		    if (fabs(ang) > (90.0 + qcutoff) * M_PI / 180.0) {
/*
   printf("Unacceptable element %lf\n", ang * 180.0 / M_PI);
 */
			form = 0;
			break;
		    }
		    if (fabs(ang) < (90.0 - qcutoff) * M_PI / 180.0) {
/*
   printf("Unacceptable element %lf\n", ang * 180.0 / M_PI);
 */
			form = 0;
			break;
		    }
		}		/* end for k */
		if (form) {
		    if (grid[gridno].icon[i].type == 3 && grid[gridno].icon[eladj].type == 3) {
			grid[gridno].ellist[i] = 1;
			grid[gridno].ellist[eladj] = 1;
			for (k = 0; k < 4; k++) {
			    x1 = grid[gridno].xord[nn[k]];
			    x2 = grid[gridno].xord[nn[(k + 1) % 4]];
			    y1 = grid[gridno].yord[nn[k]];
			    y2 = grid[gridno].yord[nn[(k + 1) % 4]];
			    my_move2(x1, y1);
			    my_draw2(x2, y2);
			}	/* end for k */
			add_quad_element(gridno, n1, n2, n3, n4);
		    } else {
/*
   adjacent element was not a triangle
   printf("Major failure in type\n");
 */
		    }
		    goto bustout;
		}		/* end if form */
	    }			/* end if */
	}			/* end for j */
      bustout:;
    }				/* end for i */
    free_edgelist(gridno, adjel);
    for (j = 0; j < 3; j++) {
	free(adjels[j]);
    }
    compact_grid(gridno);
}

void create_triquad(int gridno)
{
    int i, j, k;
    int nd, n1, n2, n3, n4, nn[4];
    int eladj, addn;
    double x1, y1, x2, y2, d1, d2;
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type != 3) {
	    errwin("Grid is not composed of strictly linear triangles, operation cancelled");
	    return;
	}
    }
    getassoc_elements(gridno);
    nd = grid[gridno].nmnp;
    for (i = 0; i < grid[gridno].nmel; i++) {
	if (grid[gridno].icon[i].type == 6) {
	    continue;
	}
	grid[gridno].icon[i].type = 6;
	for (j = 0; j < 3; j++) {
	    if ((eladj = adjels[j][i]) >= 0) {
		for (k = 0; k < 3; k++) {
		    if (i == adjels[k][eladj]) {
			break;
		    }
		}
		adjels[j][i] = -2;
		adjels[k][eladj] = -2;
		n1 = grid[gridno].icon[i].nl[j];
		n2 = grid[gridno].icon[i].nl[(j + 1) % 3];
		setcolor(j + 2);
		x1 = grid[gridno].xord[n1];
		x2 = grid[gridno].xord[n2];
		y1 = grid[gridno].yord[n1];
		y2 = grid[gridno].yord[n2];
		d1 = grid[gridno].depth[n1];
		d2 = grid[gridno].depth[n2];
		my_circlefilled((x1 + x2) * 0.5, (y1 + y2) * 0.5);
		add_node(gridno, (x1 + x2) * 0.5, (y1 + y2) * 0.5, (d1 + d2) * 0.5);
		switch (j) {
		case 0:
		    grid[gridno].icon[i].nl[3] = nd;
		    break;
		case 1:
		    grid[gridno].icon[i].nl[4] = nd;
		    break;
		case 2:
		    grid[gridno].icon[i].nl[5] = nd;
		    break;
		}
		switch (k) {	/* set element entry for adjacent element */
		case 0:
		    grid[gridno].icon[eladj].nl[3] = nd;
		    break;
		case 1:
		    grid[gridno].icon[eladj].nl[4] = nd;
		    break;
		case 2:
		    grid[gridno].icon[eladj].nl[5] = nd;
		    break;
		}
		nd++;
	    }
	    /* end if */
	    else if (eladj == -1) {	/* no adjacent element */
		adjels[j][i] = -2;
		n1 = grid[gridno].icon[i].nl[j];
		n2 = grid[gridno].icon[i].nl[(j + 1) % 3];
		setcolor(1);
		x1 = grid[gridno].xord[n1];
		x2 = grid[gridno].xord[n2];
		y1 = grid[gridno].yord[n1];
		y2 = grid[gridno].yord[n2];
		d1 = grid[gridno].depth[n1];
		d2 = grid[gridno].depth[n2];
		my_circlefilled((x1 + x2) * 0.5, (y1 + y2) * 0.5);
		add_node(gridno, (x1 + x2) * 0.5, (y1 + y2) * 0.5, (d1 + d2) * 0.5);
		switch (j) {
		case 0:
		    grid[gridno].icon[i].nl[3] = nd;
		    break;
		case 1:
		    grid[gridno].icon[i].nl[4] = nd;
		    break;
		case 2:
		    grid[gridno].icon[i].nl[5] = nd;
		    break;
		}
		nd++;
	    }
	}			/* end for j */
    }				/* end for i */
    free_edgelist(gridno, adjel);
    for (j = 0; j < 3; j++) {
	free(adjels[j]);
    }
}

/*
 * get al list of all elements associated with a give node 
 */
void getnode_els(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[3], v, nmnp;
    struct edgelist *e;

    nmnp = grid[gridno].nmnp;
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) malloc(sizeof(int));
	e[i].list[0] = 0;
    }

    for (i = 0; i < grid[gridno].nmel; i++) {
	el[0] = grid[gridno].icon[i].nl[0];
	el[1] = grid[gridno].icon[i].nl[1];
	el[2] = grid[gridno].icon[i].nl[2];
	el[3] = grid[gridno].icon[i].nl[3];

	for (j = 0; j < grid[gridno].icon[i].nn; j++) {
	    n1 = el[j];
	    cnt1 = e[n1].list[0];
	    e[n1].list = (int *) realloc(e[n1].list, (cnt1 + 2) * sizeof(int));
	    e[n1].list[cnt1 + 1] = i;
	    e[n1].list[0] = cnt1 + 1;
	}
    }
    connecl = e;
}

/*
 * return the elements associated with a given node
 */
void returnnode_els(int node, int *nels, int *els)
{
    int i;
    for (i = 0; i < connecl[node].list[0]; i++) {
	els[i] = connecl[node].list[i + 1];
    }
    *nels = connecl[node].list[0];
}

void freenode_els(gridno)
{
    int i;
    for (i = 0; i < grid[gridno].nmnp; i++) {
	free(connecl[i].list);
    }
    free(connecl);
}

void markelements(int gridno)
{
    int i, j, nn, nels, els[200];
    getnode_els(gridno);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	returnnode_els(i, &nels, els);
	if (nels > 0) {
	    for (j = 0; j < nels; j++) {
		if (grid[gridno].icon[els[j]].nn == 4) {
		    grid[gridno].nlist[i] = 1;
		    break;
		}
	    }
	}
    }
    freenode_els(gridno);
}

/*
 * mark nodes within a region
 */
void markregion(int gridno)
{
    int i;
    double x, y;
    for (i = 0; i < grid[gridno].nmnp; i++) {
	x = grid[gridno].xord[i];
	y = grid[gridno].yord[i];
	if (region_flag && !inregion(regionx, regiony, nregion, x, y)) {
	    grid[gridno].nlist[i] = 1;
	}
    }
}

void adjustnode(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[3], v, nmnp;
    double x, y;

    /*
       if (!grid_is_triangles(gridno)) {
       errwin("Grid must be compose entirely of linear elements");
       return;
       }
     */
    getadjnodes(gridno);
    nmnp = grid[gridno].nmnp;
    markboundaries(gridno);
    markelements(gridno);
    for (i = 0; i < nmnp; i++) {
	if (grid[gridno].nlist[i] == 0) {
	    x = y = 0.0;
	    for (j = 0; j < adjn[i].list[0]; j++) {
		x = x + grid[gridno].xord[adjn[i].list[j + 1]];
		y = y + grid[gridno].yord[adjn[i].list[j + 1]];
	    }
	    grid[gridno].xord[i] = x / adjn[i].list[0];
	    grid[gridno].yord[i] = y / adjn[i].list[0];
	}
    }
    free_edgelist(gridno, adjn);
    do_drawgrid();
}

void adjustnode_region(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[3], v, nmnp;
    double x, y;

    getadjnodes(gridno);
    nmnp = grid[gridno].nmnp;
    markboundaries(gridno);
    markelements(gridno);
    markregion(gridno);
    for (i = 0; i < nmnp; i++) {
	if (grid[gridno].nlist[i] == 0) {
	    x = y = 0.0;
	    for (j = 0; j < adjn[i].list[0]; j++) {
		x = x + grid[gridno].xord[adjn[i].list[j + 1]];
		y = y + grid[gridno].yord[adjn[i].list[j + 1]];
	    }
	    grid[gridno].xord[i] = x / adjn[i].list[0];
	    grid[gridno].yord[i] = y / adjn[i].list[0];
	}
    }
    free_edgelist(gridno, adjn);
    do_drawgrid();
}

void getadjnodes(int gridno)
{
    int i, j, f1, f2, cnt1, cnt2, n1, n2, n3, el[4], v, nmnp, nnodes;
    struct edgelist *e;

    nmnp = grid[gridno].nmnp;
    e = (struct edgelist *) malloc(nmnp * sizeof(struct edgelist));
    if (e == NULL) {
	errwin("Can't malloc edgelist");
	return;
    }
    for (i = 0; i < nmnp; i++) {
	e[i].list = (int *) calloc(1, sizeof(int));
	e[i].list[0] = 0;
    }

    for (i = 0; i < grid[gridno].nmel; i++) {
	nnodes = grid[gridno].icon[i].nn;
	el[0] = grid[gridno].icon[i].nl[0];
	el[1] = grid[gridno].icon[i].nl[1];
	el[2] = grid[gridno].icon[i].nl[2];
	if (nnodes == 4) {
	    el[3] = grid[gridno].icon[i].nl[3];
	}

	for (j = 0; j < nnodes; j++) {
	    n1 = el[j];
	    n2 = el[(j + 1) % nnodes];

	    cnt1 = e[n1].list[0];
	    e[n1].list = (int *) realloc(e[n1].list, (cnt1 + 2) * sizeof(int));
	    e[n1].list[cnt1 + 1] = n2;
	    e[n1].list[0] = cnt1 + 1;

	    cnt2 = e[n2].list[0];
	    e[n2].list = (int *) realloc(e[n2].list, (cnt2 + 2) * sizeof(int));
	    e[n2].list[cnt2 + 1] = n1;
	    e[n2].list[0] = cnt2 + 1;
	}
    }
    adjn = e;
}

void dofast(void)
{
    int cnt1, i, n1, n2, nmnp = grid[curgrid].nmnp;
    double x, y, x1, y1;

/*
   getedgelist(0);
 */
    for (i = 0; i < nmnp; i++) {
	x1 = grid[curgrid].xord[i];
	y1 = grid[curgrid].yord[i];
	cnt1 = 0;
	while (edgel[i].list[2 * cnt1] != -1) {
	    my_move2(x1, y1);
	    n2 = edgel[i].list[2 * cnt1];
	    x = grid[curgrid].xord[n2];
	    y = grid[curgrid].yord[n2];
	    my_draw2(x, y);
	    cnt1++;
	}
    }
/*
   free_edgelist(0, edgel);
 */
}

#undef dmin
#undef dmax

void get_gridminmax(int gridno)
{
    int i;

    grid[gridno].xmin = grid[gridno].xord[0];
    grid[gridno].ymin = grid[gridno].yord[0];
    grid[gridno].xmax = grid[gridno].xord[0];
    grid[gridno].ymax = grid[gridno].yord[0];
    grid[gridno].dmin = grid[gridno].depth[0];
    grid[gridno].dmax = grid[gridno].depth[0];
    for (i = 1; i < grid[gridno].nmnp; i++) {
	if (grid[gridno].xmin > grid[gridno].xord[i])
	    grid[gridno].xmin = grid[gridno].xord[i];
	if (grid[gridno].ymin > grid[gridno].yord[i])
	    grid[gridno].ymin = grid[gridno].yord[i];
	if (grid[gridno].xmax < grid[gridno].xord[i])
	    grid[gridno].xmax = grid[gridno].xord[i];
	if (grid[gridno].ymax < grid[gridno].yord[i])
	    grid[gridno].ymax = grid[gridno].yord[i];
	if (grid[gridno].dmin > grid[gridno].depth[i])
	    grid[gridno].dmin = grid[gridno].depth[i];
	if (grid[gridno].dmax < grid[gridno].depth[i])
	    grid[gridno].dmax = grid[gridno].depth[i];
    }
}

int gengrid2(int gridno, double x1, double y1, double x2, double y2)
{
    char buf[256];
    int i, j, itmp;
    double xmin, xmax, ymin, ymax, dmin, dmax;
    Free_grid(gridno);
    strcpy(grid[gridno].alphid, "Generated grid");
    grid[gridno].nmel = 2;
    grid[gridno].nmnp = 4;
    if (!Allocate_grid(gridno, grid[gridno].nmel, grid[gridno].nmnp)) {
	errwin("Can't allocate memory for grid");
	Free_grid(gridno);
	return 0;
    }
    if (x1 > x2) {
	fswap(&x1, &x2);
    }
    if (y1 > y2) {
	fswap(&y1, &y2);
    }
    xmin = grid[gridno].xmin = x1;
    xmax = grid[gridno].xmax = x2;
    ymin = grid[gridno].ymin = y1;
    ymax = grid[gridno].ymax = y2;
    dmin = grid[gridno].dmin = 1.0;
    dmax = grid[gridno].dmax = 1.0;
    grid[gridno].xord[0] = xmin;
    grid[gridno].xord[1] = xmax;
    grid[gridno].xord[2] = xmax;
    grid[gridno].xord[3] = xmin;
    grid[gridno].yord[0] = ymin;
    grid[gridno].yord[1] = ymin;
    grid[gridno].yord[2] = ymax;
    grid[gridno].yord[3] = ymax;
    grid[gridno].depth[0] = dmin;
    grid[gridno].depth[1] = dmin;
    grid[gridno].depth[2] = dmin;
    grid[gridno].depth[3] = dmin;
    grid[gridno].icon[0].type = 3;
    grid[gridno].icon[0].nn = 3;
    grid[gridno].icon[0].ngeom = 3;
    grid[gridno].icon[0].nl[0] = 0;
    grid[gridno].icon[0].nl[1] = 1;
    grid[gridno].icon[0].nl[2] = 3;
    grid[gridno].icon[1].type = 3;
    grid[gridno].icon[1].nn = 3;
    grid[gridno].icon[1].ngeom = 3;
    grid[gridno].icon[1].nl[0] = 1;
    grid[gridno].icon[1].nl[1] = 2;
    grid[gridno].icon[1].nl[2] = 3;
    return 1;
}

void compute_grad(int gridno, double *gx, double *gy)
{
    int i, j, n[3];
    double a, b, c, d;
    for (i = 0; i < grid[gridno].nmel; i++) {
	for (j = 0; j < 3; j++) {
	    n[j] = grid[gridno].icon[i].nl[j];
	}
	a = b = c = d = 0.0;
	for (j = 0; j < 3; j++) {
	    a += grid[gridno].yord[n[j]] * (grid[gridno].depth[n[(j + 1) % 3]]
					    - grid[gridno].depth[n[(j + 2) % 3]]);
	    b += grid[gridno].depth[n[j]] * (grid[gridno].xord[n[(j + 1) % 3]] - grid[gridno].xord[n[(j + 2) % 3]]);
	    c += grid[gridno].xord[n[j]] * (grid[gridno].yord[n[(j + 1) % 3]] - grid[gridno].yord[n[(j + 2) % 3]]);
	}
	d = -a * grid[gridno].xord[n[0]] - b * grid[gridno].yord[n[0]] - c * grid[gridno].depth[n[0]];
	if (c != 0.0) {
	    gx[i] = a / c;
	    gy[i] = b / c;
	} else {
	    printf("c = 0.0\n");
	}
    }
}

int convert_to_triangles(int gridno)
{
    int i, j, k, ok = 0, ne;
    int nd, n1, n2, n3, n4, nn[4];
    int eladj, addn;
    double x1, y1, x2, y2, d1, d2;
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
	if (grid[gridno].icon[i].type != 3) {
	    ok = 1;
	}
    }
    if (!ok) {
	errwin("Grid is already composed of linear triangles");
	return 0;
    }
    nd = grid[gridno].nmnp;
    ne = grid[gridno].nmel;
    for (i = 0; i < ne; i++) {
	if (grid[gridno].icon[i].type == 6) {
	    grid[gridno].icon[i].type = 3;
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	} else if (grid[gridno].icon[i].type == 4) {
	    grid[gridno].ellist[i] = 1;
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    n4 = grid[gridno].icon[i].nl[3];
	    add_element(gridno, n1, n2, n3);
	    add_element(gridno, n1, n3, n4);
	} else if (grid[gridno].icon[i].type == 8) {
	    grid[gridno].ellist[i] = 1;
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    n4 = grid[gridno].icon[i].nl[3];
	    add_element(gridno, n1, n2, n3);
	    add_element(gridno, n1, n3, n4);
	}
    }
    compact_grid(gridno);
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
    }
    return 1;
}

void drop_middle_nodes(int gridno)
{
    int i, j, k, ok = 0, ne;
    int nd, n1, n2, n3, n4, nn[4];
    int eladj, addn;
    double x1, y1, x2, y2, d1, d2;
    for (i = 0; i < grid[gridno].nmel; i++) {
	grid[gridno].ellist[i] = 0;
	if (grid[gridno].icon[i].type != 3) {
	    ok = 1;
	}
    }
    if (!ok) {
	errwin("Grid is already composed of linear triangles");
	return;
    }
    nd = grid[gridno].nmnp;
    ne = grid[gridno].nmel;
    for (i = 0; i < ne; i++) {
	if (grid[gridno].icon[i].type == 1) {
	    grid[gridno].icon[i].type = 2;
	    grid[gridno].icon[i].nn = 2;
	    grid[gridno].icon[i].ngeom = 2;
	    continue;
	}
	if (grid[gridno].icon[i].type == 6) {
	    grid[gridno].icon[i].type = 3;
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	    continue;
	}
	if (grid[gridno].icon[i].type == 8) {
	    grid[gridno].icon[i].type = 4;
	    grid[gridno].icon[i].nn = 4;
	    grid[gridno].icon[i].ngeom = 4;
	    continue;
	}
	if (grid[gridno].icon[i].type == 5) {
	    grid[gridno].icon[i].type = 3;
	    grid[gridno].icon[i].nn = 3;
	    grid[gridno].icon[i].ngeom = 3;
	    grid[gridno].icon[i].nl[1] = grid[gridno].icon[i].nl[2];
	    grid[gridno].icon[i].nl[2] = grid[gridno].icon[i].nl[3];
	    continue;
	}
    }
    compact_grid(gridno);
}

typedef struct edge_struct {
    int n1, n2;
    int els[2];
    int pos[2];
} Edge;

static int nedges;
static Edge *edges;

int edgeinlist(int gridno, int n1, int n2)
{
    int i;
    for (i = 0; i < nedges; i++) {
	if ((edges[i].n1 == n1 && edges[i].n2 == n2)
	    || (edges[i].n2 == n1 && edges[i].n1 == n2)) {
	    return i;
	}
    }
    return -1;
}

int addedge(int gridno, int n1, int n2, int e, int pos)
{
    int edge = edgeinlist(gridno, n1, n2);
    if (edge == -1) {		/* not in list, so add it */
	if (nedges) {
	    edges = (Edge *) realloc(edges, (nedges + 1) * sizeof(Edge));
	} else {
	    edges = (Edge *) calloc(1, sizeof(Edge));
	}
	edges[nedges].n1 = n1;
	edges[nedges].n2 = n2;
	edges[nedges].els[0] = e;
	edges[nedges].els[1] = -1;
	edges[nedges].pos[0] = pos;
	edges[nedges].pos[1] = -1;
	nedges++;
	return 1;
    } else {			/* already in list */
	edges[edge].els[1] = e;
	edges[edge].pos[1] = pos;
	return 0;
    }
}

int getedge(int gridno, int n1, int n2)
{
    return 0;
}

int addmiddle(int gridno, int n1, int n2)
{
    return 0;
}

int getmiddle(int gridno, int n1, int n2)
{
    return 0;
}

void create_edgelist(int gridno)
{
    int i, j, n1, n2;
    nedges = 0;
    for (i = 0; i < grid[gridno].nmel; i++) {
	switch (grid[gridno].icon[i].type) {
	case 0:
	    break;
	case 1:
	    break;
	case 2:
	    break;
	case 3:
	    for (j = 0; j < 3; j++) {
		n1 = grid[gridno].icon[i].nl[j];
		n2 = grid[gridno].icon[i].nl[(j + 1) % 3];
		addedge(gridno, n1, n2, i, j);
	    }
	    break;
	case 4:
	    for (j = 0; j < 4; j++) {
		n1 = grid[gridno].icon[i].nl[j];
		n2 = grid[gridno].icon[i].nl[(j + 1) % 4];
		addedge(gridno, n1, n2, i, j);
	    }
	    break;
	case 5:
	    break;
	case 6:
	    break;
	case 8:
	    break;
	}
    }
}

void delete_edgelist(void)
{
    if (nedges) {
	free(edges);
    }
    edges = (Edge *) NULL;
    nedges = 0;
}

void adj_middle_nodes(int gridno)
{
    int i, j, k, ok = 0, ne, e, savej;
    int nd, n1, n2, n3, n4, n5, nn[4], pos;
    int eladj, addn, second, el;
    double x1, y1, x2, y2, d1, d2;
    ne = grid[gridno].nmel;
    for (i = 0; i < ne; i++) {
	el = i;
	switch (grid[gridno].icon[el].type) {
	case 0:
	    break;
	case 1:
	    n1 = grid[gridno].icon[el].nl[0];
	    n2 = grid[gridno].icon[el].nl[1];
	    n3 = grid[gridno].icon[el].nl[2];
	    x1 = grid[gridno].xord[n1];
	    x2 = grid[gridno].xord[n2];
	    y1 = grid[gridno].yord[n1];
	    y2 = grid[gridno].yord[n2];
	    d1 = grid[gridno].depth[n1];
	    d2 = grid[gridno].depth[n2];
	    grid[gridno].xord[n3] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n3] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n3] = (d1 + d2) * 0.5;
	    break;
	case 2:
	    break;
	case 3:
	    break;
	case 4:
	    break;
	case 5:
	    n1 = grid[gridno].icon[el].nl[0];
	    n2 = grid[gridno].icon[el].nl[2];
	    n3 = grid[gridno].icon[el].nl[1];	/* middle node */
	    n4 = grid[gridno].icon[el].nl[3];
	    n5 = grid[gridno].icon[el].nl[4];	/* middle node */
	    x1 = grid[gridno].xord[n1];
	    x2 = grid[gridno].xord[n2];
	    y1 = grid[gridno].yord[n1];
	    y2 = grid[gridno].yord[n2];
	    d1 = grid[gridno].depth[n1];
	    d2 = grid[gridno].depth[n2];
	    grid[gridno].xord[n3] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n3] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n3] = (d1 + d2) * 0.5;

	    x1 = grid[gridno].xord[n3];
	    x2 = grid[gridno].xord[n4];
	    y1 = grid[gridno].yord[n3];
	    y2 = grid[gridno].yord[n4];
	    d1 = grid[gridno].depth[n3];
	    d2 = grid[gridno].depth[n4];
	    grid[gridno].xord[n5] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n5] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n5] = (d1 + d2) * 0.5;
	    break;
	case 6:
	    for (j = 0; j < 3; j++) {
		n1 = grid[gridno].icon[el].nl[j];
		n2 = grid[gridno].icon[el].nl[(j + 1) % 3];
		n3 = grid[gridno].icon[el].nl[j + 3];
		x1 = grid[gridno].xord[n1];
		x2 = grid[gridno].xord[n2];
		y1 = grid[gridno].yord[n1];
		y2 = grid[gridno].yord[n2];
		d1 = grid[gridno].depth[n1];
		d2 = grid[gridno].depth[n2];
		grid[gridno].xord[n3] = (x1 + x2) * 0.5;
		grid[gridno].yord[n3] = (y1 + y2) * 0.5;
		grid[gridno].depth[n3] = (d1 + d2) * 0.5;
	    }
	    break;
	case 8:
	    for (j = 0; j < 4; j++) {
		n1 = grid[gridno].icon[el].nl[j];
		n2 = grid[gridno].icon[el].nl[(j + 1) % 4];
		n3 = grid[gridno].icon[el].nl[j + 4];
		x1 = grid[gridno].xord[n1];
		x2 = grid[gridno].xord[n2];
		y1 = grid[gridno].yord[n1];
		y2 = grid[gridno].yord[n2];
		d1 = grid[gridno].depth[n1];
		d2 = grid[gridno].depth[n2];
		grid[gridno].xord[n3] = (x1 + x2) * 0.5;
		grid[gridno].yord[n3] = (y1 + y2) * 0.5;
		grid[gridno].depth[n3] = (d1 + d2) * 0.5;
	    }
	    break;
	}
    }
}

/*
 * add middle nodes to a grid
 */
void add_middle_nodes(int gridno)
{
    int i, j, k, ok = 0, ne, e, savej;
    int nd, n1, n2, n3, n4, n5, nn[4], pos;
    int eladj, addn, second, el;
    double x1, y1, x2, y2, d1, d2;
    nd = grid[gridno].nmnp;
    ne = grid[gridno].nmel;
    create_edgelist(gridno);
    for (i = 0; i < nedges; i++) {
	n1 = edges[i].n1;
	n2 = edges[i].n2;
	for (k = 0; k < 2; k++) {
	    if ((el = edges[i].els[k]) >= 0) {
		switch (grid[gridno].icon[el].type) {
		case 0:
		    break;
		case 1:
		    break;
		case 2:
		    break;
		case 3:
		    pos = edges[i].pos[k];
		    if (k == 0) {
			x1 = grid[gridno].xord[n1];
			x2 = grid[gridno].xord[n2];
			y1 = grid[gridno].yord[n1];
			y2 = grid[gridno].yord[n2];
			d1 = grid[gridno].depth[n1];
			d2 = grid[gridno].depth[n2];
			e = add_node(gridno, (x1 + x2) * 0.5, (y1 + y2) * 0.5, (d1 + d2) * 0.5);
		    }
		    grid[gridno].icon[el].nl[pos + 3] = e;
		    break;
		case 4:
		    pos = edges[i].pos[k];
		    if (k == 0) {
			x1 = grid[gridno].xord[n1];
			x2 = grid[gridno].xord[n2];
			y1 = grid[gridno].yord[n1];
			y2 = grid[gridno].yord[n2];
			d1 = grid[gridno].depth[n1];
			d2 = grid[gridno].depth[n2];
			e = add_node(gridno, (x1 + x2) * 0.5, (y1 + y2) * 0.5, (d1 + d2) * 0.5);
		    }
		    if (e == 0) {
			printf("Got 1\n");
		    }
		    grid[gridno].icon[el].nl[pos + 4] = e;
		    break;
		case 5:
		    break;
		case 6:
		    for (j = 0; j < 3; j++) {
			n1 = grid[gridno].icon[el].nl[j];
			n2 = grid[gridno].icon[el].nl[(j + 1) % 3];
			n3 = grid[gridno].icon[el].nl[j + 3];
			x1 = grid[gridno].xord[n1];
			x2 = grid[gridno].xord[n2];
			y1 = grid[gridno].yord[n1];
			y2 = grid[gridno].yord[n2];
			d1 = grid[gridno].depth[n1];
			d2 = grid[gridno].depth[n2];
			grid[gridno].xord[n3] = (x1 + x2) * 0.5;
			grid[gridno].yord[n3] = (y1 + y2) * 0.5;
			grid[gridno].depth[n3] = (d1 + d2) * 0.5;
		    }
		    break;
		case 8:
		    for (j = 0; j < 4; j++) {
			n1 = grid[gridno].icon[el].nl[j];
			n2 = grid[gridno].icon[el].nl[(j + 1) % 4];
			n3 = grid[gridno].icon[el].nl[j + 4];
			x1 = grid[gridno].xord[n1];
			x2 = grid[gridno].xord[n2];
			y1 = grid[gridno].yord[n1];
			y2 = grid[gridno].yord[n2];
			d1 = grid[gridno].depth[n1];
			d2 = grid[gridno].depth[n2];
			grid[gridno].xord[n3] = (x1 + x2) * 0.5;
			grid[gridno].yord[n3] = (y1 + y2) * 0.5;
			grid[gridno].depth[n3] = (d1 + d2) * 0.5;
		    }
		    break;
		}
	    }
	}
    }
    for (i = 0; i < grid[gridno].nmel; i++) {
	switch (grid[gridno].icon[i].type) {
	case 0:
	    break;
	case 1:
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    n3 = grid[gridno].icon[i].nl[2];
	    x1 = grid[gridno].xord[n1];
	    x2 = grid[gridno].xord[n2];
	    y1 = grid[gridno].yord[n1];
	    y2 = grid[gridno].yord[n2];
	    d1 = grid[gridno].depth[n1];
	    d2 = grid[gridno].depth[n2];
	    grid[gridno].xord[n3] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n3] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n3] = (d1 + d2) * 0.5;
	    break;
	case 2:
	    grid[gridno].icon[i].type = 1;
	    grid[gridno].icon[i].ngeom = 2;
	    grid[gridno].icon[i].nn = 3;
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[1];
	    x1 = grid[gridno].xord[n1];
	    x2 = grid[gridno].xord[n2];
	    y1 = grid[gridno].yord[n1];
	    y2 = grid[gridno].yord[n2];
	    d1 = grid[gridno].depth[n1];
	    d2 = grid[gridno].depth[n2];
	    e = add_node(gridno, (x1 + x2) * 0.5, (y1 + y2) * 0.5, (d1 + d2) * 0.5);
	    grid[gridno].icon[i].nl[2] = e;
	    break;
	case 3:
	    grid[gridno].icon[i].type = 6;
	    grid[gridno].icon[i].ngeom = 3;
	    grid[gridno].icon[i].nn = 6;
	    break;
	case 4:
	    grid[gridno].icon[i].type = 8;
	    grid[gridno].icon[i].ngeom = 4;
	    grid[gridno].icon[i].nn = 8;
	    break;
	case 5:
	    n1 = grid[gridno].icon[i].nl[0];
	    n2 = grid[gridno].icon[i].nl[2];
	    n3 = grid[gridno].icon[i].nl[1];	/* middle node */
	    n4 = grid[gridno].icon[i].nl[3];
	    n5 = grid[gridno].icon[i].nl[4];	/* middle node */
	    x1 = grid[gridno].xord[n1];
	    x2 = grid[gridno].xord[n2];
	    y1 = grid[gridno].yord[n1];
	    y2 = grid[gridno].yord[n2];
	    d1 = grid[gridno].depth[n1];
	    d2 = grid[gridno].depth[n2];
	    grid[gridno].xord[n3] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n3] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n3] = (d1 + d2) * 0.5;

	    x1 = grid[gridno].xord[n3];
	    x2 = grid[gridno].xord[n4];
	    y1 = grid[gridno].yord[n3];
	    y2 = grid[gridno].yord[n4];
	    d1 = grid[gridno].depth[n3];
	    d2 = grid[gridno].depth[n4];
	    grid[gridno].xord[n5] = (x1 + x2) * 0.5;
	    grid[gridno].yord[n5] = (y1 + y2) * 0.5;
	    grid[gridno].depth[n5] = (d1 + d2) * 0.5;
	    break;
	case 6:
	    break;
	case 8:
	    break;
	}
    }
    delete_edgelist();
}

/*
 * get the bounding box for an element
 */
void get_element_bb(int gridno, int el, double *x1, double *y1, double *x2, double *y2)
{
    int i, n;
    for (i = 0; i < grid[gridno].icon[el].nn; i++) {
	n = grid[gridno].icon[el].nl[i];
	if (i) {
	    *x1 = grid[gridno].xord[n] < *x1 ? grid[gridno].xord[n] : *x1;
	    *x2 = grid[gridno].xord[n] > *x2 ? grid[gridno].xord[n] : *x2;
	    *y1 = grid[gridno].yord[n] < *y1 ? grid[gridno].yord[n] : *y1;
	    *y2 = grid[gridno].yord[n] > *y2 ? grid[gridno].yord[n] : *y2;
	} else {
	    *x1 = *x2 = grid[gridno].xord[n];
	    *y1 = *y2 = grid[gridno].yord[n];
	}
    }
}

/*
 * compute the volume of a prism
 */
double element_volume(int gridno, int elem)
{
    double area(int gridno, int elem);
    double d, v = 0.0;
    double ar = area(gridno, elem);
    int n1 = grid[gridno].icon[elem].nl[0];
    int n2 = grid[gridno].icon[elem].nl[1];
    int n3 = grid[gridno].icon[elem].nl[2];
    v = ar * (grid[gridno].depth[n1] + grid[gridno].depth[n2] + grid[gridno].depth[n3]) * 0.333333333333;
    return v;
}

/*
 * compute the volume of a prism given the coordinates of the triangle on top.
 */
double element_volume_vertices(double x1, double y1, double d1, double x2, double y2, double d2, double x3, double y3, double d3)
{
    double d, v = 0.0;
    double ar = area_from_vertices(x1, y1, x2, y2, x3, y3);
    v = ar * (d1 + d2 + d3) * 0.3333333333333333;
    return v;
}

/*
 * compute the circumcenter of a triangle, used for median axis computation
 */
void compute_circumcenter(double xk, double yk, double xl, double yl, double xm, double ym, double *xc, double *yc)
{
    double xlk, ylk, xmk, ymk, xcc, ycc, det, detinv, rlksq, rmksq;
    xlk = xl - xk;
    ylk = yl - yk;
    xmk = xm - xk;
    ymk = ym - yk;
    det = xlk * ymk - xmk * ylk;
    if (fabs(det) <= 1.0e-14) {	/* determinant is zero */
	return;
    }
    detinv = 0.5 / det;
    rlksq = xlk * xlk + ylk * ylk;
    rmksq = xmk * xmk + ymk * ymk;
    xcc = detinv * (rlksq * ymk - rmksq * ylk);
    ycc = detinv * (xlk * rmksq - xmk * rlksq);
    *xc = xcc + xk;
    *yc = ycc + yk;
}

void PrintGrid(Grid * g)
{
    int i;
    printf("triang\n");
    printf("%d %d\n", g->nmel, g->nmnp);
    for (i = 0; i < g->nmnp; i++) {
	printf("%d %lf %lf %lf\n", i + 1, g->xord[i], g->yord[i], g->depth[i]);
    }
    for (i = 0; i < g->nmel; i++) {
	printf("%d 3 %d %d %d\n", i + 1, g->icon[i].nl[0] + 1, g->icon[i].nl[1] + 1, g->icon[i].nl[2] + 1);
    }
}

int IsCCW(Grid * g, int elno)
{
    double a[2], b[2], ar;
    int n1, n2, n3;
    n1 = g->icon[elno].nl[0];
    n2 = g->icon[elno].nl[1];
    n3 = g->icon[elno].nl[2];
    a[0] = g->xord[n3] - g->xord[n2];
    a[1] = g->xord[n1] - g->xord[n3];
    b[0] = g->yord[n2] - g->yord[n3];
    b[1] = g->yord[n3] - g->yord[n1];
    ar = 0.5 * (a[1] * b[0] - a[0] * b[1]);
    if (ar <= 0) {
	return 0;
    }
    return 1;
}

double ElementArea(Grid * g, int elno)
{
    double a[2], b[2];
    int n1, n2, n3;
    n1 = g->icon[elno].nl[0];
    n2 = g->icon[elno].nl[1];
    n3 = g->icon[elno].nl[2];
    a[0] = g->xord[n3] - g->xord[n2];
    a[1] = g->xord[n1] - g->xord[n3];
    b[0] = g->yord[n2] - g->yord[n3];
    b[1] = g->yord[n3] - g->yord[n1];
    return 0.5 * (a[1] * b[0] - a[0] * b[1]);
}

double ElementVolume(Grid * g, int elno)
{
    int n1, n2, n3;
    double ar = ElementArea(g, elno);
    n1 = g->icon[elno].nl[0];
    n2 = g->icon[elno].nl[1];
    n3 = g->icon[elno].nl[2];
    return ar * (g->depth[n3] + g->depth[n2] + g->depth[n3]) * 0.3333333333333333;
}

void GetCenter(Grid * g, int elno, double *cx, double *cy)
{
    double x1, y1, x2, y2, x3, y3;
    int n1, n2, n3, n4;
    n1 = g->icon[elno].nl[0];
    x1 = g->xord[n1];
    y1 = g->yord[n1];
    n2 = g->icon[elno].nl[1];
    x2 = g->xord[n2];
    y2 = g->yord[n2];
    n3 = g->icon[elno].nl[2];
    x3 = g->xord[n3];
    y3 = g->yord[n3];
    *cx = 0.3333333333333333 * (x1 + x2 + x3);
    *cy = 0.3333333333333333 * (y1 + y2 + y3);
}

void FixAreas(Grid * g)
{
    int i;
    for (i = 0; i < g->nmel; i++) {
	if (!IsCCW(g, i)) {
	    int itmp = g->icon[i].nl[0];
	    g->icon[i].nl[0] = g->icon[i].nl[2];
	    g->icon[i].nl[2] = itmp;
	}
    }
}

void FixArea(Grid * g, int elno)
{
    if (!IsCCW(g, elno)) {
	int itmp = g->icon[elno].nl[0];
	g->icon[elno].nl[0] = g->icon[elno].nl[2];
	g->icon[elno].nl[2] = itmp;
    }
}
