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
 * utilities for managing build points
 *
 */

#ifndef lint
static char RCSid[] = "$Id: buildutils.c,v 1.4 2011/09/14 17:44:21 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "defines.h"
#include "globals.h"

double get_depth_element();

/*
 * function declarations
 */
void drawbuild(int buildno);
void drawbuildline(int buildno);
void accumulate_nearest_buildpts(int bno, double wx, double wy);
void delete_build_points(int bno, int *ibuf, int nibuf);
void build_minmax(int bno, double *xmin, double *xmax, double *ymin, double *ymax, double *dmin, double *dmax);
void load_grid(int gridno, int bno);
void find_nearest_buildpt(int bno, double x, double y, int *ind);
void load_from_grid(int gridno, int bno);
void load_centers(int gridno, int bno);
void merge_grid_to_build(int gridno, int bno);
void move_buildpt(int bno, int ind, double wx, double wy);
void clean_up_grid(int gridno);
void merge_boundary_to_grid(int gridno, int buildno, int boundno);
void orientgrid(int gridno);
void spread_rectangular(int bno, double wx1, double wy1, double wx2, double wy2);
void spread_rectangular_offset(int bno, double wx1, double wy1, double wx2, double wy2);
void spread_random(int bno, double wx1, double wy1, double wx2, double wy2);
void spread_rotated_rect(int bno, double wx1, double wy1, double wx2, double wy2, double wx3, double wy3);
void spread_rotated_rect_offset(int bno, double wx1, double wy1, double wx2, double wy2, double wx3, double wy3);
int compare_build(const void  *iptr, const void *jptr);
void elim_build_dupes(int buildno, double mindist);
void gen_circ(int bno, double xc, double yc, double r, int nr, int na, int nainc);
void do_select_build_region(void);
void delete_build_inside_region(void);
void delete_build_outside_region(void);
int check_crit(int crit, double a, double avgd, double tol);
void auto_spread(int gridno, double x1, double y1, double x2, double y2, double dx, double dy, int crit, int maxloop, double tol, int dgrid, int reg);

void DrawBuildDepths(Buildpts * b, Isolparms ip);

/*
 * Draw build points 
 */
void drawbuild(int buildno)
{
    int i;
    char buf[24];
    if (build[buildno].isol) {
	DrawBuildDepths(&build[buildno], grid[0].ip);
    } else {
	setcolor(build[buildno].color);
	for (i = 0; i < build[buildno].nbuild; i++) {
	    switch (build[buildno].sym) {
	    case 0:
		break;
	    case 1:
		drawdot(build[buildno].bx[i], build[buildno].by[i]);
		break;
	    case 2:
		box(build[buildno].bx[i], build[buildno].by[i]);
		break;
	    case 3:
		my_circle(build[buildno].bx[i], build[buildno].by[i]);
		break;
	    case 4:
		my_circle(build[buildno].bx[i], build[buildno].by[i]);
                sprintf(buf, " %d", i + 1);
                writestr(build[buildno].bx[i], build[buildno].by[i], 0, 0, buf);
		break;
	    }
	}
	flush_pending();
    }
}

void drawbuildline(int buildno)
{
    int i;
    my_move2(build[buildno].bx[0], build[buildno].by[0]);
    for (i = 1; i < build[buildno].nbuild; i++) {
	my_draw2(build[buildno].bx[i], build[buildno].by[i]);
    }
}

void accumulate_nearest_buildpts(int bno, double wx, double wy)
{
    int ind;
    double x, y;

    find_nearest_buildpt(bno, wx, wy, &ind);
    if (ind >= 0) {
	solidbox(build[bno].bx[ind], build[bno].by[ind]);
	ibuf[nibuf++] = ind;
    }
}

void delete_build_points(int bno, int *ibuf, int nibuf)
{
    int i, j, npts = 0;

    int *itmp;

    itmp = (int *) calloc(build[bno].nbuild, sizeof(int));
    if (itmp == NULL || nibuf == 0) {
	return;
    }
    for (j = 0; j < build[bno].nbuild; j++) {
	itmp[j] = 0;
    }
    for (i = 0; i < nibuf; i++) {
	itmp[ibuf[i]] = 1;
    }
    setcolor(0);
    for (j = 0; j < build[bno].nbuild; j++) {
	if (!itmp[j]) {
	    build[bno].bx[npts] = build[bno].bx[j];
	    build[bno].by[npts] = build[bno].by[j];
	    build[bno].db[npts] = build[bno].db[j];
	    npts++;
	} else {
	    if (inwin) {
		box(build[bno].bx[j], build[bno].by[j]);
	    }
	}
    }
    setcolor(1);
    free(itmp);
    if (npts != 0) {
	Reallocate_build(bno, npts);
    } else {
	Free_build(bno);
    }
}

void build_minmax(int bno, double *xmin, double *xmax, double *ymin, double *ymax, double *dmin, double *dmax)
{
    int i;
    if (build[bno].nbuild > 3) {
	*xmin = build[bno].bx[0];
	*ymin = build[bno].by[0];
	*xmax = build[bno].bx[0];
	*ymax = build[bno].by[0];
	*dmin = build[bno].db[0];
	*dmax = build[bno].db[0];
	for (i = 1; i < build[bno].nbuild; i++) {
	    if (*xmin > build[bno].bx[i])
		*xmin = build[bno].bx[i];
	    if (*ymin > build[bno].by[i])
		*ymin = build[bno].by[i];
	    if (*xmax < build[bno].bx[i])
		*xmax = build[bno].bx[i];
	    if (*ymax < build[bno].by[i])
		*ymax = build[bno].by[i];
	    if (*dmin > build[bno].db[i])
		*dmin = build[bno].db[i];
	    if (*dmax < build[bno].db[i])
		*dmax = build[bno].db[i];
	}
    }
}

/*
 * After triangulating build points load them to the grid, replacing
 * the grid's nodes and elements.
 */
void load_grid(int gridno, int bno)
{
    int i, j;

    Free_grid_only(gridno);
    Allocate_grid_only(gridno, nels, build[bno].nbuild);
    for (i = 0; i < build[bno].nbuild; i++) {
	grid[gridno].xord[i] = build[bno].bx[i];
	grid[gridno].yord[i] = build[bno].by[i];
	grid[gridno].depth[i] = build[bno].db[i];
    }
    for (i = 0; i < nels; i++) {
	grid[gridno].icon[i].type = 3;
	grid[gridno].icon[i].nn = 3;
	grid[gridno].icon[i].ngeom = 3;
	for (j = 0; j < 3; j++) {
	    grid[gridno].icon[i].nl[j] = tmptable[i].nl[j];
	}
    }
    free(tmptable);
    CleanGrid(&grid[gridno]);
    DeleteElements(&grid[gridno]);
    setlimits_grid(gridno);
}

void find_nearest_buildpt(int bno, double x, double y, int *ind)
{
    int i;
    double radius = 1e307;
    double tmp;

    *ind = -1;
    for (i = 0; i < build[bno].nbuild; i++) {
	tmp = hypot((x - build[bno].bx[i]), (y - build[bno].by[i]));
	if (tmp < radius) {
	    radius = tmp;
	    *ind = i;
	}
    }
}

/*
 * Load build points from grid gridno
 */
void load_from_grid(int gridno, int bno)
{
    int i, j;

    Free_build(bno);
    Allocate_build(bno, grid[gridno].nmnp);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	build[bno].bx[i] = grid[gridno].xord[i];
	build[bno].by[i] = grid[gridno].yord[i];
	build[bno].db[i] = grid[gridno].depth[i];
    }
}

void load_centers(int gridno, int bno)
{
    int i, j;
    double xg, yg;
    Free_build(bno);
    Allocate_build(bno, grid[gridno].nmel);
    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &xg, &yg);
	build[bno].bx[i] = xg;
	build[bno].by[i] = yg;
	build[bno].db[i] = get_depth_element(gridno, i, xg, yg);
    }
}

/*
 * Merge the grid's nodal points with the current set of build points.
 */
void merge_grid_to_build(int gridno, int bno)
{
    int i;

    init_qadd_build(0);
    for (i = 0; i < grid[gridno].nmnp; i++) {
	qadd_build(bno, grid[gridno].xord[i], grid[gridno].yord[i], grid[gridno].depth[i]);
    }
    flush_qadd_build();
}

void move_buildpt(int bno, int ind, double wx, double wy)
{
    build[bno].bx[ind] = wx;
    build[bno].by[ind] = wy;
    return;
}

void clean_up_grid(int gridno)
{
    int i;
    double x, y;

    for (i = 0; i < grid[gridno].nmel; i++) {
	get_center(gridno, i, &x, &y);
	if (!goodpoint(gridno, x, y)) {
	    grid[gridno].ellist[i] = 1;
	} else {
	    grid[gridno].ellist[i] = 0;
	}
    }
    compact_grid(gridno);
}

void merge_boundary_to_grid(int gridno, int buildno, int boundno)
{
    int i, j;
}

void orientgrid(int gridno)
{
    int i;

    for (i = 0; i < grid[gridno].nmel; i++) {
	if (!comparea(gridno, i)) {
	    iswap(&grid[gridno].icon[i].nl[0], &grid[gridno].icon[i].nl[2]);
	}
    }
}

/*
 * spread build points in a semi-automatic way
 */
extern int npts_to_place;
extern int spreadflag;
extern double spread_dx, spread_dy;

void spread_rectangular(int bno, double wx1, double wy1, double wx2, double wy2)
{
    int i, j;
    double x, y, d = 1.0;
    double dx = fabs(wx1 - wx2);
    double dy = fabs(wy1 - wy2);
    int idx, idy;

    if (spread_dx == 0.0) {
	return;
    }
    if (spread_dy == 0.0) {
	return;
    }
    idx = (int) (dx / spread_dx) + 1;
    idy = (int) (dy / spread_dy) + 1;
    if (idx * idy > 3000) {
	if (!yesno("More than 3000 build points, OK?", "Press Yes or No", "Yes", "No")) {
	    return;
	}
    }
    init_qadd_build(0);
    for (i = 0; i < idx; i++) {
	for (j = 0; j < idy; j++) {
	    qadd_build(bno, wx1 + spread_dx * i, wy1 + spread_dy * j, 1.0);
	    box(wx1 + spread_dx * i, wy1 + spread_dy * j);
	}
    }
    flush_qadd_build(0);
}

void spread_rectangular_offset(int bno, double wx1, double wy1, double wx2, double wy2)
{
    int i, j;
    double x, y, d = 1.0;
    double dx = fabs(wx1 - wx2);
    double dy = fabs(wy1 - wy2);
    int idx, idy;

    if (spread_dx == 0.0) {
	return;
    }
    if (spread_dy == 0.0) {
	return;
    }
    idx = (int) (dx / spread_dx);
    idy = (int) (dy / spread_dy);
    init_qadd_build(bno);
    for (i = 0; i <= idx; i++) {
	for (j = 0; j <= idy; j++) {
	    if (j % 2) {
		qadd_build(bno, wx1 + spread_dx * (i + 0.5), wy1 + spread_dy * j, 1.0);
		box(wx1 + spread_dx * (i + 0.5), wy1 + spread_dy * j);
	    } else {
		qadd_build(bno, wx1 + spread_dx * i, wy1 + spread_dy * j, 1.0);
		box(wx1 + spread_dx * i, wy1 + spread_dy * j);
	    }
	}
    }
    flush_qadd_build(bno);
}

void spread_random(int bno, double wx1, double wy1, double wx2, double wy2)
{
    int i;
    double x, y, d = 1.0, uniform();
    double dx = fabs(wx1 - wx2);
    double dy = fabs(wy1 - wy2);

    init_qadd_build(bno);
    for (i = 0; i < npts_to_place; i++) {
	qadd_build(bno, x = uniform() * dx + wx1, y = uniform() * dy + wy1, d);
	box(x, y);
    }
    flush_qadd_build(bno);
}

void spread_rotated_rect(int bno, double wx1, double wy1, double wx2, double wy2, double wx3, double wy3)
{
    int i, j;
    double ang, m, mm, x, y, d = 1.0, xt, yt;
    double dx = wx1 - wx2;
    double dy = wy1 - wy2;
    double xi, yi, b, bb;
    int idx, idy;
    double s, c, xtmp, ytmp;

    m = dy / dx;
    mm = -1.0 / m;
    ang = atan(m);
    c = cos(ang);
    s = sin(ang);
    b = wy1 - m * wx1;
    bb = wy3 - mm * wx3;
    xi = (bb - b) / (m - mm);
    yi = mm * xi + bb;
    dx = hypot(wx1 - xi, wy1 - yi);
    dy = hypot(wx3 - xi, wy3 - yi);
    idx = (int) (dx / spread_dx) + 1;
    idy = (int) (dy / spread_dy) + 1;
    init_qadd_build(bno);
    for (i = 0; i <= idx; i++) {
	for (j = 0; j <= idy; j++) {
	    xt = spread_dx * i;
	    yt = spread_dy * j;
	    xtmp = xt * c - yt * s;
	    ytmp = xt * s + yt * c;
	    xtmp = wx1 + xtmp;
	    ytmp = wy1 + ytmp;
	    qadd_build(bno, xtmp, ytmp, 1.0);
	    box(xtmp, ytmp);
	}
    }
    flush_qadd_build(bno);
}

void spread_rotated_rect_offset(int bno, double wx1, double wy1, double wx2, double wy2, double wx3, double wy3)
{
    int i, j;
    double ang, m, mm, x, y, xt, yt, d = 1.0;
    double dx = wx1 - wx2;
    double dy = wy1 - wy2;
    double xi, yi, b, bb;
    int idx, idy;
    double s, c, xtmp, ytmp;

    m = dy / dx;
    mm = -1.0 / m;
    ang = atan(m);
    c = cos(ang);
    s = sin(ang);
    b = wy1 - m * wx1;
    bb = wy3 - mm * wx3;
    xi = (bb - b) / (m - mm);
    yi = mm * xi + bb;
    dx = hypot(wx1 - xi, wy1 - yi);
    dy = hypot(wx3 - xi, wy3 - yi);
    idx = (int) (dx / spread_dx);
    idy = (int) (dy / spread_dy);
    init_qadd_build(bno);
    for (i = 0; i < idx; i++) {
	for (j = 0; j < idy; j++) {
	    if (j % 2) {
		xt = spread_dx * (i + 0.5);
		yt = spread_dy * j;
		xtmp = xt * c - yt * s;
		ytmp = xt * s + yt * c;
		xtmp += wx1;
		ytmp += wy1;
	    } else {
		xt = spread_dx * i;
		yt = spread_dy * j;
		xtmp = xt * c - yt * s;
		ytmp = xt * s + yt * c;
		xtmp += wx1;
		ytmp += wy1;
	    }
	    qadd_build(bno, xtmp, ytmp, 1.0);
	    box(xtmp, ytmp);
	}
    }
    flush_qadd_build(bno);
}

static int cbuild;		/* for communication between compare_build
				 * and elim_build_dupes */

int compare_build(const void *iptr, const void *jptr)
{
    int *i = (int *)iptr;
    int *j = (int *)jptr;
    if (build[cbuild].by[*i] < build[cbuild].by[*j])
	return -1;
    if (build[cbuild].by[*i] > build[cbuild].by[*j])
	return 1;
    if (build[cbuild].bx[*i] < build[cbuild].bx[*j])
	return -1;
    if (build[cbuild].bx[*i] > build[cbuild].bx[*j])
	return 1;
    return 0;
}

void elim_build_dupes(int buildno, double mindist)
{
    int i, cnt;
    double tmp;
    int *ind, *elim, ntoelim = 0;

    cbuild = buildno;
    ind = (int *) calloc(build[buildno].nbuild, sizeof(int));
    elim = (int *) calloc(build[buildno].nbuild, sizeof(int));
    for (i = 0; i < build[buildno].nbuild; i++) {
	ind[i] = i;
    }
    qsort(ind, build[buildno].nbuild, sizeof(int), compare_build);
    cnt = ind[0];
    for (i = 1; i < build[buildno].nbuild; i++) {
	tmp = hypot(build[buildno].bx[cnt] - build[buildno].bx[ind[i]],
		    build[buildno].by[cnt] - build[buildno].by[ind[i]]);
	if (tmp < mindist) {
	    elim[ntoelim] = ind[i];
	    ntoelim++;
	} else {
	    cnt = ind[i];
	}
    }
    if (ntoelim) {
	char tmpbuf[256];
	sprintf(tmpbuf, "Found %d duplicate build points using Minimum distance=%lf, delete?\n", ntoelim, mindist);
	if (yesno(tmpbuf, "Press Yes or No", "Yes", "No")) {
	    delete_build_points(buildno, elim, ntoelim);
	}
    }
    free(ind);
    free(elim);
}

void gen_circ(int bno, double xc, double yc, double r, int nr, int na, int nainc)
{
    double teta;
    int i, j;
    double dteta, dr, rr, xtmp, ytmp;

    dr = r / nr;
    init_qadd_build(bno);
    qadd_build(bno, xc, yc, 1.0);
    for (i = 0; i < nr; ++i) {
	rr = (i + 1) * dr;
	teta = 0.0;
	dteta = M_PI * 2.0 / na;
	for (j = 0; j < na; ++j) {
	    xtmp = xc + rr * cos(teta);
	    ytmp = yc + rr * sin(teta);
	    qadd_build(bno, xtmp, ytmp, 1.0);
	    teta += dteta;
	}
	na += nainc;
    }
    flush_qadd_build(bno);
}

void do_select_build_region(void)
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

void delete_build_inside_region(void)
{
    int *elim, ntoelim = 0, i;

    elim = (int *) calloc(build[curbuild].nbuild, sizeof(int));
    for (i = 0; i < build[curbuild].nbuild; i++) {
	if (inregion(regionx, regiony, nregion, build[curbuild].bx[i], build[curbuild].by[i])) {
	    elim[ntoelim] = i;
	    ntoelim++;
	}
    }
    delete_build_points(curbuild, elim, ntoelim);
    free(elim);
}

void delete_build_outside_region(void)
{
    int *elim, ntoelim = 0, i;

    elim = (int *) calloc(build[curbuild].nbuild, sizeof(int));
    for (i = 0; i < build[curbuild].nbuild; i++) {
	if (!inregion(regionx, regiony, nregion, build[curbuild].bx[i], build[curbuild].by[i])) {
	    elim[ntoelim] = i;
	    ntoelim++;
	}
    }
    delete_build_points(curbuild, elim, ntoelim);
    free(elim);
}

int check_crit(int crit, double a, double avgd, double tol)
{
    double dfact;
    /* check the criteria */
    switch (crit) {
    case 0:
	dfact = a;
	return (dfact > tol);
    case 1:
	dfact = a / avgd;
	return (dfact > tol);
    case 2:
	dfact = a / avgd;
	return (dfact > tol);
    }
    return 0;
}

void auto_spread(int gridno, double x1, double y1, double x2, double y2, double dx, double dy, int crit, int maxloop, double tol, int dgrid, int reg)
{
    int i, j, sind, si, sj, di, dj, k, l, curi, curj, ib, itmp, loop = 1;
    double *dep, wx = x2 - x1, wy = y2 - y1, d = 0.0;
    double a = dx * dy, dfact, avgx, avgy;
    double maxd = -1.0, tmp = -1.0, sumd, avgd, sumx, sumy;
    int dobound = 1, ind, nx, ny, maxi = 0, maxj = 0, notdone = 1, navgd, cnt = 0;
    int errcnt = 0, tcnt = 0;
    int mx1, mx2, my1, my2;
    double bx1, bx2, by1, by2;
    extern double cboundmin, celevmin;
    char *v, buf[256];
    time_t tp, st, et;		/* if timing needed */
    /* initialize the search routine */
    if (grid[gridno].nbounds <= 0) {
	errwin("Need to have an edit boundary, operation cancelled");
	return;
    }
/*
    printf("Auto parms: grid=%d x1=%lf y1=%lf x2=%lf y2=%lf dx=%lf dy=%lf type=%d loop=%d tol=%lf\n", gridno, x1, y1, x2, y2, dx, dy, crit, maxloop, tol);
*/
    write_mode_str("Automatically generating build points...");
    prep_find(dgrid);
    init_qadd_build(0);
    setup_abort();
    st = time(0);
    et = time(0) - st;
    /* dimensions of temporary rectangular grid */
    nx = (int) (wx / dx + 1);
    ny = (int) (wy / dy + 1);
    /* v is a visited array */
    v = (char *) calloc((nx + 1) * (ny + 1), sizeof(char));
    if (v == NULL) {
	errwin("Unable to allocate memory for flags array");
	unsetup_abort();
	return;
    }
    dep = (double *) calloc((nx + 1) * (ny + 1), sizeof(double));
    if (dep == NULL) {
	errwin("Unable to allocate memory for depth array");
	unsetup_abort();
	return;
    }
    if (dep == NULL) {
	errwin("Unable to allocate memory for element array");
	unsetup_abort();
	return;
    }
    create_timer_popup();
    et = time(0) - st;
    sprintf(buf, "Initializing auxillary grid\nStart\nElapsed time %ld seconds", et);
    set_timer_item(buf, 0.0, (double) (nx * ny), 0.0);
    write_mode_str("Initializing auxillary grid, please wait...");
    /* initialize v - a boundary must be defined */
    for (j = 0; j < ny; j++) {
	sind = j * nx;
	if ((j % 20 == 0) && check_action()) {
	    extern int cancel_flag;
	    cancel_action(0);
	    if (cancel_flag)
		goto bailout;
	}
	for (i = 0; i < nx; i++) {
	    if (goodpoint(gridno, x1 + dx * i, y1 + dy * j)) {
		if (!reg) {
		    v[i + sind] = 0;
		} else {
		    if (inregion(regionx, regiony, nregion, x1 + dx * i, y1 + dy * j)) {
			v[i + sind] = 0;
		    } else {
			v[i + sind] = 2;
		    }
		}
	    } else {
		v[i + sind] = 2;
	    }
	}
	et = time(0) - st;
	sprintf(buf, "Initializing auxillary grid\nWorking\nElapsed time %ld seconds", et);
	set_timer_item(buf, 0.0, (double) (nx * ny), (double) (i * j));

    }
    write_mode_str("Computing depths in auxillary grid, please wait...");
    create_timer_popup();
    et = time(0) - st;
    sprintf(buf, "Computing depths in temporary grid\nStart\nElapsed time %ld seconds", et);
    set_timer_item(buf, 0.0, (double) (grid[dgrid].nmel), 0.0);
    errcnt = 0;
    for (k = 0; k < grid[dgrid].nmel; k++) {
	if ((k % 20 == 0) && check_action()) {
	    extern int cancel_flag;
	    cancel_action(0);
	    if (cancel_flag)
		goto bailout;
	}
	get_bounding_box(dgrid, k, &bx1, &bx2, &by1, &by2);
	mx1 = floor((bx1 - x1) / dx);
	mx2 = ceil((bx2 - x1) / dx);
	my1 = floor((by1 - y1) / dy);
	my2 = ceil((by2 - y1) / dy);
	for (j = my1; j <= my2; j++) {
	    sind = j * nx;
	    for (i = mx1; i <= mx2; i++) {
		if (inside_element(dgrid, k, x1 + dx * i, y1 + dy * j)) {
		    if ((j < ny && j >= 0) && (i < nx && i >= 0)) {
			dep[i + sind] = get_depth_element(dgrid, k, x1 + dx * i, y1 + dy * j);
			if (dep[i + sind] < celevmin) {
			    dep[i + sind] = celevmin;
			    errcnt++;
			}
		    } else {
			printf("Assumption failed\n");
		    }
		}
	    }
	}
	if (k % 20 == 0) {
	    et = time(0) - st;
	    sprintf(buf, "Computing depths in temporary grid\nAdjusted %d elevations\nElapsed time %ld seconds", errcnt, et);
	    set_timer_item(buf, 0.0, 1.0 * grid[dgrid].nmel, 1.0 * k);
	}
    }
    if (!reg) {
	write_mode_str("Processing boundary, please wait...");
	/* visit each of the boundary points and mark them */
	create_timer_popup();
	for (i = 0; i < grid[gridno].nbounds; i++) {
	    ib = grid[gridno].boundaries[i];
	    if (boundary[ib].bactive) {
		et = time(0) - st;
		sprintf(buf, "Start\nElapsed time %ld seconds", et);
		set_timer_item(buf, 0.0, (double) (boundary[ib].nbpts), 0.0);
		errcnt = 0;
		for (j = 0; j < boundary[ib].nbpts; j++) {
		    et = time(0) - st;
		    sprintf(buf, "Processing boundary\nAdjusted %d elevations\nElapsed time %ld seconds", errcnt, et);
		    set_timer_item(buf, 0.0, (double) (boundary[ib].nbpts), (double) (j));
		    /* find the nearest temporary grid vertex */
		    maxi = (int) ((boundary[ib].boundx[j] - x1) / dx) + 1;
		    maxj = (int) ((boundary[ib].boundy[j] - y1) / dy) + 1;
		    curi = si = maxi;
		    curj = sj = maxj;
		    /* find the element this boundary point is in */
		    find_element(dgrid, boundary[ib].boundx[j], boundary[ib].boundy[j], &ind);
		    if (ind < 0) {
			setcolor(2);
			my_circlefilled(boundary[ib].boundx[j], boundary[ib].boundy[j]);
			fprintf(stderr, "Boundary point not found at init\n");
			setcolor(1);
			goto skip1;
		    }
		    /* get the depth at this boundary point */
		    maxd = get_depth_element(dgrid, ind, boundary[ib].boundx[j], boundary[ib].boundy[j]);
		    if (maxd < cboundmin) {
			maxd = cboundmin;
		    }
		    cnt = 0;
		    /* initialize accumulation variables */
		    a = dx * dy;
		    navgd = 1;
		    sumd = avgd = maxd;
		    sumx = avgx = x1 + curi * dx;
		    sumy = avgy = y1 + curj * dy;
		    /* mark the initial point as visited */
		    v[curi + curj * nx] = 1;
		    cnt++;
		    if (check_crit(crit, a, avgd, tol)) {
			goto bustout1;
		    }
		    curi--;
		    curj--;
		    di = 1;
		    dj = 1;
		    loop = 1;
		    my_move2(x1 + curi * dx, y1 + curj * dy);
		    while (loop < maxloop) {
			for (l = 0; l < di; l++) {
			    my_draw2(x1 + curi * dx, y1 + curj * dy);
			    if ((l % 10 == 0) && check_action()) {
				extern int cancel_flag;
				cancel_action(0);
				if (cancel_flag)
				    goto bailout;
			    }
			    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
				if (v[curi + curj * nx] < 2) {
				    a += dx * dy;
				    navgd++;
				    sumd += dep[curi + curj * nx];
				    avgd = sumd / navgd;
				    sumx += x1 + curi * dx;
				    sumy += y1 + curj * dy;
				    avgx = sumx / navgd;
				    avgy = sumy / navgd;
				    v[curi + curj * nx] = 1;
				    cnt++;
				    if (check_crit(crit, a, avgd, tol)) {
					goto bustout1;
				    }
				}
			    }
			    curi++;
			}
			for (l = 0; l < dj; l++) {
			    if ((l % 10 == 0) && check_action()) {
				extern int cancel_flag;
				cancel_action(0);
				if (cancel_flag)
				    goto bailout;
			    }
			    my_draw2(x1 + curi * dx, y1 + curj * dy);
			    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
				if (v[curi + curj * nx] < 2) {
				    a += dx * dy;
				    navgd++;
				    sumd += dep[curi + curj * nx];
				    avgd = sumd / navgd;
				    sumx += x1 + curi * dx;
				    sumy += y1 + curj * dy;
				    avgx = sumx / navgd;
				    avgy = sumy / navgd;
				    v[curi + curj * nx] = 1;
				    cnt++;
				    if (check_crit(crit, a, avgd, tol)) {
					goto bustout1;
				    }
				}
			    }
			    curj++;
			}
			for (l = 0; l < di; l++) {
			    if ((l % 10 == 0) && check_action()) {
				extern int cancel_flag;
				cancel_action(0);
				if (cancel_flag)
				    goto bailout;
			    }
			    my_draw2(x1 + curi * dx, y1 + curj * dy);
			    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
				if (v[curi + curj * nx] < 2) {
				    a += dx * dy;
				    navgd++;
				    sumd += dep[curi + curj * nx];
				    avgd = sumd / navgd;
				    sumx += x1 + curi * dx;
				    sumy += y1 + curj * dy;
				    avgx = sumx / navgd;
				    avgy = sumy / navgd;
				    v[curi + curj * nx] = 1;
				    dfact = sqrt(a) / avgd;
				    cnt++;
				    if (check_crit(crit, a, avgd, tol)) {
					goto bustout1;
				    }
				}
			    }
			    curi--;
			}
			for (l = 0; l < dj; l++) {
			    if ((l % 10 == 0) && check_action()) {
				extern int cancel_flag;
				cancel_action(0);
				if (cancel_flag)
				    goto bailout;
			    }
			    my_draw2(x1 + curi * dx, y1 + curj * dy);
			    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
				if (v[curi + curj * nx] < 2) {
				    a += dx * dy;
				    navgd++;
				    sumd += dep[curi + curj * nx];
				    avgd = sumd / navgd;
				    sumx += x1 + curi * dx;
				    sumy += y1 + curj * dy;
				    avgx = sumx / navgd;
				    avgy = sumy / navgd;
				    v[curi + curj * nx] = 1;
				    cnt++;
				    if (check_crit(crit, a, avgd, tol)) {
					goto bustout1;
				    }
				}
			    }
			    curj--;
			}
			loop++;
			curi = si - loop;
			curj = sj - loop;
			di = 2 * (loop - 1) + 1;
			dj = 2 * (loop - 1) + 1;
		    }
	    bustout1:;
/* add build point */
		    if (ind >= 0) {
			double dtmp = get_depth_element(dgrid, ind,
						     boundary[ib].boundx[j],
						    boundary[ib].boundy[j]);
			if (dtmp < cboundmin) {
			    dtmp = cboundmin;
			}
			qadd_build(0, boundary[ib].boundx[j], boundary[ib].boundy[j], dtmp);
		    } else {
			setcolor(2);
			my_circlefilled(boundary[ib].boundx[j], boundary[ib].boundy[j]);
			fprintf(stderr, "Boundary point not found in build\n");
			setcolor(1);
		    }
	    skip1:    ;
		}
	    }
	}
    }
    /* END of boundary */

    write_mode_str("Processing interior of the domain, please wait...");
    et = time(0) - st;
    sprintf(buf, "Processing interior of the domain\nStart...\nElapsed time %ld seconds", et);
    set_timer_item(buf, 0.0, (double) (nx * ny), 0.0);
    /* now do the interior of the edit grid */
    while (notdone) {
	notdone = 0;
	maxd = -1.0;
	tcnt = 0;
	for (j = 0; j < ny; j++) {
	    sind = j * nx;
	    for (i = 0; i < nx; i++) {
		if (v[i + sind] == 0) {
		    tmp = dep[i + sind];
		    if (tmp > maxd) {
			maxd = tmp;
			maxi = i;
			maxj = j;
		    }
		    notdone = 1;
		} else {
		    tcnt++;
		}
	    }
	}
	et = time(0) - st;
	sprintf(buf, "Processing interior of the domain\nWorking...\nElapsed time %ld seconds", et);
	set_timer_item(buf, 0.0, (double) (nx * ny), (double) (tcnt));
	curi = si = maxi;
	curj = sj = maxj;
	cnt = 0;
	a = dx * dy;
	navgd = 1;
	sumd = avgd = maxd;
	sumx = avgx = x1 + curi * dx;
	sumy = avgy = y1 + curj * dy;
	v[curi + curj * nx] = 1;
	cnt++;
	if (check_crit(crit, a, avgd, tol)) {
	    goto bustout;
	}
	curi--;
	curj--;
	di = 1;
	dj = 1;
	loop = 1;
	if (notdone) {
	    my_move2(x1 + curi * dx, y1 + curj * dy);
	    while (loop < maxloop) {
		for (l = 0; l < di; l++) {
		    if ((l % 10 == 0) && check_action()) {
			extern int cancel_flag;
			cancel_action(0);
			if (cancel_flag)
			    goto bailout;
		    }
		    my_draw2(x1 + curi * dx, y1 + curj * dy);
		    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
			if (!v[curi + curj * nx]) {
			    a += dx * dy;
			    navgd++;
			    sumd += dep[curi + curj * nx];
			    avgd = sumd / navgd;
			    sumx += x1 + curi * dx;
			    sumy += y1 + curj * dy;
			    avgx = sumx / navgd;
			    avgy = sumy / navgd;
			    v[curi + curj * nx] = 1;
			    cnt++;
			    if (check_crit(crit, a, avgd, tol)) {
				goto bustout;
			    }
			}
		    }
		    curi++;
		}
		for (l = 0; l < dj; l++) {
		    if ((l % 10 == 0) && check_action()) {
			extern int cancel_flag;
			cancel_action(0);
			if (cancel_flag)
			    goto bailout;
		    }
		    my_draw2(x1 + curi * dx, y1 + curj * dy);
		    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
			if (!v[curi + curj * nx]) {
			    a += dx * dy;
			    navgd++;
			    sumd += dep[curi + curj * nx];
			    avgd = sumd / navgd;
			    sumx += x1 + curi * dx;
			    sumy += y1 + curj * dy;
			    avgx = sumx / navgd;
			    avgy = sumy / navgd;
			    v[curi + curj * nx] = 1;
			    cnt++;
			    if (check_crit(crit, a, avgd, tol)) {
				goto bustout;
			    }
			}
		    }
		    curj++;
		}
		for (l = 0; l < di; l++) {
		    if ((l % 10 == 0) && check_action()) {
			extern int cancel_flag;
			cancel_action(0);
			if (cancel_flag)
			    goto bailout;
		    }
		    my_draw2(x1 + curi * dx, y1 + curj * dy);
		    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
			if (!v[curi + curj * nx]) {
			    a += dx * dy;
			    navgd++;
			    sumd += dep[curi + curj * nx];
			    avgd = sumd / navgd;
			    sumx += x1 + curi * dx;
			    sumy += y1 + curj * dy;
			    avgx = sumx / navgd;
			    avgy = sumy / navgd;
			    v[curi + curj * nx] = 1;
			    cnt++;
			    if (check_crit(crit, a, avgd, tol)) {
				goto bustout;
			    }
			}
		    }
		    curi--;
		}
		for (l = 0; l < dj; l++) {
		    if ((l % 10 == 0) && check_action()) {
			extern int cancel_flag;
			cancel_action(0);
			if (cancel_flag)
			    goto bailout;
		    }
		    my_draw2(x1 + curi * dx, y1 + curj * dy);
		    if ((curi >= 0 && curi < nx) && (curj >= 0 && curj < ny)) {
			if (!v[curi + curj * nx]) {
			    a += dx * dy;
			    navgd++;
			    sumd += dep[curi + curj * nx];
			    avgd = sumd / navgd;
			    sumx += x1 + curi * dx;
			    sumy += y1 + curj * dy;
			    avgx = sumx / navgd;
			    avgy = sumy / navgd;
			    v[curi + curj * nx] = 1;
			    cnt++;
			    if (check_crit(crit, a, avgd, tol)) {
				goto bustout;
			    }
			}
		    }
		    curj--;
		}
		if (cnt < (2 * loop) && loop > 1) {
		    goto skip;
		}
		loop++;
		curi = si - loop;
		curj = sj - loop;
		di = 2 * (loop - 1) + 1;
		dj = 2 * (loop - 1) + 1;
	    }
    bustout:;
	    if (notdone) {
/* add build point */
		ind = find_element3(dgrid, avgx, avgy);
		if (ind >= 0) {
		    double dtmp = get_depth_element(dgrid, ind, avgx, avgy);
		    if (dtmp < celevmin) {
			dtmp = celevmin;
		    }
		    qadd_build(0, avgx, avgy, dtmp);
		}
	    }
    skip:    ;
	}
    }
    flush_qadd_build(0);
bailout:;
    unsetup_abort();
    end_find();
    free(v);
    free(dep);
}

void do_random_placement(void)
{
    double u1, u2, u3;
    double x1, x2, y1, y2, z1, z2, dx, dy, dz, x, y, d;
    double drand();
    int cnt, i, j, el;
    int gridno = 0;
    init_qadd_build(0);
    cnt = 0;
    x1 = grid[gridno].xmin;
    y1 = grid[gridno].ymin;
    z1 = -grid[gridno].dmax;
    x2 = grid[gridno].xmax;
    y2 = grid[gridno].ymax;
    z2 = -grid[gridno].dmin;
    dx = x2 - x1;
    dy = y2 - y1;
    dz = z2 - z1;
    while (cnt < 1000) {
	u1 = drand48() * dx + x1;
	u2 = drand48() * dy + y1;
	u3 = drand48() * dz + z1;
	x = u1;
	y = u2;
	if (goodpoint(gridno, x, y)) {
	    find_element(MAXGRIDS, x, y, &el);
	    if (el >= 0) {
		d = -get_depth_element(MAXGRIDS, el, x, y);
		if (u3 <= d) {
		    my_circlefilled(x, y);
		    qadd_build(0, x, y, -d);
		    cnt++;
		}
	    }
	}
    }
    flush_qadd_build(0);
}

/*
double potential();

void emf(bx, by, n)
double bx[], by[];
int n;
{
    int k;
    for (k = 0; k < nsteps; k++) {
	sum = 0.0;
	prevsum = sum;
	track(n, px, py, u, v, au, av, dt, &sum);
	printf("%d %lf\n", k, sum);
    }
    printf("&\n");
    for (i = 0; i < n; i++) {
	printf("%lf %lf\n", px[i], py[i]);
    }
    printf("&\n");
}

track(n, x, y, u, v, au, av, dt, sum)
    double x[], y[], u[], v[], au[], av[], dt, *sum;
{
    int i, j;
    double xn, yn, dx, dy, utmp, vtmp, dt2 = dt * dt;
    for (i = 0; i < n; i++) {
	x[i] = x[i] + u[i] * dt + 0.5 * au[i] * dt2;
	y[i] = y[i] + v[i] * dt + 0.5 * av[i] * dt2;
    }
    for (i = 0; i < n; i++) {
	u[i] = u[i] + 0.5 * au[i] * dt;
	v[i] = v[i] + 0.5 * av[i] * dt;
	au[i] = 0.0;
	av[i] = 0.0;
    }
    *sum = 0.0;
    for (i = 0; i < n - 1; i++) {
	for (j = i + 1; j < n; j++) {
	    *sum = *sum + potential(i, j, &utmp, &vtmp);
	    au[i] = au[i] + utmp;
	    av[i] = av[i] + vtmp;
	    au[j] = au[j] - utmp;
	    av[j] = av[j] - vtmp;
	}
    }
    for (i = 0; i < n; i++) {
	u[i] = u[i] + 0.5 * au[i] * dt;
	v[i] = v[i] + 0.5 * av[i] * dt;
    }
}

double potential(px, py, i, j, u, v)
    int i, j;
    double *px, *py, *u, *v;
{
    double fx, fy, dist, force, ang, ri;
    dist = hypot(px[i] - px[j], py[i] - py[j]);
    force = (sigma[i] + sigma[j]) - dist;
    if (force > 10.0) {
	force = 10.0;
    } else if (force < -10.0) {
	force = -10.0;
    }
    ang = atan2(py[i] - py[j], px[i] - px[j]);
    *u = cos(ang) * force;
    *v = sin(ang) * force;
    return force;
}
*/

void WriteBuildB(Buildpts * b, char *fname)
{
    int i, npts;
    char buf[256];
    FILE *fp;
    float *ftmp;
    int magic = 40;

    if ((fp = fopen(fname, "w")) == NULL) {
	fprintf(stderr, "WriteBuild(): unable to open file\n");
	return;
    }
    npts = b->nbuild;
    ftmp = (float *) malloc(b->nbuild * sizeof(float));
    fwrite(&magic, sizeof(int), 1, fp);
    fwrite(&b->nbuild, sizeof(int), 1, fp);
    for (i = 0; i < b->nbuild; i++) {
	ftmp[i] = b->bx[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < b->nbuild; i++) {
	ftmp[i] = b->by[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    for (i = 0; i < b->nbuild; i++) {
	ftmp[i] = b->db[i];
    }
    fwrite(ftmp, sizeof(float), npts, fp);
    free(ftmp);
    fclose(fp);
}

void WriteBuildA(Buildpts * b, char *fname)
{
    int i;
    FILE *fp;
    if ((fp = fopen(fname, "w")) == NULL) {
	fprintf(stderr, "WriteBuildA(): unable to open file\n");
	return;
    }
    for (i = 0; i < b->nbuild; i++) {
	fprintf(fp, "%d %lf %lf %lf\n", i + 1, b->bx[i], b->by[i], b->db[i]);
    }
    fclose(fp);
}

Buildpts *NewBuild(int npts)
{
    Buildpts *b = (Buildpts *) malloc(sizeof(Buildpts));
    if (b == NULL) {
	return NULL;
    }
    b->nbuild = npts;
    b->buildactive = 1;
    b->bx = (double *) calloc(npts, sizeof(double));
    b->by = (double *) calloc(npts, sizeof(double));
    b->db = (double *) calloc(npts, sizeof(double));
    return b;
}

void DefaultsBuild(Buildpts * b)
{
    b->color = 1;
    b->sym = 2;
    b->size = 1.0;
    b->isol = 0;
}

void ReallocateBuild(Buildpts * b, int npts)
{
    b->nbuild = npts;
    b->bx = (double *) realloc(b->bx, npts * sizeof(double));
    b->by = (double *) realloc(b->by, npts * sizeof(double));
    b->db = (double *) realloc(b->db, npts * sizeof(double));
}

void FreeBuild(Buildpts * b)
{
    if (b != NULL) {
	if (b->bx != NULL) {
	    free(b->bx);
	}
	if (b->by != NULL) {
	    free(b->by);
	}
	if (b->db != NULL) {
	    free(b->db);
	}
	free(b);
    }
}

void AddBuildPoints(Buildpts * b, int n, double *wx, double *wy, double *d)
{
    int i;
    int npts = b->nbuild;

    if (npts) {
	b->nbuild += n;
	b->bx = (double *) realloc(b->bx, (npts + n) * sizeof(double));
	b->by = (double *) realloc(b->by, (npts + n) * sizeof(double));
	b->db = (double *) realloc(b->db, (npts + n) * sizeof(double));
	for (i = npts; i < npts + n; i++) {
	    b->bx[i] = wx[i - npts];
	    b->by[i] = wy[i - npts];
	    b->db[i] = d[i - npts];
	}
    } else {
	b->nbuild = n;
	b->bx = (double *) calloc((npts + n), sizeof(double));
	b->by = (double *) calloc((npts + n), sizeof(double));
	b->db = (double *) calloc((npts + n), sizeof(double));
	for (i = 0; i < n; i++) {
	    b->bx[i] = wx[i];
	    b->by[i] = wy[i];
	    b->db[i] = d[i];
	}
    }
}

void DrawBuildDepths(Buildpts * b, Isolparms ip)
{
    int i, k;
    double cc, ccp1, e;
    extern int mapisolconc[], mapisolbath[];
    for (i = 0; i < b->nbuild; i++) {
	e = b->db[i];
	for (k = 0; k < ip.nisol - 1; k++) {
	    cc = ip.cis[k];
	    ccp1 = ip.cis[k + 1];
	    if (e > cc && e <= ccp1) {
		setcolor(mapisolbath[k]);
		switch (b->sym) {
		case 0:
		    break;
		case 1:
		    drawdot(b->bx[i], b->by[i]);
		    break;
		case 2:
		    solidbox(b->bx[i], b->by[i]);
		    break;
		case 3:
		    my_circlefilled(b->bx[i], b->by[i]);
		    break;
		}
		break;
	    }
	}
    }
}
