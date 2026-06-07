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
 * action procs
 */

#include "symdefs.h"
#include "defines.h"
#include "globals.h"

#ifndef lint
static char RCSid[] = "$Id: drawobjs.c,v 1.13 2007/03/21 18:06:25 pturner Exp $";
#endif

double setcharsize();

/* TODO */
extern int hardcopyflag;
extern int cursortype;
extern int cursor_oldx;
extern int cursor_oldy;
extern int display_elevdepth;
int force_redraw;

void do_animat();
void freshen_display(int restart);
void display_image(void);
void drawgraph(void);
double eval_teanl_node(int flowno, int node, double t);
double compute_ela_mass(int gno, int cno, int flowno, int step, double t);
void boxplot(int gno, framep f);
void draw_string(int gno, int i);
void draw_annotation(int gno);
void domapscale(int gno);
void dovscale(int gno);
void dowscale(int gno);
void dofscale(int gno);
void dolegend(Isolparms ip, int *cmap, int nmap);
void drawtimeline(int gno, double time);
double get_current_time(void);
void draw_ref_point(int gno);
void plotone(int gno, int ino, double curtime, int zno);
void do_drawmapscale();
void do_drawclock(int gno, double time);
void drawtimeinfo(int gno, double time);
void draw_region(int r);
void draw_single_region(double *rx, double *ry, int nr);
void draw_transect(int n, double *x, double *y);

extern int slice_adcircelev;

extern int tdevice;

static int restart = 0;

void setbgcolor(int bg)
{
    bgcolor = bg;
}

void setredraw_world(void)
{
    int restart = 0, sv = no_display;
    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    no_display = 0;
    doclear = 1;
    dobackground = 1;
    display_image();
    no_display = sv;
    if (restart) {
	setirun();
    }
}

void setreset_world(void)
{
    int restart = 0, sv = no_display;
    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    reset_world(cg);
    no_display = 0;
    doclear = 1;
    dobackground = 1;
    display_image();
    no_display = sv;
    if (restart) {
	setirun();
    }
}

void do_hardcopy(void)
{
    hardcopyflag = TRUE;
    initgraphics(hdevice);
    set_up_world(cg);
    display_image();
    hardcopyflag = FALSE;
    my_doublebuffer(0);
    initgraphics(0);
    display_image();
    my_doublebuffer(1);
}

/*
 * do a single plot then quit 
 */
void do_batch_plot(void)
{
    int i, j;
    extern int drawimage_flag;
    hardcopyflag = TRUE;
    inwin = 1;
    initgraphics(7);
    set_up_world(cg);
    if (drawimage_flag) {
	drawimage();
    }
    for (i = 0; i < maxgraph; i++) {
	if (isactive_graph(i) && !g[i].hidden) {
	    plotone(i, get_current_step(), get_current_time(), -1);
	    if (debuglevel == 10) {
		printf("%d %lf\n", get_current_step(), get_current_time());
	    }
	    for (j = 0; j < MAXZOOMBOXES; j++) {
		if (display_object(i, ZOOM, j)) {
		    draw_zoom(i, j);
		    plotone(i, timeclock.curstep, timeclock.curtime, j);
		    defineworld(g[i].w.xg1, g[i].w.yg1, g[i].w.xg2, g[i].w.yg2, islogx(i), islogy(i));
		    viewport(g[i].v.xv1, g[i].v.yv1, g[i].v.xv2, g[i].v.yv2);
		}
	    }
	    draw_annotation(i);
	}
    }
    draw_annotation(-1);
    leavegraphics();
    hardcopyflag = FALSE;
}

/*
 * draw all active graphs, when graphs are drawn, draw the focus markers
 */
void drawgraph(void)
{				/* TODO go through and change all the
				 * drawgraph()s to display_image */
    display_image();
}

void display_image(void)
{
    int i, j;
    extern int drawimage_flag;
    if (hardcopyflag) {
    } else {
	do_set_window();
    }
    if (no_display) {
	return;
    }
    if (timeclock.running) {
	my_doublebuffer(1);
    } else {
	my_doublebuffer(0);
    }
    if (inwin && (auto_redraw || force_redraw)) {
	if (cursortype) {
	    cursor_oldx = cursor_oldy = -1;
	}
/*
   set_wait_cursor(NULL);
   set_right_footer("Redraw...");
 */
	if (hardcopyflag) {
	    initgraphics(hdevice);
	} else {
	    initgraphics(tdevice);
	}
	if (drawimage_flag) {
	    drawimage();
	}
	for (i = 0; i < maxgraph; i++) {
	    if (isactive_graph(i) && !g[i].hidden) {
		plotone(i, get_current_step(), get_current_time(), -1);
		if (debuglevel == 10) {
		    printf("%d %lf\n", get_current_step(), get_current_time());
		}
		for (j = 0; j < MAXZOOMBOXES; j++) {
		    if (display_object(i, ZOOM, j)) {
			draw_zoom(i, j);
			plotone(i, timeclock.curstep, timeclock.curtime, j);
			defineworld(g[i].w.xg1, g[i].w.yg1, g[i].w.xg2, g[i].w.yg2, islogx(i), islogy(i));
			viewport(g[i].v.xv1, g[i].v.yv1, g[i].v.xv2, g[i].v.yv2);
		    }
		}
		draw_annotation(i);
	    }
	}
	draw_annotation(-1);
	leavegraphics();
	set_steps_label(timeclock.curstep);
	defineworld(g[cg].w.xg1, g[cg].w.yg1, g[cg].w.xg2, g[cg].w.yg2, islogx(cg), islogy(cg));
	viewport(g[cg].v.xv1, g[cg].v.yv1, g[cg].v.xv2, g[cg].v.yv2);
	if (!timeclock.running) {
	    my_doublebuffer(0);
	    draw_focus(cg);
	}
    }
}

void plotone(int gno, int ino, double curtime, int zno)
{
    extern double mapscale;
    char buf[256];
    int i, j;
    int st;
    double x, y;
    double *c, dy = (g[gno].w.yg2 - g[gno].w.yg1);
    double x1 = g[gno].w.xg1, x2 = g[gno].w.xg2, y1 = g[gno].w.yg1, y2 = g[gno].w.yg2;
    world w;
    view v;
    framep f;
    labels lab;

    if (zno < 0) {
	get_graph_view(gno, &v);
	get_graph_labels(gno, &lab);
	if (hardcopyflag && mapscale != 1.0) {
	    set_up_mapscale(gno, mapscale);
	}
	get_graph_world(gno, &w);
	get_graph_framep(gno, &f);
	setclipping(TRUE);
	defineworld(w.xg1, w.yg1, w.xg2, w.yg2, islogx(gno), islogy(gno));
	viewport(v.xv1, v.yv1, v.xv2, v.yv2);
    } else {
	setclipping(TRUE);
    }

    st = 0;
    if (!inwin) {
	return;
    }
    setcharsize(0.8);
    setfont(4);

    if (zno < 0) {
	if (g[gno].f.active == ON) {
	    if (g[gno].f.fillbg == ON) {
		setcolor(g[gno].f.bgcolor);
		fillrectcolor(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2);
	    }
	}
    }
    if (dobackground) {
	setbgcolor(bgcolor);
	for (i = 0; i < MAXGRIDS; i++) {
	    if (object_isactive(GRID, i) && display_object(gno, BATH, i)) {
		if (grid[i].depthflag) {
		    do_isol2(i, 0, grid[i].edepth, g[gno].grid[i].ip, mapisolbath, g[gno].grid[i].ip.nisol);
		} else {
		    do_isol(i, 0, grid[i].depth, g[gno].grid[i].ip, mapisolbath, g[gno].grid[i].ip.nisol);
		}
		if (g[gno].grid[i].ip.lactive == ON) {
		    dolegend(g[gno].grid[i].ip, mapisolbath, g[gno].grid[i].ip.nisol);
		}
	    }
	}
    }
    if (zno < 0) {
	setcolor(1);
	setlinestyle(1);
	setlinewidth(1);
	drawaxes(gno);

	setcolor(1);
	setlinestyle(1);
	setlinewidth(1);
	setcharsize(0.8);
	setfont(4);
    }
    setcolor(1);
    setlinewidth(1);
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, MAG, i)) {
	    drawflowmag(gno, i, curtime);
	    if (g[gno].flowf[i].magip.lactive == ON) {
		dolegend(g[gno].flowf[i].magip, mapisolconc, g[gno].flowf[i].magip.nisol);
	    }
	}
    }
    setcolor(1);
    setlinewidth(1);
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, WETDRY, i)) {
	    drawwetdry(gno, i, curtime);
	}
    }
    setcolor(1);
    setlinewidth(1);
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, ELEV, i)) {
	    drawelevations(gno, i, curtime);
	    if (g[gno].flowf[i].elevip.lactive == ON) {
		dolegend(g[gno].flowf[i].elevip, mapisolconc, g[gno].flowf[i].elevip.nisol);
	    }
	}
    }
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, ELEV, i) && display_object(gno, ADCIRC, ELEV, i)) {
	    int j;
	    double *tmp = (double *) malloc(flowt[i].g.nmnp * sizeof(double));
	    if (ino >= flowt[i].nsteps)
		continue;
	    if (tmp != NULL) {
		if (g[gno].flowt[i].display_elevdepth == ON) {
		    for (j = 0; j < flowt[i].g.nmnp; j++) {
			tmp[j] = flowt[i].f[ino].e[j] + flowt[i].g.depth[j];
		    }
		} else {
		    for (j = 0; j < flowt[i].g.nmnp; j++) {
			tmp[j] = flowt[i].f[ino].e[j];
			/*printf("%d %lf\n", j, tmp[j]); */
		    }
		}
		do_grid_isol(&flowt[i].g, 0, tmp, g[gno].flowt[i].elevip, mapisolconc, g[gno].flowt[i].elevip.nisol);
		free(tmp);
	    }
	    if (g[gno].flowt[i].elevip.lactive == ON) {
		dolegend(g[gno].flowt[i].elevip, mapisolconc, g[gno].flowt[i].elevip.nisol);
	    }
	}
    }
    setcolor(1);
    setlinewidth(1);
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, FLOW, i) && display_object(gno, ADCIRC, MAG, i)) {
	    if (ino >= flowt[i].nsteps)
		continue;
	    drawflowtmag(gno, i, ino);
	    if (g[gno].flowt[i].magip.lactive == ON) {
		dolegend(g[gno].flowt[i].magip, mapisolconc, g[gno].flowt[i].magip.nisol);
	    }
	}
    }
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, ELEV, i) && display_object(gno, ADCIRC, ELEV, MAXP, i)) {
	    if (compute_adcirc_maxelev(i, 0)) {
		double *tmp;
		tmp = (double *) malloc(flowt[i].g.nmnp * sizeof(double));
		if (ino >= flowt[i].nsteps)
		    continue;
		if (tmp != NULL) {
		    for (j = 0; j < flowt[i].g.nmnp; j++) {
			tmp[j] = flowt[i].global_emax[j] + g[gno].flowt[i].d;
		    }
		    do_grid_isol(&flowt[i].g, 0, tmp, g[gno].flowt[i].maxelevip, mapisolconc, g[gno].flowt[i].elevip.nisol);
		    if (g[gno].flowt[i].maxelevip.lactive == ON) {
			dolegend(g[gno].flowt[i].maxelevip, mapisolconc, g[gno].flowt[i].maxelevip.nisol);
		    }
		    free(tmp);
		}
	    }
	}
	if (object_isactive(ADCIRC, ELEV, i) && display_object(gno, ADCIRC, ELEV, NODES, i)) {
	    if (compute_adcirc_maxelev(i, 0)) {
		if (ino >= flowt[i].nsteps)
		    continue;
		for (j = 0; j < flowt[i].g.nmnp; j++) {
		    sprintf(buf, "%lf", flowt[i].global_emax[j]);
		    x = flowt[i].g.xord[j];
		    y = flowt[i].g.yord[j];
		    writestr(x, y, 0, 0, buf);
		}
	    }
	}
    }

    setcolor(1);
    setlinewidth(1);

    setcolor(1);
    setlinewidth(1);
    for (i = 0; i < MAXGRIDS; i++) {
	if (object_isactive(GRID, i)) {
	    setcolor(1);
	    setlinestyle(1);
	    setlinewidth(1);
	    if (g[gno].grid[i].display_gridf == ON) {
		drawgrid_filled(i, g[gno].grid[i].p.fillcol);
	    }
	    if (display_object(gno, GRID, GRID, i)) {
		drawgrid(gno, i, 1.0);
	    }
	    if (g[gno].grid[i].display_courant == ON) {
		drawgrid_cour_filled(gno, i);
	    }
	    if (g[gno].grid[i].display_dimw == ON) {
		drawgrid_dimw_filled(i);
	    }
	    if (g[gno].grid[i].display_nodes == ON) {
		drawgridnodes(i);
	    }
	    if (g[gno].grid[i].display_elements == ON) {
		drawgridelems(i);
	    }
	    if (g[gno].grid[i].display_depths == ON) {
		drawnodedepths(i);
	    }
	    if (display_object(gno, GRID, BOUNDARY, i)) {
		draw_boundary(gno, i);
	    }
	}
    }
    setcolor(1);
    setlinestyle(1);
    setlinewidth(1);

    /* draw custom boundaries (not attached to any grid) */
    for (i = 0; i < MAXBOUNDS; i++) {
	if (object_isactive(BOUNDARY, i) && display_object(gno, BOUNDARY, i)) {
	    draw_boundary2(gno, i);
	}
    }

    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, PHASE, i)) {
	    if (g[gno].flowf[i].freq < flowf[i].nfreq) {
		do_isol(flowf[i].grid, 0, flowf[i].elphase[g[gno].flowf[i].freq], g[gno].flowf[i].phaseip, mapisolconc, g[gno].flowf[i].phaseip.nisol);
		if (g[gno].flowf[i].phaseip.lactive == ON) {
		    dolegend(g[gno].flowf[i].phaseip, mapisolconc, g[gno].flowf[i].phaseip.nisol);
		}
	    }
	}
    }
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, AMP, i)) {
	    if (g[gno].flowf[i].freq < flowf[i].nfreq) {
		do_isol(flowf[i].grid, 0, flowf[i].elamp[g[gno].flowf[i].freq], g[gno].flowf[i].ampip, mapisolconc, g[gno].flowf[i].ampip.nisol);
		if (g[gno].flowf[i].ampip.lactive == ON) {
		    dolegend(g[gno].flowf[i].ampip, mapisolconc, g[gno].flowf[i].ampip.nisol);
		}
	    }
	}
    }
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i) && display_object(gno, TEANL, FLOW, i)) {
	    setcolor(g[gno].flowf[i].p.color);
	    if (g[gno].flowf[i].display == CENTER) {
		drawflowcenter(gno, i, curtime);
	    } else if (g[gno].flowf[i].display == NODES) {
		drawflownodes(gno, i, curtime);
	    } else if (g[gno].flowf[i].display == ELLIPSE) {
		drawflowf_ellipse(gno, i, ino);
	    }
	}
    }
    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, FLOW, i) && display_object(gno, ADCIRC, FLOW, i)) {
	    if (ino >= flowt[i].nsteps)
		continue;
	    setcolor(g[gno].flowt[i].p.color);
	    if (g[gno].flowt[i].display == CENTER) {
		drawflowtcenter(gno, i, ino);
	    } else if (g[gno].flowt[i].display == NODES) {
		drawflowtnodes(gno, i, ino);
	    } else if (g[gno].flowt[i].display == ELLIPSE) {
		drawflowt_ellipse(gno, i, ino);
	    }
	}
    }
    for (i = 0; i < nsta; i++) {
	if (object_isactive(STATION, i) && display_object(gno, STATION, i)) {
	    DrawStation(gno, i);
	}
    }
    for (i = 0; i < ntidestat; i++) {
	if (object_isactive(TIDESTATION, i) && display_object(gno, TIDESTATION, i)) {
	    drawtidestation(gno, i, ino, curtime);
	}
    }
    if (zno < 0) {
	for (i = 0; i < MAXTEANL; i++) {
	    if (object_isactive(TEANL, i) && display_object(gno, TEANL, MARKERS, i)) {
		setcolor(g[gno].flowf[i].p.color);
		setcolor(1);
		draw_elevationmarkers(gno, ino, i, curtime, 0);
	    }
	}
	for (i = 0; i < MAXADCIRC; i++) {
	    if (object_isactive(ADCIRC, ELEV, i) && display_object(gno, ADCIRC, MARKERS, i)) {
		if (ino >= flowt[i].nsteps)
		    continue;
		setcolor(g[gno].flowt[i].p.color);
		setcolor(1);
		draw_elevationmarkers(gno, ino, i, curtime, 1);
	    }
	}
	for (i = 0; i < MAXTRANS; i++) {
	    if (object_isactive(TRANSECT, i) && display_object(gno, TRANSECT, i)) {
		draw_trans(gno, ino, i, curtime, 1);
	    }
	}
	for (i = 0; i < MAXADCIRC3D; i++) {
	    if (object_isactive(ADCIRC3DFLOW, i) && display_object(gno, ADCIRC3DFLOW, FLOW, i)) {
		setcolor(1);
		drawflow3d(gno, i, ino, curtime);
	    }
	}
	if (g[gno].salip.lactive == ON) {
	    dolegend(g[gno].salip, mapisolconc, g[gno].salip.nisol);
	}
	if (g[gno].velmagip.lactive == ON) {
	    dolegend(g[gno].velmagip, mapisolconc, g[gno].velmagip.nisol);
	}
    }
    setcolor(1);
    setlinewidth(1);

    if (zno < 0) {
	for (i = 0; i < MAXSLICES; i++) {
	    if (display_object(gno, SLICE, i)) {
		drawslice(gno, i, ino, curtime);
	    }
	}
	for (i = 0; i < MAXVELHIST; i++) {
	    if (object_isactive(HISTORY, FLOW, i) && display_object(gno, HISTORY, FLOW, i)) {
		setcolor(g[gno].flowh[i].p.color);
		drawflowh(gno, i, ino);
	    }
	}

	setcolor(1);
	setlinewidth(1);
    }
    draw_init_drogues(gno);	/* placing drogues with the mouse */
    for (i = 0; i < MAXPATHLINES; i++) {
	if (object_isactive(DROGUES, i) && display_object(gno, DROGUES, i)) {
	    drawdrogues(gno, i, ino);
	}
    }
    for (i = 0; i < MAXTRACK; i++) {
	if (object_isactive(TRACK, i) && display_object(gno, TRACK, i)) {
	    drawtrack(gno, i, ino);
	}
    }

    setcolor(1);
    setlinewidth(1);

    if (zno < 0) {
	for (i = 0; i < MAXHISTMARKERS; i++) {
	    if (display_object(gno, HISTORY, 0, i)) {
		drawhist(gno, i, ino, curtime);
	    }
	}
	setcolor(1);
	setlinewidth(1);
	if (display_object(gno, TIDALCLOCK)) {
	    do_drawclock(gno, curtime);
	}
	if (display_object(gno, TIMEINFO)) {
	    drawtimeinfo(gno, curtime);
	}
	if (display_object(gno, MAPSCALE)) {
	    domapscale(gno);
	}
	if (display_object(gno, VSCALE)) {
	    dovscale(gno);
	}
	if (display_object(gno, WSCALE)) {
	    dowscale(gno);
	}
	if (display_object(gno, FLOW)) {
	    dofscale(gno);
	}
	if (display_object(gno, TIMELINE)) {
	    drawtimeline(gno, curtime);
	}
	for (i = 0; i < MAXTRANS; i++) {
	    if (trans[i].display == ON && trans[i].transgno == gno && trans[i].npts) {
		setcolor(1);
		setlinewidth(2);
		draw_transect(trans[i].npts, trans[i].x, trans[i].y);
	    }
	}
	setcolor(1);
	setlinestyle(1);
	setlinewidth(1);

	if (g[gno].pointset) {
	    draw_ref_point(gno);
	}
	if (g[gno].f.active == ON) {
	    boxplot(gno, g[gno].f);
	}
	setcolor(1);
	setlinestyle(1);
	setlinewidth(1);
	if (timestamp.active == ON) {
	    double xtmp, ytmp;
	    set_timestamp();
	    setlinewidth(timestamp.linew);
	    setcolor(timestamp.color);
	    setcharsize(timestamp.charsize);
	    setfont(timestamp.font);
	    view2world(timestamp.x, timestamp.y, &xtmp, &ytmp);
	    writestr(xtmp, ytmp, timestamp.rot, timestamp.just, timestamp.s);
	}
	if (region_flag) {
	    setcolor(1);
	    draw_single_region(regionx, regiony, nregion);
	}
	if (lab.title.s[0]) {
	    setlinewidth(lab.title.linew);
	    setcolor(lab.title.color);
	    setcharsize(lab.title.charsize);
	    setfont(lab.title.font);
	    drawtitle(lab.title.s, 0);
	}
	if (lab.stitle.s[0]) {
	    setlinewidth(lab.stitle.linew);
	    setcolor(lab.stitle.color);
	    setcharsize(lab.stitle.charsize);
	    setfont(lab.stitle.font);
	    drawtitle(lab.stitle.s, 1);
	}
/*
   for (i = 0; i < MAXREGION; i++) {
   if (rg[i].active == ON && rg[i].linkto[gno]) {
   draw_region(i);
   }
   }
 */
    }
    if (!hardcopyflag) {
    } else {
	g[gno].w.xg1 = x1;
	g[gno].w.xg2 = x2;
	g[gno].w.yg1 = y1;
	g[gno].w.yg2 = y2;
    }
}

void boxplot(int gno, framep f)
{
    int c, s, wi;

    c = setcolor(f.color);
    s = setlinestyle(f.lines);
    wi = setlinewidth(f.linew);
    if (f.type == 0) {
	rect(g[gno].w.xg1, g[gno].w.yg1, g[gno].w.xg2, g[gno].w.yg2);
    } else {
	my_move2(g[gno].w.xg1, g[gno].w.yg1);
	my_draw2(g[gno].w.xg2, g[gno].w.yg1);
	my_move2(g[gno].w.xg1, g[gno].w.yg1);
	my_draw2(g[gno].w.xg1, g[gno].w.yg2);
    }
    setcolor(c);
    setlinestyle(s);
    setlinewidth(wi);
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

/*
 * draw annotative text
 */
void draw_string(int gno, int i)
{
    double xtmp1, ytmp1;
    int f, c, w, loctype;
    double s;
    plotstr pstr;

    get_graph_string(i, &pstr);
    if (debuglevel == 5) {
	printf("String %d %s\n", i, pstr.s);
    }
    if (gno != -2) {
	if (pstr.loctype == WORLD && pstr.gno != gno) {
	    return;
	}
	if (pstr.loctype == VIEW && gno != -1) {
	    return;
	}
    }
    switch (pstr.symloc) {
    case RIGHT:
	loctype = 0;
	break;
    case LEFT:
	loctype = 1;
	break;
    case ABOVE:
	loctype = 2;
	break;
    case BELOW:
	loctype = 3;
	break;
    default:
	loctype = 0;
	break;
    }
    if (strlen(pstr.s) && (pstr.charsize > 0.0) && (pstr.active == ON)) {
	c = setcolor(pstr.color);
	w = setlinewidth(pstr.linew);
	s = setcharsize(pstr.charsize);
	f = setfont(pstr.font);
	if (pstr.loctype == WORLD) {
	    if (pstr.sym) {

		writesymstr(pstr.x, pstr.y, pstr.rot, pstr.just, pstr.s, pstr.sym, pstr.symcolor, pstr.symfill, loctype, pstr.symsize);
	    } else {
		writestr(pstr.x, pstr.y, pstr.rot, pstr.just, pstr.s);
	    }
	} else {
	    view2world(pstr.x, pstr.y, &xtmp1, &ytmp1);
	    if (pstr.sym) {

		writesymstr(xtmp1, ytmp1, pstr.rot, pstr.just, pstr.s, pstr.sym, pstr.symcolor, pstr.symfill, loctype, pstr.symsize);
	    } else {
		writestr(xtmp1, ytmp1, pstr.rot, pstr.just, pstr.s);
	    }
	}
	(void) setcolor(c);
	(void) setlinewidth(w);
	(void) setcharsize(s);
	(void) setfont(f);
    }
}

void draw_ref_point(int gno)
{
    drawpolysym(&g[gno].dsx, &g[gno].dsy, 1, SYM_CIRCLE, 0, 0, 1.0);
    drawpolysym(&g[gno].dsx, &g[gno].dsy, 1, SYM_PLUS, 0, 0, 1.0);
    drawpolysym(&g[gno].dsx, &g[gno].dsy, 1, SYM_PLUS, 0, 0, 1.0);
}

void draw_annotation(int gno)
{
    int i;

    setclipping(0);		/* shut down clipping for strings, boxes,
				 * lines, and legends */
    if (debuglevel == 5) {
	printf("Boxes\n");
    }
    for (i = 0; i < MAXBOXES; i++) {
	if (isactive_box(i)) {
	    draw_box(gno, i);
	}
    }
    if (debuglevel == 5) {
	printf("Lines\n");
    }
    for (i = 0; i < MAXLINES; i++) {
	if (isactive_line(i)) {
	    draw_line(gno, i);
	}
    }
    if (debuglevel == 5) {
	printf("Strings\n");
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

void domapscale(int gno)
{
    map_scale m;
    double x, y;
    char buf[50];
    int units;
    double unitfac;
    m = g[gno].mapscale;
    units = m.units;
    unitfac = m.unitfac;
    switch (units) {
    case MM:
	sprintf(buf, "%.1lf mm", m.len);
	break;
    case CM:
	sprintf(buf, "%.1lf cm", m.len);
	break;
    case M:
	sprintf(buf, "%.1lf m", m.len);
	break;
    case KM:
	sprintf(buf, "%.1lf km", m.len);
	break;
    }
    if (m.loctype == WORLD) {
	x = m.x;
	y = m.y;
    } else {
	view2world(m.x, m.y, &x, &y);
    }
    setcolor(m.p.color);
    drawmapscale(x, y, m.len * m.unitfac, buf);
}

void dovscale(int gno)
{
    double x, y;
    velocity_scale v;
    int units;
    double unitfac;
    char buf[50];
    v = g[gno].vl;
    units = v.units;
    unitfac = v.unitfac;
    switch (units) {
    case MM:
	sprintf(buf, "%.1lf mm/s", v.len);
	break;
    case CM:
	sprintf(buf, "%.1lf cm/s", v.len);
	break;
    case M:
	sprintf(buf, "%.1lf m/s", v.len);
	break;
    case KM:
	sprintf(buf, "%.1lf km/s", v.len);
	break;
    }
    if (v.loctype == WORLD) {
	x = v.x;
	y = v.y;
    } else {
	view2world(v.x, v.y, &x, &y);
    }
    setcolor(v.p.color);
    drawvellegend(x, y, v.len, 0.0, unitfac * v.scale, buf);
}

void dowscale(int gno)
{
    double x, y;
    velocity_scale v;
    int units;
    double unitfac;
    char buf[50];
    v = g[gno].wl;
    units = v.units;
    unitfac = v.unitfac;
    switch (units) {
    case MM:
	sprintf(buf, "%.1lf mm/s", v.len);
	break;
    case CM:
	sprintf(buf, "%.1lf cm/s", v.len);
	break;
    case M:
	sprintf(buf, "%.1lf m/s", v.len);
	break;
    case KM:
	sprintf(buf, "%.1lf km/s", v.len);
	break;
    }
    if (v.loctype == WORLD) {
	x = v.x;
	y = v.y;
    } else {
	view2world(v.x, v.y, &x, &y);
    }
    setcolor(v.p.color);
    drawvellegend(x, y, v.len, 0.0, unitfac * v.scale, buf);
}

void dofscale(int gno)
{
    double x, y;
    velocity_scale v;
    int units;
    double unitfac;
    char buf[50];
    v = g[gno].fl;
    units = v.units;
    unitfac = v.unitfac;
    switch (units) {
    case MM:
	sprintf(buf, "%.1lf mm^3/s", v.len);
	break;
    case CM:
	sprintf(buf, "%.1lf cm^3/s", v.len);
	break;
    case M:
	sprintf(buf, "%.1lf m^3/s", v.len);
	break;
    case KM:
	sprintf(buf, "%.1lf km^3/s", v.len);
	break;
    }
    if (v.loctype == WORLD) {
	x = v.x;
	y = v.y;
    } else {
	view2world(v.x, v.y, &x, &y);
    }
    setcolor(v.p.color);
    drawfluxlegend(x, y, 0.0, -v.len, unitfac * v.scale, buf);
}

void dolegend(Isolparms ip, int *cmap, int nmap)
{
    double s, x, y;
    int f;
    switch (ip.p.format) {
    case DECIMAL:
	f = 0;
	break;
    case EXPONENTIAL:
	f = 1;
	break;
    case GENERAL:
	f = 2;
	break;
    }
    if (ip.loctype == WORLD) {
	x = ip.x;
	y = ip.y;
    } else {
	view2world(ip.x, ip.y, &x, &y);
    }
    s = setcharsize(ip.p.charsize);
    if (ip.layout == HORIZONTAL) {
        drawisollegend(1, ip.frame == ON, ip.framecol, x, y, ip.cis, ip.nisol,
                       ip.xlen, ip.ylen, ip.xgap, ip.ygap, ip.p.charsize,
                       f, ip.p.prec, ip.p.color,
                       cmap, nmap, 0);
    } else {
        drawisollegend(0, ip.frame == ON, ip.framecol, x, y, ip.cis, ip.nisol,
                       ip.xlen, ip.ylen, ip.xgap, ip.ygap, ip.p.charsize,
                       f, ip.p.prec, ip.p.color,
                       cmap, nmap, 0);
    }
    (void) setcharsize(s);
}

void drawtimeline(int gno, double time)
{
    double x, y;
    time_line t;
    char buf[50];
    int f;
    t = g[gno].timeline;
    switch (t.p.format) {
    case DECIMAL:
	f = 0;
	break;
    case EXPONENTIAL:
	f = 1;
	break;
    case GENERAL:
	f = 2;
	break;
    }
    if (t.loctype == WORLD) {
	x = t.x;
	y = t.y;
    } else {
	view2world(t.x, t.y, &x, &y);
    }
    setcolor(t.p.color);
    timeline(time, t.start, t.stop, t.step, t.units, t.p.prec, x, y, t.len, 0, t.c1, t.c2, t.c3);
}

void drawtimeinfo(int gno, double runtime)
{
    double x, y;
    time_info t;
    struct tm *tmp;
    time_t ts;
    int len;
    buf[256];
    t = g[gno].timeinfo;
    if (t.loctype == WORLD) {
	x = t.x;
	y = t.y;
    } else {
	view2world(t.x, t.y, &x, &y);
    }
    if (t.display == 0) {
	ts = (time_t) t.time + (time_t) runtime;
	tmp = gmtime(&ts);
	len = strftime(buf, 255, t.format, tmp);
	if (len == 0) {
	    strcpy(t.s, ctime(&ts));
	} else {
	    strcpy(t.s, buf);
	}
    } else {
	int d, h, m, s;
	double dd, hh, mm, ss;
	dd = runtime / 86400.0;	/* convert to days */
	d = (int) dd;		/* just the days */
	hh = (dd - d) * 24.0;	/* convert fractional days to hours */
	h = (int) hh;		/* just the hours */
	mm = (hh - h) * 60.0;	/* convert fractional hours to minutes */
	m = (int) mm;		/* Just the minutes */
        ss = (mm - m) * 60.0;
        s = (int) ss;
	sprintf(t.s, "Time = %02d:%02d:%02d:%02d (DD:HH:MM:SS)", d, h, m, s);
    }
    if (strlen(t.s) > 0) {
	int c, w, f;
	double s;
	c = setcolor(t.color);
	w = setlinewidth(t.linew);
	s = setcharsize(t.charsize);
	f = setfont(t.font);
	writestrbox(x, y, t.rot, t.just, t.s);
	setcolor(c);
	setlinewidth(w);
	setcharsize(s);
	setfont(f);
    }
}

void do_drawclock(int gno, double time)
{
    double x, y;
    tidal_clock t;
    t = g[gno].tidalclock;
    if (t.loctype == WORLD) {
	x = t.x;
	y = t.y;
    } else {
	view2world(t.x, t.y, &x, &y);
    }
    drawclock(x, y, t.total_time, time, t.p.color, t.p.fillcol);
}

void draw_single_region(double *rx, double *ry, int nr)
{
    int i;

    if (nr >= 3) {
	setlinewidth(3);
	my_move2(rx[0], ry[0]);
	for (i = 0; i < nr; i++) {
	    my_draw2(rx[i], ry[i]);
	}
	my_draw2(rx[0], ry[0]);
	setlinewidth(1);
    }
}

void draw_transect(int n, double *x, double *y)
{
    int i;

    if (n > 2) {
	my_move2(x[0], y[0]);
	for (i = 0; i < n; i++) {
	    my_draw2(x[i], y[i]);
	}
	drawpolysym(x, y, n, SYM_CIRCLE, 0, 0, 0.6);
    }
}
