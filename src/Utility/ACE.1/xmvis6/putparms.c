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
 * write a parameter file
 *
 */

#ifndef lint
static char RCSid[] =
    "$Id: putparms.c,v 1.7 2007/09/06 22:47:46 pturner Exp $";

#endif

#include "defines.h"
#include "globals.h"

#define MAXAXES 6

extern char *get_format_types(int format);

/*
 * Function prototypes
 */
void put_annotation(int gno, FILE * pp, int imbed);
void put_region(int gno, FILE * pp, int imbed);
void WriteCMap(FILE * fp);
static int write_DisplayBoundary(DisplayBoundary * d, FILE * fout);
static int write_DisplayFlow(DisplayFlow d, int flowtype, int flowno,
			     FILE * pp);
static int write_Display3dFlow(int n, Display3dFlow d, FILE * pp);
static int write_DisplayScalar3d(DisplayScalar3d * d, FILE * fout);
static int write_DisplayScalar2d(DisplayScalar2d * d, FILE * fout);
static int write_DisplayParticles(DisplayParticles * d, FILE * fout);
static int write_Velocity_marker(Velocity_marker * d, FILE * fout);
static int write_Conc_marker(Conc_marker * d, FILE * fout);
static int write_Hist_marker(Hist_marker d, int hno, char *prefix,
			     FILE * fout);
static int write_north_indicator(north_indicator * d, FILE * fout);
static int write_DisplaySlice(DisplaySlice d, int sno, char *prefix,
			      FILE * pp);
static int write_DisplayFlux(DisplayFlux * d, FILE * fout);
static int write_graph(int gno, FILE * pp, int imbed);

static int write_Zoom_box(Zoom_box d, int zno, char *prefix, FILE * pp)
{
    fprintf(pp, "%s with zoombox %d\n", prefix, zno);
    fprintf(pp, "%s zoombox %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s zoombox display %s\n", prefix, on_or_off(d.display));
    fprintf(pp, "%s zoombox display marker %s\n", prefix,
	    on_or_off(d.display_marker));
    fprintf(pp, "%s zoombox color %d\n", prefix, d.p.color);
    fprintf(pp, "%s zoombox fill color %d\n", prefix, d.p.fillcol);
    fprintf(pp, "%s zoombox linewidth %d\n", prefix, d.p.linew);
    fprintf(pp, "%s zoombox loctype %s\n", prefix, w_or_v(d.loctype));
    fprintf(pp, "%s zoombox loc %lf, %lf\n", prefix, d.locx, d.locy);
    fprintf(pp, "%s zoombox xy %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s zoombox world %lf, %lf, %lf, %lf\n", prefix, d.wx1,
	    d.wy1, d.wx2, d.wy2);
    fprintf(pp, "%s zoombox view %lf, %lf\n", prefix, d.vx, d.vy);
    fprintf(pp, "%s zoombox zoom %lf\n", prefix, d.expand);
    fprintf(pp, "%s zoombox attach %d\n", prefix, d.attach);
    fprintf(pp, "%s zoombox prec %d, %d\n", prefix, d.precx, d.precy);
}

static int write_time_info(time_info s, char *prefix, FILE * pp)
{
    fprintf(pp, "%s timeinfo %s\n", prefix, on_or_off(s.active));
    fprintf(pp, "%s timeinfo loctype %s\n", prefix, w_or_v(s.loctype));
    fprintf(pp, "%s timeinfo %.12lg, %.12lg\n", prefix, s.x, s.y);
    fprintf(pp, "%s timeinfo linewidth %d\n", prefix, s.linew);
    fprintf(pp, "%s timeinfo color %d\n", prefix, s.color);
    fprintf(pp, "%s timeinfo rot %d\n", prefix, s.rot);
    fprintf(pp, "%s timeinfo font %d\n", prefix, s.font);
    fprintf(pp, "%s timeinfo just %d\n", prefix, s.just);
    fprintf(pp, "%s timeinfo char size %lf\n", prefix, s.charsize);
    fprintf(pp, "%s timeinfo format \"%s\"\n", prefix, s.format);
    fprintf(pp, "%s timeinfo start \"%s\"\n", prefix, s.start);
}

static int write_tidal_clock(tidal_clock d, char *prefix, FILE * pp)
{
    fprintf(pp, "%s tidalclock %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s tidalclock color %d\n", prefix, d.p.color);
    fprintf(pp, "%s tidalclock fill color %d\n", prefix, d.p.fillcol);
    fprintf(pp, "%s tidalclock loctype %s\n", prefix, w_or_v(d.loctype));
    fprintf(pp, "%s tidalclock loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s tidalclock total time %lf\n", prefix, d.total_time);
}

static int write_mapscale(map_scale d, char *prefix, FILE * pp)
{
    char s[256];
    fprintf(pp, "%s mapscale %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s mapscale color %d\n", prefix, d.p.color);
    fprintf(pp, "%s mapscale scale %lf\n", prefix, d.scale);
    fprintf(pp, "%s mapscale length %lf\n", prefix, d.len);
    fprintf(pp, "%s mapscale loc %lf, %lf\n", prefix, d.x, d.y);
    switch (d.units) {
    case MM:
	strcpy(s, "mm");
	break;
    case CM:
	strcpy(s, "cm");
	break;
    case M:
	strcpy(s, "m");
	break;
    case KM:
	strcpy(s, "km");
	break;
    }
    fprintf(pp, "%s mapscale units %s\n", prefix, s);
}

static int write_velocity_scale(velocity_scale d, char *prefix, FILE * pp)
{
    char s[256];
    fprintf(pp, "%s vscale %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s vscale color %d\n", prefix, d.p.color);
    fprintf(pp, "%s vscale scale %lf\n", prefix, d.scale);
    fprintf(pp, "%s vscale length %lf\n", prefix, d.len);
    fprintf(pp, "%s vscale loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s vscale loctype %s\n", prefix, w_or_v(d.loctype));
    switch (d.units) {
    case MM:
	strcpy(s, "mm");
	break;
    case CM:
	strcpy(s, "cm");
	break;
    case M:
	strcpy(s, "m");
	break;
    case KM:
	strcpy(s, "km");
	break;
    }
    fprintf(pp, "%s vscale units %s\n", prefix, s);
}

static int write_wind_scale(velocity_scale d, char *prefix, FILE * pp)
{
    char s[256];
    fprintf(pp, "%s wscale %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s wscale color %d\n", prefix, d.p.color);
    fprintf(pp, "%s wscale scale %lf\n", prefix, d.scale);
    fprintf(pp, "%s wscale length %lf\n", prefix, d.len);
    fprintf(pp, "%s wscale loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s wscale loctype %s\n", prefix, w_or_v(d.loctype));
    switch (d.units) {
    case MM:
	strcpy(s, "mm");
	break;
    case CM:
	strcpy(s, "cm");
	break;
    case M:
	strcpy(s, "m");
	break;
    case KM:
	strcpy(s, "km");
	break;
    }
    fprintf(pp, "%s wscale units %s\n", prefix, s);
}

static int write_flux_scale(flux_scale d, char *prefix, FILE * pp)
{
    char s[256];
    fprintf(pp, "%s flux scale %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s flux scale color %d\n", prefix, d.p.color);
    fprintf(pp, "%s flux scale scale %lf\n", prefix, d.scale);
    fprintf(pp, "%s flux scale length %lf\n", prefix, d.len);
    fprintf(pp, "%s flux scale loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s flux scale loctype %s\n", prefix, w_or_v(d.loctype));
    switch (d.units) {
    case MM:
	strcpy(s, "mm");
	break;
    case CM:
	strcpy(s, "cm");
	break;
    case M:
	strcpy(s, "m");
	break;
    case KM:
	strcpy(s, "km");
	break;
    }
    fprintf(pp, "%s flux scale units %s\n", prefix, s);
}

static int write_time_line(time_line d, char *prefix, FILE * pp)
{
    fprintf(pp, "%s timeline %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s timeline prec %d\n", prefix, d.p.prec);
    fprintf(pp, "%s timeline color %d\n", prefix, d.c1);
    fprintf(pp, "%s timeline fill color %d\n", prefix, d.c2);
    fprintf(pp, "%s timeline units %d\n", prefix, d.units);
    fprintf(pp, "%s timeline start %lf\n", prefix, d.start);
    fprintf(pp, "%s timeline stop %lf\n", prefix, d.stop);
    fprintf(pp, "%s timeline step %lf\n", prefix, d.step);
    fprintf(pp, "%s timeline width %d\n", prefix, d.width);
    fprintf(pp, "%s timeline length %d\n", prefix, d.len);
    fprintf(pp, "%s timeline loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s timeline loctype %s\n", prefix, w_or_v(d.loctype));
}

static int write_Elevmarker(int elno, Elevmarker d, char *ftype, FILE * pp)
{
    fprintf(pp, "%s elevmarker %d\n", ftype, elno);
    fprintf(pp, "%s elevmarker active %s\n", ftype, on_or_off(d.active));
    fprintf(pp, "%s elevmarker display %s\n", ftype, on_or_off(d.display));
    fprintf(pp, "%s elevmarker node %d\n", ftype, d.node);
    fprintf(pp, "%s elevmarker xy %lf, %lf\n", ftype, d.x, d.y);
    fprintf(pp, "%s elevmarker loctype %s\n", ftype, w_or_v(d.loctype));
    fprintf(pp, "%s elevmarker loc %lf, %lf\n", ftype, d.locx, d.locy);
    fprintf(pp, "%s elevmarker minmax %lf, %lf\n", ftype, d.emin, d.emax);
}

void putparms(int gno, char *fname)
{
    int i, j, k, cnt, ming, maxg;
    FILE *pp;
    if ((pp = fopen(fname, "w")) == NULL) {
	sprintf(buf, "Can't open parameter file %s", fname);
	errwin(buf);
    }
    write_graph(gno, pp, 0);
    WriteCMap(pp);
    fclose(pp);
}

char *getfilltype(int fillusing)
{
    switch (fillusing) {
    case COLOR:
	return "color";
	break;
    case PATTERN:
	return "pattern";
	break;
    case NONE:
	return "none";
	break;
    }
    return "none";
}

/*
 * Write type Props version 0
 */
static int write_Props(Props d, char *prefix, FILE * fout)
{
    fprintf(fout, "%s prop color %d\n", prefix, d.color);
    fprintf(fout, "%s prop linewidth %d\n", prefix, d.linew);
    fprintf(fout, "%s prop linestyle %d\n", prefix, d.lines);
    fprintf(fout, "%s prop font %d\n", prefix, d.font);
    fprintf(fout, "%s prop format %s\n", prefix,
	    get_format_types(d.format));
    fprintf(fout, "%s prop prec %d\n", prefix, d.prec);
    /*fprintf(fout, "%s prop points %d\n", prefix, d.points); */
    fprintf(fout, "%s prop char size %lf\n", prefix, d.charsize);
    fprintf(fout, "%s prop symbol %d\n", prefix, d.symbol);
    fprintf(fout, "%s prop symbol size %d\n", prefix, d.symsize);
    fprintf(fout, "%s prop fill %s\n", prefix, on_or_off(d.fill));
    fprintf(fout, "%s prop fill %s\n", prefix, getfilltype(d.fillusing));
    fprintf(fout, "%s prop fill color %d\n", prefix, d.fillcol);
    fprintf(fout, "%s prop fill pattern %d\n", prefix, d.fillpat);
    fprintf(fout, "%s prop arrow %d\n", prefix, d.arrow);
    fprintf(fout, "%s prop arrow type %d\n", prefix, d.atype);
    fprintf(fout, "%s prop arrow size %lf\n", prefix, d.asize);
}

/* write isoline settings */
static int write_Isolparms(Isolparms d, char *prefix, FILE * pp)
{
    int i;
    char s[256];
    fprintf(pp, "%s isolines %d\n", prefix, d.nisol);
    fprintf(pp, "%s isolines type %d\n", prefix, d.type);
    fprintf(pp, "%s isolines set type %d\n", prefix, d.isoltype);
    fprintf(pp, "%s isolines legend %s\n", prefix, on_or_off(d.lactive));
    fprintf(pp, "%s isolines legend loctype %s\n", prefix,
	    w_or_v(d.loctype));
    fprintf(pp, "%s isolines legend layout %s\n", prefix,
	    d.layout == HORIZONTAL ? "horizontal" : "vertical");
    fprintf(pp, "%s isolines legend %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s isolines legend size %lf, %lf\n", prefix, d.xlen,
	    d.ylen);
    fprintf(pp, "%s isolines legend hgap %lf, %lf\n", prefix, d.xgap,
	    d.ygap);
    fprintf(pp, "%s isolines start %lf step %lf\n", prefix, d.cis[0],
	    d.cint);
    sprintf(s, "%s isolines", prefix);
    write_Props(d.p, s, pp);
    for (i = 0; i < d.nisol; i++) {
	fprintf(pp, "%s isoline %d, %lf\n", prefix, i, d.cis[i]);
    }
    for (i = 0; i < d.nisol; i++) {
	fprintf(pp, "%s isoline %d color %d\n", prefix, i, d.color[i]);
	fprintf(pp, "%s isoline %d linewidth %d\n", prefix, i, d.linew[i]);
	fprintf(pp, "%s isoline %d linestyle %d\n", prefix, i, d.lines[i]);
    }
}

static int write_DisplayTideStation(DisplayTideStation * d, FILE * fout)
{
}

static int write_Hist_marker(Hist_marker d, int hno, char *prefix,
			     FILE * pp)
{
    int i;
    fprintf(pp, "%s with histbox %d\n", prefix, hno);
    fprintf(pp, "%s histbox display %s\n", prefix, on_or_off(d.display));
    fprintf(pp, "%s histbox display marker %s\n", prefix,
	    on_or_off(d.display_marker));
    fprintf(pp, "%s histbox color %d\n", prefix, d.p.color);
    fprintf(pp, "%s histbox fill color %d\n", prefix, d.p.fillcol);
    fprintf(pp, "%s histbox linewidth %d\n", prefix, d.p.linew);
    fprintf(pp, "%s histbox loctype %s\n", prefix, w_or_v(d.loctype));
    fprintf(pp, "%s histbox loc %lf, %lf\n", prefix, d.locx, d.locy);
    fprintf(pp, "%s histbox xy %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s histbox tick %lf, %lf\n", prefix, d.xtickm, d.ytickm);
    fprintf(pp, "%s histbox world %lf, %lf, %lf, %lf\n", prefix, d.wx1,
	    d.wy1, d.wx2, d.wy2);
    fprintf(pp, "%s histbox view %lf, %lf\n", prefix, d.vx, d.vy);
    fprintf(pp, "%s histbox attach %d\n", prefix, d.attach);
    fprintf(pp, "%s histbox prec %d, %d\n", prefix, d.precx, d.precy);
    if (d.thist) {
	fprintf(pp, "%s histbox display history true\n", prefix);
	fprintf(pp, "%s histbox display history color %d\n", prefix,
		d.ap[i].color);
    }
    for (i = 0; i < MAXADCIRC; i++) {
	if (d.adcirc[i]) {
	    fprintf(pp, "%s histbox display adcirc %d true\n", prefix, i);
	    fprintf(pp, "%s histbox display adcirc %d color %d\n", prefix,
		    i, d.ap[i].color);
	}
    }
}

static int write_Display3dFlow(int n, Display3dFlow d, FILE * pp)
{
    fprintf(pp, "with elcirc marker %d\n", n);
    fprintf(pp, "elcirc marker display %s\n", on_or_off(d.display));
    fprintf(pp, "elcirc marker display marker %s\n",
	    on_or_off(d.display_marker));
    fprintf(pp, "elcirc marker color %d\n", d.p.color);
    fprintf(pp, "elcirc marker fill color %d\n", d.p.fillcol);
    fprintf(pp, "elcirc marker linewidth %d\n", d.p.linew);
    fprintf(pp, "elcirc marker loctype %s\n", w_or_v(d.loctype));
    fprintf(pp, "elcirc marker loc %lf, %lf\n", d.locx, d.locy);
    fprintf(pp, "elcirc marker xy %lf, %lf\n", d.x, d.y);
    fprintf(pp, "elcirc marker world %lf, %lf, %lf, %lf\n", d.wx1, d.wy1,
	    d.wx2, d.wy2);
    fprintf(pp, "elcirc marker view %lf, %lf\n", d.vx, d.vy);
    fprintf(pp, "elcirc marker attach %d\n", d.attach);
    fprintf(pp, "elcirc marker prec %d, %d\n", d.precx, d.precy);
}

static int write_DisplayScalar3d(DisplayScalar3d * d, FILE * fout)
{
}
static int write_DisplayScalar2d(DisplayScalar2d * d, FILE * fout)
{
}
static int write_DisplayParticles(DisplayParticles * d, FILE * fout)
{
}
static int write_Velocity_marker(Velocity_marker * d, FILE * fout)
{
}
static int write_Conc_marker(Conc_marker * d, FILE * fout)
{
}

static int write_north_indicator(north_indicator * d, FILE * fout)
{
}

static int write_DisplaySlice(DisplaySlice d, int sno, char *prefix,
			      FILE * pp)
{
    fprintf(pp, "%s with slice %d\n", prefix, sno);
    fprintf(pp, "%s slice %s\n", prefix, on_or_off(d.active));
    fprintf(pp, "%s slice loc %lf, %lf\n", prefix, d.x, d.y);
    fprintf(pp, "%s slice loctype %s\n", prefix, w_or_v(d.loctype));
    fprintf(pp, "%s slice world %lf, %lf, %lf, %lf\n", d.wx1, d.wy1, d.wx2,
	    d.wy2);
    fprintf(pp, "%s slice color %d\n", d.p.color);
    fprintf(pp, "%s slice fill color %d\n", d.p.fillcol);
    fprintf(pp, "%s slice linewidth %d\n", d.p.linew);
    fprintf(pp, "%s slice world %lf, %lf, %lf, %lf\n", d.wx1, d.wy1, d.wx2,
	    d.wy2);
    fprintf(pp, "%s slice view %lf, %lf\n", d.vx, d.vy);
    fprintf(pp, "%s slice attach %d\n", d.attach);
    fprintf(pp, "%s slice prec %d, %d\n", d.precx, d.precy);
}

static int write_DisplayTransect(DisplayTransect d, int sno, char *prefix,
				 FILE * pp)
{
    fprintf(pp, "%s set elcirc transect %d\n", prefix, sno);
    write_Isolparms(d.ip, "elcirc transect", pp);
}

static int write_DisplayFlux(DisplayFlux * d, FILE * fout)
{
}

static int write_DisplayGrid(DisplayGrid d, int gridno, FILE * pp)
{
    fprintf(pp, "with grid %d\n", gridno);
    fprintf(pp, "grid display %s\n", on_or_off(d.display));
    fprintf(pp, "grid boundary display %s\n",
	    on_or_off(d.display_boundary));
    fprintf(pp, "grid nodes display %s\n", on_or_off(d.display_nodes));
    fprintf(pp, "grid elements display %s\n",
	    on_or_off(d.display_elements));
    fprintf(pp, "grid depth display %s\n", on_or_off(d.display_depths));
    fprintf(pp, "grid fill %s\n", on_or_off(d.display_gridf));
    write_Props(d.p, "grid", pp);
    write_Props(d.bp, "grid boundary", pp);
    fprintf(pp, "grid bath display %s\n", on_or_off(d.display_bath));
    if (d.display_bath == ON) {
	write_Isolparms(d.ip, "grid bath", pp);
    }
}

static int write_graph(int gno, FILE * pp, int imbed)
{
    int i, j, k, kk, ming, maxg;
    int ps, pt, gh, gl, gt, fx, fy, px, py;
    double dsx, dsy;
    char imbedstr[2], tmpstr1[128], tmpstr2[128];
    framep f;
    labels lab;
    tickmarks t;
    world w;
    view v;

    if (imbed) {
	strcpy(imbedstr, "@");
    } else {
	imbedstr[0] = 0;
    }
    fprintf(pp, "# ACE/vis parameter file\n");
    fprintf(pp, "#\n");
    fprintf(pp, "#\n");
    fprintf(pp, "%spage %d\n", imbedstr, (int) (scrollper * 100));
    fprintf(pp, "%spage inout %d\n", imbedstr, (int) (shexper * 100));
    fprintf(pp, "%slink page %s\n", imbedstr,
	    scrolling_islinked ? "on" : "off");
    put_annotation(gno, pp, imbed);
    put_region(gno, pp, imbed);
    if (gno == -1) {
	maxg = maxgraph - 1;
	ming = 0;
    } else {
	maxg = gno;
	ming = gno;
    }
    for (k = ming; k <= maxg; k++) {
	if (isactive_graph(k)) {
	    gno = k;
	    gh = g[gno].hidden;
	    gl = g[gno].label;
	    gt = g[gno].type;
	    ps = g[gno].pointset;
	    pt = g[gno].pt_type;
	    dsx = g[gno].dsx;
	    dsy = g[gno].dsy;
	    fx = g[gno].fx;
	    fy = g[gno].fy;
	    px = g[gno].px;
	    py = g[gno].py;

	    fprintf(pp, "%swith g%1d\n", imbedstr, gno);

	    fprintf(pp, "%sg%1d %s\n", imbedstr, gno,
		    on_or_off(g[gno].active));
	    fprintf(pp, "%sg%1d label %s\n", imbedstr, gno, on_or_off(gl));
	    fprintf(pp, "%sg%1d hidden %s\n", imbedstr, gno,
		    gh ? "true" : "false");
	    fprintf(pp, "%sg%1d type %s\n", imbedstr, gno,
		    graph_types(g[gno].type, 1));
	    fprintf(pp, "%sg%1d autoscale type %s\n", imbedstr, gno,
		    g[gno].auto_type == AUTO ? "AUTO" : "SPEC");
	    fprintf(pp, "%sg%1d fixedpoint %s\n", imbedstr, gno,
		    on_or_off(ps));
	    fprintf(pp, "%sg%1d fixedpoint type %d\n", imbedstr, gno, pt);
	    fprintf(pp, "%sg%1d fixedpoint xy %lf, %lf\n", imbedstr, gno,
		    dsx, dsy);
	    strcpy(tmpstr1, get_format_types(fx));
	    strcpy(tmpstr2, get_format_types(fy));
	    fprintf(pp, "%sg%1d fixedpoint format %s %s\n", imbedstr, gno,
		    tmpstr1, tmpstr2);
	    fprintf(pp, "%sg%1d fixedpoint prec %d, %d\n", imbedstr, gno,
		    px, py);

	    get_graph_world(gno, &w);
	    fprintf(pp, "%s    world xmin %.12lg\n", imbedstr, w.xg1);
	    fprintf(pp, "%s    world xmax %.12lg\n", imbedstr, w.xg2);
	    fprintf(pp, "%s    world ymin %.12lg\n", imbedstr, w.yg1);
	    fprintf(pp, "%s    world ymax %.12lg\n", imbedstr, w.yg2);

	    for (i = 0; i < g[gno].ws_top; i++) {
		fprintf(pp,
			"%s    stack world %.9lg, %.9lg, %.9lg, %.9lg tick %lg, %lg, %lg, %lg\n",
			imbedstr, g[gno].ws[i].w.xg1, g[gno].ws[i].w.xg2,
			g[gno].ws[i].w.yg1, g[gno].ws[i].w.yg2,
			g[gno].ws[i].t[0].xg1, g[gno].ws[i].t[0].xg2,
			g[gno].ws[i].t[0].yg1, g[gno].ws[i].t[0].yg2);
	    }

	    get_graph_view(gno, &v);
	    fprintf(pp, "%s    view xmin %lf\n", imbedstr, v.xv1);
	    fprintf(pp, "%s    view xmax %lf\n", imbedstr, v.xv2);
	    fprintf(pp, "%s    view ymin %lf\n", imbedstr, v.yv1);
	    fprintf(pp, "%s    view ymax %lf\n", imbedstr, v.yv2);

	    get_graph_labels(gno, &lab);
	    fprintf(pp, "%s    title \"%s\"\n", imbedstr, lab.title.s);
	    fprintf(pp, "%s    title font %d\n", imbedstr, lab.title.font);
	    fprintf(pp, "%s    title size %lf\n", imbedstr,
		    lab.title.charsize);
	    fprintf(pp, "%s    title color %d\n", imbedstr,
		    lab.title.color);
	    fprintf(pp, "%s    title linewidth %d\n", imbedstr,
		    lab.title.linew);
	    fprintf(pp, "%s    subtitle \"%s\"\n", imbedstr, lab.stitle.s);
	    fprintf(pp, "%s    subtitle font %d\n", imbedstr,
		    lab.stitle.font);
	    fprintf(pp, "%s    subtitle size %lf\n", imbedstr,
		    lab.stitle.charsize);
	    fprintf(pp, "%s    subtitle color %d\n", imbedstr,
		    lab.stitle.color);
	    fprintf(pp, "%s    subtitle linewidth %d\n", imbedstr,
		    lab.title.linew);

	    for (i = 0; i < MAXAXES; i++) {
		switch (i) {
		case 0:
		    get_graph_tickmarks(gno, &t, X_AXIS);
		    if (t.active == OFF) {
			fprintf(pp, "%s    xaxis off\n", imbedstr);
			continue;
		    }
		    sprintf(buf, "%s    xaxis ", imbedstr);
		    break;
		case 1:
		    get_graph_tickmarks(gno, &t, Y_AXIS);
		    if (t.active == OFF) {
			fprintf(pp, "%s    yaxis off\n", imbedstr);
			continue;
		    }
		    sprintf(buf, "%s    yaxis ", imbedstr);
		    break;
		case 2:
		    get_graph_tickmarks(gno, &t, ZX_AXIS);
		    if (t.active == OFF) {
			fprintf(pp, "%s    zeroxaxis off\n", imbedstr);
			continue;
		    }
		    sprintf(buf, "%s    zeroxaxis ", imbedstr);
		    break;
		case 3:
		    get_graph_tickmarks(gno, &t, ZY_AXIS);
		    if (t.active == OFF) {
			fprintf(pp, "%s    zeroyaxis off\n", imbedstr);
			continue;
		    }
		    sprintf(buf, "%s    zeroyaxis ", imbedstr);
		    break;
		}

		fprintf(pp, "%s tick %s\n", buf, on_or_off(t.active));
		fprintf(pp, "%s tick major %.12lg\n", buf, t.tmajor);
		fprintf(pp, "%s tick minor %.12lg\n", buf, t.tminor);
		fprintf(pp, "%s tick offsetx %lf\n", buf, t.offsx);
		fprintf(pp, "%s tick offsety %lf\n", buf, t.offsy);
/* DEFUNCT
		fprintf(pp, "%s tick alt %s\n", buf, on_or_off(t.alt));
		fprintf(pp, "%s tick min %.12lg\n", buf, t.tmin);
		fprintf(pp, "%s tick max %.12lg\n", buf, t.tmax);
*/

		fprintf(pp, "%s label \"%s\"\n", buf, t.label.s);
		if (t.label_layout == PERP) {
		    fprintf(pp, "%s label layout perp\n", buf);
		} else {
		    fprintf(pp, "%s label layout para\n", buf);
		}
		if (t.label_place == AUTO) {
		    fprintf(pp, "%s label place auto\n", buf);
		} else {
		    fprintf(pp, "%s label place spec\n", buf);
		}
		fprintf(pp, "%s label char size %lf\n", buf,
			t.label.charsize);
		fprintf(pp, "%s label font %d\n", buf, t.label.font);
		fprintf(pp, "%s label color %d\n", buf, t.label.color);
		fprintf(pp, "%s label linewidth %d\n", buf, t.label.linew);

		fprintf(pp, "%s ticklabel %s\n", buf,
			on_or_off(t.tl_flag));
		if (t.tl_type == AUTO) {
		    fprintf(pp, "%s ticklabel type auto\n", buf);
		} else {
		    fprintf(pp, "%s ticklabel type spec\n", buf);
		}
		fprintf(pp, "%s ticklabel prec %d\n", buf, t.tl_prec);
		fprintf(pp, "%s ticklabel format %s\n", buf,
			get_format_types(t.tl_format));
		fprintf(pp, "%s ticklabel append \"%s\"\n", buf,
			t.tl_appstr);
		fprintf(pp, "%s ticklabel prepend \"%s\"\n", buf,
			t.tl_prestr);
		switch (t.tl_layout) {
		case HORIZONTAL:
		    fprintf(pp, "%s ticklabel layout horizontal\n", buf);
		    break;
		case VERTICAL:
		    fprintf(pp, "%s ticklabel layout vertical\n", buf);
		    break;
		case SPEC:
		    fprintf(pp, "%s ticklabel layout spec\n", buf);
		    fprintf(pp, "%s ticklabel angle %d\n", buf,
			    t.tl_angle);
		    break;
		}
		fprintf(pp, "%s ticklabel skip %d\n", buf, t.tl_skip);
		fprintf(pp, "%s ticklabel stagger %d\n", buf,
			t.tl_staggered);
		switch (t.tl_op) {
		case TOP:
		    fprintf(pp, "%s ticklabel op top\n", buf);
		    break;
		case BOTTOM:
		    fprintf(pp, "%s ticklabel op bottom\n", buf);
		    break;
		case LEFT:
		    fprintf(pp, "%s ticklabel op left\n", buf);
		    break;
		case RIGHT:
		    fprintf(pp, "%s ticklabel op right\n", buf);
		    break;
		case BOTH:
		    fprintf(pp, "%s ticklabel op both\n", buf);
		    break;
		}
		switch (t.tl_sign) {
		case NORMAL:
		    fprintf(pp, "%s ticklabel sign normal\n", buf);
		    break;
		case ABSOLUTE:
		    fprintf(pp, "%s ticklabel sign absolute\n", buf);
		    break;
		case NEGATE:
		    fprintf(pp, "%s ticklabel sign negate\n", buf);
		    break;
		}
		fprintf(pp, "%s ticklabel start type %s\n", buf,
			t.tl_starttype == AUTO ? "auto" : "spec");
		fprintf(pp, "%s ticklabel start %lf\n", buf, t.tl_start);
		fprintf(pp, "%s ticklabel stop type %s\n", buf,
			t.tl_stoptype == AUTO ? "auto" : "spec");
		fprintf(pp, "%s ticklabel stop %lf\n", buf, t.tl_stop);
		fprintf(pp, "%s ticklabel char size %lf\n", buf,
			t.tl_charsize);
		fprintf(pp, "%s ticklabel font %d\n", buf, t.tl_font);
		fprintf(pp, "%s ticklabel color %d\n", buf, t.tl_color);
		fprintf(pp, "%s ticklabel linewidth %d\n", buf,
			t.tl_linew);

		fprintf(pp, "%s tick major %s\n", buf,
			on_or_off(t.t_flag));
		fprintf(pp, "%s tick minor %s\n", buf,
			on_or_off(t.t_mflag));
		fprintf(pp, "%s tick default %d\n", buf, t.t_num);
		switch (t.t_inout) {
		case IN:
		    fprintf(pp, "%s tick in\n", buf);
		    break;
		case OUT:
		    fprintf(pp, "%s tick out\n", buf);
		    break;
		case BOTH:
		    fprintf(pp, "%s tick both\n", buf);
		    break;
		}
		fprintf(pp, "%s tick major color %d\n", buf, t.t_color);
		fprintf(pp, "%s tick major linewidth %d\n", buf,
			t.t_linew);
		fprintf(pp, "%s tick major linestyle %d\n", buf,
			t.t_lines);
		fprintf(pp, "%s tick minor color %d\n", buf, t.t_mcolor);
		fprintf(pp, "%s tick minor linewidth %d\n", buf,
			t.t_mlinew);
		fprintf(pp, "%s tick minor linestyle %d\n", buf,
			t.t_mlines);
		fprintf(pp, "%s tick log %s\n", buf, on_or_off(t.t_log));
		fprintf(pp, "%s tick size %lf\n", buf, t.t_size);
		fprintf(pp, "%s tick minor size %lf\n", buf, t.t_msize);
		fprintf(pp, "%s bar %s\n", buf, on_or_off(t.t_drawbar));
		fprintf(pp, "%s bar color %d\n", buf, t.t_drawbarcolor);
		fprintf(pp, "%s bar linestyle %d\n", buf,
			t.t_drawbarlines);
		fprintf(pp, "%s bar linewidth %d\n", buf,
			t.t_drawbarlinew);
		fprintf(pp, "%s tick major grid %s\n", buf,
			on_or_off(t.t_gridflag));
		fprintf(pp, "%s tick minor grid %s\n", buf,
			on_or_off(t.t_mgridflag));
		switch (t.t_op) {
		case TOP:
		    fprintf(pp, "%s tick op top\n", buf);
		    break;
		case BOTTOM:
		    fprintf(pp, "%s tick op bottom\n", buf);
		    break;
		case LEFT:
		    fprintf(pp, "%s tick op left\n", buf);
		    break;
		case RIGHT:
		    fprintf(pp, "%s tick op right\n", buf);
		    break;
		case BOTH:
		    fprintf(pp, "%s tick op both\n", buf);
		    break;
		}
		if (t.t_type == AUTO) {
		    fprintf(pp, "%s tick type auto\n", buf);
		} else {
		    fprintf(pp, "%s tick type spec\n", buf);
		}
		fprintf(pp, "%s tick spec %d\n", buf, t.t_spec);
		for (j = 0; j < t.t_spec; j++) {
		    fprintf(pp, "%s tick %d, %lg\n", buf, j,
			    t.t_specloc[j]);
		    fprintf(pp, "%s ticklabel %d, \"%s\"\n", buf, j,
			    t.t_speclab[j].s);
		}
	    }

	    get_graph_framep(gno, &f);
	    fprintf(pp, "%s    frame %s\n", imbedstr, on_or_off(f.active));
	    fprintf(pp, "%s    frame type %d\n", imbedstr, f.type);
	    fprintf(pp, "%s    frame linestyle %d\n", imbedstr, f.lines);
	    fprintf(pp, "%s    frame linewidth %d\n", imbedstr, f.linew);
	    fprintf(pp, "%s    frame color %d\n", imbedstr, f.color);
	    fprintf(pp, "%s    frame fill %s\n", imbedstr,
		    on_or_off(f.fillbg));
	    fprintf(pp, "%s    frame background color %d\n", imbedstr,
		    f.bgcolor);
	    if (g[gno].mapscale.active == ON) {
		write_mapscale(g[gno].mapscale, "    ", pp);
	    } else {
		fprintf(pp, "%s    mapscale off\n", imbedstr);
	    }
	    if (g[gno].tidalclock.active == ON) {
		write_tidal_clock(g[gno].tidalclock, "    ", pp);
	    } else {
		fprintf(pp, "%s    tidalclock off\n", imbedstr);
	    }
	    if (g[gno].timeline.active == ON) {
		write_time_line(g[gno].timeline, "    ", pp);
	    } else {
		fprintf(pp, "%s    timeline off\n", imbedstr);
	    }
	    if (g[gno].timeinfo.active == ON) {
		write_time_info(g[gno].timeinfo, "    ", pp);
	    } else {
		fprintf(pp, "%s    timeinfo off\n", imbedstr);
	    }
	    write_velocity_scale(g[gno].vl, "    ", pp);
	    write_wind_scale(g[gno].wl, "    ", pp);

/*
	    if (g[gno].vl.active == ON) {
	        write_velocity_scale(g[gno].vl, "    ", pp);
	    } else {
		fprintf(pp, "%s    vscale off\n", imbedstr);
	    }
*/

	    for (i = 0; i < MAXHISTMARKERS; i++) {
		if (g[gno].hbox[i].display == ON) {
		    write_Hist_marker(g[gno].hbox[i], i, "", pp);
		}
	    }

	    for (i = 0; i < MAXZOOMBOXES; i++) {
		if (g[gno].zbox[i].display == ON) {
		    write_Zoom_box(g[gno].zbox[i], i, "", pp);
		}
	    }
	    for (i = 0; i < MAXSLICES; i++) {
		if (g[gno].sbox[i].display == ON) {
		    write_DisplaySlice(g[gno].sbox[i], i, "", pp);
		}
	    }
	    for (i = 0; i < MAXGRIDS; i++) {
		if (grid[i].active == ON) {
		    write_DisplayGrid(g[gno].grid[i], i, pp);
		}
	    }
	    for (i = 0; i < MAXTRANS; i++) {
		if (trans[i].active == ON) {
		    write_DisplayTransect(g[gno].trans[i], i, "", pp);
		}
	    }
	    for (i = 0; i < MAXTEANL; i++) {
	    }
	    for (i = 0; i < MAXADCIRC; i++) {
		if (flowt[i].active == ON) {
		    write_DisplayFlow(g[gno].flowt[i], ADCIRC, i, pp);
		}
	    }
	    for (i = 0; i < MAXADCIRC3D; i++) {
		if (g[gno].flow3d[i].display == ON) {
		    write_Display3dFlow(i, g[gno].flow3d[i], pp);
		}
	    }
	    write_Isolparms(g[gno].salip, "scalar", pp);
	    write_Isolparms(g[gno].velmagip, "magnitude", pp);
	}
    }
}

void put_annotation(int gno, FILE * pp, int imbed)
{
    int i;
    boxtype b;
    linetype l;
    plotstr s;
    char imbedstr[2];

    if (imbed) {
	strcpy(imbedstr, "@");
    } else {
	imbedstr[0] = 0;
    }
    for (i = 0; i < MAXBOXES; i++) {
	get_graph_box(i, &b);
	if (b.active == ON) {
	    fprintf(pp, "%swith box\n", imbedstr);
	    fprintf(pp, "%s    box on\n", imbedstr);
	    fprintf(pp, "%s    box loctype %s\n", imbedstr,
		    w_or_v(b.loctype));
	    if (b.loctype == WORLD) {
		fprintf(pp, "%s    box g%1d\n", imbedstr, b.gno);
	    }
	    fprintf(pp, "%s    box %.12lg, %.12lg, %.12lg, %.12lg\n",
		    imbedstr, b.x1, b.y1, b.x2, b.y2);
	    fprintf(pp, "%s    box linestyle %d\n", imbedstr, b.lines);
	    fprintf(pp, "%s    box linewidth %d\n", imbedstr, b.linew);
	    fprintf(pp, "%s    box color %d\n", imbedstr, b.color);
	    switch (b.fill) {
	    case NONE:
		fprintf(pp, "%s    box fill none\n", imbedstr);
		break;
	    case COLOR:
		fprintf(pp, "%s    box fill color\n", imbedstr);
		break;
	    case PATTERN:
		fprintf(pp, "%s    box fill pattern\n", imbedstr);
		break;
	    }
	    fprintf(pp, "%s    box fill color %d\n", imbedstr,
		    b.fillcolor);
	    fprintf(pp, "%s    box fill pattern %d\n", imbedstr,
		    b.fillpattern);
	    fprintf(pp, "%sbox def\n", imbedstr);
	}
    }

    for (i = 0; i < MAXLINES; i++) {
	get_graph_line(i, &l);
	if (l.active == ON) {
	    fprintf(pp, "%swith line\n", imbedstr);
	    fprintf(pp, "%s    line on\n", imbedstr);
	    fprintf(pp, "%s    line loctype %s\n", imbedstr,
		    w_or_v(l.loctype));
	    if (l.loctype == WORLD) {
		fprintf(pp, "%s    line g%1d\n", imbedstr, l.gno);
	    }
	    fprintf(pp, "%s    line %.12lg, %.12lg, %.12lg, %.12lg\n",
		    imbedstr, l.x1, l.y1, l.x2, l.y2);
	    fprintf(pp, "%s    line linewidth %d\n", imbedstr, l.linew);
	    fprintf(pp, "%s    line linestyle %d\n", imbedstr, l.lines);
	    fprintf(pp, "%s    line color %d\n", imbedstr, l.color);
	    fprintf(pp, "%s    line arrow %d\n", imbedstr, l.arrow);
	    fprintf(pp, "%s    line arrow size %lf\n", imbedstr, l.asize);
	    fprintf(pp, "%s    line arrow type %d\n", imbedstr, l.atype);
	    fprintf(pp, "%sline def\n", imbedstr);
	}
    }

    for (i = 0; i < MAXSTR; i++) {
	get_graph_string(i, &s);
	if (s.active == ON && s.s[0]) {
	    fprintf(pp, "%swith string\n", imbedstr);
	    fprintf(pp, "%s    string on\n", imbedstr);
	    fprintf(pp, "%s    string loctype %s\n", imbedstr,
		    w_or_v(s.loctype));
	    if (s.loctype == WORLD) {
		fprintf(pp, "%s    string g%1d\n", imbedstr, s.gno);
	    }
	    fprintf(pp, "%s    string %.12lg, %.12lg\n", imbedstr, s.x,
		    s.y);
	    fprintf(pp, "%s    string linewidth %d\n", imbedstr, s.linew);
	    fprintf(pp, "%s    string color %d\n", imbedstr, s.color);
	    fprintf(pp, "%s    string rot %d\n", imbedstr, s.rot);
	    fprintf(pp, "%s    string font %d\n", imbedstr, s.font);
	    fprintf(pp, "%s    string just %d\n", imbedstr, s.just);
	    fprintf(pp, "%s    string char size %lf\n", imbedstr,
		    s.charsize);
	    fprintf(pp, "%sstring def \"%s\"\n", imbedstr, s.s);
	}
    }
}

void put_region(int gno, FILE * pp, int imbed)
{
    int i, j;
    char imbedstr[2];

    if (imbed) {
	strcpy(imbedstr, "@");
    } else {
	imbedstr[0] = 0;
    }
    for (i = 0; i < MAXREGION; i++) {
	if (rg[i].active == ON) {
	    fprintf(pp, "%sr%1d ON\n", imbedstr, i);
	    switch (rg[i].type) {
	    case ABOVE:
		fprintf(pp, "%sr%1d type above\n", imbedstr, i);
		break;
	    case BELOW:
		fprintf(pp, "%sr%1d type below\n", imbedstr, i);
		break;
	    case LEFT:
		fprintf(pp, "%sr%1d type left\n", imbedstr, i);
		break;
	    case RIGHT:
		fprintf(pp, "%sr%1d type right\n", imbedstr, i);
		break;
	    case POLYI:
		fprintf(pp, "%sr%1d type polyi\n", imbedstr, i);
		break;
	    case POLYO:
		fprintf(pp, "%sr%1d type polyo\n", imbedstr, i);
		break;
	    }
	    fprintf(pp, "%sr%1d linestyle %d\n", imbedstr, i, rg[i].lines);
	    fprintf(pp, "%sr%1d linewidth %d\n", imbedstr, i, rg[i].linew);
	    fprintf(pp, "%sr%1d color %d\n", imbedstr, i, rg[i].color);
	    if (rg[i].type != POLYI && rg[i].type != POLYO) {
		fprintf(pp, "%sr%1d line %.12lg, %.12lg, %.12lg, %.12lg\n",
			imbedstr, i, rg[i].x1, rg[i].y1, rg[i].x2,
			rg[i].y2);
	    } else {
		if (rg[i].x != NULL) {
		    for (j = 0; j < rg[i].n; j++) {
			fprintf(pp, "%sr%1d xy %.12lg, %.12lg\n", imbedstr,
				i, rg[i].x[j], rg[i].y[j]);
		    }
		}
	    }
	    for (j = 0; j < maxgraph; j++) {
		if (rg[i].linkto[j] == TRUE) {
		    fprintf(pp, "%slink r%1d to g%1d\n", imbedstr, i, j);
		}
	    }
	}
    }
}

static int write_DisplayFlow(DisplayFlow d, int flowtype, int flowno,
			     FILE * pp)
{
    int i, j;
    char s[256], ftype[256];

    if (flowtype == ADCIRC) {
	strcpy(ftype, "adcirc");
    } else if (flowtype == TEANL) {
	strcpy(ftype, "teanl");
    }
    if (d.display == OFF) {
	strcpy(s, "off");
    } else if (d.display == CENTER) {
	strcpy(s, "center");
    } else if (d.display == NODES) {
	strcpy(s, "nodes");
    }
    fprintf(pp, "with %s %d\n", ftype, flowno);
    fprintf(pp, "%s display %s\n", ftype, s);
    fprintf(pp, "%s display elev %s\n", ftype, on_or_off(d.display_elev));
    fprintf(pp, "%s display elev depth %s\n", ftype,
	    on_or_off(d.display_elevdepth));
    fprintf(pp, "%s display elev value %s\n", ftype,
	    on_or_off(d.display_maxelevval));
    fprintf(pp, "%s display elev max %s\n", ftype,
	    on_or_off(d.display_maxelev));
    fprintf(pp, "%s display elev markers %s\n", ftype,
	    on_or_off(d.display_elevmarkers));
    fprintf(pp, "%s display elev amp %s\n", ftype,
	    on_or_off(d.display_amp));
    fprintf(pp, "%s display elev phase %s\n", ftype,
	    on_or_off(d.display_phase));
    fprintf(pp, "%s display flow mag %s\n", ftype,
	    on_or_off(d.display_mag));
    fprintf(pp, "%s display inundation %s\n", ftype,
	    on_or_off(d.display_inun));
    fprintf(pp, "%s display flow wind %s\n", ftype,
	    on_or_off(d.display_wind));
    fprintf(pp, "%s color %d\n", ftype, d.p.color);
    if (d.display_elev == ON) {
	sprintf(s, "%s elev", ftype);
	write_Isolparms(d.elevip, s, pp);
    }
    if (d.display_maxelev == ON) {
	sprintf(s, "%s elev max", ftype);
	write_Isolparms(d.maxelevip, s, pp);
    }
    if (d.display_amp == ON) {
	sprintf(s, "%s amp", ftype);
	write_Isolparms(d.ampip, s, pp);
    }
    if (d.display_phase == ON) {
	sprintf(s, "%s phase", ftype);
	write_Isolparms(d.phaseip, s, pp);
    }
    if (d.display_mag == ON) {
	sprintf(s, "%s flow mag", ftype);
	write_Isolparms(d.magip, s, pp);
    }
    fprintf(pp, "%s sample flow %s\n", ftype, on_or_off(d.sample));
    if (d.sample == ON) {
	int i;
	for (i = 0; i < d.nsamples; i++) {
	    fprintf(pp, "%s sample flow node %d\n", ftype,
		    d.samples[i] + 1);
	}
    }
}
