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
 * Draw ticks and tick labels
 */

#ifndef lint
static char RCSid[] = "$Id: drawticks.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "globals.h"
#include "externs.h"

extern double devcharsize;
double setcharsize();
int stringextenty(double size, char *string);

void drawxtickmarks(int gno, int caxis);
void drawytickmarks(int gno, int caxis);
void drawxticklabels(int gno, int caxis);
void drawyticklabels(int gno, int caxis);
void drawxaxisbar(int gno, int caxis);
void drawyaxisbar(int gno, int caxis);
void drawyminorticklabels(int gno, int caxis);
void drawxminorticklabels(int gno, int caxis);

void drawaxes(int gno)
{
    int i;

    for (i = 0; i < 6; i++) {
	if (g[gno].t[i].active == ON) {
	    if (g[gno].t[i].t_drawbar == ON) {
		if (i % 2 == 0) {
		    drawxaxisbar(gno, i);
		} else {
		    drawyaxisbar(gno, i);
		}
	    }
	    if (g[gno].t[i].t_flag == ON) {
		if (i % 2 == 0) {
		    drawxtickmarks(gno, i);
		} else {
		    drawytickmarks(gno, i);
		}
	    }
	    if (g[gno].t[i].tl_flag == ON) {
		if (i % 2 == 0) {
		    drawxticklabels(gno, i);
		} else {
		    drawyticklabels(gno, i);
		}
	    }
	    if (g[gno].t[i].mtl_flag == ON) {
		if (i % 2 == 0) {
		    drawxminorticklabels(gno, i);
		} else {
		    drawyminorticklabels(gno, i);
		}
	    }
	}
    }
}

void drawxaxisbar(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot;
    double delb = 0.0, delt = 0.0, vyb, vyt, vx, ofx, ofb, oft;
    double start, stop;

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);

    /* set start and end for axis bar */
    start = w.xg1;
    stop = w.xg2;
    if (caxis == ZX_AXIS) {
	top = 0.0;
	bot = 0.0;
    } else {
	top = w.yg2;
	bot = w.yg1;
    }

    /* if offsx or offsy != 0.0 then adjust */
    vyb = v.yv1;
    vyt = v.yv2;
    vx = v.xv1;
    delb = 0.0;
    delt = 0.0;
    if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	vyb = v.yv1 - t.offsx;
	view2world(vx, vyb, &ofx, &ofb);
	delb = bot - ofb;
    }
    if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	vyt = v.yv2 + t.offsy;
	view2world(vx, vyt, &ofx, &oft);
	delt = oft - top;
    }
    setcolor(t.t_drawbarcolor);
    setlinestyle(t.t_drawbarlines);
    setlinewidth(t.t_drawbarlinew);
    setclipping(0);
    if (t.t_op == BOTTOM || t.t_op == BOTH) {
	my_move2(start, bot - delb);
	my_draw2(stop, bot - delb);
    }
    if (t.t_op == TOP || t.t_op == BOTH) {
	my_move2(start, top + delt);
	my_draw2(stop, top + delt);
    }
    setclipping(1);
}

void drawyaxisbar(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot;
    double delr = 0.0, dell = 0.0, vxl, vxr, vy, ofy, ofr, ofl;
    double start, stop;

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);

    start = w.yg1;
    stop = w.yg2;
    if (caxis == ZY_AXIS) {
	top = 0.0;
	bot = 0.0;
    } else {
	top = w.xg2;
	bot = w.xg1;
    }

/* if offsx or offsy != 0.0 then adjust */
    vxl = v.xv1;
    vxr = v.xv2;
    vy = v.yv1;
    dell = 0.0;
    delr = 0.0;
    if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	vxl = v.xv1 - t.offsx;
	view2world(vxl, vy, &ofl, &ofy);
	dell = bot - ofl;
    }
    if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	vxr = v.xv2 + t.offsy;
	view2world(vxr, vy, &ofr, &ofy);
	delr = ofr - top;
    }
    setcolor(t.t_drawbarcolor);
    setlinestyle(t.t_drawbarlines);
    setlinewidth(t.t_drawbarlinew);
    setclipping(0);
    if (t.t_op == LEFT || t.t_op == BOTH) {
	my_move2(bot - dell, start);
	my_draw2(bot - dell, stop);
    }
    if (t.t_op == RIGHT || t.t_op == BOTH) {
	my_move2(top + delr, start);
	my_draw2(top + delr, stop);
    }
    setclipping(1);

}

/* create format string */
void create_ticklabel(int form, int prec, double loc, char *s)
{
    char format[256];
    char tmpbuf[256];
    double tmp;
    int m, d, y, h, mm;
    double sec;
    int isign = 0;
    extern char *months[];
    extern char *monthl[];
    extern char *dayofweekstrs[];
    extern char *dayofweekstrl[];

    strcpy(format, "%.*lf");
    switch (form) {
    case DECIMAL:
	sprintf(s, format, prec, loc);
	tmp = atof(s);		/* fix reverse axes problem when loc == -0.0 */
	if (tmp == 0.0) {
	    strcpy(format, "%.*lf");
	    loc = 0.0;
	    sprintf(s, format, prec, loc);
	}
	break;
    case EXPONENTIAL:
	strcpy(format, "%.*le");
	sprintf(s, format, prec, loc);
	tmp = atof(s);		/* fix reverse axes problem when loc == -0.0 */
	if (tmp == 0.0) {
	    strcpy(format, "%.*le");
	    loc = 0.0;
	    sprintf(s, format, prec, loc);
	}
	break;
    case POWER:
	if (loc < 0.0) {
	    isign = 1;
	    loc = log10(-loc);
	    strcpy(format, "-10\\S%.*lf\\N");
	} else if (loc == 0.0) {
	    sprintf(format, "%.*lf", prec, 0.0);
	} else {
	    loc = log10(loc);
	    strcpy(format, "10\\S%.*lf\\N");
	}
	sprintf(s, format, prec, loc);
	break;
    case GENERAL:		/* precision is ignored */
	strcpy(format, "%.*lg");
	sprintf(s, format, prec, loc);
	tmp = atof(s);
	if (tmp == 0.0) {
	    strcpy(format, "%lg");
	    loc = 0.0;
	    sprintf(s, format, loc);
	}
	break;
    case DDMMYY:
	strcpy(format, "%d-%d-%d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	if (y >= 1900 && y < 2000) {
	    y -= 1900;
	}
	sprintf(s, format, d, m, y);
	break;
    case MMDDYY:
	strcpy(format, "%d-%d-%d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	if (y >= 1900 && y < 2000) {
	    y -= 1900;
	}
	sprintf(s, format, m, d, y);
	break;
    case MMYY:
	strcpy(format, "%d-%d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	if (y >= 1900 && y < 2000) {
	    y -= 1900;
	}
	sprintf(s, format, m, y);
	break;
    case MMDD:
	strcpy(format, "%d-%d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, m, d);
	break;
    case MONTHDAY:
	strcpy(format, "%s-%d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, months[m - 1], d);
	break;
    case DAYMONTH:
	strcpy(format, "%d-%s");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, d, months[m - 1]);
	break;
    case MONTHS:
	strcpy(format, "%s");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, months[m - 1]);
	break;
    case MONTHL:
	strcpy(format, "%s");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, monthl[m - 1]);
	break;
    case DAYOFWEEKS:
	strcpy(format, "%s");
	sprintf(s, format, dayofweekstrs[dayofweek(loc)]);
	break;
    case DAYOFWEEKL:
	strcpy(format, "%s");
	sprintf(s, format, dayofweekstrl[dayofweek(loc)]);
	break;
    case DAYOFYEAR:
	strcpy(format, "%d");
	sprintf(s, format, getndays(loc));
	break;
    case HMS:
	strcpy(format, "%02d:%02d:%02d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	sprintf(s, format, h, mm, (int) sec);
	break;
    case MMDDHMS:
	strcpy(format, "%d-%d %02d:%02d:%02d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	if (y >= 1900 && y < 2000) {
	    y -= 1900;
	}
	sprintf(s, format, m, d, h, mm, (int) sec);
	break;
    case MMDDYYHMS:
	strcpy(format, "%d-%d-%d %02d:%02d:%02d");
	calcdate(loc, &m, &d, &y, &h, &mm, &sec);
	if (y >= 1900 && y < 2000) {
	    y -= 1900;
	}
	sprintf(s, format, m, d, y, h, mm, (int) sec);
	break;
    case DEGREESLON:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%.*lfW");
	} else if (loc > 0.0) {
	    strcpy(format, "%.*lfE");
	} else {
	    strcpy(format, "0");
	}
	sprintf(s, format, prec, loc);
	break;
    case DEGREESMMLON:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%d %.*lf' W");
	} else if (loc > 0.0) {
	    strcpy(format, "%d %.*lf' E");
	} else {
	    strcpy(format, "0 0'");
	}
	y = loc;
	sec = (loc - y) * 60.0;
	sprintf(s, format, y, prec, sec);
	break;
    case DEGREESMMSSLON:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%d %d' %.*lf\" W");
	} else if (loc > 0.0) {
	    strcpy(format, "%d %d' %.*lf\" E");
	} else {
	    strcpy(format, "0 0' 0\"");
	}
	y = loc;
	sec = (loc - y) * 3600.0;
	m = sec / 60.0;
	sec = (sec - m * 60);
	sprintf(s, format, y, m, prec, sec);
	break;
    case MMSSLON:
	break;
    case DEGREESLAT:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%.*lfS");
	} else if (loc > 0.0) {
	    strcpy(format, "%.*lfN");
	} else {
	    strcpy(format, "0");
	}
	sprintf(s, format, prec, loc);
	break;
    case DEGREESMMLAT:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%d %.*lf' S");
	} else if (loc > 0.0) {
	    strcpy(format, "%d %.*lf' N");
	} else {
	    strcpy(format, "0 0'");
	}
	y = loc;
	sec = (loc - y) * 60.0;
	sprintf(s, format, y, prec, sec);
	break;
    case DEGREESMMSSLAT:
	if (loc < 0.0) {
	    loc *= -1.0;
	    strcpy(format, "%d %d' %.*lf\" S");
	} else if (loc > 0.0) {
	    strcpy(format, "%d %d' %.*lf\" N");
	} else {
	    strcpy(format, "0 0' 0\"");
	}
	y = loc;
	sec = (loc - y) * 3600.0;
	m = sec / 60.0;
	sec = (sec - m * 60);
	sprintf(s, format, y, m, prec, sec);
	break;
    case MMSSLAT:
	break;
    default:
	sprintf(s, format, prec, loc);
	break;
    }
}

void drawxticklabels(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    char s[256];
    double loc, x, y, top, bot;
    int i = 0, ifudge, iy, nticks, bjust = 2, tjust = 2;
    int tang;
    double delb = 0.0, delt = 0.0, vyb, vyt, vx, ofx, ofb, oft;
    double start, stop;
    double dx, dt, ttmp;	/* for alternate mapping */

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    t.tl_skip++;
    ifudge = (int) (t.tl_vgap * stringextenty(t.tl_charsize * devcharsize, "Ny"));

    /* set start and end for ticks */
    start = w.xg1;
    stop = w.xg2;
    if (t.alt == ON) {
	dt = (t.tmax - t.tmin);
	dx = (stop - start);
	ttmp = t.tmajor / dt * dx;
    } else {
	ttmp = t.tmajor;
    }
    if (caxis == ZX_AXIS) {
	top = 0.0;
	bot = 0.0;
    } else {
	top = w.yg2;
	bot = w.yg1;
    }
    if (t.tl_starttype == SPEC) {
	if (t.alt == ON) {
	    start = (t.tl_start - t.tmin) / dt * dx + w.xg1;
	} else {
	    if (t.tl_start >= w.xg1) {
		start = t.tl_start;
	    }
	}
    }
    if (t.tl_stoptype == SPEC) {
	if (t.alt == ON) {
	    stop = (t.tl_stop - t.tmin) / dt * dx + w.xg1;
	} else {
	    if (t.tl_stop <= w.xg2) {
		stop = t.tl_stop;
	    }
	}
    }
    /* how many ticklabels */
    if (islogx(gno)) {
	nticks = (int) fabs((log10(stop) - log10(start)) / ttmp) + 1;
    } else {
	nticks = (int) fabs((stop - start) / ttmp) + 1;
    }

    /* if offsx or offsy != 0.0 then adjust */
    vyb = v.yv1;
    vyt = v.yv2;
    vx = v.xv1;
    delb = 0.0;
    delt = 0.0;
    if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	vyb = v.yv1 - t.offsx;
	view2world(vx, vyb, &ofx, &ofb);
	delb = bot - ofb;
    }
    if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	vyt = v.yv2 + t.offsy;
	view2world(vx, vyt, &ofx, &oft);
	delt = oft - top;
    }
    setcharsize(t.tl_charsize);
    setfont(t.tl_font);
    setcolor(t.tl_color);
    setlinestyle(1);
    setlinewidth(t.tl_linew);

    if (t.tl_type == SPEC) {
	nticks = t.t_spec;
    }
    /* tick label justification */
    switch (t.tl_layout) {
    case HORIZONTAL:
	tang = 0;
	break;
    case VERTICAL:
	tang = 90;
	bjust = 1;
	tjust = 0;
	break;
    case SPEC:
	tang = t.tl_angle;
	if (tang > 0 && tang < 180) {
	    bjust = 1;
	    tjust = 0;
	} else if (tang > 180 && tang < 360) {
	    bjust = 0;
	    tjust = 1;
	} else if (tang == 0 || tang == 180 || tang == 360) {
	    bjust = 2;
	    tjust = 2;
	}
	break;
    }
    if (nticks > 500) {
	errwin("Too many ticks ( > 500 )");
	return;
    }
/* draw tick labels */
    for (i = 0; i < nticks; i += t.tl_skip) {
	if (t.t_type == SPEC) {
	    loc = t.t_specloc[i];
	} else {
	    if (t.alt == ON) {
		loc = ((start - w.xg1) / dx * dt + t.tmin) + i * t.tmajor;
	    } else {
		if (islogx(gno)) {
		    ttmp = ceil(log10(start));
		    loc = pow(10.0, ttmp + i * t.tmajor);
		} else {
		    loc = start + i * t.tmajor;
		}
	    }
	}
	if (!(t.alt == ON)) {
	    if (loc < w.xg1 || loc > w.xg2) {
		continue;
	    }
	}
	if (t.alt == ON) {
	    x = (loc - t.tmin) / dt * dx + w.xg1;
	} else {
	    x = loc;
	}

	if (t.tl_type == SPEC) {
	    strcpy(s, t.t_speclab[i].s);
	} else {
	    switch (t.tl_sign) {
	    case NORMAL:
		break;
	    case ABSOLUTE:
		loc = fabs(loc);
		break;
	    case NEGATE:
		loc = -loc;
		break;
	    }
	    create_ticklabel(t.tl_format, t.tl_prec, loc, s);
	}
	if (t.tl_prestr[0]) {
	    char tmpbuf[300];

	    strcpy(tmpbuf, t.tl_prestr);
	    strcat(tmpbuf, s);
	    strcpy(s, tmpbuf);
	}
	if (t.tl_appstr[0]) {
	    strcat(s, t.tl_appstr);
	}
/*
	iy = ifudge * ((t.t_inout == BOTH || t.t_inout == OUT) ? 2 : 1);
*/
	iy = ifudge;
	if (t.tl_op == BOTTOM || t.tl_op == BOTH) {
	    (*devwritestr) ((*devconvx) (x), (*devconvy) (bot - delb) - iy * (1 + i % (t.tl_staggered + 1)), tang, s, bjust, 1);
	}
	if (t.tl_op == TOP || t.tl_op == BOTH) {
	    (*devwritestr) ((*devconvx) (x), (*devconvy) (top + delt) + iy * (1 + i % (t.tl_staggered + 1)), tang, s, tjust, 0);
	}
    }

/* TODO axis label layout */
    if (t.label.s[0]) {
	int lang1 = 0, ljust1 = 2;
	int lang2 = 0, ljust2 = 2;
	double tmp = ((*devconvx) (g[gno].w.xg2) + (*devconvx) (g[gno].w.xg1)) / 2.0;

	switch (t.label_layout) {
	case PARA:
	    if (t.tl_op == BOTTOM || t.tl_op == BOTH) {
		lang1 = 0;
		ljust1 = 2;
	    }
	    if (t.tl_op == TOP || t.tl_op == BOTH) {
		lang2 = 0;
		ljust2 = 2;
	    }
	    break;
	case PERP:
	    if (t.tl_op == BOTTOM || t.tl_op == BOTH) {
		lang1 = 90;
		ljust1 = 1;
	    }
	    if (t.tl_op == TOP || t.tl_op == BOTH) {
		lang2 = 90;
		ljust2 = 0;
	    }
	    break;
	}

	setcharsize(t.label.charsize);
	setfont(t.label.font);
	setcolor(t.label.color);
	setlinestyle(1);
	setlinewidth(t.label.linew);
	if (t.tl_op == BOTTOM || t.tl_op == BOTH) {
	    (*devwritestr) ((int) (tmp), (*devconvy) (g[gno].w.yg1 - delb) - 3 * iy, lang1, t.label.s, ljust1, 0);
	}
	if (t.tl_op == TOP || t.tl_op == BOTH) {
	    (*devwritestr) ((int) (tmp), (*devconvy) (g[gno].w.yg2 + delt) + 3 * iy, lang2, t.label.s, ljust2, 0);
	}
    }
}

void drawyticklabels(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    char s[256];
    double loc, x, y, top, bot;
    int i = 0, ifudge, ix, nticks;
    int ljust = 1, rjust = 0;
    int tang, ylabpos = 0;
    double delr = 0.0, dell = 0.0, vxl, vxr, vy, ofy, ofr, ofl;
    double start, stop;
    double dy, dt, ttmp;	/* for alternate mapping */

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    t.tl_skip++;
    ifudge = (int) (t.tl_hgap * stringextentx(t.tl_charsize * devcharsize, "M"));

    /* set start and end for ticks */
    start = w.yg1;
    stop = w.yg2;
    if (t.alt == ON) {
	dt = (t.tmax - t.tmin);
	dy = (stop - start);
	ttmp = t.tmajor / dt * dy;
    } else {
	ttmp = t.tmajor;
    }
    if (caxis == ZY_AXIS) {
	top = 0.0;
	bot = 0.0;
    } else {
	top = w.xg2;
	bot = w.xg1;
    }
    if (t.tl_starttype == SPEC) {
	if (t.alt == ON) {
	    start = (t.tl_start - t.tmin) / dt * dy + w.yg1;
	} else {
	    if (t.tl_start >= w.yg1) {
		start = t.tl_start;
	    }
	}
    }
    if (t.tl_stoptype == SPEC) {
	if (t.alt == ON) {
	    stop = (t.tl_stop - t.tmin) / dt * dy + w.yg1;
	} else {
	    if (t.tl_stop <= w.yg2) {
		stop = t.tl_stop;
	    }
	}
    }
    /* how many ticks and which direction to increment */
    if (islogy(gno)) {
	nticks = (int) fabs((log10(stop) - log10(start)) / ttmp) + 2;
    } else {
	nticks = (int) fabs((stop - start) / ttmp) + 1;
    }
    t.tmajor = fabs(t.tmajor);

/* if offsx or offsy != 0.0 then adjust */
    vxl = v.xv1;
    vxr = v.xv2;
    vy = v.yv1;
    dell = 0.0;
    delr = 0.0;
    if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	vxl = v.xv1 - t.offsx;
	view2world(vxl, vy, &ofl, &ofy);
	dell = bot - ofl;
    }
    if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	vxr = v.xv2 + t.offsy;
	view2world(vxr, vy, &ofr, &ofy);
	delr = ofr - top;
    }
    setcharsize(t.tl_charsize);
    setfont(t.tl_font);
    setcolor(t.tl_color);
    setlinestyle(1);
    setlinewidth(t.tl_linew);

    if (t.tl_type == SPEC) {
	nticks = t.t_spec;
    }
    /* tick label justification */
    switch (t.tl_layout) {
    case HORIZONTAL:
	tang = 0;
	ljust = 1;
	rjust = 0;
	break;
    case VERTICAL:
	tang = 90;
	ljust = 2;
	rjust = 2;
	break;
    case SPEC:
	tang = t.tl_angle;
	if (tang > 0 && tang < 90) {
	    ljust = 1;
	    rjust = 0;
	} else if (tang > 90 && tang < 270) {
	    ljust = 0;
	    rjust = 1;
	} else if (tang == 90 || tang == 270) {
	    ljust = 2;
	    rjust = 2;
	} else if (tang == 0 || tang == 360) {
	    ljust = 1;
	    rjust = 0;
	}
	break;
    }
    if (nticks > 500) {
	errwin("Too many ticks ( > 500 )");
	return;
    }
/* draw tick labels */
    for (i = 0; i < nticks; i += t.tl_skip) {
	if (t.t_type == SPEC) {
	    loc = t.t_specloc[i];
	} else {
	    if (t.alt == ON) {
		loc = ((start - w.yg1) / dy * dt + t.tmin) + i * t.tmajor;
	    } else {
		if (islogy(gno)) {
		    ttmp = ceil(log10(start));
		    loc = pow(10.0, ttmp + i * t.tmajor);
		} else {
		    loc = start + i * t.tmajor;
		}
	    }
	}
	if (!(t.alt == ON)) {
	    if (loc < w.yg1 || loc > w.yg2) {
		continue;
	    }
	}
	if (t.alt == ON) {
	    y = (loc - t.tmin) / dt * dy + w.yg1;
	} else {
	    y = loc;
	}

	if (t.tl_type == SPEC) {
	    strcpy(s, t.t_speclab[i].s);
	} else {
	    switch (t.tl_sign) {
	    case NORMAL:
		break;
	    case ABSOLUTE:
		loc = fabs(loc);
		break;
	    case NEGATE:
		loc = -loc;
		break;
	    }
	    create_ticklabel(t.tl_format, t.tl_prec, loc, s);
	}
	if (t.tl_prestr[0]) {
	    char tmpbuf[300];

	    strcpy(tmpbuf, t.tl_prestr);
	    strcat(tmpbuf, s);
	    strcpy(s, tmpbuf);
	}
	if (t.tl_appstr[0]) {
	    strcat(s, t.tl_appstr);
	}
	if (t.t_inout == BOTH || t.t_inout == OUT) {
	    ix = ifudge * 2;
	} else {
	    ix = ifudge;
	}
	if (t.tl_op == LEFT || t.tl_op == BOTH) {
	    if (tang == 90) {
		ylabpos = 2 * ix;
		(*devwritestr) ((*devconvx) (bot - dell) - ix, (*devconvy) (y), tang, s, ljust, 1);
	    } else {
		ix = ifudge + stringextentx(t.tl_charsize * devcharsize, s);
		if (ix > ylabpos) {
		    ylabpos = ix;
		}
		(*devwritestr) ((*devconvx) (bot - dell) - ifudge, (*devconvy) (y), tang, s, ljust, 1);
	    }
	}
	if (t.tl_op == RIGHT || t.tl_op == BOTH) {
	    if (tang == 90) {
		ylabpos = 2 * ix;
		(*devwritestr) ((*devconvx) (top + delr) + ix, (*devconvy) (y), tang, s, rjust, 0);
	    } else {
		ix = ifudge + stringextentx(t.tl_charsize * devcharsize, s);
		if (ix > ylabpos) {
		    ylabpos = ix;
		}
		(*devwritestr) ((*devconvx) (top + delr) + ifudge, (*devconvy) (y), tang, s, rjust, 1);
	    }
	}
    }
/* TODO axis label layout */
    if (t.label.s[0]) {
	double tmp;
	int lang1 = 0, ljust1 = 2;
	int lang2 = 0, ljust2 = 2;

	tmp = ((*devconvy) (g[gno].w.yg2) + (*devconvy) (g[gno].w.yg1)) / 2.0;
	switch (t.label_layout) {
	case PARA:
	    if (t.tl_op == LEFT || t.tl_op == BOTH) {
		lang1 = 90;
		ljust1 = 2;
	    }
	    if (t.tl_op == RIGHT || t.tl_op == BOTH) {
		lang2 = 270;
		ljust2 = 2;
	    }
	    break;
	case PERP:
	    if (t.tl_op == LEFT || t.tl_op == BOTH) {
		lang1 = 0;
		ljust1 = 1;
	    }
	    if (t.tl_op == RIGHT || t.tl_op == BOTH) {
		lang2 = 0;
		ljust2 = 0;
	    }
	    break;
	}

	setcharsize(t.label.charsize);
	setfont(t.label.font);
	setcolor(t.label.color);
	setlinestyle(1);
	setlinewidth(t.label.linew);
	if (t.tl_op == LEFT || t.tl_op == BOTH) {
	    (*devwritestr) ((*devconvx) (g[gno].w.xg1 - dell) - (int) (1.4 * ylabpos), (int) (tmp), lang1, t.label.s, ljust1, 0);
	}
	if (t.tl_op == RIGHT || t.tl_op == BOTH) {
	    (*devwritestr) ((*devconvx) (g[gno].w.xg2 + delr) + (int) (1.4 * ylabpos), (int) (tmp), lang2, t.label.s, ljust2, 0);
	}
    }
}

/*
 */
void drawxtickmarks(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot, s, start, stop, step;
    double delb = 0.0, delt = 0.0, vyb, vyt, vx, ofx, ofb, oft;
    double dx, dt, ttmp1, ttmp2;	/* for alternate mapping */

    int axis, nticks, i, j, k, logloop, kstart;
    double ticl[2];

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    axis = caxis % 2;		/* 0 if an axis along x, 1 if along y */
    ticl[0] = t.t_size;
    ticl[1] = t.t_msize;
    for (i = 1; i >= 0; i--) {
	if (i == 1) {
	    setcolor(t.t_mcolor);
	    setlinestyle(t.t_mlines);
	    setlinewidth(t.t_mlinew);
	}
	if (i == 0) {
	    setcolor(t.t_color);
	    setlinestyle(t.t_lines);
	    setlinewidth(t.t_linew);
	}
	setticksize(ticl[i], 1.0);
	start = w.xg1;
	stop = w.xg2;
	if (t.tl_starttype == SPEC) {
	    if (i == 0) {
		if (t.tl_start >= w.xg1 && t.tl_start <= w.xg2) {
		    start = t.tl_start;
		}
	    } else {
		if (t.tl_start >= w.xg1 && t.tl_start <= w.xg2) {
		    start = t.tl_start - (int) ((t.tl_start - w.xg1) / t.tminor) * t.tminor;
		}
	    }
	}
	if (t.tl_stoptype == SPEC) {
	    if (i == 0) {
		if (t.tl_stop <= w.xg2 && t.tl_stop >= w.xg1) {
		    stop = t.tl_stop;
		}
	    } else {
		if (t.tl_stop <= w.xg2 && t.tl_stop >= w.xg1) {
		    stop = t.tl_stop + (int) ((w.xg2 - t.tl_stop) / t.tminor) * t.tminor;
		}
	    }
	}
	switch (caxis) {
	case X_AXIS:
	case XA_AXIS:
	    top = w.yg2;
	    bot = w.yg1;
	    break;
	case ZX_AXIS:
	    top = 0.0;
	    bot = 0.0;
	    break;
	}
	if (t.alt == ON) {
	    dt = (t.tmax - t.tmin);
	    dx = (stop - start);
	    ttmp1 = t.tmajor / dt * dx;
	    ttmp2 = t.tminor / dt * dx;
	} else {
	    ttmp1 = t.tmajor;
	    ttmp2 = t.tminor;
	}

	if (i == 0) {
	    step = ttmp1;
	} else {
	    step = ttmp2;
	}

	/* if offsx or offsy != 0.0 then adjust  */
	vyb = v.yv1;
	vyt = v.yv2;
	vx = v.xv1;
	delb = 0.0;
	delt = 0.0;
	if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	    vyb = v.yv1 - t.offsx;
	    view2world(vx, vyb, &ofx, &ofb);
	    delb = bot - ofb;
	}
	if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	    vyt = v.yv2 + t.offsy;
	    view2world(vx, vyt, &ofx, &oft);
	    delt = oft - top;
	}
	/* determine the number of tick marks */
	if (t.t_type == SPEC) {
	    nticks = t.t_spec;
	} else if (islogx(gno) && i == 0) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / step) + 2;
	} else if (islogx(gno) && i == 1) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / ttmp1) + 2;
	} else {
	    nticks = (int) fabs((stop - start) / step) + 1;
	}

	if (islogx(gno)) {
	    s = pow(10.0, floor(log10(start)));
	} else {
	    s = start;
	}

	if (nticks > 500) {
	    errwin("Too many ticks ( > 500 )");
	    return;
	}
	if (t.t_op == BOTTOM || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.xg1 || start > w.xg2) {
			continue;
		    }
		} else {
		    if (islogx(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogx(gno) && i == 1) {
			logloop = t.tminor;
			s = pow(10.0, floor(log10(s)));
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}
		if (logloop != 0) {
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogx(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip1;
			    }
			}
			if (start < w.xg1 || start > w.xg2) {
			    continue;
			}
			if (t.t_gridflag == ON && i == 0) {
			    my_move2(start, bot);
			    my_draw2(start, top);
			} else if (t.t_mgridflag == ON && i == 1) {
			    my_move2(start, bot);
			    my_draw2(start, top);
			} else {
			    switch (t.t_inout) {
			    case IN:
				drawtic(start, bot - delb, axis, 0);	/* draw up from y1 */
				break;
			    case OUT:
				drawtic(start, bot - delb, axis, 1);	/* draw down from y2 */
				break;
			    case BOTH:
				drawtic(start, bot - delb, axis, 0);	/* draw up from y1 */
				drawtic(start, bot - delb, axis, 1);	/* draw down from y2 */
				break;
			    }
			}
		    }
		}
	    }
	}
      skip1:;

	if (t.t_gridflag == ON && i == 0) {
	    continue;
	} else if (t.t_mgridflag == ON && i == 1) {
	    continue;
	}
	start = s;
	if (t.t_op == TOP || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.xg1 || start > w.xg2) {
			continue;
		    }
		} else {
		    if (islogx(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogx(gno) && i == 1) {
			logloop = t.tminor;
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}
		if (logloop != 0) {
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogx(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip2;
			    }
			}
			if (start < w.xg1 || start > w.xg2) {
			    continue;
			}
			switch (t.t_inout) {
			case IN:
			    drawtic(start, top + delt, axis, 1);	/* draw down from y2 */
			    break;
			case OUT:
			    drawtic(start, top + delt, axis, 0);	/* draw up from y2 */
			    break;
			case BOTH:
			    drawtic(start, top + delt, axis, 1);	/* draw down from y2 */
			    drawtic(start, top + delt, axis, 0);	/* draw up from y2 */
			    break;
			}
		    }
		}
	    }
	}
      skip2:;
    }
}

/*
 */
void drawytickmarks(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot, s, start, stop, step;

    double delr = 0.0, dell = 0.0, vxl, vxr, vy, ofy, ofr, ofl;

    double dx, dt, ttmp1, ttmp2;	/* for alternate mapping */

    int axis, nticks, i, j, k, logloop, kstart;
    double ticl[2];

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    axis = caxis % 2;		/* 0 if an axis along x, 1 if along y */
    ticl[0] = t.t_size;
    ticl[1] = t.t_msize;
    for (i = 1; i >= 0; i--) {
	if (i == 1) {
	    setcolor(t.t_mcolor);
	    setlinestyle(t.t_mlines);
	    setlinewidth(t.t_mlinew);
	}
	if (i == 0) {
	    setcolor(t.t_color);
	    setlinestyle(t.t_lines);
	    setlinewidth(t.t_linew);
	}
	setticksize(1.0, ticl[i]);
	start = w.yg1;
	stop = w.yg2;
	if (t.tl_starttype == SPEC) {
	    if (i == 0) {
		if (t.tl_start >= w.yg1 && t.tl_start <= w.yg2) {
		    start = t.tl_start;
		}
	    } else {
		if (t.tl_start >= w.yg1 && t.tl_start <= w.yg2) {
		    start = t.tl_start - (int) ((t.tl_start - w.yg1) / t.tminor) * t.tminor;
		}
	    }
	}
	if (t.tl_stoptype == SPEC) {
	    if (i == 0) {
		if (t.tl_stop <= w.yg2 && t.tl_stop >= w.yg1) {
		    stop = t.tl_stop;
		}
	    } else {
		if (t.tl_stop <= w.yg2 && t.tl_stop >= w.yg1) {
		    stop = t.tl_stop + (int) ((w.yg2 - t.tl_stop) / t.tminor) * t.tminor;
		}
	    }
	}
	switch (caxis) {
	case Y_AXIS:
	case YA_AXIS:
	    top = w.xg2;
	    bot = w.xg1;
	    break;
	case ZY_AXIS:
	    top = 0.0;
	    bot = 0.0;
	    break;
	}

	if (t.alt == ON) {
	    dt = (t.tmax - t.tmin);
	    dx = (stop - start);
	    ttmp1 = t.tmajor / dt * dx;
	    ttmp2 = t.tminor / dt * dx;
	} else {
	    ttmp1 = t.tmajor;
	    ttmp2 = t.tminor;
	}

	if (i == 0) {
	    step = ttmp1;
	} else {
	    step = ttmp2;
	}

	/* if offsx or offsy != 0.0 then adjust  */
	vxl = v.xv1;
	vxr = v.xv2;
	vy = v.yv1;
	dell = 0.0;
	delr = 0.0;
	if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	    vxl = v.xv1 - t.offsx;
	    view2world(vxl, vy, &ofl, &ofy);
	    dell = bot - ofl;
	}
	if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	    vxr = v.xv2 + t.offsy;
	    view2world(vxr, vy, &ofr, &ofy);
	    delr = ofr - top;
	}
	/* determine the number of tick marks */
	if (t.t_type == SPEC) {
	    nticks = t.t_spec;
	} else if (islogy(gno) && i == 0) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / step) + 2;
	} else if (islogy(gno) && i == 1) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / ttmp1) + 2;
	} else {
	    nticks = (int) fabs((stop - start) / step) + 1;
	}

	if (islogy(gno)) {
	    s = pow(10.0, floor(log10(start)));
	} else {
	    s = start;
	}
	if (nticks > 500) {
	    errwin("Too many ticks ( > 500 )");
	    return;
	}
	if (t.t_op == LEFT || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.yg1 || start > w.yg2) {
			continue;
		    }
		} else {
		    if (islogy(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogy(gno) && i == 1) {
			logloop = t.tminor;
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}

		if (logloop != 0) {
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogy(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip1;
			    }
			}
			if (start < w.yg1 || start > w.yg2) {
			    continue;
			}
			if (t.t_gridflag == ON && i == 0) {
			    my_move2(bot, start);
			    my_draw2(top, start);
			} else if (t.t_mgridflag == ON && i == 1) {
			    my_move2(bot, start);
			    my_draw2(top, start);
			} else {
			    switch (t.t_inout) {
			    case IN:
				drawtic(bot - dell, start, axis, 0);	/* draw up from y1 */
				break;
			    case OUT:
				drawtic(bot - dell, start, axis, 1);	/* draw down from y2 */
				break;
			    case BOTH:
				drawtic(bot - dell, start, axis, 0);	/* draw up from y1 */
				drawtic(bot - dell, start, axis, 1);	/* draw down from y2 */
				break;

			    }
			}
		    }
		}
	    }
	}
      skip1:;

	if (t.t_gridflag == ON && i == 0) {
	    continue;
	} else if (t.t_mgridflag == ON && i == 1) {
	    continue;
	}
	start = s;
	if (t.t_op == RIGHT || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.yg1 || start > w.yg2) {
			continue;
		    }
		} else {
		    if (islogy(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogy(gno) && i == 1) {
			logloop = t.tminor;
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}
		if (logloop != 0) {
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogy(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip2;
			    }
			}
			if (start < w.yg1 || start > w.yg2) {
			    continue;
			}
			switch (t.t_inout) {
			case IN:
			    drawtic(top + delr, start, axis, 1);	/* draw down from y2 */
			    break;
			case OUT:
			    drawtic(top + delr, start, axis, 0);	/* draw up from y2 */
			    break;
			case BOTH:
			    drawtic(top + delr, start, axis, 1);	/* draw down from y2 */
			    drawtic(top + delr, start, axis, 0);	/* draw up from y2 */
			    break;

			}
		    }
		}
	    }
	}
      skip2:;
    }
}

int check_nticks(int gno, int axis, double gmin, double gmax, double tm, int maxnt)
{
    return 1;
}

/*
 */
void drawyminorticklabels(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot, s, start, stop, step;

    double delr = 0.0, dell = 0.0, vxl, vxr, vy, ofy, ofr, ofl;

    double dx, dt, ttmp1, ttmp2;	/* for alternate mapping */

    int axis, nticks, i, j, k, logloop, kstart;
    double ticl[2];
    int ifudge;
    char ls[256];

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    ifudge = (int) (t.tl_hgap * stringextentx(t.tl_charsize * devcharsize, "M"));
    axis = caxis % 2;		/* 0 if an axis along x, 1 if along y */
    setcharsize(t.tl_charsize);
    setfont(t.tl_font);
    setcolor(t.tl_color);
    setlinestyle(1);
    setlinewidth(t.tl_linew);
    for (i = 1; i >= 0; i--) {
	start = w.yg1;
	stop = w.yg2;
	if (t.tl_starttype == SPEC) {
	    if (i == 0) {
		if (t.tl_start >= w.yg1 && t.tl_start <= w.yg2) {
		    start = t.tl_start;
		}
	    } else {
		if (t.tl_start >= w.yg1 && t.tl_start <= w.yg2) {
		    start = t.tl_start - (int) ((t.tl_start - w.yg1) / t.tminor) * t.tminor;
		}
	    }
	}
	if (t.tl_stoptype == SPEC) {
	    if (i == 0) {
		if (t.tl_stop <= w.yg2 && t.tl_stop >= w.yg1) {
		    stop = t.tl_stop;
		}
	    } else {
		if (t.tl_stop <= w.yg2 && t.tl_stop >= w.yg1) {
		    stop = t.tl_stop + (int) ((w.yg2 - t.tl_stop) / t.tminor) * t.tminor;
		}
	    }
	}
	switch (caxis) {
	case Y_AXIS:
	case YA_AXIS:
	    top = w.xg2;
	    bot = w.xg1;
	    break;
	case ZY_AXIS:
	    top = 0.0;
	    bot = 0.0;
	    break;
	}

	if (t.alt == ON) {
	    dt = (t.tmax - t.tmin);
	    dx = (stop - start);
	    ttmp1 = t.tmajor / dt * dx;
	    ttmp2 = t.tminor / dt * dx;
	} else {
	    ttmp1 = t.tmajor;
	    ttmp2 = t.tminor;
	}

	if (i == 0) {
	    step = ttmp1;
	} else {
	    step = ttmp2;
	}

	/* if offsx or offsy != 0.0 then adjust  */
	vxl = v.xv1;
	vxr = v.xv2;
	vy = v.yv1;
	dell = 0.0;
	delr = 0.0;
	if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	    vxl = v.xv1 - t.offsx;
	    view2world(vxl, vy, &ofl, &ofy);
	    dell = bot - ofl;
	}
	if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	    vxr = v.xv2 + t.offsy;
	    view2world(vxr, vy, &ofr, &ofy);
	    delr = ofr - top;
	}
	/* determine the number of tick marks */
	if (t.t_type == SPEC) {
	    nticks = t.t_spec;
	} else if (islogy(gno) && i == 0) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / step) + 2;
	} else if (islogy(gno) && i == 1) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / ttmp1) + 2;
	} else {
	    nticks = (int) fabs((stop - start) / step) + 1;
	}

	if (islogy(gno)) {
	    s = pow(10.0, floor(log10(start)));
	} else {
	    s = start;
	}
	if (nticks > 500) {
	    errwin("Too many ticks ( > 500 )");
	    return;
	}
	if (t.t_op == LEFT || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.yg1 || start > w.yg2) {
			continue;
		    }
		} else {
		    if (islogy(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogy(gno) && i == 1) {
			logloop = t.tminor;
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}

		if (logloop != 0) {
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogy(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip1;
			    }
			}
			if (start < w.yg1 || start > w.yg2) {
			    continue;
			}
			if (i == 1) {
			    int ix, tang = 0, ylabpos, ljust = 1, rjust = 0;
			    double loc = start;
			    double y = start;
			    switch (t.tl_sign) {
			    case NORMAL:
				break;
			    case ABSOLUTE:
				loc = fabs(loc);
				break;
			    case NEGATE:
				loc = -loc;
				break;
			    }
			    create_ticklabel(t.tl_format, t.tl_prec, loc, ls);
			    if (t.t_inout == BOTH || t.t_inout == OUT) {
				ix = ifudge * 2;
			    } else {
				ix = ifudge;
			    }
			    if (t.tl_op == LEFT || t.tl_op == BOTH) {
				if (tang == 90) {
				    ylabpos = 2 * ix;
				    (*devwritestr) ((*devconvx) (bot - dell) - ix, (*devconvy) (y), tang, ls, ljust, 1);
				} else {
				    ix = ifudge + stringextentx(t.tl_charsize * devcharsize, ls);
				    if (ix > ylabpos) {
					ylabpos = ix;
				    }
				    (*devwritestr) ((*devconvx) (bot - dell) - ifudge, (*devconvy) (y), tang, ls, ljust, 1);
				}
			    }
			    if (t.tl_op == RIGHT || t.tl_op == BOTH) {
				if (tang == 90) {
				    ylabpos = 2 * ix;
				    (*devwritestr) ((*devconvx) (top + delr) + ix, (*devconvy) (y), tang, s, rjust, 0);
				} else {
				    ix = ifudge + stringextentx(t.tl_charsize * devcharsize, s);
				    if (ix > ylabpos) {
					ylabpos = ix;
				    }
				    (*devwritestr) ((*devconvx) (top + delr) + ifudge, (*devconvy) (y), tang, s, rjust, 1);
				}
			    }
			}
		    }
		}
	    }
	}
      skip1:;
    }
}

/*
 */
void drawxminorticklabels(int gno, int caxis)
{
    tickmarks t;
    world w;
    view v;
    double top, bot, s, start, stop, step;
    double delb = 0.0, delt = 0.0, vyb, vyt, vx, ofx, ofb, oft;
    double dx, dt, ttmp1, ttmp2;	/* for alternate mapping */

    int axis, nticks, i, j, k, logloop, kstart;
    double ticl[2];

    get_graph_tickmarks(gno, &t, caxis);
    get_graph_world(gno, &w);
    get_graph_view(gno, &v);
    axis = caxis % 2;		/* 0 if an axis along x, 1 if along y */
    setcharsize(t.tl_charsize);
    setfont(t.tl_font);
    setcolor(t.tl_color);
    setlinestyle(1);
    setlinewidth(t.tl_linew);
    for (i = 1; i >= 0; i--) {
	start = w.xg1;
	stop = w.xg2;
	if (t.tl_starttype == SPEC) {
	    if (i == 0) {
		if (t.tl_start >= w.xg1 && t.tl_start <= w.xg2) {
		    start = t.tl_start;
		}
	    } else {
		if (t.tl_start >= w.xg1 && t.tl_start <= w.xg2) {
		    start = t.tl_start - (int) ((t.tl_start - w.xg1) / t.tminor) * t.tminor;
		}
	    }
	}
	if (t.tl_stoptype == SPEC) {
	    if (i == 0) {
		if (t.tl_stop <= w.xg2 && t.tl_stop >= w.xg1) {
		    stop = t.tl_stop;
		}
	    } else {
		if (t.tl_stop <= w.xg2 && t.tl_stop >= w.xg1) {
		    stop = t.tl_stop + (int) ((w.xg2 - t.tl_stop) / t.tminor) * t.tminor;
		}
	    }
	}
	switch (caxis) {
	case X_AXIS:
	case XA_AXIS:
	    top = w.yg2;
	    bot = w.yg1;
	    break;
	case ZX_AXIS:
	    top = 0.0;
	    bot = 0.0;
	    break;
	}
	if (t.alt == ON) {
	    dt = (t.tmax - t.tmin);
	    dx = (stop - start);
	    ttmp1 = t.tmajor / dt * dx;
	    ttmp2 = t.tminor / dt * dx;
	} else {
	    ttmp1 = t.tmajor;
	    ttmp2 = t.tminor;
	}

	if (i == 0) {
	    step = ttmp1;
	} else {
	    step = ttmp2;
	}

	/* if offsx or offsy != 0.0 then adjust  */
	vyb = v.yv1;
	vyt = v.yv2;
	vx = v.xv1;
	delb = 0.0;
	delt = 0.0;
	if (fabs(t.offsx) > 0.001 && fabs(t.offsx) < 1.0) {
	    vyb = v.yv1 - t.offsx;
	    view2world(vx, vyb, &ofx, &ofb);
	    delb = bot - ofb;
	}
	if (fabs(t.offsy) > 0.001 && fabs(t.offsy) < 1.0) {
	    vyt = v.yv2 + t.offsy;
	    view2world(vx, vyt, &ofx, &oft);
	    delt = oft - top;
	}
	/* determine the number of tick marks */
	if (t.t_type == SPEC) {
	    nticks = t.t_spec;
	} else if (islogx(gno) && i == 0) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / step) + 2;
	} else if (islogx(gno) && i == 1) {
	    nticks = (int) fabs((log10(stop) - log10(start)) / ttmp1) + 2;
	} else {
	    nticks = (int) fabs((stop - start) / step) + 1;
	}

	if (islogx(gno)) {
	    s = pow(10.0, floor(log10(start)));
	} else {
	    s = start;
	}

	if (nticks > 500) {
	    errwin("Too many ticks ( > 500 )");
	    return;
	}
	if (t.t_op == BOTTOM || t.t_op == BOTH) {
	    for (j = 0; j < nticks && start <= stop; j++) {
		logloop = 5;
		kstart = 0;
		if (t.t_type == SPEC) {
		    if (t.alt == ON) {
			start = (t.t_specloc[j] - t.tmin) / dt * dx + s;
		    } else {
			start = t.t_specloc[j];
		    }
		    if (start < w.xg1 || start > w.xg2) {
			continue;
		    }
		} else {
		    if (islogx(gno) && i == 0) {
			start = pow(10.0, log10(s) + j * t.tmajor);
		    } else if (islogx(gno) && i == 1) {
			logloop = t.tminor;
			s = pow(10.0, floor(log10(s)));
			start = pow(10.0, log10(s) + j * t.tmajor);
			kstart = 0;
		    } else {
			start = s + j * step;
		    }
		}
		if (logloop != 0) {
		    int ifudge, iy, bjust = 2, tjust = 2;
		    int tang = 0;
		    char ls[256];
		    double x, loc;
		    ifudge = (int) (t.tl_vgap * stringextenty(t.tl_charsize * devcharsize, "Ny"));
		    for (k = kstart; k < 10 / logloop - ((10 % logloop) ? 0 : 1); k++) {
			if (islogx(gno) && i == 1) {
			    start = ((k + 1) * logloop) * pow(10.0, log10(s) + j * t.tmajor);
			    if (start > stop) {
				goto skip1;
			    }
			}
			if (start < w.xg1 || start > w.xg2) {
			    continue;
			}
			x = loc = start;
			switch (t.tl_sign) {
			case NORMAL:
			    break;
			case ABSOLUTE:
			    loc = fabs(loc);
			    break;
			case NEGATE:
			    loc = -loc;
			    break;
			}
			create_ticklabel(t.tl_format, t.tl_prec, loc, ls);
			iy = ifudge;
			if (t.tl_op == BOTTOM || t.tl_op == BOTH) {
			    (*devwritestr) ((*devconvx) (x), (*devconvy) (bot - delb) - iy * (1 + i % (t.tl_staggered + 1)), tang, ls, bjust, 1);
			}
			if (t.tl_op == TOP || t.tl_op == BOTH) {
			    (*devwritestr) ((*devconvx) (x), (*devconvy) (top + delt) + iy * (1 + i % (t.tl_staggered + 1)), tang, ls, tjust, 0);
			}
		    }
		}
	    }
	}
      skip1:;
    }
}
