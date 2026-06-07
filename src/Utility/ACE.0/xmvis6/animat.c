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
 * animation routines
 */

#include "motifinc.h"
#include "defines.h"
#include "globals.h"
#include "symdefs.h"

#ifndef lint
static char RCSid[] = "$Id: animat.c,v 1.9 2007/01/16 23:32:58 pturner Exp $";

#endif

/* TODO */

static int prevstep = 0;
static double prevtime = 0.0;

void freshen_display(int restart);
void display_image(void);
void drawgraph(void);
double get_current_time(void);

static int restart = 0;

#include <Xm/Xm.h>

extern XtAppContext app_con;

XtIntervalId timer_id = 0;

static long timeint = 300;
extern int anim_state;

void do_ianimat(void);
void setirun(void);
void setistop(void);
void setistep(Widget w, int cd);
void set_next_image(void);

int isrunning(void)
{
    return timeclock.running;
}

void do_ianimat(void)
{
    display_image();
    set_next_image();
    if (timer_id) {
	timer_id = XtAppAddTimeOut(app_con, timeint, (XtTimerCallbackProc) do_ianimat, NULL);
    }
}

void set_wrap(int w)
{
    timeclock.wrap = w;
}

int get_wrap(void)
{
    return timeclock.wrap;
}

void goto_step(int step)
{
    if (step < timeclock.nsteps && step >= 0) {
	setistop();
	timeclock.curstep = step;
	timeclock.curtime = get_current_time();
	display_image();
    } else {
	errwin("Not that many time steps");
    }
}

void set_display(int s)
{
}

void setirun(void)
{
    if (timeclock.nsteps == 0) {
	timeclock.running = 0;
	return;
    }
    timeclock.running = 1;
    anim_state = 1;
    timer_id = XtAppAddTimeOut(app_con, timeint, (XtTimerCallbackProc) do_ianimat, NULL);
}

void setistop(void)
{
    anim_state = 0;
    timeclock.running = 0;
    timeclock.curstep = prevstep;
    timeclock.curtime = prevtime;
    if (timer_id) {
	XtRemoveTimeOut(timer_id);
    }
    timer_id = 0;
}

void setistep(Widget w, int cd)
{
    timeclock.dir = cd;
    if (timeclock.running) {
	setistop();
	set_next_image();
    }
    set_next_image();
    display_image();
}

/*
 * generate plots in batch mode, sequence of gifs
 */
void batchrunsteps(int start, int stop, int skip)
{
    extern int geoflag;
    extern char geoprintstr[];
    int i;
    int save;
    char savestr[1024];
    setistop();
    timeclock.running = 1;
    save = ptofile;
    strcpy(savestr, printstr);
    if (skip <= 0) {
	fprintf(stderr, "batchrunsteps(); skip <= 0\n");
	return;
    }
    for (i = start - 1; i < stop; i += skip) {
	if (i >= timeclock.nsteps || i < 0) {
	    break;
	}
	timeclock.curstep = (i % timeclock.nsteps);
	timeclock.curtime = get_current_time();
	ptofile = 1;
	sprintf(printstr, "%s%09.0lf.gif", batchprefix, get_current_time());
	if (geoflag) {
	    sprintf(geoprintstr, "%s%09.0lf.gfw", batchprefix, get_current_time());
	}
	do_batch_plot();
    }
    timeclock.running = 0;
    ptofile = save;
    strcpy(printstr, savestr);
}

void batchrunstep(int step)
{
    extern int geoflag;
    extern char geoprintstr[];
    int i;
    int save;
    char savestr[1024];
    setistop();
    timeclock.running = 1;
    save = ptofile;
    strcpy(savestr, printstr);
    ptofile = 1;
    sprintf(printstr, "%s%09.0lf.gif", batchprefix, 0.0);
    if (geoflag) {
	sprintf(geoprintstr, "%s%09.0lf.gfw", batchprefix, 0.0);
    }
    do_batch_plot();
    timeclock.running = 0;
    ptofile = save;
    strcpy(printstr, savestr);
}

void runsteps(int start, int stop)
{
    int i;

    setistop();
    timeclock.running = 1;
    for (i = start - 1; i < stop; i++) {
	if (i >= timeclock.nsteps) {
	    break;
	}
	timeclock.curstep = (i % timeclock.nsteps);
	timeclock.curtime = get_current_time();
	display_image();
    }
    timeclock.running = 0;
}

void create_fli(char *fname, int start, int stop, int res)
{
    int i, nf = 0;
    extern int save_fli;

    setistop();
    display_image();
    timeclock.running = 1;
    for (i = start - 1; i < stop; i++) {
	nf++;
    }
    save_fli = 1;
    initialize_fli(fname, nf, res);
    initmake_fli();
    for (i = start - 1; i < stop; i++) {
	if (i >= timeclock.nsteps) {
	    break;
	}
	timeclock.curstep = (i % timeclock.nsteps);
	timeclock.curtime = get_current_time();
	display_image();
	add_image(i - start + 2);	/* start with image counting from 1 */
    }
    save_fli = 0;
    timeclock.running = 0;
    close_fli();
}

void set_restart(void)
{
    if (timeclock.running) {
	restart = 1;
	setistop();
    } else {
	restart = 0;
    }
}

void freshen_display(int restart)
{
    doclear = 1;
    dobackground = 1;
    display_image();
    if (restart) {
	setirun();
    }
}

int setupanimation(void)
{
    int n, i = 0;
    int restart = 0;
    double current_time;

    if (!inwin) {
	return;
    }
    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    timeclock.curstep = get_current_step();
    current_time = get_current_time();
    doclear = 1;
    display_image();
    if (restart) {
	setirun();
    }
}

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

int autoscale_clock(void)
{
    double start = BIG, stop = MBIG, step = 1.0;
    int i;
    int nsteps = -1;
    int didone = 0, dtype = 0, dind;

    for (i = 0; i < MAXPATHLINES; i++) {
	if (object_isactive(DROGUES, i)) {
	    start = min(drogues[i].start, start);
	    stop = max(drogues[i].stop, stop);
	    if (nsteps < drogues[i].nsteps) {
		nsteps = drogues[i].nsteps;
		didone = 1;
		dind = i;
		dtype = DROGUES;
	    }
	}
    }

    for (i = 0; i < MAXHISTMARKERS; i++) {
    }

    for (i = 0; i < MAXELA; i++) {
	if (object_isactive(ELA, i)) {
	    start = min(elaconc[i].start, start);
	    stop = max(elaconc[i].stop, stop);
	    if (nsteps < elaconc[i].nsteps) {
		nsteps = elaconc[i].nsteps;
		didone = 1;
		dind = i;
		dtype = ELA;
	    }
	}
    }

/*
    for (i = 0; i < MAXTEANL; i++) {
	if (object_isactive(TEANL, i)) {
	    start = min(flowf[i].start, start);
	    stop = max(flowf[i].stop, stop);
	    if (nsteps < flowf[i].nsteps) {
		nsteps = flowf[i].nsteps;
	        didone = 1;
	        dind = i;
	        dtype = TEANL;
	    }
	}
    }
*/

    for (i = 0; i < MAXADCIRC; i++) {
	if (object_isactive(ADCIRC, FLOW, i) || object_isactive(ADCIRC, ELEV, i)) {
	    start = min(flowt[i].start, start);
	    stop = max(flowt[i].stop, stop);
	    if (nsteps < flowt[i].nsteps) {
		nsteps = flowt[i].nsteps;
		didone = 1;
		dind = i;
		dtype = ADCIRC;
	    }
	}
    }

    if (didone) {
	switch (dtype) {
	case ADCIRC:
	    break;
	case ELA:
	    break;
	case DROGUES:
	    break;
	}
	timeclock.start = start;
	timeclock.stop = stop;
	set_clock(0, start, stop, step, nsteps);
	if (timeclock.nsteps != 0) {
	    timeclock.step = (timeclock.start - timeclock.stop) / timeclock.nsteps;
	} else {
	}
	timeclock.curtime = start;
	timeclock.nsteps = nsteps;
	timeclock.curstep = 0;
    } else {
	timeclock.start = 0.0;
	timeclock.stop = 0.0;
	timeclock.step = 0.0;
	timeclock.curtime = 0.0;
	timeclock.nsteps = 0;
	timeclock.curstep = 0;
    }
}

void load_clock(int obj, int w, int w2)
{
    int i, j;
    sprintf(statusstr, "In load_clock(): %d %d %d", obj, w, w2);
    writelogfile(statusstr);
    switch (obj) {
    case GRID:
	if (gridt[w].active == ON) {
	    for (j = 0; j < gridt[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = gridt[w].grids[j].time;
		} else {
		}
	    }
	}
    case ELA:
	if (object_isactive(ELA, w)) {
	    for (j = 0; j < elaconc[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = elaconc[w].data[j].time;
		} else {
		}
	    }
	}
	break;
    case TRANSECT:
	if (object_isactive(TRANSECT, w)) {
	    for (j = 0; j < trans[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = trans[w].t[j];
		} else {
		}
	    }
	}
	break;
    case ADCIRC:
	if (object_isactive(ADCIRC, FLOW, w) || object_isactive(ADCIRC, ELEV, w)) {
	    for (j = 0; j < flowt[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = flowt[w].f[j].time;
		} else {
		}
	    }
	}
	break;
    case ADCIRC3DFLOW:
	if (object_isactive(ADCIRC3DFLOW, w)) {
	    for (j = 0; j < adc3d[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = adc3d[w].time[j];
		} else {
		}
	    }
	}
	break;
    case HISTORY:
	if (object_isactive(HISTORY, FLOW, w2)) {
	    for (j = 0; j < flowh[w2].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = flowh[w2].f[j].time;
		} else {
		}
	    }
	}
	break;
    case DROGUES:
	if (object_isactive(DROGUES, w)) {
	    for (j = 0; j < drogues[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = drogues[w].p[j].time;
		} else {
		}
	    }
	}
	break;
    case TRACK:
	if (object_isactive(TRACK, w)) {
	    for (j = 0; j < track[w].nsteps; j++) {
		if (j < timeclock.nsteps) {
		    timeclock.t[j] = track[w].time[j];
		} else {
		}
	    }
	}
	break;
    }
    timeclock.curstep = 0;
    sprintf(statusstr, "Done load_clock(): %d %d %d", obj, w, w2);
    writelogfile(statusstr);
}

int get_current_step(void)
{
    return timeclock.curstep;
}

double get_current_time(void)
{
    if (timeclock.t == NULL) {
	return 0.0;
    }
    return timeclock.t[timeclock.curstep];
}

void set_next_image(void)
{
    if (timeclock.nsteps == 0 || timeclock.curstep > timeclock.nsteps) {
	timeclock.curstep = timeclock.startstep;
	timeclock.curtime = timeclock.start;
	return;
    }
    prevstep = timeclock.curstep;
    prevtime = timeclock.curtime;
    timeclock.curstep += timeclock.dir;
    if (timeclock.curstep < 0) {
	if (timeclock.wrap) {
	    timeclock.curstep = timeclock.nsteps - 1;
	    timeclock.curtime = timeclock.stop;
	} else {
	    setistop();
	    timeclock.curstep = 0;
	    timeclock.curtime = timeclock.start;
	}
    } else if (timeclock.curstep == timeclock.nsteps) {
	if (timeclock.wrap) {
	    timeclock.curstep = 0;
	    timeclock.curtime = timeclock.start;
	} else {
	    setistop();
	    timeclock.curstep = timeclock.nsteps - 1;
	    timeclock.curtime = timeclock.stop;
	}
    }
    timeclock.curtime = get_current_time();
}

void setrewind(void)
{
    int restart = 0;

    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    prevstep = timeclock.curstep = 0;
    prevtime = timeclock.curtime = timeclock.start;
    doclear = 1;
    display_image();
}

void setreverse(void)
{
    int restart = 0;

    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    timeclock.dir = -1;
    if (restart) {
	setirun();
    }
}

void setfaster(void)
{
    timeint = (timeint - 50) < 100 ? 100 : timeint - 50;
}

void setslower(void)
{
    timeint += 50;
}

void setfastforward(void)
{
    if (timeclock.running) {
	setistop();
    }
    prevstep = timeclock.curstep = timeclock.nsteps - 1;
    prevtime = timeclock.curtime = timeclock.stop;
    doclear = 1;
    display_image();
}

void setfastreverse(void)
{
}

void setforward(void)
{
    int restart = 0;

    if (timeclock.running) {
	restart = 1;
	setistop();
    }
    timeclock.dir = 1;
    if (restart) {
	setirun();
    }
}

void setstreams(void)
{
}

void do_setreverse(void)
{
    setistop();
    setreverse();
    setirun();
}

void do_setfastforward(void)
{
    setfastforward();
}

void do_setfastreverse(void)
{
}

void do_run(void)
{
    setistop();
    setirun();
}

void do_setforward(void)
{
    setistop();
    setforward();
    setirun();
}
