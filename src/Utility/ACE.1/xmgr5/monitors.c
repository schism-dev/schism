/* $Id: monitors.c,v 1.2 2004/07/07 03:11:19 pturner Exp $

 * Monitors
 *
 */

#include <stdio.h>
#ifndef WIN32
#include <sys/param.h>
#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <fcntl.h>
#ifndef WIN32
#include <sys/time.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>

#include <Xm/XmAll.h>

#include "globals.h"
#include "motifinc.h"
#include "noxprotos.h"

typedef struct _Monitor {
    int type;
    int active;
    char fname[1024];
    int gno;
    int setno;
    int load;
    int prev;
    int cur;
    int fd;
    int filesize;
    int begin;
} Monitor;

static Monitor monitors[MAXMONITORS];

void processMonitors(void);
void defaultMonitor(int i);
void createMonitor(int mon, char *fname, int gno, int setno, int load);
void createQuickMonitor(char *fname);

static int fid;
static XtIntervalId tim = 0;
static XtTimerCallbackProc timercb(XtPointer cdp, XtIntervalId * id);
extern int timer_delay;
extern XtAppContext app_con;

static XtTimerCallbackProc timercb(XtPointer cdp, XtIntervalId * id)
{
    processMonitors();
    tim = XtAppAddTimeOut(app_con, timer_delay, (XtTimerCallbackProc) timercb, NULL);
    return 0;
}

void settimer(void)
{
    tim = XtAppAddTimeOut(app_con, timer_delay, (XtTimerCallbackProc) timercb, NULL);
}

void stoptimer(void)
{
    if (tim) {
	XtRemoveTimeOut(tim);
    }
    tim = 0;
}

static Monitor monitor[MAXMONITORS];

/*
 * Monitor a growing netcdf file
 */
void createMonitor(int mon, char *fname, int gno, int setno, int load)
{
    defaultMonitor(mon);
    monitors[mon].active = ON;
    monitors[mon].gno = gno;
    monitors[mon].setno = setno;
    strcpy(monitors[mon].fname, fname);
    monitors[mon].load = load;
}

void createQuickMonitor(char *fname)
{
    int mon = 0;
    defaultMonitor(mon);
    monitors[mon].active = ON;
    monitors[mon].gno = cg;
    monitors[mon].setno = 0;
    strcpy(monitors[mon].fname, fname);
    monitors[mon].load = 0;
}

void defaultMonitor(int i)
{
    monitors[i].active = OFF;
    monitors[i].type = 0;
    monitors[i].gno = -1;
    monitors[i].setno = -1;
    monitors[i].fname[0] = 0;
    monitors[i].load = 0;
    monitors[i].prev = -1;
    monitors[i].cur = 0;
}

void killMonitor(int i)
{
    defaultMonitor(i);
}

void startMonitor(int i)
{
    monitors[i].active = ON;
}

void stopMonitor(int i)
{
    monitors[i].active = OFF;
}

int newMonitor(void)
{
    int i;
    for (i = 0; i < MAXMONITORS; i++) {
	if (monitor[i].active == OFF) {
	    defaultMonitor(i);
	    return i;
	}
    }
    return -1;
}

void initMonitor(void)
{
    int i;
    for (i = 0; i < MAXMONITORS; i++) {
	defaultMonitor(i);
    }
}

#define BUFSIZE  512

void processMonitors(void)
{
    int i, cnt, doredraw = 0;
    for (i = 0; i < MAXMONITORS; i++) {
	if (monitors[i].active == ON) {
	    FILE *fp;
	    double *x, *y;
	    if ((fp = fopen(monitors[i].fname, "r")) == NULL) {
		continue;
	    }
	    if (isactive(monitors[i].gno, monitors[i].setno)) {
		softkillset(monitors[i].gno, monitors[i].setno);
	    }
	    x = (double *) calloc(BUFSIZE, sizeof(double));
	    y = (double *) calloc(BUFSIZE, sizeof(double));
	    if (x == NULL || y == NULL) {
		errwin("Insufficient memory for set");
		cxfree(x);
		cxfree(y);
		killset(monitors[i].gno, monitors[i].setno);
		continue;
	    }
	    cnt = 0;
	    while (fgets(buf, 255, fp) != NULL) {
		if (buf[0] == '#') {
		    continue;
		}
		sscanf(buf, "%lf %lf", &x[cnt], &y[cnt]);
		cnt++;
		if (cnt % BUFSIZE == 0) {
		    x = (double *) realloc(x, (cnt + BUFSIZE) * sizeof(double));
		    y = (double *) realloc(y, (cnt + BUFSIZE) * sizeof(double));
		}
	    }
	    x = (double *) realloc(x, cnt * sizeof(double));
	    y = (double *) realloc(y, cnt * sizeof(double));
	    activateset(monitors[i].gno, monitors[i].setno);
	    settype(monitors[i].gno, monitors[i].setno, XY);
	    setcol(monitors[i].gno, x, monitors[i].setno, cnt, 0);
	    setcol(monitors[i].gno, y, monitors[i].setno, cnt, 1);
	    setcomment(monitors[i].gno, monitors[i].setno, monitors[i].fname);
	    updatesetminmax(monitors[i].gno, monitors[i].setno);
	    doredraw = 1;
	}
    }
    if (doredraw) {
	drawgraph();
    }
}
