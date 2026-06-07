/* $Id: visutils.c,v 1.1.1.1 2003/07/21 16:18:42 pturner Exp $
 *
 * utilities for 1d visualization code
 *
 */

#include <stdio.h>
#include <math.h>
#include <Xm/Xm.h>

extern XtWorkProcId wd;

#include "defines.h"
#include "globals.h"
#include "vis.h"

#define MAXVIS 20

vis_struct vis[MAXVIS];

static int vis_mode = 0;

kill_vis(vno)
    int vno;
{
    softkillset(vis[vno].set1);
    vis[vno].x = NULL;
    vis[vno].y = NULL;
    vis[vno].z = NULL;
    if (vis[vno].xf != NULL) {
	free(vis[vno].xf);
    }
    if (vis[vno].yf != NULL) {
	free(vis[vno].yf);
    }
    if (vis[vno].zf != NULL) {
	free(vis[vno].zf);
    }
    vis[vno].active = OFF;
    vis[vno].skip = 1;
}

set_vismode(m)
{
    vis_mode = m;
}

int init_vis()
{
    int i;
    for (i = 0; i < MAXVIS; i++) {
	vis[i].active = OFF;
	vis[i].skip = 1;
    }
}

int nextvis()
{
    int i;
    for (i = 0; i < MAXVIS; i++) {
	if (vis[i].active == OFF) {
	    return i;
	}
    }
    return -1;
}

check_vistype(type)
    int type;
{
    return (type == 17 || type == 117 ||
	    type == 18 || type == 118 ||
	    type == 19 ||
	    type == 20 || type == 120);
}

static char timebuf[256];
static int time_str = -1;
double stx = 0.8, sty = 0.8;
extern int vis_time;

int read_vis()
{
    int i, str;
    int retval = 0;
    double wx, wy;
    for (i = 0; i < MAXVIS; i++) {
	if (vis[i].active == ON) {
	    switch (vis[i].src) {
	    case DISK:
	    case PIPE:
		retval |= readbin_vis(i, vis_mode);
		break;
	    case MEM:
		retval |= readmem_vis(i, vis_mode);
		break;
	    }
	}
    }
    view2world(stx, sty, &wx, &wy);
    if (time_str == -1) {
	time_str = define_string(timebuf, wx, wy);
    }
    strcpy(pstr[time_str].s, timebuf);

    if (vis_time) {
	pstr[time_str].active = ON;
    } else {
	pstr[time_str].active = OFF;
    }

    drawgraph();
/*
    kill_string(str);
    if (first_pic) {
        sleep(5);
        first_pic = 0;
    }
*/
    return retval;
}

int readmem_vis(vno, mode)
{
    int i, gno = vis[vno].gno, set1 = vis[vno].set1;
    double time, x1, y1, x2, y2, sum;
    float ftime;
    static int first = 1;
    int npts, type;
    switch (mode) {
    case 0:
	vis[vno].fp = NULL;
	vis[vno].npts = npts;
	if ((vis[vno].x = (double *) calloc(npts, sizeof(double))) == NULL) {
	    vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for x");
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	if ((vis[vno].y = (double *) calloc(npts, sizeof(double))) == NULL) {
	    cxfree(vis[vno].x);
	    vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for y");
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	softkillset(gno, set1);
	activateset(gno, set1);
	settype(gno, set1, XY);

	setcol(gno, vis[vno].x, set1, npts, 0);
	setcol(gno, vis[vno].y, set1, npts, 1);

	setcomment(vis[vno].gno, vis[vno].set1, vis[vno].fname);
	if (first) {
	    first = 0;
	}
	if ((fread(&time, sizeof(double), 1, vis[vno].fp) == 1) &&
	    fread(vis[vno].y, sizeof(double), vis[vno].npts, vis[vno].fp) == npts) {
	    updatesetminmax(gno, set1);
	    getsetminmax(gno, set1, &x1, &x2, &y1, &y2);
	    vis[vno].nsteps++;
	    if (!(vis[vno].nsteps % vis[vno].skip)) {
		sprintf(timebuf, "%.0lf s", time);
		my_doublebuffer(TRUE);
	    }
	    return 0;
	} else {
	    vis_error(vis[vno].fp, vis[vno].src, NULL);
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	break;
    case 1:
	vis_error(vis[vno].fp, vis[vno].src, NULL);
	vis[vno].fp = NULL;
	return TRUE;
	break;
    case 2:
	return 0;
	break;
    }
}

int readbin_vis(vno, mode)
    int vno, mode;
{
    int i, j, nsteps, gno = vis[vno].gno, set1 = vis[vno].set1;
    double time, x1, y1, x2, y2, sum;
    float ftime;
    static int first = 1;
    int npts, nsets, type;
    if (vis[vno].fp == NULL && mode == 0) {
	vis[vno].nsteps = 0;
	if (vis[vno].src == PIPE) {
	    vis[vno].fp = popen(vis[vno].fname, "r");
	} else if (vis[vno].src == DISK) {
	    vis[vno].fp = fopen(vis[vno].fname, "r");
	}
	if (fread(&type, sizeof(int), 1, vis[vno].fp) != 1) {
	    vis_error(vis[vno].fp, vis[vno].src, "Error reading file");
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	if (!check_vistype(type)) {
	    vis_error(vis[vno].fp, vis[vno].src, "Incorrect file type for visualization");
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	vis[vno].ftype = type;
	if (type == 20) {
	    if (fread(&nsteps, sizeof(int), 1, vis[vno].fp) != 1) {
		vis_error(vis[vno].fp, vis[vno].src, "Error reading file");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    if (fread(&nsets, sizeof(int), 1, vis[vno].fp) != 1) {
		vis_error(vis[vno].fp, vis[vno].src, "Error reading file");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    if ((vis[vno].xm = (double **) calloc(nsets, sizeof(double *))) == NULL) {
		vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for **xm");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    if ((vis[vno].ym = (double **) calloc(nsets, sizeof(double *))) == NULL) {
		vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for **ym");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    for (i = 0; i < nsets; i++) {
		if (fread(&npts, sizeof(int), 1, vis[vno].fp) != 1) {
		    vis_error(vis[vno].fp, vis[vno].src, "Error reading file");
		    vis[vno].fp = NULL;
		    wd = NULL;
		    return TRUE;
		}
		vis[vno].slen[i] = npts;
		vis[vno].sets[i] = i;
		if ((vis[vno].xm[i] = (double *) calloc(npts, sizeof(double))) == NULL) {
		    vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for *xm");
		    vis[vno].fp = NULL;
		    wd = NULL;
		    return TRUE;
		}
		if ((vis[vno].ym[i] = (double *) calloc(npts, sizeof(double))) == NULL) {
		    vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for *ym");
		    vis[vno].fp = NULL;
		    wd = NULL;
		    return TRUE;
		}
		softkillset(gno, i);
		activateset(gno, i);
		settype(gno, i, XY);

		setcol(gno, vis[vno].xm[i], i, npts, 0);
		setcol(gno, vis[vno].ym[i], i, npts, 1);

		setcomment(vis[vno].gno, vis[vno].sets[i], vis[vno].fname);
	    if (fread(vis[vno].xm[i], sizeof(double), npts, vis[vno].fp) == npts &&
		fread(vis[vno].ym[i], sizeof(double), npts, vis[vno].fp) == npts) {
		updatesetminmax(gno, vis[vno].sets[i]);
		getsetminmax(gno, vis[vno].sets[i], &x1, &x2, &y1, &y2);
	    } else {
		vis_error(vis[vno].fp, vis[vno].src, NULL);
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }

	    }
	    vis[vno].nsets = nsets;
	    if (first) {
		first = 0;
	    }
	} else {
	    if (fread(&npts, sizeof(int), 1, vis[vno].fp) != 1) {
		vis_error(vis[vno].fp, vis[vno].src, "Error reading file");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    vis[vno].npts = npts;

	    if ((vis[vno].x = (double *) calloc(npts, sizeof(double))) == NULL) {
		vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for x");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    switch (type) {
	    case 18:
	    case 118:
		if (fread(vis[vno].x, sizeof(double), npts, vis[vno].fp) != npts) {
		    vis_error(vis[vno].fp, vis[vno].src, "Read wrong number of nodal locations");
		    vis[vno].fp = NULL;
		    cxfree(vis[vno].x);
		    wd = NULL;
		    return TRUE;
		}
		break;
	    case 17:
	    case 117:
		vis[vno].xf = (float *) calloc(npts, sizeof(float));
		vis[vno].yf = (float *) calloc(npts, sizeof(float));
		if (fread(vis[vno].xf, sizeof(float), npts, vis[vno].fp) != npts) {
		    vis_error(vis[vno].fp, vis[vno].src, "Read wrong number of nodal locations");
		    vis[vno].fp = NULL;
		    cxfree(vis[vno].x);
		    cxfree(vis[vno].xf);
		    wd = NULL;
		    return TRUE;
		}
		for (i = 0; i < npts; i++) {
		    vis[vno].x[i] = vis[vno].xf[i];
		}
		break;
	    case 19:
	    case 119:
		vis[vno].xf = (float *) calloc(npts, sizeof(float));
		vis[vno].yf = (float *) calloc(npts, sizeof(float));
		break;
	    }
	    vis[vno].y = (double *) calloc(npts, sizeof(double));
	    if (vis[vno].y == NULL) {
		cxfree(vis[vno].x);
		vis_error(vis[vno].fp, vis[vno].src, "Insufficient memory to allocate for y");
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	    softkillset(gno, set1);
	    activateset(gno, set1);
	    settype(gno, set1, XY);

	    setcol(gno, vis[vno].x, set1, npts, 0);
	    setcol(gno, vis[vno].y, set1, npts, 1);

	    setcomment(vis[vno].gno, vis[vno].set1, vis[vno].fname);
	    if (first) {
		first = 0;
	    }
	}
    } else if (mode == 1) {
	vis_error(vis[vno].fp, vis[vno].src, NULL);
	vis[vno].fp = NULL;
	return TRUE;
    } else if (mode == 2) {
	return 0;
    }
    npts = vis[vno].npts;
    switch (vis[vno].ftype) {
    case 18:
    case 118:
	if ((fread(&time, sizeof(double), 1, vis[vno].fp) == 1) && fread(vis[vno].y, sizeof(double), vis[vno].npts, vis[vno].fp) == npts) {
	    updatesetminmax(gno, set1);
	    getsetminmax(gno, set1, &x1, &x2, &y1, &y2);
	    vis[vno].nsteps++;
	    if (!(vis[vno].nsteps % vis[vno].skip)) {
		sprintf(timebuf, "%.0lf s", time);
		my_doublebuffer(TRUE);
	    }
	    return 0;
	} else {
	    vis_error(vis[vno].fp, vis[vno].src, NULL);
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	break;
    case 17:
    case 117:
	if ((fread(&ftime, sizeof(float), 1, vis[vno].fp) == 1) && fread(vis[vno].yf, sizeof(float), npts, vis[vno].fp) == npts) {
	    time = ftime;
	    for (i = 0; i < npts; i++) {
		vis[vno].y[i] = vis[vno].yf[i];
	    }
	    updatesetminmax(gno, set1);
	    getsetminmax(gno, set1, &x1, &x2, &y1, &y2);
	    vis[vno].nsteps++;
	    if (!(vis[vno].nsteps % vis[vno].skip)) {
		sprintf(timebuf, "%.0lf s", time);
		my_doublebuffer(TRUE);
	    }
	    return 0;
	} else {
	    vis_error(vis[vno].fp, vis[vno].src, NULL);
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	break;
    case 19:
    case 119:
	if ((fread(&ftime, sizeof(float), 1, vis[vno].fp) == 1) &&
	    fread(vis[vno].xf, sizeof(float), npts, vis[vno].fp) == npts &&
	    fread(vis[vno].yf, sizeof(float), npts, vis[vno].fp) == npts) {
	    time = ftime;
	    for (i = 0; i < npts; i++) {
		vis[vno].x[i] = vis[vno].xf[i];
		vis[vno].y[i] = vis[vno].yf[i];
	    }
	    updatesetminmax(gno, set1);
	    getsetminmax(gno, set1, &x1, &x2, &y1, &y2);
	    vis[vno].nsteps++;
	    if (!(vis[vno].nsteps % vis[vno].skip)) {
		sprintf(timebuf, "%.0lf s", time);
		my_doublebuffer(TRUE);
	    }
	    return 0;
	} else {
	    vis_error(vis[vno].fp, vis[vno].src, NULL);
	    vis[vno].fp = NULL;
	    wd = NULL;
	    return TRUE;
	}
	break;
    case 20:
    case 120:
	for (j = 0; j < vis[vno].nsets; j++) {
	    if (fread(&npts, sizeof(int), 1, vis[vno].fp) != 1) {
	        vis_error(vis[vno].fp, vis[vno].src, NULL);
	        vis[vno].fp = NULL;
	        wd = NULL;
	        return TRUE;
	    }
	    npts = vis[vno].slen[j];
	    if (fread(vis[vno].xm[j], sizeof(double), npts, vis[vno].fp) == npts &&
		fread(vis[vno].ym[j], sizeof(double), npts, vis[vno].fp) == npts) {
		updatesetminmax(gno, vis[vno].sets[j]);
		getsetminmax(gno, vis[vno].sets[j], &x1, &x2, &y1, &y2);
		vis[vno].nsteps++;
		if (!(vis[vno].nsteps % vis[vno].skip)) {
		    sprintf(timebuf, "%.0lf s", time);
		    my_doublebuffer(TRUE);
		}
	    } else {
		vis_error(vis[vno].fp, vis[vno].src, NULL);
		vis[vno].fp = NULL;
		wd = NULL;
		return TRUE;
	    }
	}
	return 0;
	break;
    }
}
