/* $Id: files.c,v 1.2 2004/05/30 16:04:02 pturner Exp $

 * read data files
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "globals.h"

int readrectgrid(int gno, char *fname, FILE * fp);
int readfegrid(int gno, char *fname, FILE * fp);

#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)

#include "netcdf.h"

#endif

#define MAXERR 50
#define MAX_LINE_LEN 512
/*
 * number of doubles to allocate for each call to realloc
 */
#define BUFSIZE  512

int realtime = 0;
int change_gno;			/* if the graph number changes on read in */
static int cur_gno;		/* if the graph number changes on read in */
int change_type;		/* current set type */
static int cur_type;		/* current set type */

static int readerror = 0;	/* number of errors */
static int readline = 0;	/* line number in file */
static int readfile = 0;	/* number of file read, not used */

int getdata(int gno, char *fn, int src, int type)
{
    FILE *fp = NULL;
    int retval;
#ifdef WIN32
    struct _stat statb;
#else
    struct stat statb;
#endif

    switch (src) {
    case DISK:
	/* check to make sure this is a file and not a dir */
#ifdef WIN32
	/* if (_stat(fn, &statb)) { */
#else
	if (stat(fn, &statb)) {
	    sprintf(buf, "Can't stat file %s", fn);
	    errwin(buf);
	    return 0;
	}
#endif
#ifdef WIN32
	/*if (!(statb.st_mode & _S_IFREG)) { */
#else
	if (!S_ISREG(statb.st_mode)) {
	    sprintf(buf, "File %s is not a regular file: %d %ld", fn, statb.st_mode, statb
		    .st_size);
	    errwin(buf);
	    return 0;
	}
#endif
	fp = fopen(fn, "r");
	readline = 0;
	break;
    case PIPE:
	fp = (FILE *) popen(fn, "r");
	readline = 0;
	break;
    case 2:
	fp = stdin;
	readline = 0;
	break;
    }
    if (fp == NULL) {
	sprintf(buf, "Can't open file %s", fn);
	errwin(buf);
	return 0;
    }
    cur_gno = gno;
    change_type = cur_type = type;
    retval = -1;
    while (retval == -1) {
	retval = 0;
	switch (cur_type) {
	case XY:
	    retval = readxy(cur_gno, fn, fp, 0);
	    break;
	case NXY:
	    retval = readnxy(cur_gno, fn, fp);
	    break;
	case CTD:
	    retval = readCTD(cur_gno, fn, fp);
	    break;
	case IHL:
	    retval = readihl(cur_gno, fn, fp);
	    break;
	case BP:
	    retval = readbp(cur_gno, fn, fp);
	    break;
	case BIN:
	    retval = readbinary(cur_gno, fn, fp);
	    break;
	case XYDX:
	case XYDY:
	case XYDXDX:
	case XYDYDY:
	case XYDXDY:
	case XYZ:
	case XYRT:
	case XYHILO:
	case XYBOXPLOT:
	case XYUV:
	case XYBOX:
	    retval = readxxyy(cur_gno, fn, fp, cur_type);
	    break;
	case XYSTRING:
	    retval = readxystring(cur_gno, fn, fp);
	    break;
	case RECTGRID:
	    retval = readrectgrid(cur_gno, fn, fp);
	    break;
	case FEGRID:
	    retval = readfegrid(cur_gno, fn, fp);
	    break;
	case BLOCK:
	    retval = readblockdata(cur_gno, fn, fp);
	    break;
	}
    }
    if (src == PIPE) {
	pclose(fp);
    } else {
	if (fp != stdin) {	/* leave stdin open */
	    fclose(fp);
	}
    }
    update_status_popup();
    return retval;
}

int getdata_step(int gno, char *fn, int src, int type)
{
    static FILE *fp;
    int retval;

    if (fp == NULL) {
	switch (src) {
	case DISK:
	    fp = fopen(fn, "r");
	    break;
	case PIPE:
	    fp = (FILE *) popen(fn, "r");
	    break;
	case 2:
	    fp = stdin;
	    break;
	case 3:
	    if (fp) {
		if (src == PIPE) {
		    pclose(fp);
		} else {
		    if (fp != stdin) {	/* leave stdin open */
			fclose(fp);
		    }
		}
	    }
	    fp = NULL;
	    return 0;		/* no break */
	}
    }
    if (fp == NULL) {
	sprintf(buf, "Can't open file %s", fn);
	errwin(buf);
	fp = NULL;
	return 0;
    }
    cur_gno = gno;
    change_type = cur_type = type;
    retval = -1;
    while (retval == -1) {
	retval = 0;
	switch (cur_type) {
	case XY:
	    retval = readxy(cur_gno, fn, fp, 1);
	    break;
	case NXY:
	    retval = readnxy(cur_gno, fn, fp);
	    break;
	case CTD:
	    retval = readCTD(cur_gno, fn, fp);
	    break;
	case IHL:
	    retval = readihl(cur_gno, fn, fp);
	    break;
	case BP:
	    retval = readbp(cur_gno, fn, fp);
	    break;
	case BIN:
	    retval = readbinary(cur_gno, fn, fp);
	    break;
	case XYDX:
	case XYDY:
	case XYDXDX:
	case XYDYDY:
	case XYDXDY:
	case XYZ:
	case XYRT:
	case XYHILO:
	case XYBOXPLOT:
	case XYBOX:
	    retval = readxxyy(cur_gno, fn, fp, cur_type);
	    break;
	case XYSTRING:
	    retval = readxystring(cur_gno, fn, fp);
	    break;
	case RECTGRID:
	    retval = readrectgrid(cur_gno, fn, fp);
	    break;
	case FEGRID:
	    retval = readfegrid(cur_gno, fn, fp);
	    break;
	case BLOCK:
	    retval = readblockdata(cur_gno, fn, fp);
	    break;
	}
    }
    if (retval != -2) {		/* means it returned because a single set was
				 * read */
	if (src == PIPE) {
	    pclose(fp);
	} else {
	    if (fp != stdin) {	/* leave stdin open */
		fclose(fp);
	    }
	}
    }
    return retval;
}

/*
 * read file type 0
 */
int readxy(int gno, char *fn, FILE * fp, int readone)
{
    int i = 0, ll, j, pstat, readset = 0, ptype, retval = 0;
    double *x, *y;

    x = (double *) malloc(BUFSIZE * sizeof(double));
    y = (double *) malloc(BUFSIZE * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	return (0);
    }
    while (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	readline++;
	ll = strlen(buf);
	if ((ll > 0) && (buf[ll - 1] != '\n')) {	/* must have a newline
							 * char at end of line */
	    readerror++;
	    fprintf(stderr, "No newline at line #%1d: %s\n", readline, buf);
	    if (readerror > MAXERR) {
		if (yesno("Lots of errors, abort?", NULL, NULL, NULL)) {
		    cxfree(x);
		    cxfree(y);
		    return (0);
		} else {
		    readerror = 0;
		}
	    }
	    continue;
	}
	if (buf[0] == '#') {
	    continue;
	}
/* commented out for 3.01pl2 PJT 12-10-94
 * This fragment would skip blank lines, some users
 * use blank lines for set separators like gnuplot.
 * uncomment for original behavior.
 if (strlen(buf) < 2) {
 continue;
 }
 */
	if (buf[0] == '@') {
	    change_gno = -1;
	    change_type = cur_type;
	    read_param(buf + 1);
	    if (change_gno >= 0) {
		cur_gno = gno = change_gno;
	    }
	    if (change_type != cur_type) {
		cur_type = change_type;
		retval = -1;
		break;		/* exit this module and store any set */
	    }
	    continue;
	}
	convertchar(buf);
	/* count the number of items scanned */
	if ((pstat = sscanf(buf, "%lf %lf", &x[i], &y[i])) >= 1) {
	    /* supply x if missing (y winds up in x) */
	    if (pstat == 1) {
		y[i] = x[i];
		x[i] = i;
	    }
	    if (realtime == 1 && inwin) {
		drawpolysym(&x[i], &y[i], 1, 3, 0, 0, 1.0);
	    }
	    /* got x and y so increment */
	    i++;
	    if (i % BUFSIZE == 0) {
		x = (double *) realloc(x, (i + BUFSIZE) * sizeof(double));
		y = (double *) realloc(y, (i + BUFSIZE) * sizeof(double));
	    }
	} else {
	    if (i != 0) {
		if ((j = nextset(gno)) == -1) {
		    cxfree(x);
		    cxfree(y);
		    return (readset);
		}
		activateset(gno, j);
		settype(gno, j, XY);
		setcol(gno, x, j, i, 0);
		setcol(gno, y, j, i, 1);
		setcomment(gno, j, fn);
		updatesetminmax(gno, j);
		if (realtime == 2 && inwin) {
		    drawsetxy(gno, g[gno].p[j], j);
		}
		readset++;
	    } else {
		readerror++;
		fprintf(stderr, "Error at line #%1d: %s", readline, buf);
		if (readerror > MAXERR) {
		    if (yesno("Lots of errors, abort?", NULL, NULL, NULL)) {
			cxfree(x);
			cxfree(y);
			return (0);
		    } else {
			readerror = 0;
		    }
		}
	    }
	    i = 0;
	    x = (double *) malloc(BUFSIZE * sizeof(double));
	    y = (double *) malloc(BUFSIZE * sizeof(double));
	    if (x == NULL || y == NULL) {
		errwin("Insufficient memory for set");
		cxfree(x);
		cxfree(y);
		return (readset);
	    }
	    if (readone) {
		return (-2);
	    }
	}
    }
    if (i != 0) {
	if ((j = nextset(gno)) == -1) {
	    cxfree(x);
	    cxfree(y);
	    return (readset);
	}
	activateset(gno, j);
	settype(gno, j, XY);
	setcol(gno, x, j, i, 0);
	setcol(gno, y, j, i, 1);
	setcomment(gno, j, fn);
	updatesetminmax(gno, j);
	if (realtime == 2 && inwin) {
	    /*
	     * TODO ??? drawsetxy(g[gno].p[j]);
	     */
	}
	readset++;
    } else {
	cxfree(x);
	cxfree(y);
    }
    if (retval == -1) {
	return retval;
    } else {
	return readset;
    }
}

/*
 * read the first set found in a file to set setno
 */
int read_set_fromfile(int gno, int setno, char *fn, int src)
{
    FILE *fp = NULL;
#ifdef WIN32
    struct _stat statb;
#else
    struct stat statb;
#endif
    int readline;
    int i = 0, j, pstat, readset = 0, ptype, retval = 0;
    double *x, *y;

    switch (src) {
    case DISK:
	/* check to make sure this is a file and not a dir */
#ifdef WIN32
	/* if (_stat(fn, &statb)) { */
#else
	if (stat(fn, &statb)) {
	    sprintf(buf, "Can't stat file %s", fn);
	    errwin(buf);
	    return 0;
	}
#endif
#ifdef WIN32
	/*if (!(statb.st_mode & _S_IFREG)) { */
#else
	if (!S_ISREG(statb.st_mode)) {
	    sprintf(buf, "File %s is not a regular file: %d %ld", fn, statb.st_mode, statb
		    .st_size);
	    errwin(buf);
	    return 0;
	}
#endif
	fp = fopen(fn, "r");
	readline = 0;
	break;
    case PIPE:
	fp = (FILE *) popen(fn, "r");
	readline = 0;
	break;
    case 2:
	fp = stdin;
	readline = 0;
	break;
    }
    if (fp == NULL) {
	sprintf(buf, "Can't open file %s", fn);
	errwin(buf);
	return 0;
    }
    softkillset(gno, setno);
    x = (double *) malloc(BUFSIZE * sizeof(double));
    y = (double *) malloc(BUFSIZE * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	goto breakout;
    }
    while (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	readline++;
	if (buf[strlen(buf) - 1] != '\n') {	/* must have a newline char
						 * at end of line */
	    readerror++;
	    fprintf(stderr, "No newline at line #%1d: %s", readline, buf);
	    continue;
	}
	if (buf[0] == '#') {
	    continue;
	}
	if (buf[0] == '@') {
	    continue;
	}
	convertchar(buf);
	/* count the number of items scanned */
	if ((pstat = sscanf(buf, "%lf %lf", &x[i], &y[i])) >= 1) {
	    /* supply x if missing (y winds up in x) */
	    if (pstat == 1) {
		y[i] = x[i];
		x[i] = i;
	    }
	    i++;
	    if (i % BUFSIZE == 0) {
		x = (double *) realloc(x, (i + BUFSIZE) * sizeof(double));
		y = (double *) realloc(y, (i + BUFSIZE) * sizeof(double));
	    }
	}
    }
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, i, 0);
    setcol(gno, y, setno, i, 1);
    setcomment(gno, setno, fn);
    updatesetminmax(gno, setno);
    retval = 1;

  breakout:;

    if (src == PIPE) {
	pclose(fp);
    } else {
	if (fp != stdin) {	/* leave stdin open */
	    fclose(fp);
	}
    }
    return retval;
}

/*
 * read IHL format
 */
int readihl(int gno, char *fn, FILE * fp)
{
    int i, j, pstat, npts;
    double *x, *y, tmp;

    i = 0;
    pstat = 0;
    if ((j = nextset(gno)) == -1) {
	return 0;
    }
    if (fgets(buf, MAX_LINE_LEN, fp) == NULL) {
	errwin("Can't read from file");
	killset(gno, j);
	return 0;
    }
    readline++;
    pstat = sscanf(buf, "%d", &npts);
    if (npts == 0) {
	errwin("Number of points = 0");
	killset(gno, j);
	return 0;
    }
    activateset(gno, j);
    settype(gno, j, XY);
    setlength(gno, j, npts);
    setcomment(gno, j, fn);
    x = getx(gno, j);
    y = gety(gno, j);
    for (i = 0; i < npts; i++) {
	if (fgets(buf, MAX_LINE_LEN, fp) == NULL) {
	    errwin("Premature EOF");
	    updatesetminmax(gno, j);
	    return 1;
	}
	readline++;
	convertchar(buf);
	pstat = sscanf(buf, "%lf %lf %lf", &tmp, &x[i], &y[i]);
    }
    updatesetminmax(gno, j);
    return 1;
}

/*
 * read BP format
 */
int readbp(int gno, char *fn, FILE * fp)
{
    int i, j, pstat, npts;
    double *x, *y, tmp;

    i = 0;
    pstat = 0;
    if ((j = nextset(gno)) == -1) {
	return 0;
    }
/* skip line */
    if (fgets(buf, MAX_LINE_LEN, fp) == NULL) {
	errwin("Can't read from file");
	killset(gno, j);
	return 0;
    }
    readline++;
    if (fgets(buf, MAX_LINE_LEN, fp) == NULL) {
	errwin("Can't read from file");
	killset(gno, j);
	return 0;
    }
    readline++;
    pstat = sscanf(buf, "%d", &npts);
    if (npts == 0) {
	errwin("Number of points = 0");
	killset(gno, j);
	return 0;
    }
    activateset(gno, j);
    settype(gno, j, XY);
    setlength(gno, j, npts);
    setcomment(gno, j, fn);
    x = getx(gno, j);
    y = gety(gno, j);
    for (i = 0; i < npts; i++) {
	if (fgets(buf, MAX_LINE_LEN, fp) == NULL) {
	    errwin("Premature EOF");
	    updatesetminmax(gno, j);
	    return 1;
	}
	readline++;
	convertchar(buf);
	pstat = sscanf(buf, "%lf %lf %lf", &tmp, &x[i], &y[i]);
    }
    updatesetminmax(gno, j);
    return 1;
}

/*
 * read x1 y1 y2 ... y30 formatted files
 * note that the maximum number of sets is 30
 */
#define MAXSETN 60

int readnxy(int gno, char *fn, FILE * fp)
{
    int i, j, pstat, rcnt, cnt, scnt[MAXSETN], setn[MAXSETN], ptype,
     retval = 0;
    double *x[MAXSETN], *y[MAXSETN], xval, yr[MAXSETN];
    char *s, buf[1024], tmpbuf[1024];
    int readerror = 0;
    int do_restart = 0;

/* if more than one set of nxy data is in the file,
 * leap to here after each is read - the goto is at the
 * bottom of this module.
 */
  restart:;

    i = 0;
    pstat = 0;
    cnt = 0;
    while ((fgets(buf, MAX_LINE_LEN, fp) != NULL) && ((buf[0] == '#') || (buf[0] == '@'))) {
	readline++;
	if (buf[0] == '@') {
	    change_gno = -1;
	    read_param(buf + 1);
	    if (change_gno >= 0) {
		cur_gno = gno = change_gno;
	    }
	}
    }
    convertchar(buf);

    /*
     * count the columns
     */
    strcpy(tmpbuf, buf);
    s = tmpbuf;
    while ((s = strtok(s, " \t\n")) != NULL) {
	cnt++;
	s = NULL;
    }
    if (cnt > maxplot) {
	errwin("Maximum number of columns exceeded, reading first 31");
	cnt = 31;
    }
    s = buf;
    s = strtok(s, " \t\n");
    if (s == NULL) {
	errwin("Read ended by a blank line at or near the beginning of file");
	return 0;
    }
    pstat = sscanf(s, "%lf", &xval);
    if (pstat == 0) {
	errwin("Read ended, non-numeric found on line at or near beginning of file");
	return 0;
    }
    s = NULL;
    for (j = 0; j < cnt - 1; j++) {
	s = strtok(s, " \t\n");
	if (s == NULL) {
	    yr[j] = 0.0;
	    errwin("Number of items in column incorrect");
	} else {
	    yr[j] = atof(s);
	}
	s = NULL;
    }
    if (cnt > 1) {
	for (i = 0; i < cnt - 1; i++) {
	    if ((setn[i] = nextset(gno)) == -1) {
		for (j = 0; j < i; j++) {
		    killset(gno, setn[j]);
		}
		return 0;
	    }
	    activateset(gno, setn[i]);
	    settype(gno, setn[i], XY);
	    x[i] = (double *) malloc(BUFSIZE * sizeof(double));
	    y[i] = (double *) malloc(BUFSIZE * sizeof(double));
	    if (x[i] == NULL || y[i] == NULL) {
		errwin("Insufficient memory for set");
		cxfree(x[i]);
		cxfree(y[i]);
		for (j = 0; j < i + 1; j++) {
		    killset(gno, setn[j]);
		}
		return (0);
	    }
	    *(x[i]) = xval;
	    *(y[i]) = yr[i];
	    scnt[i] = 1;
	}
	while (!do_restart && (fgets(buf, MAX_LINE_LEN, fp) != NULL)) {
	    readline++;
	    if (buf[0] == '#') {
		continue;
	    }
	    if (strlen(buf) < 2) {
		continue;
	    }
	    if (buf[0] == '@') {
		change_gno = -1;
		change_type = cur_type;
		read_param(buf + 1);
		if (change_gno >= 0) {
		    cur_gno = gno = change_gno;
		}
		if (change_type != cur_type) {
		    cur_type = change_type;
		    retval = -1;
		    break;	/* exit this module and store any set */
		}
		continue;
	    }
	    convertchar(buf);
	    s = buf;
	    s = strtok(s, " \t\n");
	    if (s == NULL) {
		continue;
	    }
/* check for set separator */
	    pstat = sscanf(s, "%lf", &xval);
	    if (pstat == 0) {
		do_restart = 1;
		continue;
	    } else {
		s = NULL;
		for (j = 0; j < cnt - 1; j++) {
		    s = strtok(s, " \t\n");
		    if (s == NULL) {
			yr[j] = 0.0;
			errwin("Number of items in column incorrect");
		    } else {
			yr[j] = atof(s);
		    }
		    s = NULL;
		}
		for (i = 0; i < cnt - 1; i++) {
		    *(x[i] + scnt[i]) = xval;
		    *(y[i] + scnt[i]) = yr[i];
		    scnt[i]++;
		    if (scnt[i] % BUFSIZE == 0) {
			x[i] = (double *) realloc(x[i], (scnt[i] + BUFSIZE) * sizeof(double));
			y[i] = (double *) realloc(y[i], (scnt[i] + BUFSIZE) * sizeof(double));
		    }
		}
	    }
	}
	for (i = 0; i < cnt - 1; i++) {
	    setcol(gno, x[i], setn[i], scnt[i], 0);
	    setcol(gno, y[i], setn[i], scnt[i], 1);
	    setcomment(gno, setn[i], fn);
	    updatesetminmax(gno, setn[i]);
	}
	if (!do_restart) {
	    if (retval == -1) {
		return retval;
	    } else {
		return 1;
	    }
	} else {
	    do_restart = 0;
	    goto restart;
	}
    }
    return 0;
}

int readbinary(int gno, char *fn, FILE * fp)
{
    int i, j, type, setn, nsets = 0, npts, n;
    double *x, *y;
    float *xf, *yf;

/*
   fread(&type, sizeof(int), 1, fp);
   fread(&type, sizeof(int), 1, fp);
   if (type != 32) {
   sprintf(buf, "Bad magic: have %d, need 32", type);
   errwin(buf);
   return 0;
   }
 */
    fread(&nsets, sizeof(int), 1, fp);
    if (nsets > g[gno].maxplot) {
	sprintf(buf, "Not enough sets: have %d, need %d", g[gno].maxplot, nsets);
	errwin(buf);
	return 0;
    }
    for (i = 0; i < nsets; i++) {
	fread(&npts, sizeof(int), 1, fp);
	if (npts > 0) {
	    x = (double *) malloc(npts * sizeof(double));
	    if (x == NULL) {
		errwin("Can't malloc in readbinary");
		return 0;
	    }
	    y = (double *) malloc(npts * sizeof(double));
	    if (y == NULL) {
		errwin("Can't malloc in readbinary");
		cxfree(x);
		return 0;
	    }
	    xf = (float *) malloc(npts * sizeof(float));
	    if (xf == NULL) {
		errwin("Can't malloc in readbinary");
		return 0;
	    }
	    yf = (float *) malloc(npts * sizeof(float));
	    if (yf == NULL) {
		errwin("Can't malloc in readbinary");
		cxfree(xf);
		return 0;
	    }
	    fread(xf, sizeof(float), npts, fp);
	    fread(yf, sizeof(float), npts, fp);
	    for (j = 0; j < npts; j++) {
		x[j] = xf[j];
		y[j] = yf[j];
	    }
	    free(xf);
	    free(yf);
	    if ((setn = nextset(gno)) == -1) {
		cxfree(x);
		cxfree(y);
		return 0;
	    }
	    activateset(gno, setn);
	    settype(gno, setn, XY);
	    setcol(gno, x, setn, npts, 0);
	    setcol(gno, y, setn, npts, 1);
	    setcomment(gno, setn, fn);
	    updatesetminmax(gno, setn);
	}
    }
    return 1;
}

int readxystring(void)
{
    return 0;
}

/*
 * read file types using dx and/or dy
 */
int readxxyy(int gno, char *fn, FILE * fp, int type)
{
    int i = 0, j, pstat, readset = 0, ptype, retval = 0;
    double *x, *y, *dx, *dy, *dz, *dw;
    double xtmp, ytmp, dxtmp, dytmp, dztmp, dwtmp;

    x = y = dx = dy = dz = dw = NULL;
    x = (double *) malloc(BUFSIZE * sizeof(double));
    y = (double *) malloc(BUFSIZE * sizeof(double));
    switch (type) {
    case XYZ:
    case XYRT:
    case XYDX:
    case XYDY:
	dx = (double *) malloc(BUFSIZE * sizeof(double));
	break;
    case XYDXDX:
    case XYDYDY:
    case XYDXDY:
    case XYUV:
	dx = (double *) malloc(BUFSIZE * sizeof(double));
	dy = (double *) malloc(BUFSIZE * sizeof(double));
	break;
    case XYHILO:
    case XYBOX:
	dx = (double *) malloc(BUFSIZE * sizeof(double));
	dy = (double *) malloc(BUFSIZE * sizeof(double));
	dz = (double *) malloc(BUFSIZE * sizeof(double));
	break;
    case XYBOXPLOT:
	dx = (double *) malloc(BUFSIZE * sizeof(double));
	dy = (double *) malloc(BUFSIZE * sizeof(double));
	dz = (double *) malloc(BUFSIZE * sizeof(double));
	dw = (double *) malloc(BUFSIZE * sizeof(double));
	break;
    default:
	dx = (double *) malloc(BUFSIZE * sizeof(double));
	dy = (double *) malloc(BUFSIZE * sizeof(double));
	break;
    }
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	cxfree(dx);
	cxfree(dy);
	cxfree(dz);
	cxfree(dw);
	return (0);
    }
    while (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	readline++;
	if (buf[0] == '#') {
	    continue;
	}
	if (strlen(buf) < 2) {
	    continue;
	}
	if (buf[0] == '@') {
	    change_gno = -1;
	    change_type = cur_type;
	    read_param(buf + 1);
	    if (change_gno >= 0) {
		cur_gno = gno = change_gno;
	    }
	    if (change_type != cur_type) {
		if (change_type != cur_type) {
		    cur_type = change_type;
		    retval = -1;
		    break;	/* exit this module and store any set */
		}
	    }
	    continue;
	}
	convertchar(buf);
	/* count the number of items scanned */
	if ((pstat = sscanf(buf, "%lf %lf %lf %lf %lf %lf", &xtmp, &ytmp, &dxtmp, &dytmp, &dztmp, &dwtmp)) >= 1) {
	    /* got x and y so increment */
	    x[i] = xtmp;
	    y[i] = ytmp;
	    if (type == XYDX || type == XYDY || type == XYZ || type == XYRT) {
		dx[i] = dxtmp;
	    } else if (type == XYHILO || type == XYBOX) {
		dx[i] = dxtmp;
		dy[i] = dytmp;
		dz[i] = dztmp;
	    } else if (type == XYBOXPLOT) {
		dx[i] = dxtmp;
		dy[i] = dytmp;
		dz[i] = dztmp;
		dw[i] = dwtmp;
	    } else {
		dx[i] = dxtmp;
		dy[i] = dytmp;
	    }
	    i++;
	    if (i % BUFSIZE == 0) {
		x = (double *) realloc(x, (i + BUFSIZE) * sizeof(double));
		y = (double *) realloc(y, (i + BUFSIZE) * sizeof(double));
		switch (type) {
		case XYDX:
		case XYDY:
		case XYZ:
		case XYRT:
		    dx = (double *) realloc(dx, (i + BUFSIZE) * sizeof(double));
		    break;
		case XYDXDX:
		case XYDYDY:
		case XYDXDY:
		case XYUV:
		    dx = (double *) realloc(dx, (i + BUFSIZE) * sizeof(double));
		    dy = (double *) realloc(dy, (i + BUFSIZE) * sizeof(double));
		    break;
		case XYHILO:
		case XYBOX:
		    dx = (double *) realloc(dx, (i + BUFSIZE) * sizeof(double));
		    dy = (double *) realloc(dy, (i + BUFSIZE) * sizeof(double));
		    dz = (double *) realloc(dz, (i + BUFSIZE) * sizeof(double));
		    break;
		case XYBOXPLOT:
		    dx = (double *) realloc(dx, (i + BUFSIZE) * sizeof(double));
		    dy = (double *) realloc(dy, (i + BUFSIZE) * sizeof(double));
		    dz = (double *) realloc(dz, (i + BUFSIZE) * sizeof(double));
		    dw = (double *) realloc(dz, (i + BUFSIZE) * sizeof(double));
		    break;
		default:
		    dx = (double *) realloc(dx, (i + BUFSIZE) * sizeof(double));
		    dy = (double *) realloc(dy, (i + BUFSIZE) * sizeof(double));
		    break;
		}
	    }
	} else {
	    if (i != 0) {
		if ((j = nextset(gno)) == -1) {
		    cxfree(x);
		    cxfree(y);
		    cxfree(dx);
		    cxfree(dy);
		    cxfree(dz);
		    cxfree(dw);
		    return readset;
		}
		activateset(gno, j);
		settype(gno, j, type);
		setcol(gno, x, j, i, 0);
		setcol(gno, y, j, i, 1);
		setcol(gno, dx, j, i, 2);
		setcol(gno, dy, j, i, 3);
		setcol(gno, dz, j, i, 4);
		setcol(gno, dw, j, i, 5);
		setcomment(gno, j, fn);
		updatesetminmax(gno, j);
		readset++;
	    } else {
		readerror++;
		fprintf(stderr, "Error at line #%1d: %s", readline, buf);
		if (readerror > MAXERR) {
		    if (yesno("Lots of errors, abort?", NULL, NULL, NULL)) {
			cxfree(x);
			cxfree(y);
			cxfree(dx);
			cxfree(dy);
			cxfree(dz);
			cxfree(dw);
			return (0);
		    } else {
			readerror = 0;
		    }
		}
	    }
	    i = 0;
	    x = (double *) malloc(BUFSIZE * sizeof(double));
	    y = (double *) malloc(BUFSIZE * sizeof(double));
	    switch (type) {
	    case XYDX:
	    case XYZ:
	    case XYRT:
	    case XYDY:
		dx = (double *) malloc(BUFSIZE * sizeof(double));
		break;
	    case XYDXDX:
	    case XYDYDY:
	    case XYDXDY:
	    case XYUV:
		dx = (double *) malloc(BUFSIZE * sizeof(double));
		dy = (double *) malloc(BUFSIZE * sizeof(double));
		break;
	    case XYHILO:
	    case XYBOX:
		dx = (double *) malloc(BUFSIZE * sizeof(double));
		dy = (double *) malloc(BUFSIZE * sizeof(double));
		dz = (double *) malloc(BUFSIZE * sizeof(double));
		break;
	    case XYBOXPLOT:
		dx = (double *) malloc(BUFSIZE * sizeof(double));
		dy = (double *) malloc(BUFSIZE * sizeof(double));
		dz = (double *) malloc(BUFSIZE * sizeof(double));
		dw = (double *) malloc(BUFSIZE * sizeof(double));
		break;
	    default:
		dx = (double *) malloc(BUFSIZE * sizeof(double));
		dy = (double *) malloc(BUFSIZE * sizeof(double));
		break;
	    }
	    if (x == NULL || y == NULL) {
		errwin("Insufficient memory for set");
		cxfree(x);
		cxfree(y);
		cxfree(dx);
		cxfree(dy);
		cxfree(dz);
		cxfree(dw);
		killset(gno, j);
		return (readset);
	    }
	}
    }
    if (i != 0) {
	if ((j = nextset(gno)) == -1) {
	    cxfree(x);
	    cxfree(y);
	    cxfree(dx);
	    cxfree(dy);
	    cxfree(dz);
	    cxfree(dw);
	    return readset;
	}
	activateset(gno, j);
	settype(gno, j, type);
	setcol(gno, x, j, i, 0);
	setcol(gno, y, j, i, 1);
	setcol(gno, dx, j, i, 2);
	setcol(gno, dy, j, i, 3);
	setcol(gno, dz, j, i, 4);
	setcol(gno, dw, j, i, 5);
	setcomment(gno, j, fn);
	updatesetminmax(gno, j);
	readset++;
    } else {
	cxfree(x);
	cxfree(y);
	cxfree(dx);
	cxfree(dy);
	cxfree(dz);
	cxfree(dw);
    }
    if (retval == -1) {
	return retval;
    } else {
	return readset;
    }
}

int readrectgrid(int gno, char *fname, FILE * fp)
{
    return 0;
}

int readfegrid(int gno, char *fname, FILE * fp)
{
    int nmnp, nmel, nn, n1, n2, n3, n4, itmp;
    int i = 0, j, retval = 0;
    double *x, *y, *z;
    int *elements[4];
    int *nnodes;

    x = y = z = NULL;
    if (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
/* skip a line of alpha */
    } else {
	return retval;
    }
    if (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	sscanf(buf, "%d %d", &nmel, &nmnp);
    } else {
	return retval;
    }
    x = (double *) malloc(nmnp * sizeof(double));
    if (x == NULL) {
	return retval;
    }
    y = (double *) malloc(nmnp * sizeof(double));
    if (y == NULL) {
	free(x);
	return retval;
    }
    z = (double *) malloc(nmnp * sizeof(double));
    if (z == NULL) {
	free(x);
	free(y);
	return retval;
    }
    for (i = 0; i < nmnp; i++) {
	if (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	    sscanf(buf, "%d %lf %lf %lf", &itmp, &x[i], &y[i], &z[i]);
	} else {
	    free(x);
	    free(y);
	    free(z);
	    return retval;
	}
    }

    nnodes = (int *) malloc(nmel * sizeof(int));
    if (nnodes == NULL) {
    }
    elements[0] = (int *) malloc(nmel * sizeof(int));
    if (elements[0] == NULL) {
    }
    elements[1] = (int *) malloc(nmel * sizeof(int));
    if (elements[1] == NULL) {
    }
    elements[2] = (int *) malloc(nmel * sizeof(int));
    if (elements[2] == NULL) {
    }
    elements[3] = (int *) malloc(nmel * sizeof(int));
    if (elements[3] == NULL) {
    }
    for (i = 0; i < nmel; i++) {
	if (fgets(buf, MAX_LINE_LEN, fp) != NULL) {
	    sscanf(buf, "%d %d %d %d %d %d", &itmp, &nn, &n1, &n2, &n3, &n4);
	} else {
	    free(x);
	    free(y);
	    free(z);
	    return retval;
	}

	nnodes[i] = nn;
	elements[0][i] = n1 - 1;
	elements[1][i] = n2 - 1;
	elements[2][i] = n3 - 1;
	elements[3][i] = n4 - 1;
    }
    if ((j = nextset(gno)) == -1) {
	cxfree(x);
	cxfree(y);
	cxfree(z);
	return retval;
    }
    activateset(gno, j);
    settype(gno, j, FEGRID);
    setlength(gno, j, nmnp);
    setcol(gno, x, j, nmnp, 0);
    setcol(gno, y, j, nmnp, 1);
    setcol(gno, z, j, nmnp, 2);
    setcomment(gno, j, fname);
    updatesetminmax(gno, j);
    g[gno].p[j].nelem = nmel;
    g[gno].p[j].nnodes = nnodes;
    g[gno].p[j].elements[0] = elements[0];
    g[gno].p[j].elements[1] = elements[1];
    g[gno].p[j].elements[2] = elements[2];
    g[gno].p[j].elements[3] = elements[3];

    return 1;
}

double *blockdata[MAXPLOT];
int blocklen;
int blockncols;

/*
 * read block data
 */
int readblockdata(int gno, char *fn, FILE * fp)
{
    int i = 0, j, k, gotcol, ncols, pstat, readset = 0, ptype, retval = 0;
    int first = 1, readerror = 0;
    double *data[MAXPLOT];
    char tmpbuf[1024], *s, tbuf[256];
    int linecount = 0;

    i = 0;
    pstat = 0;
    while ((s = fgets(buf, MAX_LINE_LEN, fp)) != NULL) {
	readline++;
	linecount++;
	if (buf[0] == '#') {
	    continue;
	}
	if (buf[0] == '@') {
	    read_param(buf + 1);
	    continue;
	}
	if ((int) strlen(buf) > 1) {
	    convertchar(buf);
	    if (first) {	/* count the number of columns */
		ncols = 0;
		strcpy(tmpbuf, buf);
		s = tmpbuf;
		while (*s == ' ' || *s == '\t' || *s == '\n')
		    s++;
		while ((s = strtok(s, " \t\n")) != NULL) {
		    ncols++;
		    s = NULL;
		}
		if (ncols < 1 || ncols > MAXPLOT) {
		    errwin("Column count incorrect");
		    return 0;
		}
		for (j = 0; j < MAXPLOT; j++) {
		    cxfree(blockdata[j]);
		    blockdata[j] = (double *) NULL;
		}
		for (j = 0; j < ncols; j++) {
		    data[j] = (double *) malloc(BUFSIZE * sizeof(double));
		    if (data[j] == NULL) {
			errwin("Insufficient memory for block data");
			for (k = 0; k < j; k++) {
			    cxfree(data[k]);
			}
			return 0;
		    }
		}
		first = 0;
	    }
	    s = buf;
	    while (*s == ' ' || *s == '\t' || *s == '\n')
		s++;
	    for (j = 0; j < ncols; j++) {
		s = strtok(s, " \t\n");
		if (s == NULL) {
		    data[j][i] = 0.0;
		    sprintf(tbuf, "Number of items in column incorrect at line %d, line skipped", linecount);
		    errwin(tbuf);
		    readerror++;
		    if (readerror > MAXERR) {
			if (yesno("Lots of errors, abort?", NULL, NULL, NULL)) {
			    for (k = 0; k < ncols; k++) {
				cxfree(data[k]);
			    }
			    return (0);
			} else {
			    readerror = 0;
			}
		    }
		    /* skip the rest */
		    goto bustout;
		} else {
		    data[j][i] = atof(s);
		}
		s = NULL;
	    }
	    i++;
	    if (i % BUFSIZE == 0) {
		for (j = 0; j < ncols; j++) {
		    data[j] = (double *) realloc(data[j], (i + BUFSIZE) * sizeof(double));
		    if (data[j] == NULL) {
			errwin("Insufficient memory for block data");
			for (k = 0; k < j; k++) {
			    cxfree(data[k]);
			}
			return 0;
		    }
		}
	    }
	}
      bustout:;
    }
    for (j = 0; j < ncols; j++) {
	blockdata[j] = data[j];
    }
    blocklen = i;
    blockncols = ncols;
    return 1;
}

void create_set_fromblock(int gno, int type, char *cols)
{
    int i;
    int setno, graphno;
    int d1, cx, cy, c1, c2, c3, c4;
    double *tx, *ty, *t2, *t3, *t4, *t5;
    int nc, coli[MAXPLOT];
    char *s, buf[256];
    strcpy(buf, cols);
    s = buf;
    nc = 0;
    while ((s = strtok(s, ":")) != NULL) {
	coli[nc] = atoi(s);
	coli[nc]--;
	nc++;
	s = NULL;
    }
    if (nc == 0) {
	errwin("No columns scanned in column string");
	return;
    }
    for (i = 0; i < nc; i++) {
	if (coli[i] < 0 || coli[i] >= blockncols) {
	    errwin("Incorrect column specification");
	    return;
	}
    }

    cx = coli[0];
    cy = coli[1];
    if (cx >= blockncols) {
	errwin("Column for X exceeds the number of columns in block data");
	return;
    }
    if (cy >= blockncols) {
	errwin("Column for Y exceeds the number of columns in block data");
	return;
    }
    switch (type) {
    case XY:
	break;
    case XYRT:
    case XYDX:
    case XYDY:
    case XYZ:
	c1 = coli[2];
	if (c1 >= blockncols) {
	    errwin("Column for E1 exceeds the number of columns in block data");
	    return;
	}
	break;
    case XYDXDX:
    case XYDYDY:
    case XYDXDY:
	c1 = coli[2];
	c2 = coli[3];
	if (c1 >= blockncols) {
	    errwin("Column for E1 exceeds the number of columns in block data");
	    return;
	}
	if (c2 >= blockncols) {
	    errwin("Column for E2 exceeds the number of columns in block data");
	    return;
	}
	break;
    case XYHILO:
    case XYBOX:
	c1 = coli[2];
	c2 = coli[3];
	c3 = coli[4];
	if (c1 >= blockncols) {
	    errwin("Column for E1 exceeds the number of columns in block data");
	    return;
	}
	if (c2 >= blockncols) {
	    errwin("Column for E2 exceeds the number of columns in block data");
	    return;
	}
	if (c3 >= blockncols) {
	    errwin("Column for E3 exceeds the number of columns in block data");
	    return;
	}
	break;
    case XYBOXPLOT:
	c1 = coli[2];
	c2 = coli[3];
	c3 = coli[4];
	c4 = coli[5];
	if (c1 >= blockncols) {
	    errwin("Column for E1 exceeds the number of columns in block data");
	    return;
	}
	if (c2 >= blockncols) {
	    errwin("Column for E2 exceeds the number of columns in block data");
	    return;
	}
	if (c3 >= blockncols) {
	    errwin("Column for E3 exceeds the number of columns in block data");
	    return;
	}
	if (c4 >= blockncols) {
	    errwin("Column for E4 exceeds the number of columns in block data");
	    return;
	}
	break;
    }
    setno = -1;
    graphno = -1;

    if (graphno == -1) {
	graphno = cg;
    }
    if (setno == -1) {
	setno = nextset(graphno);
    }
    if (setno == -1) {
	return;
    }
    if (g[graphno].active == OFF) {
	set_graph_active(graphno);
    }
    activateset(graphno, setno);
    settype(graphno, setno, type);

    tx = (double *) malloc(blocklen * sizeof(double));
    ty = (double *) malloc(blocklen * sizeof(double));
    for (i = 0; i < blocklen; i++) {
	tx[i] = blockdata[cx][i];
	ty[i] = blockdata[cy][i];
    }
    setcol(graphno, tx, setno, blocklen, 0);
    setcol(graphno, ty, setno, blocklen, 1);

    switch (type) {
    case XY:
	sprintf(buf, "Cols %d %d", cx + 1, cy + 1);
	break;
    case XYRT:
    case XYDX:
    case XYDY:
    case XYZ:
	sprintf(buf, "Cols %d %d %d", cx + 1, cy + 1, c1 + 1);
	t2 = (double *) malloc(blocklen * sizeof(double));
	for (i = 0; i < blocklen; i++) {
	    t2[i] = blockdata[c1][i];
	}
	setcol(graphno, t2, setno, blocklen, 2);
	break;
    case XYDXDX:
    case XYDYDY:
    case XYDXDY:
	sprintf(buf, "Cols %d %d %d %d", cx + 1, cy + 1, c1 + 1, c2 + 1);
	t2 = (double *) malloc(blocklen * sizeof(double));
	t3 = (double *) malloc(blocklen * sizeof(double));
	for (i = 0; i < blocklen; i++) {
	    t2[i] = blockdata[c1][i];
	    t3[i] = blockdata[c2][i];
	}
	setcol(graphno, t2, setno, blocklen, 2);
	setcol(graphno, t3, setno, blocklen, 3);
	break;
    case XYHILO:
    case XYBOX:
	sprintf(buf, "Cols %d %d %d %d %d", cx + 1, cy + 1, c1 + 1, c2 + 1, c3 + 1);
	t2 = (double *) malloc(blocklen * sizeof(double));
	t3 = (double *) malloc(blocklen * sizeof(double));
	t4 = (double *) malloc(blocklen * sizeof(double));
	for (i = 0; i < blocklen; i++) {
	    t2[i] = blockdata[c1][i];
	    t3[i] = blockdata[c2][i];
	    t4[i] = blockdata[c3][i];
	}
	setcol(graphno, t2, setno, blocklen, 2);
	setcol(graphno, t3, setno, blocklen, 3);
	setcol(graphno, t4, setno, blocklen, 4);
	break;
    case XYBOXPLOT:
	sprintf(buf, "Cols %d %d %d %d %d %d", cx + 1, cy + 1, c1 + 1, c2 + 1, c3 + 1, c4 + 1);
	t2 = (double *) malloc(blocklen * sizeof(double));
	t3 = (double *) malloc(blocklen * sizeof(double));
	t4 = (double *) malloc(blocklen * sizeof(double));
	t5 = (double *) malloc(blocklen * sizeof(double));
	for (i = 0; i < blocklen; i++) {
	    t2[i] = blockdata[c1][i];
	    t3[i] = blockdata[c2][i];
	    t4[i] = blockdata[c3][i];
	    t5[i] = blockdata[c4][i];
	}
	setcol(graphno, t2, setno, blocklen, 2);
	setcol(graphno, t3, setno, blocklen, 3);
	setcol(graphno, t4, setno, blocklen, 4);
	setcol(graphno, t5, setno, blocklen, 5);
	break;
    }

    setcomment(graphno, setno, buf);
    updatesetminmax(graphno, setno);
    update_status_popup();
    drawgraph();
}

#if defined(HAVE_NETCDF) || defined(HAVE_MFHDF)

void check_netcdf_err(const int stat, const int line, const char *file)
{
    if (stat != NC_NOERR) {
	(void) fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    }
}


/*
 * read a variable from netcdf file into a set in graph gno
 * xvar and yvar are the names for x, y in the netcdf file resp.
 * return 0 on fail, return 1 if success.
 *
 * if xvar == NULL, then load the index of the point to x
 *
 */
int readnetcdf(int gno,
	       int setno,
	       char *netcdfname,
	       char *xvar,
	       char *yvar,
	       int nstart,
	       int nstop,
	       int nstride,
	       int index2d)
{
    int cdfid;			/* netCDF id */
    int ndims, nvars, ngatts, recdim;
    int err;
    int i, n, retval = 0;
    double *x, *y;
    float *xf, *yf;
    short *xs, *ys;
    unsigned char *xc, *yc;
    int *xl, *yl;

    /* dimension id for unlimited dimension */
    int udim;

    /* variable ids */
    int x_id, y_id;

    /* variable shapes */
    int dims[2];
    size_t start[2];
    size_t count[2];
    size_t stride[2];

    nc_type xdatatype;
    nc_type ydatatype;
    int xndims, xdim[10], xnatts;
    int yndims, ydim[10], ynatts;
    size_t nx, ny;

    long size;
    char name[256];

    int ncstatus;
/*
   printf("readnetcdf(int gno = %d, int setno = %d, char *netcdfname = %s, char *xvar = %s, char *yvar = %s, int nstart = %d, int nstop = %d, int nstride = %d, int index2d = %d)\n", gno, setno, netcdfname, xvar, yvar, nstart, nstop, nstride, index2d);
 */

/*
 * get a set if on entry setno == -1, if setno=-1, then fail
 */
    if (setno == -1) {
	if ((setno = nextset(gno)) == -1) {
	    return 0;
	}
    } else {
	if (isactive(cg, setno)) {
	    killset(gno, setno);
	}
    }
/*
 * open the netcdf file and locate the variable to read
 */
    if ((ncstatus = nc__open(netcdfname, NC_NOWRITE, NC_SIZEHINT_DEFAULT, &cdfid)) != 0) {
	errwin("Can't open file.");
	return 0;
    }
    if (xvar != NULL) {
	if ((ncstatus = nc_inq_varid(cdfid, xvar, &x_id)) != NC_NOERR) {
	    char ebuf[256];
	    sprintf(ebuf, "readnetcdf(): No such variable %s for X", xvar);
	    errwin(ebuf);
	    return 0;
	}
	nc_inq_var(cdfid, x_id, NULL, &xdatatype, &xndims, xdim, &xnatts);
	nc_inq_dim(cdfid, xdim[0], NULL, &nx);
	if (xndims != 1) {
	    errwin("Number of dimensions for X must be 1.");
	    return 0;
	}
    }
    if ((ncstatus = nc_inq_varid(cdfid, yvar, &y_id)) != NC_NOERR) {
	char ebuf[256];
	sprintf(ebuf, "readnetcdf(): No such variable %s for Y", yvar);
	errwin(ebuf);
	return 0;
    }
    nc_inq_var(cdfid, y_id, NULL, &ydatatype, &yndims, ydim, &ynatts);
    nc_inq_dim(cdfid, ydim[0], NULL, &ny);
    if (yndims > 2) {
	errwin("Number of dimensions for Y must be 1 or 2.");
	return 0;
    }
    if (xvar != NULL) {
	n = nx < ny ? nx : ny;
    } else {
	n = ny;
    }
    if (n <= 0) {
	errwin("Length of dimension == 0.");
	return 0;
    }
    start[0] = nstart;
    count[0] = n;
    if (nstart >= 0 && nstop > 0 && nstop > nstart && (nstop - nstart + 1) < n) {
	n = nstop - nstart + 1;
	start[0] = nstart;
	count[0] = n;
    } else {
	start[0] = 0;
	count[0] = n;
    }
/*
 * allocate for this set
 */
    x = (double *) malloc(n * sizeof(double));
    y = (double *) malloc(n * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	nc_close(cdfid);
	return 0;
    }
/*
 * read the variables from the netcdf file
 */
    if (xvar != NULL) {
/* TODO should check for other data types here */
/* TODO should check for NULL on the mallocs() */
/* TODO making assumptions about the sizes of shorts and longs */
	switch (xdatatype) {
	case NC_CHAR:
	    xc = (char *) malloc(sizeof(char) * n);
	    ncstatus = nc_get_vara_text(cdfid, x_id, start, count, (void *) xc);
	    if (ncstatus != NC_NOERR) {
		fprintf(stderr, "Error reading unsigned char in x\n");
	    }
	    for (i = 0; i < n; i++) {
		x[i] = xc[i];
	    }
	    free(xc);
	    break;
	case NC_BYTE:
	    xc = (char *) malloc(sizeof(char) * n);
	    ncstatus = nc_get_vara_uchar(cdfid, x_id, start, count, (void *) xc);
	    if (ncstatus != NC_NOERR) {
		fprintf(stderr, "Error reading unsigned char in x\n");
	    }
	    for (i = 0; i < n; i++) {
		x[i] = xc[i];
	    }
	    free(xc);
	    break;
	case NC_SHORT:
	    xs = (short *) malloc(n * sizeof(short));
	    nc_get_vara_short(cdfid, x_id, start, count, (void *) xs);
	    for (i = 0; i < n; i++) {
		x[i] = xs[i];
	    }
	    free(xs);
	    break;
	case NC_INT:
	    xl = (int *) malloc(n * sizeof(int));
	    nc_get_vara_long(cdfid, x_id, start, count, (void *) xl);
	    for (i = 0; i < n; i++) {
		x[i] = xl[i];
	    }
	    free(xl);
	    break;
	case NC_FLOAT:
	    xf = (float *) malloc(n * sizeof(float));
	    nc_get_vara_float(cdfid, x_id, start, count, (void *) xf);
	    for (i = 0; i < n; i++) {
		x[i] = xf[i];
	    }
	    free(xf);
	    break;
	case NC_DOUBLE:
	    nc_get_vara_double(cdfid, x_id, start, count, (void *) x);
	    break;
	default:
	    errwin("Data type not supported");
	    cxfree(x);
	    cxfree(y);
	    nc_close(cdfid);
	    return 0;
	}
    } else {			/* just load index */
	for (i = 0; i < n; i++) {
	    x[i] = i + 1;
	}
    }
    switch (ydatatype) {
    case NC_CHAR:
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	yc = (unsigned char *) malloc(sizeof(unsigned char) * n);
	ncstatus = nc_get_vara_text(cdfid, y_id, start, count, yc);
	check_netcdf_err(ncstatus, __LINE__, __FILE__);
	for (i = 0; i < n; i++) {
	    y[i] = yc[i];
	}
	free(yc);
	break;
    case NC_BYTE:
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	yc = (unsigned char *) malloc(sizeof(unsigned char) * n);
	ncstatus = nc_get_vara_uchar(cdfid, y_id, start, count, yc);
	check_netcdf_err(ncstatus, __LINE__, __FILE__);
	for (i = 0; i < n; i++) {
	    y[i] = yc[i];
	}
	free(yc);
	break;
    case NC_SHORT:
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	ys = (short *) malloc(n * sizeof(short));
	ncstatus = nc_get_vara_short(cdfid, y_id, start, count, (void *) ys);
	check_netcdf_err(ncstatus, __LINE__, __FILE__);
	for (i = 0; i < n; i++) {
	    y[i] = ys[i];
	}
	free(ys);
	break;
    case NC_INT:
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	yl = (int *) malloc(n * sizeof(int));
	ncstatus = nc_get_vara_long(cdfid, y_id, start, count, (void *) yl);
	check_netcdf_err(ncstatus, __LINE__, __FILE__);
	for (i = 0; i < n; i++) {
	    y[i] = yl[i];
	}
	free(yl);
	break;
    case NC_FLOAT:
/* TODO should check for NULL here */
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	yf = (float *) malloc(n * sizeof(float));
	nc_get_vara_float(cdfid, y_id, start, count, (void *) yf);
	for (i = 0; i < n; i++) {
	    y[i] = yf[i];
	}
	free(yf);
	break;
    case NC_DOUBLE:
	if (yndims == 2) {
	    start[1] = index2d;
	    count[1] = 1;
	}
	nc_get_vara_double(cdfid, y_id, start, count, (void *) y);
	break;
    default:
	errwin("Data type not supported");
	cxfree(x);
	cxfree(y);
	return 0;
    }
    nc_close(cdfid);

/*
 * initialize stuff for the newly created set
 */
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, n, 0);
    setcol(gno, y, setno, n, 1);

    sprintf(buf, "File %s x = %s y = %s", netcdfname, xvar == NULL ? "Index" : xvar, yvar);
    setcomment(gno, setno, buf);
    updatesetminmax(gno, setno);
    return 1;
}

/*
 * read a 2d variable from netcdf file with each column going to
 * a different set.
 */
int netcdfreadADP(int gno, int setno,
		  char *netcdfname,
		  char *varname,
		  int nstep)
{
    int cdfid;			/* netCDF id */
    int ndims, nvars, ngatts, recdim;
    int err;
    int i, n, retval = 0;
    short ncells;
    double *x, *y;
    float *xf, *yf;
    short *xs, *ys;
    char *xc, *yc;
    int *xl, *yl;

    /* dimension id for unlimited dimension */
    int udim;

    /* variable ids */
    int varid, timeid, ncellsid;

    /* variable shapes */
    int dims[2];
    size_t start[2];
    size_t count[2];
    size_t stride[2];

    nc_type datatype;
    int dim[10], natts;
    size_t nx, ny;

    size_t size;
    char name[256];

    int ncstatus;

/*
 * open the netcdf file and locate the variable to read
 */
    if ((ncstatus = nc__open(netcdfname, NC_NOWRITE, NC_SIZEHINT_DEFAULT, &cdfid)) != 0) {
	errwin("Can't open file.");
	return 0;
    }
/*
 * get a set if on entry setno == -1, if setno=-1, then fail
 */
    if (setno == -1) {
	if ((setno = nextset(gno)) == -1) {
	    return 0;
	}
    } else {
	if (isactive(cg, setno)) {
	    killset(gno, setno);
	}
    }
    if ((ncstatus = nc_inq_varid(cdfid, "time", &timeid)) != NC_NOERR) {
	char ebuf[256];
	sprintf(ebuf, "Error reading time varid");
	errwin(ebuf);
	return 0;
    }
    if ((ncstatus = nc_inq_varid(cdfid, "ncells", &ncellsid)) != NC_NOERR) {
	char ebuf[256];
	sprintf(ebuf, "Error reading ncells varid");
	errwin(ebuf);
	return 0;
    }
    if (varname != NULL) {
	if ((ncstatus = nc_inq_varid(cdfid, varname, &varid)) != NC_NOERR) {
	    char ebuf[256];
	    sprintf(ebuf, "readnetcdf(): No such variable %s", varname);
	    errwin(ebuf);
	    return 0;
	}
	nc_inq_var(cdfid, varid, NULL, &datatype, &ndims, dim, &natts);
	nc_inq_dim(cdfid, dim[0], NULL, &ny);
	if (ndims != 2) {
	    errwin("Number of dimensions for variable must be 2.");
	    return 0;
	}
    }
    if (ny <= 0) {
	errwin("Length of time dimension == 0.");
	return 0;
    }
    start[0] = nstep;
    count[0] = 1;
    ncstatus = nc_get_var1_short(cdfid, ncellsid, (size_t *) start, (void *) &ncells);
    start[1] = 0;
    count[1] = ncells;
    n = ncells;
/*
 * allocate for this set
 */
    x = (double *) malloc(ncells * sizeof(double));
    y = (double *) malloc(ncells * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	nc_close(cdfid);
	return 0;
    }
    for (i = 0; i < ncells; i++) {
	x[i] = i + 1;
    }
/*
 * read the variables from the netcdf file
 */
/* TODO should check for other data types here */
/* TODO should check for NULL on the mallocs() */
/* TODO making assumptions about the sizes of shorts and longs */
    switch (datatype) {
    case NC_CHAR:
	xc = (char *) malloc(sizeof(char) * n);
	nc_get_vara_text(cdfid, varid, start, count, (void *) xc);
	for (i = 0; i < n; i++) {
	    y[i] = xc[i];
	}
	free(xc);
	break;
    case NC_BYTE:
	xc = (char *) malloc(sizeof(char) * n);
	nc_get_vara_uchar(cdfid, varid, start, count, (void *) xc);
	for (i = 0; i < n; i++) {
	    y[i] = xc[i];
	}
	free(xc);
	break;
    case NC_SHORT:
	xs = (short *) malloc(n * sizeof(short));
	nc_get_vara_short(cdfid, varid, start, count, (void *) xs);
	for (i = 0; i < n; i++) {
	    y[i] = xs[i];
	}
	free(xs);
	break;
    case NC_INT:
	xl = (int *) malloc(n * sizeof(int));
	nc_get_vara_int(cdfid, varid, start, count, (void *) xl);
	for (i = 0; i < n; i++) {
	    y[i] = xl[i];
	}
	free(xl);
	break;
    case NC_FLOAT:
	xf = (float *) malloc(n * sizeof(float));
	nc_get_vara_float(cdfid, varid, start, count, (void *) xf);
	for (i = 0; i < n; i++) {
	    y[i] = xf[i];
	}
	free(xf);
	break;
    case NC_DOUBLE:
	nc_get_vara_double(cdfid, varid, start, count, (void *) x);
	break;
    default:
	errwin("Data type not supported");
	cxfree(x);
	cxfree(y);
	nc_close(cdfid);
	return 0;
    }
    nc_close(cdfid);

/*
 * initialize stuff for the newly created set
 */
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, ncells, 0);
    setcol(gno, y, setno, ncells, 1);

    sprintf(buf, "File %s x = cell number y = %s", netcdfname, varname);
    setcomment(gno, setno, buf);
    updatesetminmax(gno, setno);
    return 1;
}

/*
 * write out sets in netcdf format
 */
void writesets_netcdf(char *fn)
{
    char buf[1024];
    char varname[256];
    char dim1name[256];
    char dim2name[256];
    int i, j, k, n, ncols;
    int *gno_id, *setno_id;
    int dimids[2];
    size_t count[2], start[2];
    int ncid, ncstatus, varid;
    double *x, *y, *dx, *dy, *dz, *dw;

    if (!fn[0]) {
	errwin("Define file name first");
	return;
    }
    if (fexists(fn)) {
	return;
    }
    ncstatus = nc_create(fn, NC_WRITE, &ncid);
    if (ncstatus != 0) {
	errwin("Error creating netcdf file");
	return;
    }
    for (k = 0; k < maxgraph; k++) {
	if (isactive_graph(k)) {
	    for (j = 0; j < g[k].maxplot; j++) {
		if (isactive(k, j)) {
		    x = getx(k, j);
		    y = gety(k, j);
		    n = getsetlength(k, j);
		    ncols = getncols(k, j);
		    sprintf(varname, "g%03d_s%05d", k, j);
		    sprintf(dim1name, "g%03d_s%05d_dim1", k, j);
		    sprintf(dim2name, "g%03d_s%05d_dim2", k, j);
		    ncstatus = nc_def_dim(ncid, dim1name, n, &dimids[0]);
		    ncstatus = nc_def_dim(ncid, dim2name, ncols, &dimids[1]);
		    ncstatus = nc_def_var(ncid, varname, NC_DOUBLE, 2, dimids, &varid);
		}
	    }
	}
    }
    ncstatus = nc_enddef(ncid);
    for (k = 0; k < maxgraph; k++) {
	if (isactive_graph(k)) {
	    for (j = 0; j < g[k].maxplot; j++) {
		if (isactive(k, j)) {
		    sprintf(varname, "g%03d_s%05d", k, j);
		    n = getsetlength(k, j);
		    ncols = getncols(k, j);
		    ncstatus = nc_inq_varid(ncid, varname, &varid);

		    start[0] = 0;
		    start[1] = 0;
		    count[0] = n;
		    count[1] = 1;
		    for (i = 0; i < ncols; i++) {
			x = getcol(k, j, i);
			start[1] = i;
			nc_put_vara_double(ncid, varid, start, count, x);
		    }
		}
	    }
	}
    }
    ncstatus = nc_close(ncid);
}

#endif				/* HAVE_NETCDF */

#ifdef HAVE_HDF
/* the netCDF interface is OK for now, no need to implement this */
#endif				/* HAVE_HDF */

/*
 * get a line from a file ended by a newline character
 */
int my_getline(char *buf, int maxlen, int fd)
{
    int n, count = 0;
    char cc;
    while ((n = read(fd, &cc, 1)) > 0) {
	if (cc == '\n') {
	    buf[count++] = cc;
	    buf[count] = 0;
	    return count;
	} else {
	    buf[count++] = cc;
	}
    }
    return -1;
}

int readCTD(int gno, char *fn, FILE * fp)
{
    int i;
    double year, ndays, day, hhmm, h, m, seconds, c, t, s, p;
    int j, pstat, rcnt, cnt, scnt[MAXSETN], setn[MAXSETN], ptype, retval = 0;
    double *x[MAXSETN], *y[MAXSETN], xval, yr[MAXSETN];
    char buf[1024], tmpbuf[1024];
    int readerror = 0;
    cnt = 4;
    for (i = 0; i < cnt; i++) {
	if ((setn[i] = nextset(i)) == -1) {
	    for (j = 0; j < i; j++) {
		killset(i, setn[j]);
	    }
	    return 0;
	}
	activateset(i, setn[i]);
	settype(i, setn[i], XY);
	x[i] = (double *) malloc(BUFSIZE * sizeof(double));
	y[i] = (double *) malloc(BUFSIZE * sizeof(double));
	if (x[i] == NULL || y[i] == NULL) {
	    errwin("Insufficient memory for set");
	    cxfree(x[i]);
	    cxfree(y[i]);
	    for (j = 0; j < i + 1; j++) {
		killset(j, setn[j]);
	    }
	    return (0);
	}
	scnt[i] = 0;
    }
    define_arrange(4, 1, 2, 0.0, 0.0, 0.2, 0.1, 0.7, 0.2, 0);

    while ((fgets(buf, MAX_LINE_LEN, fp) != NULL)) {
	readline++;
	if (buf[0] == '#') {
	    continue;
	}
	if (strlen(buf) < 2) {
	    continue;
	}
	convertchar(buf);
	sscanf(buf, "%lf %lf %lf %lf %lf %lf %lf %lf",
	       &year, &day, &hhmm, &seconds, &c, &s, &t, &p);
	h = (int) (hhmm / 100.0);
	m = hhmm - h * 100.0;
	day = day + h / 24.0 + m / (24.0 * 60.0) + seconds / (24.0 * 3600.0);
	ndays = 0.0;
	for (i = 1996; i < (int) year; i++) {
	    ndays = ndays + ((i % 4) ? 365.0 : 366.0);
	}
	yr[0] = c;
	yr[1] = s;
	yr[2] = t;
	yr[3] = p;
	for (i = 0; i < cnt; i++) {
	    *(x[i] + scnt[i]) = ndays + day;
	    *(y[i] + scnt[i]) = yr[i];
	    scnt[i]++;
	    if (scnt[i] % BUFSIZE == 0) {
		x[i] = (double *) realloc(x[i], (scnt[i] + BUFSIZE) * sizeof(double));
		y[i] = (double *) realloc(y[i], (scnt[i] + BUFSIZE) * sizeof(double));
	    }
	}
    }
    for (i = 0; i < cnt; i++) {
	setcol(i, x[i], setn[i], scnt[i], 0);
	setcol(i, y[i], setn[i], scnt[i], 1);
	setcomment(i, setn[i], fn);
	updatesetminmax(i, setn[i]);
    }
    return 0;
}
