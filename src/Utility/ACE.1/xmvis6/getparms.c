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
 * Read a parameter file
 */

#ifndef lint
static char RCSid[] = "$Id: getparms.c,v 1.2 2003/07/24 15:23:45 pturner Exp $";
#endif

#include <stdio.h>
#include "defines.h"
#include "globals.h"

static char readbuf[1024];

int getparms(int gno, char *plfile)
{
    int linecount = 0, icheck, ptype, errpos = 0, errcnt = 0;
    char s[256];
    FILE *pp;
    double a, b, c, d, x, y;

    sprintf(statusstr, "getparms(): %d %s", gno, plfile);
    writelogfile(statusstr);

    if ((pp = fopen(plfile, "r")) == NULL) {
	sprintf(readbuf, "Can't open parameter file %s", plfile);
	errwin(readbuf);
	return 0;
    } else {
	errcnt = 0;
	while (fgets(readbuf, 511, pp) != NULL) {
	    linecount++;
	    if (readbuf[0] == '#') {
		continue;
	    }
	    if (strlen(readbuf) <= 1) {
		continue;
	    }
	    lowtoupper(readbuf);
	    if (debuglevel == 1) {
		printf(readbuf);
	    }
	    errpos = 0;
	    scanner(readbuf, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
	    if (errpos) {
		printf("Error at line %d: %s\n", linecount, readbuf);
		errcnt++;
		if (errcnt > 5) {
		    if (yesno("Lots of errors, abort?", "Press YES or NO", "YES", "NO")) {
			fclose(pp);
			return 0;
		    } else {
			errcnt = 0;
		    }
		}
	    }
	}
	fclose(pp);
    }
    sprintf(statusstr, "getparms(): completed read");
    writelogfile(statusstr);
    return 1;
}

void read_param(char *pbuf)
{
    int icheck, ptype, errpos = 0;
    double a, b, c, d, x, y;

    if (pbuf[0] == '#') {
	return;
    }
    lowtoupper(pbuf);
    scanner(pbuf, &x, &y, 1, &a, &b, &c, &d, 1, 0, 0, &errpos);
}

int checkptr(void *ptr, char *s)
{
    char buf[4096];
    if (ptr == NULL) {
	sprintf(buf, "Error: NULL pointer, possibly missing 'with' statement. Use 'with grid N' or similar\nCommand: %s\n", s);
	errwin(buf);
	return 0;
    } else {
	return 1;
    }
}
