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
 * input error checking, fexists()
 *
 */

#ifndef lint
static char RCSid[] = "$Id: io.c,v 1.2 2003/07/24 15:44:05 pturner Exp $";
#endif

#include <stdio.h>
#include <pwd.h>
#include <sys/param.h>
#include <sys/types.h>
#include <sys/stat.h>

static char readbuf[80];

int ibounds(int x, int lower, int upper, char *name)
{
    int test;

    test = ((x >= lower) && (x <= upper));
    if (!test) {
	sprintf(readbuf, " in %s : parameter must be in (%d , %d)", name, lower, upper);
	errwin(readbuf);
    }
    return (test);
}

int fbounds(double x, double lower, double upper, char *name)
{
    int test;

    test = ((x >= lower) && (x <= upper));
    if (!test) {
	sprintf(readbuf, "In %s : parameter must be in [%lf, %lf]", name, lower, upper);
	errwin(readbuf);
    }
    return (test);
}

int fexists(char *to)
{
    int fold;
    char tbuf[128];
    struct stat stto;

    fold = open(to, 0);
    if (stat(to, &stto) >= 0) {
	sprintf(tbuf, "%s exists, replace?", to);
	if (!yesno(tbuf, "", "Yes", "No")) {
	    close(fold);
	    return (1);
	}
	close(fold);
	return (0);
    }
    close(fold);
    return (0);
}

int isfile(char *s)
{
    char buf[256];
    struct stat statb;

    /* check to make sure this is a file and not a dir */
    if (stat(s, &statb)) {
	return 0;
    }
    if (!S_ISREG(statb.st_mode)) {
	return 0;
    }
    return 1;
}
