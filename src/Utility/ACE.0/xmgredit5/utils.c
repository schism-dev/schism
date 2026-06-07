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
 * misc utilities
 *
 */

#ifndef lint
static char RCSid[] = "$Id: utils.c,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

double hypot(double x, double y);

/*
 * compute the mins and maxes of a vector x
 */
void minmax(double *x, int n, double *xmin, double *xmax, int *imin, int *imax)
{
    int i;

    *xmin = x[0];
    *xmax = x[0];
    *imin = 1;
    *imax = 1;
    for (i = 1; i < n; i++) {
	if (x[i] < *xmin) {
	    *xmin = x[i];
	    *imin = i + 1;
	}
	if (x[i] > *xmax) {
	    *xmax = x[i];
	    *imax = i + 1;
	}
    }
}

/*
 * free and check for NULL pointer
 */
void cxfree(void *ptr)
{
    if (ptr != NULL) {
	free(ptr);
    }
}

/*
 * swap doubles and ints
 */
void fswap(double *x, double *y)
{
    double tmp;

    tmp = (*x);
    *x = (*y);
    *y = tmp;
}

void iswap(int *x, int *y)
{
    int tmp;

    tmp = (*x);
    *x = (*y);
    *y = tmp;
}

void cswap(char *x, char *y)
{
    char tmp;

    tmp = *x;
    *x = *y;
    *y = tmp;
}

/*
 * convert a string from lower to upper case
 */
void lowtoupper(char *s)
{
    int i;

    for (i = 0; i < strlen(s); i++) {
	if (s[i] >= 'a' && s[i] <= 'z')
	    s[i] = s[i] - ' ';
    }
}

/*
 * remove all that fortran nastiness
 */
void convertchar(char *s)
{
    while (*s++) {
	if (*s == ',')
	    *s = ' ';
	if (*s == 'D' || *s == 'd')
	    *s = 'e';
    }
}

/*
 * log base 2
 */
int ilog2(int n)
{
    int i = 0;
    int n1 = n;

    while (n1 >>= 1)
	i++;
    if (1 << i != n)
	return -1;
    else
	return i;
}

/*
 * compute the area bounded by the polygon (xi,yi)
 */
double comp_area(int n, double *x, double *y)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < n; i++) {
	sum = sum + x[i] * y[(i + 1) % n] - y[i] * x[(i + 1) % n];
    }
    return sum * 0.5;
}

/*
 * compute the perimeter bounded by the polygon (xi,yi)
 */
double comp_perimeter(int n, double *x, double *y)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < n - 1; i++) {
	sum = sum + hypot(x[i] - x[(i + 1) % n], y[i] - y[(i + 1) % n]);
    }
    return sum;
}

int argmatch(char *s1, char *s2, int atleast)
{
    int l1 = strlen(s1);
    int l2 = strlen(s2);

    if (l1 < atleast) {
	return 0;
    }
    if (l1 > l2) {
	return 0;
    }
    return (strncmp(s1, s2, l1) == 0);
}

double uniform(void)
{
    double drand48(void);
    return drand48();
}

double expt(double a, register int n), nicenum(double x, int round);

/*
 * nicenum: find a "nice" number approximately equal to x
 * round if round=1, ceil if round=0
 */

double nicenum(double x, int round)
{
    int exp;
    double f, y;

    exp = floor(log10(x));
    f = x / expt(10., exp);	/* fraction between 1 and 10 */
    if (round)
	if (f < 1.5)
	    y = 1.;
	else if (f < 3.)
	    y = 2.;
	else if (f < 7.)
	    y = 5.;
	else
	    y = 10.;
    else if (f <= 1.)
	y = 1.;
    else if (f <= 2.)
	y = 2.;
    else if (f <= 5.)
	y = 5.;
    else
	y = 10.;
    return y * expt(10., exp);
}

/*
 * expt(a,n)=a^n for integer n
 * roundoff errors in pow were causing problems, so I wrote my own
 */

double expt(double a, register int n)
{
    double x;

    x = 1.;
    if (n > 0)
	for (; n > 0; n--)
	    x *= a;
    else
	for (; n < 0; n++)
	    x /= a;
    return x;
}

extern int swapBytes;

void revbytes(unsigned char *ptr, int len)
{
    int i, iadj = len - 1;
    unsigned char ctmp;
    for (i = 0; i < len / 2; i++) {
	ctmp = ptr[i];
	ptr[i] = ptr[iadj - i];
	ptr[iadj - i] = ctmp;
    }
}

int read_double(double *d, int n, FILE * fp)
{
    int err, type, i;
    err = fread(d, sizeof(double), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(double));
	}
    }
    return err == n ? 0 : err;
}

int write_double(double *d, int n, FILE * fp)
{
    int err, type, i;
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(double));
	}
    }
    err = fwrite(d, sizeof(double), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(double));
	}
    }
    return err == n ? 0 : err;
}

int read_int(int *d, int n, FILE * fp)
{
    int err, type, i;
    err = fread(d, sizeof(int), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(int));
	}
    }
    return err == n ? 0 : err;
}

int write_int(int *d, int n, FILE * fp)
{
    int err, type, i;
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(int));
	}
    }
    err = fwrite(d, sizeof(int), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(int));
	}
    }
    return err == n ? 0 : err;
}

int read_char(char *d, int n, FILE * fp)
{
    int err, type;
    err = fread(d, sizeof(char), n, fp);
    return err == n ? 0 : err;
}

int write_char(char *d, int n, FILE * fp)
{
    int err, type;
    err = fwrite(d, sizeof(char), n, fp);
    return err == n ? 0 : err;
}

int read_short(short *d, int n, FILE * fp)
{
    int err, type, i;
    err = fread(d, sizeof(short), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(short));
	}
    }
    return err == n ? 0 : err;
}

int write_short(short *d, int n, FILE * fp)
{
    int err, type, i;
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(short));
	}
    }
    err = fwrite(d, sizeof(short), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(short));
	}
    }
    return err == n ? 0 : err;
}

int read_float(float *d, int n, FILE * fp)
{
    int err, type, i;
    err = fread(d, sizeof(float), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(float));
	}
    }
    return err == n ? 0 : err;
}

int write_float(float *d, int n, FILE * fp)
{
    int err, type, i;
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(float));
	}
    }
    err = fwrite(d, sizeof(float), n, fp);
    if (swapBytes) {
	for (i = 0; i < n; i++) {
	    revbytes((unsigned char *) &d[i], sizeof(float));
	}
    }
    return err == n ? 0 : err;
}
