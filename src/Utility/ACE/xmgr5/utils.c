/* $Id: utils.c,v 1.2 2004/06/02 17:33:42 pturner Exp $
 *
 * misc utilities
 *
 * Contents:
 *
 * void cxfree() - cfree and check for NULL pointer
 * void fswap()  - swap doubles
 * void iswap()  - swap ints
 * void lowtoupper() - convert a string to upper case
 * void convertchar() - remove commas and Fortran D format
 * int ilog2() - integer log base 2, for the fft routine
 * double comp_area() - compute the area of a polygon
 * double comp_perimeter() - compute the perimeter
 * double fmin(), fmax()
 * Julian date routines
 *
 */

#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

/*
 * cfree and check for NULL pointer
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

    tmp = *x;
    *x = *y;
    *y = tmp;
}

void iswap(int *x, int *y)
{
    int tmp;

    tmp = *x;
    *x = *y;
    *y = tmp;
}

int isoneof(int c, char *s)
{
    while (*s) {
	if (c == *s) {
	    return 1;
	} else {
	    s++;
	}
    }
    return 0;
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

/*
 * convert a string from lower to upper case
 * leaving quoted strings alone
 */
void lowtoupper(char *s)
{
    int i, quoteon = 0;

    for (i = 0; i < strlen(s); i++) {
	if (s[i] == '"') {
	    if (!quoteon) {
		quoteon = 1;
	    } else {
		quoteon = 0;
	    }
	}
	if (s[i] >= 'a' && s[i] <= 'z' && !quoteon) {
	    s[i] -= ' ';
	}
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

double my_hypot(double x, double y)
{
    return sqrt(x * x + y * y);
}

/*
 * compute the perimeter bounded by the polygon (xi,yi)
 */
double comp_perimeter(int n, double *x, double *y)
{
    int i;
    double sum = 0.0;

    for (i = 0; i < n - 1; i++) {
	sum = sum + my_hypot(x[i] - x[(i + 1) % n], y[i] - y[(i + 1) % n]);
    }
    return sum;
}

/* should be macros */
double fmin(double x, double y)
{
    return (x < y ? x : y);
}

double fmax(double x, double y)
{
    return (x > y ? x : y);
}

/*
 * Time and date routines
 */

char *dayofweekstrs[] =
    { "Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat" };
char *dayofweekstrl[] =
    { "Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
"Saturday" };
char *months[] =
    { "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct",
"Nov", "Dec" };
char *monthl[] = { "January", "February", "March", "April", "May", "June",
    "July", "August", "September", "October", "November", "December"
};

static int days1[] =
    { 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365 };
static int days2[] =
    { 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366 };

double julday(int mon, int day, int year, int h, int mi, double se);

int corietime = 0;
double coriestart = 0;

void SetCorieTime(int b)
{
    corietime = b;
    if (b) {
/*
	putenv("TZ=PST");
	julday(12, 31, 1995, 0, 0, 0);
*/
	coriestart = 9495.0;
    } else {
	coriestart = 0.0;
    }
}

void SetStartTime(int b, char *t)
{
    int y, m, d, h, mi, s;
    corietime = b;
    if (b) {
	sscanf(t, "%2d%2d%2d%2d%2d%2d", &y, &m, &d, &h, &mi, &s);
	coriestart = julday(m, d, y, h, mi, s);
    } else {
	coriestart = 0.0;
    }
}

int GetCorieTime(void)
{
    return corietime;
}

/* printf("%12.8lf\n", jd = julday(1, 2, 1996, 0, 0, 0)); */

/*
 *  * return the Julian day + hms as a real number
 *   
*/
double julday(int mon, int day, int year, int h, int mi, double se)
{
    double frac = (int) se - se;
    struct tm t;
    time_t j;
    t.tm_sec = se;
    t.tm_min = mi;
    t.tm_hour = h;
    t.tm_mday = day;
    t.tm_mon = mon - 1;
    t.tm_year = year - 1900;
    j = mktime(&t);
    return ((double) j / 86400.0 + coriestart);
}

/*
 * 
*/
void calcdate(double jd, int *m, int *d, int *y, int *h, int *mi,
	      double *sec)
{
    struct tm *t;
    time_t seconds = (jd + coriestart) * 86400.0;
    t = gmtime(&seconds);
    *sec = t->tm_sec;
    *mi = t->tm_min;
    *h = t->tm_hour;
    *d = t->tm_mday;
    *m = t->tm_mon + 1;
    *y = t->tm_year + 1900;
}

int dayofweek(double j)
{
    return (int) (j + 1) % 7;
}

int leapyear(int year)
{
    if (year % 4 == 0) {
	return (1);
    } else {
	return (0);
    }
}

/*
   get the month and day given the number of days
   from the beginning of the year 'yr'
*/
void getmoday(int days, int yr, int *mo, int *da)
{
    int i;

    if (leapyear(yr)) {
	for (i = 0; i < 13; i++) {
	    if (days <= days2[i]) {
		*mo = i;
		*da = (days - days2[i - 1]);
		goto out1;
	    }
	}
    } else {
	for (i = 0; i < 13; i++) {
	    if (days <= days1[i]) {
		*mo = i;
		*da = (days - days1[i - 1]);
		goto out1;
	    }
	}
    }
  out1:;
}

/*
   return the number of days from the beginning of the year 'yr'
*/
int getndays(double j)
{
    int m, d, y, hh, mm;
    double ss;

    calcdate(j, &m, &d, &y, &hh, &mm, &ss);
    if (leapyear(y)) {
	return days2[m - 1] + d;
    } else {
	return days1[m - 1] + d;
    }
}

/*
   return hms
*/
int gethms(double j, int *h, int *m, double *s)
{
    double rem = j - (long) j;

    *h = (int) (rem * 24);
    *m = (int) ((rem * 24 - *h) * 60);
    *s = (int) (((rem * 24 - *h) * 60 - *m) * 60);
    *h = (*h + 12) % 24;
}

/*
 * strip special chars from a string
 */
void stripspecial(char *s, char *cs)
{
    int i, slen = strlen(s), curcnt = 0;

    for (i = 0; i < slen; i++) {
	if (s[i] == '\\' && isdigit(s[i + 1])) {
	    i++;
	} else if (s[i] == '\\' && isoneof(s[i + 1], "cCbxsSNuU+-")) {
	    i++;
	} else if (s[i] == '\\' && s[i + 1] == '\\') {
	    i++;
	} else {
	    cs[curcnt++] = s[i];
	}
    }
    cs[curcnt] = 0;
}
