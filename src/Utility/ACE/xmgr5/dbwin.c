/* $Id: dbwin.c,v 1.10 2008/10/08 22:53:19 pturner Exp $
 *
 * DB file access
 *
 */

#include <stdio.h>
#ifndef WIN32
#include <sys/param.h>
#endif

#ifdef PGSQL

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

#include <libpq-fe.h>

#include "globals.h"
#include "motifinc.h"
#include "noxprotos.h"

#define OVERWRITE_SET 0
#define APPEND_SET 1
#define MAX_ARRAY 500

/*
 * a file popup and a text item
 */
typedef struct {
    Widget file;
    Widget browser;
} Browser;

/***************************************************************
 *
 * Read DB data section.
 *
 ***************************************************************/

typedef struct {
    Widget top;
    Widget *setChoice;
    Widget *stationChoice;
    Widget *instrumentIDChoice;
    Widget *varChoice;
    Widget *yearStartChoice;
    Widget *monthStartChoice;
    Widget *dayStartChoice;
    Widget *hourStartChoice;
    Widget *minuteStartChoice;
    Widget *secondStartChoice;
    Widget *yearEndChoice;
    Widget *monthEndChoice;
    Widget *dayEndChoice;
    Widget *hourEndChoice;
    Widget *minuteEndChoice;
    Widget *secondEndChoice;
    Widget autoscale;
    Widget stride;
    Browser b;
} DBUI;

DBUI dbui;

int GetSites(void)
{
    int i, j, nFields;
    char buf[2048];
    char *dbName;
    char *tableName;
    PGconn *conn;
    PGresult *res;
    dbName = "telemetry";
    conn = PQconnectdb("dbname=telemetry host=amb104.ccalmr.ogi.edu user=nobody");
    if (PQstatus(conn) == CONNECTION_BAD) {
	fprintf(stderr, "Connection to database '%s' failed.\n", dbName);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQfinish(conn);
	return 1;
    }
    res = PQexec(conn, "BEGIN");
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "BEGIN command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	PQfinish(conn);
	return 1;
    }
    PQclear(res);
    res = PQexec(conn, "DECLARE mycursor CURSOR FOR SELECT site FROM sites where site like '%adp' order by site");
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "DECLARE CURSOR command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	return 1;
    }
    PQclear(res);
    res = PQexec(conn, "FETCH ALL in mycursor");
    if (!res || PQresultStatus(res) != PGRES_TUPLES_OK) {
	fprintf(stderr, "FETCH ALL command didn't return tuples properly\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	return 1;
    }

    /* first, print out the attribute names */
    nFields = PQnfields(res);
    for (i = 0; i < nFields; i++) {
	printf("%-15s", PQfname(res, i));
    }

    for (i = 0; i < PQntuples(res); i++) {
	for (j = 0; j < nFields; j++) {
	    printf("%-15s", PQgetvalue(res, i, j));
	}
	printf("\n");
    }

    PQclear(res);
    /* close the cursor */
    res = PQexec(conn, "CLOSE mycursor");
    PQclear(res);

    /* commit the transaction */
    res = PQexec(conn, "COMMIT");

    PQclear(res);
    /* close the connection to the database and cleanup */
    PQfinish(conn);
}

void createADPdb(Widget w, XtPointer client_data, XtPointer call_data)
{
    GetSites();
}

char dbhost[1024];
char dbase[1024];
char dbuser[1024];
char dbpasswd[1024];

int SetDBHost(char *host, char *db);
int SetDBUser(char *user, char *passwd);
int ReadDB(int gno, int setno, char *site, char *iid, char *xvar, char *yvar, double cday1, double cday2);
int ReadDBADP(int gno, int setno, char *site, char *iid, char *xvar, char *yvar, int bin, double cday1, double cday2);
int ReadDBSQL(int gno, int setno, char *sql);

int SetDBDefaults(void)
{
    strcpy(dbhost, "cdb02.stccmop.org");
    strcpy(dbase, "cmop");
    strcpy(dbuser, "reader");
    strcpy(dbpasswd, "");
}

int SetDBHost(char *host, char *db)
{
    int len1 = strlen(host);
    int len2 = strlen(db);
    strcpy(dbhost, host);
    strcpy(dbase, db);
    return 0;
}

int SetDBUser(char *user, char *passwd)
{
    strcpy(dbuser, user);
    strcpy(dbpasswd, passwd);
    return 0;
}

int ReadDB(int gno, int setno, char *site, char *iid, char *xvar, char *yvar, double cday1, double cday2)
{
    int i, j, nFields, nItems;
    char buf[2048];
    char *dbName;
    char *tableName;
    PGconn *conn;
    PGresult *res;
    double *x, *y;
    sprintf(buf, "dbname=%s host=%s user=nobody", dbase, dbhost);
    conn = PQconnectdb(buf);
    if (PQstatus(conn) == CONNECTION_BAD) {
	fprintf(stderr, "Connection to database '%s' failed.\n", dbName);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQfinish(conn);
	return 1;
    }
    res = PQexec(conn, "BEGIN");
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "BEGIN command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	PQfinish(conn);
	return 1;
    }
    PQclear(res);
    if (strcmp(iid, "NULL") == 0) {
	sprintf(buf, "DECLARE mycursor CURSOR FOR SELECT %s, %s FROM %s WHERE (corieday >= %lf AND corieday <= %lf) ORDER by %s", xvar, yvar, site, cday1, cday2, xvar);
    } else {
	sprintf(buf, "DECLARE mycursor CURSOR FOR SELECT %s, %s FROM %s WHERE (instrumentid = \'%s\' AND %s > -9999.0 AND corieday >= %lf AND corieday <= %lf) ORDER by %s", xvar, yvar, site, iid, yvar, cday1, cday2, xvar);
    }
    if (debuglevel == 11) {
        printf("%s\n", buf);
    }
    res = PQexec(conn, buf);
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "DECLARE CURSOR command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
        PQfinish(conn);
	return 1;
    }
    PQclear(res);
    res = PQexec(conn, "FETCH ALL in mycursor");
    if (!res || PQresultStatus(res) != PGRES_TUPLES_OK) {
	fprintf(stderr, "FETCH ALL command didn't return tuples properly\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
        PQfinish(conn);
	return 1;
    }

    /* first, print out the attribute names */
    nFields = PQnfields(res);
    nItems = PQntuples(res);
    if (nFields != 2 || nItems <= 0) {
	PQclear(res);
        PQfinish(conn);
	return 1;
    }
    if (setno == -1) {
	setno = nextset(gno);
    } else {
        if (isactive(gno, setno)) {
	    killset(gno, setno);
        }
    }
/*
 * allocate for this set
 */
    x = (double *) malloc(nItems * sizeof(double));
    y = (double *) malloc(nItems * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	return 1;
    }
    for (i = 0; i < nItems; i++) {
	x[i] = atof(PQgetvalue(res, i, 0));
	y[i] = atof(PQgetvalue(res, i, 1));
	//printf("%.3lf %.3lf\n", atof(PQgetvalue(res, i, 0)), atof(PQgetvalue(res, i, 1)));
    }

    PQclear(res);
    /* close the cursor */
    res = PQexec(conn, "CLOSE mycursor");
    PQclear(res);

    /* commit the transaction */
    res = PQexec(conn, "COMMIT");

    PQclear(res);
    /* close the connection to the database and cleanup */
    PQfinish(conn);
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, nItems, 0);
    setcol(gno, y, setno, nItems, 1);

    sprintf(buf, "%s %s", site, iid);
    setcomment(gno, setno, buf);
    updatesetminmax(gno, setno);

    return 0;
}

int ReadDBADP(int gno, int setno, char *site, char *iid, char *xvar, char *yvar, int bin, double cday1, double cday2)
{
    int i, j, nFields, nItems, nids;
    char buf[2048];
    char *dbName;
    char *tableName;
    PGconn *conn;
    PGresult *res;
    double *x, *y;
    //sprintf(buf, "dbname=%s host=%s user=corie password=db4ccalmr", dbase, dbhost);
    sprintf(buf, "dbname=%s host=%s user=%s password=%s", dbase, dbhost, dbuser, dbpasswd);
    conn = PQconnectdb(buf);
    if (PQstatus(conn) == CONNECTION_BAD) {
	fprintf(stderr, "Connection to database '%s' failed.\n", dbName);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQfinish(conn);
	return 1;
    }
    res = PQexec(conn, "BEGIN");
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "BEGIN command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	PQfinish(conn);
	return 1;
    }
    PQclear(res);
    sprintf(buf, "DECLARE mycursor CURSOR FOR select %s, %s from %s, %sdata where (instrumentid = \'%s\' AND corieday >= %lf AND corieday <= %lf and bin = %d and %s.id = %sdata.recordid) order by corieday",
	xvar, yvar, site, site, iid, cday1, cday2, bin, site, site);
    if (debuglevel == 11) {
        printf("%s\n", buf);
    }
    res = PQexec(conn, buf);
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "DECLARE CURSOR command failed: %s\n", buf);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	return 1;
    }
    PQclear(res);
    res = PQexec(conn, "FETCH ALL in mycursor");
    if (!res || PQresultStatus(res) != PGRES_TUPLES_OK) {
	fprintf(stderr, "FETCH ALL command didn't return tuples properly\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	return 1;
    }
    nids = PQntuples(res);
/*
 * allocate for this set
 */
    x = (double *) malloc(nids * sizeof(double));
    y = (double *) malloc(nids * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	return 1;
    }
    for (i = 0; i < nids; i++) {
	x[i] = atof(PQgetvalue(res, i, 0));
	y[i] = atof(PQgetvalue(res, i, 1));
        /*printf("%lf %lf\n", x[i], y[i]);*/
    }
    PQclear(res);

    /* close the cursor */
    res = PQexec(conn, "CLOSE mycursor");
    PQclear(res);

        /* commit the transaction */
    res = PQexec(conn, "COMMIT");
    PQclear(res);
    /* close the connection to the database and cleanup */
    PQfinish(conn);

    if (isactive(gno, setno)) {
	killset(gno, setno);
    }
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, nids, 0);
    setcol(gno, y, setno, nids, 1);

    sprintf(buf, "Site %s, IID %s (x, y) = (%s, %s)", site, iid, xvar, yvar);
    setcomment(gno, setno, buf);
    updatesetminmax(gno, setno);

    return 0;
}

int ReadDBSQL(int gno, int setno, char *sql)
{
    int i, j, nFields, nItems;
    char buf[2048];
    char *dbName;
    char *tableName;
    PGconn *conn;
    PGresult *res;
    double *x, *y;
    //sprintf(buf, "dbname=%s host=%s user=nobody", dbase, dbhost);
    sprintf(buf, "dbname=%s host=%s user=%s password=%s", dbase, dbhost, dbuser, dbpasswd);
    conn = PQconnectdb(buf);
    if (PQstatus(conn) == CONNECTION_BAD) {
	fprintf(stderr, "Connection to database '%s' failed.\n", dbase);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQfinish(conn);
	return 1;
    }
    res = PQexec(conn, "BEGIN");
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "BEGIN command failed\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
	PQfinish(conn);
	return 1;
    }
    PQclear(res);
    sprintf(buf, "DECLARE mycursor CURSOR for %s\n", sql);
    if (debuglevel == 11) {
        printf("%s\n", buf);
    }
    res = PQexec(conn, buf);
    if (!res || PQresultStatus(res) != PGRES_COMMAND_OK) {
	fprintf(stderr, "DECLARE CURSOR command failed: %s\n", buf);
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
        PQfinish(conn);
	return 1;
    }
    PQclear(res);
    res = PQexec(conn, "FETCH ALL in mycursor");
    if (!res || PQresultStatus(res) != PGRES_TUPLES_OK) {
	fprintf(stderr, "FETCH ALL command didn't return tuples properly\n");
	fprintf(stderr, "%s", PQerrorMessage(conn));
	PQclear(res);
        PQfinish(conn);
	return 1;
    }

    /* first, print out the attribute names */
    nFields = PQnfields(res);
    nItems = PQntuples(res);
    if (nFields != 2 || nItems <= 0) {
	PQclear(res);
        PQfinish(conn);
	return 1;
    }
    if (gno == -1) {
	gno = cg;
    }
    if (setno == -1) {
	setno = nextset(gno);
    } else {
        if (isactive(gno, setno)) {
	    killset(gno, setno);
        }
    }
/*
 * allocate for this set
 */
    x = (double *) malloc(nItems * sizeof(double));
    y = (double *) malloc(nItems * sizeof(double));
    if (x == NULL || y == NULL) {
	errwin("Insufficient memory for set");
	cxfree(x);
	cxfree(y);
	return 1;
    }
    for (i = 0; i < nItems; i++) {
	x[i] = atof(PQgetvalue(res, i, 0));
	y[i] = atof(PQgetvalue(res, i, 1));
	//printf("%.3lf %.3lf\n", atof(PQgetvalue(res, i, 0)), atof(PQgetvalue(res, i, 1)));
    }

    PQclear(res);
    /* close the cursor */
    res = PQexec(conn, "CLOSE mycursor");
    PQclear(res);

    /* commit the transaction */
    res = PQexec(conn, "COMMIT");

    PQclear(res);
    /* close the connection to the database and cleanup */
    PQfinish(conn);
    activateset(gno, setno);
    settype(gno, setno, XY);
    setcol(gno, x, setno, nItems, 0);
    setcol(gno, y, setno, nItems, 1);

    setcomment(gno, setno, sql);
    updatesetminmax(gno, setno);

    return 0;
}

#endif /* PGSQL */
