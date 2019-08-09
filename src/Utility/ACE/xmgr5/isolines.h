/*
 * $Id: isolines.h,v 1.1.1.1 2003/07/21 16:18:41 pturner Exp $
 * $Source: /home/workspace/ccalmr/src/ace/xmgr5/isolines.h,v $
 *
 * Isolines definition module
 *
 */

#define MAXISOLINES 50

typedef struct _Props {
    int color;
    int linew;
    int lines;
    int font;
    int format;
    int prec;
    int points;
    double charsize;
    int symbol;
    int symsize;
    int fill;
    int fillusing;
    int fillcol;
    int fillpat;
    int arrow;
    int atype;
    double asize;
} Props;
 
typedef struct isolparms {
    int active;
    int type;
    Props p;
    int nisol;
    int nsplits;
    int marker;
    int markstep;
    int numbs;
    int loctype;
    double x;
    double y;
    double xlen;                   /* in character widths */
    double ylen;                   /* in character widths */
    double xgap;                   /* in character widths */
    double ygap;                   /* in character widths */
    int lactive;                /* legend active */
    int layout;                 /* legend HORIZONTAL or VERTICAL */
    int llabels;                /* legend labels on or off */
    int isoltype;
    double cmin;
    double cmax;
    double cint;
    double cis[MAXISOLINES];
    int color[MAXISOLINES];
    int linew[MAXISOLINES];
    int lines[MAXISOLINES];
    int writeflag;
    char wname[256];
    int writelevel[MAXISOLINES];
} Isolparms;
