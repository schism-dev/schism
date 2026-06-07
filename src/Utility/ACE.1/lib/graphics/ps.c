/* $Id: ps.c,v 1.3 2003/10/17 02:38:47 pturner Exp $
 *
 * Copyright 1990-2003 Oregon Health and Science University
 *
 * driver for postscript printer
 *
 * courtesy of:
 *
 * Jim Hudgens
 * hudgens@ray.met.fsu.edu
 *
 * Further modifications by,
 * Ole Holm Nielsen
 * ohnielse@ltf.dth.dk
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "externs.h"

extern char version[];
extern double charsize;
extern double devcharsize;
extern int ptofile;
extern char printstr[];

void stripspecial(char *s, char *cs);
int isoneof(int c, char *s);

int epsflag = 1;

#define PORTRAIT 0
#define LANDSCAPE 1

static void putps();

/*
 * printer control string
 */
#ifndef PS_PRSTR
char ps_prstr[128] = "lp -c ";

#else
char ps_prstr[128] = PS_PRSTR;

#endif

/* postscript page at scale = 0.25 */

/*
 * the following defines are tuned for our HP-LJ-III
 * and may need adjustment for other printers
 */
#define PSXMIN 75
#define PSXMAX 2375
#define PSYMIN 75
#define PSYMAX 3095
#define DXPS 2300
#define DYPS 3020
#define CHARS 1.8

#define MINCOLOR 0
#define MAXCOLOR 30
#define PSMAXPAT 30
#define MINLINEWIDTH 0
#define MAXLINEWIDTH 9
#define MAXLINESTYLE 14

/*
#define PORTRAIT 0
#define LANDSCAPE 1
*/

static int psxmin = PSXMIN;
static int psxmax = PSXMAX;
static int psymin = PSYMIN;
static int psymax = PSYMAX;
static int psdx = DXPS;
static int psdy = DYPS;
static int pscolor = -1;
static int pslinewidth = -1;
static int psdmode;
static int pspattern = 0;
static int psfont = 0;
static double pscharsize = 1.7;
static int pslinestyle;
static char *fname;
static int orientflag = PORTRAIT;

/* reduced from 1000 to 800 PJT */

#define MAXPATHLEN 800		/* MAXPATHLEN points in a path between
				 * strokes */
static int pathlength = 0;
double xconv(), yconv();
static FILE *psout;

static int prevx = 99999, prevy = 99999, prevmode;

static void stroke()
{
    if (pathlength) {
	fprintf(psout, "stroke\n");
	prevx = 99999;
	prevy = 99999;
	pathlength = 0;
    }
}

int pssetmode(int mode)
{
    static char tbuf[128];
    char sysbuf[128];
    if (mode % 2) {
	if (!ptofile) {
	    strcpy(tbuf, "/tmp/ACEgrXXXXXX");
	    mkstemp(tbuf);
	    fname = tbuf;
	} else {
	    fname = printstr;
	}
	if ((psout = fopen(fname, "w")) == NULL) {
	    return 0;
	}
    }
    switch (mode) {

    case 3:			/* EPS   portrait */
	orientflag = PORTRAIT;
	pscharsize = CHARS;
	if (epsflag) {
	    psxmin = PSXMIN;
	    psymin = PSYMIN;
	    if (devwidth == 0) {
		devwidth = PSXMAX - PSXMIN;
	    }
	    if (devheight == 0) {
		devheight = PSYMAX - PSYMIN;
	    }
	    if (((double) devheight) / devwidth > ((double) (PSYMAX - PSYMIN)) / (PSXMAX - PSXMIN)) {
		psymax = PSYMAX;
		psxmax = (PSYMAX - PSYMIN) * devwidth / devheight + PSXMIN;
		psdy = PSYMAX - PSYMIN;
		psdx = (PSYMAX - PSYMIN) * devwidth / devheight;
	    } else {
		psxmax = PSXMAX;
		psymax = (PSXMAX - PSXMIN) * devheight / devwidth + PSYMIN;
		psdx = PSXMAX - PSXMIN;
		psdy = (PSXMAX - PSXMIN) * devheight / devwidth;
	    }

	    devwidth = psdx;
	    devheight = psdy;
	    devwidthmm = 25.4 * devwidth / 300.0;
	    devheightmm = 25.4 * devheight / 300.0;

	} else {
	    psxmin = PSXMIN;
	    psxmax = PSXMAX;
	    psymin = PSYMIN;
	    psymax = PSYMAX;
	    psdx = DXPS;
	    psdy = DYPS;
	    devwidth = DXPS;
	    devheight = DYPS;
	    devwidthmm = 204;
	    devheightmm = 266;
	}
	devoffsx = 0;
	devoffsy = 0;
	break;

    case 1:			/* EPS landscape */
	orientflag = LANDSCAPE;
	pscharsize = CHARS;
	if (epsflag) {
	    psxmin = PSYMIN;
	    psymin = PSXMIN;
	    if (devwidth == 0) {
		devwidth = PSYMAX - PSYMIN;
	    }
	    if (devheight == 0) {
		devheight = PSXMAX - PSXMIN;
	    }
	    if (((double) devheight) / devwidth > ((double) (PSXMAX - PSXMIN)) / (PSYMAX - PSYMIN)) {
		psymax = PSXMAX;
		psxmax = (PSXMAX - PSXMIN) * devwidth / devheight + PSYMIN;
		psdy = PSXMAX - PSXMIN;
		psdx = (PSXMAX - PSXMIN) * devwidth / devheight;
	    } else {
		psxmax = PSYMAX;
		psymax = (PSYMAX - PSYMIN) * devheight / devwidth + PSXMIN;
		psdx = PSYMAX - PSYMIN;
		psdy = (PSYMAX - PSYMIN) * devheight / devwidth;
	    }
	    devwidth = psdx;
	    devheight = psdy;
	    devwidthmm = 25.4 * devwidth / 300.0;
	    devheightmm = 25.4 * devheight / 300.0;;
	} else {
	    psxmin = PSYMIN;
	    psxmax = PSYMAX;
	    psymin = PSXMIN;
	    psymax = PSXMAX;
	    psdx = DYPS;
	    psdy = DXPS;
	    devwidth = DYPS;
	    devheight = DXPS;
	    devwidthmm = 266;
	    devheightmm = 204;
	}
	devoffsx = 0;
	devoffsy = 0;
	break;

    case 2:
    case 4:
	stroke();
	fprintf(psout, "showpage\n");
	fprintf(psout, "%%%%Trailer\n");
	fclose(psout);
	if (!ptofile) {
	    sprintf(sysbuf, "%s %s", ps_prstr, fname);
	    system(sysbuf);
	    sleep(2);		/* some systems are having problems with the following unlink() */
	    unlink(fname);
	}
	orientflag = PORTRAIT;
	break;
    }
    return mode;
}



void drawps(x2, y2, mode)
int x2, y2, mode;
{
    int xtmp, ytmp;

    if (x2 < 0 || y2 < 0) {	/* Eliminate garbage on output */
	return;
    }
    xtmp = x2;
    ytmp = y2;

    if (mode) {
	if (prevmode && xtmp == prevx && ytmp == prevy) {
	    return;		/* previous mode was draw and points are the
				 * same */
	}
	fprintf(psout, "%d %d l\n", xtmp, ytmp);	/* lineto */
    } else {
	/* Avoid excessive moveto's generated by grtool */
	if (xtmp == prevx && ytmp == prevy) {
	    return;
	}
	fprintf(psout, "%d %d m\n", xtmp, ytmp);	/* moveto */
    }
    pathlength++;
    prevx = xtmp;
    prevy = ytmp;

    /*
     * Printers have some maximum number of points in a path. See PostScript
     * Language Reference Manual (Red book), p. 261. Hence the fix that
     * follows 
     */

    prevmode = mode;
    if (pathlength > MAXPATHLEN) {
	stroke();
	prevmode = 0;
	fprintf(psout, "%d %d m\n", xtmp, ytmp);	/* moveto */
    }
}

int xconvps(x)
double x;
{
    return ((int) (psxmin + psdx * xconv(x)));
}

int yconvps(y)
double y;
{
    return ((int) (psymin + psdy * yconv(y)));
}

int pssetcolor(c)
int c;
{
    extern unsigned char red[], green[], blue[];

    stroke();
    if (c != pscolor) {
	if (c >= 0) {
	    fprintf(psout, "%lf %lf %lf setrgbcolor\n", red[c] / 255.0, green[c] / 255.0, blue[c] / 255.0);
	    pscolor = 1;
	}
    }
    pscolor = c;
    return c;
}

int pssetlinewidth(c)
int c;
{
    stroke();
    if (c != pslinewidth) {
	c = c % (MAXLINEWIDTH + 1);
	if (c == 1)
	    fprintf(psout, "1 setlinewidth\n");
	else
	    fprintf(psout, "%d setlinewidth\n", (int) (3.5 * c + 0.51));
    }
    pslinewidth = c;
    return c;
}

void psdrawtic(x, y, dir, updown)
int x, y, dir, updown;
{
    switch (dir) {
    case 0:
	switch (updown) {
	case 0:
	    drawps(x, y, 0);
	    drawps(x, y + devxticl, 1);
	    break;
	case 1:
	    drawps(x, y, 0);
	    drawps(x, y - devxticl, 1);
	    break;
	}
	break;
    case 1:
	switch (updown) {
	case 0:
	    drawps(x, y, 0);
	    drawps(x + devyticl, y, 1);
	    break;
	case 1:
	    drawps(x, y, 0);
	    drawps(x - devyticl, y, 1);
	    break;
	}
	break;
    }
}

int pssetlinestyle(style)
int style;
{
    stroke();
    if (style == pslinestyle) {
	return (pslinestyle);
    }
    switch (style) {
    case 1:			/* solid */
	fprintf(psout, "[] 0 setdash\n");
	break;
    case 2:			/* dotted */
	fprintf(psout, "[4 8] 0 setdash\n");
	break;
    case 3:			/* long dash */
	fprintf(psout, "[20 20] 0 setdash\n");
	break;
    case 4:			/* short dash */
	fprintf(psout, "[40 20] 0 setdash\n");
	break;
    case 5:			/* dot-dashed */
	fprintf(psout, "[40 20 12 20] 0 setdash\n");
	break;
    }
    return (pslinestyle = style);
}

char pscurfont[1024] = "/Times-Roman findfont \n60 scalefont\n setfont";
int psfontsize = 60;

void pssetfont(n)
int n;
{
    if (psfont == n) {
	return;
    }
    switch (n) {
    case 0:
	sprintf(pscurfont, "/Times-Roman findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 1:
	sprintf(pscurfont, "/Times-Bold findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 2:
	sprintf(pscurfont, "/Times-Italic findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 3:
	sprintf(pscurfont, "/Times-BoldItalic findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 4:
	sprintf(pscurfont, "/Helvetica findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 5:
	sprintf(pscurfont, "/Helvetica-Bold findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 6:
	sprintf(pscurfont, "/Helvetica-Oblique findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 7:
	sprintf(pscurfont, "/Helvetica-BoldOblique findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 8:
	sprintf(pscurfont, "/Symbol findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 9:
	sprintf(pscurfont, "/Symbol findfont \n%d scalefont\n setfont", psfontsize);
	break;
    case 10:
	sprintf(pscurfont, "/Symbol findfont \n%d scalefont\n setfont", psfontsize);
	break;
    }
    fprintf(psout, "%s\n", pscurfont);
    hselectfont(psfont = n);
}

void pssetfontsize(size)
double size;
{
    int sf = psfont;

    psfontsize = (int) (size * 60);
    psfont = -1;
    pssetfont(sf);
}

static void escape_paren(s)
char *s;
{
    char t[256];
    int i, cnt = 0;
    for (i = 0; i < strlen(s); i++) {
	if (s[i] == '(' || s[i] == ')') {
	    t[cnt++] = '\\';
	}
	t[cnt++] = s[i];
    }
    t[cnt] = 0;
    strcpy(s, t);
}

void dispstrps(x, y, rot, s, just, fudge)
int x, y, rot, just, fudge;
char *s;
{
    char tmpstr[256];
    int i, cnt = 0;

    stroke();
    if (psfontsize == 0 || s == NULL || strlen(s) == 0) {
	return;
    }
    fprintf(psout, "%d %d m\n", x, y);
    fprintf(psout, "gsave\n");
    fprintf(psout, "%d %d translate\n", x, y);
    fprintf(psout, "%d rotate\n", rot);
    if (fudge) {
	fprintf(psout, "0 %d  m\n", -psfontsize / 3);
    } else {
	fprintf(psout, "0 0  m\n");
    }
    switch (just) {
    case 0:
	break;
    case 1:
	stripspecial(s, tmpstr);
	escape_paren(tmpstr);
	fprintf(psout, "(%s) RJ\n", tmpstr);
	break;
    case 2:
	stripspecial(s, tmpstr);
	escape_paren(tmpstr);
	fprintf(psout, "(%s) CS\n", tmpstr);
	break;
    }
    putps(s);
    fprintf(psout, "grestore\n");
    fprintf(psout, "newpath\n");
}

static void putps(s)
char *s;
{
    int i, cf, slen = strlen(s), curcnt = 0;
    int underline = 0, offset = 0;
    double saves = psfontsize / 60.0, scale = psfontsize / 60.0;
    char curstr[256];
    int upperset = 0;
    int symfont = 0;

    if (psfont == 9) {
	symfont = 1;
	upperset = 0x80;
    } else {
	symfont = 0;
	upperset = 0;
    }
    for (i = 0; i < slen; i++) {
	if (s[i] == '-' && isdigit(s[i + 1])) {
	    /*s[i] = 0261; */
	} else if (s[i] == '\\' && isdigit(s[i + 1])) {
	    curstr[curcnt] = 0;
	    if (curcnt >= 1) {
		fprintf(psout, "(%s) show\n", curstr);
	    }
	    curcnt = 0;
	    if (symfont) {
		symfont = 0;
		upperset = 0;
	    }
	    pssetfont(s[i + 1] - '0');
	    if (psfont == 9) {
		symfont = 1;
		upperset = 0x80;
	    }
	    i++;
	    continue;
	} else if (s[i] == '(' || s[i] == ')') {
	    curstr[curcnt++] = '\\';
	} else if (s[i] == '\\' && isoneof(s[i + 1], "cCbxsSNuU+-")) {
	    switch (s[i + 1]) {
	    case 'x':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		if (symfont == 0) {
		    symfont = 1;
		    upperset = 0x80;
		}
		pssetfont(10);
		i++;
		break;
	    case 's':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		pssetfontsize(scale = 0.6 * saves);
		offset -= psfontsize / 2;
		fprintf(psout, "0 %d rmoveto\n", -(psfontsize / 2));
		i++;
		break;
	    case 'S':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		pssetfontsize(scale = 0.6 * saves);
		offset += psfontsize;
		fprintf(psout, "0 %d rmoveto\n", psfontsize);
		i++;
		break;
	    case 'N':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		scale = saves;
		pssetfontsize(scale);
		fprintf(psout, "0 %d rmoveto\n", -offset);
		offset = 0;
/*
		fprintf(psout, "0 %d rmoveto\n", psfontsize);
*/
		i++;
		break;
	    case 'b':
		i++;
		break;
	    case 'c':
		upperset = 0x80;
		i++;
		break;
	    case 'C':
		upperset = 0;
		i++;
		break;
	    case 'u':
		underline = 1;
		i++;
		break;
	    case 'U':
		underline = 0;
		i++;
		break;
	    case '-':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		scale -= 0.2;
		if (scale < 0.2) {
		    scale = 0.2;
		}
		pssetfontsize(scale);
		i++;
		break;
	    case '+':
		curstr[curcnt] = 0;
		if (curcnt >= 1) {
		    fprintf(psout, "(%s) show\n", curstr);
		}
		curcnt = 0;
		scale += 0.2;
		pssetfontsize(scale);
		i++;
		break;
	    }
	    continue;
	} else if (s[i] == '\\' && s[i + 1] == '\\') {
	    curstr[curcnt++] = '\\';
	    curstr[curcnt++] = s[i];
	    i++;
	    continue;
	}
	curstr[curcnt++] = s[i] + upperset;
    }
    curstr[curcnt] = 0;
    fprintf(psout, "(%s) show\n", curstr);
}

int pssetpat(k)
int k;
{
    stroke();
    if (k > PSMAXPAT) {
	k = PSMAXPAT;
    } else if (k < 0) {
	k = 0;
	fprintf(psout, "0.0 setgray\n");
    }
    return (pspattern = k);
}

void psfill(n, px, py)
int n;
int px[], py[];

{
    int i;

    stroke();
    drawps(px[0], py[0], 0);
    for (i = 1; i < n; i++) {
	drawps(px[i], py[i], 1);
    }
    fprintf(psout, "closepath\n");
    fprintf(psout, "%lf setgray\n", 1.0 - pspattern / (double) PSMAXPAT);
    fprintf(psout, "gsave eofill grestore\n");
    stroke();
    fprintf(psout, "0 setgray\n");
}

void psfillcolor(n, px, py)
int n;
int px[], py[];

{
    int i;

    stroke();
    drawps(px[0], py[0], 0);
    for (i = 1; i < n; i++) {
	drawps(px[i], py[i], 1);
    }
    fprintf(psout, "closepath\n");
    fprintf(psout, "gsave eofill grestore\n");
    stroke();
}

void psdrawarc(x, y, r)
int x, y, r;
{
    stroke();
    fprintf(psout, "%d %d %d %d %d arc\n", x, y, r, 0, 360);
    fprintf(psout, "stroke\n");
}

void psfillarc(x, y, r)
int x, y, r;
{
    stroke();
    fprintf(psout, "%d %d %d %d %d arc\n", x, y, r, 0, 360);
    fprintf(psout, "gsave fill grestore\n");
    fprintf(psout, "stroke\n");
}

void psdrawellipse(x, y, xm, ym)
int x, y, xm, ym;
{
    double scalex = (double) xm / (double) ym, scaley = 1.0;

    stroke();
    fprintf(psout, "gsave\n");
    fprintf(psout, "%lf %lf scale\n", scalex, scaley);
    fprintf(psout, "%d %d %d %d %d arc\n", (int) (x * 1.0 / scalex), y, ym, 0, 360);
    fprintf(psout, "stroke\n");
    fprintf(psout, "grestore\n");
}

void psfillellipse(x, y, xm, ym)
int x, y, xm, ym;
{
    double scalex = (double) xm / (double) ym, scaley = 1.0;

    stroke();
    fprintf(psout, "gsave\n");
    fprintf(psout, "%lf %lf scale\n", (double) xm, (double) xm / ym);
    fprintf(psout, "%d %d %d %d %d arc\n", x, y, ym, 0, 360);
    fprintf(psout, "gsave fill grestore\n");
    fprintf(psout, "stroke\n");
    fprintf(psout, "grestore\n");
}

void psleavegraphics()
{
    pssetmode(psdmode + 1);
}

/*           postscript initialization routine  */
int psinitgraphics(dmode)
int dmode;
{

    psdmode = dmode;
    if (!pssetmode(psdmode)) {
	return -1;
    }
    devconvx = xconvps;
    devconvy = yconvps;
    vector = drawps;
    devwritestr = dispstrps;
    devsetcolor = pssetcolor;
    devsetfont = pssetfont;
    devsetline = pssetlinestyle;
    devsetlinew = pssetlinewidth;
    devdrawtic = psdrawtic;
    devsetpat = pssetpat;
    devdrawarc = psdrawarc;
    devfillarc = psfillarc;
    devfill = psfill;
    devfillcolor = psfillcolor;
    devdrawellipse = psdrawellipse;
    devfillellipse = psfillellipse;
    devleavegraphics = psleavegraphics;
    devcharsize = pscharsize;
    devsymsize = 20;
    devxticl = 20;
    devyticl = 20;
    devarrowlength = 20;

    fprintf(psout, "%%!PostScript\n");
    if (epsflag) {
	if (dmode == 1) {
	    fprintf(psout, "%%%%BoundingBox: %d %d %d %d\n", psxmin * 72 / 300, psymin * 72 / 300, psdy * 72 / 300, psdx * 72 / 300);
	} else {
	    fprintf(psout, "%%%%BoundingBox: %d %d %d %d\n", psxmin * 72 / 300, psymin * 72 / 300, psdx * 72 / 300, psdy * 72 / 300);
	}
    }
    fprintf(psout, "%%%%Creator: ACE/graphics\n");
    fprintf(psout, "%%%%EndComments\n");

    fprintf(psout, "/m {moveto} bind def\n");
    fprintf(psout, "/l {lineto} bind def\n");
    fprintf(psout, "/RJ {\n");
    fprintf(psout, " stringwidth neg exch neg exch\n");
    fprintf(psout, " rmoveto\n");
    fprintf(psout, "} bind def\n");

    fprintf(psout, "/CS {\n");
    fprintf(psout, " stringwidth\n");
    fprintf(psout, " 2 div neg exch 2 div neg exch\n");
    fprintf(psout, " rmoveto\n");
    fprintf(psout, "} bind def\n");

    fprintf(psout, "0.24 0.24 scale\n");
    fprintf(psout, "1 setlinecap\n");
/*
 * rotate if in landscape mode
 */
    if (dmode == 1) {
	fprintf(psout, "%d 0 translate\n", 2 * psymin + psdy);
	fprintf(psout, "90 rotate\n");
    }
    pssetcolor(1);
    pssetlinewidth(1);
    setlinestyle(0);
    psfont = -1;
    setfont(2);
    return 0;
}
