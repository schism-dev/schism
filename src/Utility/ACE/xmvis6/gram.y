/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

%{
/*
 * 
 * evaluate expressions, commands, parameter files
 * 
 */

#define GRAMMAR

#ifndef lint
static char RCSid[] = "$Id: gram.y,v 1.31 2008/04/10 18:29:08 pturner Exp $";
#endif

#include <stdio.h>
#include <math.h>
#include <ctype.h>

#include "defines.h"
#include "globals.h"

typedef struct _symtab_entry {
    char *s;
    int type;
} symtab_entry;

#ifndef M_PI
#     define M_PI  3.14159265358979323846
#endif

#ifndef TRUE
#     define TRUE 1
#endif

#ifndef FALSE
#     define FALSE 0
#endif
/* for LINUX
char *gettxt(char *t, char *s) { return *s; }
*/

int elcirc_maxlevels;
int elcirc_gridno;
int elcirc_flowno;
int elcircmarker = 0;
int curg = 0;

int maxboxes = MAXBOXES;
int maxlines = MAXLINES;
int maxstring = MAXSTR;

double result, resx, resy;	/* return value if expression */
double nonl_parms[10];

double drand48(void);
long lrand48(void);

double rnorm(double mean, double sdev), fx(double x), normp(double b, double *s);
void yyerror(char *s);

static int interr;

static double *freelist[100]; 	/* temporary vectors */
static int fcnt;		/* number allocated */

int naxis = 0;	/* current axis */
int curline, curbox, curstring, curleg, curobject;

int gotbatch, gotparams, gotread; /* these guys attempt to avoid reentrancy problems */
int readtype, readsrc;
extern char batchfile[];
char paramfile[256], readfile[256];

static char f_string[512];	/* buffer for string to parse */
static int pos = 0;
static double *aa, *bb, *cc, *dd, *xx, *yy;
static int setindex, lxy, ls;
static int setsetno;
static int whichgraph;
static int whichset;

extern int change_gno;
extern int change_type;

int checkptr(void *ptr, char *buf);
static Isolparms *setisol = NULL; /* pointer to current Isolparms struct */
static DisplayFlow *setflow = NULL; /* pointer to current Isolparms struct */
static Props *setprops = NULL; /* pointer to current Isolparms struct */
static Zoom_box *setzoombox = NULL;
static Hist_marker *sethistbox = NULL;
static DisplaySlice *setslice = NULL;
static Elevmarker *setelevmarker = NULL;
static DisplayGrid *setgrid = NULL;
static Transect *settrans = NULL;
static ADCIRC3D *setadc3d = NULL;
static Display3dFlow *setflow3d = NULL;
static DisplayParticles *setdrogs = NULL;

/* may add these later TODO
*/

%}

%union {
    double val;
    int ival;
    double *ptr;
    int func;
    int pset;
    char *str;
}

%token <ptr> VAR 
%token <ptr> X
%token <ptr> Y

%token <str> CHRSTR

%token <val> FITPARM
%token <val> NUMBER

%token <func> ABS 
%token <func> ACOS
%token <func> ASIN
%token <func> ATAN
%token <func> ATAN2
%token <func> CEIL 
%token <func> COS
%token <func> DEG
%token <func> DX
%token <func> DY
%token <func> ERF
%token <func> ERFC
%token <func> EXP
%token <func> FLOOR 
%token <func> HYPOT
%token <func> INDEX
%token <func> INT
%token <func> IRAND
%token <func> LGAMMA 
%token <func> LN
%token <func> LOG 
%token <func> LOGISTIC 
%token <func> MAXP
%token <func> MINP
%token <func> MINMAX
%token <func> MOD 
%token <func> NORM
%token <func> NORMP
%token <func> PI 
%token <func> RAD
%token <func> RAND
%token <func> SETNO
%token <func> SIN
%token <func> SQR 
%token <func> SQRT
%token <func> TAN 

%token <ival> INUM

%token <pset> ABORT
%token <pset> ABOVE
%token <pset> ABSOLUTE
%token <pset> ACTIVATE
%token <pset> ACTIVE
%token <pset> ADCIRC
%token <pset> ADCIRC3DFLOW
%token <pset> ALL
%token <pset> ALT
%token <pset> ALTERNATE
%token <pset> ALTXAXIS
%token <pset> ALTYAXIS
%token <pset> AMP
%token <pset> ANGLE
%token <pset> ANNOTATE
%token <pset> APPEND
%token <pset> AREA
%token <pset> ARROW
%token <pset> ASCEND
%token <pset> AT
%token <pset> ATTACH
%token <pset> AUTO
%token <pset> AUTOSCALE
%token <pset> AUTOTICKS
%token <pset> AVERAGE
%token <pset> AVG
%token <pset> AXES
%token <pset> AXIS
%token <pset> BACKBUFFER
%token <pset> BACKGROUND
%token <pset> BAR
%token <pset> BATCH
%token <pset> BATH
%token <pset> BATHYMETRY
%token <pset> COURANT
%token <pset> BELOW
%token <pset> BIN
%token <pset> BINARY
%token <pset> BOTH
%token <pset> BOTTOM
%token <pset> BOUNDARY
%token <pset> BOX
%token <pset> CELLS
%token <pset> CENTER
%token <pset> CH3D
%token <pset> CHAR
%token <pset> CHDIR
%token <pset> CIRCLE
%token <pset> CLEAR
%token <pset> CLICK
%token <pset> CLOCK
%token <pset> CLOSE 
%token <pset> CM
%token <pset> CMAP
%token <pset> COLOR
%token <pset> COLORMAP
%token <pset> COMMENT
%token <pset> CONC
%token <pset> CONCENTRATION
%token <pset> CONCENTRATIONS
%token <pset> COPY
%token <pset> CROSS
%token <pset> CYCLE
%token <pset> DAYMONTH
%token <pset> DAYOFWEEKL
%token <pset> DAYOFWEEKS
%token <pset> DAYOFYEAR
%token <pset> DAYS
%token <pset> DDMMYY
%token <pset> DDMONTHSYYHHMMSS
%token <pset> DECIMAL
%token <pset> DEF
%token <pset> DEFAULT
%token <pset> DEGREESLAT
%token <pset> DEGREESLON
%token <pset> DEGREESMMLAT
%token <pset> DEGREESMMLON
%token <pset> DEGREESMMSSLAT
%token <pset> DEGREESMMSSLON
%token <pset> DELAYP
%token <pset> DELETE
%token <pset> DEPTH
%token <pset> DEPTHS
%token <pset> DESCEND
%token <pset> DEVICE
%token <pset> DEVXY
%token <pset> DFT
%token <pset> DT
%token <pset> DIAMOND
%token <pset> DIFFERENCE
%token <pset> DISK
%token <pset> DISPLAY
%token <pset> DOT
%token <pset> DOUBLEBUFFER
%token <pset> DOWN
%token <pset> DRAW2
%token <pset> DROGUE
%token <pset> DROGUES
%token <pset> DRY
%token <pset> DXDX
%token <pset> DXP
%token <pset> DYDY
%token <pset> DYP
%token <pset> ECHO
%token <pset> EDIT
%token <pset> ELA
%token <pset> ELCIRC
%token <pset> ELEMENT
%token <pset> ELEMENTS
%token <pset> ELEV
%token <pset> ELEVATION
%token <pset> ELEVATIONS
%token <pset> ELEVMARKER
%token <pset> ELLIPSE
%token <pset> ELLIPSES
%token <pset> ELLIPSEZ
%token <pset> ELSE
%token <pset> END
%token <pset> ERRORBAR
%token <pset> EXIT
%token <pset> EXPAND
%token <pset> EXPONENTIAL
%token <pset> FACTOR
%token <pset> FALSEP
%token <pset> FAST
%token <pset> FEET
%token <pset> FFT
%token <pset> FILEP
%token <pset> FILL
%token <pset> FIND
%token <pset> FIXEDPOINT
%token <pset> FLOW
%token <pset> FLUSH
%token <pset> FLUX
%token <pset> FOCUS
%token <pset> FOLLOWS
%token <pset> FONTP
%token <pset> FOREGROUND
%token <pset> FORMAT
%token <pset> FORT14
%token <pset> FORT63
%token <pset> FORT64
%token <pset> FORWARD
%token <pset> FRAMEP
%token <pset> FREQ
%token <pset> FRONTBUFFER
%token <pset> GENERAL
%token <pset> GETP
%token <pset> GOTO
%token <pset> GRAPH
%token <pset> GRAPHNO
%token <pset> GRAPHS
%token <pset> GRAPHTYPE
%token <pset> GRID
%token <pset> HARDCOPY
%token <pset> HBAR
%token <pset> HELP
%token <pset> HGAP
%token <pset> HIDDEN
%token <pset> HISTBOX
%token <pset> HISTO
%token <pset> HISTORY
%token <pset> HMS
%token <pset> HORIZONTAL
%token <pset> HOURS
%token <pset> HPGLL
%token <pset> HPGLP
%token <pset> IF
%token <pset> IGNORE
%token <pset> IHL
%token <pset> IMAGE
%token <pset> IMAGES
%token <pset> IN
%token <pset> INCLUDE
%token <pset> INFO
%token <pset> INIT
%token <pset> INITGRAPHICS
%token <pset> INOUT
%token <pset> INTEGRATE
%token <pset> INTERP
%token <pset> INUNDATION
%token <pset> INVDFT
%token <pset> INVFFT
%token <pset> ISOLINE
%token <pset> ISOLINES
%token <pset> JUST
%token <pset> KILL
%token <pset> KM
%token <pset> LABEL
%token <pset> LAYOUT
%token <pset> LEAVE
%token <pset> LEAVEGRAPHICS
%token <pset> LEFT
%token <pset> LEGEND
%token <pset> LENGTH
%token <pset> LEVEL
%token <pset> LEVELS
%token <pset> LIMITS
%token <pset> LINE
%token <pset> LINES
%token <pset> LINESTYLE
%token <pset> LINETO
%token <pset> LINEW
%token <pset> LINEWIDTH
%token <pset> LINK
%token <pset> LOAD
%token <pset> LOC
%token <pset> LOCATE
%token <pset> LOCATOR
%token <pset> LOCTYPE
%token <pset> LOGX
%token <pset> LOGXY
%token <pset> LOGY
%token <pset> M
%token <pset> MAG
%token <pset> MAGNITUDE
%token <pset> MAJOR
%token <pset> MAPSCALE
%token <pset> MARKER
%token <pset> MARKERS
%token <pset> MAXLEVELS
%token <pset> METHOD
%token <pset> MIFL
%token <pset> MIFP
%token <pset> MILES
%token <pset> MINOR
%token <pset> MINUTES
%token <pset> MISSINGP
%token <pset> MM
%token <pset> MMDD
%token <pset> MMDDHMS
%token <pset> MMDDYY
%token <pset> MMDDYYHMS
%token <pset> MMSSLAT
%token <pset> MMSSLON
%token <pset> MMYY
%token <pset> MONTHDAY
%token <pset> MONTHL
%token <pset> MONTHS
%token <pset> MOVE
%token <pset> MOVE2
%token <pset> MOVETO
%token <pset> NEGATE
%token <pset> NO
%token <pset> NODE
%token <pset> NODES
%token <pset> NONE
%token <pset> NORMAL
%token <pset> NORTH
%token <pset> NXY
%token <pset> OFF
%token <pset> OFFSETX
%token <pset> OFFSETY
%token <pset> ON
%token <pset> OP
%token <pset> OPEN
%token <pset> ORIENT
%token <pset> OUT
%token <pset> PAGE
%token <pset> PARA
%token <pset> PARALLEL
%token <pset> PARAMETERS
%token <pset> PARAMS
%token <pset> PARMS
%token <pset> PATTERN
%token <pset> PER
%token <pset> PERIMETER
%token <pset> PERP
%token <pset> PERPENDICULAR
%token <pset> PHASE
%token <pset> PIE
%token <pset> PIPE
%token <pset> PLACE
%token <pset> PLAN
%token <pset> PLUS
%token <pset> POINT
%token <pset> POLAR
%token <pset> POLY
%token <pset> POLYI
%token <pset> POLYO
%token <pset> POP
%token <pset> POWER
%token <pset> PREC
%token <pset> PREFIX
%token <pset> PREPEND
%token <pset> PRINT
%token <pset> PROFILE
%token <pset> PROP
%token <pset> PS
%token <pset> PSCOLORL
%token <pset> PSCOLORP
%token <pset> PSMONOL
%token <pset> PSMONOP
%token <pset> PUSH
%token <pset> PUTP
%token <pset> QUIT
%token <pset> READ
%token <pset> READBIN
%token <pset> REDRAW
%token <pset> REGION
%token <pset> REGIONS
%token <pset> REGNUM
%token <pset> REGRESS
%token <pset> REMOVE
%token <pset> RENDER
%token <pset> REPORT
%token <pset> RESET
%token <pset> REVERSE
%token <pset> REWIND
%token <pset> RIGHT
%token <pset> RISER
%token <pset> ROT
%token <pset> RUN
%token <pset> SALINITY
%token <pset> SAMPLE
%token <pset> SAVE
%token <pset> SCALAR
%token <pset> SCALE
%token <pset> SCIENTIFIC
%token <pset> SECONDS
%token <pset> SET
%token <pset> SETS
%token <pset> SHOW
%token <pset> SHRINK
%token <pset> SIGMA
%token <pset> SIGN
%token <pset> SIZE
%token <pset> SKIP
%token <pset> SLAB
%token <pset> SLEEP
%token <pset> SLICE
%token <pset> SOURCE
%token <pset> SPEC
%token <pset> SPECIFIED
%token <pset> SPECTRUM
%token <pset> SPLITS
%token <pset> SQUARE
%token <pset> STACK
%token <pset> STACKEDBAR
%token <pset> STACKEDHBAR
%token <pset> STACKEDLINE
%token <pset> STAGGER
%token <pset> STAR
%token <pset> START
%token <pset> STARTSTEP
%token <pset> STARTTYPE
%token <pset> STATION
%token <pset> STATUS
%token <pset> STEP
%token <pset> STOP
%token <pset> STREAMLINES
%token <pset> STRING
%token <pset> STRINGS
%token <pset> SUBTITLE
%token <pset> SURFACE
%token <pset> SWAPBUFFER
%token <pset> SYMBOL
%token <pset> SYSTEM
%token <pset> TEANL
%token <pset> TEXT
%token <pset> TICK
%token <pset> TICKLABEL
%token <pset> TICKMARKS
%token <pset> TICKP
%token <pset> TIDALCLOCK
%token <pset> TIDESTATION
%token <pset> TIME
%token <pset> TIMEINFO
%token <pset> TIMELINE
%token <pset> TITLE
%token <pset> TO
%token <pset> TOP
%token <pset> TOTAL
%token <pset> TRACK
%token <pset> TRANSECT
%token <pset> TRIANGLE1
%token <pset> TRIANGLE2
%token <pset> TRIANGLE3
%token <pset> TRIANGLE4
%token <pset> TRUEP
%token <pset> TYPE
%token <pset> UNITS
%token <pset> UP
%token <pset> VALUE
%token <pset> VECTOR
%token <pset> VEL
%token <pset> VELMARKER
%token <pset> VELOCITY
%token <pset> VERTICAL
%token <pset> VGAP
%token <pset> VIEW
%token <pset> VSCALE
%token <pset> VX1
%token <pset> VX2
%token <pset> VY1
%token <pset> VY2
%token <pset> WEEKS
%token <pset> WET
%token <pset> WETDRY
%token <pset> WIDTH
%token <pset> WIND
%token <pset> WITH
%token <pset> WORLD
%token <pset> WRAP
%token <pset> WRITE
%token <pset> WSCALE
%token <pset> WX1
%token <pset> WX2
%token <pset> WY1
%token <pset> WY2
%token <pset> X0
%token <pset> X1
%token <pset> X2
%token <pset> X3
%token <pset> X4
%token <pset> X5
%token <pset> XAXES
%token <pset> XAXIS
%token <pset> XCOR
%token <pset> XMAX
%token <pset> XMIN
%token <pset> XY
%token <pset> XYARC
%token <pset> XYBOX
%token <pset> XYDX
%token <pset> XYDXDX
%token <pset> XYDXDY
%token <pset> XYDY
%token <pset> XYDYDY
%token <pset> XYFIXED
%token <pset> XYHILO
%token <pset> XYRT
%token <pset> XYSEG
%token <pset> XYSTRING
%token <pset> XYUV
%token <pset> XYX2Y2
%token <pset> XYXX
%token <pset> XYYY
%token <pset> XYZ
%token <pset> XYZW
%token <pset> Y0
%token <pset> Y1
%token <pset> Y2
%token <pset> Y3
%token <pset> Y4
%token <pset> Y5
%token <pset> YAXES
%token <pset> YAXIS
%token <pset> YEARS
%token <pset> YES
%token <pset> YMAX
%token <pset> YMIN
%token <pset> ZEROXAXIS
%token <pset> ZEROYAXIS
%token <pset> ZOOM
%token <pset> ZOOMBOX

%type <ptr> asgn
%type <ptr> vasgn
%type <ptr> vexpr

%type <pset> onoff
%type <pset> flowonoff
%type <pset> worldview
%type <pset> extremetype
%type <pset> torf
%type <pset> filltype
%type <pset> opchoice
%type <pset> props
%type <pset> horv
%type <pset> units
%type <pset> timeunits
%type <pset> isolines
%type <pset> formatchoice
%type <pset> ticklabelattr
%type <pset> axislabeldesc
%type <pset> prop
%type <pset> sourcetype
%type <pset> justchoice
%type <pset> inoutchoice
%type <pset> signchoice
%type <pset> direction
%type <pset> graphtype
%type <pset> tickattr


%type <val> expr
%right '='
%left OR
%left AND
%nonassoc GT LT LE GE EQ NE
%left '+' '-'
%left '*' '/' '%'
%right '^'
%right UMINUS NOT

%%

list:
	| asgn '\n' {}
	| vasgn '\n' {}
	| expr '\n' { result = $1; }
	| vexpr '\n' { result = *$1; }
	| animation '\n'
	| annotation '\n'
	| models '\n'
	| grid '\n'
	| graph '\n'
	| setaxis '\n'
	| drogues'\n'
	| elevmarker '\n'
	| histboxes '\n'
	| zoomboxes '\n'
	| slice '\n'
	| isolines '\n' { }
	| tidalclock '\n'
	| timeline '\n'
	| timeinfo '\n'
	| mapscale '\n'
	| vscale '\n'
	| wscale '\n'
	| misc '\n'
	| error '\n' { return 1; }
	;

annotation:
	| CLEAR BOX { do_clear_boxes(); }
	| WITH BOX { curbox = next_box(); }
	| WITH BOX NUMBER { curbox = (int) $3; }
	| BOX onoff { boxes[curbox].active = $2; }
	| BOX GRAPHNO { boxes[curbox].gno = $2; }
	| BOX expr ',' expr ',' expr ',' expr
	{
	    if (curbox >= 0 && curbox < maxboxes) {
		boxes[curbox].x1 = $2;
		boxes[curbox].y1 = $4;
		boxes[curbox].x2 = $6;
		boxes[curbox].y2 = $8;
	    }
	}
	| BOX LOCTYPE worldview { sysbox.loctype = $3; }
	| BOX LINESTYLE NUMBER { sysbox.lines = (int) $3; }
	| BOX LINEWIDTH NUMBER { sysbox.linew = (int) $3; }
	| BOX COLOR NUMBER { sysbox.color = (int) $3; }
	| BOX FILL filltype { sysbox.fill = $3; }
	| BOX FILL COLOR NUMBER { sysbox.fillcolor = (int) $4; }
	| BOX FILL PATTERN NUMBER { sysbox.fillpattern = (int) $4; }
	| BOX DEF
	{
	    if (curbox >= 0 && curbox < maxboxes) {
                    boxes[curbox].loctype = sysbox.loctype;
                    boxes[curbox].color = sysbox.color;
                    boxes[curbox].linew = sysbox.linew;
                    boxes[curbox].lines = sysbox.lines;
                    boxes[curbox].fill = sysbox.fill;
                    boxes[curbox].fillcolor = sysbox.fillcolor;
                    boxes[curbox].fillpattern = sysbox.fillpattern;
	    }
	}
	| WITH LINE { curline = next_line(); }
	| WITH LINE NUMBER { curline = (int) $3; }
	| CLEAR LINE { do_clear_lines(); }
	| LINE onoff { lines[curline].active = $2; }
	| LINE GRAPHNO { lines[curline].gno = $2; }
	| LINE expr ',' expr ',' expr ',' expr
	{
	    lines[curline].x1 = $2;
	    lines[curline].y1 = $4;
	    lines[curline].x2 = $6;
	    lines[curline].y2 = $8;
	}
	| LINE LOCTYPE worldview { sysline.loctype = $3; }
	| LINE LINEWIDTH NUMBER { sysline.linew = (int) $3; }
	| LINE LINESTYLE NUMBER { sysline.lines = (int) $3; }
	| LINE COLOR NUMBER { sysline.color = (int) $3; }
	| LINE ARROW NUMBER { sysline.arrow = (int) $3; }
	| LINE ARROW SIZE NUMBER { sysline.asize = $4; }
	| LINE ARROW TYPE NUMBER { sysline.atype = (int) $4; }
	| LINE DEF
	{
	    if (curline >= 0 && curline < maxlines) {
		lines[curline].lines = sysline.lines;
		lines[curline].linew = sysline.linew;
		lines[curline].color = sysline.color;
		lines[curline].arrow = sysline.arrow;
		lines[curline].asize = sysline.asize;
		lines[curline].atype = sysline.atype;
		lines[curline].loctype = sysline.loctype;
	    }
	}
	| CLEAR STRING { do_clear_text(); }
	| WITH STRING { curstring = next_string(); }
	| WITH STRING NUMBER { curstring = (int) $3; }
	| STRING onoff { pstr[curstring].active = $2; }
	| STRING GRAPHNO { pstr[curstring].gno = $2; }
	| STRING expr ',' expr
	{
	    pstr[curstring].x = $2;
	    pstr[curstring].y = $4;
	}
	| STRING LOCTYPE worldview { sysstr.loctype = $3; }
	| STRING LINEWIDTH NUMBER { sysstr.linew = (int) $3; }
	| STRING COLOR NUMBER { sysstr.color = (int) $3; }
	| STRING ROT NUMBER { sysstr.rot = (int) $3; }
	| STRING FONTP NUMBER { sysstr.font = (int) $3; }
	| STRING JUST NUMBER { sysstr.just = (int) $3; }
	| STRING SYMBOL NUMBER { sysstr.sym = (int) $3; }
	| STRING SYMBOL LOCTYPE opchoice { sysstr.symloc = (int) $4; }
	| STRING SYMBOL SIZE NUMBER { sysstr.symsize = (double) $4; }
	| STRING SYMBOL FILL NUMBER { sysstr.symfill = (int) $4; }
	| STRING SYMBOL COLOR NUMBER { sysstr.symcolor = (int) $4; }
	| STRING CHAR SIZE NUMBER { sysstr.charsize = (double) $4; }
	| STRING DEF CHRSTR
	{
	    strcpy(pstr[curstring].s, (char *) $3);
	    pstr[curstring].linew = sysstr.linew;
	    pstr[curstring].color = sysstr.color;
	    pstr[curstring].font = sysstr.font;
	    pstr[curstring].just = sysstr.just;
	    pstr[curstring].sym = sysstr.sym;
	    pstr[curstring].symloc = sysstr.symloc;
	    pstr[curstring].symfill = sysstr.symfill;
	    pstr[curstring].symcolor = sysstr.symcolor;
	    pstr[curstring].symsize = sysstr.symsize;
	    pstr[curstring].loctype = sysstr.loctype;
	    pstr[curstring].rot = sysstr.rot;
	    pstr[curstring].charsize = sysstr.charsize;
	    /*print_plotstr(pstr[curstring]);*/
	}
	;

animation:
	STEP { setistep(); }
	| REVERSE { setreverse(); }
	| REWIND { setrewind(); }
	| FORWARD { setforward(); }
	| WRAP onoff { set_wrap($2 == ON); }
	| GOTO NUMBER { goto_step((int) $2 - 1); }
	| RUN { setirun(); }
	| RUN NUMBER ',' NUMBER { runsteps((int) $2, (int) $4); }
	| BATCH RUN NUMBER ',' NUMBER { batchrunsteps((int) $3, (int) $5, 1); }
	| BATCH RUN NUMBER ',' NUMBER SKIP NUMBER { batchrunsteps((int) $3, (int) $5, (int) $7); }
	| BATCH PREFIX CHRSTR { strcpy(batchprefix, (char *) $3); }
	| STOP { setistop(); }
	;

misc:
	SYSTEM CHRSTR { system($2); }
	| RUN BATCH CHRSTR
	{
	    gotbatch = 1;
	    batchfile[0] = 0;
	    strcpy(batchfile, $3);
	}
	| CHDIR CHRSTR {
		if (chdir((char *) $2) < 0) {
			sprintf(buf, "chdir() to %s failed", (char *) $2);
			errwin(buf);
		}
	}
	| RESET { setreset_world(); }
	| REDRAW { setredraw_world(); }
	| SLEEP NUMBER { sleep((int) $2); }
	| QUIT { exit(0); }
	| ZOOM expr ',' expr ',' expr ',' expr { my_blowup($2, $4, $6, $8); }
	| EXPAND { page(page_per, 4); }
	| SHRINK { page(page_per, 5); }
	| PAGE LEFT { page(page_per, 0); }
	| PAGE RIGHT { page(page_per, 1); }
	| PAGE UP { page(page_per, 2); }
	| PAGE DOWN { page(page_per, 3); }
	| PAGE expr { page_per = $2; }
	| PAGE INOUT NUMBER { scrollinout_proc((int) $3); }
	| LINK PAGE onoff { scrolling_islinked = $3 == ON; }
	| LOG CHRSTR
	{
	    if ((logfp = fopen($2, "w")) != NULL) {
		logfile = 1;
		printf("Opened logfile %s\n", $2);
	    } else {
		logfile = 0;
		printf("Failed to open logfile %s\n", $2);
	    }
	}
	| CLOSE LOG
	{
	    if (logfp != NULL) {
		printf("Closing logfile\n");
		logfile = 0;
		fclose(logfp);
		logfp = NULL;
	    }
	}
	| INCLUDE IMAGE CHRSTR { }
	| PRINT { batchrunstep(0); }
	| ECHO CHRSTR {
	    if (inwin) { set_left_footer($2); }
	    else { printf("%s\n", $2); }
	}
	| CMAP NUMBER ',' NUMBER ',' NUMBER ',' NUMBER
	{ set_colormapdata((int) $2, (int) $4, (int) $6, (int) $8); }
	| COLORMAP NUMBER ',' NUMBER ',' NUMBER ',' NUMBER
	{ set_colormapdata((int) $2, (int) $4, (int) $6, (int) $8); }
	| COLOR NUMBER ',' NUMBER ',' NUMBER ',' NUMBER
	{ set_colormapdata((int) $2, (int) $4, (int) $6, (int) $8); }
	| SCALAR { setisol = &(g[curg].salip); } isolines
	| MAGNITUDE { setisol = &(g[curg].velmagip); } isolines
	;

models:
	WITH TEANL NUMBER { 
		setflow = &g[curg].flowf[(int) $3]; 
	}
	| TEANL flowprops { } 
	| READ TEANL {}
        | WITH ADCIRC NUMBER { 
		setflow = &g[curg].flowt[(int) $3]; 
                elcirc_flowno = (int) $3;
	}
        | ADCIRC flowprops {}
	| READ ADCIRC ELEV NUMBER GRID CHRSTR CHRSTR
          {
	      int fno = (int) $4;
              readbin_adcirc_elev(fno, (char *) $6, (char *) $7);
              set_clock(0, flowt[fno].start, 
				flowt[fno].stop, flowt[fno].step,
                                  flowt[fno].nsteps);
              load_clock(ADCIRC, fno);
          }
	| WITH ELCIRC NUMBER { 
		setflow = &g[curg].flowt[(int) $3]; 
                elcirc_flowno = (int) $3;
		g[curg].curadc3d = (int) $3; 
	}
	| ELCIRC flowprops {}
	| ELCIRC RUN NUMBER { curadc3d = (int) $3; }
	| ELCIRC GRID NUMBER CHRSTR 
	{
		ReplaceElcircGrid((int) $3, (char *) $4);
	}
	| READ ELCIRC NUMBER REGION CHRSTR 
	{ 
		ReadElcircRegionFile((int) $3, (char *) $5);
	}
	| READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER LEVEL NUMBER
	{
		ReadElcirc((int) $3, (char *) $4, (int) $12 - 1, (int) $6, (int) $8, (int) $10, 0, 0, 0);
                elcirc_flowno = (int) $3;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER LEVEL NUMBER APPEND
	{
		ReadElcirc((int) $3, (char *) $4, (int) $12 - 1, (int) $6, (int) $8, (int) $10, 0, 0, 1);
                elcirc_flowno = (int) $3;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER DEPTH NUMBER
	{
		ReadElcircDepth((int) $3, (char *) $4, (char *) NULL, (double) $12, (int) $6, (int) $8, (int) $10, 0, 0, 0);
                elcirc_flowno = (int) $3;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER DEPTH NUMBER
	{
/* read at a given depth relative to the free surface */
		ReadElcircDepthFromFreeSurface((int) $4, (char *) $5, (char *) NULL, (double) $13, (int) $7, (int) $9, (int) $11, 0, 0, 0);
                elcirc_flowno = (int) $4;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER
	{
		ReadElcircSurf((int) $4, (char *) $5, 0, (int) $7, (int) $9, (int) $11, 0, 0.0, 0);
                elcirc_flowno = (int) $4;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER APPEND
	{
		ReadElcircSurf((int) $4, (char *) $5, 0, (int) $7, (int) $9, (int) $11, 0, 0.0, 1);
                elcirc_flowno = (int) $4;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC BOTTOM NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER
	{
		ReadElcircSurf((int) $4, (char *) $5, 1, (int) $7, (int) $9, (int) $11, 0, 0.0, 0);
                elcirc_flowno = (int) $4;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| READ ELCIRC BOTTOM NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER APPEND
	{
		ReadElcircSurf((int) $4, (char *) $5, 1, (int) $7, (int) $9, (int) $11, 0, 0.0, 1);
                elcirc_flowno = (int) $4;
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
	| WITH ELCIRC TRANSECT NUMBER { 
			curtrans = (int) $4; 
			settrans = &(trans[curtrans]); 
	}
	| SET ELCIRC TRANSECT NUMBER { curtrans = (int) $4; }
	| ELCIRC TRANSECT { setisol = &(g[curg].trans[curtrans].ip); } isolines
	| ELCIRC TRANSECT FLOW CHRSTR {strcpy(settrans->uvname, (char *) $4);}
	| ELCIRC TRANSECT VERTICAL FLOW CHRSTR {strcpy(settrans->vvname, (char *) $5);}
	| ELCIRC TRANSECT FLOW NUMBER { settrans->flowno = (int) $4;}
	| ELCIRC TRANSECT SALINITY CHRSTR {strcpy(settrans->salname, (char *) $4);}
	| ELCIRC TRANSECT ELEV CHRSTR {strcpy(settrans->elevname, (char *) $4);}
	| ELCIRC TRANSECT GRAPH NUMBER { settrans->gno = (int) $4;}
	| ELCIRC TRANSECT DISPLAY GRAPH NUMBER { settrans->transgno = (int) $5;}
	| ELCIRC TRANSECT DISPLAY LINE onoff { settrans->display = (int) $5;}
	| ELCIRC TRANSECT DISPLAY onoff { g[curg].trans[curtrans].display = (int) $4;}
	| ELCIRC TRANSECT DISPLAY MAG onoff { g[curg].trans[curtrans].display_mag = (int) $5;}
	| ELCIRC TRANSECT MAXLEVELS NUMBER { }
	| ELCIRC TRANSECT START NUMBER STOP NUMBER SKIP NUMBER 
	{
		settrans->start = (int) $4;
		settrans->stop = (int) $6;
		settrans->skip = (int) $8;
	}
	| ELCIRC TRANSECT SAMPLE NUMBER { settrans->npts = (int) $4; }
	| ELCIRC TRANSECT TYPE NUMBER { settrans->transtype = (int) $4; }
	| ELCIRC TRANSECT NUMBER ',' expr ',' expr { AddTransNXY(settrans, (int) $3, (double) $5, (double) $7); }
	| ELCIRC TRANSECT NODE NUMBER ',' NUMBER { AddTransNode(settrans, (int) $4, (int) $6); }
	| ELCIRC TRANSECT REGION FILEP CHRSTR 
	{
		settrans->transtype = 0;
		strcpy(settrans->transname, (char *) $5);
	}
	| ELCIRC TRANSECT LINE NUMBER ',' expr ',' expr ',' expr ',' expr 
	{
		settrans->transtype = 1;
		settrans->npts = (int) $4;
		settrans->x1 = (double) $6;
		settrans->y1 = (double) $8;
		settrans->x2 = (double) $10;
		settrans->y2 = (double) $12;
	}
	| READ ELCIRC TRANSECT 
	{ 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| READ ELCIRC TRANSECT APPEND
	{ 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| READ ELCIRC FLOW TRANSECT 
	{ 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| READ ELCIRC FLOW TRANSECT APPEND
	{ 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| READ ELCIRC SALINITY TRANSECT 
	{ 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| READ ELCIRC SALINITY TRANSECT APPEND
	{ 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
	| WITH ELCIRC MARKER NUMBER 
	{
		elcircmarker = (int) $4;
		setadc3d = &adc3d[(int) $4];
		setflow3d = &g[curg].flow3d[(int) $4];
	}
	| ELCIRC MARKER FLOW CHRSTR 
	{
		strcpy(setadc3d->datafile, (char *) $4);
	}
	| ELCIRC MARKER MAGNITUDE CHRSTR 
	{
		strcpy(setadc3d->datafile, (char *) $4);
	}
	| ELCIRC MARKER VECTOR CHRSTR 
	{
		strcpy(setadc3d->datafile, (char *) $4);
	}
	| ELCIRC MARKER SCALAR CHRSTR 
	{
		strcpy(setadc3d->datafile, (char *) $4);
	}
	| ELCIRC MARKER SALINITY CHRSTR 
	{
		strcpy(setadc3d->datafile, (char *) $4);
	}
	| ELCIRC MARKER ELEV CHRSTR 
	{
		strcpy(setadc3d->elevfile, (char *) $4);
	}
	| READ ELCIRC MARKER NODE NUMBER START NUMBER STEP NUMBER SKIP NUMBER 
	{
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) $5;
		ReadNodeDataNew(elcircmarker, (int) $7, (int) $9, 0);
	}
	| READ ELCIRC MARKER NODE NUMBER START NUMBER STEP NUMBER SKIP NUMBER APPEND
	{
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) $5;
		ReadNodeDataNew(elcircmarker, (int) $7, (int) $9, 1);
	}
	| READ ELCIRC MARKER XY NUMBER ',' NUMBER  START NUMBER STEP NUMBER SKIP NUMBER 
	{
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) $5;
		setadc3d->y = (double) $7;
		ReadXYDataNew(elcircmarker, (int) $9, (int) $11, 0);
	}
	| READ ELCIRC MARKER XY NUMBER ',' NUMBER  START NUMBER STEP NUMBER SKIP NUMBER APPEND
	{
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) $5;
		setadc3d->y = (double) $7;
		ReadXYDataNew(elcircmarker, (int) $9, (int) $11, 1);
	}
	| ELCIRC MARKER SCALAR { setisol = &(g[curg].salip); } isolines
	| ELCIRC MARKER SALINITY { setisol = &(g[curg].salip); } isolines
	| ELCIRC MARKER MAGNITUDE { setisol = &(g[curg].velmagip); } isolines
	| ELCIRC MARKER { setprops = &(setflow3d->p); } props
	| ELCIRC MARKER PREC NUMBER ',' NUMBER
	{
		setflow3d->precx = (int) $4;
		setflow3d->precy = (int) $6;
	}
	| ELCIRC MARKER ATTACH NUMBER { setflow3d->attach = (int) $4; }
	| ELCIRC MARKER LOCTYPE worldview { setflow3d->loctype = (int) $4; }
	| ELCIRC MARKER DISPLAY MARKER onoff { setflow3d->display_marker = (int) $5; }
	| ELCIRC MARKER DISPLAY onoff { setflow3d->display = (int) $4; }
	| ELCIRC MARKER COLOR NUMBER { setflow3d->p.color = $4; }
	| ELCIRC MARKER LINEWIDTH NUMBER { setflow3d->p.linew = $4; }
	| ELCIRC MARKER FILL COLOR NUMBER { setflow3d->p.fillcol = $5; }
	| ELCIRC MARKER WORLD expr ',' expr ',' expr ',' expr
	{ 
		setflow3d->wx1 = (double) $4; 
		setflow3d->wy1 = (double) $6; 
		setflow3d->wx2 = (double) $8; 
		setflow3d->wy2 = (double) $10; 
	}
	| ELCIRC MARKER VIEW expr ',' expr
	{ 
		setflow3d->vx = (double) $4; 
		setflow3d->vy = (double) $6; 
	}
	| ELCIRC MARKER LOC expr ',' expr
	{ 
		setflow3d->locx = (double) $4; 
		setflow3d->locy = (double) $6; 
	}
	| ELCIRC MARKER XY expr ',' expr
	{ 
		setflow3d->x = (double) $4; 
		setflow3d->y = (double) $6; 
	}
	;

flowprops:
	DISPLAY flowonoff { setflow->display = $2;  }
	| DISPLAY ELEV onoff { setflow->display_elev = $3; }
	| DISPLAY ELEV DEPTH onoff { setflow->display_elevdepth = $4; }
	| DISPLAY ELEV VALUE onoff { setflow->display_maxelevval = $4; }
	| DISPLAY ELEV MAXP onoff { setflow->display_maxelev = $4; }
	| DISPLAY ELEV AMP onoff { setflow->display_amp = $4; }
	| DISPLAY ELEV PHASE onoff { setflow->display_phase = $4; }
	| DISPLAY ELEV MARKERS onoff { setflow->display_elevmarkers = $4; }
	| DISPLAY FLOW MAG onoff { setflow->display_mag = $4; }
	| DISPLAY FLOW WIND onoff { setflow->display_wind = $4; }
	| DISPLAY INUNDATION onoff { setflow->display_inun = $2; }
	| COLOR NUMBER { setflow->p.color = $2; }
	| ELEV { setisol = &(setflow->elevip); } isolines
	| ELEV MAXP { setisol = &(setflow->maxelevip); } isolines
	| AMP { setisol = &(setflow->ampip); } isolines
	| PHASE { setisol = &(setflow->phaseip); } isolines
	| FLOW MAG { setisol = &(setflow->magip); } isolines
	| FLOW FREQ NUMBER { setflow->flowfreq = $3; }
	| FREQ NUMBER { setflow->freq = $2; }
	| ELEVMARKER NUMBER { setelevmarker = &(setflow->em[(int) $2]); }
	| SAMPLE FLOW onoff { setflow->sample = (int) $3; }
	| SAMPLE FLOW READ CHRSTR { ReadSampleFlow(setflow, (char *) $4); }
	| SAMPLE FLOW MINP NUMBER { SetMinSampleFlow(elcirc_flowno, (double) $4); }
	| SAMPLE FLOW READ XY CHRSTR { ReadSampleFlowXY(setflow, (char *) $5); }
	| SAMPLE FLOW TYPE XY { setflow->samptype = XY; }
	| SAMPLE FLOW TYPE NODE { setflow->samptype = NODE; }
	| SAMPLE FLOW NODE NUMBER { AddSampleFlowNode(setflow, (int) $4 - 1); }
	| SAMPLE FLOW ELEMENT NUMBER { AddSampleFlowElem(setflow, (int) $4 - 1); }
	| SAMPLE FLOW XY expr ',' expr { AddSampleFlowXY(setflow, (double) $4, (double) $6); }
	| DELETE SAMPLE FLOW NODE NUMBER { DeleteSampleFlowNode(setflow, (int) $5 - 1); }
	| DELETE SAMPLE FLOW ELEMENT NUMBER { DeleteSampleFlowElem(setflow, (int) $5 - 1); }
	| DELETE SAMPLE FLOW XY expr ',' expr { DeleteSampleFlowXY(setflow, (double) $5, (double) $7); }
	;

elevmarker:
	WITH ELEVMARKER NUMBER { setelevmarker = &(setflow->em[(int) $3]); }
	| ELEVMARKER ACTIVE onoff { setelevmarker->active = (int) $3; }
	| ELEVMARKER TYPE onoff { setelevmarker->type = (int) $3; }
	| ELEVMARKER DISPLAY onoff { setelevmarker->display = (int) $3; }
	| ELEVMARKER { setprops = &(setelevmarker->p); } props
	| ELEVMARKER LOCTYPE worldview { setelevmarker->loctype = (int) $3; }
	| ELEVMARKER NODE NUMBER { setelevmarker->node = (int) $3; }
	| ELEVMARKER LOC expr ',' expr
	{ 
		setelevmarker->locx = (double) $3; 
		setelevmarker->locy = (double) $5; 
	}
	| ELEVMARKER MINMAX expr ',' expr
	{ 
		setelevmarker->emin = (double) $3; 
		setelevmarker->emax = (double) $5; 
	}
	;

velocitymarker:
	VELOCITY MARKER ACTIVE onoff {}
	;

grid:
	GRID NUMBER { setgrid = &g[curg].grid[(int) $2]; }
	| WITH GRID NUMBER { setgrid = &g[curg].grid[(int) $3]; }
	| GRID DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display = $3; }
	| GRID BATH DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display_bath = $4; }
	| GRID COURANT DISPLAY onoff DT NUMBER { if (checkptr(setgrid, f_string)) setgrid->display_courant = $4; }
	| GRID COURANT DISPLAY VALUE onoff DT NUMBER { if (checkptr(setgrid, f_string)) setgrid->display_courantn = $5; }
	| GRID BOUNDARY DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display_boundary = $4; }
	| GRID { if (checkptr(setgrid, f_string)) setprops = &(setgrid->p); } props
	| GRID BOUNDARY { if (checkptr(setgrid, f_string)) setprops = &(setgrid->bp); } props
	| GRID NODES DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display_nodes = $4; }
	| GRID ELEMENTS DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display_elements = $4; }
	| GRID DEPTH DISPLAY onoff { if (checkptr(setgrid, f_string)) setgrid->display_depths = $4; }
	| GRID FILL onoff { if (checkptr(setgrid, f_string)) setgrid->display_gridf = $3; }
	| GRID BATH { if (checkptr(setgrid, f_string)) setisol = &(setgrid->ip); } isolines
	| GRAPHNO AUTOSCALE 
	{ 
		autoscale_grid((int) $1, g[(int) $1].curgrid); 
		set_defaults((int) $1);
	}
	| AUTOSCALE 
	{ 
		autoscale_grid(curg, g[curg].curgrid); 
		set_defaults(curg);
	}
	| READ GRID NUMBER CHRSTR 
	{
		extern int readgridfile;
		readgrid((int) $3, (char *) $4);
		readgridfile = 1;
	}
	;

drogues:
	WITH DROGUES NUMBER { setdrogs = &g[curg].drogues[(int) $3]; }
	| DROGUES DISPLAY onoff {  setdrogs->display = (int) $3; }
	| DROGUES { setprops = &(setdrogs->p); } props
	| READ DROGUES CHRSTR 
	{
    	if (!readdrogues(curdrog, (char *) $3, -1, 0, 0)) {
        	fprintf(stderr, "Error reading file %s", (char *) $3);
    	} else {
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
    	}
	}
	| READ DROGUES CHRSTR START NUMBER STOP NUMBER SKIP NUMBER 
	{
		readdrogues(0, (char *) $3, (int) $5, (int) $7, (int) $9);
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
	}
       ;

isolines:
	ISOLINES LEGEND onoff { setisol->lactive = $3; }
	| ISOLINES LEGEND LAYOUT horv { setisol->layout = $4; }
	| ISOLINES LEGEND LABEL onoff { setisol->llabels = $3; }
        | ISOLINES LEGEND FRAMEP onoff { setisol->frame = (int) $4; }
        | ISOLINES LEGEND FRAMEP COLOR NUMBER { setisol->framecol = (int) $5; }
        | ISOLINES NUMBER { setisol->nisol = (int) $2; }
        | ISOLINES TYPE NUMBER { setisol->type = (int) $3; }
        | ISOLINES SET TYPE NUMBER { setisol->isoltype = (int) $4; }
        | ISOLINES FILL TYPE NUMBER { setisol->visflag = (int) $4; }
        | ISOLINES SAVE ON { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
        | ISOLINES SAVE OFF { 
		int i;
		setisol->writeflag = 0; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
        | ISOLINES SAVE FILEP CHRSTR { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		strncpy(setisol->wname, (char *) $4, 1023);
	}
        | ISOLINES SAVE ISOLINE NUMBER FILEP CHRSTR { 
		int i;
		setisol->writeflag = 1; 
		strncpy(setisol->wname, (char *) $6, 1023);
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
    		setisol->writelevel[(int) $4] = 1;
	}
        | ISOLINES SAVE ISOLINE NUMBER { 
    		setisol->writelevel[(int) $4] = 1;
	}
        | ISOLINES START expr STEP expr {
            setisol->cis[0] = (double) $3;
            setisol->cint = (double) $5;
	    if (setisol->isoltype == 0) {
		int i;
		for (i=1;i< 16;i++) {
		    setisol->cis[i] = setisol->cis[0] + i * setisol->cint;
		}
	    }
        } 
        | ISOLINES MINMAX expr ',' expr { setisol->cmin = $3; setisol->cmax = $5; }
        | ISOLINE expr ',' expr { setisol->cis[(int) $2] = $4; }
        | ISOLINES LEGEND LOCTYPE worldview { setisol->loctype = $4; }
        | ISOLINES LEGEND expr ',' expr {
            setisol->x = (double) $3;
            setisol->y = (double) $5;
        } 
	| ISOLINE NUMBER COLOR NUMBER { setisol->color[(int) $2] = (int) $4; }
        | ISOLINE NUMBER LINEWIDTH NUMBER { setisol->linew[(int) $2] = (int) $4; }
        | ISOLINE NUMBER LINESTYLE NUMBER { setisol->lines[(int) $2] = (int) $4; }
	| ISOLINES { setprops = &(setisol->p); }  props
        | ISOLINES LEGEND { setprops = &(setisol->p); } props
        | ISOLINES LEGEND SIZE NUMBER ',' NUMBER {
            setisol->xlen = $4;
            setisol->ylen = $6;
        }
        | ISOLINES LEGEND HGAP NUMBER ',' NUMBER {
            setisol->xgap = $4;
            setisol->ygap = $6;
        }
        ;

histboxes:
	WITH HISTBOX NUMBER { sethistbox = &g[curg].hbox[(int) $3]; }
	| HISTBOX { setprops = &(sethistbox->p); } props
	| HISTBOX PREC NUMBER ',' NUMBER
	{
		sethistbox->precx = (int) $3;
		sethistbox->precy = (int) $5;
	}
	| HISTBOX ATTACH NUMBER { sethistbox->attach = (int) $3; }
	| HISTBOX TICKP NUMBER ',' NUMBER { sethistbox->xtickm = (double ) $3; sethistbox->ytickm = (double ) $5;}
	| HISTBOX LOCTYPE worldview { sethistbox->loctype = (int) $3; }
	| HISTBOX DISPLAY MARKER onoff { sethistbox->display_marker = (int) $4; }
	| HISTBOX DISPLAY onoff { sethistbox->display = (int) $3; }
	| HISTBOX DISPLAY ADCIRC NUMBER torf { sethistbox->adcirc[(int) $4] = (int) $5 == TRUEP; }
	| HISTBOX DISPLAY ADCIRC NUMBER COLOR NUMBER { sethistbox->ap[(int) $4].color = (int) $6; }
	| HISTBOX COLOR NUMBER { sethistbox->p.color = $3; }
	| HISTBOX LINEWIDTH NUMBER { sethistbox->p.linew = $3; }
	| HISTBOX FILL COLOR NUMBER { sethistbox->p.fillcol = $4; }
	| HISTBOX READ NUMBER CHRSTR { read_hist((int) $3, TIME, (char *) $4); }
	| HISTBOX DISPLAY HISTORY torf { sethistbox->thist = (int) $4 == TRUEP; }
	| HISTBOX DISPLAY HISTORY COLOR NUMBER { sethistbox->hp.color = (int) $5; }
	| HISTBOX WORLD expr ',' expr ',' expr ',' expr
	{ 
		sethistbox->wx1 = (double) $3; 
		sethistbox->wy1 = (double) $5; 
		sethistbox->wx2 = (double) $7; 
		sethistbox->wy2 = (double) $9; 
	}
	| HISTBOX VIEW expr ',' expr
	{ 
		sethistbox->vx = (double) $3; 
		sethistbox->vy = (double) $5; 
	}
	| HISTBOX LOC expr ',' expr
	{ 
		sethistbox->locx = (double) $3; 
		sethistbox->locy = (double) $5; 
	}
	| HISTBOX XY expr ',' expr
	{ 
		sethistbox->x = (double) $3; 
		sethistbox->y = (double) $5; 
	}
	;

zoomboxes:
	WITH ZOOMBOX NUMBER { setzoombox = &g[curg].zbox[(int) $3]; }
	| ZOOMBOX { setprops = &(setzoombox->p); } props
	| ZOOMBOX PREC NUMBER ',' NUMBER
	{
		setzoombox->precx = (int) $3;
		setzoombox->precy = (int) $5;
	}
	| ZOOMBOX ATTACH NUMBER { setzoombox->attach = (int) $3; }
	| ZOOMBOX LOCTYPE worldview { setzoombox->loctype = (int) $3; }
	| ZOOMBOX DISPLAY MARKER onoff { setzoombox->display_marker = (int) $4; }
	| ZOOMBOX onoff { setzoombox->active = (int) $2; }
	| ZOOMBOX DISPLAY onoff { setzoombox->display = (int) $3; }
	| ZOOMBOX ZOOM NUMBER { setzoombox->expand = (int) $3; }
	| ZOOMBOX SCALE NUMBER { setzoombox->expand = (int) $3; }
	| ZOOMBOX COLOR NUMBER { setzoombox->p.color = $3; }
	| ZOOMBOX LINEWIDTH NUMBER { setzoombox->p.linew = $3; }
	| ZOOMBOX FILL COLOR NUMBER { setzoombox->p.fillcol = $4; }
	| ZOOMBOX WORLD expr ',' expr ',' expr ',' expr
	{ 
		setzoombox->wx1 = (double) $3; 
		setzoombox->wy1 = (double) $5; 
		setzoombox->wx2 = (double) $7; 
		setzoombox->wy2 = (double) $9; 
	}
	| ZOOMBOX VIEW expr ',' expr
	{ 
		setzoombox->vx = (double) $3; 
		setzoombox->vy = (double) $5; 
	}
	| ZOOMBOX LOC expr ',' expr
	{ 
		setzoombox->locx = (double) $3; 
		setzoombox->locy = (double) $5; 
	}
	| ZOOMBOX XY expr ',' expr
	{ 
		setzoombox->x = (double) $3; 
		setzoombox->y = (double) $5; 
	}
	;

wscale:
	WSCALE onoff { g[curg].wl.active = $2; }
	| WSCALE LENGTH NUMBER { g[curg].wl.len = $3; }
	| WSCALE SCALE NUMBER { g[curg].wl.scale = $3; }
	| WSCALE COLOR NUMBER { g[curg].wl.p.color = $3; }
	| WSCALE LOCTYPE worldview { g[curg].wl.loctype = $3; }
	| WSCALE LOC expr ',' expr
	{ 
		g[curg].wl.x = $3;
		g[curg].wl.y = $5;
	}
	| WSCALE UNITS units 
	{ 
		g[curg].wl.units = $3;
		switch (g[curg].wl.units) {
		case MM: g[curg].wl.unitfac = 0.001; break;
		case CM: g[curg].wl.unitfac = 0.01; break;
		case M: g[curg].wl.unitfac = 1.0; break;
		case KM: g[curg].wl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
	| WSCALE props { }
	;


vscale:
	VSCALE onoff { g[curg].vl.active = $2; }
	| VSCALE LENGTH NUMBER { g[curg].vl.len = $3; }
	| VSCALE SCALE NUMBER { g[curg].vl.scale = $3; }
	| VSCALE COLOR NUMBER { g[curg].vl.p.color = $3; }
	| VSCALE LOCTYPE worldview { g[curg].vl.loctype = $3; }
	| VSCALE LOC expr ',' expr
	{ 
		g[curg].vl.x = $3;
		g[curg].vl.y = $5;
	}
	| VSCALE UNITS units 
	{ 
		g[curg].vl.units = $3;
		switch (g[curg].vl.units) {
		case MM: g[curg].vl.unitfac = 0.001; break;
		case CM: g[curg].vl.unitfac = 0.01; break;
		case M: g[curg].vl.unitfac = 1.0; break;
		case KM: g[curg].vl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
	| VSCALE props { }
	;

mapscale:
	MAPSCALE onoff { g[curg].mapscale.active = $2; }
	| MAPSCALE LENGTH NUMBER { g[curg].mapscale.len = $3; }
	| MAPSCALE COLOR NUMBER { g[curg].mapscale.p.color = $3; }
	| MAPSCALE SCALE NUMBER { g[curg].mapscale.scale = $3; }
	| MAPSCALE LOCTYPE worldview { g[curg].mapscale.loctype = $3; }
	| MAPSCALE LOC expr ',' expr
	{ 
		g[curg].mapscale.x = $3;
		g[curg].mapscale.y = $5;
	}
	| MAPSCALE UNITS units 
	{ 
		g[curg].mapscale.units = $3;
		switch (g[curg].mapscale.units) {
		case MM: g[curg].mapscale.unitfac = 0.001; break;
		case CM: g[curg].mapscale.unitfac = 0.01; break;
		case M: g[curg].mapscale.unitfac = 1.0; break;
		case KM: g[curg].mapscale.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for mapscape scale\n"); break;
		}
	}
	| MAPSCALE props { }
	;

tidalclock:
	TIDALCLOCK onoff { g[curg].tidalclock.active = $2; }
	| TIDALCLOCK COLOR NUMBER { g[curg].tidalclock.p.color = $3; }
	| TIDALCLOCK FILL COLOR NUMBER { g[curg].tidalclock.p.fillcol = $4; }
	| TIDALCLOCK TOTAL TIME NUMBER { g[curg].tidalclock.total_time = $4; }
	| TIDALCLOCK LOCTYPE worldview { g[curg].tidalclock.loctype = $3; }
	| TIDALCLOCK LOC expr ',' expr
	{ 
		g[curg].tidalclock.x = $3;
		g[curg].tidalclock.y = $5;
	}
	| TIDALCLOCK props { }
	;

timeinfo:
	TIMEINFO onoff { g[curg].timeinfo.active = $2; }
	| TIMEINFO START CHRSTR 
	{ 
		strcpy(g[curg].timeinfo.start, (char *) $3); 
		time_info_start(curg);
	}
	| TIMEINFO expr ',' expr
	{
	    g[curg].timeinfo.x = $2;
	    g[curg].timeinfo.y = $4;
	}
	| TIMEINFO LOCTYPE worldview { g[curg].timeinfo.loctype = $3; }
	| TIMEINFO LINEWIDTH NUMBER { g[curg].timeinfo.linew = (int) $3; }
	| TIMEINFO COLOR NUMBER { g[curg].timeinfo.color = (int) $3; }
	| TIMEINFO ROT NUMBER { g[curg].timeinfo.rot = (int) $3; }
	| TIMEINFO FONTP NUMBER { g[curg].timeinfo.font = (int) $3; }
	| TIMEINFO JUST NUMBER { g[curg].timeinfo.just = (int) $3; }
	| TIMEINFO CHAR SIZE NUMBER { g[curg].timeinfo.charsize = (double) $4; }
	| TIMEINFO FORMAT CHRSTR { strcpy(g[curg].timeinfo.format, (char *) $3); }
	;

timeline:
	TIMELINE onoff { g[curg].timeline.active = (int) $2; }
	| TIMELINE LENGTH NUMBER { g[curg].timeline.len = (int) $3; }
	| TIMELINE WIDTH NUMBER { g[curg].timeline.width = (int) $3; }
	| TIMELINE START NUMBER { g[curg].timeline.start = (double) $3; }
	| TIMELINE STOP NUMBER { g[curg].timeline.stop = (double) $3; }
	| TIMELINE STEP NUMBER { g[curg].timeline.step = (double) $3; }
	| TIMELINE PREC NUMBER { g[curg].timeline.p.prec = (int) $3; }
	| TIMELINE UNITS NUMBER { g[curg].timeline.units = (int) $3; }
	| TIMELINE COLOR NUMBER { g[curg].timeline.c1 = g[curg].timeline.c3 = (int) $3; }
	| TIMELINE FILL COLOR NUMBER { g[curg].timeline.c2 = (int) $4; }
	| TIMELINE LOCTYPE worldview { g[curg].timeline.loctype = (int) $3; }
	| TIMELINE LOC expr ',' expr
	{ 
		g[curg].timeline.x = (double) $3;
		g[curg].timeline.y = (double) $5;
	}
	| TIMELINE props { }
	;

slice:
	WITH SLICE NUMBER { setslice = &g[curg].sbox[(int) $3]; }
	| SLICE { setprops = &(setslice->p); } props
	| SLICE PREC NUMBER ',' NUMBER
	{
		setslice->precx = (int) $3;
		setslice->precy = (int) $5;
	}
	| SLICE ATTACH NUMBER { setslice->attach = (int) $3; }
	| SLICE LOCTYPE worldview { setslice->loctype = (int) $3; }
	| SLICE DISPLAY MARKER onoff { setslice->display_marker = (int) $4; }
	| SLICE onoff { setslice->active = (int) $2; }
	| SLICE DISPLAY onoff { setslice->display = (int) $3; }
	| SLICE COLOR NUMBER { setslice->p.color = $3; }
	| SLICE LINEWIDTH NUMBER { setslice->p.linew = $3; }
	| SLICE FILL COLOR NUMBER { setslice->p.fillcol = $4; }
	| SLICE WORLD expr ',' expr ',' expr ',' expr
	{ 
		setslice->wx1 = (double) $3; 
		setslice->wy1 = (double) $5; 
		setslice->wx2 = (double) $7; 
		setslice->wy2 = (double) $9; 
	}
	| SLICE VIEW expr ',' expr
	{ 
		setslice->vx = (double) $3; 
		setslice->vy = (double) $5; 
	}
	| SLICE LOC expr ',' expr
	{ 
		setslice->locx = (double) $3; 
		setslice->locy = (double) $5; 
	}
	| SLICE XY expr ',' expr
	{ 
		setslice->x = (double) $3; 
		setslice->y = (double) $5; 
	}
	;

props:
	PROP COLOR NUMBER { setprops->color = $3; }
	| PROP LINEWIDTH NUMBER { setprops->linew = $3; }
	| PROP LINESTYLE NUMBER { setprops->lines = $3; }
	| PROP FORMAT formatchoice { setprops->format = $3; }
	| PROP FONTP NUMBER { setprops->font = $3; }
	| PROP PREC NUMBER { setprops->prec = $3; }
	| PROP CHAR SIZE NUMBER { setprops->charsize = $4; }
	| PROP SYMBOL NUMBER { setprops->symbol = $3; }
	| PROP SYMBOL SIZE NUMBER { setprops->symsize = $4; }
	| PROP FILL onoff { setprops->fill = $3; }
	| PROP FILL filltype { setprops->fillusing = $3; }
	| PROP FILL COLOR NUMBER { setprops->fillcol = $4; }
	| PROP FILL PATTERN NUMBER { setprops->fillpat = $4; }
	| PROP ARROW NUMBER { setprops->arrow = $3; }
	| PROP ARROW TYPE NUMBER { setprops->atype = $4; }
	| PROP ARROW SIZE NUMBER { setprops->asize = $4; }
	;

graph:
	WITH GRAPHNO { curg = (int) $2; }
	| WITH GRAPH NUMBER { curg = (int) $3; }
	| KILL GRAPHNO { kill_graph($2); }
	| KILL GRAPHS { kill_graph(maxgraph); }
	| LOCATOR onoff
	{
	    extern int go_locateflag;
	    go_locateflag = ($2 == ON);
	}
	| FOCUS GRAPHNO
	{
	    cg = curg = (int) $2;
	    draw_focus(curg);
	    defineworld(g[curg].w.xg1, g[curg].w.yg1, g[curg].w.xg2, g[curg].w.yg2, 
			islogx(curg), islogy(curg));
	    viewport(g[curg].v.xv1, g[curg].v.yv1, g[curg].v.xv2, g[curg].v.yv2);
	    draw_focus(curg);
	    update_all(curg);
	}
	| FOCUS onoff { draw_focus_flag = $2; }
	| FOCUS SET { focus_policy = $2; }
	| FOCUS FOLLOWS { focus_policy = $2; }
	| FOCUS CLICK { focus_policy = $2; }
	| SOURCE sourcetype { cursource = $2; }
	| PUSH { push_world(); }
	| POP { pop_world(); }
	| CYCLE { cycle_world_stack(); }
	| STACK NUMBER {
	    if ((int) $2 > 0)
		show_world_stack((int) $2 - 1);
	}
	| STACK WORLD expr ',' expr ',' expr ',' expr TICKP expr ',' expr ',' expr ',' expr
	{
	    add_world(curg, $3, $5, $7, $9, $11, $13, $15, $17);
	}
	| CLEAR STACK { clear_world_stack(); }
	| WORLD expr ',' expr ',' expr ',' expr
	{
	    g[curg].w.xg1 = $2;
	    g[curg].w.yg1 = $4;
	    g[curg].w.xg2 = $6;
	    g[curg].w.yg2 = $8;
	}
	| WORLD XMIN expr { g[curg].w.xg1 = $3; }
	| WORLD XMAX expr { g[curg].w.xg2 = $3; }
	| WORLD YMIN expr { g[curg].w.yg1 = $3; }
	| WORLD YMAX expr { g[curg].w.yg2 = $3; }
	| VIEW expr ',' expr ',' expr ',' expr
	{
	    g[curg].v.xv1 = $2;
	    g[curg].v.yv1 = $4;
	    g[curg].v.xv2 = $6;
	    g[curg].v.yv2 = $8;
	}
	| VIEW XMIN NUMBER { g[curg].v.xv1 = $3; }
	| VIEW XMAX NUMBER { g[curg].v.xv2 = $3; }
	| VIEW YMIN NUMBER { g[curg].v.yv1 = $3; }
	| VIEW YMAX NUMBER { g[curg].v.yv2 = $3; }
	| TITLE CHRSTR {
	    set_plotstr_string(&g[curg].labs.title, (char *) $2);
	}
	| TITLE FONTP NUMBER {
	    g[curg].labs.title.font = checkon(FONTP, g[curg].labs.title.font, (int) $3);
	}
	| TITLE SIZE NUMBER { g[curg].labs.title.charsize = $3; }
	| TITLE COLOR NUMBER {
	    g[curg].labs.title.color = checkon(COLOR, g[curg].labs.title.color, (int) $3);
	}
	| TITLE LINEWIDTH NUMBER
	{
	    g[curg].labs.title.linew = checkon(LINEWIDTH, g[curg].labs.title.linew, (int) $3);
	}
	| SUBTITLE CHRSTR {
	    set_plotstr_string(&g[curg].labs.stitle, (char *) $2);
	}
	| SUBTITLE FONTP NUMBER {
	    g[curg].labs.stitle.font = checkon(FONTP, g[curg].labs.stitle.font, (int) $3);
	}
	| SUBTITLE SIZE NUMBER { g[curg].labs.stitle.charsize = $3; }
	| SUBTITLE COLOR NUMBER
	{
	    g[curg].labs.stitle.color = checkon(COLOR, g[curg].labs.stitle.color, (int) $3);
	}
	| SUBTITLE LINEWIDTH NUMBER
	{
	    g[curg].labs.stitle.linew = checkon(LINEWIDTH, g[curg].labs.stitle.color, (int) $3);
	}
	| FRAMEP onoff { g[curg].f.active = $2; }
	| FRAMEP TYPE NUMBER { g[curg].f.type = (int) $3; }
	| FRAMEP LINESTYLE NUMBER { g[curg].f.lines = checkon(LINESTYLE, g[curg].f.lines, (int) $3); }
	| FRAMEP LINEWIDTH NUMBER { g[curg].f.linew = checkon(LINEWIDTH, g[curg].f.linew, (int) $3); }
	| FRAMEP COLOR NUMBER { g[curg].f.color = checkon(COLOR, g[curg].f.color, (int) $3); }
	| FRAMEP FILL onoff { g[curg].f.fillbg = $3; }
	| FRAMEP BACKGROUND COLOR NUMBER { g[curg].f.bgcolor = (int) $4; }
	| GRAPHNO onoff { g[$1].active = $2; }
	| GRAPHNO LABEL onoff { g[$1].label = $3; }
	| GRAPHNO AUTOSCALE TYPE AUTO { g[$1].auto_type = $4; }
	| GRAPHNO AUTOSCALE TYPE SPEC { g[$1].auto_type = $4; }
	| GRAPHNO AUTOSCALE torf { g[$1].parmsread = ($3 == FALSEP); }
	| GRAPHNO HIDDEN torf { g[$1].hidden = ($3 == TRUEP); }
	| GRAPHNO TYPE graphtype { g[$1].type = $3; }
	| GRAPHNO FIXEDPOINT onoff { g[$1].pointset = ($3 == ON); }
	| GRAPHNO FIXEDPOINT FORMAT formatchoice formatchoice
	{
	    g[$1].fx = $4;
	    g[$1].fy = $5;
	}
	| GRAPHNO FIXEDPOINT PREC NUMBER ',' NUMBER
	{
	    g[$1].px = $4;
	    g[$1].py = $6;
	}
	| GRAPHNO FIXEDPOINT XY expr ',' expr
	{
	    g[$1].dsx = $4;
	    g[$1].dsy = $6;
	}
	| GRAPHNO FIXEDPOINT TYPE NUMBER { g[$1].pt_type = (int) $4; }
	;

setaxis:
	axis axisfeature
	|  allaxes
	|  GRAPHS axis {}
	|  GRAPHS axis axisfeature {}
	|  GRAPHS allaxes {}
	;

axis:
	XAXIS {}
	|  YAXIS {}
	|  ALTXAXIS {}
	|  ALTYAXIS {}
	|  ZEROXAXIS {}
	|  ZEROYAXIS {}
	;

allaxes:
	AXES axesprops {}
	|  XAXES axesprops {}
	|  YAXES axesprops {}
	;

axesprops:
	onoff { set_axis_prop(whichgraph, naxis, $1, 0.0); }
	| COLOR NUMBER { set_axis_prop(whichgraph, naxis, $1, $2); }
	| LINEWIDTH NUMBER { set_axis_prop(whichgraph, naxis, $1, $2); }
	| LINESTYLE NUMBER { set_axis_prop(whichgraph, naxis, $1, $2); }
	| FONTP NUMBER { set_axis_prop(whichgraph, naxis, $1, $2); }
	| CHAR SIZE NUMBER { set_axis_prop(whichgraph, naxis, $1, $3); }
	| GRID onoff { set_axis_prop(whichgraph, naxis, $1, $2); }
	;

/*	TICKP tickdesc {}*/
axisfeature:
	TICKP tickattr {}
	|  TICKLABEL ticklabeldesc {}
	|  LABEL axislabeldesc {}
	|  BAR axisbardesc {}
	|  onoff { g[curg].t[naxis].active = $1; }
	;

tickdesc:
	tickattr {}
	|  tickdesc tickattr {}
	;

tickattr:
	onoff
	{
	    g[curg].t[naxis].t_flag = $1;
	    g[curg].t[naxis].t_mflag = $1;
	}
	| MAJOR onoff { g[curg].t[naxis].t_flag = $2; }
	| MINOR onoff { g[curg].t[naxis].t_mflag = $2; }
	| MAJOR expr { g[curg].t[naxis].tmajor = $2; }
	| MINOR expr { g[curg].t[naxis].tminor = $2; }
	| OFFSETX expr { g[curg].t[naxis].offsx = $2; }
	| OFFSETY expr { g[curg].t[naxis].offsy = $2; }
	| ALT onoff { g[curg].t[naxis].alt = $2; }
	| MINP expr { g[curg].t[naxis].tmin = $2; }
	| MAXP expr { g[curg].t[naxis].tmax = $2; }
	| DEFAULT NUMBER { g[curg].t[naxis].t_num = (int) $2; }
	| inoutchoice { g[curg].t[naxis].t_inout = $1; }
	| LOG onoff { g[curg].t[naxis].t_log = $2; }
	| SIZE NUMBER { g[curg].t[naxis].t_size = $2; }
	| MAJOR SIZE NUMBER { g[curg].t[naxis].t_size = $3; }
	| MINOR SIZE NUMBER { g[curg].t[naxis].t_msize = $3; }
	| COLOR NUMBER { g[curg].t[naxis].t_color = g[curg].t[naxis].t_mcolor = (int) $2; }
	| LINEWIDTH NUMBER { g[curg].t[naxis].t_linew = g[curg].t[naxis].t_mlinew = (int) $2; }
	| MAJOR COLOR NUMBER { g[curg].t[naxis].t_color = (int) $3; }
	| MINOR COLOR NUMBER { g[curg].t[naxis].t_mcolor = (int) $3; }
	| MAJOR LINEWIDTH NUMBER { g[curg].t[naxis].t_linew = (int) $3; }
	| MINOR LINEWIDTH NUMBER { g[curg].t[naxis].t_mlinew = (int) $3; }
	| MAJOR LINESTYLE NUMBER { g[curg].t[naxis].t_lines = (int) $3; }
	| MINOR LINESTYLE NUMBER { g[curg].t[naxis].t_mlines = (int) $3; }
	| MAJOR GRID onoff { g[curg].t[naxis].t_gridflag = $3; }
	| MINOR GRID onoff { g[curg].t[naxis].t_mgridflag = $3; }
	| OP opchoice { g[curg].t[naxis].t_op = $2; }
	| TYPE AUTO { g[curg].t[naxis].t_type = AUTO; }
	| TYPE SPEC { g[curg].t[naxis].t_type = SPEC; }
	| SPEC NUMBER { g[curg].t[naxis].t_spec = (int) $2; }
	| NUMBER ',' expr { g[curg].t[naxis].t_specloc[(int) $1] = $3; }
	;

ticklabeldesc:
	ticklabelattr {}
	|  ticklabeldesc ticklabelattr {}
	;

ticklabelattr:
	onoff { g[curg].t[naxis].tl_flag = $1; }
	| TYPE AUTO { g[curg].t[naxis].tl_type = AUTO; }
	| TYPE SPEC { g[curg].t[naxis].tl_type = SPEC; }
	| PREC NUMBER { g[curg].t[naxis].tl_prec = (int) $2; }
	| FORMAT formatchoice { g[curg].t[naxis].tl_format = $2; }
	| FORMAT NUMBER { g[curg].t[naxis].tl_format = $2; }
	| APPEND CHRSTR { strcpy(g[curg].t[naxis].tl_appstr, (char *) $2); }
	| PREPEND CHRSTR { strcpy(g[curg].t[naxis].tl_prestr, (char *) $2); }
	| LAYOUT HORIZONTAL { g[curg].t[naxis].tl_layout = HORIZONTAL; }
	| LAYOUT VERTICAL { g[curg].t[naxis].tl_layout = VERTICAL; }
	| LAYOUT SPEC { g[curg].t[naxis].tl_layout = SPEC; }
	| ANGLE NUMBER { g[curg].t[naxis].tl_angle = (int) $2; }
	| JUST justchoice { g[curg].t[naxis].tl_just = (int) $2; }
	| SKIP NUMBER { g[curg].t[naxis].tl_skip = (int) $2; }
	| STAGGER NUMBER { g[curg].t[naxis].tl_staggered = (int) $2; }
	| OP opchoice { g[curg].t[naxis].tl_op = $2; }
	| SIGN signchoice { g[curg].t[naxis].tl_sign = $2; }
	| START expr { g[curg].t[naxis].tl_start = $2; }
	| STOP expr { g[curg].t[naxis].tl_stop = $2; }
	| START TYPE SPEC { g[curg].t[naxis].tl_starttype = (int) $3; }
	| START TYPE AUTO { g[curg].t[naxis].tl_starttype = (int) $3; }
	| STOP TYPE SPEC { g[curg].t[naxis].tl_stoptype = (int) $3; }
	| STOP TYPE AUTO { g[curg].t[naxis].tl_stoptype = (int) $3; }
	| VGAP NUMBER { g[curg].t[naxis].tl_vgap = $2; }
	| HGAP NUMBER { g[curg].t[naxis].tl_hgap = $2; }
	| CHAR SIZE NUMBER { g[curg].t[naxis].tl_charsize = $3; }
	| FONTP NUMBER { g[curg].t[naxis].tl_font = (int) $2; }
	| COLOR NUMBER { g[curg].t[naxis].tl_color = (int) $2; }
	| LINEWIDTH NUMBER { g[curg].t[naxis].tl_linew = (int) $2; }
	| NUMBER ',' CHRSTR { set_plotstr_string(&g[curg].t[naxis].t_speclab[(int) $1], (char *) $3); }
	;

axislabeldesc:
	CHRSTR { set_plotstr_string(&g[curg].t[naxis].label, (char *) $1); }
	| LAYOUT PERP { g[curg].t[naxis].label_layout = PERP; }
	| LAYOUT PARA { g[curg].t[naxis].label_layout = PARA; }
	| PLACE AUTO { g[curg].t[naxis].label_place = $2; }
	| PLACE SPEC { g[curg].t[naxis].label_place = $2; }
	| PLACE NUMBER ',' NUMBER {
	    g[curg].t[naxis].label.x = $2;
	    g[curg].t[naxis].label.y = $4;
	}
	| JUST justchoice { g[curg].t[naxis].label.just = (int) $2; }
	| CHAR SIZE NUMBER { g[curg].t[naxis].label.charsize = $3; }
	| FONTP NUMBER { g[curg].t[naxis].label.font = (int) $2; }
	| COLOR NUMBER { g[curg].t[naxis].label.color = (int) $2; }
	| LINEWIDTH NUMBER { g[curg].t[naxis].label.linew = (int) $2; }
	;

axisbardesc:
	onoff { g[curg].t[naxis].t_drawbar = $1; }
	| COLOR NUMBER { g[curg].t[naxis].t_drawbarcolor = (int) $2; }
	| LINESTYLE NUMBER { g[curg].t[naxis].t_drawbarlines = (int) $2; }
	| LINEWIDTH NUMBER { g[curg].t[naxis].t_drawbarlinew = (int) $2; }
	;

prop:
        LINESTYLE { $$ = $1; }
        | LINEWIDTH { $$ = $1; }
        | FONTP { $$ = $1; }
        | COLOR { $$ = $1; }
        | SIZE { $$ = $1; }
        ;

sourcetype:
	DISK { $$ = DISK; }
	| PIPE { $$ = PIPE; }
	;

justchoice:
	RIGHT { $$ = RIGHT; }
	| LEFT { $$ = LEFT; }
	| CENTER { $$ = CENTER; }
	;

extremetype:
	MINP { $$ = MINP; }
	| MAXP { $$ = MAXP; }
	;

graphtype:
        XY { $$ = $1; }
        | LOGX { $$ = $1; }
        | LOGY { $$ = $1; }
        | LOGXY { $$ = $1; }
        | XYFIXED { $$ = XYFIXED; }
        ;

inoutchoice:
	IN { $$ = IN; }
	| OUT { $$ = OUT; }
	| BOTH { $$ = BOTH; }
	;

signchoice:
	NORMAL { $$ = NORMAL; }
	| ABSOLUTE { $$ = ABSOLUTE; }
	| NEGATE { $$ = NEGATE; }
	;

direction:
	UP { $$ = UP; }
	| DOWN { $$ = DOWN; }
	| RIGHT { $$ = RIGHT; }
	| LEFT { $$ = LEFT; }
	| IN { $$ = IN; }
	| OUT { $$ = OUT; }
	;

formatchoice:
        DECIMAL { $$ = DECIMAL; }
        | EXPONENTIAL { $$ = EXPONENTIAL; }
        | POWER { $$ = POWER; }
        | GENERAL { $$ = GENERAL; }
        | DDMMYY { $$ = DDMMYY; }
        | MMDDYY { $$ = MMDDYY; }
        | MMYY { $$ = MMYY; }
        | MMDD { $$ = MMDD; }
        | MONTHDAY { $$ = MONTHDAY; }
        | DAYMONTH { $$ = DAYMONTH; }
        | DDMONTHSYYHHMMSS { $$ = DDMONTHSYYHHMMSS; }
        | MONTHS { $$ = MONTHS; }
        | MONTHL { $$ = MONTHL; }
        | DAYOFWEEKS { $$ = DAYOFWEEKS; }
        | DAYOFWEEKL { $$ = DAYOFWEEKL; }
        | DAYOFYEAR { $$ = DAYOFYEAR; }
        | HMS { $$ = HMS; }
        | MMDDHMS { $$ = MMDDHMS; }
        | MMDDYYHMS { $$ = MMDDYYHMS; }
        | DEGREESLON { $$ = DEGREESLON; }
        | DEGREESMMLON { $$ = DEGREESMMLON; }
        | DEGREESMMSSLON { $$ = DEGREESMMSSLON; }
        | MMSSLON { $$ = MMSSLON; }
        | DEGREESLAT { $$ = DEGREESLAT; }
        | DEGREESMMLAT { $$ = DEGREESMMLAT; }
        | DEGREESMMSSLAT { $$ = DEGREESMMSSLAT; }
        | MMSSLAT { $$ = MMSSLAT; }
        ;

horv: HORIZONTAL { $$ = $1; }
      | VERTICAL { $$ = $1; }
      ;

flowonoff:
        OFF { $$ = $1; }
        | ON { $$ = $1; }
        | NODES { $$ = $1; }
        | CENTER { $$ = $1; }
        ;

onoff:
        ON { $$ = ON; }
        | OFF { $$ = OFF; }
        ;

worldview: WORLD { $$ = WORLD; }
        | VIEW { $$ = VIEW; }
        ;

extremetype:
        XMIN { $$ = XMIN; }
        | XMAX { $$ = XMAX; }
        | YMIN { $$ = YMIN; }
        | YMAX { $$ = YMAX; }
        ;

torf:
        TRUEP { $$ = TRUEP; }
        | FALSEP { $$ = FALSEP; }
        ;

filltype: PATTERN { $$ = PATTERN; }
        | COLOR { $$ = COLOR; }
        | NONE { $$ = NONE; }
        ;

opchoice: ABOVE { $$ = ABOVE; }
        | BELOW { $$ = BELOW; }
        | LEFT { $$ = LEFT; }
        | RIGHT { $$ = RIGHT; }
        | TOP { $$ = TOP; }
        | BOTTOM { $$ = BOTTOM ; }
        | BOTH { $$ = BOTH ; }
        ;

units:
	MM { $$ = MM; }
	| CM { $$ = CM; }
	| M { $$ = M; }
	| KM { $$ = KM; }
	;

timeunits:
	SECONDS { $$ = $1; }
	| MINUTES  { $$ = $1; }
	| HOURS { $$ = $1; }
	| DAYS { $$ = $1; }
	| WEEKS { $$ = $1; }
	| MONTHS { $$ = $1; }
	| YEARS { $$ = $1; }
	;

asgn:
	VAR '[' expr ']' '=' expr
	{
	    int itmp = (int) $3 - 1;
	    if (itmp >= ls) {
		yyerror("subscript out of range");
		return 1;
	    } else {
		$1[itmp] = $6;
		result = $6;
	    }
	}
	;

vasgn:
	VAR '=' vexpr
	{
	    int i;
	    for (i = 0; i < lxy; i++) {
		$1[i] = $3[i];
	    }
	    result = $3[0];
	}
	| VAR '=' expr
	{
	    int i;
	    for (i = 0; i < lxy; i++) {
		$1[i] = $3;
	    }
	    result = $3;
	}
	;

vexpr:
	VAR
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i];
	    }
	}
	| expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1;
	    }
	}
	| expr '+' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 + $3;
	    }
	}
	| vexpr '+' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] + $3[i];
	    }
	}
	| expr '+' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 + $3[i];
	    }
	}
	| vexpr '+' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] + $3;
	    }
	}
	| expr '-' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 - $3;
	    }
	}
	| vexpr '-' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] - $3[i];
	    }
	}
	| expr '-' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 - $3[i];
	    }
	}
	| vexpr '-' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] - $3;
	    }
	}
	| expr '*' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 * $3;
	    }
	}
	| vexpr '*' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] * $3[i];
	    }
	}
	| expr '*' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 * $3[i];
	    }
	}
	| vexpr '*' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] * $3;
	    }
	}
	| expr '/' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    if ($3 == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 / $3;
	    }
	}
	| vexpr '/' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		if ($3[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] / $3[i];
	    }
	}
	| expr '/' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		if ($3[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1 / $3[i];
	    }
	}
	| vexpr '/' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    if ($3 == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] / $3;
	    }
	}
	| expr '^' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = pow($1, $3);
	    }
	}
	| expr '^' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = pow($1, $3[i]);
	    }
	}
	| vexpr '^' expr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = pow($1[i], $3);
	    }
	}
	| vexpr '^' vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = pow($1[i], $3[i]);
	    }
	}
	| ABS '(' expr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = fabs($3);
	    }
	}
	| ABS '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = fabs($3[i]);
	    }
	}
	| ACOS '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = acos($3[i]);
	    }
	}
	| ASIN '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = asin($3[i]);
	    }
	}
	| ATAN '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = atan($3[i]);
	    }
	}
	| ATAN2 '(' vexpr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = atan2($3[i], $5[i]);
	    }
	}
	| CEIL '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = ceil($3[i]);
	    }
	}
	| COS '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = cos($3[i]);
	    }
	}
	| DEG
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] *= M_PI / 180.0;
	    }
	}
	| DX
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = xx[i];
	    }
	}
	| DY
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = yy[i];
	    }
	}
	| ERF '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = erf($3[i]);
	    }
	}
	| ERFC '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = erfc($3[i]);
	    }
	}
	| EXP '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = exp($3[i]);
	    }
	}
	| FLOOR '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = floor($3[i]);
	    }
	}
	| HYPOT '(' vexpr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = hypot($3[i], $5[i]);
	    }
	}
	| HYPOT '(' expr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = hypot($3, $5[i]);
	    }
	}
	| HYPOT '(' vexpr ',' expr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = hypot($3[i], $5);
	    }
	}
	| HYPOT '(' expr ',' expr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = hypot($3, $5);
	    }
	}
	| INDEX
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = i + 1;
	    }
	}
	| SETNO
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1;
	    }
	}
	| INT '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = (int) $3[i];
	    }
	}
	| IRAND '(' NUMBER ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = lrand48() % (long) ($3);
	    }
	}
	| LGAMMA '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = lgamma($3[i]);
	    }
	}
	| LN '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = log($3[i]);
	    }
	}
	| LOG '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = log10($3[i]);
	    }
	}
	| LOGISTIC '(' vexpr ',' expr ',' expr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = 1.0 / (1.0 + exp(-($3[i] - $5)/ $7));
	    }
	}
	| MAXP '(' vexpr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $3[i] >= $5[i] ? $3[i] : $5[i];
	    }
	}
	| MINP '(' vexpr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $3[i] <= $5[i] ? $3[i] : $5[i];
	    }
	}
	| MOD '(' vexpr ',' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = fmod($3[i], $5[i]);
	    }
	}
	| NORM '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = fx($3[i]);
	    }
	}
	| NORMP '(' vexpr ')'
	{
	    int i;
	    double tmp;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = normp($3[i], &tmp);
	    }
	}
	| PI
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = M_PI;
	    }
	}
	| RAD
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = M_PI / 180.0;
	    }
	}
	| RAND
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = (double) drand48();
	    }
	}
	| SIN '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = sin($3[i]);
	    }
	}
	| SQR '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $3[i] * $3[i];
	    }
	}
	| SQRT '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = sqrt($3[i]);
	    }
	}
	| TAN '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = tan($3[i]);
	    }
	}
	| vexpr GT vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] > $3[i];
	    }
	}
	| vexpr LT vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] < $3[i];
	    }
	}
	| vexpr LE vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] <= $3[i];
	    }
	}
	| vexpr GE vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] >= $3[i];
	    }
	}
	| vexpr EQ vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] == $3[i];
	    }
	}
	| vexpr NE vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] != $3[i];
	    }
	}
	| vexpr AND vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] && $3[i];
	    }
	}
	| vexpr OR vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $1[i] || $3[i];
	    }
	}
	| NOT vexpr
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = !($2[i]);
	    }
	}
	| '(' vexpr ')'
	{
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = $2[i];
	    }
	}
	| '-' vexpr %prec UMINUS {
	    int i;
	    $$ = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = $$;
	    for (i = 0; i < lxy; i++) {
		$$[i] = -$2[i];
	    }
	}
	;

expr:	NUMBER
	|  FITPARM {
	    $$ = $1;
	}
	|  VAR '[' expr ']' {
	    $$ = $1[(int) $3];
	}
	| expr '+' expr {
	    $$ = $1 + $3;
	}
	| expr '-' expr {
	    $$ = $1 - $3;
	}
	| expr '*' expr {
	    $$ = $1 * $3;
	}
	| expr '/' expr
	{
	    if ($3 != 0.0) {
		$$ = $1 / $3;
	    } else {
		yyerror("Divide by Zero");
		return 1;
	    }
	}
	| expr '%' expr {
	    $$ = fmod($1, $3);
	}
	| expr '^' expr {
	    $$ = pow($1, $3);
	}
	| ABS '(' expr ')' {
	    $$ = fabs($3);
	}
	| ACOS '(' expr ')' {
	    $$ = acos($3);
	}
	| ASIN '(' expr ')' {
	    $$ = asin($3);
	}
	| ATAN '(' expr ')' {
	    $$ = atan($3);
	}
	| ATAN2 '(' expr ',' expr ')' {
	    $$ = atan2($3, $5);
	}
	| CEIL '(' expr ')' {
	    $$ = ceil($3);
	}
	| COS '(' expr ')' {
	    $$ = cos($3);
	}
	| DEG {
	    $$ = 180.0 / M_PI;
	}
	| DX {
	    $$ = *xx;
	}
	| DY {
	    $$ = *yy;
	}
	| ERF '(' expr ')' {
	    $$ = erf($3);
	}
	| ERFC '(' expr ')' {
	    $$ = erfc($3);
	}
	| EXP '(' expr ')' {
	    $$ = exp($3);
	}
	| FLOOR '(' expr ')' {
	    $$ = floor($3);
	}
	| HYPOT '(' expr ',' expr ')' {
	    $$ = hypot($3, $5);
	}
	| GRAPHNO '.' VX1 {
	    $$ = g[$1].v.xv1;
	}
	| GRAPHNO '.' VX2 {
	    $$ = g[$1].v.xv2;
	}
	| GRAPHNO '.' VY1 {
	    $$ = g[$1].v.yv1;
	}
	| GRAPHNO '.' VY2 {
	    $$ = g[$1].v.yv2;
	}
	| GRAPHNO '.' WX1 {
	    $$ = g[$1].w.xg1;
	}
	| GRAPHNO '.' WX2 {
	    $$ = g[$1].w.xg2;
	}
	| GRAPHNO '.' WY1 {
	    $$ = g[$1].w.yg1;
	}
	| GRAPHNO '.' WY2 {
	    $$ = g[$1].w.yg2;
	}
	| VX1 {
	    $$ = g[curg].v.xv1;
	}
	| VX2 {
	    $$ = g[curg].v.xv2;
	}
	| VY1 {
	    $$ = g[curg].v.yv1;
	}
	| VY2 {
	    $$ = g[curg].v.yv2;
	}
	| WX1 {
	    $$ = g[curg].w.xg1;
	}
	| WX2 {
	    $$ = g[curg].w.xg2;
	}
	| WY1 {
	    $$ = g[curg].w.yg1;
	}
	| WY2 {
	    $$ = g[curg].w.yg2;
	}
	| INDEX {
	    $$ = setindex;
	}
	| SETNO {
	    $$ = setsetno;
	}
	| INT '(' expr ')' {
	    $$ = (long) $3;
	}
	| IRAND '(' NUMBER ')' {
	    $$ = lrand48() % (long) ($3);
	}
	| LGAMMA '(' expr ')' {
	    $$ = lgamma($3);
	}
	| LN '(' expr ')' {
	    $$ = log($3);
	}
	| LOG '(' expr ')' {
	    $$ = log10($3);
	}
	| LOGISTIC '(' expr ',' expr ',' expr ')'
	{
	    $$ = 1.0 / (1.0 + exp(-($3 - $5)/ $7));
	}
	| MAXP '(' expr ',' expr ')' {
	    $$ = $3 >= $5 ? $3 : $5;
	}
	| MINP '(' expr ',' expr ')' {
	    $$ = $3 <= $5 ? $3 : $5;
	}
	| MOD '(' expr ',' expr ')' {
	    $$ = fmod($3, $5);
	}
	| NORM '(' expr ')' {
	    $$ = fx($3);
	}
	| NORMP '(' expr ')' {
	    double tmp;
	    $$ = normp($3, &tmp);
	}
	| PI {
	    $$ = M_PI;
	}
	| RAD {
	    $$ = M_PI / 180.0;
	}
	| RAND {
	    $$ = (double) drand48();
	}
	| SIN '(' expr ')' {
	    $$ = sin($3);
	}
	| SQR '(' expr ')' {
	    $$ = pow($3, 2.0);
	}
	| SQRT '(' expr ')' {
	    $$ = sqrt($3);
	}
	| TAN '(' expr ')' {
	    $$ = tan($3);
	}
	| IF '(' expr ')' expr {
	    if ((int) $3)
		$$ = $5;
	}
	| IF '(' expr ')' expr ELSE expr {
	    if ((int) $3) {
		$$ = $5;
	    } else {
		$$ = $7;
	    }
	}
	| expr GT expr {
	    $$ = $1 > $3;
	}
	| expr LT expr {
	    $$ = $1 < $3;
	}
	| expr LE expr {
	    $$ = $1 <= $3;
	}
	| expr GE expr {
	    $$ = $1 >= $3;
	}
	| expr EQ expr {
	    $$ = $1 == $3;
	}
	| expr NE expr {
	    $$ = $1 != $3;
	}
	| expr AND expr {
	    $$ = $1 && $3;
	}
	| expr OR expr {
	    $$ = $1 || $3;
	}
	| NOT expr {
	    $$ = !($2);
	}
	| '(' expr ')' {
	    $$ = $2;
	}
	| '-' expr %prec UMINUS {
	    $$ = -$2;
	}
	;
%%

void fixupstr(char *val)
{
    int vl = strlen(val);
    lowtoupper(val);
    val[vl + 1] = 0;
    val[vl] = '\n';
}

void scanner(char *s, double *x, double *y, int len, double *a, double *b, double *c, double *d, int lenscr, int i, int setno, int *errpos)
{
    interr = 0;
    whichgraph = cg;
    whichset = setno;
    if (s[0] == '#') {
	return;
    }
    pos = 0;
    aa = a;
    bb = b;
    cc = c;
    dd = d;
    xx = x;
    yy = y;
    lxy = len;
    ls = lenscr;
    setindex = i + 1;
    curset = setsetno = setno;
    strcpy(f_string, s);
    strcpy(statusstr, s);
    fcnt = 0;
    yyparse();
    *errpos = interr;
    for (i = 0; i < fcnt; i++) {
	free(freelist[i]);
	freelist[i] = NULL;
    }
}

void runbatch(char *bfile)
{
    double x, y, a, b, c, d;
    int i, setno, errpos, lcnt = 1;
    char stext[256];
    FILE *fp;
    if (strcmp("stdin", bfile)) {
	fp = fopen(bfile, "r");
    }
    else {
	fp = stdin;
    }
    if (fp == NULL) {
        fprintf(stderr, "Error opening batch file \"%s\"\n", bfile);
        exit(1);
    }
    while(fgets(stext, 255, fp) != NULL) {
        if (stext[0] == '#') {
            continue;
        }
        lowtoupper(stext);
/* TODO check on 0, 0 here for index and setno */
        scanner(stext, &x, &y, 1, ax, bx, cx, dx, 1, 0, 0, &errpos);
        stext[0] = 0;
        if (gotparams && paramfile[0]) {
            if (!getparms(cg, paramfile)) {
            }
            gotparams = 0;
        } else if (gotread && readfile[0]) {
            if (getdata(cg, readfile, readsrc, readtype)) {
                drawgraph();
            }
            gotread = 0;
        }
    }
    if (fp != stdin) {
	fclose(fp);
    }
}

symtab_entry key[] = {
"ABORT", ABORT,
"ABOVE", ABOVE,
"ABSOLUTE", ABSOLUTE,
"ACTIVATE", ACTIVATE,
"ACTIVE", ACTIVE,
"ADCIRC", ADCIRC,
"ADCIRC3DFLOW", ADCIRC3DFLOW,
"ALL", ALL,
"ALT", ALT,
"ALTERNATE", ALTERNATE,
"ALTXAXIS", ALTXAXIS,
"ALTYAXIS", ALTYAXIS,
"AMP", AMP,
"ANGLE", ANGLE,
"ANNOTATE", ANNOTATE,
"APPEND", APPEND,
"AREA", AREA,
"ARROW", ARROW,
"ASCEND", ASCEND,
"ASCENDING", ASCEND,
"AT", AT,
"ATTACH", ATTACH,
"AUTO", AUTO,
"AUTOSCALE", AUTOSCALE,
"AUTOTICKS", AUTOTICKS,
"AVERAGE", AVERAGE,
"AVG", AVG,
"AXES", AXES,
"AXIS", AXIS,
"BACKBUFFER", BACKBUFFER,
"BACKGROUND", BACKGROUND,
"BAR", BAR,
"BATCH", BATCH,
"BATH", BATH,
"BATHYMETRY", BATHYMETRY,
"BELOW", BELOW,
"BIN", BIN,
"BINARY", BINARY,
"BOTH", BOTH,
"BOTTOM", BOTTOM,
"BOUNDARY", BOUNDARY,
"BOX", BOX,
"CELLS", CELLS,
"CENTER", CENTER,
"CH3D", CH3D,
"CHAR", CHAR,
"CHDIR", CHDIR,
"CIRCLE", CIRCLE,
"CLEAR", CLEAR,
"CLICK", CLICK,
"CLOCK", CLOCK,
"CLOSE", CLOSE,
"CM", CM,
"CMAP", CMAP,
"COLOR", COLOR,
"COLORMAP", COLORMAP,
"COMMENT", COMMENT,
"CONC", CONC,
"CONCENTRATION", CONCENTRATION,
"CONCENTRATIONS", CONCENTRATIONS,
"COPY", COPY,
"COURANT", COURANT,
"CROSS", CROSS,
"CYCLE", CYCLE,
"DAYMONTH", DAYMONTH,
"DAYOFWEEKL", DAYOFWEEKL,
"DAYOFWEEKS", DAYOFWEEKS,
"DAYOFYEAR", DAYOFYEAR,
"DDMMYY", DDMMYY,
"DECIMAL", DECIMAL,
"DEF", DEF,
"DEFAULT", DEFAULT,
"DEGREESLAT", DEGREESLAT,
"DEGREESLON", DEGREESLON,
"DEGREESMMLAT", DEGREESMMLAT,
"DEGREESMMLON", DEGREESMMLON,
"DEGREESMMSSLAT", DEGREESMMSSLAT,
"DEGREESMMSSLON", DEGREESMMSSLON,
"DELAYP", DELAYP,
"DELETE", DELETE,
"DEPTH", DEPTH,
"DEPTHS", DEPTHS,
"DESCEND", DESCEND,
"DESCENDING", DESCEND,
"DEVICE", DEVICE,
"DEVXY", DEVXY,
"DFT", DFT,
"DIAMOND", DIAMOND,
"DIFFERENCE", DIFFERENCE,
"DISK", DISK,
"DISPLAY", DISPLAY,
"DOT", DOT,
"DOUBLEBUFFER", DOUBLEBUFFER,
"DOWN", DOWN,
"DRAW2", DRAW2,
"DROGUE", DROGUE,
"DROGUES", DROGUES,
"DT", DT,
"DXDX", DXDX,
"DXP", DXP,
"DYDY", DYDY,
"DYP", DYP,
"ECHO", ECHO,
"EDIT", EDIT,
"ELA", ELA,
"ELCIRC", ELCIRC,
"ELEMENT", ELEMENT,
"ELEMENTS", ELEMENTS,
"ELEV", ELEV,
"ELEVATION", ELEVATION,
"ELEVATIONS", ELEVATIONS,
"ELEVMARKER", ELEVMARKER,
"ELSE", ELSE,
"END", END,
"ERRORBAR", ERRORBAR,
"EXIT", EXIT,
"EXPAND", EXPAND,
"EXPONENTIAL", EXPONENTIAL,
"FACTOR", FACTOR,
"FALSE", FALSEP,
"FAST", FAST,
"FEET", FEET,
"FFT", FFT,
"FILE", FILEP,
"FILL", FILL,
"FIND", FIND,
"FIXEDPOINT", FIXEDPOINT,
"FLOW", FLOW,
"FLUSH", FLUSH,
"FLUX", FLUX,
"FOCUS", FOCUS,
"FOLLOWS", FOLLOWS,
"FONT", FONTP,
"FOREGROUND", FOREGROUND,
"FORMAT", FORMAT,
"FORT14", FORT14,
"FORT63", FORT63,
"FORT64", FORT64,
"FORWARD", FORWARD,
"FRAME", FRAMEP,
"FREQ", FREQ,
"FRONTBUFFER", FRONTBUFFER,
"GENERAL", GENERAL,
"GETP", GETP,
"GOTO", GOTO,
"GRAPH", GRAPH,
"GRAPHNO", GRAPHNO,
"GRAPHS", GRAPHS,
"GRAPHTYPE", GRAPHTYPE,
"GRID", GRID,
"HARDCOPY", HARDCOPY,
"HBAR", HBAR,
"HELP", HELP,
"HGAP", HGAP,
"HIDDEN", HIDDEN,
"HISTBOX", HISTBOX,
"HISTO", HISTO,
"HISTORY", HISTORY,
"HMS", HMS,
"HORIZONTAL", HORIZONTAL,
"HPGLL", HPGLL,
"HPGLP", HPGLP,
"IF", IF,
"IGNORE", IGNORE,
"IHL", IHL,
"IMAGE", IMAGE,
"IMAGES", IMAGES,
"IN", IN,
"INCLUDE", INCLUDE,
"INFO", INFO,
"INIT", INIT,
"INITGRAPHICS", INITGRAPHICS,
"INOUT", INOUT,
"INTEGRATE", INTEGRATE,
"INTERP", INTERP,
"INUNDATION", INUNDATION,
"INVDFT", INVDFT,
"INVFFT", INVFFT,
"ISOLINE", ISOLINE,
"ISOLINES", ISOLINES,
"JUST", JUST,
"KILL", KILL,
"KM", KM,
"LABEL", LABEL,
"LAYOUT", LAYOUT,
"LEAVE", LEAVE,
"LEAVEGRAPHICS", LEAVEGRAPHICS,
"LEFT", LEFT,
"LEGEND", LEGEND,
"LENGTH", LENGTH,
"LEVEL", LEVEL,
"LEVELS", LEVELS,
"LIMITS", LIMITS,
"LINE", LINE,
"LINES", LINES,
"LINESTYLE", LINESTYLE,
"LINETO", LINETO,
"LINEW", LINEW,
"LINEWIDTH", LINEWIDTH,
"LINK", LINK,
"LOAD", LOAD,
"LOC", LOC,
"LOCATE", LOCATE,
"LOCATOR", LOCATOR,
"LOCTYPE", LOCTYPE,
"LOG", LOG,
"LOGX", LOGX,
"LOGXY", LOGXY,
"LOGY", LOGY,
"M", M,
"MAG", MAG,
"MAGNITUDE", MAGNITUDE,
"MAJOR", MAJOR,
"MAPSCALE", MAPSCALE,
"MARKER", MARKER,
"MARKERS", MARKERS,
"MAX", MAXP,
"MAXLEVELS", MAXLEVELS,
"METHOD", METHOD,
"MIFL", MIFL,
"MIFP", MIFP,
"MILES", MILES,
"MIN", MINP,
"MINOR", MINOR,
"MISSINGP", MISSINGP,
"MM", MM,
"MMDD", MMDD,
"MMDDHMS", MMDDHMS,
"MMDDYY", MMDDYY,
"MMDDYYHMS", MMDDYYHMS,
"MMSSLAT", MMSSLAT,
"MMSSLON", MMSSLON,
"MMYY", MMYY,
"MONTHDAY", MONTHDAY,
"MONTHL", MONTHL,
"MONTHS", MONTHS,
"MOVE", MOVE,
"MOVE2", MOVE2,
"MOVETO", MOVETO,
"NEGATE", NEGATE,
"NO", NO,
"NODE", NODE,
"NODES", NODES,
"NONE", NONE,
"NORMAL", NORMAL,
"NORTH", NORTH,
"NXY", NXY,
"OFF", OFF,
"OFFSETX", OFFSETX,
"OFFSETY", OFFSETY,
"ON", ON,
"OP", OP,
"OPEN", OPEN,
"ORIENT", ORIENT,
"OUT", OUT,
"PAGE", PAGE,
"PARA", PARA,
"PARALLEL", PARALLEL,
"PARAMETERS", PARAMETERS,
"PARAMS", PARAMS,
"PARMS", PARMS,
"PATTERN", PATTERN,
"PER", PER,
"PERIMETER", PERIMETER,
"PERP", PERP,
"PERPENDICULAR", PERPENDICULAR,
"PHASE", PHASE,
"PIE", PIE,
"PIPE", PIPE,
"PLACE", PLACE,
"PLUS", PLUS,
"POINT", POINT,
"POLAR", POLAR,
"POLY", POLY,
"POLYI", POLYI,
"POLYO", POLYO,
"POP", POP,
"POWER", POWER,
"PREC", PREC,
"PREFIX", PREFIX,
"PREPEND", PREPEND,
"PRINT", PRINT,
"PROP", PROP,
"PS", PS,
"PSCOLORL", PSCOLORL,
"PSCOLORP", PSCOLORP,
"PSMONOL", PSMONOL,
"PSMONOP", PSMONOP,
"PUSH", PUSH,
"PUTP", PUTP,
"QUIT", QUIT,
"READ", READ,
"READBIN", READBIN,
"REDRAW", REDRAW,
"REGION", REGION,
"REGIONS", REGIONS,
"REGNUM", REGNUM,
"REGRESS", REGRESS,
"REMOVE", REMOVE,
"RENDER", RENDER,
"REPORT", REPORT,
"RESET", RESET,
"REVERSE", REVERSE,
"REWIND", REWIND,
"RIGHT", RIGHT,
"RISER", RISER,
"ROT", ROT,
"RUN", RUN,
"SALINITY", SALINITY,
"SAMPLE", SAMPLE,
"SAVE", SAVE,
"SCALAR", SCALAR,
"SCALE", SCALE,
"SCIENTIFIC", SCIENTIFIC,
"SET", SET,
"SETS", SETS,
"SHOW", SHOW,
"SHRINK", SHRINK,
"SIGMA", SIGMA,
"SIGN", SIGN,
"SIZE", SIZE,
"SKIP", SKIP,
"SLEEP", SLEEP,
"SLICE", SLICE,
"SOURCE", SOURCE,
"SPEC", SPEC,
"SPECIFIED", SPECIFIED,
"SPECTRUM", SPECTRUM,
"SPLITS", SPLITS,
"SQUARE", SQUARE,
"STACK", STACK,
"STACKEDBAR", STACKEDBAR,
"STACKEDHBAR", STACKEDHBAR,
"STACKEDLINE", STACKEDLINE,
"STAGGER", STAGGER,
"STAR", STAR,
"START", START,
"STARTSTEP", STARTSTEP,
"STARTTYPE", STARTTYPE,
"STATUS", STATUS,
"STEP", STEP,
"STOP", STOP,
"STREAMLINES", STREAMLINES,
"STRING", STRING,
"STRINGS", STRINGS,
"SUBTITLE", SUBTITLE,
"SURFACE", SURFACE,
"SWAPBUFFER", SWAPBUFFER,
"SYMBOL", SYMBOL,
"SYSTEM", SYSTEM,
"TEANL", TEANL,
"TEXT", TEXT,
"TICK", TICKP,
"TICKLABEL", TICKLABEL,
"TICKMARKS", TICKMARKS,
"TIDALCLOCK", TIDALCLOCK,
"TIME", TIME,
"TIMEINFO", TIMEINFO,
"TIMELINE", TIMELINE,
"TITLE", TITLE,
"TO", TO,
"TOP", TOP,
"TOTAL", TOTAL,
"TRACK", TRACK,
"TRANSECT", TRANSECT,
"TRIANGLE1", TRIANGLE1,
"TRIANGLE2", TRIANGLE2,
"TRIANGLE3", TRIANGLE3,
"TRIANGLE4", TRIANGLE4,
"TRUE", TRUEP,
"TYPE", TYPE,
"UNITS", UNITS,
"UP", UP,
"VALUE", VALUE,
"VECTOR", VECTOR,
"VEL", VEL,
"VELOCITY", VELOCITY,
"VERTICAL", VERTICAL,
"VGAP", VGAP,
"VIEW", VIEW,
"VSCALE", VSCALE,
"VX1", VX1,
"VX2", VX2,
"VY1", VY1,
"VY2", VY2,
"WIDTH", WIDTH,
"WIND", WIND,
"WITH", WITH,
"WORLD", WORLD,
"WRAP", WRAP,
"WRITE", WRITE,
"WSCALE", WSCALE,
"WX1", WX1,
"WX2", WX2,
"WY1", WY1,
"WY2", WY2,
"X1", X1,
"X2", X2,
"X3", X3,
"X4", X4,
"X5", X5,
"XAXES", XAXES,
"XAXIS", XAXIS,
"XCOR", XCOR,
"XMAX", XMAX,
"XMIN", XMIN,
"XY", XY,
"XYARC", XYARC,
"XYBOX", XYBOX,
"XYDX", XYDX,
"XYDXDX", XYDXDX,
"XYDXDY", XYDXDY,
"XYDY", XYDY,
"XYDYDY", XYDYDY,
"XYFIXED", XYFIXED,
"XYHILO", XYHILO,
"XYRT", XYRT,
"XYSEG", XYSEG,
"XYSTRING", XYSTRING,
"XYUV", XYUV,
"XYX2Y2", XYX2Y2,
"XYXX", XYXX,
"XYYY", XYYY,
"XYZ", XYZ,
"XYZW", XYZW,
"Y1", Y1,
"Y2", Y2,
"Y3", Y3,
"Y4", Y4,
"Y5", Y5,
"YAXES", YAXES,
"YAXIS", YAXIS,
"YES", YES,
"YMAX", YMAX,
"YMIN", YMIN,
"ZEROXAXIS", ZEROXAXIS,
"ZEROYAXIS", ZEROYAXIS,
"ZOOM", ZOOM,
"ZOOMBOX", ZOOMBOX
};

int maxparms = sizeof(key) / sizeof(symtab_entry);
int maxfunc = sizeof(key) / sizeof(symtab_entry);

int findf(symtab_entry *key, char *s, int tlen)
{

    int low, high, mid;

    low = 0;
    high = tlen - 1;
    while (low <= high) {
	mid = (low + high) / 2;
	if (strcmp(s, key[mid].s) < 0) {
	    high = mid - 1;
	} else {
	    if (strcmp(s, key[mid].s) > 0) {
		low = mid + 1;
	    } else {
		return (mid);
	    }
	}
    }
    return (-1);
}

int getcharstr(void)
{
    if (pos >= strlen(f_string))
	 return EOF;
    return (f_string[pos++]);
}

void ungetchstr(void)
{
    if (pos > 0)
	pos--;
}

int yylex(void)
{
    int c, i;
    int found;
    static char s[256];
    char sbuf[256];

    while ((c = getcharstr()) == ' ' || c == '\t');
    if (c == EOF) {
	return (0);
    }
    if (c == '"') {
	i = 0;
	while ((c = getcharstr()) != '"' && c != EOF) {
	    if (c == '\\') {
		int ctmp;
		ctmp = getcharstr();
		if (ctmp != '"') {
		    ungetchstr();
		}
		else {
		    c = ctmp;
		}
	    }
	    s[i] = c;
	    i++;
	}
	if (c == EOF) {
	    sprintf(sbuf, "Nonterminating string\n");
	    yyerror(sbuf);
	    return 0;
	}
	s[i] = '\0';
	yylval.str = s;
	return CHRSTR;
    }
    if (c == '.' || isdigit(c)) {
	char stmp[80];
	double d;
	int i, gotdot = 0;

	i = 0;
	while (c == '.' || isdigit(c)) {
	    if (c == '.') {
		if (gotdot) {
		    yyerror("Reading number, too many dots");
	    	    return 0;
		} else {
		    gotdot = 1;
		}
	    }
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (c == 'E' || c == 'e') {
	    stmp[i++] = c;
	    c = getcharstr();
	    if (c == '+' || c == '-') {
		stmp[i++] = c;
		c = getcharstr();
	    }
	    while (isdigit(c)) {
		stmp[i++] = c;
		c = getcharstr();
	    }
	}
	if (gotdot && i == 1) {
	    ungetchstr();
	    return '.';
	}
	stmp[i] = '\0';
	ungetchstr();
	sscanf(stmp, "%lf", &d);
	yylval.val = d;
	return NUMBER;
    }
/* graphs, sets, regions resp. */
    if (c == 'G' || c == 'S' || c == 'R') {
	char stmp[80];
	double d;
	int i = 0, ctmp = c, gn, sn, rn;
	c = getcharstr();
	while (isdigit(c)) {
	    stmp[i++] = c;
	    c = getcharstr();
	}
	if (i == 0) {
	    c = ctmp;
	    ungetchstr();
	} else {
	    ungetchstr();
	    if (ctmp == 'G') {
	        stmp[i] = '\0';
		gn = atoi(stmp);
		if (gn >= 0 && gn < maxgraph) {
		    yylval.ival = gn;
		    whichgraph = gn;
		    return GRAPHNO;
		}
	    }
	}
    }
    if (isalpha(c)) {
	char *p = sbuf;
	int gno = -1, setno = -1, xy = -1, elno = -1;

	do {
	    *p++ = c;
	} while ((c = getcharstr()) != EOF && isalnum(c));
	ungetchstr();
	*p = '\0';
        if (debuglevel == 2) {
	    printf("->%s<-\n", sbuf);
	}
	if ((found = findf(key, sbuf, maxfunc)) >= 0) {
	    if (key[found].type == VAR) {
		switch (sbuf[0]) {
		case 'A':
		    yylval.ptr = aa;
		    return VAR;
		case 'B':
		    yylval.ptr = bb;
		    return VAR;
		case 'C':
		    yylval.ptr = cc;
		    return VAR;
		case 'D':
		    yylval.ptr = dd;
		    return VAR;
		}
	    }
	    else if (key[found].type == FITPARM) {
		int index = sbuf[1] - '0';
		yylval.val = nonl_parms[index];
		return FITPARM;
	    }
	    else { /* set up special cases */
		switch (key[found].type) {
		case XAXIS:
		    naxis = 0;
		    break;
		case YAXIS:
		    naxis = 1;
		    break;
		case ZEROXAXIS:
		    naxis = 2;
		    break;
		case ZEROYAXIS:
		    naxis = 3;
		    break;
		case ALTXAXIS:
		    naxis = 4;
		    break;
		case ALTYAXIS:
		    naxis = 5;
		    break;
		case AXES:
		    naxis = 6;
		    break;
		case XAXES:
		    naxis = 7;
		    break;
		case YAXES:
		    naxis = 8;
		    break;
		case GRAPHS:
		    yylval.ival = -1;
		    whichgraph = -1;
		    return GRAPHS;
		case SETS:
		    yylval.ival = -1;
		    whichset = -1;
		    return SETS;
		default:
		    break;
		}
	    }
	    yylval.func = key[found].type;
	    return key[found].type;
	} else {
	    strcat(sbuf, ": No such function or variable");
	    yyerror(sbuf);
	    return 0;
	}
    }
    switch (c) {
    case '>':
	return follow('=', GE, GT);
    case '<':
	return follow('=', LE, LT);
    case '=':
	return follow('=', EQ, '=');
    case '!':
	return follow('=', NE, NOT);
    case '|':
	return follow('|', OR, '|');
    case '&':
	return follow('&', AND, '&');
    case '\n':
	return '\n';
    default:
	return c;
    }
}

int follow(int expect, int ifyes, int ifno)
{
    int c = getcharstr();

    if (c == expect) {
	return ifyes;
    }
    ungetchstr();
    return ifno;
}

void yyerror(char *s)
{
    int i;
    char buf[256];
    sprintf(buf, "Error: %s: %s", s, f_string);
    i = strlen(buf);
    buf[i - 1] = 0;
    errwin(buf);
    interr = 1;
}

#define C1 0.1978977093962766
#define C2 0.1352915131768107

double rnorm(double mean, double sdev)
{
    double u = drand48();

    return mean + sdev * (pow(u, C2) - pow(1.0 - u, C2)) / C1;
}

double fx(double x)
{
    return 1.0 / sqrt(2.0 * M_PI) * exp(-x * x * 0.5);
}

double normp(double b, double *s)
{
    double sum, dx, a = -8.0, fx(double x);
    int i, n = 48;

    sum = fx(a) + fx(b);
    dx = (b - a) / n;
    for (i = 1; i <= ((n - 1) / 2); i++)
	sum = sum + 4.0 * fx(a + (2.0 * i - 1.0) * dx) + 2.0 * fx(a + 2.0 * i * dx);
    sum = sum + 4.0 * fx(b - dx);
    *s = fx(b);
    return sum * dx / 3.0;
}
