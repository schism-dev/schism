/* A Bison parser, made by GNU Bison 3.7.4.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2020 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30704

/* Bison version string.  */
#define YYBISON_VERSION "3.7.4"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 11 "gram.y"

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


#line 173 "y.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    VAR = 258,                     /* VAR  */
    X = 259,                       /* X  */
    Y = 260,                       /* Y  */
    CHRSTR = 261,                  /* CHRSTR  */
    FITPARM = 262,                 /* FITPARM  */
    NUMBER = 263,                  /* NUMBER  */
    ABS = 264,                     /* ABS  */
    ACOS = 265,                    /* ACOS  */
    ASIN = 266,                    /* ASIN  */
    ATAN = 267,                    /* ATAN  */
    ATAN2 = 268,                   /* ATAN2  */
    CEIL = 269,                    /* CEIL  */
    COS = 270,                     /* COS  */
    DEG = 271,                     /* DEG  */
    DX = 272,                      /* DX  */
    DY = 273,                      /* DY  */
    ERF = 274,                     /* ERF  */
    ERFC = 275,                    /* ERFC  */
    EXP = 276,                     /* EXP  */
    FLOOR = 277,                   /* FLOOR  */
    HYPOT = 278,                   /* HYPOT  */
    INDEX = 279,                   /* INDEX  */
    INT = 280,                     /* INT  */
    IRAND = 281,                   /* IRAND  */
    LGAMMA = 282,                  /* LGAMMA  */
    LN = 283,                      /* LN  */
    LOG = 284,                     /* LOG  */
    LOGISTIC = 285,                /* LOGISTIC  */
    MAXP = 286,                    /* MAXP  */
    MINP = 287,                    /* MINP  */
    MINMAX = 288,                  /* MINMAX  */
    MOD = 289,                     /* MOD  */
    NORM = 290,                    /* NORM  */
    NORMP = 291,                   /* NORMP  */
    PI = 292,                      /* PI  */
    RAD = 293,                     /* RAD  */
    RAND = 294,                    /* RAND  */
    SETNO = 295,                   /* SETNO  */
    SIN = 296,                     /* SIN  */
    SQR = 297,                     /* SQR  */
    SQRT = 298,                    /* SQRT  */
    TAN = 299,                     /* TAN  */
    INUM = 300,                    /* INUM  */
    ABORT = 301,                   /* ABORT  */
    ABOVE = 302,                   /* ABOVE  */
    ABSOLUTE = 303,                /* ABSOLUTE  */
    ACTIVATE = 304,                /* ACTIVATE  */
    ACTIVE = 305,                  /* ACTIVE  */
    ADCIRC = 306,                  /* ADCIRC  */
    ADCIRC3DFLOW = 307,            /* ADCIRC3DFLOW  */
    ALL = 308,                     /* ALL  */
    ALT = 309,                     /* ALT  */
    ALTERNATE = 310,               /* ALTERNATE  */
    ALTXAXIS = 311,                /* ALTXAXIS  */
    ALTYAXIS = 312,                /* ALTYAXIS  */
    AMP = 313,                     /* AMP  */
    ANGLE = 314,                   /* ANGLE  */
    ANNOTATE = 315,                /* ANNOTATE  */
    APPEND = 316,                  /* APPEND  */
    AREA = 317,                    /* AREA  */
    ARROW = 318,                   /* ARROW  */
    ASCEND = 319,                  /* ASCEND  */
    AT = 320,                      /* AT  */
    ATTACH = 321,                  /* ATTACH  */
    AUTO = 322,                    /* AUTO  */
    AUTOSCALE = 323,               /* AUTOSCALE  */
    AUTOTICKS = 324,               /* AUTOTICKS  */
    AVERAGE = 325,                 /* AVERAGE  */
    AVG = 326,                     /* AVG  */
    AXES = 327,                    /* AXES  */
    AXIS = 328,                    /* AXIS  */
    BACKBUFFER = 329,              /* BACKBUFFER  */
    BACKGROUND = 330,              /* BACKGROUND  */
    BAR = 331,                     /* BAR  */
    BATCH = 332,                   /* BATCH  */
    BATH = 333,                    /* BATH  */
    BATHYMETRY = 334,              /* BATHYMETRY  */
    COURANT = 335,                 /* COURANT  */
    BELOW = 336,                   /* BELOW  */
    BIN = 337,                     /* BIN  */
    BINARY = 338,                  /* BINARY  */
    BOTH = 339,                    /* BOTH  */
    BOTTOM = 340,                  /* BOTTOM  */
    BOUNDARY = 341,                /* BOUNDARY  */
    BOX = 342,                     /* BOX  */
    CELLS = 343,                   /* CELLS  */
    CENTER = 344,                  /* CENTER  */
    CH3D = 345,                    /* CH3D  */
    CHAR = 346,                    /* CHAR  */
    CHDIR = 347,                   /* CHDIR  */
    CIRCLE = 348,                  /* CIRCLE  */
    CLEAR = 349,                   /* CLEAR  */
    CLICK = 350,                   /* CLICK  */
    CLOCK = 351,                   /* CLOCK  */
    CLOSE = 352,                   /* CLOSE  */
    CM = 353,                      /* CM  */
    CMAP = 354,                    /* CMAP  */
    COLOR = 355,                   /* COLOR  */
    COLORMAP = 356,                /* COLORMAP  */
    COMMENT = 357,                 /* COMMENT  */
    CONC = 358,                    /* CONC  */
    CONCENTRATION = 359,           /* CONCENTRATION  */
    CONCENTRATIONS = 360,          /* CONCENTRATIONS  */
    COPY = 361,                    /* COPY  */
    CROSS = 362,                   /* CROSS  */
    CYCLE = 363,                   /* CYCLE  */
    DAYMONTH = 364,                /* DAYMONTH  */
    DAYOFWEEKL = 365,              /* DAYOFWEEKL  */
    DAYOFWEEKS = 366,              /* DAYOFWEEKS  */
    DAYOFYEAR = 367,               /* DAYOFYEAR  */
    DAYS = 368,                    /* DAYS  */
    DDMMYY = 369,                  /* DDMMYY  */
    DDMONTHSYYHHMMSS = 370,        /* DDMONTHSYYHHMMSS  */
    DECIMAL = 371,                 /* DECIMAL  */
    DEF = 372,                     /* DEF  */
    DEFAULT = 373,                 /* DEFAULT  */
    DEGREESLAT = 374,              /* DEGREESLAT  */
    DEGREESLON = 375,              /* DEGREESLON  */
    DEGREESMMLAT = 376,            /* DEGREESMMLAT  */
    DEGREESMMLON = 377,            /* DEGREESMMLON  */
    DEGREESMMSSLAT = 378,          /* DEGREESMMSSLAT  */
    DEGREESMMSSLON = 379,          /* DEGREESMMSSLON  */
    DELAYP = 380,                  /* DELAYP  */
    DELETE = 381,                  /* DELETE  */
    DEPTH = 382,                   /* DEPTH  */
    DEPTHS = 383,                  /* DEPTHS  */
    DESCEND = 384,                 /* DESCEND  */
    DEVICE = 385,                  /* DEVICE  */
    DEVXY = 386,                   /* DEVXY  */
    DFT = 387,                     /* DFT  */
    DT = 388,                      /* DT  */
    DIAMOND = 389,                 /* DIAMOND  */
    DIFFERENCE = 390,              /* DIFFERENCE  */
    DISK = 391,                    /* DISK  */
    DISPLAY = 392,                 /* DISPLAY  */
    DOT = 393,                     /* DOT  */
    DOUBLEBUFFER = 394,            /* DOUBLEBUFFER  */
    DOWN = 395,                    /* DOWN  */
    DRAW2 = 396,                   /* DRAW2  */
    DROGUE = 397,                  /* DROGUE  */
    DROGUES = 398,                 /* DROGUES  */
    DRY = 399,                     /* DRY  */
    DXDX = 400,                    /* DXDX  */
    DXP = 401,                     /* DXP  */
    DYDY = 402,                    /* DYDY  */
    DYP = 403,                     /* DYP  */
    ECHO = 404,                    /* ECHO  */
    EDIT = 405,                    /* EDIT  */
    ELA = 406,                     /* ELA  */
    ELCIRC = 407,                  /* ELCIRC  */
    ELEMENT = 408,                 /* ELEMENT  */
    ELEMENTS = 409,                /* ELEMENTS  */
    ELEV = 410,                    /* ELEV  */
    ELEVATION = 411,               /* ELEVATION  */
    ELEVATIONS = 412,              /* ELEVATIONS  */
    ELEVMARKER = 413,              /* ELEVMARKER  */
    ELLIPSE = 414,                 /* ELLIPSE  */
    ELLIPSES = 415,                /* ELLIPSES  */
    ELLIPSEZ = 416,                /* ELLIPSEZ  */
    ELSE = 417,                    /* ELSE  */
    END = 418,                     /* END  */
    ERRORBAR = 419,                /* ERRORBAR  */
    EXIT = 420,                    /* EXIT  */
    EXPAND = 421,                  /* EXPAND  */
    EXPONENTIAL = 422,             /* EXPONENTIAL  */
    FACTOR = 423,                  /* FACTOR  */
    FALSEP = 424,                  /* FALSEP  */
    FAST = 425,                    /* FAST  */
    FEET = 426,                    /* FEET  */
    FFT = 427,                     /* FFT  */
    FILEP = 428,                   /* FILEP  */
    FILL = 429,                    /* FILL  */
    FIND = 430,                    /* FIND  */
    FIXEDPOINT = 431,              /* FIXEDPOINT  */
    FLOW = 432,                    /* FLOW  */
    FLUSH = 433,                   /* FLUSH  */
    FLUX = 434,                    /* FLUX  */
    FOCUS = 435,                   /* FOCUS  */
    FOLLOWS = 436,                 /* FOLLOWS  */
    FONTP = 437,                   /* FONTP  */
    FOREGROUND = 438,              /* FOREGROUND  */
    FORMAT = 439,                  /* FORMAT  */
    FORT14 = 440,                  /* FORT14  */
    FORT63 = 441,                  /* FORT63  */
    FORT64 = 442,                  /* FORT64  */
    FORWARD = 443,                 /* FORWARD  */
    FRAMEP = 444,                  /* FRAMEP  */
    FREQ = 445,                    /* FREQ  */
    FRONTBUFFER = 446,             /* FRONTBUFFER  */
    GENERAL = 447,                 /* GENERAL  */
    GETP = 448,                    /* GETP  */
    GOTO = 449,                    /* GOTO  */
    GRAPH = 450,                   /* GRAPH  */
    GRAPHNO = 451,                 /* GRAPHNO  */
    GRAPHS = 452,                  /* GRAPHS  */
    GRAPHTYPE = 453,               /* GRAPHTYPE  */
    GRID = 454,                    /* GRID  */
    HARDCOPY = 455,                /* HARDCOPY  */
    HBAR = 456,                    /* HBAR  */
    HELP = 457,                    /* HELP  */
    HGAP = 458,                    /* HGAP  */
    HIDDEN = 459,                  /* HIDDEN  */
    HISTBOX = 460,                 /* HISTBOX  */
    HISTO = 461,                   /* HISTO  */
    HISTORY = 462,                 /* HISTORY  */
    HMS = 463,                     /* HMS  */
    HORIZONTAL = 464,              /* HORIZONTAL  */
    HOURS = 465,                   /* HOURS  */
    HPGLL = 466,                   /* HPGLL  */
    HPGLP = 467,                   /* HPGLP  */
    IF = 468,                      /* IF  */
    IGNORE = 469,                  /* IGNORE  */
    IHL = 470,                     /* IHL  */
    IMAGE = 471,                   /* IMAGE  */
    IMAGES = 472,                  /* IMAGES  */
    IN = 473,                      /* IN  */
    INCLUDE = 474,                 /* INCLUDE  */
    INFO = 475,                    /* INFO  */
    INIT = 476,                    /* INIT  */
    INITGRAPHICS = 477,            /* INITGRAPHICS  */
    INOUT = 478,                   /* INOUT  */
    INTEGRATE = 479,               /* INTEGRATE  */
    INTERP = 480,                  /* INTERP  */
    INUNDATION = 481,              /* INUNDATION  */
    INVDFT = 482,                  /* INVDFT  */
    INVFFT = 483,                  /* INVFFT  */
    ISOLINE = 484,                 /* ISOLINE  */
    ISOLINES = 485,                /* ISOLINES  */
    JUST = 486,                    /* JUST  */
    KILL = 487,                    /* KILL  */
    KM = 488,                      /* KM  */
    LABEL = 489,                   /* LABEL  */
    LAYOUT = 490,                  /* LAYOUT  */
    LEAVE = 491,                   /* LEAVE  */
    LEAVEGRAPHICS = 492,           /* LEAVEGRAPHICS  */
    LEFT = 493,                    /* LEFT  */
    LEGEND = 494,                  /* LEGEND  */
    LENGTH = 495,                  /* LENGTH  */
    LEVEL = 496,                   /* LEVEL  */
    LEVELS = 497,                  /* LEVELS  */
    LIMITS = 498,                  /* LIMITS  */
    LINE = 499,                    /* LINE  */
    LINES = 500,                   /* LINES  */
    LINESTYLE = 501,               /* LINESTYLE  */
    LINETO = 502,                  /* LINETO  */
    LINEW = 503,                   /* LINEW  */
    LINEWIDTH = 504,               /* LINEWIDTH  */
    LINK = 505,                    /* LINK  */
    LOAD = 506,                    /* LOAD  */
    LOC = 507,                     /* LOC  */
    LOCATE = 508,                  /* LOCATE  */
    LOCATOR = 509,                 /* LOCATOR  */
    LOCTYPE = 510,                 /* LOCTYPE  */
    LOGX = 511,                    /* LOGX  */
    LOGXY = 512,                   /* LOGXY  */
    LOGY = 513,                    /* LOGY  */
    M = 514,                       /* M  */
    MAG = 515,                     /* MAG  */
    MAGNITUDE = 516,               /* MAGNITUDE  */
    MAJOR = 517,                   /* MAJOR  */
    MAPSCALE = 518,                /* MAPSCALE  */
    MARKER = 519,                  /* MARKER  */
    MARKERS = 520,                 /* MARKERS  */
    MAXLEVELS = 521,               /* MAXLEVELS  */
    METHOD = 522,                  /* METHOD  */
    MIFL = 523,                    /* MIFL  */
    MIFP = 524,                    /* MIFP  */
    MILES = 525,                   /* MILES  */
    MINOR = 526,                   /* MINOR  */
    MINUTES = 527,                 /* MINUTES  */
    MISSINGP = 528,                /* MISSINGP  */
    MM = 529,                      /* MM  */
    MMDD = 530,                    /* MMDD  */
    MMDDHMS = 531,                 /* MMDDHMS  */
    MMDDYY = 532,                  /* MMDDYY  */
    MMDDYYHMS = 533,               /* MMDDYYHMS  */
    MMSSLAT = 534,                 /* MMSSLAT  */
    MMSSLON = 535,                 /* MMSSLON  */
    MMYY = 536,                    /* MMYY  */
    MONTHDAY = 537,                /* MONTHDAY  */
    MONTHL = 538,                  /* MONTHL  */
    MONTHS = 539,                  /* MONTHS  */
    MOVE = 540,                    /* MOVE  */
    MOVE2 = 541,                   /* MOVE2  */
    MOVETO = 542,                  /* MOVETO  */
    NEGATE = 543,                  /* NEGATE  */
    NO = 544,                      /* NO  */
    NODE = 545,                    /* NODE  */
    NODES = 546,                   /* NODES  */
    NONE = 547,                    /* NONE  */
    NORMAL = 548,                  /* NORMAL  */
    NORTH = 549,                   /* NORTH  */
    NXY = 550,                     /* NXY  */
    OFF = 551,                     /* OFF  */
    OFFSETX = 552,                 /* OFFSETX  */
    OFFSETY = 553,                 /* OFFSETY  */
    ON = 554,                      /* ON  */
    OP = 555,                      /* OP  */
    OPEN = 556,                    /* OPEN  */
    ORIENT = 557,                  /* ORIENT  */
    OUT = 558,                     /* OUT  */
    PAGE = 559,                    /* PAGE  */
    PARA = 560,                    /* PARA  */
    PARALLEL = 561,                /* PARALLEL  */
    PARAMETERS = 562,              /* PARAMETERS  */
    PARAMS = 563,                  /* PARAMS  */
    PARMS = 564,                   /* PARMS  */
    PATTERN = 565,                 /* PATTERN  */
    PER = 566,                     /* PER  */
    PERIMETER = 567,               /* PERIMETER  */
    PERP = 568,                    /* PERP  */
    PERPENDICULAR = 569,           /* PERPENDICULAR  */
    PHASE = 570,                   /* PHASE  */
    PIE = 571,                     /* PIE  */
    PIPE = 572,                    /* PIPE  */
    PLACE = 573,                   /* PLACE  */
    PLAN = 574,                    /* PLAN  */
    PLUS = 575,                    /* PLUS  */
    POINT = 576,                   /* POINT  */
    POLAR = 577,                   /* POLAR  */
    POLY = 578,                    /* POLY  */
    POLYI = 579,                   /* POLYI  */
    POLYO = 580,                   /* POLYO  */
    POP = 581,                     /* POP  */
    POWER = 582,                   /* POWER  */
    PREC = 583,                    /* PREC  */
    PREFIX = 584,                  /* PREFIX  */
    PREPEND = 585,                 /* PREPEND  */
    PRINT = 586,                   /* PRINT  */
    PROFILE = 587,                 /* PROFILE  */
    PROP = 588,                    /* PROP  */
    PS = 589,                      /* PS  */
    PSCOLORL = 590,                /* PSCOLORL  */
    PSCOLORP = 591,                /* PSCOLORP  */
    PSMONOL = 592,                 /* PSMONOL  */
    PSMONOP = 593,                 /* PSMONOP  */
    PUSH = 594,                    /* PUSH  */
    PUTP = 595,                    /* PUTP  */
    QUIT = 596,                    /* QUIT  */
    READ = 597,                    /* READ  */
    READBIN = 598,                 /* READBIN  */
    REDRAW = 599,                  /* REDRAW  */
    REGION = 600,                  /* REGION  */
    REGIONS = 601,                 /* REGIONS  */
    REGNUM = 602,                  /* REGNUM  */
    REGRESS = 603,                 /* REGRESS  */
    REMOVE = 604,                  /* REMOVE  */
    RENDER = 605,                  /* RENDER  */
    REPORT = 606,                  /* REPORT  */
    RESET = 607,                   /* RESET  */
    REVERSE = 608,                 /* REVERSE  */
    REWIND = 609,                  /* REWIND  */
    RIGHT = 610,                   /* RIGHT  */
    RISER = 611,                   /* RISER  */
    ROT = 612,                     /* ROT  */
    RUN = 613,                     /* RUN  */
    SALINITY = 614,                /* SALINITY  */
    SAMPLE = 615,                  /* SAMPLE  */
    SAVE = 616,                    /* SAVE  */
    SCALAR = 617,                  /* SCALAR  */
    SCALE = 618,                   /* SCALE  */
    SCIENTIFIC = 619,              /* SCIENTIFIC  */
    SECONDS = 620,                 /* SECONDS  */
    SET = 621,                     /* SET  */
    SETS = 622,                    /* SETS  */
    SHOW = 623,                    /* SHOW  */
    SHRINK = 624,                  /* SHRINK  */
    SIGMA = 625,                   /* SIGMA  */
    SIGN = 626,                    /* SIGN  */
    SIZE = 627,                    /* SIZE  */
    SKIP = 628,                    /* SKIP  */
    SLAB = 629,                    /* SLAB  */
    SLEEP = 630,                   /* SLEEP  */
    SLICE = 631,                   /* SLICE  */
    SOURCE = 632,                  /* SOURCE  */
    SPEC = 633,                    /* SPEC  */
    SPECIFIED = 634,               /* SPECIFIED  */
    SPECTRUM = 635,                /* SPECTRUM  */
    SPLITS = 636,                  /* SPLITS  */
    SQUARE = 637,                  /* SQUARE  */
    STACK = 638,                   /* STACK  */
    STACKEDBAR = 639,              /* STACKEDBAR  */
    STACKEDHBAR = 640,             /* STACKEDHBAR  */
    STACKEDLINE = 641,             /* STACKEDLINE  */
    STAGGER = 642,                 /* STAGGER  */
    STAR = 643,                    /* STAR  */
    START = 644,                   /* START  */
    STARTSTEP = 645,               /* STARTSTEP  */
    STARTTYPE = 646,               /* STARTTYPE  */
    STATION = 647,                 /* STATION  */
    STATUS = 648,                  /* STATUS  */
    STEP = 649,                    /* STEP  */
    STOP = 650,                    /* STOP  */
    STREAMLINES = 651,             /* STREAMLINES  */
    STRING = 652,                  /* STRING  */
    STRINGS = 653,                 /* STRINGS  */
    SUBTITLE = 654,                /* SUBTITLE  */
    SURFACE = 655,                 /* SURFACE  */
    SWAPBUFFER = 656,              /* SWAPBUFFER  */
    SYMBOL = 657,                  /* SYMBOL  */
    SYSTEM = 658,                  /* SYSTEM  */
    TEANL = 659,                   /* TEANL  */
    TEXT = 660,                    /* TEXT  */
    TICK = 661,                    /* TICK  */
    TICKLABEL = 662,               /* TICKLABEL  */
    TICKMARKS = 663,               /* TICKMARKS  */
    TICKP = 664,                   /* TICKP  */
    TIDALCLOCK = 665,              /* TIDALCLOCK  */
    TIDESTATION = 666,             /* TIDESTATION  */
    TIME = 667,                    /* TIME  */
    TIMEINFO = 668,                /* TIMEINFO  */
    TIMELINE = 669,                /* TIMELINE  */
    TITLE = 670,                   /* TITLE  */
    TO = 671,                      /* TO  */
    TOP = 672,                     /* TOP  */
    TOTAL = 673,                   /* TOTAL  */
    TRACK = 674,                   /* TRACK  */
    TRANSECT = 675,                /* TRANSECT  */
    TRIANGLE1 = 676,               /* TRIANGLE1  */
    TRIANGLE2 = 677,               /* TRIANGLE2  */
    TRIANGLE3 = 678,               /* TRIANGLE3  */
    TRIANGLE4 = 679,               /* TRIANGLE4  */
    TRUEP = 680,                   /* TRUEP  */
    TYPE = 681,                    /* TYPE  */
    UNITS = 682,                   /* UNITS  */
    UP = 683,                      /* UP  */
    VALUE = 684,                   /* VALUE  */
    VECTOR = 685,                  /* VECTOR  */
    VEL = 686,                     /* VEL  */
    VELMARKER = 687,               /* VELMARKER  */
    VELOCITY = 688,                /* VELOCITY  */
    VERTICAL = 689,                /* VERTICAL  */
    VGAP = 690,                    /* VGAP  */
    VIEW = 691,                    /* VIEW  */
    VSCALE = 692,                  /* VSCALE  */
    VX1 = 693,                     /* VX1  */
    VX2 = 694,                     /* VX2  */
    VY1 = 695,                     /* VY1  */
    VY2 = 696,                     /* VY2  */
    WEEKS = 697,                   /* WEEKS  */
    WET = 698,                     /* WET  */
    WETDRY = 699,                  /* WETDRY  */
    WIDTH = 700,                   /* WIDTH  */
    WIND = 701,                    /* WIND  */
    WITH = 702,                    /* WITH  */
    WORLD = 703,                   /* WORLD  */
    WRAP = 704,                    /* WRAP  */
    WRITE = 705,                   /* WRITE  */
    WSCALE = 706,                  /* WSCALE  */
    WX1 = 707,                     /* WX1  */
    WX2 = 708,                     /* WX2  */
    WY1 = 709,                     /* WY1  */
    WY2 = 710,                     /* WY2  */
    X0 = 711,                      /* X0  */
    X1 = 712,                      /* X1  */
    X2 = 713,                      /* X2  */
    X3 = 714,                      /* X3  */
    X4 = 715,                      /* X4  */
    X5 = 716,                      /* X5  */
    XAXES = 717,                   /* XAXES  */
    XAXIS = 718,                   /* XAXIS  */
    XCOR = 719,                    /* XCOR  */
    XMAX = 720,                    /* XMAX  */
    XMIN = 721,                    /* XMIN  */
    XY = 722,                      /* XY  */
    XYARC = 723,                   /* XYARC  */
    XYBOX = 724,                   /* XYBOX  */
    XYDX = 725,                    /* XYDX  */
    XYDXDX = 726,                  /* XYDXDX  */
    XYDXDY = 727,                  /* XYDXDY  */
    XYDY = 728,                    /* XYDY  */
    XYDYDY = 729,                  /* XYDYDY  */
    XYFIXED = 730,                 /* XYFIXED  */
    XYHILO = 731,                  /* XYHILO  */
    XYRT = 732,                    /* XYRT  */
    XYSEG = 733,                   /* XYSEG  */
    XYSTRING = 734,                /* XYSTRING  */
    XYUV = 735,                    /* XYUV  */
    XYX2Y2 = 736,                  /* XYX2Y2  */
    XYXX = 737,                    /* XYXX  */
    XYYY = 738,                    /* XYYY  */
    XYZ = 739,                     /* XYZ  */
    XYZW = 740,                    /* XYZW  */
    Y0 = 741,                      /* Y0  */
    Y1 = 742,                      /* Y1  */
    Y2 = 743,                      /* Y2  */
    Y3 = 744,                      /* Y3  */
    Y4 = 745,                      /* Y4  */
    Y5 = 746,                      /* Y5  */
    YAXES = 747,                   /* YAXES  */
    YAXIS = 748,                   /* YAXIS  */
    YEARS = 749,                   /* YEARS  */
    YES = 750,                     /* YES  */
    YMAX = 751,                    /* YMAX  */
    YMIN = 752,                    /* YMIN  */
    ZEROXAXIS = 753,               /* ZEROXAXIS  */
    ZEROYAXIS = 754,               /* ZEROYAXIS  */
    ZOOM = 755,                    /* ZOOM  */
    ZOOMBOX = 756,                 /* ZOOMBOX  */
    OR = 757,                      /* OR  */
    AND = 758,                     /* AND  */
    GT = 759,                      /* GT  */
    LT = 760,                      /* LT  */
    LE = 761,                      /* LE  */
    GE = 762,                      /* GE  */
    EQ = 763,                      /* EQ  */
    NE = 764,                      /* NE  */
    UMINUS = 765,                  /* UMINUS  */
    NOT = 766                      /* NOT  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif
/* Token kinds.  */
#define YYEMPTY -2
#define YYEOF 0
#define YYerror 256
#define YYUNDEF 257
#define VAR 258
#define X 259
#define Y 260
#define CHRSTR 261
#define FITPARM 262
#define NUMBER 263
#define ABS 264
#define ACOS 265
#define ASIN 266
#define ATAN 267
#define ATAN2 268
#define CEIL 269
#define COS 270
#define DEG 271
#define DX 272
#define DY 273
#define ERF 274
#define ERFC 275
#define EXP 276
#define FLOOR 277
#define HYPOT 278
#define INDEX 279
#define INT 280
#define IRAND 281
#define LGAMMA 282
#define LN 283
#define LOG 284
#define LOGISTIC 285
#define MAXP 286
#define MINP 287
#define MINMAX 288
#define MOD 289
#define NORM 290
#define NORMP 291
#define PI 292
#define RAD 293
#define RAND 294
#define SETNO 295
#define SIN 296
#define SQR 297
#define SQRT 298
#define TAN 299
#define INUM 300
#define ABORT 301
#define ABOVE 302
#define ABSOLUTE 303
#define ACTIVATE 304
#define ACTIVE 305
#define ADCIRC 306
#define ADCIRC3DFLOW 307
#define ALL 308
#define ALT 309
#define ALTERNATE 310
#define ALTXAXIS 311
#define ALTYAXIS 312
#define AMP 313
#define ANGLE 314
#define ANNOTATE 315
#define APPEND 316
#define AREA 317
#define ARROW 318
#define ASCEND 319
#define AT 320
#define ATTACH 321
#define AUTO 322
#define AUTOSCALE 323
#define AUTOTICKS 324
#define AVERAGE 325
#define AVG 326
#define AXES 327
#define AXIS 328
#define BACKBUFFER 329
#define BACKGROUND 330
#define BAR 331
#define BATCH 332
#define BATH 333
#define BATHYMETRY 334
#define COURANT 335
#define BELOW 336
#define BIN 337
#define BINARY 338
#define BOTH 339
#define BOTTOM 340
#define BOUNDARY 341
#define BOX 342
#define CELLS 343
#define CENTER 344
#define CH3D 345
#define CHAR 346
#define CHDIR 347
#define CIRCLE 348
#define CLEAR 349
#define CLICK 350
#define CLOCK 351
#define CLOSE 352
#define CM 353
#define CMAP 354
#define COLOR 355
#define COLORMAP 356
#define COMMENT 357
#define CONC 358
#define CONCENTRATION 359
#define CONCENTRATIONS 360
#define COPY 361
#define CROSS 362
#define CYCLE 363
#define DAYMONTH 364
#define DAYOFWEEKL 365
#define DAYOFWEEKS 366
#define DAYOFYEAR 367
#define DAYS 368
#define DDMMYY 369
#define DDMONTHSYYHHMMSS 370
#define DECIMAL 371
#define DEF 372
#define DEFAULT 373
#define DEGREESLAT 374
#define DEGREESLON 375
#define DEGREESMMLAT 376
#define DEGREESMMLON 377
#define DEGREESMMSSLAT 378
#define DEGREESMMSSLON 379
#define DELAYP 380
#define DELETE 381
#define DEPTH 382
#define DEPTHS 383
#define DESCEND 384
#define DEVICE 385
#define DEVXY 386
#define DFT 387
#define DT 388
#define DIAMOND 389
#define DIFFERENCE 390
#define DISK 391
#define DISPLAY 392
#define DOT 393
#define DOUBLEBUFFER 394
#define DOWN 395
#define DRAW2 396
#define DROGUE 397
#define DROGUES 398
#define DRY 399
#define DXDX 400
#define DXP 401
#define DYDY 402
#define DYP 403
#define ECHO 404
#define EDIT 405
#define ELA 406
#define ELCIRC 407
#define ELEMENT 408
#define ELEMENTS 409
#define ELEV 410
#define ELEVATION 411
#define ELEVATIONS 412
#define ELEVMARKER 413
#define ELLIPSE 414
#define ELLIPSES 415
#define ELLIPSEZ 416
#define ELSE 417
#define END 418
#define ERRORBAR 419
#define EXIT 420
#define EXPAND 421
#define EXPONENTIAL 422
#define FACTOR 423
#define FALSEP 424
#define FAST 425
#define FEET 426
#define FFT 427
#define FILEP 428
#define FILL 429
#define FIND 430
#define FIXEDPOINT 431
#define FLOW 432
#define FLUSH 433
#define FLUX 434
#define FOCUS 435
#define FOLLOWS 436
#define FONTP 437
#define FOREGROUND 438
#define FORMAT 439
#define FORT14 440
#define FORT63 441
#define FORT64 442
#define FORWARD 443
#define FRAMEP 444
#define FREQ 445
#define FRONTBUFFER 446
#define GENERAL 447
#define GETP 448
#define GOTO 449
#define GRAPH 450
#define GRAPHNO 451
#define GRAPHS 452
#define GRAPHTYPE 453
#define GRID 454
#define HARDCOPY 455
#define HBAR 456
#define HELP 457
#define HGAP 458
#define HIDDEN 459
#define HISTBOX 460
#define HISTO 461
#define HISTORY 462
#define HMS 463
#define HORIZONTAL 464
#define HOURS 465
#define HPGLL 466
#define HPGLP 467
#define IF 468
#define IGNORE 469
#define IHL 470
#define IMAGE 471
#define IMAGES 472
#define IN 473
#define INCLUDE 474
#define INFO 475
#define INIT 476
#define INITGRAPHICS 477
#define INOUT 478
#define INTEGRATE 479
#define INTERP 480
#define INUNDATION 481
#define INVDFT 482
#define INVFFT 483
#define ISOLINE 484
#define ISOLINES 485
#define JUST 486
#define KILL 487
#define KM 488
#define LABEL 489
#define LAYOUT 490
#define LEAVE 491
#define LEAVEGRAPHICS 492
#define LEFT 493
#define LEGEND 494
#define LENGTH 495
#define LEVEL 496
#define LEVELS 497
#define LIMITS 498
#define LINE 499
#define LINES 500
#define LINESTYLE 501
#define LINETO 502
#define LINEW 503
#define LINEWIDTH 504
#define LINK 505
#define LOAD 506
#define LOC 507
#define LOCATE 508
#define LOCATOR 509
#define LOCTYPE 510
#define LOGX 511
#define LOGXY 512
#define LOGY 513
#define M 514
#define MAG 515
#define MAGNITUDE 516
#define MAJOR 517
#define MAPSCALE 518
#define MARKER 519
#define MARKERS 520
#define MAXLEVELS 521
#define METHOD 522
#define MIFL 523
#define MIFP 524
#define MILES 525
#define MINOR 526
#define MINUTES 527
#define MISSINGP 528
#define MM 529
#define MMDD 530
#define MMDDHMS 531
#define MMDDYY 532
#define MMDDYYHMS 533
#define MMSSLAT 534
#define MMSSLON 535
#define MMYY 536
#define MONTHDAY 537
#define MONTHL 538
#define MONTHS 539
#define MOVE 540
#define MOVE2 541
#define MOVETO 542
#define NEGATE 543
#define NO 544
#define NODE 545
#define NODES 546
#define NONE 547
#define NORMAL 548
#define NORTH 549
#define NXY 550
#define OFF 551
#define OFFSETX 552
#define OFFSETY 553
#define ON 554
#define OP 555
#define OPEN 556
#define ORIENT 557
#define OUT 558
#define PAGE 559
#define PARA 560
#define PARALLEL 561
#define PARAMETERS 562
#define PARAMS 563
#define PARMS 564
#define PATTERN 565
#define PER 566
#define PERIMETER 567
#define PERP 568
#define PERPENDICULAR 569
#define PHASE 570
#define PIE 571
#define PIPE 572
#define PLACE 573
#define PLAN 574
#define PLUS 575
#define POINT 576
#define POLAR 577
#define POLY 578
#define POLYI 579
#define POLYO 580
#define POP 581
#define POWER 582
#define PREC 583
#define PREFIX 584
#define PREPEND 585
#define PRINT 586
#define PROFILE 587
#define PROP 588
#define PS 589
#define PSCOLORL 590
#define PSCOLORP 591
#define PSMONOL 592
#define PSMONOP 593
#define PUSH 594
#define PUTP 595
#define QUIT 596
#define READ 597
#define READBIN 598
#define REDRAW 599
#define REGION 600
#define REGIONS 601
#define REGNUM 602
#define REGRESS 603
#define REMOVE 604
#define RENDER 605
#define REPORT 606
#define RESET 607
#define REVERSE 608
#define REWIND 609
#define RIGHT 610
#define RISER 611
#define ROT 612
#define RUN 613
#define SALINITY 614
#define SAMPLE 615
#define SAVE 616
#define SCALAR 617
#define SCALE 618
#define SCIENTIFIC 619
#define SECONDS 620
#define SET 621
#define SETS 622
#define SHOW 623
#define SHRINK 624
#define SIGMA 625
#define SIGN 626
#define SIZE 627
#define SKIP 628
#define SLAB 629
#define SLEEP 630
#define SLICE 631
#define SOURCE 632
#define SPEC 633
#define SPECIFIED 634
#define SPECTRUM 635
#define SPLITS 636
#define SQUARE 637
#define STACK 638
#define STACKEDBAR 639
#define STACKEDHBAR 640
#define STACKEDLINE 641
#define STAGGER 642
#define STAR 643
#define START 644
#define STARTSTEP 645
#define STARTTYPE 646
#define STATION 647
#define STATUS 648
#define STEP 649
#define STOP 650
#define STREAMLINES 651
#define STRING 652
#define STRINGS 653
#define SUBTITLE 654
#define SURFACE 655
#define SWAPBUFFER 656
#define SYMBOL 657
#define SYSTEM 658
#define TEANL 659
#define TEXT 660
#define TICK 661
#define TICKLABEL 662
#define TICKMARKS 663
#define TICKP 664
#define TIDALCLOCK 665
#define TIDESTATION 666
#define TIME 667
#define TIMEINFO 668
#define TIMELINE 669
#define TITLE 670
#define TO 671
#define TOP 672
#define TOTAL 673
#define TRACK 674
#define TRANSECT 675
#define TRIANGLE1 676
#define TRIANGLE2 677
#define TRIANGLE3 678
#define TRIANGLE4 679
#define TRUEP 680
#define TYPE 681
#define UNITS 682
#define UP 683
#define VALUE 684
#define VECTOR 685
#define VEL 686
#define VELMARKER 687
#define VELOCITY 688
#define VERTICAL 689
#define VGAP 690
#define VIEW 691
#define VSCALE 692
#define VX1 693
#define VX2 694
#define VY1 695
#define VY2 696
#define WEEKS 697
#define WET 698
#define WETDRY 699
#define WIDTH 700
#define WIND 701
#define WITH 702
#define WORLD 703
#define WRAP 704
#define WRITE 705
#define WSCALE 706
#define WX1 707
#define WX2 708
#define WY1 709
#define WY2 710
#define X0 711
#define X1 712
#define X2 713
#define X3 714
#define X4 715
#define X5 716
#define XAXES 717
#define XAXIS 718
#define XCOR 719
#define XMAX 720
#define XMIN 721
#define XY 722
#define XYARC 723
#define XYBOX 724
#define XYDX 725
#define XYDXDX 726
#define XYDXDY 727
#define XYDY 728
#define XYDYDY 729
#define XYFIXED 730
#define XYHILO 731
#define XYRT 732
#define XYSEG 733
#define XYSTRING 734
#define XYUV 735
#define XYX2Y2 736
#define XYXX 737
#define XYYY 738
#define XYZ 739
#define XYZW 740
#define Y0 741
#define Y1 742
#define Y2 743
#define Y3 744
#define Y4 745
#define Y5 746
#define YAXES 747
#define YAXIS 748
#define YEARS 749
#define YES 750
#define YMAX 751
#define YMIN 752
#define ZEROXAXIS 753
#define ZEROYAXIS 754
#define ZOOM 755
#define ZOOMBOX 756
#define OR 757
#define AND 758
#define GT 759
#define LT 760
#define LE 761
#define GE 762
#define EQ 763
#define NE 764
#define UMINUS 765
#define NOT 766

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 113 "gram.y"

    double val;
    int ival;
    double *ptr;
    int func;
    int pset;
    char *str;

#line 1257 "y.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;

int yyparse (void);

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_VAR = 3,                        /* VAR  */
  YYSYMBOL_X = 4,                          /* X  */
  YYSYMBOL_Y = 5,                          /* Y  */
  YYSYMBOL_CHRSTR = 6,                     /* CHRSTR  */
  YYSYMBOL_FITPARM = 7,                    /* FITPARM  */
  YYSYMBOL_NUMBER = 8,                     /* NUMBER  */
  YYSYMBOL_ABS = 9,                        /* ABS  */
  YYSYMBOL_ACOS = 10,                      /* ACOS  */
  YYSYMBOL_ASIN = 11,                      /* ASIN  */
  YYSYMBOL_ATAN = 12,                      /* ATAN  */
  YYSYMBOL_ATAN2 = 13,                     /* ATAN2  */
  YYSYMBOL_CEIL = 14,                      /* CEIL  */
  YYSYMBOL_COS = 15,                       /* COS  */
  YYSYMBOL_DEG = 16,                       /* DEG  */
  YYSYMBOL_DX = 17,                        /* DX  */
  YYSYMBOL_DY = 18,                        /* DY  */
  YYSYMBOL_ERF = 19,                       /* ERF  */
  YYSYMBOL_ERFC = 20,                      /* ERFC  */
  YYSYMBOL_EXP = 21,                       /* EXP  */
  YYSYMBOL_FLOOR = 22,                     /* FLOOR  */
  YYSYMBOL_HYPOT = 23,                     /* HYPOT  */
  YYSYMBOL_INDEX = 24,                     /* INDEX  */
  YYSYMBOL_INT = 25,                       /* INT  */
  YYSYMBOL_IRAND = 26,                     /* IRAND  */
  YYSYMBOL_LGAMMA = 27,                    /* LGAMMA  */
  YYSYMBOL_LN = 28,                        /* LN  */
  YYSYMBOL_LOG = 29,                       /* LOG  */
  YYSYMBOL_LOGISTIC = 30,                  /* LOGISTIC  */
  YYSYMBOL_MAXP = 31,                      /* MAXP  */
  YYSYMBOL_MINP = 32,                      /* MINP  */
  YYSYMBOL_MINMAX = 33,                    /* MINMAX  */
  YYSYMBOL_MOD = 34,                       /* MOD  */
  YYSYMBOL_NORM = 35,                      /* NORM  */
  YYSYMBOL_NORMP = 36,                     /* NORMP  */
  YYSYMBOL_PI = 37,                        /* PI  */
  YYSYMBOL_RAD = 38,                       /* RAD  */
  YYSYMBOL_RAND = 39,                      /* RAND  */
  YYSYMBOL_SETNO = 40,                     /* SETNO  */
  YYSYMBOL_SIN = 41,                       /* SIN  */
  YYSYMBOL_SQR = 42,                       /* SQR  */
  YYSYMBOL_SQRT = 43,                      /* SQRT  */
  YYSYMBOL_TAN = 44,                       /* TAN  */
  YYSYMBOL_INUM = 45,                      /* INUM  */
  YYSYMBOL_ABORT = 46,                     /* ABORT  */
  YYSYMBOL_ABOVE = 47,                     /* ABOVE  */
  YYSYMBOL_ABSOLUTE = 48,                  /* ABSOLUTE  */
  YYSYMBOL_ACTIVATE = 49,                  /* ACTIVATE  */
  YYSYMBOL_ACTIVE = 50,                    /* ACTIVE  */
  YYSYMBOL_ADCIRC = 51,                    /* ADCIRC  */
  YYSYMBOL_ADCIRC3DFLOW = 52,              /* ADCIRC3DFLOW  */
  YYSYMBOL_ALL = 53,                       /* ALL  */
  YYSYMBOL_ALT = 54,                       /* ALT  */
  YYSYMBOL_ALTERNATE = 55,                 /* ALTERNATE  */
  YYSYMBOL_ALTXAXIS = 56,                  /* ALTXAXIS  */
  YYSYMBOL_ALTYAXIS = 57,                  /* ALTYAXIS  */
  YYSYMBOL_AMP = 58,                       /* AMP  */
  YYSYMBOL_ANGLE = 59,                     /* ANGLE  */
  YYSYMBOL_ANNOTATE = 60,                  /* ANNOTATE  */
  YYSYMBOL_APPEND = 61,                    /* APPEND  */
  YYSYMBOL_AREA = 62,                      /* AREA  */
  YYSYMBOL_ARROW = 63,                     /* ARROW  */
  YYSYMBOL_ASCEND = 64,                    /* ASCEND  */
  YYSYMBOL_AT = 65,                        /* AT  */
  YYSYMBOL_ATTACH = 66,                    /* ATTACH  */
  YYSYMBOL_AUTO = 67,                      /* AUTO  */
  YYSYMBOL_AUTOSCALE = 68,                 /* AUTOSCALE  */
  YYSYMBOL_AUTOTICKS = 69,                 /* AUTOTICKS  */
  YYSYMBOL_AVERAGE = 70,                   /* AVERAGE  */
  YYSYMBOL_AVG = 71,                       /* AVG  */
  YYSYMBOL_AXES = 72,                      /* AXES  */
  YYSYMBOL_AXIS = 73,                      /* AXIS  */
  YYSYMBOL_BACKBUFFER = 74,                /* BACKBUFFER  */
  YYSYMBOL_BACKGROUND = 75,                /* BACKGROUND  */
  YYSYMBOL_BAR = 76,                       /* BAR  */
  YYSYMBOL_BATCH = 77,                     /* BATCH  */
  YYSYMBOL_BATH = 78,                      /* BATH  */
  YYSYMBOL_BATHYMETRY = 79,                /* BATHYMETRY  */
  YYSYMBOL_COURANT = 80,                   /* COURANT  */
  YYSYMBOL_BELOW = 81,                     /* BELOW  */
  YYSYMBOL_BIN = 82,                       /* BIN  */
  YYSYMBOL_BINARY = 83,                    /* BINARY  */
  YYSYMBOL_BOTH = 84,                      /* BOTH  */
  YYSYMBOL_BOTTOM = 85,                    /* BOTTOM  */
  YYSYMBOL_BOUNDARY = 86,                  /* BOUNDARY  */
  YYSYMBOL_BOX = 87,                       /* BOX  */
  YYSYMBOL_CELLS = 88,                     /* CELLS  */
  YYSYMBOL_CENTER = 89,                    /* CENTER  */
  YYSYMBOL_CH3D = 90,                      /* CH3D  */
  YYSYMBOL_CHAR = 91,                      /* CHAR  */
  YYSYMBOL_CHDIR = 92,                     /* CHDIR  */
  YYSYMBOL_CIRCLE = 93,                    /* CIRCLE  */
  YYSYMBOL_CLEAR = 94,                     /* CLEAR  */
  YYSYMBOL_CLICK = 95,                     /* CLICK  */
  YYSYMBOL_CLOCK = 96,                     /* CLOCK  */
  YYSYMBOL_CLOSE = 97,                     /* CLOSE  */
  YYSYMBOL_CM = 98,                        /* CM  */
  YYSYMBOL_CMAP = 99,                      /* CMAP  */
  YYSYMBOL_COLOR = 100,                    /* COLOR  */
  YYSYMBOL_COLORMAP = 101,                 /* COLORMAP  */
  YYSYMBOL_COMMENT = 102,                  /* COMMENT  */
  YYSYMBOL_CONC = 103,                     /* CONC  */
  YYSYMBOL_CONCENTRATION = 104,            /* CONCENTRATION  */
  YYSYMBOL_CONCENTRATIONS = 105,           /* CONCENTRATIONS  */
  YYSYMBOL_COPY = 106,                     /* COPY  */
  YYSYMBOL_CROSS = 107,                    /* CROSS  */
  YYSYMBOL_CYCLE = 108,                    /* CYCLE  */
  YYSYMBOL_DAYMONTH = 109,                 /* DAYMONTH  */
  YYSYMBOL_DAYOFWEEKL = 110,               /* DAYOFWEEKL  */
  YYSYMBOL_DAYOFWEEKS = 111,               /* DAYOFWEEKS  */
  YYSYMBOL_DAYOFYEAR = 112,                /* DAYOFYEAR  */
  YYSYMBOL_DAYS = 113,                     /* DAYS  */
  YYSYMBOL_DDMMYY = 114,                   /* DDMMYY  */
  YYSYMBOL_DDMONTHSYYHHMMSS = 115,         /* DDMONTHSYYHHMMSS  */
  YYSYMBOL_DECIMAL = 116,                  /* DECIMAL  */
  YYSYMBOL_DEF = 117,                      /* DEF  */
  YYSYMBOL_DEFAULT = 118,                  /* DEFAULT  */
  YYSYMBOL_DEGREESLAT = 119,               /* DEGREESLAT  */
  YYSYMBOL_DEGREESLON = 120,               /* DEGREESLON  */
  YYSYMBOL_DEGREESMMLAT = 121,             /* DEGREESMMLAT  */
  YYSYMBOL_DEGREESMMLON = 122,             /* DEGREESMMLON  */
  YYSYMBOL_DEGREESMMSSLAT = 123,           /* DEGREESMMSSLAT  */
  YYSYMBOL_DEGREESMMSSLON = 124,           /* DEGREESMMSSLON  */
  YYSYMBOL_DELAYP = 125,                   /* DELAYP  */
  YYSYMBOL_DELETE = 126,                   /* DELETE  */
  YYSYMBOL_DEPTH = 127,                    /* DEPTH  */
  YYSYMBOL_DEPTHS = 128,                   /* DEPTHS  */
  YYSYMBOL_DESCEND = 129,                  /* DESCEND  */
  YYSYMBOL_DEVICE = 130,                   /* DEVICE  */
  YYSYMBOL_DEVXY = 131,                    /* DEVXY  */
  YYSYMBOL_DFT = 132,                      /* DFT  */
  YYSYMBOL_DT = 133,                       /* DT  */
  YYSYMBOL_DIAMOND = 134,                  /* DIAMOND  */
  YYSYMBOL_DIFFERENCE = 135,               /* DIFFERENCE  */
  YYSYMBOL_DISK = 136,                     /* DISK  */
  YYSYMBOL_DISPLAY = 137,                  /* DISPLAY  */
  YYSYMBOL_DOT = 138,                      /* DOT  */
  YYSYMBOL_DOUBLEBUFFER = 139,             /* DOUBLEBUFFER  */
  YYSYMBOL_DOWN = 140,                     /* DOWN  */
  YYSYMBOL_DRAW2 = 141,                    /* DRAW2  */
  YYSYMBOL_DROGUE = 142,                   /* DROGUE  */
  YYSYMBOL_DROGUES = 143,                  /* DROGUES  */
  YYSYMBOL_DRY = 144,                      /* DRY  */
  YYSYMBOL_DXDX = 145,                     /* DXDX  */
  YYSYMBOL_DXP = 146,                      /* DXP  */
  YYSYMBOL_DYDY = 147,                     /* DYDY  */
  YYSYMBOL_DYP = 148,                      /* DYP  */
  YYSYMBOL_ECHO = 149,                     /* ECHO  */
  YYSYMBOL_EDIT = 150,                     /* EDIT  */
  YYSYMBOL_ELA = 151,                      /* ELA  */
  YYSYMBOL_ELCIRC = 152,                   /* ELCIRC  */
  YYSYMBOL_ELEMENT = 153,                  /* ELEMENT  */
  YYSYMBOL_ELEMENTS = 154,                 /* ELEMENTS  */
  YYSYMBOL_ELEV = 155,                     /* ELEV  */
  YYSYMBOL_ELEVATION = 156,                /* ELEVATION  */
  YYSYMBOL_ELEVATIONS = 157,               /* ELEVATIONS  */
  YYSYMBOL_ELEVMARKER = 158,               /* ELEVMARKER  */
  YYSYMBOL_ELLIPSE = 159,                  /* ELLIPSE  */
  YYSYMBOL_ELLIPSES = 160,                 /* ELLIPSES  */
  YYSYMBOL_ELLIPSEZ = 161,                 /* ELLIPSEZ  */
  YYSYMBOL_ELSE = 162,                     /* ELSE  */
  YYSYMBOL_END = 163,                      /* END  */
  YYSYMBOL_ERRORBAR = 164,                 /* ERRORBAR  */
  YYSYMBOL_EXIT = 165,                     /* EXIT  */
  YYSYMBOL_EXPAND = 166,                   /* EXPAND  */
  YYSYMBOL_EXPONENTIAL = 167,              /* EXPONENTIAL  */
  YYSYMBOL_FACTOR = 168,                   /* FACTOR  */
  YYSYMBOL_FALSEP = 169,                   /* FALSEP  */
  YYSYMBOL_FAST = 170,                     /* FAST  */
  YYSYMBOL_FEET = 171,                     /* FEET  */
  YYSYMBOL_FFT = 172,                      /* FFT  */
  YYSYMBOL_FILEP = 173,                    /* FILEP  */
  YYSYMBOL_FILL = 174,                     /* FILL  */
  YYSYMBOL_FIND = 175,                     /* FIND  */
  YYSYMBOL_FIXEDPOINT = 176,               /* FIXEDPOINT  */
  YYSYMBOL_FLOW = 177,                     /* FLOW  */
  YYSYMBOL_FLUSH = 178,                    /* FLUSH  */
  YYSYMBOL_FLUX = 179,                     /* FLUX  */
  YYSYMBOL_FOCUS = 180,                    /* FOCUS  */
  YYSYMBOL_FOLLOWS = 181,                  /* FOLLOWS  */
  YYSYMBOL_FONTP = 182,                    /* FONTP  */
  YYSYMBOL_FOREGROUND = 183,               /* FOREGROUND  */
  YYSYMBOL_FORMAT = 184,                   /* FORMAT  */
  YYSYMBOL_FORT14 = 185,                   /* FORT14  */
  YYSYMBOL_FORT63 = 186,                   /* FORT63  */
  YYSYMBOL_FORT64 = 187,                   /* FORT64  */
  YYSYMBOL_FORWARD = 188,                  /* FORWARD  */
  YYSYMBOL_FRAMEP = 189,                   /* FRAMEP  */
  YYSYMBOL_FREQ = 190,                     /* FREQ  */
  YYSYMBOL_FRONTBUFFER = 191,              /* FRONTBUFFER  */
  YYSYMBOL_GENERAL = 192,                  /* GENERAL  */
  YYSYMBOL_GETP = 193,                     /* GETP  */
  YYSYMBOL_GOTO = 194,                     /* GOTO  */
  YYSYMBOL_GRAPH = 195,                    /* GRAPH  */
  YYSYMBOL_GRAPHNO = 196,                  /* GRAPHNO  */
  YYSYMBOL_GRAPHS = 197,                   /* GRAPHS  */
  YYSYMBOL_GRAPHTYPE = 198,                /* GRAPHTYPE  */
  YYSYMBOL_GRID = 199,                     /* GRID  */
  YYSYMBOL_HARDCOPY = 200,                 /* HARDCOPY  */
  YYSYMBOL_HBAR = 201,                     /* HBAR  */
  YYSYMBOL_HELP = 202,                     /* HELP  */
  YYSYMBOL_HGAP = 203,                     /* HGAP  */
  YYSYMBOL_HIDDEN = 204,                   /* HIDDEN  */
  YYSYMBOL_HISTBOX = 205,                  /* HISTBOX  */
  YYSYMBOL_HISTO = 206,                    /* HISTO  */
  YYSYMBOL_HISTORY = 207,                  /* HISTORY  */
  YYSYMBOL_HMS = 208,                      /* HMS  */
  YYSYMBOL_HORIZONTAL = 209,               /* HORIZONTAL  */
  YYSYMBOL_HOURS = 210,                    /* HOURS  */
  YYSYMBOL_HPGLL = 211,                    /* HPGLL  */
  YYSYMBOL_HPGLP = 212,                    /* HPGLP  */
  YYSYMBOL_IF = 213,                       /* IF  */
  YYSYMBOL_IGNORE = 214,                   /* IGNORE  */
  YYSYMBOL_IHL = 215,                      /* IHL  */
  YYSYMBOL_IMAGE = 216,                    /* IMAGE  */
  YYSYMBOL_IMAGES = 217,                   /* IMAGES  */
  YYSYMBOL_IN = 218,                       /* IN  */
  YYSYMBOL_INCLUDE = 219,                  /* INCLUDE  */
  YYSYMBOL_INFO = 220,                     /* INFO  */
  YYSYMBOL_INIT = 221,                     /* INIT  */
  YYSYMBOL_INITGRAPHICS = 222,             /* INITGRAPHICS  */
  YYSYMBOL_INOUT = 223,                    /* INOUT  */
  YYSYMBOL_INTEGRATE = 224,                /* INTEGRATE  */
  YYSYMBOL_INTERP = 225,                   /* INTERP  */
  YYSYMBOL_INUNDATION = 226,               /* INUNDATION  */
  YYSYMBOL_INVDFT = 227,                   /* INVDFT  */
  YYSYMBOL_INVFFT = 228,                   /* INVFFT  */
  YYSYMBOL_ISOLINE = 229,                  /* ISOLINE  */
  YYSYMBOL_ISOLINES = 230,                 /* ISOLINES  */
  YYSYMBOL_JUST = 231,                     /* JUST  */
  YYSYMBOL_KILL = 232,                     /* KILL  */
  YYSYMBOL_KM = 233,                       /* KM  */
  YYSYMBOL_LABEL = 234,                    /* LABEL  */
  YYSYMBOL_LAYOUT = 235,                   /* LAYOUT  */
  YYSYMBOL_LEAVE = 236,                    /* LEAVE  */
  YYSYMBOL_LEAVEGRAPHICS = 237,            /* LEAVEGRAPHICS  */
  YYSYMBOL_LEFT = 238,                     /* LEFT  */
  YYSYMBOL_LEGEND = 239,                   /* LEGEND  */
  YYSYMBOL_LENGTH = 240,                   /* LENGTH  */
  YYSYMBOL_LEVEL = 241,                    /* LEVEL  */
  YYSYMBOL_LEVELS = 242,                   /* LEVELS  */
  YYSYMBOL_LIMITS = 243,                   /* LIMITS  */
  YYSYMBOL_LINE = 244,                     /* LINE  */
  YYSYMBOL_LINES = 245,                    /* LINES  */
  YYSYMBOL_LINESTYLE = 246,                /* LINESTYLE  */
  YYSYMBOL_LINETO = 247,                   /* LINETO  */
  YYSYMBOL_LINEW = 248,                    /* LINEW  */
  YYSYMBOL_LINEWIDTH = 249,                /* LINEWIDTH  */
  YYSYMBOL_LINK = 250,                     /* LINK  */
  YYSYMBOL_LOAD = 251,                     /* LOAD  */
  YYSYMBOL_LOC = 252,                      /* LOC  */
  YYSYMBOL_LOCATE = 253,                   /* LOCATE  */
  YYSYMBOL_LOCATOR = 254,                  /* LOCATOR  */
  YYSYMBOL_LOCTYPE = 255,                  /* LOCTYPE  */
  YYSYMBOL_LOGX = 256,                     /* LOGX  */
  YYSYMBOL_LOGXY = 257,                    /* LOGXY  */
  YYSYMBOL_LOGY = 258,                     /* LOGY  */
  YYSYMBOL_M = 259,                        /* M  */
  YYSYMBOL_MAG = 260,                      /* MAG  */
  YYSYMBOL_MAGNITUDE = 261,                /* MAGNITUDE  */
  YYSYMBOL_MAJOR = 262,                    /* MAJOR  */
  YYSYMBOL_MAPSCALE = 263,                 /* MAPSCALE  */
  YYSYMBOL_MARKER = 264,                   /* MARKER  */
  YYSYMBOL_MARKERS = 265,                  /* MARKERS  */
  YYSYMBOL_MAXLEVELS = 266,                /* MAXLEVELS  */
  YYSYMBOL_METHOD = 267,                   /* METHOD  */
  YYSYMBOL_MIFL = 268,                     /* MIFL  */
  YYSYMBOL_MIFP = 269,                     /* MIFP  */
  YYSYMBOL_MILES = 270,                    /* MILES  */
  YYSYMBOL_MINOR = 271,                    /* MINOR  */
  YYSYMBOL_MINUTES = 272,                  /* MINUTES  */
  YYSYMBOL_MISSINGP = 273,                 /* MISSINGP  */
  YYSYMBOL_MM = 274,                       /* MM  */
  YYSYMBOL_MMDD = 275,                     /* MMDD  */
  YYSYMBOL_MMDDHMS = 276,                  /* MMDDHMS  */
  YYSYMBOL_MMDDYY = 277,                   /* MMDDYY  */
  YYSYMBOL_MMDDYYHMS = 278,                /* MMDDYYHMS  */
  YYSYMBOL_MMSSLAT = 279,                  /* MMSSLAT  */
  YYSYMBOL_MMSSLON = 280,                  /* MMSSLON  */
  YYSYMBOL_MMYY = 281,                     /* MMYY  */
  YYSYMBOL_MONTHDAY = 282,                 /* MONTHDAY  */
  YYSYMBOL_MONTHL = 283,                   /* MONTHL  */
  YYSYMBOL_MONTHS = 284,                   /* MONTHS  */
  YYSYMBOL_MOVE = 285,                     /* MOVE  */
  YYSYMBOL_MOVE2 = 286,                    /* MOVE2  */
  YYSYMBOL_MOVETO = 287,                   /* MOVETO  */
  YYSYMBOL_NEGATE = 288,                   /* NEGATE  */
  YYSYMBOL_NO = 289,                       /* NO  */
  YYSYMBOL_NODE = 290,                     /* NODE  */
  YYSYMBOL_NODES = 291,                    /* NODES  */
  YYSYMBOL_NONE = 292,                     /* NONE  */
  YYSYMBOL_NORMAL = 293,                   /* NORMAL  */
  YYSYMBOL_NORTH = 294,                    /* NORTH  */
  YYSYMBOL_NXY = 295,                      /* NXY  */
  YYSYMBOL_OFF = 296,                      /* OFF  */
  YYSYMBOL_OFFSETX = 297,                  /* OFFSETX  */
  YYSYMBOL_OFFSETY = 298,                  /* OFFSETY  */
  YYSYMBOL_ON = 299,                       /* ON  */
  YYSYMBOL_OP = 300,                       /* OP  */
  YYSYMBOL_OPEN = 301,                     /* OPEN  */
  YYSYMBOL_ORIENT = 302,                   /* ORIENT  */
  YYSYMBOL_OUT = 303,                      /* OUT  */
  YYSYMBOL_PAGE = 304,                     /* PAGE  */
  YYSYMBOL_PARA = 305,                     /* PARA  */
  YYSYMBOL_PARALLEL = 306,                 /* PARALLEL  */
  YYSYMBOL_PARAMETERS = 307,               /* PARAMETERS  */
  YYSYMBOL_PARAMS = 308,                   /* PARAMS  */
  YYSYMBOL_PARMS = 309,                    /* PARMS  */
  YYSYMBOL_PATTERN = 310,                  /* PATTERN  */
  YYSYMBOL_PER = 311,                      /* PER  */
  YYSYMBOL_PERIMETER = 312,                /* PERIMETER  */
  YYSYMBOL_PERP = 313,                     /* PERP  */
  YYSYMBOL_PERPENDICULAR = 314,            /* PERPENDICULAR  */
  YYSYMBOL_PHASE = 315,                    /* PHASE  */
  YYSYMBOL_PIE = 316,                      /* PIE  */
  YYSYMBOL_PIPE = 317,                     /* PIPE  */
  YYSYMBOL_PLACE = 318,                    /* PLACE  */
  YYSYMBOL_PLAN = 319,                     /* PLAN  */
  YYSYMBOL_PLUS = 320,                     /* PLUS  */
  YYSYMBOL_POINT = 321,                    /* POINT  */
  YYSYMBOL_POLAR = 322,                    /* POLAR  */
  YYSYMBOL_POLY = 323,                     /* POLY  */
  YYSYMBOL_POLYI = 324,                    /* POLYI  */
  YYSYMBOL_POLYO = 325,                    /* POLYO  */
  YYSYMBOL_POP = 326,                      /* POP  */
  YYSYMBOL_POWER = 327,                    /* POWER  */
  YYSYMBOL_PREC = 328,                     /* PREC  */
  YYSYMBOL_PREFIX = 329,                   /* PREFIX  */
  YYSYMBOL_PREPEND = 330,                  /* PREPEND  */
  YYSYMBOL_PRINT = 331,                    /* PRINT  */
  YYSYMBOL_PROFILE = 332,                  /* PROFILE  */
  YYSYMBOL_PROP = 333,                     /* PROP  */
  YYSYMBOL_PS = 334,                       /* PS  */
  YYSYMBOL_PSCOLORL = 335,                 /* PSCOLORL  */
  YYSYMBOL_PSCOLORP = 336,                 /* PSCOLORP  */
  YYSYMBOL_PSMONOL = 337,                  /* PSMONOL  */
  YYSYMBOL_PSMONOP = 338,                  /* PSMONOP  */
  YYSYMBOL_PUSH = 339,                     /* PUSH  */
  YYSYMBOL_PUTP = 340,                     /* PUTP  */
  YYSYMBOL_QUIT = 341,                     /* QUIT  */
  YYSYMBOL_READ = 342,                     /* READ  */
  YYSYMBOL_READBIN = 343,                  /* READBIN  */
  YYSYMBOL_REDRAW = 344,                   /* REDRAW  */
  YYSYMBOL_REGION = 345,                   /* REGION  */
  YYSYMBOL_REGIONS = 346,                  /* REGIONS  */
  YYSYMBOL_REGNUM = 347,                   /* REGNUM  */
  YYSYMBOL_REGRESS = 348,                  /* REGRESS  */
  YYSYMBOL_REMOVE = 349,                   /* REMOVE  */
  YYSYMBOL_RENDER = 350,                   /* RENDER  */
  YYSYMBOL_REPORT = 351,                   /* REPORT  */
  YYSYMBOL_RESET = 352,                    /* RESET  */
  YYSYMBOL_REVERSE = 353,                  /* REVERSE  */
  YYSYMBOL_REWIND = 354,                   /* REWIND  */
  YYSYMBOL_RIGHT = 355,                    /* RIGHT  */
  YYSYMBOL_RISER = 356,                    /* RISER  */
  YYSYMBOL_ROT = 357,                      /* ROT  */
  YYSYMBOL_RUN = 358,                      /* RUN  */
  YYSYMBOL_SALINITY = 359,                 /* SALINITY  */
  YYSYMBOL_SAMPLE = 360,                   /* SAMPLE  */
  YYSYMBOL_SAVE = 361,                     /* SAVE  */
  YYSYMBOL_SCALAR = 362,                   /* SCALAR  */
  YYSYMBOL_SCALE = 363,                    /* SCALE  */
  YYSYMBOL_SCIENTIFIC = 364,               /* SCIENTIFIC  */
  YYSYMBOL_SECONDS = 365,                  /* SECONDS  */
  YYSYMBOL_SET = 366,                      /* SET  */
  YYSYMBOL_SETS = 367,                     /* SETS  */
  YYSYMBOL_SHOW = 368,                     /* SHOW  */
  YYSYMBOL_SHRINK = 369,                   /* SHRINK  */
  YYSYMBOL_SIGMA = 370,                    /* SIGMA  */
  YYSYMBOL_SIGN = 371,                     /* SIGN  */
  YYSYMBOL_SIZE = 372,                     /* SIZE  */
  YYSYMBOL_SKIP = 373,                     /* SKIP  */
  YYSYMBOL_SLAB = 374,                     /* SLAB  */
  YYSYMBOL_SLEEP = 375,                    /* SLEEP  */
  YYSYMBOL_SLICE = 376,                    /* SLICE  */
  YYSYMBOL_SOURCE = 377,                   /* SOURCE  */
  YYSYMBOL_SPEC = 378,                     /* SPEC  */
  YYSYMBOL_SPECIFIED = 379,                /* SPECIFIED  */
  YYSYMBOL_SPECTRUM = 380,                 /* SPECTRUM  */
  YYSYMBOL_SPLITS = 381,                   /* SPLITS  */
  YYSYMBOL_SQUARE = 382,                   /* SQUARE  */
  YYSYMBOL_STACK = 383,                    /* STACK  */
  YYSYMBOL_STACKEDBAR = 384,               /* STACKEDBAR  */
  YYSYMBOL_STACKEDHBAR = 385,              /* STACKEDHBAR  */
  YYSYMBOL_STACKEDLINE = 386,              /* STACKEDLINE  */
  YYSYMBOL_STAGGER = 387,                  /* STAGGER  */
  YYSYMBOL_STAR = 388,                     /* STAR  */
  YYSYMBOL_START = 389,                    /* START  */
  YYSYMBOL_STARTSTEP = 390,                /* STARTSTEP  */
  YYSYMBOL_STARTTYPE = 391,                /* STARTTYPE  */
  YYSYMBOL_STATION = 392,                  /* STATION  */
  YYSYMBOL_STATUS = 393,                   /* STATUS  */
  YYSYMBOL_STEP = 394,                     /* STEP  */
  YYSYMBOL_STOP = 395,                     /* STOP  */
  YYSYMBOL_STREAMLINES = 396,              /* STREAMLINES  */
  YYSYMBOL_STRING = 397,                   /* STRING  */
  YYSYMBOL_STRINGS = 398,                  /* STRINGS  */
  YYSYMBOL_SUBTITLE = 399,                 /* SUBTITLE  */
  YYSYMBOL_SURFACE = 400,                  /* SURFACE  */
  YYSYMBOL_SWAPBUFFER = 401,               /* SWAPBUFFER  */
  YYSYMBOL_SYMBOL = 402,                   /* SYMBOL  */
  YYSYMBOL_SYSTEM = 403,                   /* SYSTEM  */
  YYSYMBOL_TEANL = 404,                    /* TEANL  */
  YYSYMBOL_TEXT = 405,                     /* TEXT  */
  YYSYMBOL_TICK = 406,                     /* TICK  */
  YYSYMBOL_TICKLABEL = 407,                /* TICKLABEL  */
  YYSYMBOL_TICKMARKS = 408,                /* TICKMARKS  */
  YYSYMBOL_TICKP = 409,                    /* TICKP  */
  YYSYMBOL_TIDALCLOCK = 410,               /* TIDALCLOCK  */
  YYSYMBOL_TIDESTATION = 411,              /* TIDESTATION  */
  YYSYMBOL_TIME = 412,                     /* TIME  */
  YYSYMBOL_TIMEINFO = 413,                 /* TIMEINFO  */
  YYSYMBOL_TIMELINE = 414,                 /* TIMELINE  */
  YYSYMBOL_TITLE = 415,                    /* TITLE  */
  YYSYMBOL_TO = 416,                       /* TO  */
  YYSYMBOL_TOP = 417,                      /* TOP  */
  YYSYMBOL_TOTAL = 418,                    /* TOTAL  */
  YYSYMBOL_TRACK = 419,                    /* TRACK  */
  YYSYMBOL_TRANSECT = 420,                 /* TRANSECT  */
  YYSYMBOL_TRIANGLE1 = 421,                /* TRIANGLE1  */
  YYSYMBOL_TRIANGLE2 = 422,                /* TRIANGLE2  */
  YYSYMBOL_TRIANGLE3 = 423,                /* TRIANGLE3  */
  YYSYMBOL_TRIANGLE4 = 424,                /* TRIANGLE4  */
  YYSYMBOL_TRUEP = 425,                    /* TRUEP  */
  YYSYMBOL_TYPE = 426,                     /* TYPE  */
  YYSYMBOL_UNITS = 427,                    /* UNITS  */
  YYSYMBOL_UP = 428,                       /* UP  */
  YYSYMBOL_VALUE = 429,                    /* VALUE  */
  YYSYMBOL_VECTOR = 430,                   /* VECTOR  */
  YYSYMBOL_VEL = 431,                      /* VEL  */
  YYSYMBOL_VELMARKER = 432,                /* VELMARKER  */
  YYSYMBOL_VELOCITY = 433,                 /* VELOCITY  */
  YYSYMBOL_VERTICAL = 434,                 /* VERTICAL  */
  YYSYMBOL_VGAP = 435,                     /* VGAP  */
  YYSYMBOL_VIEW = 436,                     /* VIEW  */
  YYSYMBOL_VSCALE = 437,                   /* VSCALE  */
  YYSYMBOL_VX1 = 438,                      /* VX1  */
  YYSYMBOL_VX2 = 439,                      /* VX2  */
  YYSYMBOL_VY1 = 440,                      /* VY1  */
  YYSYMBOL_VY2 = 441,                      /* VY2  */
  YYSYMBOL_WEEKS = 442,                    /* WEEKS  */
  YYSYMBOL_WET = 443,                      /* WET  */
  YYSYMBOL_WETDRY = 444,                   /* WETDRY  */
  YYSYMBOL_WIDTH = 445,                    /* WIDTH  */
  YYSYMBOL_WIND = 446,                     /* WIND  */
  YYSYMBOL_WITH = 447,                     /* WITH  */
  YYSYMBOL_WORLD = 448,                    /* WORLD  */
  YYSYMBOL_WRAP = 449,                     /* WRAP  */
  YYSYMBOL_WRITE = 450,                    /* WRITE  */
  YYSYMBOL_WSCALE = 451,                   /* WSCALE  */
  YYSYMBOL_WX1 = 452,                      /* WX1  */
  YYSYMBOL_WX2 = 453,                      /* WX2  */
  YYSYMBOL_WY1 = 454,                      /* WY1  */
  YYSYMBOL_WY2 = 455,                      /* WY2  */
  YYSYMBOL_X0 = 456,                       /* X0  */
  YYSYMBOL_X1 = 457,                       /* X1  */
  YYSYMBOL_X2 = 458,                       /* X2  */
  YYSYMBOL_X3 = 459,                       /* X3  */
  YYSYMBOL_X4 = 460,                       /* X4  */
  YYSYMBOL_X5 = 461,                       /* X5  */
  YYSYMBOL_XAXES = 462,                    /* XAXES  */
  YYSYMBOL_XAXIS = 463,                    /* XAXIS  */
  YYSYMBOL_XCOR = 464,                     /* XCOR  */
  YYSYMBOL_XMAX = 465,                     /* XMAX  */
  YYSYMBOL_XMIN = 466,                     /* XMIN  */
  YYSYMBOL_XY = 467,                       /* XY  */
  YYSYMBOL_XYARC = 468,                    /* XYARC  */
  YYSYMBOL_XYBOX = 469,                    /* XYBOX  */
  YYSYMBOL_XYDX = 470,                     /* XYDX  */
  YYSYMBOL_XYDXDX = 471,                   /* XYDXDX  */
  YYSYMBOL_XYDXDY = 472,                   /* XYDXDY  */
  YYSYMBOL_XYDY = 473,                     /* XYDY  */
  YYSYMBOL_XYDYDY = 474,                   /* XYDYDY  */
  YYSYMBOL_XYFIXED = 475,                  /* XYFIXED  */
  YYSYMBOL_XYHILO = 476,                   /* XYHILO  */
  YYSYMBOL_XYRT = 477,                     /* XYRT  */
  YYSYMBOL_XYSEG = 478,                    /* XYSEG  */
  YYSYMBOL_XYSTRING = 479,                 /* XYSTRING  */
  YYSYMBOL_XYUV = 480,                     /* XYUV  */
  YYSYMBOL_XYX2Y2 = 481,                   /* XYX2Y2  */
  YYSYMBOL_XYXX = 482,                     /* XYXX  */
  YYSYMBOL_XYYY = 483,                     /* XYYY  */
  YYSYMBOL_XYZ = 484,                      /* XYZ  */
  YYSYMBOL_XYZW = 485,                     /* XYZW  */
  YYSYMBOL_Y0 = 486,                       /* Y0  */
  YYSYMBOL_Y1 = 487,                       /* Y1  */
  YYSYMBOL_Y2 = 488,                       /* Y2  */
  YYSYMBOL_Y3 = 489,                       /* Y3  */
  YYSYMBOL_Y4 = 490,                       /* Y4  */
  YYSYMBOL_Y5 = 491,                       /* Y5  */
  YYSYMBOL_YAXES = 492,                    /* YAXES  */
  YYSYMBOL_YAXIS = 493,                    /* YAXIS  */
  YYSYMBOL_YEARS = 494,                    /* YEARS  */
  YYSYMBOL_YES = 495,                      /* YES  */
  YYSYMBOL_YMAX = 496,                     /* YMAX  */
  YYSYMBOL_YMIN = 497,                     /* YMIN  */
  YYSYMBOL_ZEROXAXIS = 498,                /* ZEROXAXIS  */
  YYSYMBOL_ZEROYAXIS = 499,                /* ZEROYAXIS  */
  YYSYMBOL_ZOOM = 500,                     /* ZOOM  */
  YYSYMBOL_ZOOMBOX = 501,                  /* ZOOMBOX  */
  YYSYMBOL_502_ = 502,                     /* '='  */
  YYSYMBOL_OR = 503,                       /* OR  */
  YYSYMBOL_AND = 504,                      /* AND  */
  YYSYMBOL_GT = 505,                       /* GT  */
  YYSYMBOL_LT = 506,                       /* LT  */
  YYSYMBOL_LE = 507,                       /* LE  */
  YYSYMBOL_GE = 508,                       /* GE  */
  YYSYMBOL_EQ = 509,                       /* EQ  */
  YYSYMBOL_NE = 510,                       /* NE  */
  YYSYMBOL_511_ = 511,                     /* '+'  */
  YYSYMBOL_512_ = 512,                     /* '-'  */
  YYSYMBOL_513_ = 513,                     /* '*'  */
  YYSYMBOL_514_ = 514,                     /* '/'  */
  YYSYMBOL_515_ = 515,                     /* '%'  */
  YYSYMBOL_516_ = 516,                     /* '^'  */
  YYSYMBOL_UMINUS = 517,                   /* UMINUS  */
  YYSYMBOL_NOT = 518,                      /* NOT  */
  YYSYMBOL_519_n_ = 519,                   /* '\n'  */
  YYSYMBOL_520_ = 520,                     /* ','  */
  YYSYMBOL_521_ = 521,                     /* '['  */
  YYSYMBOL_522_ = 522,                     /* ']'  */
  YYSYMBOL_523_ = 523,                     /* '('  */
  YYSYMBOL_524_ = 524,                     /* ')'  */
  YYSYMBOL_525_ = 525,                     /* '.'  */
  YYSYMBOL_YYACCEPT = 526,                 /* $accept  */
  YYSYMBOL_list = 527,                     /* list  */
  YYSYMBOL_annotation = 528,               /* annotation  */
  YYSYMBOL_animation = 529,                /* animation  */
  YYSYMBOL_misc = 530,                     /* misc  */
  YYSYMBOL_531_1 = 531,                    /* $@1  */
  YYSYMBOL_532_2 = 532,                    /* $@2  */
  YYSYMBOL_models = 533,                   /* models  */
  YYSYMBOL_534_3 = 534,                    /* $@3  */
  YYSYMBOL_535_4 = 535,                    /* $@4  */
  YYSYMBOL_536_5 = 536,                    /* $@5  */
  YYSYMBOL_537_6 = 537,                    /* $@6  */
  YYSYMBOL_538_7 = 538,                    /* $@7  */
  YYSYMBOL_flowprops = 539,                /* flowprops  */
  YYSYMBOL_540_8 = 540,                    /* $@8  */
  YYSYMBOL_541_9 = 541,                    /* $@9  */
  YYSYMBOL_542_10 = 542,                   /* $@10  */
  YYSYMBOL_543_11 = 543,                   /* $@11  */
  YYSYMBOL_544_12 = 544,                   /* $@12  */
  YYSYMBOL_elevmarker = 545,               /* elevmarker  */
  YYSYMBOL_546_13 = 546,                   /* $@13  */
  YYSYMBOL_grid = 547,                     /* grid  */
  YYSYMBOL_548_14 = 548,                   /* $@14  */
  YYSYMBOL_549_15 = 549,                   /* $@15  */
  YYSYMBOL_550_16 = 550,                   /* $@16  */
  YYSYMBOL_drogues = 551,                  /* drogues  */
  YYSYMBOL_552_17 = 552,                   /* $@17  */
  YYSYMBOL_isolines = 553,                 /* isolines  */
  YYSYMBOL_554_18 = 554,                   /* $@18  */
  YYSYMBOL_555_19 = 555,                   /* $@19  */
  YYSYMBOL_histboxes = 556,                /* histboxes  */
  YYSYMBOL_557_20 = 557,                   /* $@20  */
  YYSYMBOL_zoomboxes = 558,                /* zoomboxes  */
  YYSYMBOL_559_21 = 559,                   /* $@21  */
  YYSYMBOL_wscale = 560,                   /* wscale  */
  YYSYMBOL_vscale = 561,                   /* vscale  */
  YYSYMBOL_mapscale = 562,                 /* mapscale  */
  YYSYMBOL_tidalclock = 563,               /* tidalclock  */
  YYSYMBOL_timeinfo = 564,                 /* timeinfo  */
  YYSYMBOL_timeline = 565,                 /* timeline  */
  YYSYMBOL_slice = 566,                    /* slice  */
  YYSYMBOL_567_22 = 567,                   /* $@22  */
  YYSYMBOL_props = 568,                    /* props  */
  YYSYMBOL_graph = 569,                    /* graph  */
  YYSYMBOL_setaxis = 570,                  /* setaxis  */
  YYSYMBOL_axis = 571,                     /* axis  */
  YYSYMBOL_allaxes = 572,                  /* allaxes  */
  YYSYMBOL_axesprops = 573,                /* axesprops  */
  YYSYMBOL_axisfeature = 574,              /* axisfeature  */
  YYSYMBOL_tickattr = 575,                 /* tickattr  */
  YYSYMBOL_ticklabeldesc = 576,            /* ticklabeldesc  */
  YYSYMBOL_ticklabelattr = 577,            /* ticklabelattr  */
  YYSYMBOL_axislabeldesc = 578,            /* axislabeldesc  */
  YYSYMBOL_axisbardesc = 579,              /* axisbardesc  */
  YYSYMBOL_sourcetype = 580,               /* sourcetype  */
  YYSYMBOL_justchoice = 581,               /* justchoice  */
  YYSYMBOL_graphtype = 582,                /* graphtype  */
  YYSYMBOL_inoutchoice = 583,              /* inoutchoice  */
  YYSYMBOL_signchoice = 584,               /* signchoice  */
  YYSYMBOL_formatchoice = 585,             /* formatchoice  */
  YYSYMBOL_horv = 586,                     /* horv  */
  YYSYMBOL_flowonoff = 587,                /* flowonoff  */
  YYSYMBOL_onoff = 588,                    /* onoff  */
  YYSYMBOL_worldview = 589,                /* worldview  */
  YYSYMBOL_torf = 590,                     /* torf  */
  YYSYMBOL_filltype = 591,                 /* filltype  */
  YYSYMBOL_opchoice = 592,                 /* opchoice  */
  YYSYMBOL_units = 593,                    /* units  */
  YYSYMBOL_asgn = 594,                     /* asgn  */
  YYSYMBOL_vasgn = 595,                    /* vasgn  */
  YYSYMBOL_vexpr = 596,                    /* vexpr  */
  YYSYMBOL_expr = 597                      /* expr  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

#if defined __GNUC__ && ! defined __ICC && 407 <= __GNUC__ * 100 + __GNUC_MINOR__
/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                            \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  522
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   6977

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  526
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  72
/* YYNRULES -- Number of rules.  */
#define YYNRULES  797
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  1798

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   766


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_int16 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
     519,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,   515,     2,     2,
     523,   524,   513,   511,   520,   512,   525,   514,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   502,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   521,     2,   522,   516,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160,   161,   162,   163,   164,
     165,   166,   167,   168,   169,   170,   171,   172,   173,   174,
     175,   176,   177,   178,   179,   180,   181,   182,   183,   184,
     185,   186,   187,   188,   189,   190,   191,   192,   193,   194,
     195,   196,   197,   198,   199,   200,   201,   202,   203,   204,
     205,   206,   207,   208,   209,   210,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,   246,   247,   248,   249,   250,   251,   252,   253,   254,
     255,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   503,   504,   505,
     506,   507,   508,   509,   510,   517,   518
};

#if YYDEBUG
  /* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   668,   668,   669,   670,   671,   672,   673,   674,   675,
     676,   677,   678,   679,   680,   681,   682,   683,   684,   685,
     686,   687,   688,   689,   690,   691,   692,   695,   696,   697,
     698,   699,   700,   701,   710,   711,   712,   713,   714,   715,
     716,   717,   729,   730,   731,   732,   733,   734,   741,   742,
     743,   744,   745,   746,   747,   748,   760,   761,   762,   763,
     764,   765,   770,   771,   772,   773,   774,   775,   776,   777,
     778,   779,   780,   781,   782,   802,   803,   804,   805,   806,
     807,   808,   809,   810,   811,   812,   813,   817,   818,   824,
     830,   831,   832,   833,   834,   835,   836,   837,   838,   839,
     840,   841,   842,   843,   844,   854,   863,   864,   865,   869,
     871,   873,   875,   875,   876,   876,   880,   883,   884,   885,
     889,   890,   899,   904,   905,   906,   910,   914,   923,   932,
     941,   951,   960,   969,   978,   987,   991,   992,   992,   993,
     994,   995,   996,   997,   998,   999,  1000,  1001,  1002,  1003,
    1004,  1010,  1011,  1012,  1013,  1014,  1019,  1028,  1036,  1044,
    1052,  1060,  1068,  1076,  1082,  1086,  1090,  1094,  1098,  1102,
    1106,  1113,  1120,  1128,  1136,  1136,  1137,  1137,  1138,  1138,
    1139,  1139,  1140,  1145,  1146,  1147,  1148,  1149,  1150,  1151,
    1152,  1159,  1164,  1169,  1177,  1178,  1179,  1180,  1181,  1182,
    1183,  1184,  1185,  1186,  1187,  1188,  1189,  1189,  1190,  1190,
    1191,  1191,  1192,  1192,  1193,  1193,  1194,  1195,  1196,  1197,
    1198,  1199,  1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,
    1208,  1212,  1213,  1214,  1215,  1216,  1216,  1217,  1218,  1219,
    1224,  1236,  1237,  1238,  1239,  1240,  1241,  1242,  1243,  1243,
    1244,  1244,  1245,  1246,  1247,  1248,  1249,  1249,  1250,  1255,
    1260,  1269,  1270,  1271,  1271,  1272,  1283,  1294,  1295,  1296,
    1297,  1298,  1299,  1300,  1301,  1302,  1303,  1311,  1319,  1327,
    1336,  1339,  1349,  1350,  1351,  1352,  1356,  1357,  1358,  1359,
    1359,  1360,  1360,  1361,  1365,  1372,  1373,  1373,  1374,  1379,
    1380,  1381,  1382,  1383,  1384,  1385,  1386,  1387,  1388,  1389,
    1390,  1391,  1392,  1399,  1404,  1409,  1417,  1418,  1418,  1419,
    1424,  1425,  1426,  1427,  1428,  1429,  1430,  1431,  1432,  1433,
    1434,  1441,  1446,  1451,  1459,  1460,  1461,  1462,  1463,  1464,
    1469,  1480,  1485,  1486,  1487,  1488,  1489,  1490,  1495,  1506,
    1510,  1511,  1512,  1513,  1514,  1515,  1520,  1531,  1535,  1536,
    1537,  1538,  1539,  1540,  1545,  1549,  1550,  1555,  1560,  1561,
    1562,  1563,  1564,  1565,  1566,  1567,  1571,  1572,  1573,  1574,
    1575,  1576,  1577,  1578,  1579,  1580,  1581,  1582,  1587,  1591,
    1592,  1592,  1593,  1598,  1599,  1600,  1601,  1602,  1603,  1604,
    1605,  1606,  1613,  1618,  1623,  1631,  1632,  1633,  1634,  1635,
    1636,  1637,  1638,  1639,  1640,  1641,  1642,  1643,  1644,  1645,
    1646,  1650,  1651,  1652,  1653,  1654,  1659,  1669,  1670,  1671,
    1672,  1673,  1674,  1675,  1676,  1677,  1681,  1685,  1686,  1693,
    1694,  1695,  1696,  1697,  1704,  1705,  1706,  1707,  1708,  1711,
    1714,  1715,  1718,  1722,  1725,  1728,  1729,  1733,  1737,  1738,
    1739,  1740,  1741,  1742,  1743,  1744,  1745,  1746,  1747,  1748,
    1749,  1750,  1751,  1752,  1757,  1762,  1767,  1771,  1772,  1773,
    1774,  1775,  1779,  1780,  1781,  1782,  1783,  1784,  1788,  1789,
    1790,  1794,  1795,  1796,  1797,  1798,  1799,  1800,  1805,  1806,
    1807,  1808,  1809,  1818,  1823,  1824,  1825,  1826,  1827,  1828,
    1829,  1830,  1831,  1832,  1833,  1834,  1835,  1836,  1837,  1838,
    1839,  1840,  1841,  1842,  1843,  1844,  1845,  1846,  1847,  1848,
    1849,  1850,  1851,  1852,  1856,  1857,  1861,  1862,  1863,  1864,
    1865,  1866,  1867,  1868,  1869,  1870,  1871,  1872,  1873,  1874,
    1875,  1876,  1877,  1878,  1879,  1880,  1881,  1882,  1883,  1884,
    1885,  1886,  1887,  1888,  1889,  1890,  1894,  1895,  1896,  1897,
    1898,  1899,  1903,  1904,  1905,  1906,  1907,  1911,  1912,  1913,
    1914,  1926,  1927,  1931,  1932,  1933,  1942,  1943,  1944,  1945,
    1946,  1950,  1951,  1952,  1956,  1957,  1958,  1971,  1972,  1973,
    1974,  1975,  1976,  1977,  1978,  1979,  1980,  1981,  1982,  1983,
    1984,  1985,  1986,  1987,  1988,  1989,  1990,  1991,  1992,  1993,
    1994,  1995,  1996,  1997,  2000,  2001,  2005,  2006,  2007,  2008,
    2012,  2013,  2016,  2017,  2028,  2029,  2032,  2033,  2034,  2037,
    2038,  2039,  2040,  2041,  2042,  2043,  2047,  2048,  2049,  2050,
    2064,  2078,  2086,  2097,  2106,  2115,  2124,  2133,  2142,  2151,
    2160,  2169,  2178,  2187,  2196,  2205,  2214,  2223,  2236,  2251,
    2266,  2279,  2288,  2297,  2306,  2315,  2324,  2333,  2342,  2351,
    2360,  2369,  2378,  2387,  2396,  2405,  2414,  2423,  2432,  2441,
    2450,  2459,  2468,  2477,  2486,  2495,  2504,  2513,  2522,  2531,
    2540,  2549,  2558,  2567,  2576,  2585,  2594,  2604,  2613,  2622,
    2631,  2640,  2649,  2658,  2667,  2676,  2685,  2694,  2703,  2712,
    2721,  2730,  2739,  2748,  2757,  2767,  2768,  2771,  2774,  2777,
    2780,  2783,  2792,  2795,  2798,  2801,  2804,  2807,  2810,  2813,
    2816,  2819,  2822,  2825,  2828,  2831,  2834,  2837,  2840,  2843,
    2846,  2849,  2852,  2855,  2858,  2861,  2864,  2867,  2870,  2873,
    2876,  2879,  2882,  2885,  2888,  2891,  2894,  2897,  2900,  2903,
    2906,  2909,  2912,  2916,  2919,  2922,  2925,  2928,  2932,  2935,
    2938,  2941,  2944,  2947,  2950,  2953,  2957,  2964,  2967,  2970,
    2973,  2976,  2979,  2982,  2985,  2988,  2991,  2994
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "VAR", "X", "Y",
  "CHRSTR", "FITPARM", "NUMBER", "ABS", "ACOS", "ASIN", "ATAN", "ATAN2",
  "CEIL", "COS", "DEG", "DX", "DY", "ERF", "ERFC", "EXP", "FLOOR", "HYPOT",
  "INDEX", "INT", "IRAND", "LGAMMA", "LN", "LOG", "LOGISTIC", "MAXP",
  "MINP", "MINMAX", "MOD", "NORM", "NORMP", "PI", "RAD", "RAND", "SETNO",
  "SIN", "SQR", "SQRT", "TAN", "INUM", "ABORT", "ABOVE", "ABSOLUTE",
  "ACTIVATE", "ACTIVE", "ADCIRC", "ADCIRC3DFLOW", "ALL", "ALT",
  "ALTERNATE", "ALTXAXIS", "ALTYAXIS", "AMP", "ANGLE", "ANNOTATE",
  "APPEND", "AREA", "ARROW", "ASCEND", "AT", "ATTACH", "AUTO", "AUTOSCALE",
  "AUTOTICKS", "AVERAGE", "AVG", "AXES", "AXIS", "BACKBUFFER",
  "BACKGROUND", "BAR", "BATCH", "BATH", "BATHYMETRY", "COURANT", "BELOW",
  "BIN", "BINARY", "BOTH", "BOTTOM", "BOUNDARY", "BOX", "CELLS", "CENTER",
  "CH3D", "CHAR", "CHDIR", "CIRCLE", "CLEAR", "CLICK", "CLOCK", "CLOSE",
  "CM", "CMAP", "COLOR", "COLORMAP", "COMMENT", "CONC", "CONCENTRATION",
  "CONCENTRATIONS", "COPY", "CROSS", "CYCLE", "DAYMONTH", "DAYOFWEEKL",
  "DAYOFWEEKS", "DAYOFYEAR", "DAYS", "DDMMYY", "DDMONTHSYYHHMMSS",
  "DECIMAL", "DEF", "DEFAULT", "DEGREESLAT", "DEGREESLON", "DEGREESMMLAT",
  "DEGREESMMLON", "DEGREESMMSSLAT", "DEGREESMMSSLON", "DELAYP", "DELETE",
  "DEPTH", "DEPTHS", "DESCEND", "DEVICE", "DEVXY", "DFT", "DT", "DIAMOND",
  "DIFFERENCE", "DISK", "DISPLAY", "DOT", "DOUBLEBUFFER", "DOWN", "DRAW2",
  "DROGUE", "DROGUES", "DRY", "DXDX", "DXP", "DYDY", "DYP", "ECHO", "EDIT",
  "ELA", "ELCIRC", "ELEMENT", "ELEMENTS", "ELEV", "ELEVATION",
  "ELEVATIONS", "ELEVMARKER", "ELLIPSE", "ELLIPSES", "ELLIPSEZ", "ELSE",
  "END", "ERRORBAR", "EXIT", "EXPAND", "EXPONENTIAL", "FACTOR", "FALSEP",
  "FAST", "FEET", "FFT", "FILEP", "FILL", "FIND", "FIXEDPOINT", "FLOW",
  "FLUSH", "FLUX", "FOCUS", "FOLLOWS", "FONTP", "FOREGROUND", "FORMAT",
  "FORT14", "FORT63", "FORT64", "FORWARD", "FRAMEP", "FREQ", "FRONTBUFFER",
  "GENERAL", "GETP", "GOTO", "GRAPH", "GRAPHNO", "GRAPHS", "GRAPHTYPE",
  "GRID", "HARDCOPY", "HBAR", "HELP", "HGAP", "HIDDEN", "HISTBOX", "HISTO",
  "HISTORY", "HMS", "HORIZONTAL", "HOURS", "HPGLL", "HPGLP", "IF",
  "IGNORE", "IHL", "IMAGE", "IMAGES", "IN", "INCLUDE", "INFO", "INIT",
  "INITGRAPHICS", "INOUT", "INTEGRATE", "INTERP", "INUNDATION", "INVDFT",
  "INVFFT", "ISOLINE", "ISOLINES", "JUST", "KILL", "KM", "LABEL", "LAYOUT",
  "LEAVE", "LEAVEGRAPHICS", "LEFT", "LEGEND", "LENGTH", "LEVEL", "LEVELS",
  "LIMITS", "LINE", "LINES", "LINESTYLE", "LINETO", "LINEW", "LINEWIDTH",
  "LINK", "LOAD", "LOC", "LOCATE", "LOCATOR", "LOCTYPE", "LOGX", "LOGXY",
  "LOGY", "M", "MAG", "MAGNITUDE", "MAJOR", "MAPSCALE", "MARKER",
  "MARKERS", "MAXLEVELS", "METHOD", "MIFL", "MIFP", "MILES", "MINOR",
  "MINUTES", "MISSINGP", "MM", "MMDD", "MMDDHMS", "MMDDYY", "MMDDYYHMS",
  "MMSSLAT", "MMSSLON", "MMYY", "MONTHDAY", "MONTHL", "MONTHS", "MOVE",
  "MOVE2", "MOVETO", "NEGATE", "NO", "NODE", "NODES", "NONE", "NORMAL",
  "NORTH", "NXY", "OFF", "OFFSETX", "OFFSETY", "ON", "OP", "OPEN",
  "ORIENT", "OUT", "PAGE", "PARA", "PARALLEL", "PARAMETERS", "PARAMS",
  "PARMS", "PATTERN", "PER", "PERIMETER", "PERP", "PERPENDICULAR", "PHASE",
  "PIE", "PIPE", "PLACE", "PLAN", "PLUS", "POINT", "POLAR", "POLY",
  "POLYI", "POLYO", "POP", "POWER", "PREC", "PREFIX", "PREPEND", "PRINT",
  "PROFILE", "PROP", "PS", "PSCOLORL", "PSCOLORP", "PSMONOL", "PSMONOP",
  "PUSH", "PUTP", "QUIT", "READ", "READBIN", "REDRAW", "REGION", "REGIONS",
  "REGNUM", "REGRESS", "REMOVE", "RENDER", "REPORT", "RESET", "REVERSE",
  "REWIND", "RIGHT", "RISER", "ROT", "RUN", "SALINITY", "SAMPLE", "SAVE",
  "SCALAR", "SCALE", "SCIENTIFIC", "SECONDS", "SET", "SETS", "SHOW",
  "SHRINK", "SIGMA", "SIGN", "SIZE", "SKIP", "SLAB", "SLEEP", "SLICE",
  "SOURCE", "SPEC", "SPECIFIED", "SPECTRUM", "SPLITS", "SQUARE", "STACK",
  "STACKEDBAR", "STACKEDHBAR", "STACKEDLINE", "STAGGER", "STAR", "START",
  "STARTSTEP", "STARTTYPE", "STATION", "STATUS", "STEP", "STOP",
  "STREAMLINES", "STRING", "STRINGS", "SUBTITLE", "SURFACE", "SWAPBUFFER",
  "SYMBOL", "SYSTEM", "TEANL", "TEXT", "TICK", "TICKLABEL", "TICKMARKS",
  "TICKP", "TIDALCLOCK", "TIDESTATION", "TIME", "TIMEINFO", "TIMELINE",
  "TITLE", "TO", "TOP", "TOTAL", "TRACK", "TRANSECT", "TRIANGLE1",
  "TRIANGLE2", "TRIANGLE3", "TRIANGLE4", "TRUEP", "TYPE", "UNITS", "UP",
  "VALUE", "VECTOR", "VEL", "VELMARKER", "VELOCITY", "VERTICAL", "VGAP",
  "VIEW", "VSCALE", "VX1", "VX2", "VY1", "VY2", "WEEKS", "WET", "WETDRY",
  "WIDTH", "WIND", "WITH", "WORLD", "WRAP", "WRITE", "WSCALE", "WX1",
  "WX2", "WY1", "WY2", "X0", "X1", "X2", "X3", "X4", "X5", "XAXES",
  "XAXIS", "XCOR", "XMAX", "XMIN", "XY", "XYARC", "XYBOX", "XYDX",
  "XYDXDX", "XYDXDY", "XYDY", "XYDYDY", "XYFIXED", "XYHILO", "XYRT",
  "XYSEG", "XYSTRING", "XYUV", "XYX2Y2", "XYXX", "XYYY", "XYZ", "XYZW",
  "Y0", "Y1", "Y2", "Y3", "Y4", "Y5", "YAXES", "YAXIS", "YEARS", "YES",
  "YMAX", "YMIN", "ZEROXAXIS", "ZEROYAXIS", "ZOOM", "ZOOMBOX", "'='", "OR",
  "AND", "GT", "LT", "LE", "GE", "EQ", "NE", "'+'", "'-'", "'*'", "'/'",
  "'%'", "'^'", "UMINUS", "NOT", "'\\n'", "','", "'['", "']'", "'('",
  "')'", "'.'", "$accept", "list", "annotation", "animation", "misc",
  "$@1", "$@2", "models", "$@3", "$@4", "$@5", "$@6", "$@7", "flowprops",
  "$@8", "$@9", "$@10", "$@11", "$@12", "elevmarker", "$@13", "grid",
  "$@14", "$@15", "$@16", "drogues", "$@17", "isolines", "$@18", "$@19",
  "histboxes", "$@20", "zoomboxes", "$@21", "wscale", "vscale", "mapscale",
  "tidalclock", "timeinfo", "timeline", "slice", "$@22", "props", "graph",
  "setaxis", "axis", "allaxes", "axesprops", "axisfeature", "tickattr",
  "ticklabeldesc", "ticklabelattr", "axislabeldesc", "axisbardesc",
  "sourcetype", "justchoice", "graphtype", "inoutchoice", "signchoice",
  "formatchoice", "horv", "flowonoff", "onoff", "worldview", "torf",
  "filltype", "opchoice", "units", "asgn", "vasgn", "vexpr", "expr", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#ifdef YYPRINT
/* YYTOKNUM[NUM] -- (External) token number corresponding to the
   (internal) symbol number NUM (which must be that of a token).  */
static const yytype_int16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283,   284,
     285,   286,   287,   288,   289,   290,   291,   292,   293,   294,
     295,   296,   297,   298,   299,   300,   301,   302,   303,   304,
     305,   306,   307,   308,   309,   310,   311,   312,   313,   314,
     315,   316,   317,   318,   319,   320,   321,   322,   323,   324,
     325,   326,   327,   328,   329,   330,   331,   332,   333,   334,
     335,   336,   337,   338,   339,   340,   341,   342,   343,   344,
     345,   346,   347,   348,   349,   350,   351,   352,   353,   354,
     355,   356,   357,   358,   359,   360,   361,   362,   363,   364,
     365,   366,   367,   368,   369,   370,   371,   372,   373,   374,
     375,   376,   377,   378,   379,   380,   381,   382,   383,   384,
     385,   386,   387,   388,   389,   390,   391,   392,   393,   394,
     395,   396,   397,   398,   399,   400,   401,   402,   403,   404,
     405,   406,   407,   408,   409,   410,   411,   412,   413,   414,
     415,   416,   417,   418,   419,   420,   421,   422,   423,   424,
     425,   426,   427,   428,   429,   430,   431,   432,   433,   434,
     435,   436,   437,   438,   439,   440,   441,   442,   443,   444,
     445,   446,   447,   448,   449,   450,   451,   452,   453,   454,
     455,   456,   457,   458,   459,   460,   461,   462,   463,   464,
     465,   466,   467,   468,   469,   470,   471,   472,   473,   474,
     475,   476,   477,   478,   479,   480,   481,   482,   483,   484,
     485,   486,   487,   488,   489,   490,   491,   492,   493,   494,
     495,   496,   497,   498,   499,   500,   501,   502,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   517,   518,   519,   520,   521,   522,   523,   524,
     525,   526,   527,   528,   529,   530,   531,   532,   533,   534,
     535,   536,   537,   538,   539,   540,   541,   542,   543,   544,
     545,   546,   547,   548,   549,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,   562,   563,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   578,   579,   580,   581,   582,   583,   584,
     585,   586,   587,   588,   589,   590,   591,   592,   593,   594,
     595,   596,   597,   598,   599,   600,   601,   602,   603,   604,
     605,   606,   607,   608,   609,   610,   611,   612,   613,   614,
     615,   616,   617,   618,   619,   620,   621,   622,   623,   624,
     625,   626,   627,   628,   629,   630,   631,   632,   633,   634,
     635,   636,   637,   638,   639,   640,   641,   642,   643,   644,
     645,   646,   647,   648,   649,   650,   651,   652,   653,   654,
     655,   656,   657,   658,   659,   660,   661,   662,   663,   664,
     665,   666,   667,   668,   669,   670,   671,   672,   673,   674,
     675,   676,   677,   678,   679,   680,   681,   682,   683,   684,
     685,   686,   687,   688,   689,   690,   691,   692,   693,   694,
     695,   696,   697,   698,   699,   700,   701,   702,   703,   704,
     705,   706,   707,   708,   709,   710,   711,   712,   713,   714,
     715,   716,   717,   718,   719,   720,   721,   722,   723,   724,
     725,   726,   727,   728,   729,   730,   731,   732,   733,   734,
     735,   736,   737,   738,   739,   740,   741,   742,   743,   744,
     745,   746,   747,   748,   749,   750,   751,   752,   753,   754,
     755,   756,    61,   757,   758,   759,   760,   761,   762,   763,
     764,    43,    45,    42,    47,    37,    94,   765,   766,    10,
      44,    91,    93,    40,    41,    46
};
#endif

#define YYPACT_NINF (-954)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-798)

#define yytable_value_is_error(Yyn) \
  ((Yyn) == YYTABLE_NINF)

  /* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
     STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    1797,  -452,  -423,  -954,  -954,  -448,  -443,  -431,  -427,  -415,
    -407,  -401,  -402,  -379,  -375,  -372,  -368,  -319,  -312,  -310,
    -371,  -308,  -304,  -287,  -261,    14,  -254,  -241,  -235,  -226,
    -219,  -190,  -225,  -208,  -191,  -179,  -181,  -180,  -167,  -165,
     687,  -954,  -954,  -954,   676,  -269,  2451,    94,   -34,   318,
     354,   360,   364,  -954,   239,   372,   728,    60,  -954,   285,
    -954,   177,   378,   -29,   -51,   741,   -43,  -125,   187,  3898,
      18,   114,  2904,   111,  -195,  -954,    32,   504,  -954,  -954,
    -954,  -954,    -4,  -954,  -954,  -954,  -954,    53,  -954,   266,
    -954,   420,   958,  -109,    -5,  -954,  -954,  2314,    11,   425,
     687,   165,  2357,   996,    34,  2068,   403,  -954,  -954,  -954,
    -954,   558,  3586,  -195,   485,  -954,  -954,  -954,  -954,   676,
    -954,   676,  -954,  -954,  -954,  4192,   655,  4234,  4234,  4234,
     433,   -85,   -80,   -74,   -73,   -70,   -68,   -67,   -65,   -64,
     -63,   -56,   -52,   -44,   -42,   -41,   -39,   -28,   -25,   -22,
     -40,  -954,   -18,   -14,  6430,  6413,  -954,  4234,  4192,  4234,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,
    4234,  4234,   464,  4234,  4234,  -954,  4234,  4234,  4234,  4234,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  -954,   466,   119,
     571,   456,   487,  -161,   491,  -954,   323,  -954,   181,   547,
     548,  -195,   555,   570,  -954,  -954,  -954,  -954,   573,   572,
      65,    59,    64,    82,    83,    85,    87,   102,  -954,  -954,
    -954,   103,   104,   106,   107,   109,  -954,   116,   117,   118,
     123,   125,   126,   127,   130,   131,   134,   136,  -954,  -954,
    -954,  -954,   138,   153,   167,   168,   629,  -954,    54,   169,
     644,   684,  -358,  4192,  4192,  4192,  -954,  4955,  -954,  -954,
    -954,  -954,  -954,  -954,   173,   175,   176,  -195,   365,  -954,
     689,  1007,   701,   716,  -954,  4192,  -195,  -195,  4192,  -358,
     703,  -195,   365,  -954,  -954,  -954,  -954,  -954,   612,   706,
    -195,   707,   710,   711,  -954,  -954,  -134,   463,  -151,  -195,
    -137,   736,  -954,   -40,  -954,  -954,   586,   591,   595,   596,
    -195,   597,  -195,   598,   365,   730,   731,   537,   646,   739,
    4192,  -358,   743,   744,   753,  4192,  4192,  4192,   365,  4192,
     758,   -46,   169,  4973,  -954,  4192,   339,  3015,  -143,   342,
    4192,   764,   365,  -954,  -954,    -1,   767,  -954,   169,   769,
     770,  -358,  -954,  4991,  -195,  -954,    84,   771,   772,  4192,
    -358,   -32,   774,   101,  -954,  -954,  -954,   775,  -954,  -954,
    -954,  6447,   630,   782,    30,   785,  -954,   270,   789,    84,
     376,  -954,   796,   797,  -173,   708,   798,  4192,  -358,   799,
    4192,  4192,  4192,   365,  -954,  -954,  -954,  -954,  -954,  4192,
     437,   807,   810,   812,   169,   817,   823,  -358,   826,    24,
    -954,  5009,  -954,   827,   829,   830,   831,  -954,  -954,   832,
     746,  4192,  -358,   429,  -954,  -954,   472,   842,   843,   846,
     847,   848,  -358,   849,   854,  -954,  5027,   853,   763,   858,
    4192,  -358,   861,   864,   865,   866,   868,   872,  -954,  -954,
    -954,   879,   880,   882,   883,   884,   889,   892,   893,  5045,
     894,   895,  4192,  -358,   898,   101,  -954,  -954,   900,   901,
     905,     2,   908,   909,  -954,   911,   912,   913,   915,   916,
     920,   921,  4192,  4192,  4192,  4192,  5063,  -954,   923,   928,
    4192,  -358,   929,   101,  -954,  -954,  -954,  -954,  5081,   930,
     931,  -158,   840,   933,  4192,  -358,   938,   939,  4192,  4192,
    4192,   940,   365,  -954,    65,   426,  -954,   435,  -954,   446,
    2091,  1085,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,   308,   650,   835,  1061,  -954,  -954,  -954,  -954,
    4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,  4234,
    4234,  4234,  4234,  -954,  4192,  4192,  4192,  4192,  4192,  4192,
    4192,  4192,  4234,  4234,  4234,  4234,  4192,  4234,  -954,  2503,
    6461,  4915,  2311,  1251,  3373,  1403,  4387,  1715,  4409,  1762,
    1273,  5099,  4431,  1912,  4453,  2130,  4475,  2172,  4497,  2214,
    4519,  2335,  4541,  2446,  1425,  5117,  4563,  2472,   428,  4585,
    2615,  4607,  2665,  4629,  2727,  1462,  5135,  6359,  5153,  6377,
    5171,  6395,  5189,  4651,  2778,  4673,  2814,  4695,  2864,  4717,
    2886,  4739,  2968,  4761,  3243,    84,  -954,   786,  -954,     6,
    -245,  -195,  -954,  -954,  -954,  -954,  -954,    84,  -954,   957,
    -954,  -954,    84,   -10,   959,  -954,  -954,  -954,  -954,  -954,
    -954,   449,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,
    4192,  4192,  4192,  4192,  4192,  4192,   962,  4192,  4192,  4192,
    4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,  4192,
    -954,   963,  -954,   968,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  3280,  4192,  4192,  4192,  4192,  4192,  4192,   969,
     972,   973,  -954,  -954,   960,   976,   977,  -150,   982,   890,
     983,   985,  4192,  -358,   989,   990,   991,   993,   995,  4192,
    4192,  4192,   365,  -954,   471,   193,   997,   129,  1000,  1002,
    1004,  1005,   850,   998,  1012,  1013,  1017,   828,    84,  5207,
    -954,  -954,  5225,  -954,  -954,  -954,  -954,  1018,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,   -33,  -954,  2507,  1020,  1021,
    4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -195,
      84,  -233,  -195,   365,  -195,  -954,  -195,  -954,  -195,  -954,
    -954,  -954,  1022,   -84,  -195,  -954,  1023,  -954,  5243,  -954,
     495,  1027,   517,  5261,  5279,  5297,  -954,  3302,  -954,  1031,
    1032,  1033,  4192,  5315,  1034,   -35,  1038,  -195,  -176,  -358,
    1040,   365,  -954,  5333,  1043,  1042,  -954,  -954,  1048,  1060,
    -954,  -954,  -954,  1049,  1052,  -954,  -954,  -954,  -954,  4192,
    -954,  -954,  -954,  -954,  5351,  -954,     1,   690,  1055,   332,
    1056,  2507,  1057,  1059,  1063,     5,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  1064,   679,     8,  1066,   657,  -234,   658,
    1071,  1019,  1076,  1075,  -954,  -954,  1077,  -954,  -954,  -195,
    -954,  1079,  -954,  5369,  -954,   569,  5387,  5405,  5423,  -954,
    5441,  1086,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    1089,  1090,   -36,  1091,  4192,  -954,  -954,  -954,  -954,  -954,
    1092,  5459,  -954,  1093,  1094,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  4192,  -954,  1096,  -954,  5477,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  4192,  -954,  -954,  5495,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  1098,  1100,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  6447,  6447,  6447,  6447,  4192,
    -954,  -954,  5513,  -954,  -954,  -954,  4192,  -954,  -954,  -195,
    -954,  1101,  -954,  5531,  -954,   590,  -954,  5549,  5567,  5585,
    -954,  -954,  -954,  -954,  1103,  1105,  1106,  -954,  -954,  -954,
     745,  1108,  1110,   -47,  -243,  1111,     4,  -954,   600,  1113,
    1117,   752,  1119,  1120,  1347,  1121,   -47,  -128,  1122,   -36,
    1125,  1132,   -20,  1138,  1139,  3633,  3675,   -26,  1141,   835,
    -954,  -954,   620,  -195,  4192,  4192,  -195,  -954,  1143,  1144,
    -954,  1145,  3067,  3532,  4192,  4192,   -36,  -954,  1146,  1148,
     -24,  -954,  -954,  -954,  1153,  6461,  1526,  1556,  1556,  1556,
    1556,  1556,  1556,  -213,   -45,  -213,   -45,   641,  -199,   641,
    -199,   641,  -199,  2936,  1127,  2265,  2265,  2265,  2265,  2265,
    2265,  -213,   -45,  -213,   -45,   641,  -114,   641,  -106,   642,
     641,  -102,   662,  -954,   645,  -954,  -954,  -954,  -954,  -954,
    -954,  4234,  4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  4234,  4234,  -954,  -954,   651,
    -954,  -954,  -954,  -954,  -954,  -954,  4192,  4192,  4234,  4192,
    4234,  4192,  4234,  4192,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -145,  -195,  -195,
    -195,  -195,  -195,  -195,  -954,  -195,  -195,  -954,    84,  -954,
    -954,    84,  -954,  1151,  1159,  1160,    -2,  -216,  4192,  -954,
    -954,  1161,  4935,  3329,  3351,  3440,  3485,  5603,  3507,  3628,
    3654,  3776,  3798,  3851,  5621,  3887,   647,  3928,  3950,  3972,
    5639,  5657,  5675,  5693,  3994,  4016,  4038,  4060,  4082,  4104,
    -954,  -954,    36,    36,   642,   642,   642,  5711,   653,   660,
     663,  -954,  -954,  -954,  -195,  -954,  -954,  1174,  -954,  -954,
    5729,  -954,  -954,    84,   665,  -954,    84,  -954,    84,  -954,
    5747,  5765,  5783,  -954,  4192,  1184,  -195,  -195,  -954,  -954,
    -954,  -954,  -954,   673,  -954,   674,  1189,  -954,  -954,   801,
    -954,  1191,  -954,  4192,  4192,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  2507,   680,  -954,  5801,  -954,
    -954,  -195,  1072,  -954,  -954,  -954,  -954,  -954,   -81,  1196,
    -954,  -954,  -954,  4192,  1201,  -954,  1203,  4192,  4192,  4192,
    4192,  -954,  -954,  -954,  6447,  4192,  -954,  1204,  -954,   694,
    -954,  -954,  -954,  -954,  -954,   695,  -954,  4192,  -954,  1044,
    -954,  4192,  -954,  -954,  5819,  4192,  -954,  1208,  1210,  1211,
    -954,  1212,  1213,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  1215,  1026,  1218,   838,  1222,  1223,  1170,  1224,  1225,
    1173,  1229,  -954,  -954,  -954,  -954,  -954,  -954,  4192,  1230,
    4192,  4192,  4192,  4192,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  4192,  -954,
    -954,  6447,  -954,  4192,  5837,  4192,  -954,  -954,  5855,  4192,
    5873,  -954,  -954,  4192,  1232,  4192,  4192,  4192,  -954,  -954,
    -954,  1235,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,   717,  -954,  -954,  1239,  -954,  -954,  1238,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,   -23,  6447,   -17,
    6447,  -954,  -954,  -954,  -954,  4192,  -954,  6447,  6447,  -954,
    -954,  -954,  -954,  1241,  -195,  1242,  1244,  1245,  -954,  6447,
    1247,  -195,  1250,  1252,  1257,  -954,  6447,  6447,  6447,  -954,
    -954,  -954,  -954,  -954,  4192,  4783,  4145,  4805,  4187,  4827,
    4213,  5891,  5909,  4849,  4255,  4871,  4277,  4893,  4299,  1258,
    1259,  4192,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  1263,  -954,  -954,  5927,
     874,  -954,  -954,  4192,  -954,  4192,  1264,  1265,  1266,  -954,
    -954,  4192,  -954,  1267,  -954,  -954,  4192,  4192,  4192,  5945,
    -954,  -954,  -954,  4192,  1268,  -954,  1269,  -954,  6447,  6447,
    -954,  1270,  4192,  1147,  1274,  1275,  -954,  -954,  6447,  -954,
    -954,  6447,  5963,  6447,    86,  6447,  -954,  1276,  1277,  6447,
    1281,  6447,  4192,  6447,  -954,  -954,  -954,  -954,  -954,  -954,
    1282,   876,  1283,  -954,   904,  -954,   922,   761,  -954,   925,
    6447,  -954,  6447,  5981,  6447,  5999,  6447,  6447,  4192,  6447,
    4192,  6447,  4192,  6447,  -954,  6447,  6017,  6447,  -954,  1286,
    -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  6447,  -954,  -954,
    -954,  -954,  -954,   800,  4192,  4192,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  6035,  -954,  4192,  1304,  4321,  6053,
     802,   806,   808,  6447,  -954,  6447,  6071,  6447,  4192,  6089,
    -954,   917,  -954,  6447,  1308,  -954,  -954,  4192,  4192,  -954,
    -954,  -954,  6107,  1315,  1319,   935,  1323,  1325,  1326,  1331,
    4192,  4192,  6125,  6143,  6161,  4192,  -954,  4343,  4365,  4192,
    6447,  -954,  -954,  4192,  1332,  1334,  1335,  4192,  6447,  4192,
    1336,  -954,  6179,  6447,  4192,  -954,   974,  1338,   953,   955,
     964,   956,  6197,  6215,  4192,  4192,  4192,  6233,  -954,  -954,
    6447,  6447,  -954,  -954,  -954,  6251,  6269,  -954,  4192,  6447,
    1344,   981,  1348,  1354,  1355,  1357,  4192,  4192,  6447,  6447,
    6447,  4192,  4192,  4192,  6447,  -954,  1359,   999,  1001,   979,
    1003,  6447,   793,  6447,  6447,  6287,   -69,  1360,  1362,  1363,
    1367,  4192,  4192,  1369,  1371,  1228,  1327,  1011,    -9,  6305,
    6447,  -954,  1328,  -954,  -954,  1378,  -954,  1379,  4192,  -954,
    1339,  -954,  6323,  -954,  4192,  6341,  4192,  6447
};

  /* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
     Performed when YYTABLE does not specify something else to do.  Zero
     means the default is an error.  */
static const yytype_int16 yydefact[] =
{
       0,     0,   653,   726,   725,     0,     0,     0,     0,     0,
       0,     0,   683,   684,   685,     0,     0,     0,     0,     0,
     694,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   707,   708,   709,   695,     0,     0,     0,     0,
       0,   484,   485,   259,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   434,   263,     0,     0,   235,    95,     0,
      78,     0,     0,     0,     0,   248,   296,     0,     0,     0,
     289,     0,     0,     0,     0,   114,     0,     0,   433,   107,
     432,    93,     0,    91,    90,    76,    77,    81,   112,     0,
      96,     0,   390,     0,     0,    75,    86,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   757,   758,   759,
     760,     0,     0,     0,     0,   761,   762,   763,   764,     0,
     482,     0,   483,   486,   487,     0,   317,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   478,     0,     0,     0,     0,    26,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   104,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   210,     0,     0,
       0,   206,     0,     0,     0,   212,     0,   120,     0,     0,
       0,     0,     0,     0,   631,   630,   488,   491,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   741,   742,
     743,     0,     0,     0,     0,     0,   765,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   778,   779,
     780,   766,     0,     0,     0,     0,     0,    41,     0,    32,
       0,     0,     0,     0,     0,     0,    31,     0,    89,    28,
      44,   437,    56,   105,     0,     0,     0,     0,     0,   108,
       0,   180,     0,   137,   123,     0,     0,     0,     0,     0,
       0,     0,     0,   430,   429,   426,   428,   427,     0,     0,
       0,     0,     0,     0,   458,    80,   258,     0,     0,     0,
       0,     0,   465,   479,   481,   241,   256,     0,   250,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   725,     0,     0,   272,     0,     0,   291,     0,     0,
       0,     0,     0,   423,   424,     0,     0,    55,    46,     0,
       0,     0,    45,     0,     0,   425,     0,     0,     0,     0,
       0,     0,     0,     0,   357,   350,   100,     0,    97,    98,
      99,   101,     0,     0,     0,     0,   118,     0,     0,     0,
       0,    92,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   396,   581,   582,   431,   435,     0,
       0,     0,     0,     0,    60,     0,     0,     0,     0,     0,
      59,     0,   453,     0,     0,     0,     0,    87,   117,     0,
       0,     0,     0,     0,   364,   358,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   365,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   388,   376,
     448,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   349,   342,     0,    29,
       0,     0,     0,     0,   421,     0,     0,    42,     0,    57,
       0,     0,     0,     0,     0,     0,     0,    79,     0,     0,
       0,     0,     0,     0,   341,   334,   489,   490,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   323,   653,     0,   724,   654,   722,   654,
       0,     0,     1,     8,     7,    25,     9,    14,    10,    13,
      18,    15,    16,    24,    23,    22,    19,    21,    20,    17,
      11,    12,     0,     0,     0,     0,   477,   502,     3,     4,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     6,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     5,   651,
     652,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   205,     0,   629,     0,
       0,     0,   628,   626,   627,   194,   208,     0,   218,     0,
     214,   217,     0,     0,     0,   492,   495,   497,   494,   493,
      85,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      37,   637,   638,   636,    38,    35,    36,   633,   632,    34,
     797,   795,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   262,   264,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   178,     0,   176,   174,     0,     0,
       0,     0,     0,   124,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     232,   234,     0,   237,   238,   233,   236,     0,   462,   463,
     460,   461,   459,   635,   634,     0,   469,     0,     0,     0,
       0,   472,   470,   466,   587,   589,   588,   586,   590,   471,
     749,   750,   751,   752,   753,   754,   755,   756,   480,     0,
       0,     0,     0,     0,     0,   243,     0,   255,     0,   249,
     299,   306,     0,     0,     0,   303,     0,   307,     0,   301,
       0,     0,     0,     0,     0,     0,   297,     0,   106,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   267,     0,     0,     0,   277,   276,     0,     0,
     273,   290,    52,     0,     0,    51,    50,    49,    48,     0,
     103,   115,   352,   351,     0,   354,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   353,   647,   649,   648,
     646,   356,   102,     0,   265,     0,     0,     0,     0,     0,
       0,   157,     0,     0,    88,   113,     0,   393,   398,     0,
     397,     0,   399,     0,   394,     0,     0,     0,     0,   391,
       0,     0,    64,    74,    66,    67,    63,    62,    65,    68,
       0,     0,     0,     0,     0,   456,   454,   457,   455,   359,
       0,     0,   362,     0,     0,   370,   372,   375,   373,   369,
     368,   371,   366,     0,   384,     0,   377,     0,   386,   382,
     379,   381,   380,   383,   378,   451,   449,   452,   450,   445,
     444,   447,   446,     0,   345,   343,     0,   346,   344,   348,
     119,    30,   261,   122,     0,     0,   231,   422,   242,   295,
      43,   389,    58,   116,   316,   440,   439,   442,   441,     0,
     337,   335,     0,   338,   336,   340,     0,   320,   327,     0,
     324,     0,   328,     0,   321,     0,   326,     0,     0,     0,
     325,   318,   723,   796,     0,     0,     0,   501,   577,   566,
       0,     0,     0,     0,     0,     0,     0,   500,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   499,
     534,   536,     0,     0,     0,     0,     0,   593,     0,     0,
     591,     0,     0,     0,     0,     0,     0,   592,     0,     0,
       0,   498,   514,   503,   721,   654,   720,   714,   715,   716,
     717,   718,   719,   656,   654,   660,   654,   664,   654,   668,
     654,   674,   654,   794,   793,   787,   788,   789,   790,   791,
     792,   657,   654,   661,   654,   665,   654,   669,   654,   732,
     672,   654,   727,   676,   675,   677,   735,   678,   736,   679,
     737,     0,     0,   681,   739,   682,   740,   686,   744,   687,
     745,   688,   746,   689,   747,     0,     0,   696,   767,   697,
     698,   769,   699,   770,   700,   771,     0,     0,     0,     0,
       0,     0,     0,     0,   705,   776,   706,   777,   710,   781,
     711,   782,   712,   783,   713,   784,   211,     0,     0,     0,
       0,     0,     0,     0,   195,     0,     0,   204,     0,   207,
     216,     0,   213,     0,     0,     0,     0,     0,     0,   219,
     496,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      39,    40,   728,   729,   730,   731,   733,     0,     0,     0,
       0,   125,   183,   187,     0,   186,   169,     0,   164,   188,
       0,   184,   165,     0,     0,   168,     0,   167,     0,   166,
       0,     0,     0,   181,     0,     0,     0,     0,   147,   143,
     139,   141,   144,     0,   149,     0,     0,   142,   151,     0,
     152,     0,   138,     0,     0,   464,   467,   468,   606,   611,
     610,   612,   601,   607,   597,   620,   616,   621,   617,   622,
     618,   598,   600,   613,   604,   614,   602,   615,   623,   619,
     603,   605,   609,   608,   599,     0,     0,   476,     0,   244,
     257,     0,     0,   247,   251,   254,   253,   252,     0,     0,
     310,   302,   308,     0,     0,   309,     0,     0,     0,     0,
       0,   286,   288,   287,   283,     0,   275,     0,   270,     0,
     269,   624,   625,   268,   284,     0,   292,     0,   278,   280,
     274,     0,    53,    54,     0,     0,   418,     0,     0,     0,
     405,   637,   636,   414,   415,   409,   408,   407,   406,   410,
     412,     0,     0,     0,     0,     0,     0,   159,     0,     0,
     161,     0,   158,   260,    82,   136,   395,   400,     0,     0,
       0,     0,     0,     0,    73,    72,    71,   639,   640,   645,
     644,   641,   642,   643,    69,    70,    61,   360,     0,   361,
     374,   367,   385,     0,     0,     0,   163,   135,     0,     0,
       0,   322,   329,     0,     0,     0,     0,     0,   578,   579,
     580,     0,   575,   574,   585,   584,   583,   572,   568,   567,
     576,     0,   569,   570,     0,   547,   542,     0,   563,   562,
     541,   540,   560,   548,   544,   546,   545,   564,   551,   539,
     543,   595,   596,   594,   552,   549,   550,     0,   553,     0,
     554,   537,   538,   559,   535,     0,   515,   512,   511,   510,
     519,   513,   520,     0,     0,     0,     0,     0,   504,   506,
       0,     0,     0,     0,     0,   505,   507,   508,   509,   529,
     516,   532,   530,   531,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   198,   199,   196,   201,   200,   197,   202,   203,
     209,   215,   221,   226,   225,   220,     0,   224,   223,     0,
      83,   727,   734,     0,   768,     0,     0,     0,     0,   185,
     189,     0,   179,     0,   177,   175,     0,     0,     0,     0,
     145,   146,   148,     0,     0,   155,     0,   140,   240,   239,
     473,     0,     0,     0,     0,     0,   304,   311,   314,   298,
     300,   313,     0,   315,   785,   282,   271,     0,     0,   285,
       0,   281,     0,   355,   420,   419,   411,   416,   417,   413,
       0,     0,     0,   126,     0,   160,     0,     0,   162,     0,
     403,   392,   402,     0,   404,     0,   363,   387,     0,   347,
       0,   339,     0,   332,   319,   331,     0,   333,   573,     0,
     565,   561,   556,   555,   558,   557,   533,   521,   527,   525,
     523,   517,   522,   528,   526,   524,   518,   650,   680,   738,
     690,   692,   691,   693,     0,     0,   702,   773,   703,   774,
     704,   775,   229,   228,     0,   222,     0,     0,     0,     0,
       0,     0,     0,   192,   182,   191,     0,   193,     0,     0,
     154,     0,   474,   475,     0,   245,   305,     0,     0,   294,
     293,   279,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   571,     0,     0,     0,
     227,    84,   748,     0,     0,     0,     0,     0,   153,     0,
       0,   246,     0,   786,     0,   121,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   701,   772,
     230,    33,   109,   111,   110,     0,     0,   150,     0,    47,
       0,     0,     0,     0,     0,     0,     0,     0,   443,   438,
      94,     0,     0,     0,   312,   266,     0,     0,     0,     0,
       0,   401,     0,   330,   190,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   133,   170,     0,   131,     0,
     156,   129,   127,   134,   171,     0,   132,     0,     0,   128,
     172,   130,     0,   173,     0,     0,     0,   436
};

  /* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,    17,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -354,  -954,  -954,
    -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,  -954,
    -954,  -954,   -19,  -954,  -954,  1329,  1333,   -12,  1099,  -954,
    -954,   353,  -954,  -954,  -954,   369,  -954,  -954,  -954,  -837,
    -954,  -954,   697,    -3,  -297,   540,  -953,  -410,  -954,  -954,
      61,     0
};

  /* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
      -1,   130,   131,   132,   133,   379,   356,   134,   748,  1238,
    1236,  1233,   732,   197,   647,  1168,   635,   652,  1171,   135,
     282,   136,   314,   793,   790,   137,   268,   138,   342,   831,
     139,   328,   140,   512,   141,   142,   143,   144,   145,   146,
     147,   393,   364,   148,   149,   150,   151,   206,   546,  1061,
    1039,  1040,  1017,  1007,   397,  1427,   779,  1062,  1454,  1295,
    1333,   645,   207,   699,   766,   694,  1394,   871,   152,   153,
     154,  1065
};

  /* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
     positive, shift that token.  If negative, reduce the rule whose
     number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
     155,   772,   851,   398,  1525,    41,    42,   842,  1509,  1346,
     963,  1387,  1431,  1360,  1364,  1165,  1309,   412,   763,  1565,
     175,    44,  1173,   315,  1356,   885,   334,   395,  1451,   649,
     834,   856,   909,  1331,  1266,   763,   542,  1158,   875,   296,
     450,  1461,  1424,  1492,  1622,  1388,   257,   372,  1389,  1390,
    1624,   335,  1786,   259,   819,   959,  1368,   316,  1773,   857,
     208,   377,  1428,   204,  1159,  1327,   205,   156,   858,   333,
    1429,  1432,   353,   274,  1527,   159,  1448,   371,   697,   157,
     160,  1444,   424,   985,   448,   763,   835,   466,   763,   209,
     698,   889,   161,   275,   317,   494,   162,   411,   158,   650,
     258,   204,   436,  1489,   205,   459,   989,   496,   163,   497,
     276,   413,   486,  -741,  1224,   876,   164,   418,  1787,   774,
     775,   776,   165,   204,   910,   498,   205,   517,   519,   521,
     378,   318,   357,  1160,   451,  1250,  -742,  1251,   204,   373,
    -743,   205,   859,  1174,  -765,  1510,   204,   297,   374,   205,
     860,   166,   861,   836,   691,   167,   837,   580,   581,   583,
     585,   587,   589,   591,   593,   595,   597,   599,   601,   603,
     605,   607,  1774,   610,   612,   298,   614,   616,   618,   620,
     622,   624,   626,   628,   630,   632,   634,  1441,   516,   518,
     520,  1425,   336,   414,   543,   375,  1301,   277,   911,   867,
     820,  1166,  1391,   821,   168,   299,   319,   877,   396,   320,
     260,   169,   321,   170,   862,   171,   452,   863,   579,   172,
     582,   584,   586,   588,   590,   592,   594,   596,   598,   600,
     602,   604,   606,  1369,   609,   611,   173,   613,   615,   617,
     619,   621,   623,   625,   627,   629,   631,   633,  1678,   713,
    1445,  1528,   288,   700,   701,   702,   204,   337,  1332,   205,
     415,   204,   174,   756,   205,   419,   964,   204,  1452,   177,
     205,  1161,   358,  1453,   764,   749,   753,   289,   752,   912,
    1175,  1156,   178,   453,   359,   322,   204,   360,   179,   205,
    -778,   764,   765,  1169,   878,   799,   864,   180,  1172,   323,
     560,   561,   204,   562,   181,   205,  1446,  -779,  1426,   816,
     343,   344,   278,    69,    70,   279,   576,   577,   809,  1392,
     808,  1162,  1511,   841,  -780,   813,   814,   815,   204,   817,
     777,   205,  1176,   182,   868,   823,  -766,   833,   778,   420,
     839,   764,   183,   184,   764,  1267,   692,   263,   848,   261,
     280,   290,  1462,  1365,  1493,  1623,   185,   855,   186,   854,
     869,  1625,   264,   262,   693,   361,   324,   544,   265,   545,
     865,   843,   266,  1347,   899,   870,   267,  1361,   269,   338,
     283,  1393,  1433,   416,   339,   894,   295,   893,  1245,   879,
     896,   897,   898,   325,  1262,   362,   913,   300,   329,   900,
     376,  -730,   577,   330,   907,   326,   454,   340,  1004,  -731,
     577,   119,   120,  -733,   577,   354,  1177,   421,   380,   922,
     422,   921,   965,   291,   327,   844,   292,  1348,   381,   930,
     880,   417,  1351,   522,   523,  1163,  1300,  1246,   938,   524,
     937,   121,   122,   399,   341,   525,   526,   123,   124,   527,
     881,   528,   529,  1247,   530,   531,   532,  1178,  1560,   363,
     957,   204,   956,   533,   205,  1526,   284,   534,   574,   575,
     576,   577,   608,   204,   636,   535,   205,   536,   537,   637,
     538,   285,   975,   976,   977,   978,   281,   646,   983,   204,
     982,   539,   205,  1001,   540,   648,   301,   541,   361,   651,
     653,   548,   994,   460,   993,   549,  1310,   210,   997,   998,
     999,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,   176,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   705,
     706,   576,   707,   654,  1005,   655,   656,  1006,  1074,  1076,
    1078,  1080,  1082,   658,  1083,  1084,  1085,  1086,  1087,  1088,
    1089,  1090,  1092,  1094,  1096,  1098,  1099,  1101,   659,   660,
     661,   204,   663,   423,   205,   488,   662,   664,   802,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,   293,   204,   665,   666,   205,   667,   468,
     668,  1064,  1066,  1067,  1068,  1069,  1070,  1071,  1072,  1073,
    1075,  1077,  1079,  1081,   692,   669,   670,   671,   204,   672,
     673,   205,   674,  1091,  1093,  1095,  1097,   690,  1100,   675,
     676,   677,  1352,   461,   366,   469,   678,   767,   679,   680,
     681,   286,   695,   682,   683,   462,  1009,   684,   463,   685,
     638,   686,  1182,  1183,  1184,  1185,  1186,  1187,  1188,  1189,
    1190,  1191,  1192,  1193,  1194,  1195,   687,  1197,  1198,  1199,
    1200,  1201,  1202,  1203,  1204,  1205,  1206,  1207,  1208,  1209,
     688,   689,   696,   709,   301,   710,   711,   714,   361,   204,
     332,   470,   205,  1212,  1213,  1214,  1215,  1216,  1217,   733,
     471,   754,   757,  1243,   758,   760,   472,    67,   761,   762,
    1231,   499,  1230,   789,   734,   489,   639,   367,   791,  1240,
    1241,  1242,   792,   794,   796,   798,   361,   490,   800,   801,
     491,  1010,   368,   256,   803,   187,   806,   807,   640,   305,
    1011,   810,   811,   473,   474,   500,   287,   475,   294,   204,
     302,   812,   205,   476,   818,   824,   464,   198,   838,   352,
    1298,   355,   840,   365,  1304,   845,   199,   846,   847,   852,
     853,   204,   866,   872,   205,   873,   187,   188,   874,   394,
     883,   768,   501,   882,   410,   884,   886,   641,   425,   435,
     449,   804,   477,   467,   887,   888,   892,   895,   891,   901,
     487,   495,  1336,   189,  1520,   902,   903,  1521,   361,   306,
     904,   307,  1324,   513,   190,   905,  1334,   308,   188,   502,
     465,   906,  1012,   204,   908,   915,   205,   916,   917,   918,
     919,   923,   191,  1018,   924,   192,   920,   547,   492,  1344,
     925,   926,   927,   735,   189,   928,   929,   931,   200,   369,
     932,   934,   642,   935,   193,   190,   936,   643,   309,   939,
     644,   736,   940,   941,   942,   201,   943,   194,   310,  1542,
     944,  1013,  1544,   191,  1545,  1014,   192,   945,   946,   769,
     947,   948,   949,   737,  1019,   311,  1020,   950,   657,  1015,
     951,   952,   954,   955,   503,   193,   958,   504,   960,   961,
     505,   738,   493,   962,  1396,   312,   966,   967,   194,   968,
     969,   970,   202,   971,   972,   203,  1021,   270,   973,   974,
     770,   980,   370,  1401,   478,  1022,   981,   984,   987,   988,
     991,   992,   107,   108,   109,   110,   995,   996,  1000,   176,
    -797,   204,  1129,  1404,   205,   479,   115,   116,   117,   118,
     739,  -795,   480,  1157,   712,  1170,  1221,  1180,  1016,  1181,
    1196,  1210,   204,   750,   751,   205,  1211,  1218,   755,  1408,
    1219,  1220,   740,   506,  1222,  1223,  1410,   759,  1226,  1228,
    1227,  1244,   271,  1229,   771,  1232,   773,  1235,  1234,  1237,
     547,  1239,   195,  1249,  1257,  1261,   741,   795,  1252,   797,
    1253,  1566,  1254,  1255,   805,  1314,   253,  1023,   507,  1024,
    1258,  1259,   254,  1256,   382,  1260,  1265,   255,  1296,  1297,
    1308,  1312,   313,  1315,   832,  1458,  1460,  1316,  1025,  1321,
    1322,  1323,  1326,   195,  1467,  1468,  1329,   196,  1335,  1338,
    1339,   850,  1479,  1486,  1487,  1488,  1340,  1342,   383,   481,
    1343,   742,  1349,  1350,  1355,  1357,  1026,  1358,  1363,  1042,
    1027,  1359,  1362,   715,  1366,   743,   744,  1367,  1370,  1371,
    1372,   890,  1373,  1374,  1028,  1375,   272,  1377,   196,  1379,
    1043,   508,  1044,  1045,  1384,   384,   437,  1385,  1386,  1395,
    1397,  1399,  1400,   509,  1402,   745,  1406,   716,  1407,  1412,
    1414,  1418,  1496,  1419,  1420,  1046,  1422,  1421,  1423,  1430,
    1434,  1435,   510,  1436,  1437,  1498,  1500,  1438,  1439,  1442,
    1447,   204,   385,  1449,   205,  1029,  1501,  1502,  1450,  1504,
    1465,  1506,   746,  1508,   717,  1047,  1455,  1456,   273,  1463,
     747,  1470,  1471,  1472,  1490,   511,  1491,   562,   707,  1522,
    -734,  1048,   718,  1030,  1494,  1031,  -768,  1523,  1524,  1530,
     438,  1534,  1495,  1536,   780,   781,   782,   783,  1529,  1049,
    1537,   719,  1540,  1538,   720,  1543,  1497,  1499,   784,   785,
     786,   787,  1550,  1553,  1554,  1555,  1556,  1557,   990,  1503,
    1561,  1505,  1771,  1507,  1567,  1564,  1032,   386,  1033,  1569,
     387,  1570,  1576,   388,  1577,  1578,  1584,  1580,  1585,  1586,
    1587,  1588,  1034,  1589,  1035,  1590,  1591,  1592,  1593,  1594,
    1036,  1595,  1596,  1597,  1598,  1599,   439,  1619,  1601,  1008,
    1614,  1041,  1063,  1618,  1549,  1620,  1621,  1657,   440,  1627,
    1629,   441,  1630,  1631,   204,  1632,   721,   205,  1634,   722,
    1635,  1037,   723,  1558,  1559,  1636,  1652,  1653,   724,  1655,
    1038,  1684,  1660,  1661,  1662,  1664,  1670,  1671,  1672,  1050,
    1674,  1688,  1675,  1676,  1679,  1680,   389,  1681,  1683,  1783,
    1710,  1685,   204,  1686,  1696,   205,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
    1051,  1687,  1701,  1568,  1689,  -748,  1711,  1571,  1572,  1573,
    1574,  1715,  1704,  1052,   442,  1575,  1705,  1716,  1706,   361,
    1717,  1718,  1053,  1719,  1720,   725,  1164,  1579,  1167,  1721,
    1732,  1581,  1733,  1734,  1737,  1583,  1741,  1740,  1742,  1743,
    1179,  1745,  1755,  1744,  1756,  1440,  1757,   204,  1054,  1055,
     205,  1056,  1758,  1759,  1057,  1760,   726,  1766,  1775,   727,
    1776,  1777,  1767,  1769,  1768,  1778,  1770,  1781,  1600,  1782,
    1602,  1603,  1604,  1605,  1785,   443,  1790,  1791,  1784,  1789,
     444,   445,  1464,   303,   390,  1443,     0,   304,  1606,  1354,
    1793,     0,   788,  1607,     0,  1609,   391,     0,     0,  1611,
       0,     0,     0,  1613,  1225,  1615,  1616,  1617,     0,     0,
       0,     0,     0,   446,     0,   392,     0,     0,     0,     0,
       0,     0,  1248,  1058,     0,     0,     0,   728,     0,  1059,
       0,   447,     0,   729,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1341,   730,  1268,  1269,  1270,  1271,
       0,  1272,  1273,  1274,     0,  1626,  1275,  1276,  1277,  1278,
    1279,  1280,     0,     0,   731,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1299,  1060,  1302,  1303,
       0,  1305,     0,  1306,  1637,  1307,     0,     0,     0,     0,
       0,  1311,     0,     0,     0,     0,     0,     0,     0,     0,
       0,  1654,     0,     0,  1281,     0,     0,     0,     0,     0,
       0,     0,  1328,     0,  1330,     0,     0,     0,     0,     0,
       0,     0,     0,  1658,     0,  1659,     0,     0,     0,  1282,
       0,  1663,     0,     0,     0,     0,  1665,  1666,  1667,     0,
       0,     0,     0,  1669,     0,  1283,  1353,     0,     0,     0,
       0,     0,  1673,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,  1682,     0,     0,     0,  1376,     0,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,  1692,  1003,
    1693,     0,  1694,     0,     0,     0,     0,     0,     0,     0,
       0,     0,  1284,  1285,  1286,  1287,  1288,  1289,  1290,  1291,
    1292,  1293,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,  1697,  1698,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1700,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,  1708,   562,
       0,     0,     0,     0,  1294,     0,     0,  1712,  1713,     0,
       0,     0,     0,     0,     0,     0,  1411,     0,     0,     0,
    1722,  1723,     0,     0,     0,  1727,     0,     0,     0,  1730,
       0,     0,     0,  1731,     0,     0,     0,  1735,     0,  1736,
       0,     0,     0,     0,  1739,     0,     0,     0,     0,     0,
       0,     0,     0,     0,  1748,  1749,  1750,     0,     0,     0,
       0,     0,     0,     0,     0,     0,  1041,     0,  1754,     0,
    1466,     0,     0,  1469,     0,     0,  1761,  1762,     0,  1478,
    1485,  1763,  1764,  1765,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     0,     0,
       0,  1779,  1780,     0,     0,  1104,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,  1792,   562,
       0,     0,     0,  1111,  1795,     0,  1797,    -2,     1,     0,
       2,     0,     0,     0,     3,     4,     5,     6,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
       0,    29,    30,    31,    32,    33,    34,    35,    36,    37,
      38,    39,     0,     0,     0,     0,     0,     0,    40,     0,
       0,     0,     0,    41,    42,  1512,  1513,  1514,  1515,  1516,
    1517,     0,  1518,  1519,     0,    43,     0,     0,     0,    44,
       0,     0,     0,     0,    45,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    46,     0,     0,     0,     0,    47,
       0,    48,     0,     0,    49,     0,    50,    51,    52,     0,
       0,     0,     0,     0,     0,    53,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,  1539,     0,     0,     0,     0,     0,  1106,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
      54,   562,     0,  1551,  1552,  1125,    55,     0,     0,    56,
       0,     0,     0,     0,     0,    57,     0,     0,     0,     0,
       0,     0,     0,    58,     0,   550,   551,   552,   553,   554,
     555,   556,   557,   558,   559,   560,   561,    59,   562,     0,
       0,     0,  1136,     0,     0,    60,    61,     0,     0,     0,
       0,    62,     0,    63,    64,     0,    65,     0,  1563,     0,
       0,     0,    66,     0,     0,     0,     0,     0,     0,     0,
      67,     0,     0,     0,     0,     0,    68,     0,     0,     0,
       0,     0,     0,     0,     0,     0,    69,    70,     0,    71,
       0,   552,   553,   554,   555,   556,   557,   558,   559,   560,
     561,    72,   562,     0,     0,     0,     0,    73,     0,     0,
       0,    74,     0,     0,     0,     0,     0,     0,    75,     0,
      76,  -798,  -798,  -798,  -798,  -798,  -798,   558,   559,   560,
     561,   210,   562,     0,     0,     3,     4,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,    77,   235,   236,   237,   238,   239,   240,   241,   242,
     243,   244,   245,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    78,     0,     0,     0,     0,    79,     0,
       0,     0,     0,     0,     0,     0,    80,     0,    81,    82,
       0,    83,     0,     0,     0,     0,     0,     0,     0,    84,
      85,    86,     0,     0,     0,    87,     0,     0,     0,    88,
       0,     0,     0,    89,     0,     0,    90,     0,     0,     0,
       0,  1628,    91,    92,    93,     0,     0,     0,  1633,     0,
      94,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    95,    96,     0,    97,     0,    98,     0,     0,     0,
      99,   100,     0,     0,     0,     0,     0,   101,     0,     0,
     102,   103,   104,     0,     0,     0,     0,     0,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,   105,   106,   107,   108,   109,   110,  1108,
       0,     0,     0,     0,   111,   112,   113,     0,   114,   115,
     116,   117,   118,     0,     0,     0,     0,     0,     0,   119,
     120,     0,     0,     0,   332,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,    67,     0,     0,     0,     0,  1110,     0,     0,   121,
     122,     0,     0,     0,     0,   123,   124,   125,   126,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   127,
       0,     0,     0,     0,     0,   128,   -27,   210,     0,     0,
     129,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,     0,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,     0,
     210,     0,     0,     0,     3,     4,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
       0,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,     0,     0,     0,   400,     0,     0,     0,     0,
       0,     0,     0,     0,   401,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,   402,     0,     0,     0,     0,  1114,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   426,     0,
       0,     0,     0,     0,   210,     0,     0,   427,     3,     4,
     211,   212,   213,   214,   215,   216,   217,   218,   219,   220,
     221,   222,   223,   224,   225,   226,   227,   228,   229,   230,
     231,   232,   233,   234,     0,   235,   236,   237,   238,   239,
     240,   241,   242,   243,   244,   245,   403,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   107,   108,   109,   110,
     404,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     115,   116,   117,   118,     0,     0,     0,    67,     0,     0,
       0,     0,     0,   455,   456,     0,     0,     0,     0,   428,
       0,   429,     0,     0,     0,   405,     0,     0,     0,     0,
       0,   246,     0,   332,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   406,   457,   458,     0,     0,   247,   407,
      67,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     253,     0,     0,     0,     0,     0,   254,     0,   430,     0,
       0,   255,     0,     0,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   431,   562,     0,     0,
     204,     0,   432,   205,     0,  1002,  1268,  1269,  1270,  1271,
       0,  1272,  1273,  1274,     0,   248,  1275,  1276,  1277,  1278,
    1279,  1280,     0,   564,   565,   566,   567,   568,   569,   570,
     571,   572,   573,   574,   575,   576,   577,   249,     0,     0,
       0,     0,     0,   204,  1116,     0,   205,     0,     0,     0,
       0,     0,     0,     0,    67,     0,     0,     0,     0,     0,
       0,   408,     0,     0,  1281,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,     0,     0,     0,     0,     0,  1118,   250,     0,  1282,
     251,     0,     0,     0,     0,     0,   252,     0,     0,     0,
       0,     0,     0,     0,   433,  1283,   409,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,     0,  1120,     0,
       0,     0,     0,     0,     0,     0,   434,   204,     0,     0,
     205,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
    -798,  -798,  -798,  -798,  -798,  -798,   703,   704,   705,   706,
     576,   707,  1284,  1285,  1286,  1287,  1288,  1289,  1290,  1291,
    1292,  1293,     0,     0,     0,   107,   108,   109,   110,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   115,
     116,   117,   118,     0,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,   253,   562,     0,     0,
       0,     0,   254,     0,  1294,  1103,     0,   255,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1122,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   253,
       0,     0,     0,     0,     0,   254,     0,     0,     0,     0,
     255,     0,     0,     0,     0,     0,     0,     0,     0,   107,
     108,   109,   110,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   115,   116,   117,   118,   210,     0,     0,
       0,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,     0,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   253,     0,     0,     0,   345,     0,   254,
    1124,     0,     0,     0,   255,   564,   565,   566,   567,   568,
     569,   570,   571,   572,   573,   574,   575,   576,   577,     0,
       0,     0,     0,     0,     0,     0,  1128,     0,     0,     0,
       0,     0,     0,     0,   346,     0,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,   210,   562,
       0,   347,     3,     4,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,   235,
     236,   237,   238,   239,   240,   241,   242,   243,   244,   245,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     210,     0,     0,     0,     3,     4,   211,   212,   213,   214,
     215,   216,   217,   218,   219,   220,   221,   222,   223,   224,
     225,   226,   227,   228,   229,   230,   231,   232,   233,   234,
     348,   235,   236,   237,   238,   239,   240,   241,   242,   243,
     244,   245,     0,     0,     0,     0,     0,    67,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1131,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     349,     0,     0,   350,     0,     0,     0,     0,     0,   351,
       0,     0,     0,     0,     0,     0,     0,  1473,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1133,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     204,     0,     0,   205,   825,     0,     0,     0,     0,     0,
       0,   332,     0,     0,     0,     0,     0,     0,   826,     0,
       0,     0,     0,     0,     0,     0,     0,     0,    67,     0,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,     0,     0,     0,     0,     0,   827,
     828,  1135,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   332,     0,     0,  1474,     0,     0,     0,
     829,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      67,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,     0,     0,     0,     0,     0,
       0,     0,  1145,     0,     0,     0,     0,     0,     0,     0,
       0,   204,     0,  1475,   205,     0,  1476,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,     0,  1147,     0,
       0,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
       0,     0,     0,   204,     0,     0,   205,   564,   565,   566,
     567,   568,   569,   570,   571,   572,   573,   574,   575,   576,
     577,     0,     0,     0,     0,     0,     0,   830,  1149,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,     0,     0,     0,     0,     0,     0,     0,
    1151,     0,     0,     0,     0,     0,   253,     0,     0,     0,
       0,     0,   254,     0,     0,     0,     0,   255,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,  1477,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,   107,   108,   109,   110,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   115,   116,   117,
     118,   564,   565,   566,   567,   568,   569,   570,   571,   572,
     573,   574,   575,   576,   577,     0,     0,     0,     0,     0,
       0,     0,  1153,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,   107,   108,   109,   110,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   115,
     116,   117,   118,     0,     0,     0,     0,   253,     0,     0,
       0,     0,     0,   254,     0,   210,     0,     0,   255,     3,
       4,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,     0,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   245,     0,     0,   253,
       0,     0,     0,     0,     0,   254,     0,     0,     0,   210,
     255,     0,     0,     3,     4,   211,   212,   213,   214,   215,
     216,   217,   218,   219,   220,   221,   222,   223,   224,   225,
     226,   227,   228,   229,   230,   231,   232,   233,   234,     0,
     235,   236,   237,   238,   239,   240,   241,   242,   243,   244,
     245,     0,  1480,     0,     0,     0,   210,     0,     0,     0,
       3,     4,   211,   212,   213,   214,   215,   216,   217,   218,
     219,   220,   221,   222,   223,   224,   225,   226,   227,   228,
     229,   230,   231,   232,   233,   234,     0,   235,   236,   237,
     238,   239,   240,   241,   242,   243,   244,   245,   210,     0,
       0,     0,     3,     4,   211,   212,   213,   214,   215,   216,
     217,   218,   219,   220,   221,   222,   223,   224,   225,   226,
     227,   228,   229,   230,   231,   232,   233,   234,     0,   235,
     236,   237,   238,   239,   240,   241,   242,   243,   244,   245,
       0,     0,     0,     0,     0,     0,     0,     0,   332,     0,
       0,  1481,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    67,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,     0,     0,     0,     0,  1155,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,  1482,     0,
       0,  1483,   332,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,    67,
       0,     0,     0,     0,  1003,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1320,     0,   204,   332,
       0,   205,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,    67,     0,     0,     0,
       0,     0,     0,  1532,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   332,     0,     0,     0,  1106,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,    67,   562,
       0,     0,     0,     0,     0,     0,     0,  1105,     0,     0,
       0,   210,     0,     0,  1484,     3,   331,   211,   212,   213,
     214,   215,   216,   217,   218,   219,   220,   221,   222,   223,
     224,   225,   226,   227,   228,   229,   230,   231,   232,   233,
     234,     0,   235,   236,   237,   238,   239,   240,   241,   242,
     243,   244,   245,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1108,     0,     0,     0,     0,     0,
     107,   108,   109,   110,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   115,   116,   117,   118,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1110,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,   107,   108,   109,   110,     0,     0,
       0,  1114,     0,     0,     0,     0,     0,     0,   115,   116,
     117,   118,     0,     0,   253,     0,     0,     0,     0,     0,
     254,   482,   483,     0,     0,   255,     0,     0,     0,  1457,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   107,   108,   109,   110,     0,     0,     0,     0,     0,
       0,     0,   484,   485,     0,   115,   116,   117,   118,     0,
       0,     0,     0,     0,   332,     0,     0,     0,   253,     0,
       0,  1459,     0,     0,   254,     0,     0,     0,     0,   255,
       0,    67,     0,   107,   108,   109,   110,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   115,   116,   117,
     118,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,   253,     0,     0,     0,     0,
       0,   254,  1116,     0,     0,     0,   255,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1118,     0,
       0,     0,     0,     0,     0,     0,     0,   253,     0,     0,
       0,     0,     0,   254,     0,   210,     0,     0,   255,     3,
       4,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,     0,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   245,   514,     0,     0,
       0,     3,     4,     5,     6,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,   515,    26,    27,    28,     0,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,     0,     0,     0,     0,     0,     0,     0,
    1120,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,     0,     0,     0,
       0,     0,  1122,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   107,   108,   109,   110,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     115,   116,   117,   118,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,     0,     0,     0,     0,  1124,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   332,     0,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,    67,     0,     0,     0,     0,
     253,  1128,     0,     0,     0,     0,   254,     0,     0,     0,
       0,   255,     0,     0,     0,     0,     0,     0,     0,     0,
     332,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,    67,     0,     0,
       0,     0,  1131,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1133,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1135,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1145,   564,
     565,   566,   567,   568,   569,   570,   571,   703,   704,   705,
     706,   576,   707,     0,     0,     0,     0,     0,     0,     0,
    1147,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,     0,     0,     0,     0,     0,
       0,     0,  1149,   564,   565,   566,   567,   568,   569,   570,
     571,   703,   704,   705,   706,   576,   707,     0,     0,     0,
       0,     0,     0,     0,  1151,   564,   565,   566,   567,   568,
     569,   570,   571,   703,   704,   705,   706,   576,   707,     0,
       0,     0,     0,     0,     0,     0,  1153,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1155,     0,
     107,   108,   109,   110,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   115,   116,   117,   118,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1639,
       0,     0,   107,   108,   109,   110,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,   115,   116,   117,   118,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,   253,     0,     0,     0,     0,     0,
     254,  1641,     0,     0,     0,   255,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,     0,     0,     0,     0,  1643,     0,     0,
       0,     0,     0,     0,     0,     0,   127,     0,     0,     0,
       0,     0,   128,     0,     0,     0,     0,   129,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1647,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,     0,     0,     0,
       0,  1649,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,     0,
       0,     0,     0,  1651,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,     0,     0,     0,     0,  1702,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,     0,     0,     0,     0,  1728,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,     0,     0,  1729,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1107,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1109,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1113,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1115,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1117,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1119,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1121,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1123,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1127,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1130,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1132,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1134,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1144,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1146,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1148,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1150,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1152,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1154,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1638,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,     0,     0,     0,     0,  1640,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,     0,     0,     0,
       0,  1642,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,     0,
       0,     0,     0,  1646,   550,   551,   552,   553,   554,   555,
     556,   557,   558,   559,   560,   561,     0,   562,     0,     0,
       0,     0,     0,     0,     0,  1648,   550,   551,   552,   553,
     554,   555,   556,   557,   558,   559,   560,   561,     0,   562,
       0,     0,     0,     0,     0,     0,     0,  1650,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,  1102,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,     0,     0,  1531,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,   708,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,   822,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   849,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,   914,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,   933,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,   953,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,   979,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,   986,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,     0,     0,     0,  1112,
     564,   565,   566,   567,   568,   569,   570,   571,   572,   573,
     574,   575,   576,   577,     0,     0,     0,  1126,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,  1137,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,  1139,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     0,     0,
       0,  1141,   564,   565,   566,   567,   568,   569,   570,   571,
     572,   573,   574,   575,   576,   577,     0,     0,     0,  1143,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1263,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1264,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1313,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1317,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1318,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1319,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1325,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1337,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1345,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1378,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1380,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1381,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1382,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1383,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1398,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1403,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1405,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1409,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1413,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1415,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1416,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1417,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1112,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1533,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1137,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1139,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1141,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1143,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1535,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1541,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1546,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1547,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1548,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1562,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1582,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1608,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1610,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1612,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1644,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1645,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1656,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1668,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1677,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1690,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1691,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1695,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1699,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1703,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1707,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1709,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1714,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1724,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1725,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1726,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1738,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1746,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1747,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1751,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1752,   564,   565,   566,   567,   568,   569,   570,   571,
     703,   704,   705,   706,   576,   707,     0,     0,     0,  1753,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,     0,     0,     0,  1772,   564,   565,
     566,   567,   568,   569,   570,   571,   703,   704,   705,   706,
     576,   707,     0,     0,     0,  1788,   564,   565,   566,   567,
     568,   569,   570,   571,   703,   704,   705,   706,   576,   707,
       0,     0,     0,  1794,   564,   565,   566,   567,   568,   569,
     570,   571,   703,   704,   705,   706,   576,   707,     0,     0,
       0,  1796,   550,   551,   552,   553,   554,   555,   556,   557,
     558,   559,   560,   561,     0,   562,     0,     0,     0,  1138,
     550,   551,   552,   553,   554,   555,   556,   557,   558,   559,
     560,   561,     0,   562,     0,     0,     0,  1140,   550,   551,
     552,   553,   554,   555,   556,   557,   558,   559,   560,   561,
       0,   562,     0,     0,     0,  1142,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,   578,   550,   551,   552,   553,   554,   555,   556,
     557,   558,   559,   560,   561,     0,   562,     0,     0,   563,
     564,   565,   566,   567,   568,   569,   570,   571,   703,   704,
     705,   706,   576,   707,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577
};

static const yytype_int16 yycheck[] =
{
       0,   298,   356,     8,     6,    56,    57,     8,   153,     8,
       8,    47,     8,     8,     6,   260,   100,     6,   169,   100,
       6,    72,    32,    66,   861,   379,     8,   136,    48,   190,
     173,    63,     8,   209,    67,   169,    76,    31,     8,    68,
       6,    67,    89,    67,    67,    81,    46,    51,    84,    85,
      67,    33,    61,    87,   100,   465,   290,   100,   127,    91,
     329,     8,   305,   296,    58,   100,   299,   519,   100,    69,
     313,    67,    72,    56,   290,   523,  1029,    77,   436,   502,
     523,   209,   101,   493,   103,   169,   229,   106,   169,   358,
     448,   264,   523,    33,   137,   114,   523,    97,   521,   260,
       6,   296,   102,  1056,   299,   105,   264,   119,   523,   121,
      50,   100,   112,   515,   264,    85,   523,   100,   127,   256,
     257,   258,   523,   296,   100,   125,   299,   127,   128,   129,
      77,   174,   100,   127,   100,     6,   515,     8,   296,   143,
     515,   299,   174,   153,   515,   290,   296,   176,   152,   299,
     182,   523,   184,   296,   100,   523,   299,   157,   158,   159,
     160,   161,   162,   163,   164,   165,   166,   167,   168,   169,
     170,   171,   241,   173,   174,   204,   176,   177,   178,   179,
     180,   181,   182,   183,   184,   185,   186,  1024,   127,   128,
     129,   238,   174,   182,   234,   199,   429,   137,   174,    98,
     246,   446,   238,   249,   523,   234,   249,   177,   317,   252,
     244,   523,   255,   523,   246,   523,   182,   249,   157,   523,
     159,   160,   161,   162,   163,   164,   165,   166,   167,   168,
     169,   170,   171,   467,   173,   174,   523,   176,   177,   178,
     179,   180,   181,   182,   183,   184,   185,   186,   162,   268,
     378,   467,    75,   253,   254,   255,   296,   239,   434,   299,
     249,   296,   523,   282,   299,   100,   264,   296,   288,   523,
     299,   265,   240,   293,   425,   275,   279,   100,   278,   255,
     290,   635,   523,   249,   252,   328,   296,   255,   523,   299,
     515,   425,   426,   647,   264,   314,   328,   523,   652,   342,
     513,   514,   296,   516,   523,   299,   434,   515,   355,   328,
     196,   197,   252,   229,   230,   255,   515,   516,   321,   355,
     320,   315,   467,   342,   515,   325,   326,   327,   296,   329,
     467,   299,   342,   523,   233,   335,   515,   337,   475,   174,
     340,   425,   523,   523,   425,   378,   292,    29,   351,   383,
     290,   174,   378,   345,   378,   378,   523,   360,   523,   359,
     259,   378,     8,   397,   310,   333,   409,   407,     8,   409,
     402,   372,     8,   372,   393,   274,   137,   372,     6,   361,
      95,   417,   378,   372,   366,   388,     8,   387,   195,   359,
     390,   391,   392,   436,   748,   363,   372,   426,   523,   399,
     404,   515,   516,   216,   407,   448,   372,   389,   100,   515,
     516,   462,   463,   515,   516,   304,   426,   252,   152,   422,
     255,   421,   420,   246,   467,   426,   249,   426,     8,   432,
     400,     6,   100,     0,   519,   429,   790,   244,   441,   519,
     440,   492,   493,   448,   426,   519,   519,   498,   499,   519,
     420,   519,   519,   260,   519,   519,   519,   467,  1295,   427,
     463,   296,   462,   519,   299,   467,   181,   519,   513,   514,
     515,   516,     8,   296,     8,   519,   299,   519,   519,   360,
     519,   196,   482,   483,   484,   485,   426,    31,   491,   296,
     490,   519,   299,   512,   519,     8,   525,   519,   333,     8,
     177,   519,   505,   100,   504,   519,   803,     3,   508,   509,
     510,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,   523,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   513,
     514,   515,   516,   372,   246,     8,     8,   249,   558,   559,
     560,   561,   562,     8,   564,   565,   566,   567,   568,   569,
     570,   571,   572,   573,   574,   575,   576,   577,     8,     6,
       8,   296,   523,   418,   299,   100,   521,   523,    51,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   426,   296,   523,   523,   299,   523,    51,
     523,   550,   551,   552,   553,   554,   555,   556,   557,   558,
     559,   560,   561,   562,   292,   523,   523,   523,   296,   523,
     523,   299,   523,   572,   573,   574,   575,     8,   577,   523,
     523,   523,   310,   240,   140,    87,   523,   184,   523,   523,
     523,   366,     8,   523,   523,   252,     6,   523,   255,   523,
      89,   523,   662,   663,   664,   665,   666,   667,   668,   669,
     670,   671,   672,   673,   674,   675,   523,   677,   678,   679,
     680,   681,   682,   683,   684,   685,   686,   687,   688,   689,
     523,   523,     8,   520,   525,   520,   520,     8,   333,   296,
     196,   143,   299,   703,   704,   705,   706,   707,   708,     8,
     152,     8,   100,   732,     8,     8,   158,   213,     8,     8,
     723,    66,   722,   137,     8,   240,   155,   223,   137,   729,
     730,   731,   137,   137,   137,   137,   333,   252,     8,     8,
     255,    91,   238,    46,   207,    58,   100,     8,   177,     8,
     100,     8,     8,   195,   196,   100,    59,   199,    61,   296,
      63,     8,   299,   205,     6,   426,   363,    91,   426,    72,
     770,    74,     8,    76,   793,     8,   100,     8,     8,     8,
       8,   296,     8,     8,   299,   155,    58,   100,     6,    92,
     520,   328,   137,     8,    97,     6,   420,   226,   101,   102,
     103,   264,   244,   106,     8,     8,     8,     8,   100,   372,
     113,   114,   831,   126,  1168,     8,     6,  1171,   333,    78,
       8,    80,   822,   126,   137,     8,   829,    86,   100,   174,
     427,     8,   182,   296,     8,     8,   299,     8,     8,     8,
       8,   412,   155,     8,   372,   158,   100,   150,   363,   849,
       8,     8,     6,   137,   126,     8,     8,     8,   182,   355,
       6,     8,   291,   100,   177,   137,     8,   296,   127,     8,
     299,   155,     8,     8,     8,   199,     8,   190,   137,  1233,
       8,   231,  1236,   155,  1238,   235,   158,     8,     8,   426,
       8,     8,     8,   177,    59,   154,    61,     8,   201,   249,
       8,     8,     8,     8,   249,   177,     8,   252,     8,     8,
     255,   195,   427,     8,   914,   174,     8,     8,   190,     8,
       8,     8,   246,     8,     8,   249,    91,   199,     8,     8,
     467,     8,   428,   933,   376,   100,     8,     8,     8,     8,
     100,     8,   438,   439,   440,   441,     8,     8,     8,   523,
     515,   296,   524,   953,   299,   397,   452,   453,   454,   455,
     244,   515,   404,   177,   267,     8,     6,     8,   318,   520,
       8,     8,   296,   276,   277,   299,     8,     8,   281,   979,
       8,     8,   266,   328,     8,     8,   986,   290,     6,     6,
     100,   520,   264,     8,   297,     6,   299,     6,     8,     6,
     303,     6,   315,     6,     6,   177,   290,   310,     8,   312,
       8,  1308,     8,     8,   317,   520,   512,   182,   363,   184,
       8,     8,   518,   173,    66,     8,     8,   523,     8,     8,
       8,     8,   291,     6,   337,  1035,  1036,   520,   203,     8,
       8,     8,     8,   315,  1044,  1045,     8,   360,     8,     6,
       8,   354,  1052,  1053,  1054,  1055,     8,     8,   100,   501,
       8,   345,   372,     8,     8,     8,   231,     8,   389,     8,
     235,     8,     8,    66,     8,   359,   360,   420,   420,     8,
      61,   384,     6,     8,   249,     8,   358,     8,   360,   520,
      29,   436,    31,    32,     8,   137,   100,     8,     8,     8,
       8,     8,     8,   448,     8,   389,     8,   100,     8,     8,
     520,     8,  1112,     8,     8,    54,     8,   372,     8,     8,
     520,     8,   467,     6,   372,  1125,  1126,     8,     8,     8,
       8,   296,   174,     8,   299,   300,  1136,  1137,     6,  1139,
     520,  1141,   426,  1143,   137,    84,     8,     8,   420,     8,
     434,     8,     8,     8,     8,   500,     8,   516,   516,     8,
     515,   100,   155,   328,   502,   330,   515,     8,     8,     8,
     174,   524,  1111,   520,   438,   439,   440,   441,  1178,   118,
     520,   174,     8,   520,   177,   520,  1125,  1126,   452,   453,
     454,   455,     8,   520,   520,     6,   395,     6,   501,  1138,
     520,  1140,   409,  1142,     8,   133,   371,   249,   373,     8,
     252,     8,     8,   255,   520,   520,     8,   173,     8,     8,
       8,     8,   387,     8,   389,   199,     8,   389,     6,     6,
     395,    61,     8,     8,    61,     6,   240,   520,     8,   542,
       8,   544,   545,     8,  1244,     6,     8,   373,   252,     8,
       8,   255,     8,     8,   296,     8,   249,   299,     8,   252,
       8,   426,   255,  1263,  1264,     8,     8,     8,   261,     6,
     435,   395,     8,     8,     8,     8,     8,     8,     8,   218,
     133,   520,     8,     8,     8,     8,   328,     6,     6,    61,
     373,     8,   296,   389,     8,   299,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
     249,   389,     8,  1313,   389,   515,     8,  1317,  1318,  1319,
    1320,     6,   520,   262,   328,  1325,   520,     8,   520,   333,
     395,     8,   271,     8,     8,   328,   639,  1337,   641,     8,
       8,  1341,     8,     8,     8,  1345,     8,   373,   395,   394,
     653,   395,     8,   389,   373,     8,     8,   296,   297,   298,
     299,   300,     8,     8,   303,     8,   359,     8,     8,   362,
       8,     8,   373,   394,   373,     8,   373,     8,  1378,     8,
    1380,  1381,  1382,  1383,   373,   389,     8,     8,    61,    61,
     394,   395,  1039,    64,   436,  1026,    -1,    64,  1398,   859,
      61,    -1,   303,  1403,    -1,  1405,   448,    -1,    -1,  1409,
      -1,    -1,    -1,  1413,   717,  1415,  1416,  1417,    -1,    -1,
      -1,    -1,    -1,   427,    -1,   467,    -1,    -1,    -1,    -1,
      -1,    -1,   735,   372,    -1,    -1,    -1,   430,    -1,   378,
      -1,   445,    -1,   436,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   394,   448,   109,   110,   111,   112,
      -1,   114,   115,   116,    -1,  1465,   119,   120,   121,   122,
     123,   124,    -1,    -1,   467,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   789,   426,   791,   792,
      -1,   794,    -1,   796,  1494,   798,    -1,    -1,    -1,    -1,
      -1,   804,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,  1511,    -1,    -1,   167,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   825,    -1,   827,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,  1533,    -1,  1535,    -1,    -1,    -1,   192,
      -1,  1541,    -1,    -1,    -1,    -1,  1546,  1547,  1548,    -1,
      -1,    -1,    -1,  1553,    -1,   208,   859,    -1,    -1,    -1,
      -1,    -1,  1562,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,  1582,    -1,    -1,    -1,   889,    -1,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,  1608,   524,
    1610,    -1,  1612,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,  1644,  1645,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1656,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,  1668,   516,
      -1,    -1,    -1,    -1,   327,    -1,    -1,  1677,  1678,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   989,    -1,    -1,    -1,
    1690,  1691,    -1,    -1,    -1,  1695,    -1,    -1,    -1,  1699,
      -1,    -1,    -1,  1703,    -1,    -1,    -1,  1707,    -1,  1709,
      -1,    -1,    -1,    -1,  1714,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,  1724,  1725,  1726,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,  1039,    -1,  1738,    -1,
    1043,    -1,    -1,  1046,    -1,    -1,  1746,  1747,    -1,  1052,
    1053,  1751,  1752,  1753,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,  1771,  1772,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,  1788,   516,
      -1,    -1,    -1,   520,  1794,    -1,  1796,     0,     1,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      -1,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    -1,    -1,    -1,    51,    -1,
      -1,    -1,    -1,    56,    57,  1158,  1159,  1160,  1161,  1162,
    1163,    -1,  1165,  1166,    -1,    68,    -1,    -1,    -1,    72,
      -1,    -1,    -1,    -1,    77,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    87,    -1,    -1,    -1,    -1,    92,
      -1,    94,    -1,    -1,    97,    -1,    99,   100,   101,    -1,
      -1,    -1,    -1,    -1,    -1,   108,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,  1224,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     143,   516,    -1,  1246,  1247,   520,   149,    -1,    -1,   152,
      -1,    -1,    -1,    -1,    -1,   158,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   166,    -1,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   180,   516,    -1,
      -1,    -1,   520,    -1,    -1,   188,   189,    -1,    -1,    -1,
      -1,   194,    -1,   196,   197,    -1,   199,    -1,  1301,    -1,
      -1,    -1,   205,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     213,    -1,    -1,    -1,    -1,    -1,   219,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   229,   230,    -1,   232,
      -1,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   244,   516,    -1,    -1,    -1,    -1,   250,    -1,    -1,
      -1,   254,    -1,    -1,    -1,    -1,    -1,    -1,   261,    -1,
     263,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,     3,   516,    -1,    -1,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,   304,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   326,    -1,    -1,    -1,    -1,   331,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   339,    -1,   341,   342,
      -1,   344,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   352,
     353,   354,    -1,    -1,    -1,   358,    -1,    -1,    -1,   362,
      -1,    -1,    -1,   366,    -1,    -1,   369,    -1,    -1,    -1,
      -1,  1474,   375,   376,   377,    -1,    -1,    -1,  1481,    -1,
     383,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   394,   395,    -1,   397,    -1,   399,    -1,    -1,    -1,
     403,   404,    -1,    -1,    -1,    -1,    -1,   410,    -1,    -1,
     413,   414,   415,    -1,    -1,    -1,    -1,    -1,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,   436,   437,   438,   439,   440,   441,   524,
      -1,    -1,    -1,    -1,   447,   448,   449,    -1,   451,   452,
     453,   454,   455,    -1,    -1,    -1,    -1,    -1,    -1,   462,
     463,    -1,    -1,    -1,   196,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,   213,    -1,    -1,    -1,    -1,   524,    -1,    -1,   492,
     493,    -1,    -1,    -1,    -1,   498,   499,   500,   501,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,   519,     3,    -1,    -1,
     523,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
      -1,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    91,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   100,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,   117,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    91,    -1,
      -1,    -1,    -1,    -1,     3,    -1,    -1,   100,     7,     8,
       9,    10,    11,    12,    13,    14,    15,    16,    17,    18,
      19,    20,    21,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    -1,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,   182,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,
     196,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     452,   453,   454,   455,    -1,    -1,    -1,   213,    -1,    -1,
      -1,    -1,    -1,   465,   466,    -1,    -1,    -1,    -1,   182,
      -1,   184,    -1,    -1,    -1,   231,    -1,    -1,    -1,    -1,
      -1,   100,    -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   249,   496,   497,    -1,    -1,   117,   255,
     213,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     512,    -1,    -1,    -1,    -1,    -1,   518,    -1,   231,    -1,
      -1,   523,    -1,    -1,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   249,   516,    -1,    -1,
     296,    -1,   255,   299,    -1,   524,   109,   110,   111,   112,
      -1,   114,   115,   116,    -1,   174,   119,   120,   121,   122,
     123,   124,    -1,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,   196,    -1,    -1,
      -1,    -1,    -1,   296,   524,    -1,   299,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   213,    -1,    -1,    -1,    -1,    -1,
      -1,   357,    -1,    -1,   167,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   246,    -1,   192,
     249,    -1,    -1,    -1,    -1,    -1,   255,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   357,   208,   402,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   389,   296,    -1,    -1,
     299,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,   275,   276,   277,   278,   279,   280,   281,   282,
     283,   284,    -1,    -1,    -1,   438,   439,   440,   441,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,
     453,   454,   455,    -1,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   512,   516,    -1,    -1,
      -1,    -1,   518,    -1,   327,   524,    -1,   523,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,    -1,
     523,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   438,
     439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   452,   453,   454,   455,     3,    -1,    -1,
      -1,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   512,    -1,    -1,    -1,    63,    -1,   518,
     524,    -1,    -1,    -1,   523,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   100,    -1,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,     3,   516,
      -1,   117,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    -1,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
       3,    -1,    -1,    -1,     7,     8,     9,    10,    11,    12,
      13,    14,    15,    16,    17,    18,    19,    20,    21,    22,
      23,    24,    25,    26,    27,    28,    29,    30,    31,    32,
     196,    34,    35,    36,    37,    38,    39,    40,    41,    42,
      43,    44,    -1,    -1,    -1,    -1,    -1,   213,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     246,    -1,    -1,   249,    -1,    -1,    -1,    -1,    -1,   255,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   100,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     296,    -1,    -1,   299,   189,    -1,    -1,    -1,    -1,    -1,
      -1,   196,    -1,    -1,    -1,    -1,    -1,    -1,   203,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   213,    -1,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,   234,
     235,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   196,    -1,    -1,   199,    -1,    -1,    -1,
     255,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     213,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   296,    -1,   246,   299,    -1,   249,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
      -1,    -1,    -1,   296,    -1,    -1,   299,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,   372,   524,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,
      -1,    -1,   518,    -1,    -1,    -1,    -1,   523,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   372,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   438,   439,   440,   441,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,
     455,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,
     453,   454,   455,    -1,    -1,    -1,    -1,   512,    -1,    -1,
      -1,    -1,    -1,   518,    -1,     3,    -1,    -1,   523,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    -1,    -1,   512,
      -1,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,     3,
     523,    -1,    -1,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    -1,
      34,    35,    36,    37,    38,    39,    40,    41,    42,    43,
      44,    -1,   100,    -1,    -1,    -1,     3,    -1,    -1,    -1,
       7,     8,     9,    10,    11,    12,    13,    14,    15,    16,
      17,    18,    19,    20,    21,    22,    23,    24,    25,    26,
      27,    28,    29,    30,    31,    32,    -1,    34,    35,    36,
      37,    38,    39,    40,    41,    42,    43,    44,     3,    -1,
      -1,    -1,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    -1,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
      -1,   199,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   213,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   246,    -1,
      -1,   249,   196,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,   213,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,   296,   196,
      -1,   299,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,   213,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   196,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   213,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,     3,    -1,    -1,   372,     7,     8,     9,    10,    11,
      12,    13,    14,    15,    16,    17,    18,    19,    20,    21,
      22,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    -1,    34,    35,    36,    37,    38,    39,    40,    41,
      42,    43,    44,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,
     438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   452,   453,   454,   455,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   438,   439,   440,   441,    -1,    -1,
      -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,
     454,   455,    -1,    -1,   512,    -1,    -1,    -1,    -1,    -1,
     518,   465,   466,    -1,    -1,   523,    -1,    -1,    -1,   426,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   496,   497,    -1,   452,   453,   454,   455,    -1,
      -1,    -1,    -1,    -1,   196,    -1,    -1,    -1,   512,    -1,
      -1,   426,    -1,    -1,   518,    -1,    -1,    -1,    -1,   523,
      -1,   213,    -1,   438,   439,   440,   441,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,
     455,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,
      -1,   518,   524,    -1,    -1,    -1,   523,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,
      -1,    -1,    -1,   518,    -1,     3,    -1,    -1,   523,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,     3,    -1,    -1,
      -1,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   438,   439,   440,   441,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     452,   453,   454,   455,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   196,    -1,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,   213,    -1,    -1,    -1,    -1,
     512,   524,    -1,    -1,    -1,    -1,   518,    -1,    -1,    -1,
      -1,   523,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     196,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,   213,    -1,    -1,
      -1,    -1,   524,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     524,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   524,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,   515,   516,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,   507,
     508,   509,   510,   511,   512,   513,   514,   515,   516,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
     438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   452,   453,   454,   455,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,   455,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,    -1,
     518,   524,    -1,    -1,    -1,   523,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,
      -1,    -1,   518,    -1,    -1,    -1,    -1,   523,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   524,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   524,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,    -1,   516,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   524,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,    -1,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,   522,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,   522,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,   515,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,   520,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,    -1,    -1,
      -1,   520,   503,   504,   505,   506,   507,   508,   509,   510,
     511,   512,   513,   514,    -1,   516,    -1,    -1,    -1,   520,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,    -1,   516,    -1,    -1,    -1,   520,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
      -1,   516,    -1,    -1,    -1,   520,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,   519,   503,   504,   505,   506,   507,   508,   509,
     510,   511,   512,   513,   514,    -1,   516,    -1,    -1,   519,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   515,   516,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516
};

  /* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
     symbol of state STATE-NUM.  */
static const yytype_int16 yystos[] =
{
       0,     1,     3,     7,     8,     9,    10,    11,    12,    13,
      14,    15,    16,    17,    18,    19,    20,    21,    22,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      51,    56,    57,    68,    72,    77,    87,    92,    94,    97,
      99,   100,   101,   108,   143,   149,   152,   158,   166,   180,
     188,   189,   194,   196,   197,   199,   205,   213,   219,   229,
     230,   232,   244,   250,   254,   261,   263,   304,   326,   331,
     339,   341,   342,   344,   352,   353,   354,   358,   362,   366,
     369,   375,   376,   377,   383,   394,   395,   397,   399,   403,
     404,   410,   413,   414,   415,   436,   437,   438,   439,   440,
     441,   447,   448,   449,   451,   452,   453,   454,   455,   462,
     463,   492,   493,   498,   499,   500,   501,   512,   518,   523,
     527,   528,   529,   530,   533,   545,   547,   551,   553,   556,
     558,   560,   561,   562,   563,   564,   565,   566,   569,   570,
     571,   572,   594,   595,   596,   597,   519,   502,   521,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,     6,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,    58,   100,   126,
     137,   155,   158,   177,   190,   315,   360,   539,    91,   100,
     182,   199,   246,   249,   296,   299,   573,   588,   329,   358,
       3,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    34,    35,    36,    37,    38,
      39,    40,    41,    42,    43,    44,   100,   117,   174,   196,
     246,   249,   255,   512,   518,   523,   588,   597,     6,    87,
     244,   383,   397,    29,     8,     8,     8,   137,   552,     6,
     199,   264,   358,   420,   539,    33,    50,   137,   252,   255,
     290,   426,   546,    95,   181,   196,   366,   588,    75,   100,
     174,   246,   249,   426,   588,     8,    68,   176,   204,   234,
     426,   525,   588,   571,   572,     8,    78,    80,    86,   127,
     137,   154,   174,   291,   548,    66,   100,   137,   174,   249,
     252,   255,   328,   342,   409,   436,   448,   467,   557,   523,
     216,     8,   196,   597,     8,    33,   174,   239,   361,   366,
     389,   426,   554,   196,   197,    63,   100,   117,   196,   246,
     249,   255,   588,   597,   304,   588,   532,   100,   240,   252,
     255,   333,   363,   427,   568,   588,   140,   223,   238,   355,
     428,   597,    51,   143,   152,   199,   404,     8,    77,   531,
     152,     8,    66,   100,   137,   174,   249,   252,   255,   328,
     436,   448,   467,   567,   588,   136,   317,   580,     8,   448,
      91,   100,   117,   182,   196,   231,   249,   255,   357,   402,
     588,   597,     6,   100,   182,   249,   372,     6,   539,   100,
     174,   252,   255,   418,   568,   588,    91,   100,   182,   184,
     231,   249,   255,   357,   389,   588,   597,   100,   174,   240,
     252,   255,   328,   389,   394,   395,   427,   445,   568,   588,
       6,   100,   182,   249,   372,   465,   466,   496,   497,   597,
     100,   240,   252,   255,   363,   427,   568,   588,    51,    87,
     143,   152,   158,   195,   196,   199,   205,   244,   376,   397,
     404,   501,   465,   466,   496,   497,   597,   588,   100,   240,
     252,   255,   363,   427,   568,   588,   573,   573,   597,    66,
     100,   137,   174,   249,   252,   255,   328,   363,   436,   448,
     467,   500,   559,   588,     3,    29,   596,   597,   596,   597,
     596,   597,     0,   519,   519,   519,   519,   519,   519,   519,
     519,   519,   519,   519,   519,   519,   519,   519,   519,   519,
     519,   519,    76,   234,   407,   409,   574,   588,   519,   519,
     503,   504,   505,   506,   507,   508,   509,   510,   511,   512,
     513,   514,   516,   519,   503,   504,   505,   506,   507,   508,
     509,   510,   511,   512,   513,   514,   515,   516,   519,   596,
     597,   597,   596,   597,   596,   597,   596,   597,   596,   597,
     596,   597,   596,   597,   596,   597,   596,   597,   596,   597,
     596,   597,   596,   597,   596,   597,   596,   597,     8,   596,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   596,   597,   542,     8,   360,    89,   155,
     177,   226,   291,   296,   299,   587,    31,   540,     8,   190,
     260,     8,   543,   177,   372,     8,     8,   588,     8,     8,
       6,     8,   521,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
     523,   523,   523,   523,   523,   523,   523,   523,   523,   523,
       8,   100,   292,   310,   591,     8,     8,   436,   448,   589,
     597,   597,   597,   511,   512,   513,   514,   516,   520,   520,
     520,   520,   588,   568,     8,    66,   100,   137,   155,   174,
     177,   249,   252,   255,   261,   328,   359,   362,   430,   436,
     448,   467,   538,     8,     8,   137,   155,   177,   195,   244,
     266,   290,   345,   359,   360,   389,   426,   434,   534,   597,
     588,   588,   597,   589,     8,   588,   568,   100,     8,   588,
       8,     8,     8,   169,   425,   426,   590,   184,   328,   426,
     467,   588,   590,   588,   256,   257,   258,   467,   475,   582,
     438,   439,   440,   441,   452,   453,   454,   455,   574,   137,
     550,   137,   137,   549,   137,   588,   137,   588,   137,   568,
       8,     8,    51,   207,   264,   588,   100,     8,   597,   589,
       8,     8,     8,   597,   597,   597,   568,   597,     6,   100,
     246,   249,   520,   597,   426,   189,   203,   234,   235,   255,
     372,   555,   588,   597,   173,   229,   296,   299,   426,   597,
       8,   568,     8,   372,   426,     8,     8,     8,   589,   520,
     588,   553,     8,     8,   597,   589,    63,    91,   100,   174,
     182,   184,   246,   249,   328,   402,     8,    98,   233,   259,
     274,   593,     8,   155,     6,     8,    85,   177,   264,   359,
     400,   420,     8,   520,     6,   553,   420,     8,     8,   264,
     588,   100,     8,   597,   589,     8,   597,   597,   597,   568,
     597,   372,     8,     6,     8,     8,     8,   589,     8,     8,
     100,   174,   255,   372,   520,     8,     8,     8,     8,     8,
     100,   597,   589,   412,   372,     8,     8,     6,     8,     8,
     589,     8,     6,   520,     8,   100,     8,   597,   589,     8,
       8,     8,     8,     8,     8,     8,     8,     8,     8,     8,
       8,     8,     8,   520,     8,     8,   597,   589,     8,   593,
       8,     8,     8,     8,   264,   420,     8,     8,     8,     8,
       8,     8,     8,     8,     8,   597,   597,   597,   597,   520,
       8,     8,   597,   589,     8,   593,   520,     8,     8,   264,
     588,   100,     8,   597,   589,     8,     8,   597,   597,   597,
       8,   568,   524,   524,   100,   246,   249,   579,   588,     6,
      91,   100,   182,   231,   235,   249,   318,   578,     8,    59,
      61,    91,   100,   182,   184,   203,   231,   235,   249,   300,
     328,   330,   371,   373,   387,   389,   395,   426,   435,   576,
     577,   588,     8,    29,    31,    32,    54,    84,   100,   118,
     218,   249,   262,   271,   297,   298,   300,   303,   372,   378,
     426,   575,   583,   588,   596,   597,   596,   596,   596,   596,
     596,   596,   596,   596,   597,   596,   597,   596,   597,   596,
     597,   596,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   596,   597,   596,   597,   596,   597,   596,   597,   597,
     596,   597,   522,   524,   524,   524,   524,   524,   524,   524,
     524,   520,   520,   524,   524,   524,   524,   524,   524,   524,
     524,   524,   524,   524,   524,   520,   520,   524,   524,   524,
     524,   524,   524,   524,   524,   524,   520,   520,   520,   520,
     520,   520,   520,   520,   524,   524,   524,   524,   524,   524,
     524,   524,   524,   524,   524,   524,   553,   177,    31,    58,
     127,   265,   315,   429,   588,   260,   446,   588,   541,   553,
       8,   544,   553,    32,   153,   290,   342,   426,   467,   588,
       8,   520,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,     8,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
       8,     8,   597,   597,   597,   597,   597,   597,     8,     8,
       8,     6,     8,     8,   264,   588,     6,   100,     6,     8,
     597,   589,     6,   537,     8,     6,   536,     6,   535,     6,
     597,   597,   597,   568,   520,   195,   244,   260,   588,     6,
       6,     8,     8,     8,     8,     8,   173,     6,     8,     8,
       8,   177,   553,   520,   520,     8,    67,   378,   109,   110,
     111,   112,   114,   115,   116,   119,   120,   121,   122,   123,
     124,   167,   192,   208,   275,   276,   277,   278,   279,   280,
     281,   282,   283,   284,   327,   585,     8,     8,   597,   588,
     553,   429,   588,   588,   568,   588,   588,   588,     8,   100,
     590,   588,     8,   520,   520,     6,   520,   520,   520,   520,
     524,     8,     8,     8,   597,   520,     8,   100,   588,     8,
     588,   209,   434,   586,   589,     8,   568,   520,     6,     8,
       8,   394,     8,     8,   597,   520,     8,   372,   426,   372,
       8,   100,   310,   588,   591,     8,   585,     8,     8,     8,
       8,   372,     8,   389,     6,   345,     8,   420,   290,   467,
     420,     8,    61,     6,     8,     8,   588,     8,   520,   520,
     520,   520,   520,   520,     8,     8,     8,    47,    81,    84,
      85,   238,   355,   417,   592,     8,   597,     8,   520,     8,
       8,   597,     8,   520,   597,   520,     8,     8,   597,   520,
     597,   588,     8,   520,   520,   520,   520,   520,     8,     8,
       8,   372,     8,     8,    89,   238,   355,   581,   305,   313,
       8,     8,    67,   378,   520,     8,     6,   372,     8,     8,
       8,   585,     8,   581,   209,   378,   434,     8,   592,     8,
       6,    48,   288,   293,   584,     8,     8,   426,   597,   426,
     597,    67,   378,     8,   577,   520,   588,   597,   597,   588,
       8,     8,     8,   100,   199,   246,   249,   372,   588,   597,
     100,   199,   246,   249,   372,   588,   597,   597,   597,   592,
       8,     8,    67,   378,   502,   596,   597,   596,   597,   596,
     597,   597,   597,   596,   597,   596,   597,   596,   597,   153,
     290,   467,   588,   588,   588,   588,   588,   588,   588,   588,
     553,   553,     8,     8,     8,     6,   467,   290,   467,   597,
       8,   522,   524,   520,   524,   520,   520,   520,   520,   588,
       8,   520,   553,   520,   553,   553,   520,   520,   520,   597,
       8,   588,   588,   520,   520,     6,   395,     6,   597,   597,
     585,   520,   520,   588,   133,   100,   590,     8,   597,     8,
       8,   597,   597,   597,   597,   597,     8,   520,   520,   597,
     173,   597,   520,   597,     8,     8,     8,     8,     8,     8,
     199,     8,   389,     6,     6,    61,     8,     8,    61,     6,
     597,     8,   597,   597,   597,   597,   597,   597,   520,   597,
     520,   597,   520,   597,     8,   597,   597,   597,     8,   520,
       6,     8,    67,   378,    67,   378,   597,     8,   588,     8,
       8,     8,     8,   588,     8,     8,     8,   597,   524,   524,
     524,   524,   524,   524,   520,   520,   524,   524,   524,   524,
     524,   524,     8,     8,   597,     6,   520,   373,   597,   597,
       8,     8,     8,   597,     8,   597,   597,   597,   520,   597,
       8,     8,     8,   597,   133,     8,     8,   520,   162,     8,
       8,     6,   597,     6,   395,     8,   389,   389,   520,   389,
     520,   520,   597,   597,   597,   520,     8,   597,   597,   520,
     597,     8,   524,   520,   520,   520,   520,   520,   597,   520,
     373,     8,   597,   597,   520,     6,     8,   395,     8,     8,
       8,     8,   597,   597,   520,   520,   520,   597,   524,   524,
     597,   597,     8,     8,     8,   597,   597,     8,   520,   597,
     373,     8,   395,   394,   389,   395,   520,   520,   597,   597,
     597,   520,   520,   520,   597,     8,   373,     8,     8,     8,
       8,   597,   597,   597,   597,   597,     8,   373,   373,   394,
     373,   409,   520,   127,   241,     8,     8,     8,     8,   597,
     597,     8,     8,    61,    61,   373,    61,   127,   520,    61,
       8,     8,   597,    61,   520,   597,   520,   597
};

  /* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_int16 yyr1[] =
{
       0,   526,   527,   527,   527,   527,   527,   527,   527,   527,
     527,   527,   527,   527,   527,   527,   527,   527,   527,   527,
     527,   527,   527,   527,   527,   527,   527,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   528,   528,   528,   528,   528,
     528,   528,   528,   528,   528,   529,   529,   529,   529,   529,
     529,   529,   529,   529,   529,   529,   529,   530,   530,   530,
     530,   530,   530,   530,   530,   530,   530,   530,   530,   530,
     530,   530,   530,   530,   530,   530,   530,   530,   530,   530,
     530,   530,   531,   530,   532,   530,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   534,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   535,   533,   536,   533,   537,   533,
     538,   533,   533,   533,   533,   533,   533,   533,   533,   533,
     533,   533,   533,   533,   539,   539,   539,   539,   539,   539,
     539,   539,   539,   539,   539,   539,   540,   539,   541,   539,
     542,   539,   543,   539,   544,   539,   539,   539,   539,   539,
     539,   539,   539,   539,   539,   539,   539,   539,   539,   539,
     539,   545,   545,   545,   545,   546,   545,   545,   545,   545,
     545,   547,   547,   547,   547,   547,   547,   547,   548,   547,
     549,   547,   547,   547,   547,   547,   550,   547,   547,   547,
     547,   551,   551,   552,   551,   551,   551,   553,   553,   553,
     553,   553,   553,   553,   553,   553,   553,   553,   553,   553,
     553,   553,   553,   553,   553,   553,   553,   553,   553,   554,
     553,   555,   553,   553,   553,   556,   557,   556,   556,   556,
     556,   556,   556,   556,   556,   556,   556,   556,   556,   556,
     556,   556,   556,   556,   556,   556,   558,   559,   558,   558,
     558,   558,   558,   558,   558,   558,   558,   558,   558,   558,
     558,   558,   558,   558,   560,   560,   560,   560,   560,   560,
     560,   560,   561,   561,   561,   561,   561,   561,   561,   561,
     562,   562,   562,   562,   562,   562,   562,   562,   563,   563,
     563,   563,   563,   563,   563,   564,   564,   564,   564,   564,
     564,   564,   564,   564,   564,   564,   565,   565,   565,   565,
     565,   565,   565,   565,   565,   565,   565,   565,   565,   566,
     567,   566,   566,   566,   566,   566,   566,   566,   566,   566,
     566,   566,   566,   566,   566,   568,   568,   568,   568,   568,
     568,   568,   568,   568,   568,   568,   568,   568,   568,   568,
     568,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   569,   569,   569,
     569,   569,   569,   569,   569,   569,   569,   570,   570,   570,
     570,   570,   571,   571,   571,   571,   571,   571,   572,   572,
     572,   573,   573,   573,   573,   573,   573,   573,   574,   574,
     574,   574,   574,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   575,   575,   575,   575,   575,   575,
     575,   575,   575,   575,   576,   576,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   577,   577,   577,   577,
     577,   577,   577,   577,   577,   577,   578,   578,   578,   578,
     578,   578,   578,   578,   578,   578,   578,   579,   579,   579,
     579,   580,   580,   581,   581,   581,   582,   582,   582,   582,
     582,   583,   583,   583,   584,   584,   584,   585,   585,   585,
     585,   585,   585,   585,   585,   585,   585,   585,   585,   585,
     585,   585,   585,   585,   585,   585,   585,   585,   585,   585,
     585,   585,   585,   585,   586,   586,   587,   587,   587,   587,
     588,   588,   589,   589,   590,   590,   591,   591,   591,   592,
     592,   592,   592,   592,   592,   592,   593,   593,   593,   593,
     594,   595,   595,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   596,   596,   596,   596,   596,
     596,   596,   596,   596,   596,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597,   597,   597,
     597,   597,   597,   597,   597,   597,   597,   597
};

  /* YYR2[YYN] -- Number of symbols on the right hand side of rule YYN.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     0,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     0,     2,     2,
       3,     2,     2,     8,     3,     3,     3,     3,     3,     4,
       4,     2,     2,     3,     2,     2,     2,     8,     3,     3,
       3,     3,     3,     4,     4,     2,     2,     2,     3,     2,
       2,     4,     3,     3,     3,     3,     3,     3,     3,     4,
       4,     4,     4,     4,     3,     1,     1,     1,     1,     2,
       2,     1,     4,     5,     7,     3,     1,     2,     3,     2,
       1,     1,     2,     1,     8,     1,     1,     2,     2,     2,
       2,     2,     3,     3,     2,     2,     3,     1,     2,     8,
       8,     8,     0,     3,     0,     3,     3,     2,     2,     3,
       2,     7,     3,     2,     3,     4,     5,    12,    13,    12,
      13,    11,    12,    11,    12,     4,     4,     0,     4,     4,
       5,     4,     4,     4,     4,     5,     5,     4,     5,     4,
       8,     4,     4,     7,     6,     5,    12,     3,     4,     4,
       5,     4,     5,     4,     4,     4,     4,     4,     4,     4,
      11,    12,    13,    14,     0,     5,     0,     5,     0,     5,
       0,     4,     6,     4,     4,     5,     4,     4,     4,     5,
      10,     6,     6,     6,     2,     3,     4,     4,     4,     4,
       4,     4,     4,     4,     3,     2,     0,     3,     0,     4,
       0,     3,     0,     3,     0,     4,     3,     2,     2,     3,
       4,     4,     5,     4,     4,     4,     4,     6,     5,     5,
       7,     3,     3,     3,     3,     0,     3,     3,     3,     5,
       5,     2,     3,     3,     4,     6,     7,     4,     0,     3,
       0,     4,     4,     4,     4,     3,     0,     4,     2,     1,
       4,     3,     3,     0,     3,     3,     9,     3,     4,     4,
       4,     5,     2,     3,     4,     4,     3,     3,     4,     6,
       4,     5,     5,     4,     4,     5,     4,     4,     4,     0,
       3,     0,     4,     6,     6,     3,     0,     3,     5,     3,
       5,     3,     4,     3,     5,     6,     3,     3,     4,     4,
       4,     5,     9,     5,     5,     5,     3,     0,     3,     5,
       3,     3,     4,     2,     3,     3,     3,     3,     3,     4,
       9,     5,     5,     5,     2,     3,     3,     3,     3,     5,
       3,     2,     2,     3,     3,     3,     3,     5,     3,     2,
       2,     3,     3,     3,     3,     5,     3,     2,     2,     3,
       4,     4,     3,     5,     2,     2,     3,     4,     3,     3,
       3,     3,     3,     3,     4,     3,     2,     3,     3,     3,
       3,     3,     3,     3,     3,     4,     3,     5,     2,     3,
       0,     3,     5,     3,     3,     4,     2,     3,     3,     3,
       4,     9,     5,     5,     5,     3,     3,     3,     3,     3,
       3,     4,     3,     4,     3,     3,     4,     4,     3,     4,
       4,     2,     3,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     1,     1,     1,     2,    17,     2,     8,     3,
       3,     3,     3,     8,     3,     3,     3,     3,     2,     3,
       3,     3,     3,     2,     3,     3,     3,     3,     2,     3,
       3,     3,     3,     3,     4,     2,     3,     4,     4,     3,
       3,     3,     3,     5,     6,     6,     4,     2,     1,     2,
       3,     2,     1,     1,     1,     1,     1,     1,     2,     2,
       2,     1,     2,     2,     2,     2,     3,     2,     2,     2,
       2,     2,     1,     1,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     1,     2,     2,     3,     3,     2,
       2,     3,     3,     3,     3,     3,     3,     3,     3,     2,
       2,     2,     2,     3,     1,     2,     1,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     3,     3,     3,     3,     2,
       2,     3,     2,     2,     2,     3,     1,     2,     2,     2,
       2,     4,     2,     3,     2,     2,     2,     1,     2,     2,
       2,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       6,     3,     3,     1,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3,     4,     4,     4,     4,     4,
       6,     4,     4,     1,     1,     1,     4,     4,     4,     4,
       6,     6,     6,     6,     1,     1,     4,     4,     4,     4,
       4,     8,     6,     6,     6,     4,     4,     1,     1,     1,
       4,     4,     4,     4,     3,     3,     3,     3,     3,     3,
       3,     3,     2,     3,     2,     1,     1,     4,     3,     3,
       3,     3,     3,     3,     4,     4,     4,     4,     6,     4,
       4,     1,     1,     1,     4,     4,     4,     4,     6,     3,
       3,     3,     3,     3,     3,     3,     3,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     4,     4,     4,
       4,     4,     8,     6,     6,     6,     4,     4,     1,     1,
       1,     4,     4,     4,     4,     5,     7,     3,     3,     3,
       3,     3,     3,     3,     3,     2,     3,     2
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)

/* This macro is provided for backward compatibility. */
# ifndef YY_LOCATION_PRINT
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif


# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YYUSE (yyoutput);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yykind < YYNTOKENS)
    YYPRINT (yyo, yytoknum[yykind], *yyvaluep);
# endif
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YYUSE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YYUSE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    goto yyexhaustedlab;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          goto yyexhaustedlab;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */

  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 3: /* list: asgn '\n'  */
#line 669 "gram.y"
                    {}
#line 5177 "y.tab.c"
    break;

  case 4: /* list: vasgn '\n'  */
#line 670 "gram.y"
                     {}
#line 5183 "y.tab.c"
    break;

  case 5: /* list: expr '\n'  */
#line 671 "gram.y"
                    { result = (yyvsp[-1].val); }
#line 5189 "y.tab.c"
    break;

  case 6: /* list: vexpr '\n'  */
#line 672 "gram.y"
                     { result = *(yyvsp[-1].ptr); }
#line 5195 "y.tab.c"
    break;

  case 18: /* list: isolines '\n'  */
#line 684 "gram.y"
                        { }
#line 5201 "y.tab.c"
    break;

  case 26: /* list: error '\n'  */
#line 692 "gram.y"
                     { return 1; }
#line 5207 "y.tab.c"
    break;

  case 28: /* annotation: CLEAR BOX  */
#line 696 "gram.y"
                    { do_clear_boxes(); }
#line 5213 "y.tab.c"
    break;

  case 29: /* annotation: WITH BOX  */
#line 697 "gram.y"
                   { curbox = next_box(); }
#line 5219 "y.tab.c"
    break;

  case 30: /* annotation: WITH BOX NUMBER  */
#line 698 "gram.y"
                          { curbox = (int) (yyvsp[0].val); }
#line 5225 "y.tab.c"
    break;

  case 31: /* annotation: BOX onoff  */
#line 699 "gram.y"
                    { boxes[curbox].active = (yyvsp[0].pset); }
#line 5231 "y.tab.c"
    break;

  case 32: /* annotation: BOX GRAPHNO  */
#line 700 "gram.y"
                      { boxes[curbox].gno = (yyvsp[0].pset); }
#line 5237 "y.tab.c"
    break;

  case 33: /* annotation: BOX expr ',' expr ',' expr ',' expr  */
#line 702 "gram.y"
        {
	    if (curbox >= 0 && curbox < maxboxes) {
		boxes[curbox].x1 = (yyvsp[-6].val);
		boxes[curbox].y1 = (yyvsp[-4].val);
		boxes[curbox].x2 = (yyvsp[-2].val);
		boxes[curbox].y2 = (yyvsp[0].val);
	    }
	}
#line 5250 "y.tab.c"
    break;

  case 34: /* annotation: BOX LOCTYPE worldview  */
#line 710 "gram.y"
                                { sysbox.loctype = (yyvsp[0].pset); }
#line 5256 "y.tab.c"
    break;

  case 35: /* annotation: BOX LINESTYLE NUMBER  */
#line 711 "gram.y"
                               { sysbox.lines = (int) (yyvsp[0].val); }
#line 5262 "y.tab.c"
    break;

  case 36: /* annotation: BOX LINEWIDTH NUMBER  */
#line 712 "gram.y"
                               { sysbox.linew = (int) (yyvsp[0].val); }
#line 5268 "y.tab.c"
    break;

  case 37: /* annotation: BOX COLOR NUMBER  */
#line 713 "gram.y"
                           { sysbox.color = (int) (yyvsp[0].val); }
#line 5274 "y.tab.c"
    break;

  case 38: /* annotation: BOX FILL filltype  */
#line 714 "gram.y"
                            { sysbox.fill = (yyvsp[0].pset); }
#line 5280 "y.tab.c"
    break;

  case 39: /* annotation: BOX FILL COLOR NUMBER  */
#line 715 "gram.y"
                                { sysbox.fillcolor = (int) (yyvsp[0].val); }
#line 5286 "y.tab.c"
    break;

  case 40: /* annotation: BOX FILL PATTERN NUMBER  */
#line 716 "gram.y"
                                  { sysbox.fillpattern = (int) (yyvsp[0].val); }
#line 5292 "y.tab.c"
    break;

  case 41: /* annotation: BOX DEF  */
#line 718 "gram.y"
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
#line 5308 "y.tab.c"
    break;

  case 42: /* annotation: WITH LINE  */
#line 729 "gram.y"
                    { curline = next_line(); }
#line 5314 "y.tab.c"
    break;

  case 43: /* annotation: WITH LINE NUMBER  */
#line 730 "gram.y"
                           { curline = (int) (yyvsp[0].val); }
#line 5320 "y.tab.c"
    break;

  case 44: /* annotation: CLEAR LINE  */
#line 731 "gram.y"
                     { do_clear_lines(); }
#line 5326 "y.tab.c"
    break;

  case 45: /* annotation: LINE onoff  */
#line 732 "gram.y"
                     { lines[curline].active = (yyvsp[0].pset); }
#line 5332 "y.tab.c"
    break;

  case 46: /* annotation: LINE GRAPHNO  */
#line 733 "gram.y"
                       { lines[curline].gno = (yyvsp[0].pset); }
#line 5338 "y.tab.c"
    break;

  case 47: /* annotation: LINE expr ',' expr ',' expr ',' expr  */
#line 735 "gram.y"
        {
	    lines[curline].x1 = (yyvsp[-6].val);
	    lines[curline].y1 = (yyvsp[-4].val);
	    lines[curline].x2 = (yyvsp[-2].val);
	    lines[curline].y2 = (yyvsp[0].val);
	}
#line 5349 "y.tab.c"
    break;

  case 48: /* annotation: LINE LOCTYPE worldview  */
#line 741 "gram.y"
                                 { sysline.loctype = (yyvsp[0].pset); }
#line 5355 "y.tab.c"
    break;

  case 49: /* annotation: LINE LINEWIDTH NUMBER  */
#line 742 "gram.y"
                                { sysline.linew = (int) (yyvsp[0].val); }
#line 5361 "y.tab.c"
    break;

  case 50: /* annotation: LINE LINESTYLE NUMBER  */
#line 743 "gram.y"
                                { sysline.lines = (int) (yyvsp[0].val); }
#line 5367 "y.tab.c"
    break;

  case 51: /* annotation: LINE COLOR NUMBER  */
#line 744 "gram.y"
                            { sysline.color = (int) (yyvsp[0].val); }
#line 5373 "y.tab.c"
    break;

  case 52: /* annotation: LINE ARROW NUMBER  */
#line 745 "gram.y"
                            { sysline.arrow = (int) (yyvsp[0].val); }
#line 5379 "y.tab.c"
    break;

  case 53: /* annotation: LINE ARROW SIZE NUMBER  */
#line 746 "gram.y"
                                 { sysline.asize = (yyvsp[0].val); }
#line 5385 "y.tab.c"
    break;

  case 54: /* annotation: LINE ARROW TYPE NUMBER  */
#line 747 "gram.y"
                                 { sysline.atype = (int) (yyvsp[0].val); }
#line 5391 "y.tab.c"
    break;

  case 55: /* annotation: LINE DEF  */
#line 749 "gram.y"
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
#line 5407 "y.tab.c"
    break;

  case 56: /* annotation: CLEAR STRING  */
#line 760 "gram.y"
                       { do_clear_text(); }
#line 5413 "y.tab.c"
    break;

  case 57: /* annotation: WITH STRING  */
#line 761 "gram.y"
                      { curstring = next_string(); }
#line 5419 "y.tab.c"
    break;

  case 58: /* annotation: WITH STRING NUMBER  */
#line 762 "gram.y"
                             { curstring = (int) (yyvsp[0].val); }
#line 5425 "y.tab.c"
    break;

  case 59: /* annotation: STRING onoff  */
#line 763 "gram.y"
                       { pstr[curstring].active = (yyvsp[0].pset); }
#line 5431 "y.tab.c"
    break;

  case 60: /* annotation: STRING GRAPHNO  */
#line 764 "gram.y"
                         { pstr[curstring].gno = (yyvsp[0].pset); }
#line 5437 "y.tab.c"
    break;

  case 61: /* annotation: STRING expr ',' expr  */
#line 766 "gram.y"
        {
	    pstr[curstring].x = (yyvsp[-2].val);
	    pstr[curstring].y = (yyvsp[0].val);
	}
#line 5446 "y.tab.c"
    break;

  case 62: /* annotation: STRING LOCTYPE worldview  */
#line 770 "gram.y"
                                   { sysstr.loctype = (yyvsp[0].pset); }
#line 5452 "y.tab.c"
    break;

  case 63: /* annotation: STRING LINEWIDTH NUMBER  */
#line 771 "gram.y"
                                  { sysstr.linew = (int) (yyvsp[0].val); }
#line 5458 "y.tab.c"
    break;

  case 64: /* annotation: STRING COLOR NUMBER  */
#line 772 "gram.y"
                              { sysstr.color = (int) (yyvsp[0].val); }
#line 5464 "y.tab.c"
    break;

  case 65: /* annotation: STRING ROT NUMBER  */
#line 773 "gram.y"
                            { sysstr.rot = (int) (yyvsp[0].val); }
#line 5470 "y.tab.c"
    break;

  case 66: /* annotation: STRING FONTP NUMBER  */
#line 774 "gram.y"
                              { sysstr.font = (int) (yyvsp[0].val); }
#line 5476 "y.tab.c"
    break;

  case 67: /* annotation: STRING JUST NUMBER  */
#line 775 "gram.y"
                             { sysstr.just = (int) (yyvsp[0].val); }
#line 5482 "y.tab.c"
    break;

  case 68: /* annotation: STRING SYMBOL NUMBER  */
#line 776 "gram.y"
                               { sysstr.sym = (int) (yyvsp[0].val); }
#line 5488 "y.tab.c"
    break;

  case 69: /* annotation: STRING SYMBOL LOCTYPE opchoice  */
#line 777 "gram.y"
                                         { sysstr.symloc = (int) (yyvsp[0].pset); }
#line 5494 "y.tab.c"
    break;

  case 70: /* annotation: STRING SYMBOL SIZE NUMBER  */
#line 778 "gram.y"
                                    { sysstr.symsize = (double) (yyvsp[0].val); }
#line 5500 "y.tab.c"
    break;

  case 71: /* annotation: STRING SYMBOL FILL NUMBER  */
#line 779 "gram.y"
                                    { sysstr.symfill = (int) (yyvsp[0].val); }
#line 5506 "y.tab.c"
    break;

  case 72: /* annotation: STRING SYMBOL COLOR NUMBER  */
#line 780 "gram.y"
                                     { sysstr.symcolor = (int) (yyvsp[0].val); }
#line 5512 "y.tab.c"
    break;

  case 73: /* annotation: STRING CHAR SIZE NUMBER  */
#line 781 "gram.y"
                                  { sysstr.charsize = (double) (yyvsp[0].val); }
#line 5518 "y.tab.c"
    break;

  case 74: /* annotation: STRING DEF CHRSTR  */
#line 783 "gram.y"
        {
	    strcpy(pstr[curstring].s, (char *) (yyvsp[0].str));
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
#line 5539 "y.tab.c"
    break;

  case 75: /* animation: STEP  */
#line 802 "gram.y"
             { setistep(); }
#line 5545 "y.tab.c"
    break;

  case 76: /* animation: REVERSE  */
#line 803 "gram.y"
                  { setreverse(); }
#line 5551 "y.tab.c"
    break;

  case 77: /* animation: REWIND  */
#line 804 "gram.y"
                 { setrewind(); }
#line 5557 "y.tab.c"
    break;

  case 78: /* animation: FORWARD  */
#line 805 "gram.y"
                  { setforward(); }
#line 5563 "y.tab.c"
    break;

  case 79: /* animation: WRAP onoff  */
#line 806 "gram.y"
                     { set_wrap((yyvsp[0].pset) == ON); }
#line 5569 "y.tab.c"
    break;

  case 80: /* animation: GOTO NUMBER  */
#line 807 "gram.y"
                      { goto_step((int) (yyvsp[0].val) - 1); }
#line 5575 "y.tab.c"
    break;

  case 81: /* animation: RUN  */
#line 808 "gram.y"
              { setirun(); }
#line 5581 "y.tab.c"
    break;

  case 82: /* animation: RUN NUMBER ',' NUMBER  */
#line 809 "gram.y"
                                { runsteps((int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5587 "y.tab.c"
    break;

  case 83: /* animation: BATCH RUN NUMBER ',' NUMBER  */
#line 810 "gram.y"
                                      { batchrunsteps((int) (yyvsp[-2].val), (int) (yyvsp[0].val), 1); }
#line 5593 "y.tab.c"
    break;

  case 84: /* animation: BATCH RUN NUMBER ',' NUMBER SKIP NUMBER  */
#line 811 "gram.y"
                                                  { batchrunsteps((int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5599 "y.tab.c"
    break;

  case 85: /* animation: BATCH PREFIX CHRSTR  */
#line 812 "gram.y"
                              { strcpy(batchprefix, (char *) (yyvsp[0].str)); }
#line 5605 "y.tab.c"
    break;

  case 86: /* animation: STOP  */
#line 813 "gram.y"
               { setistop(); }
#line 5611 "y.tab.c"
    break;

  case 87: /* misc: SYSTEM CHRSTR  */
#line 817 "gram.y"
                      { system((yyvsp[0].str)); }
#line 5617 "y.tab.c"
    break;

  case 88: /* misc: RUN BATCH CHRSTR  */
#line 819 "gram.y"
        {
	    gotbatch = 1;
	    batchfile[0] = 0;
	    strcpy(batchfile, (yyvsp[0].str));
	}
#line 5627 "y.tab.c"
    break;

  case 89: /* misc: CHDIR CHRSTR  */
#line 824 "gram.y"
                       {
		if (chdir((char *) (yyvsp[0].str)) < 0) {
			sprintf(buf, "chdir() to %s failed", (char *) (yyvsp[0].str));
			errwin(buf);
		}
	}
#line 5638 "y.tab.c"
    break;

  case 90: /* misc: RESET  */
#line 830 "gram.y"
                { setreset_world(); }
#line 5644 "y.tab.c"
    break;

  case 91: /* misc: REDRAW  */
#line 831 "gram.y"
                 { setredraw_world(); }
#line 5650 "y.tab.c"
    break;

  case 92: /* misc: SLEEP NUMBER  */
#line 832 "gram.y"
                       { sleep((int) (yyvsp[0].val)); }
#line 5656 "y.tab.c"
    break;

  case 93: /* misc: QUIT  */
#line 833 "gram.y"
               { exit(0); }
#line 5662 "y.tab.c"
    break;

  case 94: /* misc: ZOOM expr ',' expr ',' expr ',' expr  */
#line 834 "gram.y"
                                               { my_blowup((yyvsp[-6].val), (yyvsp[-4].val), (yyvsp[-2].val), (yyvsp[0].val)); }
#line 5668 "y.tab.c"
    break;

  case 95: /* misc: EXPAND  */
#line 835 "gram.y"
                 { page(page_per, 4); }
#line 5674 "y.tab.c"
    break;

  case 96: /* misc: SHRINK  */
#line 836 "gram.y"
                 { page(page_per, 5); }
#line 5680 "y.tab.c"
    break;

  case 97: /* misc: PAGE LEFT  */
#line 837 "gram.y"
                    { page(page_per, 0); }
#line 5686 "y.tab.c"
    break;

  case 98: /* misc: PAGE RIGHT  */
#line 838 "gram.y"
                     { page(page_per, 1); }
#line 5692 "y.tab.c"
    break;

  case 99: /* misc: PAGE UP  */
#line 839 "gram.y"
                  { page(page_per, 2); }
#line 5698 "y.tab.c"
    break;

  case 100: /* misc: PAGE DOWN  */
#line 840 "gram.y"
                    { page(page_per, 3); }
#line 5704 "y.tab.c"
    break;

  case 101: /* misc: PAGE expr  */
#line 841 "gram.y"
                    { page_per = (yyvsp[0].val); }
#line 5710 "y.tab.c"
    break;

  case 102: /* misc: PAGE INOUT NUMBER  */
#line 842 "gram.y"
                            { scrollinout_proc((int) (yyvsp[0].val)); }
#line 5716 "y.tab.c"
    break;

  case 103: /* misc: LINK PAGE onoff  */
#line 843 "gram.y"
                          { scrolling_islinked = (yyvsp[0].pset) == ON; }
#line 5722 "y.tab.c"
    break;

  case 104: /* misc: LOG CHRSTR  */
#line 845 "gram.y"
        {
	    if ((logfp = fopen((yyvsp[0].str), "w")) != NULL) {
		logfile = 1;
		printf("Opened logfile %s\n", (yyvsp[0].str));
	    } else {
		logfile = 0;
		printf("Failed to open logfile %s\n", (yyvsp[0].str));
	    }
	}
#line 5736 "y.tab.c"
    break;

  case 105: /* misc: CLOSE LOG  */
#line 855 "gram.y"
        {
	    if (logfp != NULL) {
		printf("Closing logfile\n");
		logfile = 0;
		fclose(logfp);
		logfp = NULL;
	    }
	}
#line 5749 "y.tab.c"
    break;

  case 106: /* misc: INCLUDE IMAGE CHRSTR  */
#line 863 "gram.y"
                               { }
#line 5755 "y.tab.c"
    break;

  case 107: /* misc: PRINT  */
#line 864 "gram.y"
                { batchrunstep(0); }
#line 5761 "y.tab.c"
    break;

  case 108: /* misc: ECHO CHRSTR  */
#line 865 "gram.y"
                      {
	    if (inwin) { set_left_footer((yyvsp[0].str)); }
	    else { printf("%s\n", (yyvsp[0].str)); }
	}
#line 5770 "y.tab.c"
    break;

  case 109: /* misc: CMAP NUMBER ',' NUMBER ',' NUMBER ',' NUMBER  */
#line 870 "gram.y"
        { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5776 "y.tab.c"
    break;

  case 110: /* misc: COLORMAP NUMBER ',' NUMBER ',' NUMBER ',' NUMBER  */
#line 872 "gram.y"
        { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5782 "y.tab.c"
    break;

  case 111: /* misc: COLOR NUMBER ',' NUMBER ',' NUMBER ',' NUMBER  */
#line 874 "gram.y"
        { set_colormapdata((int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 5788 "y.tab.c"
    break;

  case 112: /* $@1: %empty  */
#line 875 "gram.y"
                 { setisol = &(g[curg].salip); }
#line 5794 "y.tab.c"
    break;

  case 114: /* $@2: %empty  */
#line 876 "gram.y"
                    { setisol = &(g[curg].velmagip); }
#line 5800 "y.tab.c"
    break;

  case 116: /* models: WITH TEANL NUMBER  */
#line 880 "gram.y"
                          { 
		setflow = &g[curg].flowf[(int) (yyvsp[0].val)]; 
	}
#line 5808 "y.tab.c"
    break;

  case 117: /* models: TEANL flowprops  */
#line 883 "gram.y"
                          { }
#line 5814 "y.tab.c"
    break;

  case 118: /* models: READ TEANL  */
#line 884 "gram.y"
                     {}
#line 5820 "y.tab.c"
    break;

  case 119: /* models: WITH ADCIRC NUMBER  */
#line 885 "gram.y"
                             { 
		setflow = &g[curg].flowt[(int) (yyvsp[0].val)]; 
                elcirc_flowno = (int) (yyvsp[0].val);
	}
#line 5829 "y.tab.c"
    break;

  case 120: /* models: ADCIRC flowprops  */
#line 889 "gram.y"
                           {}
#line 5835 "y.tab.c"
    break;

  case 121: /* models: READ ADCIRC ELEV NUMBER GRID CHRSTR CHRSTR  */
#line 891 "gram.y"
          {
	      int fno = (int) (yyvsp[-3].val);
              readbin_adcirc_elev(fno, (char *) (yyvsp[-1].str), (char *) (yyvsp[0].str));
              set_clock(0, flowt[fno].start, 
				flowt[fno].stop, flowt[fno].step,
                                  flowt[fno].nsteps);
              load_clock(ADCIRC, fno);
          }
#line 5848 "y.tab.c"
    break;

  case 122: /* models: WITH ELCIRC NUMBER  */
#line 899 "gram.y"
                             { 
		setflow = &g[curg].flowt[(int) (yyvsp[0].val)]; 
                elcirc_flowno = (int) (yyvsp[0].val);
		g[curg].curadc3d = (int) (yyvsp[0].val); 
	}
#line 5858 "y.tab.c"
    break;

  case 123: /* models: ELCIRC flowprops  */
#line 904 "gram.y"
                           {}
#line 5864 "y.tab.c"
    break;

  case 124: /* models: ELCIRC RUN NUMBER  */
#line 905 "gram.y"
                            { curadc3d = (int) (yyvsp[0].val); }
#line 5870 "y.tab.c"
    break;

  case 125: /* models: ELCIRC GRID NUMBER CHRSTR  */
#line 907 "gram.y"
        {
		ReplaceElcircGrid((int) (yyvsp[-1].val), (char *) (yyvsp[0].str));
	}
#line 5878 "y.tab.c"
    break;

  case 126: /* models: READ ELCIRC NUMBER REGION CHRSTR  */
#line 911 "gram.y"
        { 
		ReadElcircRegionFile((int) (yyvsp[-2].val), (char *) (yyvsp[0].str));
	}
#line 5886 "y.tab.c"
    break;

  case 127: /* models: READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER LEVEL NUMBER  */
#line 915 "gram.y"
        {
		ReadElcirc((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (int) (yyvsp[0].val) - 1, (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5899 "y.tab.c"
    break;

  case 128: /* models: READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER LEVEL NUMBER APPEND  */
#line 924 "gram.y"
        {
		ReadElcirc((int) (yyvsp[-10].val), (char *) (yyvsp[-9].str), (int) (yyvsp[-1].val) - 1, (int) (yyvsp[-7].val), (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 0, 0, 1);
                elcirc_flowno = (int) (yyvsp[-10].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5912 "y.tab.c"
    break;

  case 129: /* models: READ ELCIRC NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER DEPTH NUMBER  */
#line 933 "gram.y"
        {
		ReadElcircDepth((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (char *) NULL, (double) (yyvsp[0].val), (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5925 "y.tab.c"
    break;

  case 130: /* models: READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER DEPTH NUMBER  */
#line 942 "gram.y"
        {
/* read at a given depth relative to the free surface */
		ReadElcircDepthFromFreeSurface((int) (yyvsp[-9].val), (char *) (yyvsp[-8].str), (char *) NULL, (double) (yyvsp[0].val), (int) (yyvsp[-6].val), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[-9].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5939 "y.tab.c"
    break;

  case 131: /* models: READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER  */
#line 952 "gram.y"
        {
		ReadElcircSurf((int) (yyvsp[-7].val), (char *) (yyvsp[-6].str), 0, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[-7].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5952 "y.tab.c"
    break;

  case 132: /* models: READ ELCIRC SURFACE NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER APPEND  */
#line 961 "gram.y"
        {
		ReadElcircSurf((int) (yyvsp[-8].val), (char *) (yyvsp[-7].str), 0, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), (int) (yyvsp[-1].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[-8].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5965 "y.tab.c"
    break;

  case 133: /* models: READ ELCIRC BOTTOM NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER  */
#line 970 "gram.y"
        {
		ReadElcircSurf((int) (yyvsp[-7].val), (char *) (yyvsp[-6].str), 1, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[-7].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5978 "y.tab.c"
    break;

  case 134: /* models: READ ELCIRC BOTTOM NUMBER CHRSTR START NUMBER STOP NUMBER SKIP NUMBER APPEND  */
#line 979 "gram.y"
        {
		ReadElcircSurf((int) (yyvsp[-8].val), (char *) (yyvsp[-7].str), 1, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), (int) (yyvsp[-1].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[-8].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
#line 5991 "y.tab.c"
    break;

  case 135: /* models: WITH ELCIRC TRANSECT NUMBER  */
#line 987 "gram.y"
                                      { 
			curtrans = (int) (yyvsp[0].val); 
			settrans = &(trans[curtrans]); 
	}
#line 6000 "y.tab.c"
    break;

  case 136: /* models: SET ELCIRC TRANSECT NUMBER  */
#line 991 "gram.y"
                                     { curtrans = (int) (yyvsp[0].val); }
#line 6006 "y.tab.c"
    break;

  case 137: /* $@3: %empty  */
#line 992 "gram.y"
                          { setisol = &(g[curg].trans[curtrans].ip); }
#line 6012 "y.tab.c"
    break;

  case 139: /* models: ELCIRC TRANSECT FLOW CHRSTR  */
#line 993 "gram.y"
                                      {strcpy(settrans->uvname, (char *) (yyvsp[0].str));}
#line 6018 "y.tab.c"
    break;

  case 140: /* models: ELCIRC TRANSECT VERTICAL FLOW CHRSTR  */
#line 994 "gram.y"
                                               {strcpy(settrans->vvname, (char *) (yyvsp[0].str));}
#line 6024 "y.tab.c"
    break;

  case 141: /* models: ELCIRC TRANSECT FLOW NUMBER  */
#line 995 "gram.y"
                                      { settrans->flowno = (int) (yyvsp[0].val);}
#line 6030 "y.tab.c"
    break;

  case 142: /* models: ELCIRC TRANSECT SALINITY CHRSTR  */
#line 996 "gram.y"
                                          {strcpy(settrans->salname, (char *) (yyvsp[0].str));}
#line 6036 "y.tab.c"
    break;

  case 143: /* models: ELCIRC TRANSECT ELEV CHRSTR  */
#line 997 "gram.y"
                                      {strcpy(settrans->elevname, (char *) (yyvsp[0].str));}
#line 6042 "y.tab.c"
    break;

  case 144: /* models: ELCIRC TRANSECT GRAPH NUMBER  */
#line 998 "gram.y"
                                       { settrans->gno = (int) (yyvsp[0].val);}
#line 6048 "y.tab.c"
    break;

  case 145: /* models: ELCIRC TRANSECT DISPLAY GRAPH NUMBER  */
#line 999 "gram.y"
                                               { settrans->transgno = (int) (yyvsp[0].val);}
#line 6054 "y.tab.c"
    break;

  case 146: /* models: ELCIRC TRANSECT DISPLAY LINE onoff  */
#line 1000 "gram.y"
                                             { settrans->display = (int) (yyvsp[0].pset);}
#line 6060 "y.tab.c"
    break;

  case 147: /* models: ELCIRC TRANSECT DISPLAY onoff  */
#line 1001 "gram.y"
                                        { g[curg].trans[curtrans].display = (int) (yyvsp[0].pset);}
#line 6066 "y.tab.c"
    break;

  case 148: /* models: ELCIRC TRANSECT DISPLAY MAG onoff  */
#line 1002 "gram.y"
                                            { g[curg].trans[curtrans].display_mag = (int) (yyvsp[0].pset);}
#line 6072 "y.tab.c"
    break;

  case 149: /* models: ELCIRC TRANSECT MAXLEVELS NUMBER  */
#line 1003 "gram.y"
                                           { }
#line 6078 "y.tab.c"
    break;

  case 150: /* models: ELCIRC TRANSECT START NUMBER STOP NUMBER SKIP NUMBER  */
#line 1005 "gram.y"
        {
		settrans->start = (int) (yyvsp[-4].val);
		settrans->stop = (int) (yyvsp[-2].val);
		settrans->skip = (int) (yyvsp[0].val);
	}
#line 6088 "y.tab.c"
    break;

  case 151: /* models: ELCIRC TRANSECT SAMPLE NUMBER  */
#line 1010 "gram.y"
                                        { settrans->npts = (int) (yyvsp[0].val); }
#line 6094 "y.tab.c"
    break;

  case 152: /* models: ELCIRC TRANSECT TYPE NUMBER  */
#line 1011 "gram.y"
                                      { settrans->transtype = (int) (yyvsp[0].val); }
#line 6100 "y.tab.c"
    break;

  case 153: /* models: ELCIRC TRANSECT NUMBER ',' expr ',' expr  */
#line 1012 "gram.y"
                                                   { AddTransNXY(settrans, (int) (yyvsp[-4].val), (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 6106 "y.tab.c"
    break;

  case 154: /* models: ELCIRC TRANSECT NODE NUMBER ',' NUMBER  */
#line 1013 "gram.y"
                                                 { AddTransNode(settrans, (int) (yyvsp[-2].val), (int) (yyvsp[0].val)); }
#line 6112 "y.tab.c"
    break;

  case 155: /* models: ELCIRC TRANSECT REGION FILEP CHRSTR  */
#line 1015 "gram.y"
        {
		settrans->transtype = 0;
		strcpy(settrans->transname, (char *) (yyvsp[0].str));
	}
#line 6121 "y.tab.c"
    break;

  case 156: /* models: ELCIRC TRANSECT LINE NUMBER ',' expr ',' expr ',' expr ',' expr  */
#line 1020 "gram.y"
        {
		settrans->transtype = 1;
		settrans->npts = (int) (yyvsp[-8].val);
		settrans->x1 = (double) (yyvsp[-6].val);
		settrans->y1 = (double) (yyvsp[-4].val);
		settrans->x2 = (double) (yyvsp[-2].val);
		settrans->y2 = (double) (yyvsp[0].val);
	}
#line 6134 "y.tab.c"
    break;

  case 157: /* models: READ ELCIRC TRANSECT  */
#line 1029 "gram.y"
        { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6146 "y.tab.c"
    break;

  case 158: /* models: READ ELCIRC TRANSECT APPEND  */
#line 1037 "gram.y"
        { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6158 "y.tab.c"
    break;

  case 159: /* models: READ ELCIRC FLOW TRANSECT  */
#line 1045 "gram.y"
        { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6170 "y.tab.c"
    break;

  case 160: /* models: READ ELCIRC FLOW TRANSECT APPEND  */
#line 1053 "gram.y"
        { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6182 "y.tab.c"
    break;

  case 161: /* models: READ ELCIRC SALINITY TRANSECT  */
#line 1061 "gram.y"
        { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6194 "y.tab.c"
    break;

  case 162: /* models: READ ELCIRC SALINITY TRANSECT APPEND  */
#line 1069 "gram.y"
        { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
#line 6206 "y.tab.c"
    break;

  case 163: /* models: WITH ELCIRC MARKER NUMBER  */
#line 1077 "gram.y"
        {
		elcircmarker = (int) (yyvsp[0].val);
		setadc3d = &adc3d[(int) (yyvsp[0].val)];
		setflow3d = &g[curg].flow3d[(int) (yyvsp[0].val)];
	}
#line 6216 "y.tab.c"
    break;

  case 164: /* models: ELCIRC MARKER FLOW CHRSTR  */
#line 1083 "gram.y"
        {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 6224 "y.tab.c"
    break;

  case 165: /* models: ELCIRC MARKER MAGNITUDE CHRSTR  */
#line 1087 "gram.y"
        {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 6232 "y.tab.c"
    break;

  case 166: /* models: ELCIRC MARKER VECTOR CHRSTR  */
#line 1091 "gram.y"
        {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 6240 "y.tab.c"
    break;

  case 167: /* models: ELCIRC MARKER SCALAR CHRSTR  */
#line 1095 "gram.y"
        {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 6248 "y.tab.c"
    break;

  case 168: /* models: ELCIRC MARKER SALINITY CHRSTR  */
#line 1099 "gram.y"
        {
		strcpy(setadc3d->datafile, (char *) (yyvsp[0].str));
	}
#line 6256 "y.tab.c"
    break;

  case 169: /* models: ELCIRC MARKER ELEV CHRSTR  */
#line 1103 "gram.y"
        {
		strcpy(setadc3d->elevfile, (char *) (yyvsp[0].str));
	}
#line 6264 "y.tab.c"
    break;

  case 170: /* models: READ ELCIRC MARKER NODE NUMBER START NUMBER STEP NUMBER SKIP NUMBER  */
#line 1107 "gram.y"
        {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[-6].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0);
	}
#line 6275 "y.tab.c"
    break;

  case 171: /* models: READ ELCIRC MARKER NODE NUMBER START NUMBER STEP NUMBER SKIP NUMBER APPEND  */
#line 1114 "gram.y"
        {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[-7].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 1);
	}
#line 6286 "y.tab.c"
    break;

  case 172: /* models: READ ELCIRC MARKER XY NUMBER ',' NUMBER START NUMBER STEP NUMBER SKIP NUMBER  */
#line 1121 "gram.y"
        {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[-8].val);
		setadc3d->y = (double) (yyvsp[-6].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), 0);
	}
#line 6298 "y.tab.c"
    break;

  case 173: /* models: READ ELCIRC MARKER XY NUMBER ',' NUMBER START NUMBER STEP NUMBER SKIP NUMBER APPEND  */
#line 1129 "gram.y"
        {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[-9].val);
		setadc3d->y = (double) (yyvsp[-7].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[-5].val), (int) (yyvsp[-3].val), 1);
	}
#line 6310 "y.tab.c"
    break;

  case 174: /* $@4: %empty  */
#line 1136 "gram.y"
                               { setisol = &(g[curg].salip); }
#line 6316 "y.tab.c"
    break;

  case 176: /* $@5: %empty  */
#line 1137 "gram.y"
                                 { setisol = &(g[curg].salip); }
#line 6322 "y.tab.c"
    break;

  case 178: /* $@6: %empty  */
#line 1138 "gram.y"
                                  { setisol = &(g[curg].velmagip); }
#line 6328 "y.tab.c"
    break;

  case 180: /* $@7: %empty  */
#line 1139 "gram.y"
                        { setprops = &(setflow3d->p); }
#line 6334 "y.tab.c"
    break;

  case 182: /* models: ELCIRC MARKER PREC NUMBER ',' NUMBER  */
#line 1141 "gram.y"
        {
		setflow3d->precx = (int) (yyvsp[-2].val);
		setflow3d->precy = (int) (yyvsp[0].val);
	}
#line 6343 "y.tab.c"
    break;

  case 183: /* models: ELCIRC MARKER ATTACH NUMBER  */
#line 1145 "gram.y"
                                      { setflow3d->attach = (int) (yyvsp[0].val); }
#line 6349 "y.tab.c"
    break;

  case 184: /* models: ELCIRC MARKER LOCTYPE worldview  */
#line 1146 "gram.y"
                                          { setflow3d->loctype = (int) (yyvsp[0].pset); }
#line 6355 "y.tab.c"
    break;

  case 185: /* models: ELCIRC MARKER DISPLAY MARKER onoff  */
#line 1147 "gram.y"
                                             { setflow3d->display_marker = (int) (yyvsp[0].pset); }
#line 6361 "y.tab.c"
    break;

  case 186: /* models: ELCIRC MARKER DISPLAY onoff  */
#line 1148 "gram.y"
                                      { setflow3d->display = (int) (yyvsp[0].pset); }
#line 6367 "y.tab.c"
    break;

  case 187: /* models: ELCIRC MARKER COLOR NUMBER  */
#line 1149 "gram.y"
                                     { setflow3d->p.color = (yyvsp[0].val); }
#line 6373 "y.tab.c"
    break;

  case 188: /* models: ELCIRC MARKER LINEWIDTH NUMBER  */
#line 1150 "gram.y"
                                         { setflow3d->p.linew = (yyvsp[0].val); }
#line 6379 "y.tab.c"
    break;

  case 189: /* models: ELCIRC MARKER FILL COLOR NUMBER  */
#line 1151 "gram.y"
                                          { setflow3d->p.fillcol = (yyvsp[0].val); }
#line 6385 "y.tab.c"
    break;

  case 190: /* models: ELCIRC MARKER WORLD expr ',' expr ',' expr ',' expr  */
#line 1153 "gram.y"
        { 
		setflow3d->wx1 = (double) (yyvsp[-6].val); 
		setflow3d->wy1 = (double) (yyvsp[-4].val); 
		setflow3d->wx2 = (double) (yyvsp[-2].val); 
		setflow3d->wy2 = (double) (yyvsp[0].val); 
	}
#line 6396 "y.tab.c"
    break;

  case 191: /* models: ELCIRC MARKER VIEW expr ',' expr  */
#line 1160 "gram.y"
        { 
		setflow3d->vx = (double) (yyvsp[-2].val); 
		setflow3d->vy = (double) (yyvsp[0].val); 
	}
#line 6405 "y.tab.c"
    break;

  case 192: /* models: ELCIRC MARKER LOC expr ',' expr  */
#line 1165 "gram.y"
        { 
		setflow3d->locx = (double) (yyvsp[-2].val); 
		setflow3d->locy = (double) (yyvsp[0].val); 
	}
#line 6414 "y.tab.c"
    break;

  case 193: /* models: ELCIRC MARKER XY expr ',' expr  */
#line 1170 "gram.y"
        { 
		setflow3d->x = (double) (yyvsp[-2].val); 
		setflow3d->y = (double) (yyvsp[0].val); 
	}
#line 6423 "y.tab.c"
    break;

  case 194: /* flowprops: DISPLAY flowonoff  */
#line 1177 "gram.y"
                          { setflow->display = (yyvsp[0].pset);  }
#line 6429 "y.tab.c"
    break;

  case 195: /* flowprops: DISPLAY ELEV onoff  */
#line 1178 "gram.y"
                             { setflow->display_elev = (yyvsp[0].pset); }
#line 6435 "y.tab.c"
    break;

  case 196: /* flowprops: DISPLAY ELEV DEPTH onoff  */
#line 1179 "gram.y"
                                   { setflow->display_elevdepth = (yyvsp[0].pset); }
#line 6441 "y.tab.c"
    break;

  case 197: /* flowprops: DISPLAY ELEV VALUE onoff  */
#line 1180 "gram.y"
                                   { setflow->display_maxelevval = (yyvsp[0].pset); }
#line 6447 "y.tab.c"
    break;

  case 198: /* flowprops: DISPLAY ELEV MAXP onoff  */
#line 1181 "gram.y"
                                  { setflow->display_maxelev = (yyvsp[0].pset); }
#line 6453 "y.tab.c"
    break;

  case 199: /* flowprops: DISPLAY ELEV AMP onoff  */
#line 1182 "gram.y"
                                 { setflow->display_amp = (yyvsp[0].pset); }
#line 6459 "y.tab.c"
    break;

  case 200: /* flowprops: DISPLAY ELEV PHASE onoff  */
#line 1183 "gram.y"
                                   { setflow->display_phase = (yyvsp[0].pset); }
#line 6465 "y.tab.c"
    break;

  case 201: /* flowprops: DISPLAY ELEV MARKERS onoff  */
#line 1184 "gram.y"
                                     { setflow->display_elevmarkers = (yyvsp[0].pset); }
#line 6471 "y.tab.c"
    break;

  case 202: /* flowprops: DISPLAY FLOW MAG onoff  */
#line 1185 "gram.y"
                                 { setflow->display_mag = (yyvsp[0].pset); }
#line 6477 "y.tab.c"
    break;

  case 203: /* flowprops: DISPLAY FLOW WIND onoff  */
#line 1186 "gram.y"
                                  { setflow->display_wind = (yyvsp[0].pset); }
#line 6483 "y.tab.c"
    break;

  case 204: /* flowprops: DISPLAY INUNDATION onoff  */
#line 1187 "gram.y"
                                   { setflow->display_inun = (yyvsp[-1].pset); }
#line 6489 "y.tab.c"
    break;

  case 205: /* flowprops: COLOR NUMBER  */
#line 1188 "gram.y"
                       { setflow->p.color = (yyvsp[0].val); }
#line 6495 "y.tab.c"
    break;

  case 206: /* $@8: %empty  */
#line 1189 "gram.y"
               { setisol = &(setflow->elevip); }
#line 6501 "y.tab.c"
    break;

  case 208: /* $@9: %empty  */
#line 1190 "gram.y"
                    { setisol = &(setflow->maxelevip); }
#line 6507 "y.tab.c"
    break;

  case 210: /* $@10: %empty  */
#line 1191 "gram.y"
              { setisol = &(setflow->ampip); }
#line 6513 "y.tab.c"
    break;

  case 212: /* $@11: %empty  */
#line 1192 "gram.y"
                { setisol = &(setflow->phaseip); }
#line 6519 "y.tab.c"
    break;

  case 214: /* $@12: %empty  */
#line 1193 "gram.y"
                   { setisol = &(setflow->magip); }
#line 6525 "y.tab.c"
    break;

  case 216: /* flowprops: FLOW FREQ NUMBER  */
#line 1194 "gram.y"
                           { setflow->flowfreq = (yyvsp[0].val); }
#line 6531 "y.tab.c"
    break;

  case 217: /* flowprops: FREQ NUMBER  */
#line 1195 "gram.y"
                      { setflow->freq = (yyvsp[0].val); }
#line 6537 "y.tab.c"
    break;

  case 218: /* flowprops: ELEVMARKER NUMBER  */
#line 1196 "gram.y"
                            { setelevmarker = &(setflow->em[(int) (yyvsp[0].val)]); }
#line 6543 "y.tab.c"
    break;

  case 219: /* flowprops: SAMPLE FLOW onoff  */
#line 1197 "gram.y"
                            { setflow->sample = (int) (yyvsp[0].pset); }
#line 6549 "y.tab.c"
    break;

  case 220: /* flowprops: SAMPLE FLOW READ CHRSTR  */
#line 1198 "gram.y"
                                  { ReadSampleFlow(setflow, (char *) (yyvsp[0].str)); }
#line 6555 "y.tab.c"
    break;

  case 221: /* flowprops: SAMPLE FLOW MINP NUMBER  */
#line 1199 "gram.y"
                                  { SetMinSampleFlow(elcirc_flowno, (double) (yyvsp[0].val)); }
#line 6561 "y.tab.c"
    break;

  case 222: /* flowprops: SAMPLE FLOW READ XY CHRSTR  */
#line 1200 "gram.y"
                                     { ReadSampleFlowXY(setflow, (char *) (yyvsp[0].str)); }
#line 6567 "y.tab.c"
    break;

  case 223: /* flowprops: SAMPLE FLOW TYPE XY  */
#line 1201 "gram.y"
                              { setflow->samptype = XY; }
#line 6573 "y.tab.c"
    break;

  case 224: /* flowprops: SAMPLE FLOW TYPE NODE  */
#line 1202 "gram.y"
                                { setflow->samptype = NODE; }
#line 6579 "y.tab.c"
    break;

  case 225: /* flowprops: SAMPLE FLOW NODE NUMBER  */
#line 1203 "gram.y"
                                  { AddSampleFlowNode(setflow, (int) (yyvsp[0].val) - 1); }
#line 6585 "y.tab.c"
    break;

  case 226: /* flowprops: SAMPLE FLOW ELEMENT NUMBER  */
#line 1204 "gram.y"
                                     { AddSampleFlowElem(setflow, (int) (yyvsp[0].val) - 1); }
#line 6591 "y.tab.c"
    break;

  case 227: /* flowprops: SAMPLE FLOW XY expr ',' expr  */
#line 1205 "gram.y"
                                       { AddSampleFlowXY(setflow, (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 6597 "y.tab.c"
    break;

  case 228: /* flowprops: DELETE SAMPLE FLOW NODE NUMBER  */
#line 1206 "gram.y"
                                         { DeleteSampleFlowNode(setflow, (int) (yyvsp[0].val) - 1); }
#line 6603 "y.tab.c"
    break;

  case 229: /* flowprops: DELETE SAMPLE FLOW ELEMENT NUMBER  */
#line 1207 "gram.y"
                                            { DeleteSampleFlowElem(setflow, (int) (yyvsp[0].val) - 1); }
#line 6609 "y.tab.c"
    break;

  case 230: /* flowprops: DELETE SAMPLE FLOW XY expr ',' expr  */
#line 1208 "gram.y"
                                              { DeleteSampleFlowXY(setflow, (double) (yyvsp[-2].val), (double) (yyvsp[0].val)); }
#line 6615 "y.tab.c"
    break;

  case 231: /* elevmarker: WITH ELEVMARKER NUMBER  */
#line 1212 "gram.y"
                               { setelevmarker = &(setflow->em[(int) (yyvsp[0].val)]); }
#line 6621 "y.tab.c"
    break;

  case 232: /* elevmarker: ELEVMARKER ACTIVE onoff  */
#line 1213 "gram.y"
                                  { setelevmarker->active = (int) (yyvsp[0].pset); }
#line 6627 "y.tab.c"
    break;

  case 233: /* elevmarker: ELEVMARKER TYPE onoff  */
#line 1214 "gram.y"
                                { setelevmarker->type = (int) (yyvsp[0].pset); }
#line 6633 "y.tab.c"
    break;

  case 234: /* elevmarker: ELEVMARKER DISPLAY onoff  */
#line 1215 "gram.y"
                                   { setelevmarker->display = (int) (yyvsp[0].pset); }
#line 6639 "y.tab.c"
    break;

  case 235: /* $@13: %empty  */
#line 1216 "gram.y"
                     { setprops = &(setelevmarker->p); }
#line 6645 "y.tab.c"
    break;

  case 237: /* elevmarker: ELEVMARKER LOCTYPE worldview  */
#line 1217 "gram.y"
                                       { setelevmarker->loctype = (int) (yyvsp[0].pset); }
#line 6651 "y.tab.c"
    break;

  case 238: /* elevmarker: ELEVMARKER NODE NUMBER  */
#line 1218 "gram.y"
                                 { setelevmarker->node = (int) (yyvsp[0].val); }
#line 6657 "y.tab.c"
    break;

  case 239: /* elevmarker: ELEVMARKER LOC expr ',' expr  */
#line 1220 "gram.y"
        { 
		setelevmarker->locx = (double) (yyvsp[-2].val); 
		setelevmarker->locy = (double) (yyvsp[0].val); 
	}
#line 6666 "y.tab.c"
    break;

  case 240: /* elevmarker: ELEVMARKER MINMAX expr ',' expr  */
#line 1225 "gram.y"
        { 
		setelevmarker->emin = (double) (yyvsp[-2].val); 
		setelevmarker->emax = (double) (yyvsp[0].val); 
	}
#line 6675 "y.tab.c"
    break;

  case 241: /* grid: GRID NUMBER  */
#line 1236 "gram.y"
                    { setgrid = &g[curg].grid[(int) (yyvsp[0].val)]; }
#line 6681 "y.tab.c"
    break;

  case 242: /* grid: WITH GRID NUMBER  */
#line 1237 "gram.y"
                           { setgrid = &g[curg].grid[(int) (yyvsp[0].val)]; }
#line 6687 "y.tab.c"
    break;

  case 243: /* grid: GRID DISPLAY onoff  */
#line 1238 "gram.y"
                             { if (checkptr(setgrid, f_string)) setgrid->display = (yyvsp[0].pset); }
#line 6693 "y.tab.c"
    break;

  case 244: /* grid: GRID BATH DISPLAY onoff  */
#line 1239 "gram.y"
                                  { if (checkptr(setgrid, f_string)) setgrid->display_bath = (yyvsp[0].pset); }
#line 6699 "y.tab.c"
    break;

  case 245: /* grid: GRID COURANT DISPLAY onoff DT NUMBER  */
#line 1240 "gram.y"
                                               { if (checkptr(setgrid, f_string)) setgrid->display_courant = (yyvsp[-2].pset); }
#line 6705 "y.tab.c"
    break;

  case 246: /* grid: GRID COURANT DISPLAY VALUE onoff DT NUMBER  */
#line 1241 "gram.y"
                                                     { if (checkptr(setgrid, f_string)) setgrid->display_courantn = (yyvsp[-2].pset); }
#line 6711 "y.tab.c"
    break;

  case 247: /* grid: GRID BOUNDARY DISPLAY onoff  */
#line 1242 "gram.y"
                                      { if (checkptr(setgrid, f_string)) setgrid->display_boundary = (yyvsp[0].pset); }
#line 6717 "y.tab.c"
    break;

  case 248: /* $@14: %empty  */
#line 1243 "gram.y"
               { if (checkptr(setgrid, f_string)) setprops = &(setgrid->p); }
#line 6723 "y.tab.c"
    break;

  case 250: /* $@15: %empty  */
#line 1244 "gram.y"
                        { if (checkptr(setgrid, f_string)) setprops = &(setgrid->bp); }
#line 6729 "y.tab.c"
    break;

  case 252: /* grid: GRID NODES DISPLAY onoff  */
#line 1245 "gram.y"
                                   { if (checkptr(setgrid, f_string)) setgrid->display_nodes = (yyvsp[0].pset); }
#line 6735 "y.tab.c"
    break;

  case 253: /* grid: GRID ELEMENTS DISPLAY onoff  */
#line 1246 "gram.y"
                                      { if (checkptr(setgrid, f_string)) setgrid->display_elements = (yyvsp[0].pset); }
#line 6741 "y.tab.c"
    break;

  case 254: /* grid: GRID DEPTH DISPLAY onoff  */
#line 1247 "gram.y"
                                   { if (checkptr(setgrid, f_string)) setgrid->display_depths = (yyvsp[0].pset); }
#line 6747 "y.tab.c"
    break;

  case 255: /* grid: GRID FILL onoff  */
#line 1248 "gram.y"
                          { if (checkptr(setgrid, f_string)) setgrid->display_gridf = (yyvsp[0].pset); }
#line 6753 "y.tab.c"
    break;

  case 256: /* $@16: %empty  */
#line 1249 "gram.y"
                    { if (checkptr(setgrid, f_string)) setisol = &(setgrid->ip); }
#line 6759 "y.tab.c"
    break;

  case 258: /* grid: GRAPHNO AUTOSCALE  */
#line 1251 "gram.y"
        { 
		autoscale_grid((int) (yyvsp[-1].pset), g[(int) (yyvsp[-1].pset)].curgrid); 
		set_defaults((int) (yyvsp[-1].pset));
	}
#line 6768 "y.tab.c"
    break;

  case 259: /* grid: AUTOSCALE  */
#line 1256 "gram.y"
        { 
		autoscale_grid(curg, g[curg].curgrid); 
		set_defaults(curg);
	}
#line 6777 "y.tab.c"
    break;

  case 260: /* grid: READ GRID NUMBER CHRSTR  */
#line 1261 "gram.y"
        {
		extern int readgridfile;
		readgrid((int) (yyvsp[-1].val), (char *) (yyvsp[0].str));
		readgridfile = 1;
	}
#line 6787 "y.tab.c"
    break;

  case 261: /* drogues: WITH DROGUES NUMBER  */
#line 1269 "gram.y"
                            { setdrogs = &g[curg].drogues[(int) (yyvsp[0].val)]; }
#line 6793 "y.tab.c"
    break;

  case 262: /* drogues: DROGUES DISPLAY onoff  */
#line 1270 "gram.y"
                                {  setdrogs->display = (int) (yyvsp[0].pset); }
#line 6799 "y.tab.c"
    break;

  case 263: /* $@17: %empty  */
#line 1271 "gram.y"
                  { setprops = &(setdrogs->p); }
#line 6805 "y.tab.c"
    break;

  case 265: /* drogues: READ DROGUES CHRSTR  */
#line 1273 "gram.y"
        {
    	if (!readdrogues(curdrog, (char *) (yyvsp[0].str), -1, 0, 0)) {
        	fprintf(stderr, "Error reading file %s", (char *) (yyvsp[0].str));
    	} else {
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
    	}
	}
#line 6820 "y.tab.c"
    break;

  case 266: /* drogues: READ DROGUES CHRSTR START NUMBER STOP NUMBER SKIP NUMBER  */
#line 1284 "gram.y"
        {
		readdrogues(0, (char *) (yyvsp[-6].str), (int) (yyvsp[-4].val), (int) (yyvsp[-2].val), (int) (yyvsp[0].val));
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
	}
#line 6832 "y.tab.c"
    break;

  case 267: /* isolines: ISOLINES LEGEND onoff  */
#line 1294 "gram.y"
                              { setisol->lactive = (yyvsp[0].pset); }
#line 6838 "y.tab.c"
    break;

  case 268: /* isolines: ISOLINES LEGEND LAYOUT horv  */
#line 1295 "gram.y"
                                      { setisol->layout = (yyvsp[0].pset); }
#line 6844 "y.tab.c"
    break;

  case 269: /* isolines: ISOLINES LEGEND LABEL onoff  */
#line 1296 "gram.y"
                                      { setisol->llabels = (yyvsp[-1].pset); }
#line 6850 "y.tab.c"
    break;

  case 270: /* isolines: ISOLINES LEGEND FRAMEP onoff  */
#line 1297 "gram.y"
                                       { setisol->frame = (int) (yyvsp[0].pset); }
#line 6856 "y.tab.c"
    break;

  case 271: /* isolines: ISOLINES LEGEND FRAMEP COLOR NUMBER  */
#line 1298 "gram.y"
                                              { setisol->framecol = (int) (yyvsp[0].val); }
#line 6862 "y.tab.c"
    break;

  case 272: /* isolines: ISOLINES NUMBER  */
#line 1299 "gram.y"
                          { setisol->nisol = (int) (yyvsp[0].val); }
#line 6868 "y.tab.c"
    break;

  case 273: /* isolines: ISOLINES TYPE NUMBER  */
#line 1300 "gram.y"
                               { setisol->type = (int) (yyvsp[0].val); }
#line 6874 "y.tab.c"
    break;

  case 274: /* isolines: ISOLINES SET TYPE NUMBER  */
#line 1301 "gram.y"
                                   { setisol->isoltype = (int) (yyvsp[0].val); }
#line 6880 "y.tab.c"
    break;

  case 275: /* isolines: ISOLINES FILL TYPE NUMBER  */
#line 1302 "gram.y"
                                    { setisol->visflag = (int) (yyvsp[0].val); }
#line 6886 "y.tab.c"
    break;

  case 276: /* isolines: ISOLINES SAVE ON  */
#line 1303 "gram.y"
                           { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
#line 6899 "y.tab.c"
    break;

  case 277: /* isolines: ISOLINES SAVE OFF  */
#line 1311 "gram.y"
                            { 
		int i;
		setisol->writeflag = 0; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
#line 6912 "y.tab.c"
    break;

  case 278: /* isolines: ISOLINES SAVE FILEP CHRSTR  */
#line 1319 "gram.y"
                                     { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		strncpy(setisol->wname, (char *) (yyvsp[0].str), 1023);
	}
#line 6925 "y.tab.c"
    break;

  case 279: /* isolines: ISOLINES SAVE ISOLINE NUMBER FILEP CHRSTR  */
#line 1327 "gram.y"
                                                    { 
		int i;
		setisol->writeflag = 1; 
		strncpy(setisol->wname, (char *) (yyvsp[0].str), 1023);
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
    		setisol->writelevel[(int) (yyvsp[-2].val)] = 1;
	}
#line 6939 "y.tab.c"
    break;

  case 280: /* isolines: ISOLINES SAVE ISOLINE NUMBER  */
#line 1336 "gram.y"
                                       { 
    		setisol->writelevel[(int) (yyvsp[0].val)] = 1;
	}
#line 6947 "y.tab.c"
    break;

  case 281: /* isolines: ISOLINES START expr STEP expr  */
#line 1339 "gram.y"
                                        {
            setisol->cis[0] = (double) (yyvsp[-2].val);
            setisol->cint = (double) (yyvsp[0].val);
	    if (setisol->isoltype == 0) {
		int i;
		for (i=1;i< 16;i++) {
		    setisol->cis[i] = setisol->cis[0] + i * setisol->cint;
		}
	    }
        }
#line 6962 "y.tab.c"
    break;

  case 282: /* isolines: ISOLINES MINMAX expr ',' expr  */
#line 1349 "gram.y"
                                        { setisol->cmin = (yyvsp[-2].val); setisol->cmax = (yyvsp[0].val); }
#line 6968 "y.tab.c"
    break;

  case 283: /* isolines: ISOLINE expr ',' expr  */
#line 1350 "gram.y"
                                { setisol->cis[(int) (yyvsp[-2].val)] = (yyvsp[0].val); }
#line 6974 "y.tab.c"
    break;

  case 284: /* isolines: ISOLINES LEGEND LOCTYPE worldview  */
#line 1351 "gram.y"
                                            { setisol->loctype = (yyvsp[0].pset); }
#line 6980 "y.tab.c"
    break;

  case 285: /* isolines: ISOLINES LEGEND expr ',' expr  */
#line 1352 "gram.y"
                                        {
            setisol->x = (double) (yyvsp[-2].val);
            setisol->y = (double) (yyvsp[0].val);
        }
#line 6989 "y.tab.c"
    break;

  case 286: /* isolines: ISOLINE NUMBER COLOR NUMBER  */
#line 1356 "gram.y"
                                      { setisol->color[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 6995 "y.tab.c"
    break;

  case 287: /* isolines: ISOLINE NUMBER LINEWIDTH NUMBER  */
#line 1357 "gram.y"
                                          { setisol->linew[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 7001 "y.tab.c"
    break;

  case 288: /* isolines: ISOLINE NUMBER LINESTYLE NUMBER  */
#line 1358 "gram.y"
                                          { setisol->lines[(int) (yyvsp[-2].val)] = (int) (yyvsp[0].val); }
#line 7007 "y.tab.c"
    break;

  case 289: /* $@18: %empty  */
#line 1359 "gram.y"
                   { setprops = &(setisol->p); }
#line 7013 "y.tab.c"
    break;

  case 291: /* $@19: %empty  */
#line 1360 "gram.y"
                          { setprops = &(setisol->p); }
#line 7019 "y.tab.c"
    break;

  case 293: /* isolines: ISOLINES LEGEND SIZE NUMBER ',' NUMBER  */
#line 1361 "gram.y"
                                                 {
            setisol->xlen = (yyvsp[-2].val);
            setisol->ylen = (yyvsp[0].val);
        }
#line 7028 "y.tab.c"
    break;

  case 294: /* isolines: ISOLINES LEGEND HGAP NUMBER ',' NUMBER  */
#line 1365 "gram.y"
                                                 {
            setisol->xgap = (yyvsp[-2].val);
            setisol->ygap = (yyvsp[0].val);
        }
#line 7037 "y.tab.c"
    break;

  case 295: /* histboxes: WITH HISTBOX NUMBER  */
#line 1372 "gram.y"
                            { sethistbox = &g[curg].hbox[(int) (yyvsp[0].val)]; }
#line 7043 "y.tab.c"
    break;

  case 296: /* $@20: %empty  */
#line 1373 "gram.y"
                  { setprops = &(sethistbox->p); }
#line 7049 "y.tab.c"
    break;

  case 298: /* histboxes: HISTBOX PREC NUMBER ',' NUMBER  */
#line 1375 "gram.y"
        {
		sethistbox->precx = (int) (yyvsp[-2].val);
		sethistbox->precy = (int) (yyvsp[0].val);
	}
#line 7058 "y.tab.c"
    break;

  case 299: /* histboxes: HISTBOX ATTACH NUMBER  */
#line 1379 "gram.y"
                                { sethistbox->attach = (int) (yyvsp[0].val); }
#line 7064 "y.tab.c"
    break;

  case 300: /* histboxes: HISTBOX TICKP NUMBER ',' NUMBER  */
#line 1380 "gram.y"
                                          { sethistbox->xtickm = (double ) (yyvsp[-2].val); sethistbox->ytickm = (double ) (yyvsp[0].val);}
#line 7070 "y.tab.c"
    break;

  case 301: /* histboxes: HISTBOX LOCTYPE worldview  */
#line 1381 "gram.y"
                                    { sethistbox->loctype = (int) (yyvsp[0].pset); }
#line 7076 "y.tab.c"
    break;

  case 302: /* histboxes: HISTBOX DISPLAY MARKER onoff  */
#line 1382 "gram.y"
                                       { sethistbox->display_marker = (int) (yyvsp[0].pset); }
#line 7082 "y.tab.c"
    break;

  case 303: /* histboxes: HISTBOX DISPLAY onoff  */
#line 1383 "gram.y"
                                { sethistbox->display = (int) (yyvsp[0].pset); }
#line 7088 "y.tab.c"
    break;

  case 304: /* histboxes: HISTBOX DISPLAY ADCIRC NUMBER torf  */
#line 1384 "gram.y"
                                             { sethistbox->adcirc[(int) (yyvsp[-1].val)] = (int) (yyvsp[0].pset) == TRUEP; }
#line 7094 "y.tab.c"
    break;

  case 305: /* histboxes: HISTBOX DISPLAY ADCIRC NUMBER COLOR NUMBER  */
#line 1385 "gram.y"
                                                     { sethistbox->ap[(int) (yyvsp[-2].val)].color = (int) (yyvsp[0].val); }
#line 7100 "y.tab.c"
    break;

  case 306: /* histboxes: HISTBOX COLOR NUMBER  */
#line 1386 "gram.y"
                               { sethistbox->p.color = (yyvsp[0].val); }
#line 7106 "y.tab.c"
    break;

  case 307: /* histboxes: HISTBOX LINEWIDTH NUMBER  */
#line 1387 "gram.y"
                                   { sethistbox->p.linew = (yyvsp[0].val); }
#line 7112 "y.tab.c"
    break;

  case 308: /* histboxes: HISTBOX FILL COLOR NUMBER  */
#line 1388 "gram.y"
                                    { sethistbox->p.fillcol = (yyvsp[0].val); }
#line 7118 "y.tab.c"
    break;

  case 309: /* histboxes: HISTBOX READ NUMBER CHRSTR  */
#line 1389 "gram.y"
                                     { read_hist((int) (yyvsp[-1].val), TIME, (char *) (yyvsp[0].str)); }
#line 7124 "y.tab.c"
    break;

  case 310: /* histboxes: HISTBOX DISPLAY HISTORY torf  */
#line 1390 "gram.y"
                                       { sethistbox->thist = (int) (yyvsp[0].pset) == TRUEP; }
#line 7130 "y.tab.c"
    break;

  case 311: /* histboxes: HISTBOX DISPLAY HISTORY COLOR NUMBER  */
#line 1391 "gram.y"
                                               { sethistbox->hp.color = (int) (yyvsp[0].val); }
#line 7136 "y.tab.c"
    break;

  case 312: /* histboxes: HISTBOX WORLD expr ',' expr ',' expr ',' expr  */
#line 1393 "gram.y"
        { 
		sethistbox->wx1 = (double) (yyvsp[-6].val); 
		sethistbox->wy1 = (double) (yyvsp[-4].val); 
		sethistbox->wx2 = (double) (yyvsp[-2].val); 
		sethistbox->wy2 = (double) (yyvsp[0].val); 
	}
#line 7147 "y.tab.c"
    break;

  case 313: /* histboxes: HISTBOX VIEW expr ',' expr  */
#line 1400 "gram.y"
        { 
		sethistbox->vx = (double) (yyvsp[-2].val); 
		sethistbox->vy = (double) (yyvsp[0].val); 
	}
#line 7156 "y.tab.c"
    break;

  case 314: /* histboxes: HISTBOX LOC expr ',' expr  */
#line 1405 "gram.y"
        { 
		sethistbox->locx = (double) (yyvsp[-2].val); 
		sethistbox->locy = (double) (yyvsp[0].val); 
	}
#line 7165 "y.tab.c"
    break;

  case 315: /* histboxes: HISTBOX XY expr ',' expr  */
#line 1410 "gram.y"
        { 
		sethistbox->x = (double) (yyvsp[-2].val); 
		sethistbox->y = (double) (yyvsp[0].val); 
	}
#line 7174 "y.tab.c"
    break;

  case 316: /* zoomboxes: WITH ZOOMBOX NUMBER  */
#line 1417 "gram.y"
                            { setzoombox = &g[curg].zbox[(int) (yyvsp[0].val)]; }
#line 7180 "y.tab.c"
    break;

  case 317: /* $@21: %empty  */
#line 1418 "gram.y"
                  { setprops = &(setzoombox->p); }
#line 7186 "y.tab.c"
    break;

  case 319: /* zoomboxes: ZOOMBOX PREC NUMBER ',' NUMBER  */
#line 1420 "gram.y"
        {
		setzoombox->precx = (int) (yyvsp[-2].val);
		setzoombox->precy = (int) (yyvsp[0].val);
	}
#line 7195 "y.tab.c"
    break;

  case 320: /* zoomboxes: ZOOMBOX ATTACH NUMBER  */
#line 1424 "gram.y"
                                { setzoombox->attach = (int) (yyvsp[0].val); }
#line 7201 "y.tab.c"
    break;

  case 321: /* zoomboxes: ZOOMBOX LOCTYPE worldview  */
#line 1425 "gram.y"
                                    { setzoombox->loctype = (int) (yyvsp[0].pset); }
#line 7207 "y.tab.c"
    break;

  case 322: /* zoomboxes: ZOOMBOX DISPLAY MARKER onoff  */
#line 1426 "gram.y"
                                       { setzoombox->display_marker = (int) (yyvsp[0].pset); }
#line 7213 "y.tab.c"
    break;

  case 323: /* zoomboxes: ZOOMBOX onoff  */
#line 1427 "gram.y"
                        { setzoombox->active = (int) (yyvsp[0].pset); }
#line 7219 "y.tab.c"
    break;

  case 324: /* zoomboxes: ZOOMBOX DISPLAY onoff  */
#line 1428 "gram.y"
                                { setzoombox->display = (int) (yyvsp[0].pset); }
#line 7225 "y.tab.c"
    break;

  case 325: /* zoomboxes: ZOOMBOX ZOOM NUMBER  */
#line 1429 "gram.y"
                              { setzoombox->expand = (int) (yyvsp[0].val); }
#line 7231 "y.tab.c"
    break;

  case 326: /* zoomboxes: ZOOMBOX SCALE NUMBER  */
#line 1430 "gram.y"
                               { setzoombox->expand = (int) (yyvsp[0].val); }
#line 7237 "y.tab.c"
    break;

  case 327: /* zoomboxes: ZOOMBOX COLOR NUMBER  */
#line 1431 "gram.y"
                               { setzoombox->p.color = (yyvsp[0].val); }
#line 7243 "y.tab.c"
    break;

  case 328: /* zoomboxes: ZOOMBOX LINEWIDTH NUMBER  */
#line 1432 "gram.y"
                                   { setzoombox->p.linew = (yyvsp[0].val); }
#line 7249 "y.tab.c"
    break;

  case 329: /* zoomboxes: ZOOMBOX FILL COLOR NUMBER  */
#line 1433 "gram.y"
                                    { setzoombox->p.fillcol = (yyvsp[0].val); }
#line 7255 "y.tab.c"
    break;

  case 330: /* zoomboxes: ZOOMBOX WORLD expr ',' expr ',' expr ',' expr  */
#line 1435 "gram.y"
        { 
		setzoombox->wx1 = (double) (yyvsp[-6].val); 
		setzoombox->wy1 = (double) (yyvsp[-4].val); 
		setzoombox->wx2 = (double) (yyvsp[-2].val); 
		setzoombox->wy2 = (double) (yyvsp[0].val); 
	}
#line 7266 "y.tab.c"
    break;

  case 331: /* zoomboxes: ZOOMBOX VIEW expr ',' expr  */
#line 1442 "gram.y"
        { 
		setzoombox->vx = (double) (yyvsp[-2].val); 
		setzoombox->vy = (double) (yyvsp[0].val); 
	}
#line 7275 "y.tab.c"
    break;

  case 332: /* zoomboxes: ZOOMBOX LOC expr ',' expr  */
#line 1447 "gram.y"
        { 
		setzoombox->locx = (double) (yyvsp[-2].val); 
		setzoombox->locy = (double) (yyvsp[0].val); 
	}
#line 7284 "y.tab.c"
    break;

  case 333: /* zoomboxes: ZOOMBOX XY expr ',' expr  */
#line 1452 "gram.y"
        { 
		setzoombox->x = (double) (yyvsp[-2].val); 
		setzoombox->y = (double) (yyvsp[0].val); 
	}
#line 7293 "y.tab.c"
    break;

  case 334: /* wscale: WSCALE onoff  */
#line 1459 "gram.y"
                     { g[curg].wl.active = (yyvsp[0].pset); }
#line 7299 "y.tab.c"
    break;

  case 335: /* wscale: WSCALE LENGTH NUMBER  */
#line 1460 "gram.y"
                               { g[curg].wl.len = (yyvsp[0].val); }
#line 7305 "y.tab.c"
    break;

  case 336: /* wscale: WSCALE SCALE NUMBER  */
#line 1461 "gram.y"
                              { g[curg].wl.scale = (yyvsp[0].val); }
#line 7311 "y.tab.c"
    break;

  case 337: /* wscale: WSCALE COLOR NUMBER  */
#line 1462 "gram.y"
                              { g[curg].wl.p.color = (yyvsp[0].val); }
#line 7317 "y.tab.c"
    break;

  case 338: /* wscale: WSCALE LOCTYPE worldview  */
#line 1463 "gram.y"
                                   { g[curg].wl.loctype = (yyvsp[0].pset); }
#line 7323 "y.tab.c"
    break;

  case 339: /* wscale: WSCALE LOC expr ',' expr  */
#line 1465 "gram.y"
        { 
		g[curg].wl.x = (yyvsp[-2].val);
		g[curg].wl.y = (yyvsp[0].val);
	}
#line 7332 "y.tab.c"
    break;

  case 340: /* wscale: WSCALE UNITS units  */
#line 1470 "gram.y"
        { 
		g[curg].wl.units = (yyvsp[0].pset);
		switch (g[curg].wl.units) {
		case MM: g[curg].wl.unitfac = 0.001; break;
		case CM: g[curg].wl.unitfac = 0.01; break;
		case M: g[curg].wl.unitfac = 1.0; break;
		case KM: g[curg].wl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
#line 7347 "y.tab.c"
    break;

  case 341: /* wscale: WSCALE props  */
#line 1480 "gram.y"
                       { }
#line 7353 "y.tab.c"
    break;

  case 342: /* vscale: VSCALE onoff  */
#line 1485 "gram.y"
                     { g[curg].vl.active = (yyvsp[0].pset); }
#line 7359 "y.tab.c"
    break;

  case 343: /* vscale: VSCALE LENGTH NUMBER  */
#line 1486 "gram.y"
                               { g[curg].vl.len = (yyvsp[0].val); }
#line 7365 "y.tab.c"
    break;

  case 344: /* vscale: VSCALE SCALE NUMBER  */
#line 1487 "gram.y"
                              { g[curg].vl.scale = (yyvsp[0].val); }
#line 7371 "y.tab.c"
    break;

  case 345: /* vscale: VSCALE COLOR NUMBER  */
#line 1488 "gram.y"
                              { g[curg].vl.p.color = (yyvsp[0].val); }
#line 7377 "y.tab.c"
    break;

  case 346: /* vscale: VSCALE LOCTYPE worldview  */
#line 1489 "gram.y"
                                   { g[curg].vl.loctype = (yyvsp[0].pset); }
#line 7383 "y.tab.c"
    break;

  case 347: /* vscale: VSCALE LOC expr ',' expr  */
#line 1491 "gram.y"
        { 
		g[curg].vl.x = (yyvsp[-2].val);
		g[curg].vl.y = (yyvsp[0].val);
	}
#line 7392 "y.tab.c"
    break;

  case 348: /* vscale: VSCALE UNITS units  */
#line 1496 "gram.y"
        { 
		g[curg].vl.units = (yyvsp[0].pset);
		switch (g[curg].vl.units) {
		case MM: g[curg].vl.unitfac = 0.001; break;
		case CM: g[curg].vl.unitfac = 0.01; break;
		case M: g[curg].vl.unitfac = 1.0; break;
		case KM: g[curg].vl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
#line 7407 "y.tab.c"
    break;

  case 349: /* vscale: VSCALE props  */
#line 1506 "gram.y"
                       { }
#line 7413 "y.tab.c"
    break;

  case 350: /* mapscale: MAPSCALE onoff  */
#line 1510 "gram.y"
                       { g[curg].mapscale.active = (yyvsp[0].pset); }
#line 7419 "y.tab.c"
    break;

  case 351: /* mapscale: MAPSCALE LENGTH NUMBER  */
#line 1511 "gram.y"
                                 { g[curg].mapscale.len = (yyvsp[0].val); }
#line 7425 "y.tab.c"
    break;

  case 352: /* mapscale: MAPSCALE COLOR NUMBER  */
#line 1512 "gram.y"
                                { g[curg].mapscale.p.color = (yyvsp[0].val); }
#line 7431 "y.tab.c"
    break;

  case 353: /* mapscale: MAPSCALE SCALE NUMBER  */
#line 1513 "gram.y"
                                { g[curg].mapscale.scale = (yyvsp[0].val); }
#line 7437 "y.tab.c"
    break;

  case 354: /* mapscale: MAPSCALE LOCTYPE worldview  */
#line 1514 "gram.y"
                                     { g[curg].mapscale.loctype = (yyvsp[0].pset); }
#line 7443 "y.tab.c"
    break;

  case 355: /* mapscale: MAPSCALE LOC expr ',' expr  */
#line 1516 "gram.y"
        { 
		g[curg].mapscale.x = (yyvsp[-2].val);
		g[curg].mapscale.y = (yyvsp[0].val);
	}
#line 7452 "y.tab.c"
    break;

  case 356: /* mapscale: MAPSCALE UNITS units  */
#line 1521 "gram.y"
        { 
		g[curg].mapscale.units = (yyvsp[0].pset);
		switch (g[curg].mapscale.units) {
		case MM: g[curg].mapscale.unitfac = 0.001; break;
		case CM: g[curg].mapscale.unitfac = 0.01; break;
		case M: g[curg].mapscale.unitfac = 1.0; break;
		case KM: g[curg].mapscale.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for mapscape scale\n"); break;
		}
	}
#line 7467 "y.tab.c"
    break;

  case 357: /* mapscale: MAPSCALE props  */
#line 1531 "gram.y"
                         { }
#line 7473 "y.tab.c"
    break;

  case 358: /* tidalclock: TIDALCLOCK onoff  */
#line 1535 "gram.y"
                         { g[curg].tidalclock.active = (yyvsp[0].pset); }
#line 7479 "y.tab.c"
    break;

  case 359: /* tidalclock: TIDALCLOCK COLOR NUMBER  */
#line 1536 "gram.y"
                                  { g[curg].tidalclock.p.color = (yyvsp[0].val); }
#line 7485 "y.tab.c"
    break;

  case 360: /* tidalclock: TIDALCLOCK FILL COLOR NUMBER  */
#line 1537 "gram.y"
                                       { g[curg].tidalclock.p.fillcol = (yyvsp[0].val); }
#line 7491 "y.tab.c"
    break;

  case 361: /* tidalclock: TIDALCLOCK TOTAL TIME NUMBER  */
#line 1538 "gram.y"
                                       { g[curg].tidalclock.total_time = (yyvsp[0].val); }
#line 7497 "y.tab.c"
    break;

  case 362: /* tidalclock: TIDALCLOCK LOCTYPE worldview  */
#line 1539 "gram.y"
                                       { g[curg].tidalclock.loctype = (yyvsp[0].pset); }
#line 7503 "y.tab.c"
    break;

  case 363: /* tidalclock: TIDALCLOCK LOC expr ',' expr  */
#line 1541 "gram.y"
        { 
		g[curg].tidalclock.x = (yyvsp[-2].val);
		g[curg].tidalclock.y = (yyvsp[0].val);
	}
#line 7512 "y.tab.c"
    break;

  case 364: /* tidalclock: TIDALCLOCK props  */
#line 1545 "gram.y"
                           { }
#line 7518 "y.tab.c"
    break;

  case 365: /* timeinfo: TIMEINFO onoff  */
#line 1549 "gram.y"
                       { g[curg].timeinfo.active = (yyvsp[0].pset); }
#line 7524 "y.tab.c"
    break;

  case 366: /* timeinfo: TIMEINFO START CHRSTR  */
#line 1551 "gram.y"
        { 
		strcpy(g[curg].timeinfo.start, (char *) (yyvsp[0].str)); 
		time_info_start(curg);
	}
#line 7533 "y.tab.c"
    break;

  case 367: /* timeinfo: TIMEINFO expr ',' expr  */
#line 1556 "gram.y"
        {
	    g[curg].timeinfo.x = (yyvsp[-2].val);
	    g[curg].timeinfo.y = (yyvsp[0].val);
	}
#line 7542 "y.tab.c"
    break;

  case 368: /* timeinfo: TIMEINFO LOCTYPE worldview  */
#line 1560 "gram.y"
                                     { g[curg].timeinfo.loctype = (yyvsp[0].pset); }
#line 7548 "y.tab.c"
    break;

  case 369: /* timeinfo: TIMEINFO LINEWIDTH NUMBER  */
#line 1561 "gram.y"
                                    { g[curg].timeinfo.linew = (int) (yyvsp[0].val); }
#line 7554 "y.tab.c"
    break;

  case 370: /* timeinfo: TIMEINFO COLOR NUMBER  */
#line 1562 "gram.y"
                                { g[curg].timeinfo.color = (int) (yyvsp[0].val); }
#line 7560 "y.tab.c"
    break;

  case 371: /* timeinfo: TIMEINFO ROT NUMBER  */
#line 1563 "gram.y"
                              { g[curg].timeinfo.rot = (int) (yyvsp[0].val); }
#line 7566 "y.tab.c"
    break;

  case 372: /* timeinfo: TIMEINFO FONTP NUMBER  */
#line 1564 "gram.y"
                                { g[curg].timeinfo.font = (int) (yyvsp[0].val); }
#line 7572 "y.tab.c"
    break;

  case 373: /* timeinfo: TIMEINFO JUST NUMBER  */
#line 1565 "gram.y"
                               { g[curg].timeinfo.just = (int) (yyvsp[0].val); }
#line 7578 "y.tab.c"
    break;

  case 374: /* timeinfo: TIMEINFO CHAR SIZE NUMBER  */
#line 1566 "gram.y"
                                    { g[curg].timeinfo.charsize = (double) (yyvsp[0].val); }
#line 7584 "y.tab.c"
    break;

  case 375: /* timeinfo: TIMEINFO FORMAT CHRSTR  */
#line 1567 "gram.y"
                                 { strcpy(g[curg].timeinfo.format, (char *) (yyvsp[0].str)); }
#line 7590 "y.tab.c"
    break;

  case 376: /* timeline: TIMELINE onoff  */
#line 1571 "gram.y"
                       { g[curg].timeline.active = (int) (yyvsp[0].pset); }
#line 7596 "y.tab.c"
    break;

  case 377: /* timeline: TIMELINE LENGTH NUMBER  */
#line 1572 "gram.y"
                                 { g[curg].timeline.len = (int) (yyvsp[0].val); }
#line 7602 "y.tab.c"
    break;

  case 378: /* timeline: TIMELINE WIDTH NUMBER  */
#line 1573 "gram.y"
                                { g[curg].timeline.width = (int) (yyvsp[0].val); }
#line 7608 "y.tab.c"
    break;

  case 379: /* timeline: TIMELINE START NUMBER  */
#line 1574 "gram.y"
                                { g[curg].timeline.start = (double) (yyvsp[0].val); }
#line 7614 "y.tab.c"
    break;

  case 380: /* timeline: TIMELINE STOP NUMBER  */
#line 1575 "gram.y"
                               { g[curg].timeline.stop = (double) (yyvsp[0].val); }
#line 7620 "y.tab.c"
    break;

  case 381: /* timeline: TIMELINE STEP NUMBER  */
#line 1576 "gram.y"
                               { g[curg].timeline.step = (double) (yyvsp[0].val); }
#line 7626 "y.tab.c"
    break;

  case 382: /* timeline: TIMELINE PREC NUMBER  */
#line 1577 "gram.y"
                               { g[curg].timeline.p.prec = (int) (yyvsp[0].val); }
#line 7632 "y.tab.c"
    break;

  case 383: /* timeline: TIMELINE UNITS NUMBER  */
#line 1578 "gram.y"
                                { g[curg].timeline.units = (int) (yyvsp[0].val); }
#line 7638 "y.tab.c"
    break;

  case 384: /* timeline: TIMELINE COLOR NUMBER  */
#line 1579 "gram.y"
                                { g[curg].timeline.c1 = g[curg].timeline.c3 = (int) (yyvsp[0].val); }
#line 7644 "y.tab.c"
    break;

  case 385: /* timeline: TIMELINE FILL COLOR NUMBER  */
#line 1580 "gram.y"
                                     { g[curg].timeline.c2 = (int) (yyvsp[0].val); }
#line 7650 "y.tab.c"
    break;

  case 386: /* timeline: TIMELINE LOCTYPE worldview  */
#line 1581 "gram.y"
                                     { g[curg].timeline.loctype = (int) (yyvsp[0].pset); }
#line 7656 "y.tab.c"
    break;

  case 387: /* timeline: TIMELINE LOC expr ',' expr  */
#line 1583 "gram.y"
        { 
		g[curg].timeline.x = (double) (yyvsp[-2].val);
		g[curg].timeline.y = (double) (yyvsp[0].val);
	}
#line 7665 "y.tab.c"
    break;

  case 388: /* timeline: TIMELINE props  */
#line 1587 "gram.y"
                         { }
#line 7671 "y.tab.c"
    break;

  case 389: /* slice: WITH SLICE NUMBER  */
#line 1591 "gram.y"
                          { setslice = &g[curg].sbox[(int) (yyvsp[0].val)]; }
#line 7677 "y.tab.c"
    break;

  case 390: /* $@22: %empty  */
#line 1592 "gram.y"
                { setprops = &(setslice->p); }
#line 7683 "y.tab.c"
    break;

  case 392: /* slice: SLICE PREC NUMBER ',' NUMBER  */
#line 1594 "gram.y"
        {
		setslice->precx = (int) (yyvsp[-2].val);
		setslice->precy = (int) (yyvsp[0].val);
	}
#line 7692 "y.tab.c"
    break;

  case 393: /* slice: SLICE ATTACH NUMBER  */
#line 1598 "gram.y"
                              { setslice->attach = (int) (yyvsp[0].val); }
#line 7698 "y.tab.c"
    break;

  case 394: /* slice: SLICE LOCTYPE worldview  */
#line 1599 "gram.y"
                                  { setslice->loctype = (int) (yyvsp[0].pset); }
#line 7704 "y.tab.c"
    break;

  case 395: /* slice: SLICE DISPLAY MARKER onoff  */
#line 1600 "gram.y"
                                     { setslice->display_marker = (int) (yyvsp[0].pset); }
#line 7710 "y.tab.c"
    break;

  case 396: /* slice: SLICE onoff  */
#line 1601 "gram.y"
                      { setslice->active = (int) (yyvsp[0].pset); }
#line 7716 "y.tab.c"
    break;

  case 397: /* slice: SLICE DISPLAY onoff  */
#line 1602 "gram.y"
                              { setslice->display = (int) (yyvsp[0].pset); }
#line 7722 "y.tab.c"
    break;

  case 398: /* slice: SLICE COLOR NUMBER  */
#line 1603 "gram.y"
                             { setslice->p.color = (yyvsp[0].val); }
#line 7728 "y.tab.c"
    break;

  case 399: /* slice: SLICE LINEWIDTH NUMBER  */
#line 1604 "gram.y"
                                 { setslice->p.linew = (yyvsp[0].val); }
#line 7734 "y.tab.c"
    break;

  case 400: /* slice: SLICE FILL COLOR NUMBER  */
#line 1605 "gram.y"
                                  { setslice->p.fillcol = (yyvsp[0].val); }
#line 7740 "y.tab.c"
    break;

  case 401: /* slice: SLICE WORLD expr ',' expr ',' expr ',' expr  */
#line 1607 "gram.y"
        { 
		setslice->wx1 = (double) (yyvsp[-6].val); 
		setslice->wy1 = (double) (yyvsp[-4].val); 
		setslice->wx2 = (double) (yyvsp[-2].val); 
		setslice->wy2 = (double) (yyvsp[0].val); 
	}
#line 7751 "y.tab.c"
    break;

  case 402: /* slice: SLICE VIEW expr ',' expr  */
#line 1614 "gram.y"
        { 
		setslice->vx = (double) (yyvsp[-2].val); 
		setslice->vy = (double) (yyvsp[0].val); 
	}
#line 7760 "y.tab.c"
    break;

  case 403: /* slice: SLICE LOC expr ',' expr  */
#line 1619 "gram.y"
        { 
		setslice->locx = (double) (yyvsp[-2].val); 
		setslice->locy = (double) (yyvsp[0].val); 
	}
#line 7769 "y.tab.c"
    break;

  case 404: /* slice: SLICE XY expr ',' expr  */
#line 1624 "gram.y"
        { 
		setslice->x = (double) (yyvsp[-2].val); 
		setslice->y = (double) (yyvsp[0].val); 
	}
#line 7778 "y.tab.c"
    break;

  case 405: /* props: PROP COLOR NUMBER  */
#line 1631 "gram.y"
                          { setprops->color = (yyvsp[0].val); }
#line 7784 "y.tab.c"
    break;

  case 406: /* props: PROP LINEWIDTH NUMBER  */
#line 1632 "gram.y"
                                { setprops->linew = (yyvsp[0].val); }
#line 7790 "y.tab.c"
    break;

  case 407: /* props: PROP LINESTYLE NUMBER  */
#line 1633 "gram.y"
                                { setprops->lines = (yyvsp[0].val); }
#line 7796 "y.tab.c"
    break;

  case 408: /* props: PROP FORMAT formatchoice  */
#line 1634 "gram.y"
                                   { setprops->format = (yyvsp[0].pset); }
#line 7802 "y.tab.c"
    break;

  case 409: /* props: PROP FONTP NUMBER  */
#line 1635 "gram.y"
                            { setprops->font = (yyvsp[0].val); }
#line 7808 "y.tab.c"
    break;

  case 410: /* props: PROP PREC NUMBER  */
#line 1636 "gram.y"
                           { setprops->prec = (yyvsp[0].val); }
#line 7814 "y.tab.c"
    break;

  case 411: /* props: PROP CHAR SIZE NUMBER  */
#line 1637 "gram.y"
                                { setprops->charsize = (yyvsp[0].val); }
#line 7820 "y.tab.c"
    break;

  case 412: /* props: PROP SYMBOL NUMBER  */
#line 1638 "gram.y"
                             { setprops->symbol = (yyvsp[0].val); }
#line 7826 "y.tab.c"
    break;

  case 413: /* props: PROP SYMBOL SIZE NUMBER  */
#line 1639 "gram.y"
                                  { setprops->symsize = (yyvsp[0].val); }
#line 7832 "y.tab.c"
    break;

  case 414: /* props: PROP FILL onoff  */
#line 1640 "gram.y"
                          { setprops->fill = (yyvsp[0].pset); }
#line 7838 "y.tab.c"
    break;

  case 415: /* props: PROP FILL filltype  */
#line 1641 "gram.y"
                             { setprops->fillusing = (yyvsp[0].pset); }
#line 7844 "y.tab.c"
    break;

  case 416: /* props: PROP FILL COLOR NUMBER  */
#line 1642 "gram.y"
                                 { setprops->fillcol = (yyvsp[0].val); }
#line 7850 "y.tab.c"
    break;

  case 417: /* props: PROP FILL PATTERN NUMBER  */
#line 1643 "gram.y"
                                   { setprops->fillpat = (yyvsp[0].val); }
#line 7856 "y.tab.c"
    break;

  case 418: /* props: PROP ARROW NUMBER  */
#line 1644 "gram.y"
                            { setprops->arrow = (yyvsp[0].val); }
#line 7862 "y.tab.c"
    break;

  case 419: /* props: PROP ARROW TYPE NUMBER  */
#line 1645 "gram.y"
                                 { setprops->atype = (yyvsp[0].val); }
#line 7868 "y.tab.c"
    break;

  case 420: /* props: PROP ARROW SIZE NUMBER  */
#line 1646 "gram.y"
                                 { setprops->asize = (yyvsp[0].val); }
#line 7874 "y.tab.c"
    break;

  case 421: /* graph: WITH GRAPHNO  */
#line 1650 "gram.y"
                     { curg = (int) (yyvsp[0].pset); }
#line 7880 "y.tab.c"
    break;

  case 422: /* graph: WITH GRAPH NUMBER  */
#line 1651 "gram.y"
                            { curg = (int) (yyvsp[0].val); }
#line 7886 "y.tab.c"
    break;

  case 423: /* graph: KILL GRAPHNO  */
#line 1652 "gram.y"
                       { kill_graph((yyvsp[0].pset)); }
#line 7892 "y.tab.c"
    break;

  case 424: /* graph: KILL GRAPHS  */
#line 1653 "gram.y"
                      { kill_graph(maxgraph); }
#line 7898 "y.tab.c"
    break;

  case 425: /* graph: LOCATOR onoff  */
#line 1655 "gram.y"
        {
	    extern int go_locateflag;
	    go_locateflag = ((yyvsp[0].pset) == ON);
	}
#line 7907 "y.tab.c"
    break;

  case 426: /* graph: FOCUS GRAPHNO  */
#line 1660 "gram.y"
        {
	    cg = curg = (int) (yyvsp[0].pset);
	    draw_focus(curg);
	    defineworld(g[curg].w.xg1, g[curg].w.yg1, g[curg].w.xg2, g[curg].w.yg2, 
			islogx(curg), islogy(curg));
	    viewport(g[curg].v.xv1, g[curg].v.yv1, g[curg].v.xv2, g[curg].v.yv2);
	    draw_focus(curg);
	    update_all(curg);
	}
#line 7921 "y.tab.c"
    break;

  case 427: /* graph: FOCUS onoff  */
#line 1669 "gram.y"
                      { draw_focus_flag = (yyvsp[0].pset); }
#line 7927 "y.tab.c"
    break;

  case 428: /* graph: FOCUS SET  */
#line 1670 "gram.y"
                    { focus_policy = (yyvsp[0].pset); }
#line 7933 "y.tab.c"
    break;

  case 429: /* graph: FOCUS FOLLOWS  */
#line 1671 "gram.y"
                        { focus_policy = (yyvsp[0].pset); }
#line 7939 "y.tab.c"
    break;

  case 430: /* graph: FOCUS CLICK  */
#line 1672 "gram.y"
                      { focus_policy = (yyvsp[0].pset); }
#line 7945 "y.tab.c"
    break;

  case 431: /* graph: SOURCE sourcetype  */
#line 1673 "gram.y"
                            { cursource = (yyvsp[0].pset); }
#line 7951 "y.tab.c"
    break;

  case 432: /* graph: PUSH  */
#line 1674 "gram.y"
               { push_world(); }
#line 7957 "y.tab.c"
    break;

  case 433: /* graph: POP  */
#line 1675 "gram.y"
              { pop_world(); }
#line 7963 "y.tab.c"
    break;

  case 434: /* graph: CYCLE  */
#line 1676 "gram.y"
                { cycle_world_stack(); }
#line 7969 "y.tab.c"
    break;

  case 435: /* graph: STACK NUMBER  */
#line 1677 "gram.y"
                       {
	    if ((int) (yyvsp[0].val) > 0)
		show_world_stack((int) (yyvsp[0].val) - 1);
	}
#line 7978 "y.tab.c"
    break;

  case 436: /* graph: STACK WORLD expr ',' expr ',' expr ',' expr TICKP expr ',' expr ',' expr ',' expr  */
#line 1682 "gram.y"
        {
	    add_world(curg, (yyvsp[-14].val), (yyvsp[-12].val), (yyvsp[-10].val), (yyvsp[-8].val), (yyvsp[-6].val), (yyvsp[-4].val), (yyvsp[-2].val), (yyvsp[0].val));
	}
#line 7986 "y.tab.c"
    break;

  case 437: /* graph: CLEAR STACK  */
#line 1685 "gram.y"
                      { clear_world_stack(); }
#line 7992 "y.tab.c"
    break;

  case 438: /* graph: WORLD expr ',' expr ',' expr ',' expr  */
#line 1687 "gram.y"
        {
	    g[curg].w.xg1 = (yyvsp[-6].val);
	    g[curg].w.yg1 = (yyvsp[-4].val);
	    g[curg].w.xg2 = (yyvsp[-2].val);
	    g[curg].w.yg2 = (yyvsp[0].val);
	}
#line 8003 "y.tab.c"
    break;

  case 439: /* graph: WORLD XMIN expr  */
#line 1693 "gram.y"
                          { g[curg].w.xg1 = (yyvsp[0].val); }
#line 8009 "y.tab.c"
    break;

  case 440: /* graph: WORLD XMAX expr  */
#line 1694 "gram.y"
                          { g[curg].w.xg2 = (yyvsp[0].val); }
#line 8015 "y.tab.c"
    break;

  case 441: /* graph: WORLD YMIN expr  */
#line 1695 "gram.y"
                          { g[curg].w.yg1 = (yyvsp[0].val); }
#line 8021 "y.tab.c"
    break;

  case 442: /* graph: WORLD YMAX expr  */
#line 1696 "gram.y"
                          { g[curg].w.yg2 = (yyvsp[0].val); }
#line 8027 "y.tab.c"
    break;

  case 443: /* graph: VIEW expr ',' expr ',' expr ',' expr  */
#line 1698 "gram.y"
        {
	    g[curg].v.xv1 = (yyvsp[-6].val);
	    g[curg].v.yv1 = (yyvsp[-4].val);
	    g[curg].v.xv2 = (yyvsp[-2].val);
	    g[curg].v.yv2 = (yyvsp[0].val);
	}
#line 8038 "y.tab.c"
    break;

  case 444: /* graph: VIEW XMIN NUMBER  */
#line 1704 "gram.y"
                           { g[curg].v.xv1 = (yyvsp[0].val); }
#line 8044 "y.tab.c"
    break;

  case 445: /* graph: VIEW XMAX NUMBER  */
#line 1705 "gram.y"
                           { g[curg].v.xv2 = (yyvsp[0].val); }
#line 8050 "y.tab.c"
    break;

  case 446: /* graph: VIEW YMIN NUMBER  */
#line 1706 "gram.y"
                           { g[curg].v.yv1 = (yyvsp[0].val); }
#line 8056 "y.tab.c"
    break;

  case 447: /* graph: VIEW YMAX NUMBER  */
#line 1707 "gram.y"
                           { g[curg].v.yv2 = (yyvsp[0].val); }
#line 8062 "y.tab.c"
    break;

  case 448: /* graph: TITLE CHRSTR  */
#line 1708 "gram.y"
                       {
	    set_plotstr_string(&g[curg].labs.title, (char *) (yyvsp[0].str));
	}
#line 8070 "y.tab.c"
    break;

  case 449: /* graph: TITLE FONTP NUMBER  */
#line 1711 "gram.y"
                             {
	    g[curg].labs.title.font = checkon(FONTP, g[curg].labs.title.font, (int) (yyvsp[0].val));
	}
#line 8078 "y.tab.c"
    break;

  case 450: /* graph: TITLE SIZE NUMBER  */
#line 1714 "gram.y"
                            { g[curg].labs.title.charsize = (yyvsp[0].val); }
#line 8084 "y.tab.c"
    break;

  case 451: /* graph: TITLE COLOR NUMBER  */
#line 1715 "gram.y"
                             {
	    g[curg].labs.title.color = checkon(COLOR, g[curg].labs.title.color, (int) (yyvsp[0].val));
	}
#line 8092 "y.tab.c"
    break;

  case 452: /* graph: TITLE LINEWIDTH NUMBER  */
#line 1719 "gram.y"
        {
	    g[curg].labs.title.linew = checkon(LINEWIDTH, g[curg].labs.title.linew, (int) (yyvsp[0].val));
	}
#line 8100 "y.tab.c"
    break;

  case 453: /* graph: SUBTITLE CHRSTR  */
#line 1722 "gram.y"
                          {
	    set_plotstr_string(&g[curg].labs.stitle, (char *) (yyvsp[0].str));
	}
#line 8108 "y.tab.c"
    break;

  case 454: /* graph: SUBTITLE FONTP NUMBER  */
#line 1725 "gram.y"
                                {
	    g[curg].labs.stitle.font = checkon(FONTP, g[curg].labs.stitle.font, (int) (yyvsp[0].val));
	}
#line 8116 "y.tab.c"
    break;

  case 455: /* graph: SUBTITLE SIZE NUMBER  */
#line 1728 "gram.y"
                               { g[curg].labs.stitle.charsize = (yyvsp[0].val); }
#line 8122 "y.tab.c"
    break;

  case 456: /* graph: SUBTITLE COLOR NUMBER  */
#line 1730 "gram.y"
        {
	    g[curg].labs.stitle.color = checkon(COLOR, g[curg].labs.stitle.color, (int) (yyvsp[0].val));
	}
#line 8130 "y.tab.c"
    break;

  case 457: /* graph: SUBTITLE LINEWIDTH NUMBER  */
#line 1734 "gram.y"
        {
	    g[curg].labs.stitle.linew = checkon(LINEWIDTH, g[curg].labs.stitle.color, (int) (yyvsp[0].val));
	}
#line 8138 "y.tab.c"
    break;

  case 458: /* graph: FRAMEP onoff  */
#line 1737 "gram.y"
                       { g[curg].f.active = (yyvsp[0].pset); }
#line 8144 "y.tab.c"
    break;

  case 459: /* graph: FRAMEP TYPE NUMBER  */
#line 1738 "gram.y"
                             { g[curg].f.type = (int) (yyvsp[0].val); }
#line 8150 "y.tab.c"
    break;

  case 460: /* graph: FRAMEP LINESTYLE NUMBER  */
#line 1739 "gram.y"
                                  { g[curg].f.lines = checkon(LINESTYLE, g[curg].f.lines, (int) (yyvsp[0].val)); }
#line 8156 "y.tab.c"
    break;

  case 461: /* graph: FRAMEP LINEWIDTH NUMBER  */
#line 1740 "gram.y"
                                  { g[curg].f.linew = checkon(LINEWIDTH, g[curg].f.linew, (int) (yyvsp[0].val)); }
#line 8162 "y.tab.c"
    break;

  case 462: /* graph: FRAMEP COLOR NUMBER  */
#line 1741 "gram.y"
                              { g[curg].f.color = checkon(COLOR, g[curg].f.color, (int) (yyvsp[0].val)); }
#line 8168 "y.tab.c"
    break;

  case 463: /* graph: FRAMEP FILL onoff  */
#line 1742 "gram.y"
                            { g[curg].f.fillbg = (yyvsp[0].pset); }
#line 8174 "y.tab.c"
    break;

  case 464: /* graph: FRAMEP BACKGROUND COLOR NUMBER  */
#line 1743 "gram.y"
                                         { g[curg].f.bgcolor = (int) (yyvsp[0].val); }
#line 8180 "y.tab.c"
    break;

  case 465: /* graph: GRAPHNO onoff  */
#line 1744 "gram.y"
                        { g[(yyvsp[-1].pset)].active = (yyvsp[0].pset); }
#line 8186 "y.tab.c"
    break;

  case 466: /* graph: GRAPHNO LABEL onoff  */
#line 1745 "gram.y"
                              { g[(yyvsp[-2].pset)].label = (yyvsp[0].pset); }
#line 8192 "y.tab.c"
    break;

  case 467: /* graph: GRAPHNO AUTOSCALE TYPE AUTO  */
#line 1746 "gram.y"
                                      { g[(yyvsp[-3].pset)].auto_type = (yyvsp[0].pset); }
#line 8198 "y.tab.c"
    break;

  case 468: /* graph: GRAPHNO AUTOSCALE TYPE SPEC  */
#line 1747 "gram.y"
                                      { g[(yyvsp[-3].pset)].auto_type = (yyvsp[0].pset); }
#line 8204 "y.tab.c"
    break;

  case 469: /* graph: GRAPHNO AUTOSCALE torf  */
#line 1748 "gram.y"
                                 { g[(yyvsp[-2].pset)].parmsread = ((yyvsp[0].pset) == FALSEP); }
#line 8210 "y.tab.c"
    break;

  case 470: /* graph: GRAPHNO HIDDEN torf  */
#line 1749 "gram.y"
                              { g[(yyvsp[-2].pset)].hidden = ((yyvsp[0].pset) == TRUEP); }
#line 8216 "y.tab.c"
    break;

  case 471: /* graph: GRAPHNO TYPE graphtype  */
#line 1750 "gram.y"
                                 { g[(yyvsp[-2].pset)].type = (yyvsp[0].pset); }
#line 8222 "y.tab.c"
    break;

  case 472: /* graph: GRAPHNO FIXEDPOINT onoff  */
#line 1751 "gram.y"
                                   { g[(yyvsp[-2].pset)].pointset = ((yyvsp[0].pset) == ON); }
#line 8228 "y.tab.c"
    break;

  case 473: /* graph: GRAPHNO FIXEDPOINT FORMAT formatchoice formatchoice  */
#line 1753 "gram.y"
        {
	    g[(yyvsp[-4].pset)].fx = (yyvsp[-1].pset);
	    g[(yyvsp[-4].pset)].fy = (yyvsp[0].pset);
	}
#line 8237 "y.tab.c"
    break;

  case 474: /* graph: GRAPHNO FIXEDPOINT PREC NUMBER ',' NUMBER  */
#line 1758 "gram.y"
        {
	    g[(yyvsp[-5].pset)].px = (yyvsp[-2].val);
	    g[(yyvsp[-5].pset)].py = (yyvsp[0].val);
	}
#line 8246 "y.tab.c"
    break;

  case 475: /* graph: GRAPHNO FIXEDPOINT XY expr ',' expr  */
#line 1763 "gram.y"
        {
	    g[(yyvsp[-5].pset)].dsx = (yyvsp[-2].val);
	    g[(yyvsp[-5].pset)].dsy = (yyvsp[0].val);
	}
#line 8255 "y.tab.c"
    break;

  case 476: /* graph: GRAPHNO FIXEDPOINT TYPE NUMBER  */
#line 1767 "gram.y"
                                         { g[(yyvsp[-3].pset)].pt_type = (int) (yyvsp[0].val); }
#line 8261 "y.tab.c"
    break;

  case 479: /* setaxis: GRAPHS axis  */
#line 1773 "gram.y"
                       {}
#line 8267 "y.tab.c"
    break;

  case 480: /* setaxis: GRAPHS axis axisfeature  */
#line 1774 "gram.y"
                                   {}
#line 8273 "y.tab.c"
    break;

  case 481: /* setaxis: GRAPHS allaxes  */
#line 1775 "gram.y"
                          {}
#line 8279 "y.tab.c"
    break;

  case 482: /* axis: XAXIS  */
#line 1779 "gram.y"
              {}
#line 8285 "y.tab.c"
    break;

  case 483: /* axis: YAXIS  */
#line 1780 "gram.y"
                 {}
#line 8291 "y.tab.c"
    break;

  case 484: /* axis: ALTXAXIS  */
#line 1781 "gram.y"
                    {}
#line 8297 "y.tab.c"
    break;

  case 485: /* axis: ALTYAXIS  */
#line 1782 "gram.y"
                    {}
#line 8303 "y.tab.c"
    break;

  case 486: /* axis: ZEROXAXIS  */
#line 1783 "gram.y"
                     {}
#line 8309 "y.tab.c"
    break;

  case 487: /* axis: ZEROYAXIS  */
#line 1784 "gram.y"
                     {}
#line 8315 "y.tab.c"
    break;

  case 488: /* allaxes: AXES axesprops  */
#line 1788 "gram.y"
                       {}
#line 8321 "y.tab.c"
    break;

  case 489: /* allaxes: XAXES axesprops  */
#line 1789 "gram.y"
                           {}
#line 8327 "y.tab.c"
    break;

  case 490: /* allaxes: YAXES axesprops  */
#line 1790 "gram.y"
                           {}
#line 8333 "y.tab.c"
    break;

  case 491: /* axesprops: onoff  */
#line 1794 "gram.y"
              { set_axis_prop(whichgraph, naxis, (yyvsp[0].pset), 0.0); }
#line 8339 "y.tab.c"
    break;

  case 492: /* axesprops: COLOR NUMBER  */
#line 1795 "gram.y"
                       { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 8345 "y.tab.c"
    break;

  case 493: /* axesprops: LINEWIDTH NUMBER  */
#line 1796 "gram.y"
                           { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 8351 "y.tab.c"
    break;

  case 494: /* axesprops: LINESTYLE NUMBER  */
#line 1797 "gram.y"
                           { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 8357 "y.tab.c"
    break;

  case 495: /* axesprops: FONTP NUMBER  */
#line 1798 "gram.y"
                       { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].val)); }
#line 8363 "y.tab.c"
    break;

  case 496: /* axesprops: CHAR SIZE NUMBER  */
#line 1799 "gram.y"
                           { set_axis_prop(whichgraph, naxis, (yyvsp[-2].pset), (yyvsp[0].val)); }
#line 8369 "y.tab.c"
    break;

  case 497: /* axesprops: GRID onoff  */
#line 1800 "gram.y"
                     { set_axis_prop(whichgraph, naxis, (yyvsp[-1].pset), (yyvsp[0].pset)); }
#line 8375 "y.tab.c"
    break;

  case 498: /* axisfeature: TICKP tickattr  */
#line 1805 "gram.y"
                       {}
#line 8381 "y.tab.c"
    break;

  case 499: /* axisfeature: TICKLABEL ticklabeldesc  */
#line 1806 "gram.y"
                                   {}
#line 8387 "y.tab.c"
    break;

  case 500: /* axisfeature: LABEL axislabeldesc  */
#line 1807 "gram.y"
                               {}
#line 8393 "y.tab.c"
    break;

  case 501: /* axisfeature: BAR axisbardesc  */
#line 1808 "gram.y"
                           {}
#line 8399 "y.tab.c"
    break;

  case 502: /* axisfeature: onoff  */
#line 1809 "gram.y"
                 { g[curg].t[naxis].active = (yyvsp[0].pset); }
#line 8405 "y.tab.c"
    break;

  case 503: /* tickattr: onoff  */
#line 1819 "gram.y"
        {
	    g[curg].t[naxis].t_flag = (yyvsp[0].pset);
	    g[curg].t[naxis].t_mflag = (yyvsp[0].pset);
	}
#line 8414 "y.tab.c"
    break;

  case 504: /* tickattr: MAJOR onoff  */
#line 1823 "gram.y"
                      { g[curg].t[naxis].t_flag = (yyvsp[0].pset); }
#line 8420 "y.tab.c"
    break;

  case 505: /* tickattr: MINOR onoff  */
#line 1824 "gram.y"
                      { g[curg].t[naxis].t_mflag = (yyvsp[0].pset); }
#line 8426 "y.tab.c"
    break;

  case 506: /* tickattr: MAJOR expr  */
#line 1825 "gram.y"
                     { g[curg].t[naxis].tmajor = (yyvsp[0].val); }
#line 8432 "y.tab.c"
    break;

  case 507: /* tickattr: MINOR expr  */
#line 1826 "gram.y"
                     { g[curg].t[naxis].tminor = (yyvsp[0].val); }
#line 8438 "y.tab.c"
    break;

  case 508: /* tickattr: OFFSETX expr  */
#line 1827 "gram.y"
                       { g[curg].t[naxis].offsx = (yyvsp[0].val); }
#line 8444 "y.tab.c"
    break;

  case 509: /* tickattr: OFFSETY expr  */
#line 1828 "gram.y"
                       { g[curg].t[naxis].offsy = (yyvsp[0].val); }
#line 8450 "y.tab.c"
    break;

  case 510: /* tickattr: ALT onoff  */
#line 1829 "gram.y"
                    { g[curg].t[naxis].alt = (yyvsp[0].pset); }
#line 8456 "y.tab.c"
    break;

  case 511: /* tickattr: MINP expr  */
#line 1830 "gram.y"
                    { g[curg].t[naxis].tmin = (yyvsp[0].val); }
#line 8462 "y.tab.c"
    break;

  case 512: /* tickattr: MAXP expr  */
#line 1831 "gram.y"
                    { g[curg].t[naxis].tmax = (yyvsp[0].val); }
#line 8468 "y.tab.c"
    break;

  case 513: /* tickattr: DEFAULT NUMBER  */
#line 1832 "gram.y"
                         { g[curg].t[naxis].t_num = (int) (yyvsp[0].val); }
#line 8474 "y.tab.c"
    break;

  case 514: /* tickattr: inoutchoice  */
#line 1833 "gram.y"
                      { g[curg].t[naxis].t_inout = (yyvsp[0].pset); }
#line 8480 "y.tab.c"
    break;

  case 515: /* tickattr: LOG onoff  */
#line 1834 "gram.y"
                    { g[curg].t[naxis].t_log = (yyvsp[0].pset); }
#line 8486 "y.tab.c"
    break;

  case 516: /* tickattr: SIZE NUMBER  */
#line 1835 "gram.y"
                      { g[curg].t[naxis].t_size = (yyvsp[0].val); }
#line 8492 "y.tab.c"
    break;

  case 517: /* tickattr: MAJOR SIZE NUMBER  */
#line 1836 "gram.y"
                            { g[curg].t[naxis].t_size = (yyvsp[0].val); }
#line 8498 "y.tab.c"
    break;

  case 518: /* tickattr: MINOR SIZE NUMBER  */
#line 1837 "gram.y"
                            { g[curg].t[naxis].t_msize = (yyvsp[0].val); }
#line 8504 "y.tab.c"
    break;

  case 519: /* tickattr: COLOR NUMBER  */
#line 1838 "gram.y"
                       { g[curg].t[naxis].t_color = g[curg].t[naxis].t_mcolor = (int) (yyvsp[0].val); }
#line 8510 "y.tab.c"
    break;

  case 520: /* tickattr: LINEWIDTH NUMBER  */
#line 1839 "gram.y"
                           { g[curg].t[naxis].t_linew = g[curg].t[naxis].t_mlinew = (int) (yyvsp[0].val); }
#line 8516 "y.tab.c"
    break;

  case 521: /* tickattr: MAJOR COLOR NUMBER  */
#line 1840 "gram.y"
                             { g[curg].t[naxis].t_color = (int) (yyvsp[0].val); }
#line 8522 "y.tab.c"
    break;

  case 522: /* tickattr: MINOR COLOR NUMBER  */
#line 1841 "gram.y"
                             { g[curg].t[naxis].t_mcolor = (int) (yyvsp[0].val); }
#line 8528 "y.tab.c"
    break;

  case 523: /* tickattr: MAJOR LINEWIDTH NUMBER  */
#line 1842 "gram.y"
                                 { g[curg].t[naxis].t_linew = (int) (yyvsp[0].val); }
#line 8534 "y.tab.c"
    break;

  case 524: /* tickattr: MINOR LINEWIDTH NUMBER  */
#line 1843 "gram.y"
                                 { g[curg].t[naxis].t_mlinew = (int) (yyvsp[0].val); }
#line 8540 "y.tab.c"
    break;

  case 525: /* tickattr: MAJOR LINESTYLE NUMBER  */
#line 1844 "gram.y"
                                 { g[curg].t[naxis].t_lines = (int) (yyvsp[0].val); }
#line 8546 "y.tab.c"
    break;

  case 526: /* tickattr: MINOR LINESTYLE NUMBER  */
#line 1845 "gram.y"
                                 { g[curg].t[naxis].t_mlines = (int) (yyvsp[0].val); }
#line 8552 "y.tab.c"
    break;

  case 527: /* tickattr: MAJOR GRID onoff  */
#line 1846 "gram.y"
                           { g[curg].t[naxis].t_gridflag = (yyvsp[0].pset); }
#line 8558 "y.tab.c"
    break;

  case 528: /* tickattr: MINOR GRID onoff  */
#line 1847 "gram.y"
                           { g[curg].t[naxis].t_mgridflag = (yyvsp[0].pset); }
#line 8564 "y.tab.c"
    break;

  case 529: /* tickattr: OP opchoice  */
#line 1848 "gram.y"
                      { g[curg].t[naxis].t_op = (yyvsp[0].pset); }
#line 8570 "y.tab.c"
    break;

  case 530: /* tickattr: TYPE AUTO  */
#line 1849 "gram.y"
                    { g[curg].t[naxis].t_type = AUTO; }
#line 8576 "y.tab.c"
    break;

  case 531: /* tickattr: TYPE SPEC  */
#line 1850 "gram.y"
                    { g[curg].t[naxis].t_type = SPEC; }
#line 8582 "y.tab.c"
    break;

  case 532: /* tickattr: SPEC NUMBER  */
#line 1851 "gram.y"
                      { g[curg].t[naxis].t_spec = (int) (yyvsp[0].val); }
#line 8588 "y.tab.c"
    break;

  case 533: /* tickattr: NUMBER ',' expr  */
#line 1852 "gram.y"
                          { g[curg].t[naxis].t_specloc[(int) (yyvsp[-2].val)] = (yyvsp[0].val); }
#line 8594 "y.tab.c"
    break;

  case 534: /* ticklabeldesc: ticklabelattr  */
#line 1856 "gram.y"
                      {}
#line 8600 "y.tab.c"
    break;

  case 535: /* ticklabeldesc: ticklabeldesc ticklabelattr  */
#line 1857 "gram.y"
                                       {}
#line 8606 "y.tab.c"
    break;

  case 536: /* ticklabelattr: onoff  */
#line 1861 "gram.y"
              { g[curg].t[naxis].tl_flag = (yyvsp[0].pset); }
#line 8612 "y.tab.c"
    break;

  case 537: /* ticklabelattr: TYPE AUTO  */
#line 1862 "gram.y"
                    { g[curg].t[naxis].tl_type = AUTO; }
#line 8618 "y.tab.c"
    break;

  case 538: /* ticklabelattr: TYPE SPEC  */
#line 1863 "gram.y"
                    { g[curg].t[naxis].tl_type = SPEC; }
#line 8624 "y.tab.c"
    break;

  case 539: /* ticklabelattr: PREC NUMBER  */
#line 1864 "gram.y"
                      { g[curg].t[naxis].tl_prec = (int) (yyvsp[0].val); }
#line 8630 "y.tab.c"
    break;

  case 540: /* ticklabelattr: FORMAT formatchoice  */
#line 1865 "gram.y"
                              { g[curg].t[naxis].tl_format = (yyvsp[0].pset); }
#line 8636 "y.tab.c"
    break;

  case 541: /* ticklabelattr: FORMAT NUMBER  */
#line 1866 "gram.y"
                        { g[curg].t[naxis].tl_format = (yyvsp[0].val); }
#line 8642 "y.tab.c"
    break;

  case 542: /* ticklabelattr: APPEND CHRSTR  */
#line 1867 "gram.y"
                        { strcpy(g[curg].t[naxis].tl_appstr, (char *) (yyvsp[0].str)); }
#line 8648 "y.tab.c"
    break;

  case 543: /* ticklabelattr: PREPEND CHRSTR  */
#line 1868 "gram.y"
                         { strcpy(g[curg].t[naxis].tl_prestr, (char *) (yyvsp[0].str)); }
#line 8654 "y.tab.c"
    break;

  case 544: /* ticklabelattr: LAYOUT HORIZONTAL  */
#line 1869 "gram.y"
                            { g[curg].t[naxis].tl_layout = HORIZONTAL; }
#line 8660 "y.tab.c"
    break;

  case 545: /* ticklabelattr: LAYOUT VERTICAL  */
#line 1870 "gram.y"
                          { g[curg].t[naxis].tl_layout = VERTICAL; }
#line 8666 "y.tab.c"
    break;

  case 546: /* ticklabelattr: LAYOUT SPEC  */
#line 1871 "gram.y"
                      { g[curg].t[naxis].tl_layout = SPEC; }
#line 8672 "y.tab.c"
    break;

  case 547: /* ticklabelattr: ANGLE NUMBER  */
#line 1872 "gram.y"
                       { g[curg].t[naxis].tl_angle = (int) (yyvsp[0].val); }
#line 8678 "y.tab.c"
    break;

  case 548: /* ticklabelattr: JUST justchoice  */
#line 1873 "gram.y"
                          { g[curg].t[naxis].tl_just = (int) (yyvsp[0].pset); }
#line 8684 "y.tab.c"
    break;

  case 549: /* ticklabelattr: SKIP NUMBER  */
#line 1874 "gram.y"
                      { g[curg].t[naxis].tl_skip = (int) (yyvsp[0].val); }
#line 8690 "y.tab.c"
    break;

  case 550: /* ticklabelattr: STAGGER NUMBER  */
#line 1875 "gram.y"
                         { g[curg].t[naxis].tl_staggered = (int) (yyvsp[0].val); }
#line 8696 "y.tab.c"
    break;

  case 551: /* ticklabelattr: OP opchoice  */
#line 1876 "gram.y"
                      { g[curg].t[naxis].tl_op = (yyvsp[0].pset); }
#line 8702 "y.tab.c"
    break;

  case 552: /* ticklabelattr: SIGN signchoice  */
#line 1877 "gram.y"
                          { g[curg].t[naxis].tl_sign = (yyvsp[0].pset); }
#line 8708 "y.tab.c"
    break;

  case 553: /* ticklabelattr: START expr  */
#line 1878 "gram.y"
                     { g[curg].t[naxis].tl_start = (yyvsp[0].val); }
#line 8714 "y.tab.c"
    break;

  case 554: /* ticklabelattr: STOP expr  */
#line 1879 "gram.y"
                    { g[curg].t[naxis].tl_stop = (yyvsp[0].val); }
#line 8720 "y.tab.c"
    break;

  case 555: /* ticklabelattr: START TYPE SPEC  */
#line 1880 "gram.y"
                          { g[curg].t[naxis].tl_starttype = (int) (yyvsp[0].pset); }
#line 8726 "y.tab.c"
    break;

  case 556: /* ticklabelattr: START TYPE AUTO  */
#line 1881 "gram.y"
                          { g[curg].t[naxis].tl_starttype = (int) (yyvsp[0].pset); }
#line 8732 "y.tab.c"
    break;

  case 557: /* ticklabelattr: STOP TYPE SPEC  */
#line 1882 "gram.y"
                         { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[0].pset); }
#line 8738 "y.tab.c"
    break;

  case 558: /* ticklabelattr: STOP TYPE AUTO  */
#line 1883 "gram.y"
                         { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[0].pset); }
#line 8744 "y.tab.c"
    break;

  case 559: /* ticklabelattr: VGAP NUMBER  */
#line 1884 "gram.y"
                      { g[curg].t[naxis].tl_vgap = (yyvsp[0].val); }
#line 8750 "y.tab.c"
    break;

  case 560: /* ticklabelattr: HGAP NUMBER  */
#line 1885 "gram.y"
                      { g[curg].t[naxis].tl_hgap = (yyvsp[0].val); }
#line 8756 "y.tab.c"
    break;

  case 561: /* ticklabelattr: CHAR SIZE NUMBER  */
#line 1886 "gram.y"
                           { g[curg].t[naxis].tl_charsize = (yyvsp[0].val); }
#line 8762 "y.tab.c"
    break;

  case 562: /* ticklabelattr: FONTP NUMBER  */
#line 1887 "gram.y"
                       { g[curg].t[naxis].tl_font = (int) (yyvsp[0].val); }
#line 8768 "y.tab.c"
    break;

  case 563: /* ticklabelattr: COLOR NUMBER  */
#line 1888 "gram.y"
                       { g[curg].t[naxis].tl_color = (int) (yyvsp[0].val); }
#line 8774 "y.tab.c"
    break;

  case 564: /* ticklabelattr: LINEWIDTH NUMBER  */
#line 1889 "gram.y"
                           { g[curg].t[naxis].tl_linew = (int) (yyvsp[0].val); }
#line 8780 "y.tab.c"
    break;

  case 565: /* ticklabelattr: NUMBER ',' CHRSTR  */
#line 1890 "gram.y"
                            { set_plotstr_string(&g[curg].t[naxis].t_speclab[(int) (yyvsp[-2].val)], (char *) (yyvsp[0].str)); }
#line 8786 "y.tab.c"
    break;

  case 566: /* axislabeldesc: CHRSTR  */
#line 1894 "gram.y"
               { set_plotstr_string(&g[curg].t[naxis].label, (char *) (yyvsp[0].str)); }
#line 8792 "y.tab.c"
    break;

  case 567: /* axislabeldesc: LAYOUT PERP  */
#line 1895 "gram.y"
                      { g[curg].t[naxis].label_layout = PERP; }
#line 8798 "y.tab.c"
    break;

  case 568: /* axislabeldesc: LAYOUT PARA  */
#line 1896 "gram.y"
                      { g[curg].t[naxis].label_layout = PARA; }
#line 8804 "y.tab.c"
    break;

  case 569: /* axislabeldesc: PLACE AUTO  */
#line 1897 "gram.y"
                     { g[curg].t[naxis].label_place = (yyvsp[0].pset); }
#line 8810 "y.tab.c"
    break;

  case 570: /* axislabeldesc: PLACE SPEC  */
#line 1898 "gram.y"
                     { g[curg].t[naxis].label_place = (yyvsp[0].pset); }
#line 8816 "y.tab.c"
    break;

  case 571: /* axislabeldesc: PLACE NUMBER ',' NUMBER  */
#line 1899 "gram.y"
                                  {
	    g[curg].t[naxis].label.x = (yyvsp[-2].val);
	    g[curg].t[naxis].label.y = (yyvsp[0].val);
	}
#line 8825 "y.tab.c"
    break;

  case 572: /* axislabeldesc: JUST justchoice  */
#line 1903 "gram.y"
                          { g[curg].t[naxis].label.just = (int) (yyvsp[0].pset); }
#line 8831 "y.tab.c"
    break;

  case 573: /* axislabeldesc: CHAR SIZE NUMBER  */
#line 1904 "gram.y"
                           { g[curg].t[naxis].label.charsize = (yyvsp[0].val); }
#line 8837 "y.tab.c"
    break;

  case 574: /* axislabeldesc: FONTP NUMBER  */
#line 1905 "gram.y"
                       { g[curg].t[naxis].label.font = (int) (yyvsp[0].val); }
#line 8843 "y.tab.c"
    break;

  case 575: /* axislabeldesc: COLOR NUMBER  */
#line 1906 "gram.y"
                       { g[curg].t[naxis].label.color = (int) (yyvsp[0].val); }
#line 8849 "y.tab.c"
    break;

  case 576: /* axislabeldesc: LINEWIDTH NUMBER  */
#line 1907 "gram.y"
                           { g[curg].t[naxis].label.linew = (int) (yyvsp[0].val); }
#line 8855 "y.tab.c"
    break;

  case 577: /* axisbardesc: onoff  */
#line 1911 "gram.y"
              { g[curg].t[naxis].t_drawbar = (yyvsp[0].pset); }
#line 8861 "y.tab.c"
    break;

  case 578: /* axisbardesc: COLOR NUMBER  */
#line 1912 "gram.y"
                       { g[curg].t[naxis].t_drawbarcolor = (int) (yyvsp[0].val); }
#line 8867 "y.tab.c"
    break;

  case 579: /* axisbardesc: LINESTYLE NUMBER  */
#line 1913 "gram.y"
                           { g[curg].t[naxis].t_drawbarlines = (int) (yyvsp[0].val); }
#line 8873 "y.tab.c"
    break;

  case 580: /* axisbardesc: LINEWIDTH NUMBER  */
#line 1914 "gram.y"
                           { g[curg].t[naxis].t_drawbarlinew = (int) (yyvsp[0].val); }
#line 8879 "y.tab.c"
    break;

  case 581: /* sourcetype: DISK  */
#line 1926 "gram.y"
             { (yyval.pset) = DISK; }
#line 8885 "y.tab.c"
    break;

  case 582: /* sourcetype: PIPE  */
#line 1927 "gram.y"
               { (yyval.pset) = PIPE; }
#line 8891 "y.tab.c"
    break;

  case 583: /* justchoice: RIGHT  */
#line 1931 "gram.y"
              { (yyval.pset) = RIGHT; }
#line 8897 "y.tab.c"
    break;

  case 584: /* justchoice: LEFT  */
#line 1932 "gram.y"
               { (yyval.pset) = LEFT; }
#line 8903 "y.tab.c"
    break;

  case 585: /* justchoice: CENTER  */
#line 1933 "gram.y"
                 { (yyval.pset) = CENTER; }
#line 8909 "y.tab.c"
    break;

  case 586: /* graphtype: XY  */
#line 1942 "gram.y"
           { (yyval.pset) = (yyvsp[0].pset); }
#line 8915 "y.tab.c"
    break;

  case 587: /* graphtype: LOGX  */
#line 1943 "gram.y"
               { (yyval.pset) = (yyvsp[0].pset); }
#line 8921 "y.tab.c"
    break;

  case 588: /* graphtype: LOGY  */
#line 1944 "gram.y"
               { (yyval.pset) = (yyvsp[0].pset); }
#line 8927 "y.tab.c"
    break;

  case 589: /* graphtype: LOGXY  */
#line 1945 "gram.y"
                { (yyval.pset) = (yyvsp[0].pset); }
#line 8933 "y.tab.c"
    break;

  case 590: /* graphtype: XYFIXED  */
#line 1946 "gram.y"
                  { (yyval.pset) = XYFIXED; }
#line 8939 "y.tab.c"
    break;

  case 591: /* inoutchoice: IN  */
#line 1950 "gram.y"
           { (yyval.pset) = IN; }
#line 8945 "y.tab.c"
    break;

  case 592: /* inoutchoice: OUT  */
#line 1951 "gram.y"
              { (yyval.pset) = OUT; }
#line 8951 "y.tab.c"
    break;

  case 593: /* inoutchoice: BOTH  */
#line 1952 "gram.y"
               { (yyval.pset) = BOTH; }
#line 8957 "y.tab.c"
    break;

  case 594: /* signchoice: NORMAL  */
#line 1956 "gram.y"
               { (yyval.pset) = NORMAL; }
#line 8963 "y.tab.c"
    break;

  case 595: /* signchoice: ABSOLUTE  */
#line 1957 "gram.y"
                   { (yyval.pset) = ABSOLUTE; }
#line 8969 "y.tab.c"
    break;

  case 596: /* signchoice: NEGATE  */
#line 1958 "gram.y"
                 { (yyval.pset) = NEGATE; }
#line 8975 "y.tab.c"
    break;

  case 597: /* formatchoice: DECIMAL  */
#line 1971 "gram.y"
                { (yyval.pset) = DECIMAL; }
#line 8981 "y.tab.c"
    break;

  case 598: /* formatchoice: EXPONENTIAL  */
#line 1972 "gram.y"
                      { (yyval.pset) = EXPONENTIAL; }
#line 8987 "y.tab.c"
    break;

  case 599: /* formatchoice: POWER  */
#line 1973 "gram.y"
                { (yyval.pset) = POWER; }
#line 8993 "y.tab.c"
    break;

  case 600: /* formatchoice: GENERAL  */
#line 1974 "gram.y"
                  { (yyval.pset) = GENERAL; }
#line 8999 "y.tab.c"
    break;

  case 601: /* formatchoice: DDMMYY  */
#line 1975 "gram.y"
                 { (yyval.pset) = DDMMYY; }
#line 9005 "y.tab.c"
    break;

  case 602: /* formatchoice: MMDDYY  */
#line 1976 "gram.y"
                 { (yyval.pset) = MMDDYY; }
#line 9011 "y.tab.c"
    break;

  case 603: /* formatchoice: MMYY  */
#line 1977 "gram.y"
               { (yyval.pset) = MMYY; }
#line 9017 "y.tab.c"
    break;

  case 604: /* formatchoice: MMDD  */
#line 1978 "gram.y"
               { (yyval.pset) = MMDD; }
#line 9023 "y.tab.c"
    break;

  case 605: /* formatchoice: MONTHDAY  */
#line 1979 "gram.y"
                   { (yyval.pset) = MONTHDAY; }
#line 9029 "y.tab.c"
    break;

  case 606: /* formatchoice: DAYMONTH  */
#line 1980 "gram.y"
                   { (yyval.pset) = DAYMONTH; }
#line 9035 "y.tab.c"
    break;

  case 607: /* formatchoice: DDMONTHSYYHHMMSS  */
#line 1981 "gram.y"
                           { (yyval.pset) = DDMONTHSYYHHMMSS; }
#line 9041 "y.tab.c"
    break;

  case 608: /* formatchoice: MONTHS  */
#line 1982 "gram.y"
                 { (yyval.pset) = MONTHS; }
#line 9047 "y.tab.c"
    break;

  case 609: /* formatchoice: MONTHL  */
#line 1983 "gram.y"
                 { (yyval.pset) = MONTHL; }
#line 9053 "y.tab.c"
    break;

  case 610: /* formatchoice: DAYOFWEEKS  */
#line 1984 "gram.y"
                     { (yyval.pset) = DAYOFWEEKS; }
#line 9059 "y.tab.c"
    break;

  case 611: /* formatchoice: DAYOFWEEKL  */
#line 1985 "gram.y"
                     { (yyval.pset) = DAYOFWEEKL; }
#line 9065 "y.tab.c"
    break;

  case 612: /* formatchoice: DAYOFYEAR  */
#line 1986 "gram.y"
                    { (yyval.pset) = DAYOFYEAR; }
#line 9071 "y.tab.c"
    break;

  case 613: /* formatchoice: HMS  */
#line 1987 "gram.y"
              { (yyval.pset) = HMS; }
#line 9077 "y.tab.c"
    break;

  case 614: /* formatchoice: MMDDHMS  */
#line 1988 "gram.y"
                  { (yyval.pset) = MMDDHMS; }
#line 9083 "y.tab.c"
    break;

  case 615: /* formatchoice: MMDDYYHMS  */
#line 1989 "gram.y"
                    { (yyval.pset) = MMDDYYHMS; }
#line 9089 "y.tab.c"
    break;

  case 616: /* formatchoice: DEGREESLON  */
#line 1990 "gram.y"
                     { (yyval.pset) = DEGREESLON; }
#line 9095 "y.tab.c"
    break;

  case 617: /* formatchoice: DEGREESMMLON  */
#line 1991 "gram.y"
                       { (yyval.pset) = DEGREESMMLON; }
#line 9101 "y.tab.c"
    break;

  case 618: /* formatchoice: DEGREESMMSSLON  */
#line 1992 "gram.y"
                         { (yyval.pset) = DEGREESMMSSLON; }
#line 9107 "y.tab.c"
    break;

  case 619: /* formatchoice: MMSSLON  */
#line 1993 "gram.y"
                  { (yyval.pset) = MMSSLON; }
#line 9113 "y.tab.c"
    break;

  case 620: /* formatchoice: DEGREESLAT  */
#line 1994 "gram.y"
                     { (yyval.pset) = DEGREESLAT; }
#line 9119 "y.tab.c"
    break;

  case 621: /* formatchoice: DEGREESMMLAT  */
#line 1995 "gram.y"
                       { (yyval.pset) = DEGREESMMLAT; }
#line 9125 "y.tab.c"
    break;

  case 622: /* formatchoice: DEGREESMMSSLAT  */
#line 1996 "gram.y"
                         { (yyval.pset) = DEGREESMMSSLAT; }
#line 9131 "y.tab.c"
    break;

  case 623: /* formatchoice: MMSSLAT  */
#line 1997 "gram.y"
                  { (yyval.pset) = MMSSLAT; }
#line 9137 "y.tab.c"
    break;

  case 624: /* horv: HORIZONTAL  */
#line 2000 "gram.y"
                 { (yyval.pset) = (yyvsp[0].pset); }
#line 9143 "y.tab.c"
    break;

  case 625: /* horv: VERTICAL  */
#line 2001 "gram.y"
                 { (yyval.pset) = (yyvsp[0].pset); }
#line 9149 "y.tab.c"
    break;

  case 626: /* flowonoff: OFF  */
#line 2005 "gram.y"
            { (yyval.pset) = (yyvsp[0].pset); }
#line 9155 "y.tab.c"
    break;

  case 627: /* flowonoff: ON  */
#line 2006 "gram.y"
             { (yyval.pset) = (yyvsp[0].pset); }
#line 9161 "y.tab.c"
    break;

  case 628: /* flowonoff: NODES  */
#line 2007 "gram.y"
                { (yyval.pset) = (yyvsp[0].pset); }
#line 9167 "y.tab.c"
    break;

  case 629: /* flowonoff: CENTER  */
#line 2008 "gram.y"
                 { (yyval.pset) = (yyvsp[0].pset); }
#line 9173 "y.tab.c"
    break;

  case 630: /* onoff: ON  */
#line 2012 "gram.y"
           { (yyval.pset) = ON; }
#line 9179 "y.tab.c"
    break;

  case 631: /* onoff: OFF  */
#line 2013 "gram.y"
              { (yyval.pset) = OFF; }
#line 9185 "y.tab.c"
    break;

  case 632: /* worldview: WORLD  */
#line 2016 "gram.y"
                 { (yyval.pset) = WORLD; }
#line 9191 "y.tab.c"
    break;

  case 633: /* worldview: VIEW  */
#line 2017 "gram.y"
               { (yyval.pset) = VIEW; }
#line 9197 "y.tab.c"
    break;

  case 634: /* torf: TRUEP  */
#line 2028 "gram.y"
              { (yyval.pset) = TRUEP; }
#line 9203 "y.tab.c"
    break;

  case 635: /* torf: FALSEP  */
#line 2029 "gram.y"
                 { (yyval.pset) = FALSEP; }
#line 9209 "y.tab.c"
    break;

  case 636: /* filltype: PATTERN  */
#line 2032 "gram.y"
                  { (yyval.pset) = PATTERN; }
#line 9215 "y.tab.c"
    break;

  case 637: /* filltype: COLOR  */
#line 2033 "gram.y"
                { (yyval.pset) = COLOR; }
#line 9221 "y.tab.c"
    break;

  case 638: /* filltype: NONE  */
#line 2034 "gram.y"
               { (yyval.pset) = NONE; }
#line 9227 "y.tab.c"
    break;

  case 639: /* opchoice: ABOVE  */
#line 2037 "gram.y"
                { (yyval.pset) = ABOVE; }
#line 9233 "y.tab.c"
    break;

  case 640: /* opchoice: BELOW  */
#line 2038 "gram.y"
                { (yyval.pset) = BELOW; }
#line 9239 "y.tab.c"
    break;

  case 641: /* opchoice: LEFT  */
#line 2039 "gram.y"
               { (yyval.pset) = LEFT; }
#line 9245 "y.tab.c"
    break;

  case 642: /* opchoice: RIGHT  */
#line 2040 "gram.y"
                { (yyval.pset) = RIGHT; }
#line 9251 "y.tab.c"
    break;

  case 643: /* opchoice: TOP  */
#line 2041 "gram.y"
              { (yyval.pset) = TOP; }
#line 9257 "y.tab.c"
    break;

  case 644: /* opchoice: BOTTOM  */
#line 2042 "gram.y"
                 { (yyval.pset) = BOTTOM ; }
#line 9263 "y.tab.c"
    break;

  case 645: /* opchoice: BOTH  */
#line 2043 "gram.y"
               { (yyval.pset) = BOTH ; }
#line 9269 "y.tab.c"
    break;

  case 646: /* units: MM  */
#line 2047 "gram.y"
           { (yyval.pset) = MM; }
#line 9275 "y.tab.c"
    break;

  case 647: /* units: CM  */
#line 2048 "gram.y"
             { (yyval.pset) = CM; }
#line 9281 "y.tab.c"
    break;

  case 648: /* units: M  */
#line 2049 "gram.y"
            { (yyval.pset) = M; }
#line 9287 "y.tab.c"
    break;

  case 649: /* units: KM  */
#line 2050 "gram.y"
             { (yyval.pset) = KM; }
#line 9293 "y.tab.c"
    break;

  case 650: /* asgn: VAR '[' expr ']' '=' expr  */
#line 2065 "gram.y"
        {
	    int itmp = (int) (yyvsp[-3].val) - 1;
	    if (itmp >= ls) {
		yyerror("subscript out of range");
		return 1;
	    } else {
		(yyvsp[-5].ptr)[itmp] = (yyvsp[0].val);
		result = (yyvsp[0].val);
	    }
	}
#line 9308 "y.tab.c"
    break;

  case 651: /* vasgn: VAR '=' vexpr  */
#line 2079 "gram.y"
        {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[-2].ptr)[i] = (yyvsp[0].ptr)[i];
	    }
	    result = (yyvsp[0].ptr)[0];
	}
#line 9320 "y.tab.c"
    break;

  case 652: /* vasgn: VAR '=' expr  */
#line 2087 "gram.y"
        {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[-2].ptr)[i] = (yyvsp[0].val);
	    }
	    result = (yyvsp[0].val);
	}
#line 9332 "y.tab.c"
    break;

  case 653: /* vexpr: VAR  */
#line 2098 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].ptr)[i];
	    }
	}
#line 9345 "y.tab.c"
    break;

  case 654: /* vexpr: expr  */
#line 2107 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].val);
	    }
	}
#line 9358 "y.tab.c"
    break;

  case 655: /* vexpr: expr '+' expr  */
#line 2116 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) + (yyvsp[0].val);
	    }
	}
#line 9371 "y.tab.c"
    break;

  case 656: /* vexpr: vexpr '+' vexpr  */
#line 2125 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] + (yyvsp[0].ptr)[i];
	    }
	}
#line 9384 "y.tab.c"
    break;

  case 657: /* vexpr: expr '+' vexpr  */
#line 2134 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) + (yyvsp[0].ptr)[i];
	    }
	}
#line 9397 "y.tab.c"
    break;

  case 658: /* vexpr: vexpr '+' expr  */
#line 2143 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] + (yyvsp[0].val);
	    }
	}
#line 9410 "y.tab.c"
    break;

  case 659: /* vexpr: expr '-' expr  */
#line 2152 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) - (yyvsp[0].val);
	    }
	}
#line 9423 "y.tab.c"
    break;

  case 660: /* vexpr: vexpr '-' vexpr  */
#line 2161 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] - (yyvsp[0].ptr)[i];
	    }
	}
#line 9436 "y.tab.c"
    break;

  case 661: /* vexpr: expr '-' vexpr  */
#line 2170 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) - (yyvsp[0].ptr)[i];
	    }
	}
#line 9449 "y.tab.c"
    break;

  case 662: /* vexpr: vexpr '-' expr  */
#line 2179 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] - (yyvsp[0].val);
	    }
	}
#line 9462 "y.tab.c"
    break;

  case 663: /* vexpr: expr '*' expr  */
#line 2188 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) * (yyvsp[0].val);
	    }
	}
#line 9475 "y.tab.c"
    break;

  case 664: /* vexpr: vexpr '*' vexpr  */
#line 2197 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] * (yyvsp[0].ptr)[i];
	    }
	}
#line 9488 "y.tab.c"
    break;

  case 665: /* vexpr: expr '*' vexpr  */
#line 2206 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) * (yyvsp[0].ptr)[i];
	    }
	}
#line 9501 "y.tab.c"
    break;

  case 666: /* vexpr: vexpr '*' expr  */
#line 2215 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] * (yyvsp[0].val);
	    }
	}
#line 9514 "y.tab.c"
    break;

  case 667: /* vexpr: expr '/' expr  */
#line 2224 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[0].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) / (yyvsp[0].val);
	    }
	}
#line 9531 "y.tab.c"
    break;

  case 668: /* vexpr: vexpr '/' vexpr  */
#line 2237 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[0].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] / (yyvsp[0].ptr)[i];
	    }
	}
#line 9550 "y.tab.c"
    break;

  case 669: /* vexpr: expr '/' vexpr  */
#line 2252 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[0].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].val) / (yyvsp[0].ptr)[i];
	    }
	}
#line 9569 "y.tab.c"
    break;

  case 670: /* vexpr: vexpr '/' expr  */
#line 2267 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[0].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] / (yyvsp[0].val);
	    }
	}
#line 9586 "y.tab.c"
    break;

  case 671: /* vexpr: expr '^' expr  */
#line 2280 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].val), (yyvsp[0].val));
	    }
	}
#line 9599 "y.tab.c"
    break;

  case 672: /* vexpr: expr '^' vexpr  */
#line 2289 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].val), (yyvsp[0].ptr)[i]);
	    }
	}
#line 9612 "y.tab.c"
    break;

  case 673: /* vexpr: vexpr '^' expr  */
#line 2298 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].ptr)[i], (yyvsp[0].val));
	    }
	}
#line 9625 "y.tab.c"
    break;

  case 674: /* vexpr: vexpr '^' vexpr  */
#line 2307 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[-2].ptr)[i], (yyvsp[0].ptr)[i]);
	    }
	}
#line 9638 "y.tab.c"
    break;

  case 675: /* vexpr: ABS '(' expr ')'  */
#line 2316 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[-1].val));
	    }
	}
#line 9651 "y.tab.c"
    break;

  case 676: /* vexpr: ABS '(' vexpr ')'  */
#line 2325 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9664 "y.tab.c"
    break;

  case 677: /* vexpr: ACOS '(' vexpr ')'  */
#line 2334 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = acos((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9677 "y.tab.c"
    break;

  case 678: /* vexpr: ASIN '(' vexpr ')'  */
#line 2343 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = asin((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9690 "y.tab.c"
    break;

  case 679: /* vexpr: ATAN '(' vexpr ')'  */
#line 2352 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9703 "y.tab.c"
    break;

  case 680: /* vexpr: ATAN2 '(' vexpr ',' vexpr ')'  */
#line 2361 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan2((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9716 "y.tab.c"
    break;

  case 681: /* vexpr: CEIL '(' vexpr ')'  */
#line 2370 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ceil((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9729 "y.tab.c"
    break;

  case 682: /* vexpr: COS '(' vexpr ')'  */
#line 2379 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = cos((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9742 "y.tab.c"
    break;

  case 683: /* vexpr: DEG  */
#line 2388 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] *= M_PI / 180.0;
	    }
	}
#line 9755 "y.tab.c"
    break;

  case 684: /* vexpr: DX  */
#line 2397 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = xx[i];
	    }
	}
#line 9768 "y.tab.c"
    break;

  case 685: /* vexpr: DY  */
#line 2406 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = yy[i];
	    }
	}
#line 9781 "y.tab.c"
    break;

  case 686: /* vexpr: ERF '(' vexpr ')'  */
#line 2415 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erf((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9794 "y.tab.c"
    break;

  case 687: /* vexpr: ERFC '(' vexpr ')'  */
#line 2424 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erfc((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9807 "y.tab.c"
    break;

  case 688: /* vexpr: EXP '(' vexpr ')'  */
#line 2433 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = exp((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9820 "y.tab.c"
    break;

  case 689: /* vexpr: FLOOR '(' vexpr ')'  */
#line 2442 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = floor((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9833 "y.tab.c"
    break;

  case 690: /* vexpr: HYPOT '(' vexpr ',' vexpr ')'  */
#line 2451 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9846 "y.tab.c"
    break;

  case 691: /* vexpr: HYPOT '(' expr ',' vexpr ')'  */
#line 2460 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].val), (yyvsp[-1].ptr)[i]);
	    }
	}
#line 9859 "y.tab.c"
    break;

  case 692: /* vexpr: HYPOT '(' vexpr ',' expr ')'  */
#line 2469 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].ptr)[i], (yyvsp[-1].val));
	    }
	}
#line 9872 "y.tab.c"
    break;

  case 693: /* vexpr: HYPOT '(' expr ',' expr ')'  */
#line 2478 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[-3].val), (yyvsp[-1].val));
	    }
	}
#line 9885 "y.tab.c"
    break;

  case 694: /* vexpr: INDEX  */
#line 2487 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = i + 1;
	    }
	}
#line 9898 "y.tab.c"
    break;

  case 695: /* vexpr: SETNO  */
#line 2496 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[0].func);
	    }
	}
#line 9911 "y.tab.c"
    break;

  case 696: /* vexpr: INT '(' vexpr ')'  */
#line 2505 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (int) (yyvsp[-1].ptr)[i];
	    }
	}
#line 9924 "y.tab.c"
    break;

  case 697: /* vexpr: IRAND '(' NUMBER ')'  */
#line 2514 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lrand48() % (long) ((yyvsp[-1].val));
	    }
	}
#line 9937 "y.tab.c"
    break;

  case 698: /* vexpr: LGAMMA '(' vexpr ')'  */
#line 2523 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lgamma((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9950 "y.tab.c"
    break;

  case 699: /* vexpr: LN '(' vexpr ')'  */
#line 2532 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9963 "y.tab.c"
    break;

  case 700: /* vexpr: LOG '(' vexpr ')'  */
#line 2541 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log10((yyvsp[-1].ptr)[i]);
	    }
	}
#line 9976 "y.tab.c"
    break;

  case 701: /* vexpr: LOGISTIC '(' vexpr ',' expr ',' expr ')'  */
#line 2550 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = 1.0 / (1.0 + exp(-((yyvsp[-5].ptr)[i] - (yyvsp[-3].val))/ (yyvsp[-1].val)));
	    }
	}
#line 9989 "y.tab.c"
    break;

  case 702: /* vexpr: MAXP '(' vexpr ',' vexpr ')'  */
#line 2559 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-3].ptr)[i] >= (yyvsp[-1].ptr)[i] ? (yyvsp[-3].ptr)[i] : (yyvsp[-1].ptr)[i];
	    }
	}
#line 10002 "y.tab.c"
    break;

  case 703: /* vexpr: MINP '(' vexpr ',' vexpr ')'  */
#line 2568 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-3].ptr)[i] <= (yyvsp[-1].ptr)[i] ? (yyvsp[-3].ptr)[i] : (yyvsp[-1].ptr)[i];
	    }
	}
#line 10015 "y.tab.c"
    break;

  case 704: /* vexpr: MOD '(' vexpr ',' vexpr ')'  */
#line 2577 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fmod((yyvsp[-3].ptr)[i], (yyvsp[-1].ptr)[i]);
	    }
	}
#line 10028 "y.tab.c"
    break;

  case 705: /* vexpr: NORM '(' vexpr ')'  */
#line 2586 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fx((yyvsp[-1].ptr)[i]);
	    }
	}
#line 10041 "y.tab.c"
    break;

  case 706: /* vexpr: NORMP '(' vexpr ')'  */
#line 2595 "gram.y"
        {
	    int i;
	    double tmp;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = normp((yyvsp[-1].ptr)[i], &tmp);
	    }
	}
#line 10055 "y.tab.c"
    break;

  case 707: /* vexpr: PI  */
#line 2605 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI;
	    }
	}
#line 10068 "y.tab.c"
    break;

  case 708: /* vexpr: RAD  */
#line 2614 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI / 180.0;
	    }
	}
#line 10081 "y.tab.c"
    break;

  case 709: /* vexpr: RAND  */
#line 2623 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (double) drand48();
	    }
	}
#line 10094 "y.tab.c"
    break;

  case 710: /* vexpr: SIN '(' vexpr ')'  */
#line 2632 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sin((yyvsp[-1].ptr)[i]);
	    }
	}
#line 10107 "y.tab.c"
    break;

  case 711: /* vexpr: SQR '(' vexpr ')'  */
#line 2641 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-1].ptr)[i] * (yyvsp[-1].ptr)[i];
	    }
	}
#line 10120 "y.tab.c"
    break;

  case 712: /* vexpr: SQRT '(' vexpr ')'  */
#line 2650 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sqrt((yyvsp[-1].ptr)[i]);
	    }
	}
#line 10133 "y.tab.c"
    break;

  case 713: /* vexpr: TAN '(' vexpr ')'  */
#line 2659 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = tan((yyvsp[-1].ptr)[i]);
	    }
	}
#line 10146 "y.tab.c"
    break;

  case 714: /* vexpr: vexpr GT vexpr  */
#line 2668 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] > (yyvsp[0].ptr)[i];
	    }
	}
#line 10159 "y.tab.c"
    break;

  case 715: /* vexpr: vexpr LT vexpr  */
#line 2677 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] < (yyvsp[0].ptr)[i];
	    }
	}
#line 10172 "y.tab.c"
    break;

  case 716: /* vexpr: vexpr LE vexpr  */
#line 2686 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] <= (yyvsp[0].ptr)[i];
	    }
	}
#line 10185 "y.tab.c"
    break;

  case 717: /* vexpr: vexpr GE vexpr  */
#line 2695 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] >= (yyvsp[0].ptr)[i];
	    }
	}
#line 10198 "y.tab.c"
    break;

  case 718: /* vexpr: vexpr EQ vexpr  */
#line 2704 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] == (yyvsp[0].ptr)[i];
	    }
	}
#line 10211 "y.tab.c"
    break;

  case 719: /* vexpr: vexpr NE vexpr  */
#line 2713 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] != (yyvsp[0].ptr)[i];
	    }
	}
#line 10224 "y.tab.c"
    break;

  case 720: /* vexpr: vexpr AND vexpr  */
#line 2722 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] && (yyvsp[0].ptr)[i];
	    }
	}
#line 10237 "y.tab.c"
    break;

  case 721: /* vexpr: vexpr OR vexpr  */
#line 2731 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-2].ptr)[i] || (yyvsp[0].ptr)[i];
	    }
	}
#line 10250 "y.tab.c"
    break;

  case 722: /* vexpr: NOT vexpr  */
#line 2740 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = !((yyvsp[0].ptr)[i]);
	    }
	}
#line 10263 "y.tab.c"
    break;

  case 723: /* vexpr: '(' vexpr ')'  */
#line 2749 "gram.y"
        {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[-1].ptr)[i];
	    }
	}
#line 10276 "y.tab.c"
    break;

  case 724: /* vexpr: '-' vexpr  */
#line 2757 "gram.y"
                                 {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = -(yyvsp[0].ptr)[i];
	    }
	}
#line 10289 "y.tab.c"
    break;

  case 726: /* expr: FITPARM  */
#line 2768 "gram.y"
                   {
	    (yyval.val) = (yyvsp[0].val);
	}
#line 10297 "y.tab.c"
    break;

  case 727: /* expr: VAR '[' expr ']'  */
#line 2771 "gram.y"
                            {
	    (yyval.val) = (yyvsp[-3].ptr)[(int) (yyvsp[-1].val)];
	}
#line 10305 "y.tab.c"
    break;

  case 728: /* expr: expr '+' expr  */
#line 2774 "gram.y"
                        {
	    (yyval.val) = (yyvsp[-2].val) + (yyvsp[0].val);
	}
#line 10313 "y.tab.c"
    break;

  case 729: /* expr: expr '-' expr  */
#line 2777 "gram.y"
                        {
	    (yyval.val) = (yyvsp[-2].val) - (yyvsp[0].val);
	}
#line 10321 "y.tab.c"
    break;

  case 730: /* expr: expr '*' expr  */
#line 2780 "gram.y"
                        {
	    (yyval.val) = (yyvsp[-2].val) * (yyvsp[0].val);
	}
#line 10329 "y.tab.c"
    break;

  case 731: /* expr: expr '/' expr  */
#line 2784 "gram.y"
        {
	    if ((yyvsp[0].val) != 0.0) {
		(yyval.val) = (yyvsp[-2].val) / (yyvsp[0].val);
	    } else {
		yyerror("Divide by Zero");
		return 1;
	    }
	}
#line 10342 "y.tab.c"
    break;

  case 732: /* expr: expr '%' expr  */
#line 2792 "gram.y"
                        {
	    (yyval.val) = fmod((yyvsp[-2].val), (yyvsp[0].val));
	}
#line 10350 "y.tab.c"
    break;

  case 733: /* expr: expr '^' expr  */
#line 2795 "gram.y"
                        {
	    (yyval.val) = pow((yyvsp[-2].val), (yyvsp[0].val));
	}
#line 10358 "y.tab.c"
    break;

  case 734: /* expr: ABS '(' expr ')'  */
#line 2798 "gram.y"
                           {
	    (yyval.val) = fabs((yyvsp[-1].val));
	}
#line 10366 "y.tab.c"
    break;

  case 735: /* expr: ACOS '(' expr ')'  */
#line 2801 "gram.y"
                            {
	    (yyval.val) = acos((yyvsp[-1].val));
	}
#line 10374 "y.tab.c"
    break;

  case 736: /* expr: ASIN '(' expr ')'  */
#line 2804 "gram.y"
                            {
	    (yyval.val) = asin((yyvsp[-1].val));
	}
#line 10382 "y.tab.c"
    break;

  case 737: /* expr: ATAN '(' expr ')'  */
#line 2807 "gram.y"
                            {
	    (yyval.val) = atan((yyvsp[-1].val));
	}
#line 10390 "y.tab.c"
    break;

  case 738: /* expr: ATAN2 '(' expr ',' expr ')'  */
#line 2810 "gram.y"
                                      {
	    (yyval.val) = atan2((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 10398 "y.tab.c"
    break;

  case 739: /* expr: CEIL '(' expr ')'  */
#line 2813 "gram.y"
                            {
	    (yyval.val) = ceil((yyvsp[-1].val));
	}
#line 10406 "y.tab.c"
    break;

  case 740: /* expr: COS '(' expr ')'  */
#line 2816 "gram.y"
                           {
	    (yyval.val) = cos((yyvsp[-1].val));
	}
#line 10414 "y.tab.c"
    break;

  case 741: /* expr: DEG  */
#line 2819 "gram.y"
              {
	    (yyval.val) = 180.0 / M_PI;
	}
#line 10422 "y.tab.c"
    break;

  case 742: /* expr: DX  */
#line 2822 "gram.y"
             {
	    (yyval.val) = *xx;
	}
#line 10430 "y.tab.c"
    break;

  case 743: /* expr: DY  */
#line 2825 "gram.y"
             {
	    (yyval.val) = *yy;
	}
#line 10438 "y.tab.c"
    break;

  case 744: /* expr: ERF '(' expr ')'  */
#line 2828 "gram.y"
                           {
	    (yyval.val) = erf((yyvsp[-1].val));
	}
#line 10446 "y.tab.c"
    break;

  case 745: /* expr: ERFC '(' expr ')'  */
#line 2831 "gram.y"
                            {
	    (yyval.val) = erfc((yyvsp[-1].val));
	}
#line 10454 "y.tab.c"
    break;

  case 746: /* expr: EXP '(' expr ')'  */
#line 2834 "gram.y"
                           {
	    (yyval.val) = exp((yyvsp[-1].val));
	}
#line 10462 "y.tab.c"
    break;

  case 747: /* expr: FLOOR '(' expr ')'  */
#line 2837 "gram.y"
                             {
	    (yyval.val) = floor((yyvsp[-1].val));
	}
#line 10470 "y.tab.c"
    break;

  case 748: /* expr: HYPOT '(' expr ',' expr ')'  */
#line 2840 "gram.y"
                                      {
	    (yyval.val) = hypot((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 10478 "y.tab.c"
    break;

  case 749: /* expr: GRAPHNO '.' VX1  */
#line 2843 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.xv1;
	}
#line 10486 "y.tab.c"
    break;

  case 750: /* expr: GRAPHNO '.' VX2  */
#line 2846 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.xv2;
	}
#line 10494 "y.tab.c"
    break;

  case 751: /* expr: GRAPHNO '.' VY1  */
#line 2849 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.yv1;
	}
#line 10502 "y.tab.c"
    break;

  case 752: /* expr: GRAPHNO '.' VY2  */
#line 2852 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].v.yv2;
	}
#line 10510 "y.tab.c"
    break;

  case 753: /* expr: GRAPHNO '.' WX1  */
#line 2855 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.xg1;
	}
#line 10518 "y.tab.c"
    break;

  case 754: /* expr: GRAPHNO '.' WX2  */
#line 2858 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.xg2;
	}
#line 10526 "y.tab.c"
    break;

  case 755: /* expr: GRAPHNO '.' WY1  */
#line 2861 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.yg1;
	}
#line 10534 "y.tab.c"
    break;

  case 756: /* expr: GRAPHNO '.' WY2  */
#line 2864 "gram.y"
                          {
	    (yyval.val) = g[(yyvsp[-2].pset)].w.yg2;
	}
#line 10542 "y.tab.c"
    break;

  case 757: /* expr: VX1  */
#line 2867 "gram.y"
              {
	    (yyval.val) = g[curg].v.xv1;
	}
#line 10550 "y.tab.c"
    break;

  case 758: /* expr: VX2  */
#line 2870 "gram.y"
              {
	    (yyval.val) = g[curg].v.xv2;
	}
#line 10558 "y.tab.c"
    break;

  case 759: /* expr: VY1  */
#line 2873 "gram.y"
              {
	    (yyval.val) = g[curg].v.yv1;
	}
#line 10566 "y.tab.c"
    break;

  case 760: /* expr: VY2  */
#line 2876 "gram.y"
              {
	    (yyval.val) = g[curg].v.yv2;
	}
#line 10574 "y.tab.c"
    break;

  case 761: /* expr: WX1  */
#line 2879 "gram.y"
              {
	    (yyval.val) = g[curg].w.xg1;
	}
#line 10582 "y.tab.c"
    break;

  case 762: /* expr: WX2  */
#line 2882 "gram.y"
              {
	    (yyval.val) = g[curg].w.xg2;
	}
#line 10590 "y.tab.c"
    break;

  case 763: /* expr: WY1  */
#line 2885 "gram.y"
              {
	    (yyval.val) = g[curg].w.yg1;
	}
#line 10598 "y.tab.c"
    break;

  case 764: /* expr: WY2  */
#line 2888 "gram.y"
              {
	    (yyval.val) = g[curg].w.yg2;
	}
#line 10606 "y.tab.c"
    break;

  case 765: /* expr: INDEX  */
#line 2891 "gram.y"
                {
	    (yyval.val) = setindex;
	}
#line 10614 "y.tab.c"
    break;

  case 766: /* expr: SETNO  */
#line 2894 "gram.y"
                {
	    (yyval.val) = setsetno;
	}
#line 10622 "y.tab.c"
    break;

  case 767: /* expr: INT '(' expr ')'  */
#line 2897 "gram.y"
                           {
	    (yyval.val) = (long) (yyvsp[-1].val);
	}
#line 10630 "y.tab.c"
    break;

  case 768: /* expr: IRAND '(' NUMBER ')'  */
#line 2900 "gram.y"
                               {
	    (yyval.val) = lrand48() % (long) ((yyvsp[-1].val));
	}
#line 10638 "y.tab.c"
    break;

  case 769: /* expr: LGAMMA '(' expr ')'  */
#line 2903 "gram.y"
                              {
	    (yyval.val) = lgamma((yyvsp[-1].val));
	}
#line 10646 "y.tab.c"
    break;

  case 770: /* expr: LN '(' expr ')'  */
#line 2906 "gram.y"
                          {
	    (yyval.val) = log((yyvsp[-1].val));
	}
#line 10654 "y.tab.c"
    break;

  case 771: /* expr: LOG '(' expr ')'  */
#line 2909 "gram.y"
                           {
	    (yyval.val) = log10((yyvsp[-1].val));
	}
#line 10662 "y.tab.c"
    break;

  case 772: /* expr: LOGISTIC '(' expr ',' expr ',' expr ')'  */
#line 2913 "gram.y"
        {
	    (yyval.val) = 1.0 / (1.0 + exp(-((yyvsp[-5].val) - (yyvsp[-3].val))/ (yyvsp[-1].val)));
	}
#line 10670 "y.tab.c"
    break;

  case 773: /* expr: MAXP '(' expr ',' expr ')'  */
#line 2916 "gram.y"
                                     {
	    (yyval.val) = (yyvsp[-3].val) >= (yyvsp[-1].val) ? (yyvsp[-3].val) : (yyvsp[-1].val);
	}
#line 10678 "y.tab.c"
    break;

  case 774: /* expr: MINP '(' expr ',' expr ')'  */
#line 2919 "gram.y"
                                     {
	    (yyval.val) = (yyvsp[-3].val) <= (yyvsp[-1].val) ? (yyvsp[-3].val) : (yyvsp[-1].val);
	}
#line 10686 "y.tab.c"
    break;

  case 775: /* expr: MOD '(' expr ',' expr ')'  */
#line 2922 "gram.y"
                                    {
	    (yyval.val) = fmod((yyvsp[-3].val), (yyvsp[-1].val));
	}
#line 10694 "y.tab.c"
    break;

  case 776: /* expr: NORM '(' expr ')'  */
#line 2925 "gram.y"
                            {
	    (yyval.val) = fx((yyvsp[-1].val));
	}
#line 10702 "y.tab.c"
    break;

  case 777: /* expr: NORMP '(' expr ')'  */
#line 2928 "gram.y"
                             {
	    double tmp;
	    (yyval.val) = normp((yyvsp[-1].val), &tmp);
	}
#line 10711 "y.tab.c"
    break;

  case 778: /* expr: PI  */
#line 2932 "gram.y"
             {
	    (yyval.val) = M_PI;
	}
#line 10719 "y.tab.c"
    break;

  case 779: /* expr: RAD  */
#line 2935 "gram.y"
              {
	    (yyval.val) = M_PI / 180.0;
	}
#line 10727 "y.tab.c"
    break;

  case 780: /* expr: RAND  */
#line 2938 "gram.y"
               {
	    (yyval.val) = (double) drand48();
	}
#line 10735 "y.tab.c"
    break;

  case 781: /* expr: SIN '(' expr ')'  */
#line 2941 "gram.y"
                           {
	    (yyval.val) = sin((yyvsp[-1].val));
	}
#line 10743 "y.tab.c"
    break;

  case 782: /* expr: SQR '(' expr ')'  */
#line 2944 "gram.y"
                           {
	    (yyval.val) = pow((yyvsp[-1].val), 2.0);
	}
#line 10751 "y.tab.c"
    break;

  case 783: /* expr: SQRT '(' expr ')'  */
#line 2947 "gram.y"
                            {
	    (yyval.val) = sqrt((yyvsp[-1].val));
	}
#line 10759 "y.tab.c"
    break;

  case 784: /* expr: TAN '(' expr ')'  */
#line 2950 "gram.y"
                           {
	    (yyval.val) = tan((yyvsp[-1].val));
	}
#line 10767 "y.tab.c"
    break;

  case 785: /* expr: IF '(' expr ')' expr  */
#line 2953 "gram.y"
                               {
	    if ((int) (yyvsp[-2].val))
		(yyval.val) = (yyvsp[0].val);
	}
#line 10776 "y.tab.c"
    break;

  case 786: /* expr: IF '(' expr ')' expr ELSE expr  */
#line 2957 "gram.y"
                                         {
	    if ((int) (yyvsp[-4].val)) {
		(yyval.val) = (yyvsp[-2].val);
	    } else {
		(yyval.val) = (yyvsp[0].val);
	    }
	}
#line 10788 "y.tab.c"
    break;

  case 787: /* expr: expr GT expr  */
#line 2964 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) > (yyvsp[0].val);
	}
#line 10796 "y.tab.c"
    break;

  case 788: /* expr: expr LT expr  */
#line 2967 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) < (yyvsp[0].val);
	}
#line 10804 "y.tab.c"
    break;

  case 789: /* expr: expr LE expr  */
#line 2970 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) <= (yyvsp[0].val);
	}
#line 10812 "y.tab.c"
    break;

  case 790: /* expr: expr GE expr  */
#line 2973 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) >= (yyvsp[0].val);
	}
#line 10820 "y.tab.c"
    break;

  case 791: /* expr: expr EQ expr  */
#line 2976 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) == (yyvsp[0].val);
	}
#line 10828 "y.tab.c"
    break;

  case 792: /* expr: expr NE expr  */
#line 2979 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) != (yyvsp[0].val);
	}
#line 10836 "y.tab.c"
    break;

  case 793: /* expr: expr AND expr  */
#line 2982 "gram.y"
                        {
	    (yyval.val) = (yyvsp[-2].val) && (yyvsp[0].val);
	}
#line 10844 "y.tab.c"
    break;

  case 794: /* expr: expr OR expr  */
#line 2985 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-2].val) || (yyvsp[0].val);
	}
#line 10852 "y.tab.c"
    break;

  case 795: /* expr: NOT expr  */
#line 2988 "gram.y"
                   {
	    (yyval.val) = !((yyvsp[0].val));
	}
#line 10860 "y.tab.c"
    break;

  case 796: /* expr: '(' expr ')'  */
#line 2991 "gram.y"
                       {
	    (yyval.val) = (yyvsp[-1].val);
	}
#line 10868 "y.tab.c"
    break;

  case 797: /* expr: '-' expr  */
#line 2994 "gram.y"
                                {
	    (yyval.val) = -(yyvsp[0].val);
	}
#line 10876 "y.tab.c"
    break;


#line 10880 "y.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;


#if !defined yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturn;
#endif


/*-------------------------------------------------------.
| yyreturn -- parsing is finished, clean up and return.  |
`-------------------------------------------------------*/
yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 2998 "gram.y"


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
