/* A Bison parser, made by GNU Bison 2.7.  */

/* Bison implementation for Yacc-like parsers in C
   
      Copyright (C) 1984, 1989-1990, 2000-2012 Free Software Foundation, Inc.
   
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

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.7"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* Copy the first part of user declarations.  */
/* Line 371 of yacc.c  */
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


/* Line 371 of yacc.c  */
#line 170 "y.tab.c"

# ifndef YY_NULL
#  if defined __cplusplus && 201103L <= __cplusplus
#   define YY_NULL nullptr
#  else
#   define YY_NULL 0
#  endif
# endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* In a future release of Bison, this section will be replaced
   by #include "y.tab.h".  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     VAR = 258,
     X = 259,
     Y = 260,
     CHRSTR = 261,
     FITPARM = 262,
     NUMBER = 263,
     ABS = 264,
     ACOS = 265,
     ASIN = 266,
     ATAN = 267,
     ATAN2 = 268,
     CEIL = 269,
     COS = 270,
     DEG = 271,
     DX = 272,
     DY = 273,
     ERF = 274,
     ERFC = 275,
     EXP = 276,
     FLOOR = 277,
     HYPOT = 278,
     INDEX = 279,
     INT = 280,
     IRAND = 281,
     LGAMMA = 282,
     LN = 283,
     LOG = 284,
     LOGISTIC = 285,
     MAXP = 286,
     MINP = 287,
     MINMAX = 288,
     MOD = 289,
     NORM = 290,
     NORMP = 291,
     PI = 292,
     RAD = 293,
     RAND = 294,
     SETNO = 295,
     SIN = 296,
     SQR = 297,
     SQRT = 298,
     TAN = 299,
     INUM = 300,
     ABORT = 301,
     ABOVE = 302,
     ABSOLUTE = 303,
     ACTIVATE = 304,
     ACTIVE = 305,
     ADCIRC = 306,
     ADCIRC3DFLOW = 307,
     ALL = 308,
     ALT = 309,
     ALTERNATE = 310,
     ALTXAXIS = 311,
     ALTYAXIS = 312,
     AMP = 313,
     ANGLE = 314,
     ANNOTATE = 315,
     APPEND = 316,
     AREA = 317,
     ARROW = 318,
     ASCEND = 319,
     AT = 320,
     ATTACH = 321,
     AUTO = 322,
     AUTOSCALE = 323,
     AUTOTICKS = 324,
     AVERAGE = 325,
     AVG = 326,
     AXES = 327,
     AXIS = 328,
     BACKBUFFER = 329,
     BACKGROUND = 330,
     BAR = 331,
     BATCH = 332,
     BATH = 333,
     BATHYMETRY = 334,
     COURANT = 335,
     BELOW = 336,
     BIN = 337,
     BINARY = 338,
     BOTH = 339,
     BOTTOM = 340,
     BOUNDARY = 341,
     BOX = 342,
     CELLS = 343,
     CENTER = 344,
     CH3D = 345,
     CHAR = 346,
     CHDIR = 347,
     CIRCLE = 348,
     CLEAR = 349,
     CLICK = 350,
     CLOCK = 351,
     CLOSE = 352,
     CM = 353,
     CMAP = 354,
     COLOR = 355,
     COLORMAP = 356,
     COMMENT = 357,
     CONC = 358,
     CONCENTRATION = 359,
     CONCENTRATIONS = 360,
     COPY = 361,
     CROSS = 362,
     CYCLE = 363,
     DAYMONTH = 364,
     DAYOFWEEKL = 365,
     DAYOFWEEKS = 366,
     DAYOFYEAR = 367,
     DAYS = 368,
     DDMMYY = 369,
     DDMONTHSYYHHMMSS = 370,
     DECIMAL = 371,
     DEF = 372,
     DEFAULT = 373,
     DEGREESLAT = 374,
     DEGREESLON = 375,
     DEGREESMMLAT = 376,
     DEGREESMMLON = 377,
     DEGREESMMSSLAT = 378,
     DEGREESMMSSLON = 379,
     DELAYP = 380,
     DELETE = 381,
     DEPTH = 382,
     DEPTHS = 383,
     DESCEND = 384,
     DEVICE = 385,
     DEVXY = 386,
     DFT = 387,
     DT = 388,
     DIAMOND = 389,
     DIFFERENCE = 390,
     DISK = 391,
     DISPLAY = 392,
     DOT = 393,
     DOUBLEBUFFER = 394,
     DOWN = 395,
     DRAW2 = 396,
     DROGUE = 397,
     DROGUES = 398,
     DRY = 399,
     DXDX = 400,
     DXP = 401,
     DYDY = 402,
     DYP = 403,
     ECHO = 404,
     EDIT = 405,
     ELA = 406,
     ELCIRC = 407,
     ELEMENT = 408,
     ELEMENTS = 409,
     ELEV = 410,
     ELEVATION = 411,
     ELEVATIONS = 412,
     ELEVMARKER = 413,
     ELLIPSE = 414,
     ELLIPSES = 415,
     ELLIPSEZ = 416,
     ELSE = 417,
     END = 418,
     ERRORBAR = 419,
     EXIT = 420,
     EXPAND = 421,
     EXPONENTIAL = 422,
     FACTOR = 423,
     FALSEP = 424,
     FAST = 425,
     FEET = 426,
     FFT = 427,
     FILEP = 428,
     FILL = 429,
     FIND = 430,
     FIXEDPOINT = 431,
     FLOW = 432,
     FLUSH = 433,
     FLUX = 434,
     FOCUS = 435,
     FOLLOWS = 436,
     FONTP = 437,
     FOREGROUND = 438,
     FORMAT = 439,
     FORT14 = 440,
     FORT63 = 441,
     FORT64 = 442,
     FORWARD = 443,
     FRAMEP = 444,
     FREQ = 445,
     FRONTBUFFER = 446,
     GENERAL = 447,
     GETP = 448,
     GOTO = 449,
     GRAPH = 450,
     GRAPHNO = 451,
     GRAPHS = 452,
     GRAPHTYPE = 453,
     GRID = 454,
     HARDCOPY = 455,
     HBAR = 456,
     HELP = 457,
     HGAP = 458,
     HIDDEN = 459,
     HISTBOX = 460,
     HISTO = 461,
     HISTORY = 462,
     HMS = 463,
     HORIZONTAL = 464,
     HOURS = 465,
     HPGLL = 466,
     HPGLP = 467,
     IF = 468,
     IGNORE = 469,
     IHL = 470,
     IMAGE = 471,
     IMAGES = 472,
     IN = 473,
     INCLUDE = 474,
     INFO = 475,
     INIT = 476,
     INITGRAPHICS = 477,
     INOUT = 478,
     INTEGRATE = 479,
     INTERP = 480,
     INUNDATION = 481,
     INVDFT = 482,
     INVFFT = 483,
     ISOLINE = 484,
     ISOLINES = 485,
     JUST = 486,
     KILL = 487,
     KM = 488,
     LABEL = 489,
     LAYOUT = 490,
     LEAVE = 491,
     LEAVEGRAPHICS = 492,
     LEFT = 493,
     LEGEND = 494,
     LENGTH = 495,
     LEVEL = 496,
     LEVELS = 497,
     LIMITS = 498,
     LINE = 499,
     LINES = 500,
     LINESTYLE = 501,
     LINETO = 502,
     LINEW = 503,
     LINEWIDTH = 504,
     LINK = 505,
     LOAD = 506,
     LOC = 507,
     LOCATE = 508,
     LOCATOR = 509,
     LOCTYPE = 510,
     LOGX = 511,
     LOGXY = 512,
     LOGY = 513,
     M = 514,
     MAG = 515,
     MAGNITUDE = 516,
     MAJOR = 517,
     MAPSCALE = 518,
     MARKER = 519,
     MARKERS = 520,
     MAXLEVELS = 521,
     METHOD = 522,
     MIFL = 523,
     MIFP = 524,
     MILES = 525,
     MINOR = 526,
     MINUTES = 527,
     MISSINGP = 528,
     MM = 529,
     MMDD = 530,
     MMDDHMS = 531,
     MMDDYY = 532,
     MMDDYYHMS = 533,
     MMSSLAT = 534,
     MMSSLON = 535,
     MMYY = 536,
     MONTHDAY = 537,
     MONTHL = 538,
     MONTHS = 539,
     MOVE = 540,
     MOVE2 = 541,
     MOVETO = 542,
     NEGATE = 543,
     NO = 544,
     NODE = 545,
     NODES = 546,
     NONE = 547,
     NORMAL = 548,
     NORTH = 549,
     NXY = 550,
     OFF = 551,
     OFFSETX = 552,
     OFFSETY = 553,
     ON = 554,
     OP = 555,
     OPEN = 556,
     ORIENT = 557,
     OUT = 558,
     PAGE = 559,
     PARA = 560,
     PARALLEL = 561,
     PARAMETERS = 562,
     PARAMS = 563,
     PARMS = 564,
     PATTERN = 565,
     PER = 566,
     PERIMETER = 567,
     PERP = 568,
     PERPENDICULAR = 569,
     PHASE = 570,
     PIE = 571,
     PIPE = 572,
     PLACE = 573,
     PLAN = 574,
     PLUS = 575,
     POINT = 576,
     POLAR = 577,
     POLY = 578,
     POLYI = 579,
     POLYO = 580,
     POP = 581,
     POWER = 582,
     PREC = 583,
     PREFIX = 584,
     PREPEND = 585,
     PRINT = 586,
     PROFILE = 587,
     PROP = 588,
     PS = 589,
     PSCOLORL = 590,
     PSCOLORP = 591,
     PSMONOL = 592,
     PSMONOP = 593,
     PUSH = 594,
     PUTP = 595,
     QUIT = 596,
     READ = 597,
     READBIN = 598,
     REDRAW = 599,
     REGION = 600,
     REGIONS = 601,
     REGNUM = 602,
     REGRESS = 603,
     REMOVE = 604,
     RENDER = 605,
     REPORT = 606,
     RESET = 607,
     REVERSE = 608,
     REWIND = 609,
     RIGHT = 610,
     RISER = 611,
     ROT = 612,
     RUN = 613,
     SALINITY = 614,
     SAMPLE = 615,
     SAVE = 616,
     SCALAR = 617,
     SCALE = 618,
     SCIENTIFIC = 619,
     SECONDS = 620,
     SET = 621,
     SETS = 622,
     SHOW = 623,
     SHRINK = 624,
     SIGMA = 625,
     SIGN = 626,
     SIZE = 627,
     SKIP = 628,
     SLAB = 629,
     SLEEP = 630,
     SLICE = 631,
     SOURCE = 632,
     SPEC = 633,
     SPECIFIED = 634,
     SPECTRUM = 635,
     SPLITS = 636,
     SQUARE = 637,
     STACK = 638,
     STACKEDBAR = 639,
     STACKEDHBAR = 640,
     STACKEDLINE = 641,
     STAGGER = 642,
     STAR = 643,
     START = 644,
     STARTSTEP = 645,
     STARTTYPE = 646,
     STATION = 647,
     STATUS = 648,
     STEP = 649,
     STOP = 650,
     STREAMLINES = 651,
     STRING = 652,
     STRINGS = 653,
     SUBTITLE = 654,
     SURFACE = 655,
     SWAPBUFFER = 656,
     SYMBOL = 657,
     SYSTEM = 658,
     TEANL = 659,
     TEXT = 660,
     TICK = 661,
     TICKLABEL = 662,
     TICKMARKS = 663,
     TICKP = 664,
     TIDALCLOCK = 665,
     TIDESTATION = 666,
     TIME = 667,
     TIMEINFO = 668,
     TIMELINE = 669,
     TITLE = 670,
     TO = 671,
     TOP = 672,
     TOTAL = 673,
     TRACK = 674,
     TRANSECT = 675,
     TRIANGLE1 = 676,
     TRIANGLE2 = 677,
     TRIANGLE3 = 678,
     TRIANGLE4 = 679,
     TRUEP = 680,
     TYPE = 681,
     UNITS = 682,
     UP = 683,
     VALUE = 684,
     VECTOR = 685,
     VEL = 686,
     VELMARKER = 687,
     VELOCITY = 688,
     VERTICAL = 689,
     VGAP = 690,
     VIEW = 691,
     VSCALE = 692,
     VX1 = 693,
     VX2 = 694,
     VY1 = 695,
     VY2 = 696,
     WEEKS = 697,
     WET = 698,
     WETDRY = 699,
     WIDTH = 700,
     WIND = 701,
     WITH = 702,
     WORLD = 703,
     WRAP = 704,
     WRITE = 705,
     WSCALE = 706,
     WX1 = 707,
     WX2 = 708,
     WY1 = 709,
     WY2 = 710,
     X0 = 711,
     X1 = 712,
     X2 = 713,
     X3 = 714,
     X4 = 715,
     X5 = 716,
     XAXES = 717,
     XAXIS = 718,
     XCOR = 719,
     XMAX = 720,
     XMIN = 721,
     XY = 722,
     XYARC = 723,
     XYBOX = 724,
     XYDX = 725,
     XYDXDX = 726,
     XYDXDY = 727,
     XYDY = 728,
     XYDYDY = 729,
     XYFIXED = 730,
     XYHILO = 731,
     XYRT = 732,
     XYSEG = 733,
     XYSTRING = 734,
     XYUV = 735,
     XYX2Y2 = 736,
     XYXX = 737,
     XYYY = 738,
     XYZ = 739,
     XYZW = 740,
     Y0 = 741,
     Y1 = 742,
     Y2 = 743,
     Y3 = 744,
     Y4 = 745,
     Y5 = 746,
     YAXES = 747,
     YAXIS = 748,
     YEARS = 749,
     YES = 750,
     YMAX = 751,
     YMIN = 752,
     ZEROXAXIS = 753,
     ZEROYAXIS = 754,
     ZOOM = 755,
     ZOOMBOX = 756,
     OR = 757,
     AND = 758,
     NE = 759,
     EQ = 760,
     GE = 761,
     LE = 762,
     LT = 763,
     GT = 764,
     NOT = 765,
     UMINUS = 766
   };
#endif
/* Tokens.  */
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
#define NE 759
#define EQ 760
#define GE 761
#define LE 762
#define LT 763
#define GT 764
#define NOT 765
#define UMINUS 766



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{
/* Line 387 of yacc.c  */
#line 113 "gram.y"

    double val;
    int ival;
    double *ptr;
    int func;
    int pset;
    char *str;


/* Line 387 of yacc.c  */
#line 1245 "y.tab.c"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */

#endif /* !YY_YY_Y_TAB_H_INCLUDED  */

/* Copy the second part of user declarations.  */

/* Line 390 of yacc.c  */
#line 1273 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

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

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(E) ((void) (E))
#else
# define YYUSE(E) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(N) (N)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int yyi)
#else
static int
YYID (yyi)
    int yyi;
#endif
{
  return yyi;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

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
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
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
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
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
#   if ! defined malloc && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)				\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack_alloc, Stack, yysize);			\
	Stack = &yyptr->Stack_alloc;					\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, (Count) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYSIZE_T yyi;                         \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (YYID (0))
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
/* YYNRULES -- Number of states.  */
#define YYNSTATES  1798

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   766

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint16 yytranslate[] =
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
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint16 yyprhs[] =
{
       0,     0,     3,     4,     7,    10,    13,    16,    19,    22,
      25,    28,    31,    34,    37,    40,    43,    46,    49,    52,
      55,    58,    61,    64,    67,    70,    73,    76,    77,    80,
      83,    87,    90,    93,   102,   106,   110,   114,   118,   122,
     127,   132,   135,   138,   142,   145,   148,   151,   160,   164,
     168,   172,   176,   180,   185,   190,   193,   196,   199,   203,
     206,   209,   214,   218,   222,   226,   230,   234,   238,   242,
     247,   252,   257,   262,   267,   271,   273,   275,   277,   279,
     282,   285,   287,   292,   298,   306,   310,   312,   315,   319,
     322,   324,   326,   329,   331,   340,   342,   344,   347,   350,
     353,   356,   359,   363,   367,   370,   373,   377,   379,   382,
     391,   400,   409,   410,   414,   415,   419,   423,   426,   429,
     433,   436,   444,   448,   451,   455,   460,   466,   479,   493,
     506,   520,   532,   545,   557,   570,   575,   580,   581,   586,
     591,   597,   602,   607,   612,   617,   623,   629,   634,   640,
     645,   654,   659,   664,   672,   679,   685,   698,   702,   707,
     712,   718,   723,   729,   734,   739,   744,   749,   754,   759,
     764,   776,   789,   803,   818,   819,   825,   826,   832,   833,
     839,   840,   845,   852,   857,   862,   868,   873,   878,   883,
     889,   900,   907,   914,   921,   924,   928,   933,   938,   943,
     948,   953,   958,   963,   968,   972,   975,   976,   980,   981,
     986,   987,   991,   992,   996,   997,  1002,  1006,  1009,  1012,
    1016,  1021,  1026,  1032,  1037,  1042,  1047,  1052,  1059,  1065,
    1071,  1079,  1083,  1087,  1091,  1095,  1096,  1100,  1104,  1108,
    1114,  1120,  1123,  1127,  1131,  1136,  1143,  1151,  1156,  1157,
    1161,  1162,  1167,  1172,  1177,  1182,  1186,  1187,  1192,  1195,
    1197,  1202,  1206,  1210,  1211,  1215,  1219,  1229,  1233,  1238,
    1243,  1248,  1254,  1257,  1261,  1266,  1271,  1275,  1279,  1284,
    1291,  1296,  1302,  1308,  1313,  1318,  1324,  1329,  1334,  1339,
    1340,  1344,  1345,  1350,  1357,  1364,  1368,  1369,  1373,  1379,
    1383,  1389,  1393,  1398,  1402,  1408,  1415,  1419,  1423,  1428,
    1433,  1438,  1444,  1454,  1460,  1466,  1472,  1476,  1477,  1481,
    1487,  1491,  1495,  1500,  1503,  1507,  1511,  1515,  1519,  1523,
    1528,  1538,  1544,  1550,  1556,  1559,  1563,  1567,  1571,  1575,
    1581,  1585,  1588,  1591,  1595,  1599,  1603,  1607,  1613,  1617,
    1620,  1623,  1627,  1631,  1635,  1639,  1645,  1649,  1652,  1655,
    1659,  1664,  1669,  1673,  1679,  1682,  1685,  1689,  1694,  1698,
    1702,  1706,  1710,  1714,  1718,  1723,  1727,  1730,  1734,  1738,
    1742,  1746,  1750,  1754,  1758,  1762,  1767,  1771,  1777,  1780,
    1784,  1785,  1789,  1795,  1799,  1803,  1808,  1811,  1815,  1819,
    1823,  1828,  1838,  1844,  1850,  1856,  1860,  1864,  1868,  1872,
    1876,  1880,  1885,  1889,  1894,  1898,  1902,  1907,  1912,  1916,
    1921,  1926,  1929,  1933,  1936,  1939,  1942,  1945,  1948,  1951,
    1954,  1957,  1960,  1962,  1964,  1966,  1969,  1987,  1990,  1999,
    2003,  2007,  2011,  2015,  2024,  2028,  2032,  2036,  2040,  2043,
    2047,  2051,  2055,  2059,  2062,  2066,  2070,  2074,  2078,  2081,
    2085,  2089,  2093,  2097,  2101,  2106,  2109,  2113,  2118,  2123,
    2127,  2131,  2135,  2139,  2145,  2152,  2159,  2164,  2167,  2169,
    2172,  2176,  2179,  2181,  2183,  2185,  2187,  2189,  2191,  2194,
    2197,  2200,  2202,  2205,  2208,  2211,  2214,  2218,  2221,  2224,
    2227,  2230,  2233,  2235,  2237,  2240,  2243,  2246,  2249,  2252,
    2255,  2258,  2261,  2264,  2267,  2269,  2272,  2275,  2279,  2283,
    2286,  2289,  2293,  2297,  2301,  2305,  2309,  2313,  2317,  2321,
    2324,  2327,  2330,  2333,  2337,  2339,  2342,  2344,  2347,  2350,
    2353,  2356,  2359,  2362,  2365,  2368,  2371,  2374,  2377,  2380,
    2383,  2386,  2389,  2392,  2395,  2398,  2402,  2406,  2410,  2414,
    2417,  2420,  2424,  2427,  2430,  2433,  2437,  2439,  2442,  2445,
    2448,  2451,  2456,  2459,  2463,  2466,  2469,  2472,  2474,  2477,
    2480,  2483,  2485,  2487,  2489,  2491,  2493,  2495,  2497,  2499,
    2501,  2503,  2505,  2507,  2509,  2511,  2513,  2515,  2517,  2519,
    2521,  2523,  2525,  2527,  2529,  2531,  2533,  2535,  2537,  2539,
    2541,  2543,  2545,  2547,  2549,  2551,  2553,  2555,  2557,  2559,
    2561,  2563,  2565,  2567,  2569,  2571,  2573,  2575,  2577,  2579,
    2581,  2583,  2585,  2587,  2589,  2591,  2593,  2595,  2597,  2599,
    2601,  2603,  2605,  2607,  2609,  2611,  2613,  2615,  2617,  2619,
    2621,  2628,  2632,  2636,  2638,  2640,  2644,  2648,  2652,  2656,
    2660,  2664,  2668,  2672,  2676,  2680,  2684,  2688,  2692,  2696,
    2700,  2704,  2708,  2712,  2716,  2720,  2725,  2730,  2735,  2740,
    2745,  2752,  2757,  2762,  2764,  2766,  2768,  2773,  2778,  2783,
    2788,  2795,  2802,  2809,  2816,  2818,  2820,  2825,  2830,  2835,
    2840,  2845,  2854,  2861,  2868,  2875,  2880,  2885,  2887,  2889,
    2891,  2896,  2901,  2906,  2911,  2915,  2919,  2923,  2927,  2931,
    2935,  2939,  2943,  2946,  2950,  2953,  2955,  2957,  2962,  2966,
    2970,  2974,  2978,  2982,  2986,  2991,  2996,  3001,  3006,  3013,
    3018,  3023,  3025,  3027,  3029,  3034,  3039,  3044,  3049,  3056,
    3060,  3064,  3068,  3072,  3076,  3080,  3084,  3088,  3090,  3092,
    3094,  3096,  3098,  3100,  3102,  3104,  3106,  3108,  3113,  3118,
    3123,  3128,  3133,  3142,  3149,  3156,  3163,  3168,  3173,  3175,
    3177,  3179,  3184,  3189,  3194,  3199,  3205,  3213,  3217,  3221,
    3225,  3229,  3233,  3237,  3241,  3245,  3248,  3252
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int16 yyrhs[] =
{
     527,     0,    -1,    -1,   594,   519,    -1,   595,   519,    -1,
     597,   519,    -1,   596,   519,    -1,   529,   519,    -1,   528,
     519,    -1,   533,   519,    -1,   547,   519,    -1,   569,   519,
      -1,   570,   519,    -1,   551,   519,    -1,   545,   519,    -1,
     556,   519,    -1,   558,   519,    -1,   566,   519,    -1,   553,
     519,    -1,   563,   519,    -1,   565,   519,    -1,   564,   519,
      -1,   562,   519,    -1,   561,   519,    -1,   560,   519,    -1,
     530,   519,    -1,     1,   519,    -1,    -1,    94,    87,    -1,
     447,    87,    -1,   447,    87,     8,    -1,    87,   588,    -1,
      87,   196,    -1,    87,   597,   520,   597,   520,   597,   520,
     597,    -1,    87,   255,   589,    -1,    87,   246,     8,    -1,
      87,   249,     8,    -1,    87,   100,     8,    -1,    87,   174,
     591,    -1,    87,   174,   100,     8,    -1,    87,   174,   310,
       8,    -1,    87,   117,    -1,   447,   244,    -1,   447,   244,
       8,    -1,    94,   244,    -1,   244,   588,    -1,   244,   196,
      -1,   244,   597,   520,   597,   520,   597,   520,   597,    -1,
     244,   255,   589,    -1,   244,   249,     8,    -1,   244,   246,
       8,    -1,   244,   100,     8,    -1,   244,    63,     8,    -1,
     244,    63,   372,     8,    -1,   244,    63,   426,     8,    -1,
     244,   117,    -1,    94,   397,    -1,   447,   397,    -1,   447,
     397,     8,    -1,   397,   588,    -1,   397,   196,    -1,   397,
     597,   520,   597,    -1,   397,   255,   589,    -1,   397,   249,
       8,    -1,   397,   100,     8,    -1,   397,   357,     8,    -1,
     397,   182,     8,    -1,   397,   231,     8,    -1,   397,   402,
       8,    -1,   397,   402,   255,   592,    -1,   397,   402,   372,
       8,    -1,   397,   402,   174,     8,    -1,   397,   402,   100,
       8,    -1,   397,    91,   372,     8,    -1,   397,   117,     6,
      -1,   394,    -1,   353,    -1,   354,    -1,   188,    -1,   449,
     588,    -1,   194,     8,    -1,   358,    -1,   358,     8,   520,
       8,    -1,    77,   358,     8,   520,     8,    -1,    77,   358,
       8,   520,     8,   373,     8,    -1,    77,   329,     6,    -1,
     395,    -1,   403,     6,    -1,   358,    77,     6,    -1,    92,
       6,    -1,   352,    -1,   344,    -1,   375,     8,    -1,   341,
      -1,   500,   597,   520,   597,   520,   597,   520,   597,    -1,
     166,    -1,   369,    -1,   304,   238,    -1,   304,   355,    -1,
     304,   428,    -1,   304,   140,    -1,   304,   597,    -1,   304,
     223,     8,    -1,   250,   304,   588,    -1,    29,     6,    -1,
      97,    29,    -1,   219,   216,     6,    -1,   331,    -1,   149,
       6,    -1,    99,     8,   520,     8,   520,     8,   520,     8,
      -1,   101,     8,   520,     8,   520,     8,   520,     8,    -1,
     100,     8,   520,     8,   520,     8,   520,     8,    -1,    -1,
     362,   531,   553,    -1,    -1,   261,   532,   553,    -1,   447,
     404,     8,    -1,   404,   539,    -1,   342,   404,    -1,   447,
      51,     8,    -1,    51,   539,    -1,   342,    51,   155,     8,
     199,     6,     6,    -1,   447,   152,     8,    -1,   152,   539,
      -1,   152,   358,     8,    -1,   152,   199,     8,     6,    -1,
     342,   152,     8,   345,     6,    -1,   342,   152,     8,     6,
     389,     8,   395,     8,   373,     8,   241,     8,    -1,   342,
     152,     8,     6,   389,     8,   395,     8,   373,     8,   241,
       8,    61,    -1,   342,   152,     8,     6,   389,     8,   395,
       8,   373,     8,   127,     8,    -1,   342,   152,   400,     8,
       6,   389,     8,   395,     8,   373,     8,   127,     8,    -1,
     342,   152,   400,     8,     6,   389,     8,   395,     8,   373,
       8,    -1,   342,   152,   400,     8,     6,   389,     8,   395,
       8,   373,     8,    61,    -1,   342,   152,    85,     8,     6,
     389,     8,   395,     8,   373,     8,    -1,   342,   152,    85,
       8,     6,   389,     8,   395,     8,   373,     8,    61,    -1,
     447,   152,   420,     8,    -1,   366,   152,   420,     8,    -1,
      -1,   152,   420,   534,   553,    -1,   152,   420,   177,     6,
      -1,   152,   420,   434,   177,     6,    -1,   152,   420,   177,
       8,    -1,   152,   420,   359,     6,    -1,   152,   420,   155,
       6,    -1,   152,   420,   195,     8,    -1,   152,   420,   137,
     195,     8,    -1,   152,   420,   137,   244,   588,    -1,   152,
     420,   137,   588,    -1,   152,   420,   137,   260,   588,    -1,
     152,   420,   266,     8,    -1,   152,   420,   389,     8,   395,
       8,   373,     8,    -1,   152,   420,   360,     8,    -1,   152,
     420,   426,     8,    -1,   152,   420,     8,   520,   597,   520,
     597,    -1,   152,   420,   290,     8,   520,     8,    -1,   152,
     420,   345,   173,     6,    -1,   152,   420,   244,     8,   520,
     597,   520,   597,   520,   597,   520,   597,    -1,   342,   152,
     420,    -1,   342,   152,   420,    61,    -1,   342,   152,   177,
     420,    -1,   342,   152,   177,   420,    61,    -1,   342,   152,
     359,   420,    -1,   342,   152,   359,   420,    61,    -1,   447,
     152,   264,     8,    -1,   152,   264,   177,     6,    -1,   152,
     264,   261,     6,    -1,   152,   264,   430,     6,    -1,   152,
     264,   362,     6,    -1,   152,   264,   359,     6,    -1,   152,
     264,   155,     6,    -1,   342,   152,   264,   290,     8,   389,
       8,   394,     8,   373,     8,    -1,   342,   152,   264,   290,
       8,   389,     8,   394,     8,   373,     8,    61,    -1,   342,
     152,   264,   467,     8,   520,     8,   389,     8,   394,     8,
     373,     8,    -1,   342,   152,   264,   467,     8,   520,     8,
     389,     8,   394,     8,   373,     8,    61,    -1,    -1,   152,
     264,   362,   535,   553,    -1,    -1,   152,   264,   359,   536,
     553,    -1,    -1,   152,   264,   261,   537,   553,    -1,    -1,
     152,   264,   538,   568,    -1,   152,   264,   328,     8,   520,
       8,    -1,   152,   264,    66,     8,    -1,   152,   264,   255,
     589,    -1,   152,   264,   137,   264,   588,    -1,   152,   264,
     137,   588,    -1,   152,   264,   100,     8,    -1,   152,   264,
     249,     8,    -1,   152,   264,   174,   100,     8,    -1,   152,
     264,   448,   597,   520,   597,   520,   597,   520,   597,    -1,
     152,   264,   436,   597,   520,   597,    -1,   152,   264,   252,
     597,   520,   597,    -1,   152,   264,   467,   597,   520,   597,
      -1,   137,   587,    -1,   137,   155,   588,    -1,   137,   155,
     127,   588,    -1,   137,   155,   429,   588,    -1,   137,   155,
      31,   588,    -1,   137,   155,    58,   588,    -1,   137,   155,
     315,   588,    -1,   137,   155,   265,   588,    -1,   137,   177,
     260,   588,    -1,   137,   177,   446,   588,    -1,   137,   226,
     588,    -1,   100,     8,    -1,    -1,   155,   540,   553,    -1,
      -1,   155,    31,   541,   553,    -1,    -1,    58,   542,   553,
      -1,    -1,   315,   543,   553,    -1,    -1,   177,   260,   544,
     553,    -1,   177,   190,     8,    -1,   190,     8,    -1,   158,
       8,    -1,   360,   177,   588,    -1,   360,   177,   342,     6,
      -1,   360,   177,    32,     8,    -1,   360,   177,   342,   467,
       6,    -1,   360,   177,   426,   467,    -1,   360,   177,   426,
     290,    -1,   360,   177,   290,     8,    -1,   360,   177,   153,
       8,    -1,   360,   177,   467,   597,   520,   597,    -1,   126,
     360,   177,   290,     8,    -1,   126,   360,   177,   153,     8,
      -1,   126,   360,   177,   467,   597,   520,   597,    -1,   447,
     158,     8,    -1,   158,    50,   588,    -1,   158,   426,   588,
      -1,   158,   137,   588,    -1,    -1,   158,   546,   568,    -1,
     158,   255,   589,    -1,   158,   290,     8,    -1,   158,   252,
     597,   520,   597,    -1,   158,    33,   597,   520,   597,    -1,
     199,     8,    -1,   447,   199,     8,    -1,   199,   137,   588,
      -1,   199,    78,   137,   588,    -1,   199,    80,   137,   588,
     133,     8,    -1,   199,    80,   137,   429,   588,   133,     8,
      -1,   199,    86,   137,   588,    -1,    -1,   199,   548,   568,
      -1,    -1,   199,    86,   549,   568,    -1,   199,   291,   137,
     588,    -1,   199,   154,   137,   588,    -1,   199,   127,   137,
     588,    -1,   199,   174,   588,    -1,    -1,   199,    78,   550,
     553,    -1,   196,    68,    -1,    68,    -1,   342,   199,     8,
       6,    -1,   447,   143,     8,    -1,   143,   137,   588,    -1,
      -1,   143,   552,   568,    -1,   342,   143,     6,    -1,   342,
     143,     6,   389,     8,   395,     8,   373,     8,    -1,   230,
     239,   588,    -1,   230,   239,   235,   586,    -1,   230,   239,
     234,   588,    -1,   230,   239,   189,   588,    -1,   230,   239,
     189,   100,     8,    -1,   230,     8,    -1,   230,   426,     8,
      -1,   230,   366,   426,     8,    -1,   230,   174,   426,     8,
      -1,   230,   361,   299,    -1,   230,   361,   296,    -1,   230,
     361,   173,     6,    -1,   230,   361,   229,     8,   173,     6,
      -1,   230,   361,   229,     8,    -1,   230,   389,   597,   394,
     597,    -1,   230,    33,   597,   520,   597,    -1,   229,   597,
     520,   597,    -1,   230,   239,   255,   589,    -1,   230,   239,
     597,   520,   597,    -1,   229,     8,   100,     8,    -1,   229,
       8,   249,     8,    -1,   229,     8,   246,     8,    -1,    -1,
     230,   554,   568,    -1,    -1,   230,   239,   555,   568,    -1,
     230,   239,   372,     8,   520,     8,    -1,   230,   239,   203,
       8,   520,     8,    -1,   447,   205,     8,    -1,    -1,   205,
     557,   568,    -1,   205,   328,     8,   520,     8,    -1,   205,
      66,     8,    -1,   205,   409,     8,   520,     8,    -1,   205,
     255,   589,    -1,   205,   137,   264,   588,    -1,   205,   137,
     588,    -1,   205,   137,    51,     8,   590,    -1,   205,   137,
      51,     8,   100,     8,    -1,   205,   100,     8,    -1,   205,
     249,     8,    -1,   205,   174,   100,     8,    -1,   205,   342,
       8,     6,    -1,   205,   137,   207,   590,    -1,   205,   137,
     207,   100,     8,    -1,   205,   448,   597,   520,   597,   520,
     597,   520,   597,    -1,   205,   436,   597,   520,   597,    -1,
     205,   252,   597,   520,   597,    -1,   205,   467,   597,   520,
     597,    -1,   447,   501,     8,    -1,    -1,   501,   559,   568,
      -1,   501,   328,     8,   520,     8,    -1,   501,    66,     8,
      -1,   501,   255,   589,    -1,   501,   137,   264,   588,    -1,
     501,   588,    -1,   501,   137,   588,    -1,   501,   500,     8,
      -1,   501,   363,     8,    -1,   501,   100,     8,    -1,   501,
     249,     8,    -1,   501,   174,   100,     8,    -1,   501,   448,
     597,   520,   597,   520,   597,   520,   597,    -1,   501,   436,
     597,   520,   597,    -1,   501,   252,   597,   520,   597,    -1,
     501,   467,   597,   520,   597,    -1,   451,   588,    -1,   451,
     240,     8,    -1,   451,   363,     8,    -1,   451,   100,     8,
      -1,   451,   255,   589,    -1,   451,   252,   597,   520,   597,
      -1,   451,   427,   593,    -1,   451,   568,    -1,   437,   588,
      -1,   437,   240,     8,    -1,   437,   363,     8,    -1,   437,
     100,     8,    -1,   437,   255,   589,    -1,   437,   252,   597,
     520,   597,    -1,   437,   427,   593,    -1,   437,   568,    -1,
     263,   588,    -1,   263,   240,     8,    -1,   263,   100,     8,
      -1,   263,   363,     8,    -1,   263,   255,   589,    -1,   263,
     252,   597,   520,   597,    -1,   263,   427,   593,    -1,   263,
     568,    -1,   410,   588,    -1,   410,   100,     8,    -1,   410,
     174,   100,     8,    -1,   410,   418,   412,     8,    -1,   410,
     255,   589,    -1,   410,   252,   597,   520,   597,    -1,   410,
     568,    -1,   413,   588,    -1,   413,   389,     6,    -1,   413,
     597,   520,   597,    -1,   413,   255,   589,    -1,   413,   249,
       8,    -1,   413,   100,     8,    -1,   413,   357,     8,    -1,
     413,   182,     8,    -1,   413,   231,     8,    -1,   413,    91,
     372,     8,    -1,   413,   184,     6,    -1,   414,   588,    -1,
     414,   240,     8,    -1,   414,   445,     8,    -1,   414,   389,
       8,    -1,   414,   395,     8,    -1,   414,   394,     8,    -1,
     414,   328,     8,    -1,   414,   427,     8,    -1,   414,   100,
       8,    -1,   414,   174,   100,     8,    -1,   414,   255,   589,
      -1,   414,   252,   597,   520,   597,    -1,   414,   568,    -1,
     447,   376,     8,    -1,    -1,   376,   567,   568,    -1,   376,
     328,     8,   520,     8,    -1,   376,    66,     8,    -1,   376,
     255,   589,    -1,   376,   137,   264,   588,    -1,   376,   588,
      -1,   376,   137,   588,    -1,   376,   100,     8,    -1,   376,
     249,     8,    -1,   376,   174,   100,     8,    -1,   376,   448,
     597,   520,   597,   520,   597,   520,   597,    -1,   376,   436,
     597,   520,   597,    -1,   376,   252,   597,   520,   597,    -1,
     376,   467,   597,   520,   597,    -1,   333,   100,     8,    -1,
     333,   249,     8,    -1,   333,   246,     8,    -1,   333,   184,
     585,    -1,   333,   182,     8,    -1,   333,   328,     8,    -1,
     333,    91,   372,     8,    -1,   333,   402,     8,    -1,   333,
     402,   372,     8,    -1,   333,   174,   588,    -1,   333,   174,
     591,    -1,   333,   174,   100,     8,    -1,   333,   174,   310,
       8,    -1,   333,    63,     8,    -1,   333,    63,   426,     8,
      -1,   333,    63,   372,     8,    -1,   447,   196,    -1,   447,
     195,     8,    -1,   232,   196,    -1,   232,   197,    -1,   254,
     588,    -1,   180,   196,    -1,   180,   588,    -1,   180,   366,
      -1,   180,   181,    -1,   180,    95,    -1,   377,   580,    -1,
     339,    -1,   326,    -1,   108,    -1,   383,     8,    -1,   383,
     448,   597,   520,   597,   520,   597,   520,   597,   409,   597,
     520,   597,   520,   597,   520,   597,    -1,    94,   383,    -1,
     448,   597,   520,   597,   520,   597,   520,   597,    -1,   448,
     466,   597,    -1,   448,   465,   597,    -1,   448,   497,   597,
      -1,   448,   496,   597,    -1,   436,   597,   520,   597,   520,
     597,   520,   597,    -1,   436,   466,     8,    -1,   436,   465,
       8,    -1,   436,   497,     8,    -1,   436,   496,     8,    -1,
     415,     6,    -1,   415,   182,     8,    -1,   415,   372,     8,
      -1,   415,   100,     8,    -1,   415,   249,     8,    -1,   399,
       6,    -1,   399,   182,     8,    -1,   399,   372,     8,    -1,
     399,   100,     8,    -1,   399,   249,     8,    -1,   189,   588,
      -1,   189,   426,     8,    -1,   189,   246,     8,    -1,   189,
     249,     8,    -1,   189,   100,     8,    -1,   189,   174,   588,
      -1,   189,    75,   100,     8,    -1,   196,   588,    -1,   196,
     234,   588,    -1,   196,    68,   426,    67,    -1,   196,    68,
     426,   378,    -1,   196,    68,   590,    -1,   196,   204,   590,
      -1,   196,   426,   582,    -1,   196,   176,   588,    -1,   196,
     176,   184,   585,   585,    -1,   196,   176,   328,     8,   520,
       8,    -1,   196,   176,   467,   597,   520,   597,    -1,   196,
     176,   426,     8,    -1,   571,   574,    -1,   572,    -1,   197,
     571,    -1,   197,   571,   574,    -1,   197,   572,    -1,   463,
      -1,   493,    -1,    56,    -1,    57,    -1,   498,    -1,   499,
      -1,    72,   573,    -1,   462,   573,    -1,   492,   573,    -1,
     588,    -1,   100,     8,    -1,   249,     8,    -1,   246,     8,
      -1,   182,     8,    -1,    91,   372,     8,    -1,   199,   588,
      -1,   409,   575,    -1,   407,   576,    -1,   234,   578,    -1,
      76,   579,    -1,   588,    -1,   588,    -1,   262,   588,    -1,
     271,   588,    -1,   262,   597,    -1,   271,   597,    -1,   297,
     597,    -1,   298,   597,    -1,    54,   588,    -1,    32,   597,
      -1,    31,   597,    -1,   118,     8,    -1,   583,    -1,    29,
     588,    -1,   372,     8,    -1,   262,   372,     8,    -1,   271,
     372,     8,    -1,   100,     8,    -1,   249,     8,    -1,   262,
     100,     8,    -1,   271,   100,     8,    -1,   262,   249,     8,
      -1,   271,   249,     8,    -1,   262,   246,     8,    -1,   271,
     246,     8,    -1,   262,   199,   588,    -1,   271,   199,   588,
      -1,   300,   592,    -1,   426,    67,    -1,   426,   378,    -1,
     378,     8,    -1,     8,   520,   597,    -1,   577,    -1,   576,
     577,    -1,   588,    -1,   426,    67,    -1,   426,   378,    -1,
     328,     8,    -1,   184,   585,    -1,   184,     8,    -1,    61,
       6,    -1,   330,     6,    -1,   235,   209,    -1,   235,   434,
      -1,   235,   378,    -1,    59,     8,    -1,   231,   581,    -1,
     373,     8,    -1,   387,     8,    -1,   300,   592,    -1,   371,
     584,    -1,   389,   597,    -1,   395,   597,    -1,   389,   426,
     378,    -1,   389,   426,    67,    -1,   395,   426,   378,    -1,
     395,   426,    67,    -1,   435,     8,    -1,   203,     8,    -1,
      91,   372,     8,    -1,   182,     8,    -1,   100,     8,    -1,
     249,     8,    -1,     8,   520,     6,    -1,     6,    -1,   235,
     313,    -1,   235,   305,    -1,   318,    67,    -1,   318,   378,
      -1,   318,     8,   520,     8,    -1,   231,   581,    -1,    91,
     372,     8,    -1,   182,     8,    -1,   100,     8,    -1,   249,
       8,    -1,   588,    -1,   100,     8,    -1,   246,     8,    -1,
     249,     8,    -1,   136,    -1,   317,    -1,   355,    -1,   238,
      -1,    89,    -1,   467,    -1,   256,    -1,   258,    -1,   257,
      -1,   475,    -1,   218,    -1,   303,    -1,    84,    -1,   293,
      -1,    48,    -1,   288,    -1,   116,    -1,   167,    -1,   327,
      -1,   192,    -1,   114,    -1,   277,    -1,   281,    -1,   275,
      -1,   282,    -1,   109,    -1,   115,    -1,   284,    -1,   283,
      -1,   111,    -1,   110,    -1,   112,    -1,   208,    -1,   276,
      -1,   278,    -1,   120,    -1,   122,    -1,   124,    -1,   280,
      -1,   119,    -1,   121,    -1,   123,    -1,   279,    -1,   209,
      -1,   434,    -1,   296,    -1,   299,    -1,   291,    -1,    89,
      -1,   299,    -1,   296,    -1,   448,    -1,   436,    -1,   425,
      -1,   169,    -1,   310,    -1,   100,    -1,   292,    -1,    47,
      -1,    81,    -1,   238,    -1,   355,    -1,   417,    -1,    85,
      -1,    84,    -1,   274,    -1,    98,    -1,   259,    -1,   233,
      -1,     3,   521,   597,   522,   502,   597,    -1,     3,   502,
     596,    -1,     3,   502,   597,    -1,     3,    -1,   597,    -1,
     597,   511,   597,    -1,   596,   511,   596,    -1,   597,   511,
     596,    -1,   596,   511,   597,    -1,   597,   512,   597,    -1,
     596,   512,   596,    -1,   597,   512,   596,    -1,   596,   512,
     597,    -1,   597,   513,   597,    -1,   596,   513,   596,    -1,
     597,   513,   596,    -1,   596,   513,   597,    -1,   597,   514,
     597,    -1,   596,   514,   596,    -1,   597,   514,   596,    -1,
     596,   514,   597,    -1,   597,   516,   597,    -1,   597,   516,
     596,    -1,   596,   516,   597,    -1,   596,   516,   596,    -1,
       9,   523,   597,   524,    -1,     9,   523,   596,   524,    -1,
      10,   523,   596,   524,    -1,    11,   523,   596,   524,    -1,
      12,   523,   596,   524,    -1,    13,   523,   596,   520,   596,
     524,    -1,    14,   523,   596,   524,    -1,    15,   523,   596,
     524,    -1,    16,    -1,    17,    -1,    18,    -1,    19,   523,
     596,   524,    -1,    20,   523,   596,   524,    -1,    21,   523,
     596,   524,    -1,    22,   523,   596,   524,    -1,    23,   523,
     596,   520,   596,   524,    -1,    23,   523,   597,   520,   596,
     524,    -1,    23,   523,   596,   520,   597,   524,    -1,    23,
     523,   597,   520,   597,   524,    -1,    24,    -1,    40,    -1,
      25,   523,   596,   524,    -1,    26,   523,     8,   524,    -1,
      27,   523,   596,   524,    -1,    28,   523,   596,   524,    -1,
      29,   523,   596,   524,    -1,    30,   523,   596,   520,   597,
     520,   597,   524,    -1,    31,   523,   596,   520,   596,   524,
      -1,    32,   523,   596,   520,   596,   524,    -1,    34,   523,
     596,   520,   596,   524,    -1,    35,   523,   596,   524,    -1,
      36,   523,   596,   524,    -1,    37,    -1,    38,    -1,    39,
      -1,    41,   523,   596,   524,    -1,    42,   523,   596,   524,
      -1,    43,   523,   596,   524,    -1,    44,   523,   596,   524,
      -1,   596,   510,   596,    -1,   596,   509,   596,    -1,   596,
     508,   596,    -1,   596,   507,   596,    -1,   596,   506,   596,
      -1,   596,   505,   596,    -1,   596,   504,   596,    -1,   596,
     503,   596,    -1,   517,   596,    -1,   523,   596,   524,    -1,
     512,   596,    -1,     8,    -1,     7,    -1,     3,   521,   597,
     522,    -1,   597,   511,   597,    -1,   597,   512,   597,    -1,
     597,   513,   597,    -1,   597,   514,   597,    -1,   597,   515,
     597,    -1,   597,   516,   597,    -1,     9,   523,   597,   524,
      -1,    10,   523,   597,   524,    -1,    11,   523,   597,   524,
      -1,    12,   523,   597,   524,    -1,    13,   523,   597,   520,
     597,   524,    -1,    14,   523,   597,   524,    -1,    15,   523,
     597,   524,    -1,    16,    -1,    17,    -1,    18,    -1,    19,
     523,   597,   524,    -1,    20,   523,   597,   524,    -1,    21,
     523,   597,   524,    -1,    22,   523,   597,   524,    -1,    23,
     523,   597,   520,   597,   524,    -1,   196,   525,   438,    -1,
     196,   525,   439,    -1,   196,   525,   440,    -1,   196,   525,
     441,    -1,   196,   525,   452,    -1,   196,   525,   453,    -1,
     196,   525,   454,    -1,   196,   525,   455,    -1,   438,    -1,
     439,    -1,   440,    -1,   441,    -1,   452,    -1,   453,    -1,
     454,    -1,   455,    -1,    24,    -1,    40,    -1,    25,   523,
     597,   524,    -1,    26,   523,     8,   524,    -1,    27,   523,
     597,   524,    -1,    28,   523,   597,   524,    -1,    29,   523,
     597,   524,    -1,    30,   523,   597,   520,   597,   520,   597,
     524,    -1,    31,   523,   597,   520,   597,   524,    -1,    32,
     523,   597,   520,   597,   524,    -1,    34,   523,   597,   520,
     597,   524,    -1,    35,   523,   597,   524,    -1,    36,   523,
     597,   524,    -1,    37,    -1,    38,    -1,    39,    -1,    41,
     523,   597,   524,    -1,    42,   523,   597,   524,    -1,    43,
     523,   597,   524,    -1,    44,   523,   597,   524,    -1,   213,
     523,   597,   524,   597,    -1,   213,   523,   597,   524,   597,
     162,   597,    -1,   597,   510,   597,    -1,   597,   509,   597,
      -1,   597,   508,   597,    -1,   597,   507,   597,    -1,   597,
     506,   597,    -1,   597,   505,   597,    -1,   597,   504,   597,
      -1,   597,   503,   597,    -1,   517,   597,    -1,   523,   597,
     524,    -1,   512,   597,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
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

#if YYDEBUG || YYERROR_VERBOSE || 0
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "VAR", "X", "Y", "CHRSTR", "FITPARM",
  "NUMBER", "ABS", "ACOS", "ASIN", "ATAN", "ATAN2", "CEIL", "COS", "DEG",
  "DX", "DY", "ERF", "ERFC", "EXP", "FLOOR", "HYPOT", "INDEX", "INT",
  "IRAND", "LGAMMA", "LN", "LOG", "LOGISTIC", "MAXP", "MINP", "MINMAX",
  "MOD", "NORM", "NORMP", "PI", "RAD", "RAND", "SETNO", "SIN", "SQR",
  "SQRT", "TAN", "INUM", "ABORT", "ABOVE", "ABSOLUTE", "ACTIVATE",
  "ACTIVE", "ADCIRC", "ADCIRC3DFLOW", "ALL", "ALT", "ALTERNATE",
  "ALTXAXIS", "ALTYAXIS", "AMP", "ANGLE", "ANNOTATE", "APPEND", "AREA",
  "ARROW", "ASCEND", "AT", "ATTACH", "AUTO", "AUTOSCALE", "AUTOTICKS",
  "AVERAGE", "AVG", "AXES", "AXIS", "BACKBUFFER", "BACKGROUND", "BAR",
  "BATCH", "BATH", "BATHYMETRY", "COURANT", "BELOW", "BIN", "BINARY",
  "BOTH", "BOTTOM", "BOUNDARY", "BOX", "CELLS", "CENTER", "CH3D", "CHAR",
  "CHDIR", "CIRCLE", "CLEAR", "CLICK", "CLOCK", "CLOSE", "CM", "CMAP",
  "COLOR", "COLORMAP", "COMMENT", "CONC", "CONCENTRATION",
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
  "AND", "NE", "EQ", "GE", "LE", "LT", "GT", "'+'", "'-'", "'*'", "'/'",
  "'%'", "'^'", "NOT", "UMINUS", "'\\n'", "','", "'['", "']'", "'('",
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
  "filltype", "opchoice", "units", "asgn", "vasgn", "vexpr", "expr", YY_NULL
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
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
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint16 yyr1[] =
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

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
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

/* YYDEFACT[STATE-NAME] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint16 yydefact[] =
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
       0,   498,   514,   503,   721,   654,   720,   719,   718,   717,
     716,   715,   714,   656,   654,   660,   654,   664,   654,   668,
     654,   674,   654,   794,   793,   792,   791,   790,   789,   788,
     787,   657,   654,   661,   654,   665,   654,   669,   654,   732,
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

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -954
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
    1004,  1005,   850,   998,  1012,  1014,  1017,   828,    84,  5207,
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

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -798
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
    1258,   254,  1259,  1256,   382,  1260,  1265,   255,  1296,  1297,
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
       0,     0,     0,     0,   128,     0,   -27,   210,     0,     0,
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
     253,     0,     0,     0,     0,   254,     0,     0,   430,     0,
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
       0,   254,     0,     0,  1294,  1103,     0,   255,   564,   565,
     566,   567,   568,   569,   570,   571,   572,   573,   574,   575,
     576,   577,     0,     0,     0,     0,     0,     0,     0,  1122,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   253,
       0,     0,     0,     0,   254,     0,     0,     0,     0,     0,
     255,     0,     0,     0,     0,     0,     0,     0,     0,   107,
     108,   109,   110,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   115,   116,   117,   118,   210,     0,     0,
       0,     3,     4,   211,   212,   213,   214,   215,   216,   217,
     218,   219,   220,   221,   222,   223,   224,   225,   226,   227,
     228,   229,   230,   231,   232,   233,   234,     0,   235,   236,
     237,   238,   239,   240,   241,   242,   243,   244,   245,   564,
     565,   566,   567,   568,   569,   570,   571,   572,   573,   574,
     575,   576,   577,   253,     0,     0,     0,   345,   254,     0,
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
       0,   254,     0,     0,     0,     0,     0,   255,     0,     0,
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
       0,     0,   254,     0,     0,   210,     0,     0,   255,     3,
       4,   211,   212,   213,   214,   215,   216,   217,   218,   219,
     220,   221,   222,   223,   224,   225,   226,   227,   228,   229,
     230,   231,   232,   233,   234,     0,   235,   236,   237,   238,
     239,   240,   241,   242,   243,   244,   245,     0,     0,   253,
       0,     0,     0,     0,   254,     0,     0,     0,     0,   210,
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
     117,   118,     0,     0,   253,     0,     0,     0,     0,   254,
       0,   482,   483,     0,     0,   255,     0,     0,     0,  1457,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   107,   108,   109,   110,     0,     0,     0,     0,     0,
       0,     0,   484,   485,     0,   115,   116,   117,   118,     0,
       0,     0,     0,     0,   332,     0,     0,     0,   253,     0,
       0,  1459,     0,   254,     0,     0,     0,     0,     0,   255,
       0,    67,     0,   107,   108,   109,   110,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,   115,   116,   117,
     118,   564,   565,   566,   567,   568,   569,   570,   571,   703,
     704,   705,   706,   576,   707,   253,     0,     0,     0,     0,
     254,     0,  1116,     0,     0,     0,   255,   564,   565,   566,
     567,   568,   569,   570,   571,   703,   704,   705,   706,   576,
     707,     0,     0,     0,     0,     0,     0,     0,  1118,     0,
       0,     0,     0,     0,     0,     0,     0,   253,     0,     0,
       0,     0,   254,     0,     0,   210,     0,     0,   255,     3,
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
     253,  1128,     0,     0,     0,   254,     0,     0,     0,     0,
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
     574,   575,   576,   577,   253,     0,     0,     0,     0,   254,
       0,  1641,     0,     0,     0,   255,   564,   565,   566,   567,
     568,   569,   570,   571,   572,   573,   574,   575,   576,   577,
       0,     0,     0,     0,     0,     0,     0,  1643,     0,     0,
       0,     0,     0,     0,     0,     0,   127,     0,     0,     0,
       0,   128,     0,     0,     0,     0,     0,   129,   564,   565,
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

#define yypact_value_is_default(Yystate) \
  (!!((Yystate) == (-954)))

#define yytable_value_is_error(Yytable_value) \
  (!!((Yytable_value) == (-798)))

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
       8,   517,     8,   173,    66,     8,     8,   523,     8,     8,
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
      -1,    -1,    -1,    -1,   517,    -1,   519,     3,    -1,    -1,
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
     512,    -1,    -1,    -1,    -1,   517,    -1,    -1,   231,    -1,
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
      -1,   517,    -1,    -1,   327,   524,    -1,   523,   503,   504,
     505,   506,   507,   508,   509,   510,   511,   512,   513,   514,
     515,   516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,
      -1,    -1,    -1,    -1,   517,    -1,    -1,    -1,    -1,    -1,
     523,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   438,
     439,   440,   441,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   452,   453,   454,   455,     3,    -1,    -1,
      -1,     7,     8,     9,    10,    11,    12,    13,    14,    15,
      16,    17,    18,    19,    20,    21,    22,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    -1,    34,    35,
      36,    37,    38,    39,    40,    41,    42,    43,    44,   503,
     504,   505,   506,   507,   508,   509,   510,   511,   512,   513,
     514,   515,   516,   512,    -1,    -1,    -1,    63,   517,    -1,
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
      -1,   517,    -1,    -1,    -1,    -1,    -1,   523,    -1,    -1,
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
      -1,    -1,   517,    -1,    -1,     3,    -1,    -1,   523,     7,
       8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
      18,    19,    20,    21,    22,    23,    24,    25,    26,    27,
      28,    29,    30,    31,    32,    -1,    34,    35,    36,    37,
      38,    39,    40,    41,    42,    43,    44,    -1,    -1,   512,
      -1,    -1,    -1,    -1,   517,    -1,    -1,    -1,    -1,     3,
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
     454,   455,    -1,    -1,   512,    -1,    -1,    -1,    -1,   517,
      -1,   465,   466,    -1,    -1,   523,    -1,    -1,    -1,   426,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   438,   439,   440,   441,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   496,   497,    -1,   452,   453,   454,   455,    -1,
      -1,    -1,    -1,    -1,   196,    -1,    -1,    -1,   512,    -1,
      -1,   426,    -1,   517,    -1,    -1,    -1,    -1,    -1,   523,
      -1,   213,    -1,   438,   439,   440,   441,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   452,   453,   454,
     455,   503,   504,   505,   506,   507,   508,   509,   510,   511,
     512,   513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,
     517,    -1,   524,    -1,    -1,    -1,   523,   503,   504,   505,
     506,   507,   508,   509,   510,   511,   512,   513,   514,   515,
     516,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,
      -1,    -1,   517,    -1,    -1,     3,    -1,    -1,   523,     7,
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
     512,   524,    -1,    -1,    -1,   517,    -1,    -1,    -1,    -1,
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
     513,   514,   515,   516,   512,    -1,    -1,    -1,    -1,   517,
      -1,   524,    -1,    -1,    -1,   523,   503,   504,   505,   506,
     507,   508,   509,   510,   511,   512,   513,   514,   515,   516,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   524,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,   512,    -1,    -1,    -1,
      -1,   517,    -1,    -1,    -1,    -1,    -1,   523,   503,   504,
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
static const yytype_uint16 yystos[] =
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
     463,   492,   493,   498,   499,   500,   501,   512,   517,   523,
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
     246,   249,   255,   512,   517,   523,   588,   597,     6,    87,
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

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  However,
   YYFAIL appears to be in use.  Nevertheless, it is formally deprecated
   in Bison 2.4.2's NEWS entry, where a plan to phase it out is
   discussed.  */

#define YYFAIL		goto yyerrlab
#if defined YYFAIL
  /* This is here to suppress warnings from the GCC cpp's
     -Wunused-macros.  Normally we don't worry about that warning, but
     some users do, and we want to make it easy for users to remove
     YYFAIL uses, which will produce warnings from Bison 2.5.  */
#endif

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                  \
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
      YYERROR;							\
    }								\
while (YYID (0))

/* Error token number */
#define YYTERROR	1
#define YYERRCODE	256


/* This macro is provided for backward compatibility. */
#ifndef YY_LOCATION_PRINT
# define YY_LOCATION_PRINT(File, Loc) ((void) 0)
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */
#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  FILE *yyo = yyoutput;
  YYUSE (yyo);
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
        break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *yybottom, yytype_int16 *yytop)
#else
static void
yy_stack_print (yybottom, yytop)
    yytype_int16 *yybottom;
    yytype_int16 *yytop;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
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


#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into *YYMSG, which is of size *YYMSG_ALLOC, an error message
   about the unexpected token YYTOKEN for the state stack whose top is
   YYSSP.

   Return 0 if *YYMSG was successfully written.  Return 1 if *YYMSG is
   not large enough to hold the message.  In that case, also set
   *YYMSG_ALLOC to the required number of bytes.  Return 2 if the
   required number of bytes is too large to store.  */
static int
yysyntax_error (YYSIZE_T *yymsg_alloc, char **yymsg,
                yytype_int16 *yyssp, int yytoken)
{
  YYSIZE_T yysize0 = yytnamerr (YY_NULL, yytname[yytoken]);
  YYSIZE_T yysize = yysize0;
  enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
  /* Internationalized format string. */
  const char *yyformat = YY_NULL;
  /* Arguments of yyformat. */
  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
  /* Number of reported tokens (one for the "unexpected", one per
     "expected"). */
  int yycount = 0;

  /* There are many possibilities here to consider:
     - Assume YYFAIL is not used.  It's too flawed to consider.  See
       <http://lists.gnu.org/archive/html/bison-patches/2009-12/msg00024.html>
       for details.  YYERROR is fine as it does not invoke this
       function.
     - If this state is a consistent state with a default action, then
       the only way this function was invoked is if the default action
       is an error action.  In that case, don't check for expected
       tokens because there are none.
     - The only way there can be no lookahead present (in yychar) is if
       this state is a consistent state with a default action.  Thus,
       detecting the absence of a lookahead is sufficient to determine
       that there is no unexpected or expected token to report.  In that
       case, just report a simple "syntax error".
     - Don't assume there isn't a lookahead just because this state is a
       consistent state with a default action.  There might have been a
       previous inconsistent state, consistent state with a non-default
       action, or user semantic action that manipulated yychar.
     - Of course, the expected token list depends on states to have
       correct lookahead information, and it depends on the parser not
       to perform extra reductions after fetching a lookahead from the
       scanner and before detecting a syntax error.  Thus, state merging
       (from LALR or IELR) and default reductions corrupt the expected
       token list.  However, the list is correct for canonical LR with
       one exception: it will still contain any token that will not be
       accepted due to an error action in a later state.
  */
  if (yytoken != YYEMPTY)
    {
      int yyn = yypact[*yyssp];
      yyarg[yycount++] = yytname[yytoken];
      if (!yypact_value_is_default (yyn))
        {
          /* Start YYX at -YYN if negative to avoid negative indexes in
             YYCHECK.  In other words, skip the first -YYN actions for
             this state because they are default actions.  */
          int yyxbegin = yyn < 0 ? -yyn : 0;
          /* Stay within bounds of both yycheck and yytname.  */
          int yychecklim = YYLAST - yyn + 1;
          int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
          int yyx;

          for (yyx = yyxbegin; yyx < yyxend; ++yyx)
            if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR
                && !yytable_value_is_error (yytable[yyx + yyn]))
              {
                if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                  {
                    yycount = 1;
                    yysize = yysize0;
                    break;
                  }
                yyarg[yycount++] = yytname[yyx];
                {
                  YYSIZE_T yysize1 = yysize + yytnamerr (YY_NULL, yytname[yyx]);
                  if (! (yysize <= yysize1
                         && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
                    return 2;
                  yysize = yysize1;
                }
              }
        }
    }

  switch (yycount)
    {
# define YYCASE_(N, S)                      \
      case N:                               \
        yyformat = S;                       \
      break
      YYCASE_(0, YY_("syntax error"));
      YYCASE_(1, YY_("syntax error, unexpected %s"));
      YYCASE_(2, YY_("syntax error, unexpected %s, expecting %s"));
      YYCASE_(3, YY_("syntax error, unexpected %s, expecting %s or %s"));
      YYCASE_(4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
      YYCASE_(5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
# undef YYCASE_
    }

  {
    YYSIZE_T yysize1 = yysize + yystrlen (yyformat);
    if (! (yysize <= yysize1 && yysize1 <= YYSTACK_ALLOC_MAXIMUM))
      return 2;
    yysize = yysize1;
  }

  if (*yymsg_alloc < yysize)
    {
      *yymsg_alloc = 2 * yysize;
      if (! (yysize <= *yymsg_alloc
             && *yymsg_alloc <= YYSTACK_ALLOC_MAXIMUM))
        *yymsg_alloc = YYSTACK_ALLOC_MAXIMUM;
      return 1;
    }

  /* Avoid sprintf, as that infringes on the user's name space.
     Don't have undefined behavior even if the translation
     produced a string with the wrong number of "%s"s.  */
  {
    char *yyp = *yymsg;
    int yyi = 0;
    while ((*yyp = *yyformat) != '\0')
      if (*yyp == '%' && yyformat[1] == 's' && yyi < yycount)
        {
          yyp += yytnamerr (yyp, yyarg[yyi++]);
          yyformat += 2;
        }
      else
        {
          yyp++;
          yyformat++;
        }
  }
  return 0;
}
#endif /* YYERROR_VERBOSE */

/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
        break;
    }
}




/* The lookahead symbol.  */
int yychar;


#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval YY_INITIAL_VALUE(yyval_default);

/* Number of syntax errors so far.  */
int yynerrs;


/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
    int yystate;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;

    /* The stacks and their tools:
       `yyss': related to states.
       `yyvs': related to semantic values.

       Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack.  */
    yytype_int16 yyssa[YYINITDEPTH];
    yytype_int16 *yyss;
    yytype_int16 *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs;
    YYSTYPE *yyvsp;

    YYSIZE_T yystacksize;

  int yyn;
  int yyresult;
  /* Lookahead token as an internal (translated) token number.  */
  int yytoken = 0;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;

#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  yyssp = yyss = yyssa;
  yyvsp = yyvs = yyvsa;
  yystacksize = YYINITDEPTH;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY; /* Cause a token to be read.  */
  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;

	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),
		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss_alloc, yyss);
	YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

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

  /* YYCHAR is either YYEMPTY or YYEOF or a valid lookahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
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

  /* Discard the shifted token.  */
  yychar = YYEMPTY;

  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

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
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 3:
/* Line 1792 of yacc.c  */
#line 669 "gram.y"
    {}
    break;

  case 4:
/* Line 1792 of yacc.c  */
#line 670 "gram.y"
    {}
    break;

  case 5:
/* Line 1792 of yacc.c  */
#line 671 "gram.y"
    { result = (yyvsp[(1) - (2)].val); }
    break;

  case 6:
/* Line 1792 of yacc.c  */
#line 672 "gram.y"
    { result = *(yyvsp[(1) - (2)].ptr); }
    break;

  case 18:
/* Line 1792 of yacc.c  */
#line 684 "gram.y"
    { }
    break;

  case 26:
/* Line 1792 of yacc.c  */
#line 692 "gram.y"
    { return 1; }
    break;

  case 28:
/* Line 1792 of yacc.c  */
#line 696 "gram.y"
    { do_clear_boxes(); }
    break;

  case 29:
/* Line 1792 of yacc.c  */
#line 697 "gram.y"
    { curbox = next_box(); }
    break;

  case 30:
/* Line 1792 of yacc.c  */
#line 698 "gram.y"
    { curbox = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 31:
/* Line 1792 of yacc.c  */
#line 699 "gram.y"
    { boxes[curbox].active = (yyvsp[(2) - (2)].pset); }
    break;

  case 32:
/* Line 1792 of yacc.c  */
#line 700 "gram.y"
    { boxes[curbox].gno = (yyvsp[(2) - (2)].pset); }
    break;

  case 33:
/* Line 1792 of yacc.c  */
#line 702 "gram.y"
    {
	    if (curbox >= 0 && curbox < maxboxes) {
		boxes[curbox].x1 = (yyvsp[(2) - (8)].val);
		boxes[curbox].y1 = (yyvsp[(4) - (8)].val);
		boxes[curbox].x2 = (yyvsp[(6) - (8)].val);
		boxes[curbox].y2 = (yyvsp[(8) - (8)].val);
	    }
	}
    break;

  case 34:
/* Line 1792 of yacc.c  */
#line 710 "gram.y"
    { sysbox.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 35:
/* Line 1792 of yacc.c  */
#line 711 "gram.y"
    { sysbox.lines = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 36:
/* Line 1792 of yacc.c  */
#line 712 "gram.y"
    { sysbox.linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 37:
/* Line 1792 of yacc.c  */
#line 713 "gram.y"
    { sysbox.color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 38:
/* Line 1792 of yacc.c  */
#line 714 "gram.y"
    { sysbox.fill = (yyvsp[(3) - (3)].pset); }
    break;

  case 39:
/* Line 1792 of yacc.c  */
#line 715 "gram.y"
    { sysbox.fillcolor = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 40:
/* Line 1792 of yacc.c  */
#line 716 "gram.y"
    { sysbox.fillpattern = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 41:
/* Line 1792 of yacc.c  */
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
    break;

  case 42:
/* Line 1792 of yacc.c  */
#line 729 "gram.y"
    { curline = next_line(); }
    break;

  case 43:
/* Line 1792 of yacc.c  */
#line 730 "gram.y"
    { curline = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 44:
/* Line 1792 of yacc.c  */
#line 731 "gram.y"
    { do_clear_lines(); }
    break;

  case 45:
/* Line 1792 of yacc.c  */
#line 732 "gram.y"
    { lines[curline].active = (yyvsp[(2) - (2)].pset); }
    break;

  case 46:
/* Line 1792 of yacc.c  */
#line 733 "gram.y"
    { lines[curline].gno = (yyvsp[(2) - (2)].pset); }
    break;

  case 47:
/* Line 1792 of yacc.c  */
#line 735 "gram.y"
    {
	    lines[curline].x1 = (yyvsp[(2) - (8)].val);
	    lines[curline].y1 = (yyvsp[(4) - (8)].val);
	    lines[curline].x2 = (yyvsp[(6) - (8)].val);
	    lines[curline].y2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 48:
/* Line 1792 of yacc.c  */
#line 741 "gram.y"
    { sysline.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 49:
/* Line 1792 of yacc.c  */
#line 742 "gram.y"
    { sysline.linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 50:
/* Line 1792 of yacc.c  */
#line 743 "gram.y"
    { sysline.lines = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 51:
/* Line 1792 of yacc.c  */
#line 744 "gram.y"
    { sysline.color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 52:
/* Line 1792 of yacc.c  */
#line 745 "gram.y"
    { sysline.arrow = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 53:
/* Line 1792 of yacc.c  */
#line 746 "gram.y"
    { sysline.asize = (yyvsp[(4) - (4)].val); }
    break;

  case 54:
/* Line 1792 of yacc.c  */
#line 747 "gram.y"
    { sysline.atype = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 55:
/* Line 1792 of yacc.c  */
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
    break;

  case 56:
/* Line 1792 of yacc.c  */
#line 760 "gram.y"
    { do_clear_text(); }
    break;

  case 57:
/* Line 1792 of yacc.c  */
#line 761 "gram.y"
    { curstring = next_string(); }
    break;

  case 58:
/* Line 1792 of yacc.c  */
#line 762 "gram.y"
    { curstring = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 59:
/* Line 1792 of yacc.c  */
#line 763 "gram.y"
    { pstr[curstring].active = (yyvsp[(2) - (2)].pset); }
    break;

  case 60:
/* Line 1792 of yacc.c  */
#line 764 "gram.y"
    { pstr[curstring].gno = (yyvsp[(2) - (2)].pset); }
    break;

  case 61:
/* Line 1792 of yacc.c  */
#line 766 "gram.y"
    {
	    pstr[curstring].x = (yyvsp[(2) - (4)].val);
	    pstr[curstring].y = (yyvsp[(4) - (4)].val);
	}
    break;

  case 62:
/* Line 1792 of yacc.c  */
#line 770 "gram.y"
    { sysstr.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 63:
/* Line 1792 of yacc.c  */
#line 771 "gram.y"
    { sysstr.linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 64:
/* Line 1792 of yacc.c  */
#line 772 "gram.y"
    { sysstr.color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 65:
/* Line 1792 of yacc.c  */
#line 773 "gram.y"
    { sysstr.rot = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 66:
/* Line 1792 of yacc.c  */
#line 774 "gram.y"
    { sysstr.font = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 67:
/* Line 1792 of yacc.c  */
#line 775 "gram.y"
    { sysstr.just = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 68:
/* Line 1792 of yacc.c  */
#line 776 "gram.y"
    { sysstr.sym = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 69:
/* Line 1792 of yacc.c  */
#line 777 "gram.y"
    { sysstr.symloc = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 70:
/* Line 1792 of yacc.c  */
#line 778 "gram.y"
    { sysstr.symsize = (double) (yyvsp[(4) - (4)].val); }
    break;

  case 71:
/* Line 1792 of yacc.c  */
#line 779 "gram.y"
    { sysstr.symfill = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 72:
/* Line 1792 of yacc.c  */
#line 780 "gram.y"
    { sysstr.symcolor = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 73:
/* Line 1792 of yacc.c  */
#line 781 "gram.y"
    { sysstr.charsize = (double) (yyvsp[(4) - (4)].val); }
    break;

  case 74:
/* Line 1792 of yacc.c  */
#line 783 "gram.y"
    {
	    strcpy(pstr[curstring].s, (char *) (yyvsp[(3) - (3)].str));
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
    break;

  case 75:
/* Line 1792 of yacc.c  */
#line 802 "gram.y"
    { setistep(); }
    break;

  case 76:
/* Line 1792 of yacc.c  */
#line 803 "gram.y"
    { setreverse(); }
    break;

  case 77:
/* Line 1792 of yacc.c  */
#line 804 "gram.y"
    { setrewind(); }
    break;

  case 78:
/* Line 1792 of yacc.c  */
#line 805 "gram.y"
    { setforward(); }
    break;

  case 79:
/* Line 1792 of yacc.c  */
#line 806 "gram.y"
    { set_wrap((yyvsp[(2) - (2)].pset) == ON); }
    break;

  case 80:
/* Line 1792 of yacc.c  */
#line 807 "gram.y"
    { goto_step((int) (yyvsp[(2) - (2)].val) - 1); }
    break;

  case 81:
/* Line 1792 of yacc.c  */
#line 808 "gram.y"
    { setirun(); }
    break;

  case 82:
/* Line 1792 of yacc.c  */
#line 809 "gram.y"
    { runsteps((int) (yyvsp[(2) - (4)].val), (int) (yyvsp[(4) - (4)].val)); }
    break;

  case 83:
/* Line 1792 of yacc.c  */
#line 810 "gram.y"
    { batchrunsteps((int) (yyvsp[(3) - (5)].val), (int) (yyvsp[(5) - (5)].val), 1); }
    break;

  case 84:
/* Line 1792 of yacc.c  */
#line 811 "gram.y"
    { batchrunsteps((int) (yyvsp[(3) - (7)].val), (int) (yyvsp[(5) - (7)].val), (int) (yyvsp[(7) - (7)].val)); }
    break;

  case 85:
/* Line 1792 of yacc.c  */
#line 812 "gram.y"
    { strcpy(batchprefix, (char *) (yyvsp[(3) - (3)].str)); }
    break;

  case 86:
/* Line 1792 of yacc.c  */
#line 813 "gram.y"
    { setistop(); }
    break;

  case 87:
/* Line 1792 of yacc.c  */
#line 817 "gram.y"
    { system((yyvsp[(2) - (2)].str)); }
    break;

  case 88:
/* Line 1792 of yacc.c  */
#line 819 "gram.y"
    {
	    gotbatch = 1;
	    batchfile[0] = 0;
	    strcpy(batchfile, (yyvsp[(3) - (3)].str));
	}
    break;

  case 89:
/* Line 1792 of yacc.c  */
#line 824 "gram.y"
    {
		if (chdir((char *) (yyvsp[(2) - (2)].str)) < 0) {
			sprintf(buf, "chdir() to %s failed", (char *) (yyvsp[(2) - (2)].str));
			errwin(buf);
		}
	}
    break;

  case 90:
/* Line 1792 of yacc.c  */
#line 830 "gram.y"
    { setreset_world(); }
    break;

  case 91:
/* Line 1792 of yacc.c  */
#line 831 "gram.y"
    { setredraw_world(); }
    break;

  case 92:
/* Line 1792 of yacc.c  */
#line 832 "gram.y"
    { sleep((int) (yyvsp[(2) - (2)].val)); }
    break;

  case 93:
/* Line 1792 of yacc.c  */
#line 833 "gram.y"
    { exit(0); }
    break;

  case 94:
/* Line 1792 of yacc.c  */
#line 834 "gram.y"
    { my_blowup((yyvsp[(2) - (8)].val), (yyvsp[(4) - (8)].val), (yyvsp[(6) - (8)].val), (yyvsp[(8) - (8)].val)); }
    break;

  case 95:
/* Line 1792 of yacc.c  */
#line 835 "gram.y"
    { page(page_per, 4); }
    break;

  case 96:
/* Line 1792 of yacc.c  */
#line 836 "gram.y"
    { page(page_per, 5); }
    break;

  case 97:
/* Line 1792 of yacc.c  */
#line 837 "gram.y"
    { page(page_per, 0); }
    break;

  case 98:
/* Line 1792 of yacc.c  */
#line 838 "gram.y"
    { page(page_per, 1); }
    break;

  case 99:
/* Line 1792 of yacc.c  */
#line 839 "gram.y"
    { page(page_per, 2); }
    break;

  case 100:
/* Line 1792 of yacc.c  */
#line 840 "gram.y"
    { page(page_per, 3); }
    break;

  case 101:
/* Line 1792 of yacc.c  */
#line 841 "gram.y"
    { page_per = (yyvsp[(2) - (2)].val); }
    break;

  case 102:
/* Line 1792 of yacc.c  */
#line 842 "gram.y"
    { scrollinout_proc((int) (yyvsp[(3) - (3)].val)); }
    break;

  case 103:
/* Line 1792 of yacc.c  */
#line 843 "gram.y"
    { scrolling_islinked = (yyvsp[(3) - (3)].pset) == ON; }
    break;

  case 104:
/* Line 1792 of yacc.c  */
#line 845 "gram.y"
    {
	    if ((logfp = fopen((yyvsp[(2) - (2)].str), "w")) != NULL) {
		logfile = 1;
		printf("Opened logfile %s\n", (yyvsp[(2) - (2)].str));
	    } else {
		logfile = 0;
		printf("Failed to open logfile %s\n", (yyvsp[(2) - (2)].str));
	    }
	}
    break;

  case 105:
/* Line 1792 of yacc.c  */
#line 855 "gram.y"
    {
	    if (logfp != NULL) {
		printf("Closing logfile\n");
		logfile = 0;
		fclose(logfp);
		logfp = NULL;
	    }
	}
    break;

  case 106:
/* Line 1792 of yacc.c  */
#line 863 "gram.y"
    { }
    break;

  case 107:
/* Line 1792 of yacc.c  */
#line 864 "gram.y"
    { batchrunstep(0); }
    break;

  case 108:
/* Line 1792 of yacc.c  */
#line 865 "gram.y"
    {
	    if (inwin) { set_left_footer((yyvsp[(2) - (2)].str)); }
	    else { printf("%s\n", (yyvsp[(2) - (2)].str)); }
	}
    break;

  case 109:
/* Line 1792 of yacc.c  */
#line 870 "gram.y"
    { set_colormapdata((int) (yyvsp[(2) - (8)].val), (int) (yyvsp[(4) - (8)].val), (int) (yyvsp[(6) - (8)].val), (int) (yyvsp[(8) - (8)].val)); }
    break;

  case 110:
/* Line 1792 of yacc.c  */
#line 872 "gram.y"
    { set_colormapdata((int) (yyvsp[(2) - (8)].val), (int) (yyvsp[(4) - (8)].val), (int) (yyvsp[(6) - (8)].val), (int) (yyvsp[(8) - (8)].val)); }
    break;

  case 111:
/* Line 1792 of yacc.c  */
#line 874 "gram.y"
    { set_colormapdata((int) (yyvsp[(2) - (8)].val), (int) (yyvsp[(4) - (8)].val), (int) (yyvsp[(6) - (8)].val), (int) (yyvsp[(8) - (8)].val)); }
    break;

  case 112:
/* Line 1792 of yacc.c  */
#line 875 "gram.y"
    { setisol = &(g[curg].salip); }
    break;

  case 114:
/* Line 1792 of yacc.c  */
#line 876 "gram.y"
    { setisol = &(g[curg].velmagip); }
    break;

  case 116:
/* Line 1792 of yacc.c  */
#line 880 "gram.y"
    { 
		setflow = &g[curg].flowf[(int) (yyvsp[(3) - (3)].val)]; 
	}
    break;

  case 117:
/* Line 1792 of yacc.c  */
#line 883 "gram.y"
    { }
    break;

  case 118:
/* Line 1792 of yacc.c  */
#line 884 "gram.y"
    {}
    break;

  case 119:
/* Line 1792 of yacc.c  */
#line 885 "gram.y"
    { 
		setflow = &g[curg].flowt[(int) (yyvsp[(3) - (3)].val)]; 
                elcirc_flowno = (int) (yyvsp[(3) - (3)].val);
	}
    break;

  case 120:
/* Line 1792 of yacc.c  */
#line 889 "gram.y"
    {}
    break;

  case 121:
/* Line 1792 of yacc.c  */
#line 891 "gram.y"
    {
	      int fno = (int) (yyvsp[(4) - (7)].val);
              readbin_adcirc_elev(fno, (char *) (yyvsp[(6) - (7)].str), (char *) (yyvsp[(7) - (7)].str));
              set_clock(0, flowt[fno].start, 
				flowt[fno].stop, flowt[fno].step,
                                  flowt[fno].nsteps);
              load_clock(ADCIRC, fno);
          }
    break;

  case 122:
/* Line 1792 of yacc.c  */
#line 899 "gram.y"
    { 
		setflow = &g[curg].flowt[(int) (yyvsp[(3) - (3)].val)]; 
                elcirc_flowno = (int) (yyvsp[(3) - (3)].val);
		g[curg].curadc3d = (int) (yyvsp[(3) - (3)].val); 
	}
    break;

  case 123:
/* Line 1792 of yacc.c  */
#line 904 "gram.y"
    {}
    break;

  case 124:
/* Line 1792 of yacc.c  */
#line 905 "gram.y"
    { curadc3d = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 125:
/* Line 1792 of yacc.c  */
#line 907 "gram.y"
    {
		ReplaceElcircGrid((int) (yyvsp[(3) - (4)].val), (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 126:
/* Line 1792 of yacc.c  */
#line 911 "gram.y"
    { 
		ReadElcircRegionFile((int) (yyvsp[(3) - (5)].val), (char *) (yyvsp[(5) - (5)].str));
	}
    break;

  case 127:
/* Line 1792 of yacc.c  */
#line 915 "gram.y"
    {
		ReadElcirc((int) (yyvsp[(3) - (12)].val), (char *) (yyvsp[(4) - (12)].str), (int) (yyvsp[(12) - (12)].val) - 1, (int) (yyvsp[(6) - (12)].val), (int) (yyvsp[(8) - (12)].val), (int) (yyvsp[(10) - (12)].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[(3) - (12)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 128:
/* Line 1792 of yacc.c  */
#line 924 "gram.y"
    {
		ReadElcirc((int) (yyvsp[(3) - (13)].val), (char *) (yyvsp[(4) - (13)].str), (int) (yyvsp[(12) - (13)].val) - 1, (int) (yyvsp[(6) - (13)].val), (int) (yyvsp[(8) - (13)].val), (int) (yyvsp[(10) - (13)].val), 0, 0, 1);
                elcirc_flowno = (int) (yyvsp[(3) - (13)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 129:
/* Line 1792 of yacc.c  */
#line 933 "gram.y"
    {
		ReadElcircDepth((int) (yyvsp[(3) - (12)].val), (char *) (yyvsp[(4) - (12)].str), (char *) NULL, (double) (yyvsp[(12) - (12)].val), (int) (yyvsp[(6) - (12)].val), (int) (yyvsp[(8) - (12)].val), (int) (yyvsp[(10) - (12)].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[(3) - (12)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 130:
/* Line 1792 of yacc.c  */
#line 942 "gram.y"
    {
/* read at a given depth relative to the free surface */
		ReadElcircDepthFromFreeSurface((int) (yyvsp[(4) - (13)].val), (char *) (yyvsp[(5) - (13)].str), (char *) NULL, (double) (yyvsp[(13) - (13)].val), (int) (yyvsp[(7) - (13)].val), (int) (yyvsp[(9) - (13)].val), (int) (yyvsp[(11) - (13)].val), 0, 0, 0);
                elcirc_flowno = (int) (yyvsp[(4) - (13)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 131:
/* Line 1792 of yacc.c  */
#line 952 "gram.y"
    {
		ReadElcircSurf((int) (yyvsp[(4) - (11)].val), (char *) (yyvsp[(5) - (11)].str), 0, (int) (yyvsp[(7) - (11)].val), (int) (yyvsp[(9) - (11)].val), (int) (yyvsp[(11) - (11)].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[(4) - (11)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 132:
/* Line 1792 of yacc.c  */
#line 961 "gram.y"
    {
		ReadElcircSurf((int) (yyvsp[(4) - (12)].val), (char *) (yyvsp[(5) - (12)].str), 0, (int) (yyvsp[(7) - (12)].val), (int) (yyvsp[(9) - (12)].val), (int) (yyvsp[(11) - (12)].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[(4) - (12)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 133:
/* Line 1792 of yacc.c  */
#line 970 "gram.y"
    {
		ReadElcircSurf((int) (yyvsp[(4) - (11)].val), (char *) (yyvsp[(5) - (11)].str), 1, (int) (yyvsp[(7) - (11)].val), (int) (yyvsp[(9) - (11)].val), (int) (yyvsp[(11) - (11)].val), 0, 0.0, 0);
                elcirc_flowno = (int) (yyvsp[(4) - (11)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 134:
/* Line 1792 of yacc.c  */
#line 979 "gram.y"
    {
		ReadElcircSurf((int) (yyvsp[(4) - (12)].val), (char *) (yyvsp[(5) - (12)].str), 1, (int) (yyvsp[(7) - (12)].val), (int) (yyvsp[(9) - (12)].val), (int) (yyvsp[(11) - (12)].val), 0, 0.0, 1);
                elcirc_flowno = (int) (yyvsp[(4) - (12)].val);
                set_clock(0, flowt[elcirc_flowno].start, 
				flowt[elcirc_flowno].stop, flowt[elcirc_flowno].step,
                                  flowt[elcirc_flowno].nsteps);
                load_clock(ADCIRC, elcirc_flowno);
	}
    break;

  case 135:
/* Line 1792 of yacc.c  */
#line 987 "gram.y"
    { 
			curtrans = (int) (yyvsp[(4) - (4)].val); 
			settrans = &(trans[curtrans]); 
	}
    break;

  case 136:
/* Line 1792 of yacc.c  */
#line 991 "gram.y"
    { curtrans = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 137:
/* Line 1792 of yacc.c  */
#line 992 "gram.y"
    { setisol = &(g[curg].trans[curtrans].ip); }
    break;

  case 139:
/* Line 1792 of yacc.c  */
#line 993 "gram.y"
    {strcpy(settrans->uvname, (char *) (yyvsp[(4) - (4)].str));}
    break;

  case 140:
/* Line 1792 of yacc.c  */
#line 994 "gram.y"
    {strcpy(settrans->vvname, (char *) (yyvsp[(5) - (5)].str));}
    break;

  case 141:
/* Line 1792 of yacc.c  */
#line 995 "gram.y"
    { settrans->flowno = (int) (yyvsp[(4) - (4)].val);}
    break;

  case 142:
/* Line 1792 of yacc.c  */
#line 996 "gram.y"
    {strcpy(settrans->salname, (char *) (yyvsp[(4) - (4)].str));}
    break;

  case 143:
/* Line 1792 of yacc.c  */
#line 997 "gram.y"
    {strcpy(settrans->elevname, (char *) (yyvsp[(4) - (4)].str));}
    break;

  case 144:
/* Line 1792 of yacc.c  */
#line 998 "gram.y"
    { settrans->gno = (int) (yyvsp[(4) - (4)].val);}
    break;

  case 145:
/* Line 1792 of yacc.c  */
#line 999 "gram.y"
    { settrans->transgno = (int) (yyvsp[(5) - (5)].val);}
    break;

  case 146:
/* Line 1792 of yacc.c  */
#line 1000 "gram.y"
    { settrans->display = (int) (yyvsp[(5) - (5)].pset);}
    break;

  case 147:
/* Line 1792 of yacc.c  */
#line 1001 "gram.y"
    { g[curg].trans[curtrans].display = (int) (yyvsp[(4) - (4)].pset);}
    break;

  case 148:
/* Line 1792 of yacc.c  */
#line 1002 "gram.y"
    { g[curg].trans[curtrans].display_mag = (int) (yyvsp[(5) - (5)].pset);}
    break;

  case 149:
/* Line 1792 of yacc.c  */
#line 1003 "gram.y"
    { }
    break;

  case 150:
/* Line 1792 of yacc.c  */
#line 1005 "gram.y"
    {
		settrans->start = (int) (yyvsp[(4) - (8)].val);
		settrans->stop = (int) (yyvsp[(6) - (8)].val);
		settrans->skip = (int) (yyvsp[(8) - (8)].val);
	}
    break;

  case 151:
/* Line 1792 of yacc.c  */
#line 1010 "gram.y"
    { settrans->npts = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 152:
/* Line 1792 of yacc.c  */
#line 1011 "gram.y"
    { settrans->transtype = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 153:
/* Line 1792 of yacc.c  */
#line 1012 "gram.y"
    { AddTransNXY(settrans, (int) (yyvsp[(3) - (7)].val), (double) (yyvsp[(5) - (7)].val), (double) (yyvsp[(7) - (7)].val)); }
    break;

  case 154:
/* Line 1792 of yacc.c  */
#line 1013 "gram.y"
    { AddTransNode(settrans, (int) (yyvsp[(4) - (6)].val), (int) (yyvsp[(6) - (6)].val)); }
    break;

  case 155:
/* Line 1792 of yacc.c  */
#line 1015 "gram.y"
    {
		settrans->transtype = 0;
		strcpy(settrans->transname, (char *) (yyvsp[(5) - (5)].str));
	}
    break;

  case 156:
/* Line 1792 of yacc.c  */
#line 1020 "gram.y"
    {
		settrans->transtype = 1;
		settrans->npts = (int) (yyvsp[(4) - (12)].val);
		settrans->x1 = (double) (yyvsp[(6) - (12)].val);
		settrans->y1 = (double) (yyvsp[(8) - (12)].val);
		settrans->x2 = (double) (yyvsp[(10) - (12)].val);
		settrans->y2 = (double) (yyvsp[(12) - (12)].val);
	}
    break;

  case 157:
/* Line 1792 of yacc.c  */
#line 1029 "gram.y"
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 158:
/* Line 1792 of yacc.c  */
#line 1037 "gram.y"
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 159:
/* Line 1792 of yacc.c  */
#line 1045 "gram.y"
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 160:
/* Line 1792 of yacc.c  */
#line 1053 "gram.y"
    { 
		settrans->type = VECTOR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 161:
/* Line 1792 of yacc.c  */
#line 1061 "gram.y"
    { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 0);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 162:
/* Line 1792 of yacc.c  */
#line 1069 "gram.y"
    { 
		settrans->type = SCALAR;
		ReadNewTrans(settrans, 1);
                trans[curtrans].active = ON;
                set_clock(0, trans[curtrans].tstart, trans[curtrans].tstop, trans[curtrans].tstep, trans[curtrans].nsteps);
                load_clock(TRANSECT, curtrans);
	}
    break;

  case 163:
/* Line 1792 of yacc.c  */
#line 1077 "gram.y"
    {
		elcircmarker = (int) (yyvsp[(4) - (4)].val);
		setadc3d = &adc3d[(int) (yyvsp[(4) - (4)].val)];
		setflow3d = &g[curg].flow3d[(int) (yyvsp[(4) - (4)].val)];
	}
    break;

  case 164:
/* Line 1792 of yacc.c  */
#line 1083 "gram.y"
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 165:
/* Line 1792 of yacc.c  */
#line 1087 "gram.y"
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 166:
/* Line 1792 of yacc.c  */
#line 1091 "gram.y"
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 167:
/* Line 1792 of yacc.c  */
#line 1095 "gram.y"
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 168:
/* Line 1792 of yacc.c  */
#line 1099 "gram.y"
    {
		strcpy(setadc3d->datafile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 169:
/* Line 1792 of yacc.c  */
#line 1103 "gram.y"
    {
		strcpy(setadc3d->elevfile, (char *) (yyvsp[(4) - (4)].str));
	}
    break;

  case 170:
/* Line 1792 of yacc.c  */
#line 1107 "gram.y"
    {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[(5) - (11)].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[(7) - (11)].val), (int) (yyvsp[(9) - (11)].val), 0);
	}
    break;

  case 171:
/* Line 1792 of yacc.c  */
#line 1114 "gram.y"
    {
		setadc3d->loctype = NODE;
		setadc3d->loctype = 1;
		setadc3d->node = (int) (yyvsp[(5) - (12)].val);
		ReadNodeDataNew(elcircmarker, (int) (yyvsp[(7) - (12)].val), (int) (yyvsp[(9) - (12)].val), 1);
	}
    break;

  case 172:
/* Line 1792 of yacc.c  */
#line 1121 "gram.y"
    {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[(5) - (13)].val);
		setadc3d->y = (double) (yyvsp[(7) - (13)].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[(9) - (13)].val), (int) (yyvsp[(11) - (13)].val), 0);
	}
    break;

  case 173:
/* Line 1792 of yacc.c  */
#line 1129 "gram.y"
    {
		setadc3d->loctype = XY;
		setadc3d->loctype = 0;
		setadc3d->x = (double) (yyvsp[(5) - (14)].val);
		setadc3d->y = (double) (yyvsp[(7) - (14)].val);
		ReadXYDataNew(elcircmarker, (int) (yyvsp[(9) - (14)].val), (int) (yyvsp[(11) - (14)].val), 1);
	}
    break;

  case 174:
/* Line 1792 of yacc.c  */
#line 1136 "gram.y"
    { setisol = &(g[curg].salip); }
    break;

  case 176:
/* Line 1792 of yacc.c  */
#line 1137 "gram.y"
    { setisol = &(g[curg].salip); }
    break;

  case 178:
/* Line 1792 of yacc.c  */
#line 1138 "gram.y"
    { setisol = &(g[curg].velmagip); }
    break;

  case 180:
/* Line 1792 of yacc.c  */
#line 1139 "gram.y"
    { setprops = &(setflow3d->p); }
    break;

  case 182:
/* Line 1792 of yacc.c  */
#line 1141 "gram.y"
    {
		setflow3d->precx = (int) (yyvsp[(4) - (6)].val);
		setflow3d->precy = (int) (yyvsp[(6) - (6)].val);
	}
    break;

  case 183:
/* Line 1792 of yacc.c  */
#line 1145 "gram.y"
    { setflow3d->attach = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 184:
/* Line 1792 of yacc.c  */
#line 1146 "gram.y"
    { setflow3d->loctype = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 185:
/* Line 1792 of yacc.c  */
#line 1147 "gram.y"
    { setflow3d->display_marker = (int) (yyvsp[(5) - (5)].pset); }
    break;

  case 186:
/* Line 1792 of yacc.c  */
#line 1148 "gram.y"
    { setflow3d->display = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 187:
/* Line 1792 of yacc.c  */
#line 1149 "gram.y"
    { setflow3d->p.color = (yyvsp[(4) - (4)].val); }
    break;

  case 188:
/* Line 1792 of yacc.c  */
#line 1150 "gram.y"
    { setflow3d->p.linew = (yyvsp[(4) - (4)].val); }
    break;

  case 189:
/* Line 1792 of yacc.c  */
#line 1151 "gram.y"
    { setflow3d->p.fillcol = (yyvsp[(5) - (5)].val); }
    break;

  case 190:
/* Line 1792 of yacc.c  */
#line 1153 "gram.y"
    { 
		setflow3d->wx1 = (double) (yyvsp[(4) - (10)].val); 
		setflow3d->wy1 = (double) (yyvsp[(6) - (10)].val); 
		setflow3d->wx2 = (double) (yyvsp[(8) - (10)].val); 
		setflow3d->wy2 = (double) (yyvsp[(10) - (10)].val); 
	}
    break;

  case 191:
/* Line 1792 of yacc.c  */
#line 1160 "gram.y"
    { 
		setflow3d->vx = (double) (yyvsp[(4) - (6)].val); 
		setflow3d->vy = (double) (yyvsp[(6) - (6)].val); 
	}
    break;

  case 192:
/* Line 1792 of yacc.c  */
#line 1165 "gram.y"
    { 
		setflow3d->locx = (double) (yyvsp[(4) - (6)].val); 
		setflow3d->locy = (double) (yyvsp[(6) - (6)].val); 
	}
    break;

  case 193:
/* Line 1792 of yacc.c  */
#line 1170 "gram.y"
    { 
		setflow3d->x = (double) (yyvsp[(4) - (6)].val); 
		setflow3d->y = (double) (yyvsp[(6) - (6)].val); 
	}
    break;

  case 194:
/* Line 1792 of yacc.c  */
#line 1177 "gram.y"
    { setflow->display = (yyvsp[(2) - (2)].pset);  }
    break;

  case 195:
/* Line 1792 of yacc.c  */
#line 1178 "gram.y"
    { setflow->display_elev = (yyvsp[(3) - (3)].pset); }
    break;

  case 196:
/* Line 1792 of yacc.c  */
#line 1179 "gram.y"
    { setflow->display_elevdepth = (yyvsp[(4) - (4)].pset); }
    break;

  case 197:
/* Line 1792 of yacc.c  */
#line 1180 "gram.y"
    { setflow->display_maxelevval = (yyvsp[(4) - (4)].pset); }
    break;

  case 198:
/* Line 1792 of yacc.c  */
#line 1181 "gram.y"
    { setflow->display_maxelev = (yyvsp[(4) - (4)].pset); }
    break;

  case 199:
/* Line 1792 of yacc.c  */
#line 1182 "gram.y"
    { setflow->display_amp = (yyvsp[(4) - (4)].pset); }
    break;

  case 200:
/* Line 1792 of yacc.c  */
#line 1183 "gram.y"
    { setflow->display_phase = (yyvsp[(4) - (4)].pset); }
    break;

  case 201:
/* Line 1792 of yacc.c  */
#line 1184 "gram.y"
    { setflow->display_elevmarkers = (yyvsp[(4) - (4)].pset); }
    break;

  case 202:
/* Line 1792 of yacc.c  */
#line 1185 "gram.y"
    { setflow->display_mag = (yyvsp[(4) - (4)].pset); }
    break;

  case 203:
/* Line 1792 of yacc.c  */
#line 1186 "gram.y"
    { setflow->display_wind = (yyvsp[(4) - (4)].pset); }
    break;

  case 204:
/* Line 1792 of yacc.c  */
#line 1187 "gram.y"
    { setflow->display_inun = (yyvsp[(2) - (3)].pset); }
    break;

  case 205:
/* Line 1792 of yacc.c  */
#line 1188 "gram.y"
    { setflow->p.color = (yyvsp[(2) - (2)].val); }
    break;

  case 206:
/* Line 1792 of yacc.c  */
#line 1189 "gram.y"
    { setisol = &(setflow->elevip); }
    break;

  case 208:
/* Line 1792 of yacc.c  */
#line 1190 "gram.y"
    { setisol = &(setflow->maxelevip); }
    break;

  case 210:
/* Line 1792 of yacc.c  */
#line 1191 "gram.y"
    { setisol = &(setflow->ampip); }
    break;

  case 212:
/* Line 1792 of yacc.c  */
#line 1192 "gram.y"
    { setisol = &(setflow->phaseip); }
    break;

  case 214:
/* Line 1792 of yacc.c  */
#line 1193 "gram.y"
    { setisol = &(setflow->magip); }
    break;

  case 216:
/* Line 1792 of yacc.c  */
#line 1194 "gram.y"
    { setflow->flowfreq = (yyvsp[(3) - (3)].val); }
    break;

  case 217:
/* Line 1792 of yacc.c  */
#line 1195 "gram.y"
    { setflow->freq = (yyvsp[(2) - (2)].val); }
    break;

  case 218:
/* Line 1792 of yacc.c  */
#line 1196 "gram.y"
    { setelevmarker = &(setflow->em[(int) (yyvsp[(2) - (2)].val)]); }
    break;

  case 219:
/* Line 1792 of yacc.c  */
#line 1197 "gram.y"
    { setflow->sample = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 220:
/* Line 1792 of yacc.c  */
#line 1198 "gram.y"
    { ReadSampleFlow(setflow, (char *) (yyvsp[(4) - (4)].str)); }
    break;

  case 221:
/* Line 1792 of yacc.c  */
#line 1199 "gram.y"
    { SetMinSampleFlow(elcirc_flowno, (double) (yyvsp[(4) - (4)].val)); }
    break;

  case 222:
/* Line 1792 of yacc.c  */
#line 1200 "gram.y"
    { ReadSampleFlowXY(setflow, (char *) (yyvsp[(5) - (5)].str)); }
    break;

  case 223:
/* Line 1792 of yacc.c  */
#line 1201 "gram.y"
    { setflow->samptype = XY; }
    break;

  case 224:
/* Line 1792 of yacc.c  */
#line 1202 "gram.y"
    { setflow->samptype = NODE; }
    break;

  case 225:
/* Line 1792 of yacc.c  */
#line 1203 "gram.y"
    { AddSampleFlowNode(setflow, (int) (yyvsp[(4) - (4)].val) - 1); }
    break;

  case 226:
/* Line 1792 of yacc.c  */
#line 1204 "gram.y"
    { AddSampleFlowElem(setflow, (int) (yyvsp[(4) - (4)].val) - 1); }
    break;

  case 227:
/* Line 1792 of yacc.c  */
#line 1205 "gram.y"
    { AddSampleFlowXY(setflow, (double) (yyvsp[(4) - (6)].val), (double) (yyvsp[(6) - (6)].val)); }
    break;

  case 228:
/* Line 1792 of yacc.c  */
#line 1206 "gram.y"
    { DeleteSampleFlowNode(setflow, (int) (yyvsp[(5) - (5)].val) - 1); }
    break;

  case 229:
/* Line 1792 of yacc.c  */
#line 1207 "gram.y"
    { DeleteSampleFlowElem(setflow, (int) (yyvsp[(5) - (5)].val) - 1); }
    break;

  case 230:
/* Line 1792 of yacc.c  */
#line 1208 "gram.y"
    { DeleteSampleFlowXY(setflow, (double) (yyvsp[(5) - (7)].val), (double) (yyvsp[(7) - (7)].val)); }
    break;

  case 231:
/* Line 1792 of yacc.c  */
#line 1212 "gram.y"
    { setelevmarker = &(setflow->em[(int) (yyvsp[(3) - (3)].val)]); }
    break;

  case 232:
/* Line 1792 of yacc.c  */
#line 1213 "gram.y"
    { setelevmarker->active = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 233:
/* Line 1792 of yacc.c  */
#line 1214 "gram.y"
    { setelevmarker->type = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 234:
/* Line 1792 of yacc.c  */
#line 1215 "gram.y"
    { setelevmarker->display = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 235:
/* Line 1792 of yacc.c  */
#line 1216 "gram.y"
    { setprops = &(setelevmarker->p); }
    break;

  case 237:
/* Line 1792 of yacc.c  */
#line 1217 "gram.y"
    { setelevmarker->loctype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 238:
/* Line 1792 of yacc.c  */
#line 1218 "gram.y"
    { setelevmarker->node = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 239:
/* Line 1792 of yacc.c  */
#line 1220 "gram.y"
    { 
		setelevmarker->locx = (double) (yyvsp[(3) - (5)].val); 
		setelevmarker->locy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 240:
/* Line 1792 of yacc.c  */
#line 1225 "gram.y"
    { 
		setelevmarker->emin = (double) (yyvsp[(3) - (5)].val); 
		setelevmarker->emax = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 241:
/* Line 1792 of yacc.c  */
#line 1236 "gram.y"
    { setgrid = &g[curg].grid[(int) (yyvsp[(2) - (2)].val)]; }
    break;

  case 242:
/* Line 1792 of yacc.c  */
#line 1237 "gram.y"
    { setgrid = &g[curg].grid[(int) (yyvsp[(3) - (3)].val)]; }
    break;

  case 243:
/* Line 1792 of yacc.c  */
#line 1238 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display = (yyvsp[(3) - (3)].pset); }
    break;

  case 244:
/* Line 1792 of yacc.c  */
#line 1239 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_bath = (yyvsp[(4) - (4)].pset); }
    break;

  case 245:
/* Line 1792 of yacc.c  */
#line 1240 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_courant = (yyvsp[(4) - (6)].pset); }
    break;

  case 246:
/* Line 1792 of yacc.c  */
#line 1241 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_courantn = (yyvsp[(5) - (7)].pset); }
    break;

  case 247:
/* Line 1792 of yacc.c  */
#line 1242 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_boundary = (yyvsp[(4) - (4)].pset); }
    break;

  case 248:
/* Line 1792 of yacc.c  */
#line 1243 "gram.y"
    { if (checkptr(setgrid, f_string)) setprops = &(setgrid->p); }
    break;

  case 250:
/* Line 1792 of yacc.c  */
#line 1244 "gram.y"
    { if (checkptr(setgrid, f_string)) setprops = &(setgrid->bp); }
    break;

  case 252:
/* Line 1792 of yacc.c  */
#line 1245 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_nodes = (yyvsp[(4) - (4)].pset); }
    break;

  case 253:
/* Line 1792 of yacc.c  */
#line 1246 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_elements = (yyvsp[(4) - (4)].pset); }
    break;

  case 254:
/* Line 1792 of yacc.c  */
#line 1247 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_depths = (yyvsp[(4) - (4)].pset); }
    break;

  case 255:
/* Line 1792 of yacc.c  */
#line 1248 "gram.y"
    { if (checkptr(setgrid, f_string)) setgrid->display_gridf = (yyvsp[(3) - (3)].pset); }
    break;

  case 256:
/* Line 1792 of yacc.c  */
#line 1249 "gram.y"
    { if (checkptr(setgrid, f_string)) setisol = &(setgrid->ip); }
    break;

  case 258:
/* Line 1792 of yacc.c  */
#line 1251 "gram.y"
    { 
		autoscale_grid((int) (yyvsp[(1) - (2)].pset), g[(int) (yyvsp[(1) - (2)].pset)].curgrid); 
		set_defaults((int) (yyvsp[(1) - (2)].pset));
	}
    break;

  case 259:
/* Line 1792 of yacc.c  */
#line 1256 "gram.y"
    { 
		autoscale_grid(curg, g[curg].curgrid); 
		set_defaults(curg);
	}
    break;

  case 260:
/* Line 1792 of yacc.c  */
#line 1261 "gram.y"
    {
		extern int readgridfile;
		readgrid((int) (yyvsp[(3) - (4)].val), (char *) (yyvsp[(4) - (4)].str));
		readgridfile = 1;
	}
    break;

  case 261:
/* Line 1792 of yacc.c  */
#line 1269 "gram.y"
    { setdrogs = &g[curg].drogues[(int) (yyvsp[(3) - (3)].val)]; }
    break;

  case 262:
/* Line 1792 of yacc.c  */
#line 1270 "gram.y"
    {  setdrogs->display = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 263:
/* Line 1792 of yacc.c  */
#line 1271 "gram.y"
    { setprops = &(setdrogs->p); }
    break;

  case 265:
/* Line 1792 of yacc.c  */
#line 1273 "gram.y"
    {
    	if (!readdrogues(curdrog, (char *) (yyvsp[(3) - (3)].str), -1, 0, 0)) {
        	fprintf(stderr, "Error reading file %s", (char *) (yyvsp[(3) - (3)].str));
    	} else {
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
    	}
	}
    break;

  case 266:
/* Line 1792 of yacc.c  */
#line 1284 "gram.y"
    {
		readdrogues(0, (char *) (yyvsp[(3) - (9)].str), (int) (yyvsp[(5) - (9)].val), (int) (yyvsp[(7) - (9)].val), (int) (yyvsp[(9) - (9)].val));
        	set_clock(0, drogues[curdrog].start, drogues[curdrog].stop,
                  drogues[curdrog].step,
                  drogues[curdrog].nsteps);
        	load_clock(DROGUES, curdrog);
	}
    break;

  case 267:
/* Line 1792 of yacc.c  */
#line 1294 "gram.y"
    { setisol->lactive = (yyvsp[(3) - (3)].pset); }
    break;

  case 268:
/* Line 1792 of yacc.c  */
#line 1295 "gram.y"
    { setisol->layout = (yyvsp[(4) - (4)].pset); }
    break;

  case 269:
/* Line 1792 of yacc.c  */
#line 1296 "gram.y"
    { setisol->llabels = (yyvsp[(3) - (4)].pset); }
    break;

  case 270:
/* Line 1792 of yacc.c  */
#line 1297 "gram.y"
    { setisol->frame = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 271:
/* Line 1792 of yacc.c  */
#line 1298 "gram.y"
    { setisol->framecol = (int) (yyvsp[(5) - (5)].val); }
    break;

  case 272:
/* Line 1792 of yacc.c  */
#line 1299 "gram.y"
    { setisol->nisol = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 273:
/* Line 1792 of yacc.c  */
#line 1300 "gram.y"
    { setisol->type = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 274:
/* Line 1792 of yacc.c  */
#line 1301 "gram.y"
    { setisol->isoltype = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 275:
/* Line 1792 of yacc.c  */
#line 1302 "gram.y"
    { setisol->visflag = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 276:
/* Line 1792 of yacc.c  */
#line 1303 "gram.y"
    { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
    break;

  case 277:
/* Line 1792 of yacc.c  */
#line 1311 "gram.y"
    { 
		int i;
		setisol->writeflag = 0; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		setisol->wname[0] = 0;
	}
    break;

  case 278:
/* Line 1792 of yacc.c  */
#line 1319 "gram.y"
    { 
		int i;
		setisol->writeflag = 1; 
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
		strncpy(setisol->wname, (char *) (yyvsp[(4) - (4)].str), 1023);
	}
    break;

  case 279:
/* Line 1792 of yacc.c  */
#line 1327 "gram.y"
    { 
		int i;
		setisol->writeflag = 1; 
		strncpy(setisol->wname, (char *) (yyvsp[(6) - (6)].str), 1023);
		for (i=0;i<MAXISOLINES;i++) {
    		    setisol->writelevel[i] = 0;
		}
    		setisol->writelevel[(int) (yyvsp[(4) - (6)].val)] = 1;
	}
    break;

  case 280:
/* Line 1792 of yacc.c  */
#line 1336 "gram.y"
    { 
    		setisol->writelevel[(int) (yyvsp[(4) - (4)].val)] = 1;
	}
    break;

  case 281:
/* Line 1792 of yacc.c  */
#line 1339 "gram.y"
    {
            setisol->cis[0] = (double) (yyvsp[(3) - (5)].val);
            setisol->cint = (double) (yyvsp[(5) - (5)].val);
	    if (setisol->isoltype == 0) {
		int i;
		for (i=1;i< 16;i++) {
		    setisol->cis[i] = setisol->cis[0] + i * setisol->cint;
		}
	    }
        }
    break;

  case 282:
/* Line 1792 of yacc.c  */
#line 1349 "gram.y"
    { setisol->cmin = (yyvsp[(3) - (5)].val); setisol->cmax = (yyvsp[(5) - (5)].val); }
    break;

  case 283:
/* Line 1792 of yacc.c  */
#line 1350 "gram.y"
    { setisol->cis[(int) (yyvsp[(2) - (4)].val)] = (yyvsp[(4) - (4)].val); }
    break;

  case 284:
/* Line 1792 of yacc.c  */
#line 1351 "gram.y"
    { setisol->loctype = (yyvsp[(4) - (4)].pset); }
    break;

  case 285:
/* Line 1792 of yacc.c  */
#line 1352 "gram.y"
    {
            setisol->x = (double) (yyvsp[(3) - (5)].val);
            setisol->y = (double) (yyvsp[(5) - (5)].val);
        }
    break;

  case 286:
/* Line 1792 of yacc.c  */
#line 1356 "gram.y"
    { setisol->color[(int) (yyvsp[(2) - (4)].val)] = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 287:
/* Line 1792 of yacc.c  */
#line 1357 "gram.y"
    { setisol->linew[(int) (yyvsp[(2) - (4)].val)] = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 288:
/* Line 1792 of yacc.c  */
#line 1358 "gram.y"
    { setisol->lines[(int) (yyvsp[(2) - (4)].val)] = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 289:
/* Line 1792 of yacc.c  */
#line 1359 "gram.y"
    { setprops = &(setisol->p); }
    break;

  case 291:
/* Line 1792 of yacc.c  */
#line 1360 "gram.y"
    { setprops = &(setisol->p); }
    break;

  case 293:
/* Line 1792 of yacc.c  */
#line 1361 "gram.y"
    {
            setisol->xlen = (yyvsp[(4) - (6)].val);
            setisol->ylen = (yyvsp[(6) - (6)].val);
        }
    break;

  case 294:
/* Line 1792 of yacc.c  */
#line 1365 "gram.y"
    {
            setisol->xgap = (yyvsp[(4) - (6)].val);
            setisol->ygap = (yyvsp[(6) - (6)].val);
        }
    break;

  case 295:
/* Line 1792 of yacc.c  */
#line 1372 "gram.y"
    { sethistbox = &g[curg].hbox[(int) (yyvsp[(3) - (3)].val)]; }
    break;

  case 296:
/* Line 1792 of yacc.c  */
#line 1373 "gram.y"
    { setprops = &(sethistbox->p); }
    break;

  case 298:
/* Line 1792 of yacc.c  */
#line 1375 "gram.y"
    {
		sethistbox->precx = (int) (yyvsp[(3) - (5)].val);
		sethistbox->precy = (int) (yyvsp[(5) - (5)].val);
	}
    break;

  case 299:
/* Line 1792 of yacc.c  */
#line 1379 "gram.y"
    { sethistbox->attach = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 300:
/* Line 1792 of yacc.c  */
#line 1380 "gram.y"
    { sethistbox->xtickm = (double ) (yyvsp[(3) - (5)].val); sethistbox->ytickm = (double ) (yyvsp[(5) - (5)].val);}
    break;

  case 301:
/* Line 1792 of yacc.c  */
#line 1381 "gram.y"
    { sethistbox->loctype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 302:
/* Line 1792 of yacc.c  */
#line 1382 "gram.y"
    { sethistbox->display_marker = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 303:
/* Line 1792 of yacc.c  */
#line 1383 "gram.y"
    { sethistbox->display = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 304:
/* Line 1792 of yacc.c  */
#line 1384 "gram.y"
    { sethistbox->adcirc[(int) (yyvsp[(4) - (5)].val)] = (int) (yyvsp[(5) - (5)].pset) == TRUEP; }
    break;

  case 305:
/* Line 1792 of yacc.c  */
#line 1385 "gram.y"
    { sethistbox->ap[(int) (yyvsp[(4) - (6)].val)].color = (int) (yyvsp[(6) - (6)].val); }
    break;

  case 306:
/* Line 1792 of yacc.c  */
#line 1386 "gram.y"
    { sethistbox->p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 307:
/* Line 1792 of yacc.c  */
#line 1387 "gram.y"
    { sethistbox->p.linew = (yyvsp[(3) - (3)].val); }
    break;

  case 308:
/* Line 1792 of yacc.c  */
#line 1388 "gram.y"
    { sethistbox->p.fillcol = (yyvsp[(4) - (4)].val); }
    break;

  case 309:
/* Line 1792 of yacc.c  */
#line 1389 "gram.y"
    { read_hist((int) (yyvsp[(3) - (4)].val), TIME, (char *) (yyvsp[(4) - (4)].str)); }
    break;

  case 310:
/* Line 1792 of yacc.c  */
#line 1390 "gram.y"
    { sethistbox->thist = (int) (yyvsp[(4) - (4)].pset) == TRUEP; }
    break;

  case 311:
/* Line 1792 of yacc.c  */
#line 1391 "gram.y"
    { sethistbox->hp.color = (int) (yyvsp[(5) - (5)].val); }
    break;

  case 312:
/* Line 1792 of yacc.c  */
#line 1393 "gram.y"
    { 
		sethistbox->wx1 = (double) (yyvsp[(3) - (9)].val); 
		sethistbox->wy1 = (double) (yyvsp[(5) - (9)].val); 
		sethistbox->wx2 = (double) (yyvsp[(7) - (9)].val); 
		sethistbox->wy2 = (double) (yyvsp[(9) - (9)].val); 
	}
    break;

  case 313:
/* Line 1792 of yacc.c  */
#line 1400 "gram.y"
    { 
		sethistbox->vx = (double) (yyvsp[(3) - (5)].val); 
		sethistbox->vy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 314:
/* Line 1792 of yacc.c  */
#line 1405 "gram.y"
    { 
		sethistbox->locx = (double) (yyvsp[(3) - (5)].val); 
		sethistbox->locy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 315:
/* Line 1792 of yacc.c  */
#line 1410 "gram.y"
    { 
		sethistbox->x = (double) (yyvsp[(3) - (5)].val); 
		sethistbox->y = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 316:
/* Line 1792 of yacc.c  */
#line 1417 "gram.y"
    { setzoombox = &g[curg].zbox[(int) (yyvsp[(3) - (3)].val)]; }
    break;

  case 317:
/* Line 1792 of yacc.c  */
#line 1418 "gram.y"
    { setprops = &(setzoombox->p); }
    break;

  case 319:
/* Line 1792 of yacc.c  */
#line 1420 "gram.y"
    {
		setzoombox->precx = (int) (yyvsp[(3) - (5)].val);
		setzoombox->precy = (int) (yyvsp[(5) - (5)].val);
	}
    break;

  case 320:
/* Line 1792 of yacc.c  */
#line 1424 "gram.y"
    { setzoombox->attach = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 321:
/* Line 1792 of yacc.c  */
#line 1425 "gram.y"
    { setzoombox->loctype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 322:
/* Line 1792 of yacc.c  */
#line 1426 "gram.y"
    { setzoombox->display_marker = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 323:
/* Line 1792 of yacc.c  */
#line 1427 "gram.y"
    { setzoombox->active = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 324:
/* Line 1792 of yacc.c  */
#line 1428 "gram.y"
    { setzoombox->display = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 325:
/* Line 1792 of yacc.c  */
#line 1429 "gram.y"
    { setzoombox->expand = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 326:
/* Line 1792 of yacc.c  */
#line 1430 "gram.y"
    { setzoombox->expand = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 327:
/* Line 1792 of yacc.c  */
#line 1431 "gram.y"
    { setzoombox->p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 328:
/* Line 1792 of yacc.c  */
#line 1432 "gram.y"
    { setzoombox->p.linew = (yyvsp[(3) - (3)].val); }
    break;

  case 329:
/* Line 1792 of yacc.c  */
#line 1433 "gram.y"
    { setzoombox->p.fillcol = (yyvsp[(4) - (4)].val); }
    break;

  case 330:
/* Line 1792 of yacc.c  */
#line 1435 "gram.y"
    { 
		setzoombox->wx1 = (double) (yyvsp[(3) - (9)].val); 
		setzoombox->wy1 = (double) (yyvsp[(5) - (9)].val); 
		setzoombox->wx2 = (double) (yyvsp[(7) - (9)].val); 
		setzoombox->wy2 = (double) (yyvsp[(9) - (9)].val); 
	}
    break;

  case 331:
/* Line 1792 of yacc.c  */
#line 1442 "gram.y"
    { 
		setzoombox->vx = (double) (yyvsp[(3) - (5)].val); 
		setzoombox->vy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 332:
/* Line 1792 of yacc.c  */
#line 1447 "gram.y"
    { 
		setzoombox->locx = (double) (yyvsp[(3) - (5)].val); 
		setzoombox->locy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 333:
/* Line 1792 of yacc.c  */
#line 1452 "gram.y"
    { 
		setzoombox->x = (double) (yyvsp[(3) - (5)].val); 
		setzoombox->y = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 334:
/* Line 1792 of yacc.c  */
#line 1459 "gram.y"
    { g[curg].wl.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 335:
/* Line 1792 of yacc.c  */
#line 1460 "gram.y"
    { g[curg].wl.len = (yyvsp[(3) - (3)].val); }
    break;

  case 336:
/* Line 1792 of yacc.c  */
#line 1461 "gram.y"
    { g[curg].wl.scale = (yyvsp[(3) - (3)].val); }
    break;

  case 337:
/* Line 1792 of yacc.c  */
#line 1462 "gram.y"
    { g[curg].wl.p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 338:
/* Line 1792 of yacc.c  */
#line 1463 "gram.y"
    { g[curg].wl.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 339:
/* Line 1792 of yacc.c  */
#line 1465 "gram.y"
    { 
		g[curg].wl.x = (yyvsp[(3) - (5)].val);
		g[curg].wl.y = (yyvsp[(5) - (5)].val);
	}
    break;

  case 340:
/* Line 1792 of yacc.c  */
#line 1470 "gram.y"
    { 
		g[curg].wl.units = (yyvsp[(3) - (3)].pset);
		switch (g[curg].wl.units) {
		case MM: g[curg].wl.unitfac = 0.001; break;
		case CM: g[curg].wl.unitfac = 0.01; break;
		case M: g[curg].wl.unitfac = 1.0; break;
		case KM: g[curg].wl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
    break;

  case 341:
/* Line 1792 of yacc.c  */
#line 1480 "gram.y"
    { }
    break;

  case 342:
/* Line 1792 of yacc.c  */
#line 1485 "gram.y"
    { g[curg].vl.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 343:
/* Line 1792 of yacc.c  */
#line 1486 "gram.y"
    { g[curg].vl.len = (yyvsp[(3) - (3)].val); }
    break;

  case 344:
/* Line 1792 of yacc.c  */
#line 1487 "gram.y"
    { g[curg].vl.scale = (yyvsp[(3) - (3)].val); }
    break;

  case 345:
/* Line 1792 of yacc.c  */
#line 1488 "gram.y"
    { g[curg].vl.p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 346:
/* Line 1792 of yacc.c  */
#line 1489 "gram.y"
    { g[curg].vl.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 347:
/* Line 1792 of yacc.c  */
#line 1491 "gram.y"
    { 
		g[curg].vl.x = (yyvsp[(3) - (5)].val);
		g[curg].vl.y = (yyvsp[(5) - (5)].val);
	}
    break;

  case 348:
/* Line 1792 of yacc.c  */
#line 1496 "gram.y"
    { 
		g[curg].vl.units = (yyvsp[(3) - (3)].pset);
		switch (g[curg].vl.units) {
		case MM: g[curg].vl.unitfac = 0.001; break;
		case CM: g[curg].vl.unitfac = 0.01; break;
		case M: g[curg].vl.unitfac = 1.0; break;
		case KM: g[curg].vl.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for velocity scale\n"); break;
		}
	}
    break;

  case 349:
/* Line 1792 of yacc.c  */
#line 1506 "gram.y"
    { }
    break;

  case 350:
/* Line 1792 of yacc.c  */
#line 1510 "gram.y"
    { g[curg].mapscale.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 351:
/* Line 1792 of yacc.c  */
#line 1511 "gram.y"
    { g[curg].mapscale.len = (yyvsp[(3) - (3)].val); }
    break;

  case 352:
/* Line 1792 of yacc.c  */
#line 1512 "gram.y"
    { g[curg].mapscale.p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 353:
/* Line 1792 of yacc.c  */
#line 1513 "gram.y"
    { g[curg].mapscale.scale = (yyvsp[(3) - (3)].val); }
    break;

  case 354:
/* Line 1792 of yacc.c  */
#line 1514 "gram.y"
    { g[curg].mapscale.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 355:
/* Line 1792 of yacc.c  */
#line 1516 "gram.y"
    { 
		g[curg].mapscale.x = (yyvsp[(3) - (5)].val);
		g[curg].mapscale.y = (yyvsp[(5) - (5)].val);
	}
    break;

  case 356:
/* Line 1792 of yacc.c  */
#line 1521 "gram.y"
    { 
		g[curg].mapscale.units = (yyvsp[(3) - (3)].pset);
		switch (g[curg].mapscale.units) {
		case MM: g[curg].mapscale.unitfac = 0.001; break;
		case CM: g[curg].mapscale.unitfac = 0.01; break;
		case M: g[curg].mapscale.unitfac = 1.0; break;
		case KM: g[curg].mapscale.unitfac = 1000.0; break;
		default: fprintf(stderr, "Unknown units for mapscape scale\n"); break;
		}
	}
    break;

  case 357:
/* Line 1792 of yacc.c  */
#line 1531 "gram.y"
    { }
    break;

  case 358:
/* Line 1792 of yacc.c  */
#line 1535 "gram.y"
    { g[curg].tidalclock.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 359:
/* Line 1792 of yacc.c  */
#line 1536 "gram.y"
    { g[curg].tidalclock.p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 360:
/* Line 1792 of yacc.c  */
#line 1537 "gram.y"
    { g[curg].tidalclock.p.fillcol = (yyvsp[(4) - (4)].val); }
    break;

  case 361:
/* Line 1792 of yacc.c  */
#line 1538 "gram.y"
    { g[curg].tidalclock.total_time = (yyvsp[(4) - (4)].val); }
    break;

  case 362:
/* Line 1792 of yacc.c  */
#line 1539 "gram.y"
    { g[curg].tidalclock.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 363:
/* Line 1792 of yacc.c  */
#line 1541 "gram.y"
    { 
		g[curg].tidalclock.x = (yyvsp[(3) - (5)].val);
		g[curg].tidalclock.y = (yyvsp[(5) - (5)].val);
	}
    break;

  case 364:
/* Line 1792 of yacc.c  */
#line 1545 "gram.y"
    { }
    break;

  case 365:
/* Line 1792 of yacc.c  */
#line 1549 "gram.y"
    { g[curg].timeinfo.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 366:
/* Line 1792 of yacc.c  */
#line 1551 "gram.y"
    { 
		strcpy(g[curg].timeinfo.start, (char *) (yyvsp[(3) - (3)].str)); 
		time_info_start(curg);
	}
    break;

  case 367:
/* Line 1792 of yacc.c  */
#line 1556 "gram.y"
    {
	    g[curg].timeinfo.x = (yyvsp[(2) - (4)].val);
	    g[curg].timeinfo.y = (yyvsp[(4) - (4)].val);
	}
    break;

  case 368:
/* Line 1792 of yacc.c  */
#line 1560 "gram.y"
    { g[curg].timeinfo.loctype = (yyvsp[(3) - (3)].pset); }
    break;

  case 369:
/* Line 1792 of yacc.c  */
#line 1561 "gram.y"
    { g[curg].timeinfo.linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 370:
/* Line 1792 of yacc.c  */
#line 1562 "gram.y"
    { g[curg].timeinfo.color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 371:
/* Line 1792 of yacc.c  */
#line 1563 "gram.y"
    { g[curg].timeinfo.rot = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 372:
/* Line 1792 of yacc.c  */
#line 1564 "gram.y"
    { g[curg].timeinfo.font = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 373:
/* Line 1792 of yacc.c  */
#line 1565 "gram.y"
    { g[curg].timeinfo.just = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 374:
/* Line 1792 of yacc.c  */
#line 1566 "gram.y"
    { g[curg].timeinfo.charsize = (double) (yyvsp[(4) - (4)].val); }
    break;

  case 375:
/* Line 1792 of yacc.c  */
#line 1567 "gram.y"
    { strcpy(g[curg].timeinfo.format, (char *) (yyvsp[(3) - (3)].str)); }
    break;

  case 376:
/* Line 1792 of yacc.c  */
#line 1571 "gram.y"
    { g[curg].timeline.active = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 377:
/* Line 1792 of yacc.c  */
#line 1572 "gram.y"
    { g[curg].timeline.len = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 378:
/* Line 1792 of yacc.c  */
#line 1573 "gram.y"
    { g[curg].timeline.width = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 379:
/* Line 1792 of yacc.c  */
#line 1574 "gram.y"
    { g[curg].timeline.start = (double) (yyvsp[(3) - (3)].val); }
    break;

  case 380:
/* Line 1792 of yacc.c  */
#line 1575 "gram.y"
    { g[curg].timeline.stop = (double) (yyvsp[(3) - (3)].val); }
    break;

  case 381:
/* Line 1792 of yacc.c  */
#line 1576 "gram.y"
    { g[curg].timeline.step = (double) (yyvsp[(3) - (3)].val); }
    break;

  case 382:
/* Line 1792 of yacc.c  */
#line 1577 "gram.y"
    { g[curg].timeline.p.prec = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 383:
/* Line 1792 of yacc.c  */
#line 1578 "gram.y"
    { g[curg].timeline.units = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 384:
/* Line 1792 of yacc.c  */
#line 1579 "gram.y"
    { g[curg].timeline.c1 = g[curg].timeline.c3 = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 385:
/* Line 1792 of yacc.c  */
#line 1580 "gram.y"
    { g[curg].timeline.c2 = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 386:
/* Line 1792 of yacc.c  */
#line 1581 "gram.y"
    { g[curg].timeline.loctype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 387:
/* Line 1792 of yacc.c  */
#line 1583 "gram.y"
    { 
		g[curg].timeline.x = (double) (yyvsp[(3) - (5)].val);
		g[curg].timeline.y = (double) (yyvsp[(5) - (5)].val);
	}
    break;

  case 388:
/* Line 1792 of yacc.c  */
#line 1587 "gram.y"
    { }
    break;

  case 389:
/* Line 1792 of yacc.c  */
#line 1591 "gram.y"
    { setslice = &g[curg].sbox[(int) (yyvsp[(3) - (3)].val)]; }
    break;

  case 390:
/* Line 1792 of yacc.c  */
#line 1592 "gram.y"
    { setprops = &(setslice->p); }
    break;

  case 392:
/* Line 1792 of yacc.c  */
#line 1594 "gram.y"
    {
		setslice->precx = (int) (yyvsp[(3) - (5)].val);
		setslice->precy = (int) (yyvsp[(5) - (5)].val);
	}
    break;

  case 393:
/* Line 1792 of yacc.c  */
#line 1598 "gram.y"
    { setslice->attach = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 394:
/* Line 1792 of yacc.c  */
#line 1599 "gram.y"
    { setslice->loctype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 395:
/* Line 1792 of yacc.c  */
#line 1600 "gram.y"
    { setslice->display_marker = (int) (yyvsp[(4) - (4)].pset); }
    break;

  case 396:
/* Line 1792 of yacc.c  */
#line 1601 "gram.y"
    { setslice->active = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 397:
/* Line 1792 of yacc.c  */
#line 1602 "gram.y"
    { setslice->display = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 398:
/* Line 1792 of yacc.c  */
#line 1603 "gram.y"
    { setslice->p.color = (yyvsp[(3) - (3)].val); }
    break;

  case 399:
/* Line 1792 of yacc.c  */
#line 1604 "gram.y"
    { setslice->p.linew = (yyvsp[(3) - (3)].val); }
    break;

  case 400:
/* Line 1792 of yacc.c  */
#line 1605 "gram.y"
    { setslice->p.fillcol = (yyvsp[(4) - (4)].val); }
    break;

  case 401:
/* Line 1792 of yacc.c  */
#line 1607 "gram.y"
    { 
		setslice->wx1 = (double) (yyvsp[(3) - (9)].val); 
		setslice->wy1 = (double) (yyvsp[(5) - (9)].val); 
		setslice->wx2 = (double) (yyvsp[(7) - (9)].val); 
		setslice->wy2 = (double) (yyvsp[(9) - (9)].val); 
	}
    break;

  case 402:
/* Line 1792 of yacc.c  */
#line 1614 "gram.y"
    { 
		setslice->vx = (double) (yyvsp[(3) - (5)].val); 
		setslice->vy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 403:
/* Line 1792 of yacc.c  */
#line 1619 "gram.y"
    { 
		setslice->locx = (double) (yyvsp[(3) - (5)].val); 
		setslice->locy = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 404:
/* Line 1792 of yacc.c  */
#line 1624 "gram.y"
    { 
		setslice->x = (double) (yyvsp[(3) - (5)].val); 
		setslice->y = (double) (yyvsp[(5) - (5)].val); 
	}
    break;

  case 405:
/* Line 1792 of yacc.c  */
#line 1631 "gram.y"
    { setprops->color = (yyvsp[(3) - (3)].val); }
    break;

  case 406:
/* Line 1792 of yacc.c  */
#line 1632 "gram.y"
    { setprops->linew = (yyvsp[(3) - (3)].val); }
    break;

  case 407:
/* Line 1792 of yacc.c  */
#line 1633 "gram.y"
    { setprops->lines = (yyvsp[(3) - (3)].val); }
    break;

  case 408:
/* Line 1792 of yacc.c  */
#line 1634 "gram.y"
    { setprops->format = (yyvsp[(3) - (3)].pset); }
    break;

  case 409:
/* Line 1792 of yacc.c  */
#line 1635 "gram.y"
    { setprops->font = (yyvsp[(3) - (3)].val); }
    break;

  case 410:
/* Line 1792 of yacc.c  */
#line 1636 "gram.y"
    { setprops->prec = (yyvsp[(3) - (3)].val); }
    break;

  case 411:
/* Line 1792 of yacc.c  */
#line 1637 "gram.y"
    { setprops->charsize = (yyvsp[(4) - (4)].val); }
    break;

  case 412:
/* Line 1792 of yacc.c  */
#line 1638 "gram.y"
    { setprops->symbol = (yyvsp[(3) - (3)].val); }
    break;

  case 413:
/* Line 1792 of yacc.c  */
#line 1639 "gram.y"
    { setprops->symsize = (yyvsp[(4) - (4)].val); }
    break;

  case 414:
/* Line 1792 of yacc.c  */
#line 1640 "gram.y"
    { setprops->fill = (yyvsp[(3) - (3)].pset); }
    break;

  case 415:
/* Line 1792 of yacc.c  */
#line 1641 "gram.y"
    { setprops->fillusing = (yyvsp[(3) - (3)].pset); }
    break;

  case 416:
/* Line 1792 of yacc.c  */
#line 1642 "gram.y"
    { setprops->fillcol = (yyvsp[(4) - (4)].val); }
    break;

  case 417:
/* Line 1792 of yacc.c  */
#line 1643 "gram.y"
    { setprops->fillpat = (yyvsp[(4) - (4)].val); }
    break;

  case 418:
/* Line 1792 of yacc.c  */
#line 1644 "gram.y"
    { setprops->arrow = (yyvsp[(3) - (3)].val); }
    break;

  case 419:
/* Line 1792 of yacc.c  */
#line 1645 "gram.y"
    { setprops->atype = (yyvsp[(4) - (4)].val); }
    break;

  case 420:
/* Line 1792 of yacc.c  */
#line 1646 "gram.y"
    { setprops->asize = (yyvsp[(4) - (4)].val); }
    break;

  case 421:
/* Line 1792 of yacc.c  */
#line 1650 "gram.y"
    { curg = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 422:
/* Line 1792 of yacc.c  */
#line 1651 "gram.y"
    { curg = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 423:
/* Line 1792 of yacc.c  */
#line 1652 "gram.y"
    { kill_graph((yyvsp[(2) - (2)].pset)); }
    break;

  case 424:
/* Line 1792 of yacc.c  */
#line 1653 "gram.y"
    { kill_graph(maxgraph); }
    break;

  case 425:
/* Line 1792 of yacc.c  */
#line 1655 "gram.y"
    {
	    extern int go_locateflag;
	    go_locateflag = ((yyvsp[(2) - (2)].pset) == ON);
	}
    break;

  case 426:
/* Line 1792 of yacc.c  */
#line 1660 "gram.y"
    {
	    cg = curg = (int) (yyvsp[(2) - (2)].pset);
	    draw_focus(curg);
	    defineworld(g[curg].w.xg1, g[curg].w.yg1, g[curg].w.xg2, g[curg].w.yg2, 
			islogx(curg), islogy(curg));
	    viewport(g[curg].v.xv1, g[curg].v.yv1, g[curg].v.xv2, g[curg].v.yv2);
	    draw_focus(curg);
	    update_all(curg);
	}
    break;

  case 427:
/* Line 1792 of yacc.c  */
#line 1669 "gram.y"
    { draw_focus_flag = (yyvsp[(2) - (2)].pset); }
    break;

  case 428:
/* Line 1792 of yacc.c  */
#line 1670 "gram.y"
    { focus_policy = (yyvsp[(2) - (2)].pset); }
    break;

  case 429:
/* Line 1792 of yacc.c  */
#line 1671 "gram.y"
    { focus_policy = (yyvsp[(2) - (2)].pset); }
    break;

  case 430:
/* Line 1792 of yacc.c  */
#line 1672 "gram.y"
    { focus_policy = (yyvsp[(2) - (2)].pset); }
    break;

  case 431:
/* Line 1792 of yacc.c  */
#line 1673 "gram.y"
    { cursource = (yyvsp[(2) - (2)].pset); }
    break;

  case 432:
/* Line 1792 of yacc.c  */
#line 1674 "gram.y"
    { push_world(); }
    break;

  case 433:
/* Line 1792 of yacc.c  */
#line 1675 "gram.y"
    { pop_world(); }
    break;

  case 434:
/* Line 1792 of yacc.c  */
#line 1676 "gram.y"
    { cycle_world_stack(); }
    break;

  case 435:
/* Line 1792 of yacc.c  */
#line 1677 "gram.y"
    {
	    if ((int) (yyvsp[(2) - (2)].val) > 0)
		show_world_stack((int) (yyvsp[(2) - (2)].val) - 1);
	}
    break;

  case 436:
/* Line 1792 of yacc.c  */
#line 1682 "gram.y"
    {
	    add_world(curg, (yyvsp[(3) - (17)].val), (yyvsp[(5) - (17)].val), (yyvsp[(7) - (17)].val), (yyvsp[(9) - (17)].val), (yyvsp[(11) - (17)].val), (yyvsp[(13) - (17)].val), (yyvsp[(15) - (17)].val), (yyvsp[(17) - (17)].val));
	}
    break;

  case 437:
/* Line 1792 of yacc.c  */
#line 1685 "gram.y"
    { clear_world_stack(); }
    break;

  case 438:
/* Line 1792 of yacc.c  */
#line 1687 "gram.y"
    {
	    g[curg].w.xg1 = (yyvsp[(2) - (8)].val);
	    g[curg].w.yg1 = (yyvsp[(4) - (8)].val);
	    g[curg].w.xg2 = (yyvsp[(6) - (8)].val);
	    g[curg].w.yg2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 439:
/* Line 1792 of yacc.c  */
#line 1693 "gram.y"
    { g[curg].w.xg1 = (yyvsp[(3) - (3)].val); }
    break;

  case 440:
/* Line 1792 of yacc.c  */
#line 1694 "gram.y"
    { g[curg].w.xg2 = (yyvsp[(3) - (3)].val); }
    break;

  case 441:
/* Line 1792 of yacc.c  */
#line 1695 "gram.y"
    { g[curg].w.yg1 = (yyvsp[(3) - (3)].val); }
    break;

  case 442:
/* Line 1792 of yacc.c  */
#line 1696 "gram.y"
    { g[curg].w.yg2 = (yyvsp[(3) - (3)].val); }
    break;

  case 443:
/* Line 1792 of yacc.c  */
#line 1698 "gram.y"
    {
	    g[curg].v.xv1 = (yyvsp[(2) - (8)].val);
	    g[curg].v.yv1 = (yyvsp[(4) - (8)].val);
	    g[curg].v.xv2 = (yyvsp[(6) - (8)].val);
	    g[curg].v.yv2 = (yyvsp[(8) - (8)].val);
	}
    break;

  case 444:
/* Line 1792 of yacc.c  */
#line 1704 "gram.y"
    { g[curg].v.xv1 = (yyvsp[(3) - (3)].val); }
    break;

  case 445:
/* Line 1792 of yacc.c  */
#line 1705 "gram.y"
    { g[curg].v.xv2 = (yyvsp[(3) - (3)].val); }
    break;

  case 446:
/* Line 1792 of yacc.c  */
#line 1706 "gram.y"
    { g[curg].v.yv1 = (yyvsp[(3) - (3)].val); }
    break;

  case 447:
/* Line 1792 of yacc.c  */
#line 1707 "gram.y"
    { g[curg].v.yv2 = (yyvsp[(3) - (3)].val); }
    break;

  case 448:
/* Line 1792 of yacc.c  */
#line 1708 "gram.y"
    {
	    set_plotstr_string(&g[curg].labs.title, (char *) (yyvsp[(2) - (2)].str));
	}
    break;

  case 449:
/* Line 1792 of yacc.c  */
#line 1711 "gram.y"
    {
	    g[curg].labs.title.font = checkon(FONTP, g[curg].labs.title.font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 450:
/* Line 1792 of yacc.c  */
#line 1714 "gram.y"
    { g[curg].labs.title.charsize = (yyvsp[(3) - (3)].val); }
    break;

  case 451:
/* Line 1792 of yacc.c  */
#line 1715 "gram.y"
    {
	    g[curg].labs.title.color = checkon(COLOR, g[curg].labs.title.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 452:
/* Line 1792 of yacc.c  */
#line 1719 "gram.y"
    {
	    g[curg].labs.title.linew = checkon(LINEWIDTH, g[curg].labs.title.linew, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 453:
/* Line 1792 of yacc.c  */
#line 1722 "gram.y"
    {
	    set_plotstr_string(&g[curg].labs.stitle, (char *) (yyvsp[(2) - (2)].str));
	}
    break;

  case 454:
/* Line 1792 of yacc.c  */
#line 1725 "gram.y"
    {
	    g[curg].labs.stitle.font = checkon(FONTP, g[curg].labs.stitle.font, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 455:
/* Line 1792 of yacc.c  */
#line 1728 "gram.y"
    { g[curg].labs.stitle.charsize = (yyvsp[(3) - (3)].val); }
    break;

  case 456:
/* Line 1792 of yacc.c  */
#line 1730 "gram.y"
    {
	    g[curg].labs.stitle.color = checkon(COLOR, g[curg].labs.stitle.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 457:
/* Line 1792 of yacc.c  */
#line 1734 "gram.y"
    {
	    g[curg].labs.stitle.linew = checkon(LINEWIDTH, g[curg].labs.stitle.color, (int) (yyvsp[(3) - (3)].val));
	}
    break;

  case 458:
/* Line 1792 of yacc.c  */
#line 1737 "gram.y"
    { g[curg].f.active = (yyvsp[(2) - (2)].pset); }
    break;

  case 459:
/* Line 1792 of yacc.c  */
#line 1738 "gram.y"
    { g[curg].f.type = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 460:
/* Line 1792 of yacc.c  */
#line 1739 "gram.y"
    { g[curg].f.lines = checkon(LINESTYLE, g[curg].f.lines, (int) (yyvsp[(3) - (3)].val)); }
    break;

  case 461:
/* Line 1792 of yacc.c  */
#line 1740 "gram.y"
    { g[curg].f.linew = checkon(LINEWIDTH, g[curg].f.linew, (int) (yyvsp[(3) - (3)].val)); }
    break;

  case 462:
/* Line 1792 of yacc.c  */
#line 1741 "gram.y"
    { g[curg].f.color = checkon(COLOR, g[curg].f.color, (int) (yyvsp[(3) - (3)].val)); }
    break;

  case 463:
/* Line 1792 of yacc.c  */
#line 1742 "gram.y"
    { g[curg].f.fillbg = (yyvsp[(3) - (3)].pset); }
    break;

  case 464:
/* Line 1792 of yacc.c  */
#line 1743 "gram.y"
    { g[curg].f.bgcolor = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 465:
/* Line 1792 of yacc.c  */
#line 1744 "gram.y"
    { g[(yyvsp[(1) - (2)].pset)].active = (yyvsp[(2) - (2)].pset); }
    break;

  case 466:
/* Line 1792 of yacc.c  */
#line 1745 "gram.y"
    { g[(yyvsp[(1) - (3)].pset)].label = (yyvsp[(3) - (3)].pset); }
    break;

  case 467:
/* Line 1792 of yacc.c  */
#line 1746 "gram.y"
    { g[(yyvsp[(1) - (4)].pset)].auto_type = (yyvsp[(4) - (4)].pset); }
    break;

  case 468:
/* Line 1792 of yacc.c  */
#line 1747 "gram.y"
    { g[(yyvsp[(1) - (4)].pset)].auto_type = (yyvsp[(4) - (4)].pset); }
    break;

  case 469:
/* Line 1792 of yacc.c  */
#line 1748 "gram.y"
    { g[(yyvsp[(1) - (3)].pset)].parmsread = ((yyvsp[(3) - (3)].pset) == FALSEP); }
    break;

  case 470:
/* Line 1792 of yacc.c  */
#line 1749 "gram.y"
    { g[(yyvsp[(1) - (3)].pset)].hidden = ((yyvsp[(3) - (3)].pset) == TRUEP); }
    break;

  case 471:
/* Line 1792 of yacc.c  */
#line 1750 "gram.y"
    { g[(yyvsp[(1) - (3)].pset)].type = (yyvsp[(3) - (3)].pset); }
    break;

  case 472:
/* Line 1792 of yacc.c  */
#line 1751 "gram.y"
    { g[(yyvsp[(1) - (3)].pset)].pointset = ((yyvsp[(3) - (3)].pset) == ON); }
    break;

  case 473:
/* Line 1792 of yacc.c  */
#line 1753 "gram.y"
    {
	    g[(yyvsp[(1) - (5)].pset)].fx = (yyvsp[(4) - (5)].pset);
	    g[(yyvsp[(1) - (5)].pset)].fy = (yyvsp[(5) - (5)].pset);
	}
    break;

  case 474:
/* Line 1792 of yacc.c  */
#line 1758 "gram.y"
    {
	    g[(yyvsp[(1) - (6)].pset)].px = (yyvsp[(4) - (6)].val);
	    g[(yyvsp[(1) - (6)].pset)].py = (yyvsp[(6) - (6)].val);
	}
    break;

  case 475:
/* Line 1792 of yacc.c  */
#line 1763 "gram.y"
    {
	    g[(yyvsp[(1) - (6)].pset)].dsx = (yyvsp[(4) - (6)].val);
	    g[(yyvsp[(1) - (6)].pset)].dsy = (yyvsp[(6) - (6)].val);
	}
    break;

  case 476:
/* Line 1792 of yacc.c  */
#line 1767 "gram.y"
    { g[(yyvsp[(1) - (4)].pset)].pt_type = (int) (yyvsp[(4) - (4)].val); }
    break;

  case 479:
/* Line 1792 of yacc.c  */
#line 1773 "gram.y"
    {}
    break;

  case 480:
/* Line 1792 of yacc.c  */
#line 1774 "gram.y"
    {}
    break;

  case 481:
/* Line 1792 of yacc.c  */
#line 1775 "gram.y"
    {}
    break;

  case 482:
/* Line 1792 of yacc.c  */
#line 1779 "gram.y"
    {}
    break;

  case 483:
/* Line 1792 of yacc.c  */
#line 1780 "gram.y"
    {}
    break;

  case 484:
/* Line 1792 of yacc.c  */
#line 1781 "gram.y"
    {}
    break;

  case 485:
/* Line 1792 of yacc.c  */
#line 1782 "gram.y"
    {}
    break;

  case 486:
/* Line 1792 of yacc.c  */
#line 1783 "gram.y"
    {}
    break;

  case 487:
/* Line 1792 of yacc.c  */
#line 1784 "gram.y"
    {}
    break;

  case 488:
/* Line 1792 of yacc.c  */
#line 1788 "gram.y"
    {}
    break;

  case 489:
/* Line 1792 of yacc.c  */
#line 1789 "gram.y"
    {}
    break;

  case 490:
/* Line 1792 of yacc.c  */
#line 1790 "gram.y"
    {}
    break;

  case 491:
/* Line 1792 of yacc.c  */
#line 1794 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (1)].pset), 0.0); }
    break;

  case 492:
/* Line 1792 of yacc.c  */
#line 1795 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val)); }
    break;

  case 493:
/* Line 1792 of yacc.c  */
#line 1796 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val)); }
    break;

  case 494:
/* Line 1792 of yacc.c  */
#line 1797 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val)); }
    break;

  case 495:
/* Line 1792 of yacc.c  */
#line 1798 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].val)); }
    break;

  case 496:
/* Line 1792 of yacc.c  */
#line 1799 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (3)].pset), (yyvsp[(3) - (3)].val)); }
    break;

  case 497:
/* Line 1792 of yacc.c  */
#line 1800 "gram.y"
    { set_axis_prop(whichgraph, naxis, (yyvsp[(1) - (2)].pset), (yyvsp[(2) - (2)].pset)); }
    break;

  case 498:
/* Line 1792 of yacc.c  */
#line 1805 "gram.y"
    {}
    break;

  case 499:
/* Line 1792 of yacc.c  */
#line 1806 "gram.y"
    {}
    break;

  case 500:
/* Line 1792 of yacc.c  */
#line 1807 "gram.y"
    {}
    break;

  case 501:
/* Line 1792 of yacc.c  */
#line 1808 "gram.y"
    {}
    break;

  case 502:
/* Line 1792 of yacc.c  */
#line 1809 "gram.y"
    { g[curg].t[naxis].active = (yyvsp[(1) - (1)].pset); }
    break;

  case 503:
/* Line 1792 of yacc.c  */
#line 1819 "gram.y"
    {
	    g[curg].t[naxis].t_flag = (yyvsp[(1) - (1)].pset);
	    g[curg].t[naxis].t_mflag = (yyvsp[(1) - (1)].pset);
	}
    break;

  case 504:
/* Line 1792 of yacc.c  */
#line 1823 "gram.y"
    { g[curg].t[naxis].t_flag = (yyvsp[(2) - (2)].pset); }
    break;

  case 505:
/* Line 1792 of yacc.c  */
#line 1824 "gram.y"
    { g[curg].t[naxis].t_mflag = (yyvsp[(2) - (2)].pset); }
    break;

  case 506:
/* Line 1792 of yacc.c  */
#line 1825 "gram.y"
    { g[curg].t[naxis].tmajor = (yyvsp[(2) - (2)].val); }
    break;

  case 507:
/* Line 1792 of yacc.c  */
#line 1826 "gram.y"
    { g[curg].t[naxis].tminor = (yyvsp[(2) - (2)].val); }
    break;

  case 508:
/* Line 1792 of yacc.c  */
#line 1827 "gram.y"
    { g[curg].t[naxis].offsx = (yyvsp[(2) - (2)].val); }
    break;

  case 509:
/* Line 1792 of yacc.c  */
#line 1828 "gram.y"
    { g[curg].t[naxis].offsy = (yyvsp[(2) - (2)].val); }
    break;

  case 510:
/* Line 1792 of yacc.c  */
#line 1829 "gram.y"
    { g[curg].t[naxis].alt = (yyvsp[(2) - (2)].pset); }
    break;

  case 511:
/* Line 1792 of yacc.c  */
#line 1830 "gram.y"
    { g[curg].t[naxis].tmin = (yyvsp[(2) - (2)].val); }
    break;

  case 512:
/* Line 1792 of yacc.c  */
#line 1831 "gram.y"
    { g[curg].t[naxis].tmax = (yyvsp[(2) - (2)].val); }
    break;

  case 513:
/* Line 1792 of yacc.c  */
#line 1832 "gram.y"
    { g[curg].t[naxis].t_num = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 514:
/* Line 1792 of yacc.c  */
#line 1833 "gram.y"
    { g[curg].t[naxis].t_inout = (yyvsp[(1) - (1)].pset); }
    break;

  case 515:
/* Line 1792 of yacc.c  */
#line 1834 "gram.y"
    { g[curg].t[naxis].t_log = (yyvsp[(2) - (2)].pset); }
    break;

  case 516:
/* Line 1792 of yacc.c  */
#line 1835 "gram.y"
    { g[curg].t[naxis].t_size = (yyvsp[(2) - (2)].val); }
    break;

  case 517:
/* Line 1792 of yacc.c  */
#line 1836 "gram.y"
    { g[curg].t[naxis].t_size = (yyvsp[(3) - (3)].val); }
    break;

  case 518:
/* Line 1792 of yacc.c  */
#line 1837 "gram.y"
    { g[curg].t[naxis].t_msize = (yyvsp[(3) - (3)].val); }
    break;

  case 519:
/* Line 1792 of yacc.c  */
#line 1838 "gram.y"
    { g[curg].t[naxis].t_color = g[curg].t[naxis].t_mcolor = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 520:
/* Line 1792 of yacc.c  */
#line 1839 "gram.y"
    { g[curg].t[naxis].t_linew = g[curg].t[naxis].t_mlinew = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 521:
/* Line 1792 of yacc.c  */
#line 1840 "gram.y"
    { g[curg].t[naxis].t_color = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 522:
/* Line 1792 of yacc.c  */
#line 1841 "gram.y"
    { g[curg].t[naxis].t_mcolor = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 523:
/* Line 1792 of yacc.c  */
#line 1842 "gram.y"
    { g[curg].t[naxis].t_linew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 524:
/* Line 1792 of yacc.c  */
#line 1843 "gram.y"
    { g[curg].t[naxis].t_mlinew = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 525:
/* Line 1792 of yacc.c  */
#line 1844 "gram.y"
    { g[curg].t[naxis].t_lines = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 526:
/* Line 1792 of yacc.c  */
#line 1845 "gram.y"
    { g[curg].t[naxis].t_mlines = (int) (yyvsp[(3) - (3)].val); }
    break;

  case 527:
/* Line 1792 of yacc.c  */
#line 1846 "gram.y"
    { g[curg].t[naxis].t_gridflag = (yyvsp[(3) - (3)].pset); }
    break;

  case 528:
/* Line 1792 of yacc.c  */
#line 1847 "gram.y"
    { g[curg].t[naxis].t_mgridflag = (yyvsp[(3) - (3)].pset); }
    break;

  case 529:
/* Line 1792 of yacc.c  */
#line 1848 "gram.y"
    { g[curg].t[naxis].t_op = (yyvsp[(2) - (2)].pset); }
    break;

  case 530:
/* Line 1792 of yacc.c  */
#line 1849 "gram.y"
    { g[curg].t[naxis].t_type = AUTO; }
    break;

  case 531:
/* Line 1792 of yacc.c  */
#line 1850 "gram.y"
    { g[curg].t[naxis].t_type = SPEC; }
    break;

  case 532:
/* Line 1792 of yacc.c  */
#line 1851 "gram.y"
    { g[curg].t[naxis].t_spec = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 533:
/* Line 1792 of yacc.c  */
#line 1852 "gram.y"
    { g[curg].t[naxis].t_specloc[(int) (yyvsp[(1) - (3)].val)] = (yyvsp[(3) - (3)].val); }
    break;

  case 534:
/* Line 1792 of yacc.c  */
#line 1856 "gram.y"
    {}
    break;

  case 535:
/* Line 1792 of yacc.c  */
#line 1857 "gram.y"
    {}
    break;

  case 536:
/* Line 1792 of yacc.c  */
#line 1861 "gram.y"
    { g[curg].t[naxis].tl_flag = (yyvsp[(1) - (1)].pset); }
    break;

  case 537:
/* Line 1792 of yacc.c  */
#line 1862 "gram.y"
    { g[curg].t[naxis].tl_type = AUTO; }
    break;

  case 538:
/* Line 1792 of yacc.c  */
#line 1863 "gram.y"
    { g[curg].t[naxis].tl_type = SPEC; }
    break;

  case 539:
/* Line 1792 of yacc.c  */
#line 1864 "gram.y"
    { g[curg].t[naxis].tl_prec = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 540:
/* Line 1792 of yacc.c  */
#line 1865 "gram.y"
    { g[curg].t[naxis].tl_format = (yyvsp[(2) - (2)].pset); }
    break;

  case 541:
/* Line 1792 of yacc.c  */
#line 1866 "gram.y"
    { g[curg].t[naxis].tl_format = (yyvsp[(2) - (2)].val); }
    break;

  case 542:
/* Line 1792 of yacc.c  */
#line 1867 "gram.y"
    { strcpy(g[curg].t[naxis].tl_appstr, (char *) (yyvsp[(2) - (2)].str)); }
    break;

  case 543:
/* Line 1792 of yacc.c  */
#line 1868 "gram.y"
    { strcpy(g[curg].t[naxis].tl_prestr, (char *) (yyvsp[(2) - (2)].str)); }
    break;

  case 544:
/* Line 1792 of yacc.c  */
#line 1869 "gram.y"
    { g[curg].t[naxis].tl_layout = HORIZONTAL; }
    break;

  case 545:
/* Line 1792 of yacc.c  */
#line 1870 "gram.y"
    { g[curg].t[naxis].tl_layout = VERTICAL; }
    break;

  case 546:
/* Line 1792 of yacc.c  */
#line 1871 "gram.y"
    { g[curg].t[naxis].tl_layout = SPEC; }
    break;

  case 547:
/* Line 1792 of yacc.c  */
#line 1872 "gram.y"
    { g[curg].t[naxis].tl_angle = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 548:
/* Line 1792 of yacc.c  */
#line 1873 "gram.y"
    { g[curg].t[naxis].tl_just = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 549:
/* Line 1792 of yacc.c  */
#line 1874 "gram.y"
    { g[curg].t[naxis].tl_skip = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 550:
/* Line 1792 of yacc.c  */
#line 1875 "gram.y"
    { g[curg].t[naxis].tl_staggered = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 551:
/* Line 1792 of yacc.c  */
#line 1876 "gram.y"
    { g[curg].t[naxis].tl_op = (yyvsp[(2) - (2)].pset); }
    break;

  case 552:
/* Line 1792 of yacc.c  */
#line 1877 "gram.y"
    { g[curg].t[naxis].tl_sign = (yyvsp[(2) - (2)].pset); }
    break;

  case 553:
/* Line 1792 of yacc.c  */
#line 1878 "gram.y"
    { g[curg].t[naxis].tl_start = (yyvsp[(2) - (2)].val); }
    break;

  case 554:
/* Line 1792 of yacc.c  */
#line 1879 "gram.y"
    { g[curg].t[naxis].tl_stop = (yyvsp[(2) - (2)].val); }
    break;

  case 555:
/* Line 1792 of yacc.c  */
#line 1880 "gram.y"
    { g[curg].t[naxis].tl_starttype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 556:
/* Line 1792 of yacc.c  */
#line 1881 "gram.y"
    { g[curg].t[naxis].tl_starttype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 557:
/* Line 1792 of yacc.c  */
#line 1882 "gram.y"
    { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 558:
/* Line 1792 of yacc.c  */
#line 1883 "gram.y"
    { g[curg].t[naxis].tl_stoptype = (int) (yyvsp[(3) - (3)].pset); }
    break;

  case 559:
/* Line 1792 of yacc.c  */
#line 1884 "gram.y"
    { g[curg].t[naxis].tl_vgap = (yyvsp[(2) - (2)].val); }
    break;

  case 560:
/* Line 1792 of yacc.c  */
#line 1885 "gram.y"
    { g[curg].t[naxis].tl_hgap = (yyvsp[(2) - (2)].val); }
    break;

  case 561:
/* Line 1792 of yacc.c  */
#line 1886 "gram.y"
    { g[curg].t[naxis].tl_charsize = (yyvsp[(3) - (3)].val); }
    break;

  case 562:
/* Line 1792 of yacc.c  */
#line 1887 "gram.y"
    { g[curg].t[naxis].tl_font = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 563:
/* Line 1792 of yacc.c  */
#line 1888 "gram.y"
    { g[curg].t[naxis].tl_color = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 564:
/* Line 1792 of yacc.c  */
#line 1889 "gram.y"
    { g[curg].t[naxis].tl_linew = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 565:
/* Line 1792 of yacc.c  */
#line 1890 "gram.y"
    { set_plotstr_string(&g[curg].t[naxis].t_speclab[(int) (yyvsp[(1) - (3)].val)], (char *) (yyvsp[(3) - (3)].str)); }
    break;

  case 566:
/* Line 1792 of yacc.c  */
#line 1894 "gram.y"
    { set_plotstr_string(&g[curg].t[naxis].label, (char *) (yyvsp[(1) - (1)].str)); }
    break;

  case 567:
/* Line 1792 of yacc.c  */
#line 1895 "gram.y"
    { g[curg].t[naxis].label_layout = PERP; }
    break;

  case 568:
/* Line 1792 of yacc.c  */
#line 1896 "gram.y"
    { g[curg].t[naxis].label_layout = PARA; }
    break;

  case 569:
/* Line 1792 of yacc.c  */
#line 1897 "gram.y"
    { g[curg].t[naxis].label_place = (yyvsp[(2) - (2)].pset); }
    break;

  case 570:
/* Line 1792 of yacc.c  */
#line 1898 "gram.y"
    { g[curg].t[naxis].label_place = (yyvsp[(2) - (2)].pset); }
    break;

  case 571:
/* Line 1792 of yacc.c  */
#line 1899 "gram.y"
    {
	    g[curg].t[naxis].label.x = (yyvsp[(2) - (4)].val);
	    g[curg].t[naxis].label.y = (yyvsp[(4) - (4)].val);
	}
    break;

  case 572:
/* Line 1792 of yacc.c  */
#line 1903 "gram.y"
    { g[curg].t[naxis].label.just = (int) (yyvsp[(2) - (2)].pset); }
    break;

  case 573:
/* Line 1792 of yacc.c  */
#line 1904 "gram.y"
    { g[curg].t[naxis].label.charsize = (yyvsp[(3) - (3)].val); }
    break;

  case 574:
/* Line 1792 of yacc.c  */
#line 1905 "gram.y"
    { g[curg].t[naxis].label.font = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 575:
/* Line 1792 of yacc.c  */
#line 1906 "gram.y"
    { g[curg].t[naxis].label.color = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 576:
/* Line 1792 of yacc.c  */
#line 1907 "gram.y"
    { g[curg].t[naxis].label.linew = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 577:
/* Line 1792 of yacc.c  */
#line 1911 "gram.y"
    { g[curg].t[naxis].t_drawbar = (yyvsp[(1) - (1)].pset); }
    break;

  case 578:
/* Line 1792 of yacc.c  */
#line 1912 "gram.y"
    { g[curg].t[naxis].t_drawbarcolor = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 579:
/* Line 1792 of yacc.c  */
#line 1913 "gram.y"
    { g[curg].t[naxis].t_drawbarlines = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 580:
/* Line 1792 of yacc.c  */
#line 1914 "gram.y"
    { g[curg].t[naxis].t_drawbarlinew = (int) (yyvsp[(2) - (2)].val); }
    break;

  case 581:
/* Line 1792 of yacc.c  */
#line 1926 "gram.y"
    { (yyval.pset) = DISK; }
    break;

  case 582:
/* Line 1792 of yacc.c  */
#line 1927 "gram.y"
    { (yyval.pset) = PIPE; }
    break;

  case 583:
/* Line 1792 of yacc.c  */
#line 1931 "gram.y"
    { (yyval.pset) = RIGHT; }
    break;

  case 584:
/* Line 1792 of yacc.c  */
#line 1932 "gram.y"
    { (yyval.pset) = LEFT; }
    break;

  case 585:
/* Line 1792 of yacc.c  */
#line 1933 "gram.y"
    { (yyval.pset) = CENTER; }
    break;

  case 586:
/* Line 1792 of yacc.c  */
#line 1942 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 587:
/* Line 1792 of yacc.c  */
#line 1943 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 588:
/* Line 1792 of yacc.c  */
#line 1944 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 589:
/* Line 1792 of yacc.c  */
#line 1945 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 590:
/* Line 1792 of yacc.c  */
#line 1946 "gram.y"
    { (yyval.pset) = XYFIXED; }
    break;

  case 591:
/* Line 1792 of yacc.c  */
#line 1950 "gram.y"
    { (yyval.pset) = IN; }
    break;

  case 592:
/* Line 1792 of yacc.c  */
#line 1951 "gram.y"
    { (yyval.pset) = OUT; }
    break;

  case 593:
/* Line 1792 of yacc.c  */
#line 1952 "gram.y"
    { (yyval.pset) = BOTH; }
    break;

  case 594:
/* Line 1792 of yacc.c  */
#line 1956 "gram.y"
    { (yyval.pset) = NORMAL; }
    break;

  case 595:
/* Line 1792 of yacc.c  */
#line 1957 "gram.y"
    { (yyval.pset) = ABSOLUTE; }
    break;

  case 596:
/* Line 1792 of yacc.c  */
#line 1958 "gram.y"
    { (yyval.pset) = NEGATE; }
    break;

  case 597:
/* Line 1792 of yacc.c  */
#line 1971 "gram.y"
    { (yyval.pset) = DECIMAL; }
    break;

  case 598:
/* Line 1792 of yacc.c  */
#line 1972 "gram.y"
    { (yyval.pset) = EXPONENTIAL; }
    break;

  case 599:
/* Line 1792 of yacc.c  */
#line 1973 "gram.y"
    { (yyval.pset) = POWER; }
    break;

  case 600:
/* Line 1792 of yacc.c  */
#line 1974 "gram.y"
    { (yyval.pset) = GENERAL; }
    break;

  case 601:
/* Line 1792 of yacc.c  */
#line 1975 "gram.y"
    { (yyval.pset) = DDMMYY; }
    break;

  case 602:
/* Line 1792 of yacc.c  */
#line 1976 "gram.y"
    { (yyval.pset) = MMDDYY; }
    break;

  case 603:
/* Line 1792 of yacc.c  */
#line 1977 "gram.y"
    { (yyval.pset) = MMYY; }
    break;

  case 604:
/* Line 1792 of yacc.c  */
#line 1978 "gram.y"
    { (yyval.pset) = MMDD; }
    break;

  case 605:
/* Line 1792 of yacc.c  */
#line 1979 "gram.y"
    { (yyval.pset) = MONTHDAY; }
    break;

  case 606:
/* Line 1792 of yacc.c  */
#line 1980 "gram.y"
    { (yyval.pset) = DAYMONTH; }
    break;

  case 607:
/* Line 1792 of yacc.c  */
#line 1981 "gram.y"
    { (yyval.pset) = DDMONTHSYYHHMMSS; }
    break;

  case 608:
/* Line 1792 of yacc.c  */
#line 1982 "gram.y"
    { (yyval.pset) = MONTHS; }
    break;

  case 609:
/* Line 1792 of yacc.c  */
#line 1983 "gram.y"
    { (yyval.pset) = MONTHL; }
    break;

  case 610:
/* Line 1792 of yacc.c  */
#line 1984 "gram.y"
    { (yyval.pset) = DAYOFWEEKS; }
    break;

  case 611:
/* Line 1792 of yacc.c  */
#line 1985 "gram.y"
    { (yyval.pset) = DAYOFWEEKL; }
    break;

  case 612:
/* Line 1792 of yacc.c  */
#line 1986 "gram.y"
    { (yyval.pset) = DAYOFYEAR; }
    break;

  case 613:
/* Line 1792 of yacc.c  */
#line 1987 "gram.y"
    { (yyval.pset) = HMS; }
    break;

  case 614:
/* Line 1792 of yacc.c  */
#line 1988 "gram.y"
    { (yyval.pset) = MMDDHMS; }
    break;

  case 615:
/* Line 1792 of yacc.c  */
#line 1989 "gram.y"
    { (yyval.pset) = MMDDYYHMS; }
    break;

  case 616:
/* Line 1792 of yacc.c  */
#line 1990 "gram.y"
    { (yyval.pset) = DEGREESLON; }
    break;

  case 617:
/* Line 1792 of yacc.c  */
#line 1991 "gram.y"
    { (yyval.pset) = DEGREESMMLON; }
    break;

  case 618:
/* Line 1792 of yacc.c  */
#line 1992 "gram.y"
    { (yyval.pset) = DEGREESMMSSLON; }
    break;

  case 619:
/* Line 1792 of yacc.c  */
#line 1993 "gram.y"
    { (yyval.pset) = MMSSLON; }
    break;

  case 620:
/* Line 1792 of yacc.c  */
#line 1994 "gram.y"
    { (yyval.pset) = DEGREESLAT; }
    break;

  case 621:
/* Line 1792 of yacc.c  */
#line 1995 "gram.y"
    { (yyval.pset) = DEGREESMMLAT; }
    break;

  case 622:
/* Line 1792 of yacc.c  */
#line 1996 "gram.y"
    { (yyval.pset) = DEGREESMMSSLAT; }
    break;

  case 623:
/* Line 1792 of yacc.c  */
#line 1997 "gram.y"
    { (yyval.pset) = MMSSLAT; }
    break;

  case 624:
/* Line 1792 of yacc.c  */
#line 2000 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 625:
/* Line 1792 of yacc.c  */
#line 2001 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 626:
/* Line 1792 of yacc.c  */
#line 2005 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 627:
/* Line 1792 of yacc.c  */
#line 2006 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 628:
/* Line 1792 of yacc.c  */
#line 2007 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 629:
/* Line 1792 of yacc.c  */
#line 2008 "gram.y"
    { (yyval.pset) = (yyvsp[(1) - (1)].pset); }
    break;

  case 630:
/* Line 1792 of yacc.c  */
#line 2012 "gram.y"
    { (yyval.pset) = ON; }
    break;

  case 631:
/* Line 1792 of yacc.c  */
#line 2013 "gram.y"
    { (yyval.pset) = OFF; }
    break;

  case 632:
/* Line 1792 of yacc.c  */
#line 2016 "gram.y"
    { (yyval.pset) = WORLD; }
    break;

  case 633:
/* Line 1792 of yacc.c  */
#line 2017 "gram.y"
    { (yyval.pset) = VIEW; }
    break;

  case 634:
/* Line 1792 of yacc.c  */
#line 2028 "gram.y"
    { (yyval.pset) = TRUEP; }
    break;

  case 635:
/* Line 1792 of yacc.c  */
#line 2029 "gram.y"
    { (yyval.pset) = FALSEP; }
    break;

  case 636:
/* Line 1792 of yacc.c  */
#line 2032 "gram.y"
    { (yyval.pset) = PATTERN; }
    break;

  case 637:
/* Line 1792 of yacc.c  */
#line 2033 "gram.y"
    { (yyval.pset) = COLOR; }
    break;

  case 638:
/* Line 1792 of yacc.c  */
#line 2034 "gram.y"
    { (yyval.pset) = NONE; }
    break;

  case 639:
/* Line 1792 of yacc.c  */
#line 2037 "gram.y"
    { (yyval.pset) = ABOVE; }
    break;

  case 640:
/* Line 1792 of yacc.c  */
#line 2038 "gram.y"
    { (yyval.pset) = BELOW; }
    break;

  case 641:
/* Line 1792 of yacc.c  */
#line 2039 "gram.y"
    { (yyval.pset) = LEFT; }
    break;

  case 642:
/* Line 1792 of yacc.c  */
#line 2040 "gram.y"
    { (yyval.pset) = RIGHT; }
    break;

  case 643:
/* Line 1792 of yacc.c  */
#line 2041 "gram.y"
    { (yyval.pset) = TOP; }
    break;

  case 644:
/* Line 1792 of yacc.c  */
#line 2042 "gram.y"
    { (yyval.pset) = BOTTOM ; }
    break;

  case 645:
/* Line 1792 of yacc.c  */
#line 2043 "gram.y"
    { (yyval.pset) = BOTH ; }
    break;

  case 646:
/* Line 1792 of yacc.c  */
#line 2047 "gram.y"
    { (yyval.pset) = MM; }
    break;

  case 647:
/* Line 1792 of yacc.c  */
#line 2048 "gram.y"
    { (yyval.pset) = CM; }
    break;

  case 648:
/* Line 1792 of yacc.c  */
#line 2049 "gram.y"
    { (yyval.pset) = M; }
    break;

  case 649:
/* Line 1792 of yacc.c  */
#line 2050 "gram.y"
    { (yyval.pset) = KM; }
    break;

  case 650:
/* Line 1792 of yacc.c  */
#line 2065 "gram.y"
    {
	    int itmp = (int) (yyvsp[(3) - (6)].val) - 1;
	    if (itmp >= ls) {
		yyerror("subscript out of range");
		return 1;
	    } else {
		(yyvsp[(1) - (6)].ptr)[itmp] = (yyvsp[(6) - (6)].val);
		result = (yyvsp[(6) - (6)].val);
	    }
	}
    break;

  case 651:
/* Line 1792 of yacc.c  */
#line 2079 "gram.y"
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[(1) - (3)].ptr)[i] = (yyvsp[(3) - (3)].ptr)[i];
	    }
	    result = (yyvsp[(3) - (3)].ptr)[0];
	}
    break;

  case 652:
/* Line 1792 of yacc.c  */
#line 2087 "gram.y"
    {
	    int i;
	    for (i = 0; i < lxy; i++) {
		(yyvsp[(1) - (3)].ptr)[i] = (yyvsp[(3) - (3)].val);
	    }
	    result = (yyvsp[(3) - (3)].val);
	}
    break;

  case 653:
/* Line 1792 of yacc.c  */
#line 2098 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].ptr)[i];
	    }
	}
    break;

  case 654:
/* Line 1792 of yacc.c  */
#line 2107 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].val);
	    }
	}
    break;

  case 655:
/* Line 1792 of yacc.c  */
#line 2116 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 656:
/* Line 1792 of yacc.c  */
#line 2125 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] + (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 657:
/* Line 1792 of yacc.c  */
#line 2134 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 658:
/* Line 1792 of yacc.c  */
#line 2143 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] + (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 659:
/* Line 1792 of yacc.c  */
#line 2152 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 660:
/* Line 1792 of yacc.c  */
#line 2161 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] - (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 661:
/* Line 1792 of yacc.c  */
#line 2170 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 662:
/* Line 1792 of yacc.c  */
#line 2179 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] - (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 663:
/* Line 1792 of yacc.c  */
#line 2188 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 664:
/* Line 1792 of yacc.c  */
#line 2197 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] * (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 665:
/* Line 1792 of yacc.c  */
#line 2206 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 666:
/* Line 1792 of yacc.c  */
#line 2215 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] * (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 667:
/* Line 1792 of yacc.c  */
#line 2224 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[(3) - (3)].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 668:
/* Line 1792 of yacc.c  */
#line 2237 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[(3) - (3)].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] / (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 669:
/* Line 1792 of yacc.c  */
#line 2252 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		if ((yyvsp[(3) - (3)].ptr)[i] == 0.0) {
		    yyerror("Divide by Zero");
		    return 1;
		}
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 670:
/* Line 1792 of yacc.c  */
#line 2267 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    if ((yyvsp[(3) - (3)].val) == 0.0) {
		yyerror("Divide by Zero");
		return 1;
	    }
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] / (yyvsp[(3) - (3)].val);
	    }
	}
    break;

  case 671:
/* Line 1792 of yacc.c  */
#line 2280 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	    }
	}
    break;

  case 672:
/* Line 1792 of yacc.c  */
#line 2289 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].ptr)[i]);
	    }
	}
    break;

  case 673:
/* Line 1792 of yacc.c  */
#line 2298 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].ptr)[i], (yyvsp[(3) - (3)].val));
	    }
	}
    break;

  case 674:
/* Line 1792 of yacc.c  */
#line 2307 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = pow((yyvsp[(1) - (3)].ptr)[i], (yyvsp[(3) - (3)].ptr)[i]);
	    }
	}
    break;

  case 675:
/* Line 1792 of yacc.c  */
#line 2316 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[(3) - (4)].val));
	    }
	}
    break;

  case 676:
/* Line 1792 of yacc.c  */
#line 2325 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fabs((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 677:
/* Line 1792 of yacc.c  */
#line 2334 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = acos((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 678:
/* Line 1792 of yacc.c  */
#line 2343 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = asin((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 679:
/* Line 1792 of yacc.c  */
#line 2352 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 680:
/* Line 1792 of yacc.c  */
#line 2361 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = atan2((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 681:
/* Line 1792 of yacc.c  */
#line 2370 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = ceil((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 682:
/* Line 1792 of yacc.c  */
#line 2379 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = cos((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 683:
/* Line 1792 of yacc.c  */
#line 2388 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] *= M_PI / 180.0;
	    }
	}
    break;

  case 684:
/* Line 1792 of yacc.c  */
#line 2397 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = xx[i];
	    }
	}
    break;

  case 685:
/* Line 1792 of yacc.c  */
#line 2406 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = yy[i];
	    }
	}
    break;

  case 686:
/* Line 1792 of yacc.c  */
#line 2415 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erf((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 687:
/* Line 1792 of yacc.c  */
#line 2424 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = erfc((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 688:
/* Line 1792 of yacc.c  */
#line 2433 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = exp((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 689:
/* Line 1792 of yacc.c  */
#line 2442 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = floor((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 690:
/* Line 1792 of yacc.c  */
#line 2451 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 691:
/* Line 1792 of yacc.c  */
#line 2460 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 692:
/* Line 1792 of yacc.c  */
#line 2469 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 693:
/* Line 1792 of yacc.c  */
#line 2478 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	    }
	}
    break;

  case 694:
/* Line 1792 of yacc.c  */
#line 2487 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = i + 1;
	    }
	}
    break;

  case 695:
/* Line 1792 of yacc.c  */
#line 2496 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (1)].func);
	    }
	}
    break;

  case 696:
/* Line 1792 of yacc.c  */
#line 2505 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (int) (yyvsp[(3) - (4)].ptr)[i];
	    }
	}
    break;

  case 697:
/* Line 1792 of yacc.c  */
#line 2514 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lrand48() % (long) ((yyvsp[(3) - (4)].val));
	    }
	}
    break;

  case 698:
/* Line 1792 of yacc.c  */
#line 2523 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = lgamma((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 699:
/* Line 1792 of yacc.c  */
#line 2532 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 700:
/* Line 1792 of yacc.c  */
#line 2541 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = log10((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 701:
/* Line 1792 of yacc.c  */
#line 2550 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = 1.0 / (1.0 + exp(-((yyvsp[(3) - (8)].ptr)[i] - (yyvsp[(5) - (8)].val))/ (yyvsp[(7) - (8)].val)));
	    }
	}
    break;

  case 702:
/* Line 1792 of yacc.c  */
#line 2559 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (6)].ptr)[i] >= (yyvsp[(5) - (6)].ptr)[i] ? (yyvsp[(3) - (6)].ptr)[i] : (yyvsp[(5) - (6)].ptr)[i];
	    }
	}
    break;

  case 703:
/* Line 1792 of yacc.c  */
#line 2568 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (6)].ptr)[i] <= (yyvsp[(5) - (6)].ptr)[i] ? (yyvsp[(3) - (6)].ptr)[i] : (yyvsp[(5) - (6)].ptr)[i];
	    }
	}
    break;

  case 704:
/* Line 1792 of yacc.c  */
#line 2577 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fmod((yyvsp[(3) - (6)].ptr)[i], (yyvsp[(5) - (6)].ptr)[i]);
	    }
	}
    break;

  case 705:
/* Line 1792 of yacc.c  */
#line 2586 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = fx((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 706:
/* Line 1792 of yacc.c  */
#line 2595 "gram.y"
    {
	    int i;
	    double tmp;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = normp((yyvsp[(3) - (4)].ptr)[i], &tmp);
	    }
	}
    break;

  case 707:
/* Line 1792 of yacc.c  */
#line 2605 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI;
	    }
	}
    break;

  case 708:
/* Line 1792 of yacc.c  */
#line 2614 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = M_PI / 180.0;
	    }
	}
    break;

  case 709:
/* Line 1792 of yacc.c  */
#line 2623 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (double) drand48();
	    }
	}
    break;

  case 710:
/* Line 1792 of yacc.c  */
#line 2632 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sin((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 711:
/* Line 1792 of yacc.c  */
#line 2641 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(3) - (4)].ptr)[i] * (yyvsp[(3) - (4)].ptr)[i];
	    }
	}
    break;

  case 712:
/* Line 1792 of yacc.c  */
#line 2650 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = sqrt((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 713:
/* Line 1792 of yacc.c  */
#line 2659 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = tan((yyvsp[(3) - (4)].ptr)[i]);
	    }
	}
    break;

  case 714:
/* Line 1792 of yacc.c  */
#line 2668 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] > (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 715:
/* Line 1792 of yacc.c  */
#line 2677 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] < (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 716:
/* Line 1792 of yacc.c  */
#line 2686 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] <= (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 717:
/* Line 1792 of yacc.c  */
#line 2695 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] >= (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 718:
/* Line 1792 of yacc.c  */
#line 2704 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] == (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 719:
/* Line 1792 of yacc.c  */
#line 2713 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] != (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 720:
/* Line 1792 of yacc.c  */
#line 2722 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] && (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 721:
/* Line 1792 of yacc.c  */
#line 2731 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(1) - (3)].ptr)[i] || (yyvsp[(3) - (3)].ptr)[i];
	    }
	}
    break;

  case 722:
/* Line 1792 of yacc.c  */
#line 2740 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = !((yyvsp[(2) - (2)].ptr)[i]);
	    }
	}
    break;

  case 723:
/* Line 1792 of yacc.c  */
#line 2749 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = (yyvsp[(2) - (3)].ptr)[i];
	    }
	}
    break;

  case 724:
/* Line 1792 of yacc.c  */
#line 2757 "gram.y"
    {
	    int i;
	    (yyval.ptr) = (double *) calloc(lxy, sizeof(double));
	    freelist[fcnt++] = (yyval.ptr);
	    for (i = 0; i < lxy; i++) {
		(yyval.ptr)[i] = -(yyvsp[(2) - (2)].ptr)[i];
	    }
	}
    break;

  case 726:
/* Line 1792 of yacc.c  */
#line 2768 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (1)].val);
	}
    break;

  case 727:
/* Line 1792 of yacc.c  */
#line 2771 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (4)].ptr)[(int) (yyvsp[(3) - (4)].val)];
	}
    break;

  case 728:
/* Line 1792 of yacc.c  */
#line 2774 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) + (yyvsp[(3) - (3)].val);
	}
    break;

  case 729:
/* Line 1792 of yacc.c  */
#line 2777 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) - (yyvsp[(3) - (3)].val);
	}
    break;

  case 730:
/* Line 1792 of yacc.c  */
#line 2780 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) * (yyvsp[(3) - (3)].val);
	}
    break;

  case 731:
/* Line 1792 of yacc.c  */
#line 2784 "gram.y"
    {
	    if ((yyvsp[(3) - (3)].val) != 0.0) {
		(yyval.val) = (yyvsp[(1) - (3)].val) / (yyvsp[(3) - (3)].val);
	    } else {
		yyerror("Divide by Zero");
		return 1;
	    }
	}
    break;

  case 732:
/* Line 1792 of yacc.c  */
#line 2792 "gram.y"
    {
	    (yyval.val) = fmod((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	}
    break;

  case 733:
/* Line 1792 of yacc.c  */
#line 2795 "gram.y"
    {
	    (yyval.val) = pow((yyvsp[(1) - (3)].val), (yyvsp[(3) - (3)].val));
	}
    break;

  case 734:
/* Line 1792 of yacc.c  */
#line 2798 "gram.y"
    {
	    (yyval.val) = fabs((yyvsp[(3) - (4)].val));
	}
    break;

  case 735:
/* Line 1792 of yacc.c  */
#line 2801 "gram.y"
    {
	    (yyval.val) = acos((yyvsp[(3) - (4)].val));
	}
    break;

  case 736:
/* Line 1792 of yacc.c  */
#line 2804 "gram.y"
    {
	    (yyval.val) = asin((yyvsp[(3) - (4)].val));
	}
    break;

  case 737:
/* Line 1792 of yacc.c  */
#line 2807 "gram.y"
    {
	    (yyval.val) = atan((yyvsp[(3) - (4)].val));
	}
    break;

  case 738:
/* Line 1792 of yacc.c  */
#line 2810 "gram.y"
    {
	    (yyval.val) = atan2((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 739:
/* Line 1792 of yacc.c  */
#line 2813 "gram.y"
    {
	    (yyval.val) = ceil((yyvsp[(3) - (4)].val));
	}
    break;

  case 740:
/* Line 1792 of yacc.c  */
#line 2816 "gram.y"
    {
	    (yyval.val) = cos((yyvsp[(3) - (4)].val));
	}
    break;

  case 741:
/* Line 1792 of yacc.c  */
#line 2819 "gram.y"
    {
	    (yyval.val) = 180.0 / M_PI;
	}
    break;

  case 742:
/* Line 1792 of yacc.c  */
#line 2822 "gram.y"
    {
	    (yyval.val) = *xx;
	}
    break;

  case 743:
/* Line 1792 of yacc.c  */
#line 2825 "gram.y"
    {
	    (yyval.val) = *yy;
	}
    break;

  case 744:
/* Line 1792 of yacc.c  */
#line 2828 "gram.y"
    {
	    (yyval.val) = erf((yyvsp[(3) - (4)].val));
	}
    break;

  case 745:
/* Line 1792 of yacc.c  */
#line 2831 "gram.y"
    {
	    (yyval.val) = erfc((yyvsp[(3) - (4)].val));
	}
    break;

  case 746:
/* Line 1792 of yacc.c  */
#line 2834 "gram.y"
    {
	    (yyval.val) = exp((yyvsp[(3) - (4)].val));
	}
    break;

  case 747:
/* Line 1792 of yacc.c  */
#line 2837 "gram.y"
    {
	    (yyval.val) = floor((yyvsp[(3) - (4)].val));
	}
    break;

  case 748:
/* Line 1792 of yacc.c  */
#line 2840 "gram.y"
    {
	    (yyval.val) = hypot((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 749:
/* Line 1792 of yacc.c  */
#line 2843 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.xv1;
	}
    break;

  case 750:
/* Line 1792 of yacc.c  */
#line 2846 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.xv2;
	}
    break;

  case 751:
/* Line 1792 of yacc.c  */
#line 2849 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.yv1;
	}
    break;

  case 752:
/* Line 1792 of yacc.c  */
#line 2852 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].v.yv2;
	}
    break;

  case 753:
/* Line 1792 of yacc.c  */
#line 2855 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.xg1;
	}
    break;

  case 754:
/* Line 1792 of yacc.c  */
#line 2858 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.xg2;
	}
    break;

  case 755:
/* Line 1792 of yacc.c  */
#line 2861 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.yg1;
	}
    break;

  case 756:
/* Line 1792 of yacc.c  */
#line 2864 "gram.y"
    {
	    (yyval.val) = g[(yyvsp[(1) - (3)].pset)].w.yg2;
	}
    break;

  case 757:
/* Line 1792 of yacc.c  */
#line 2867 "gram.y"
    {
	    (yyval.val) = g[curg].v.xv1;
	}
    break;

  case 758:
/* Line 1792 of yacc.c  */
#line 2870 "gram.y"
    {
	    (yyval.val) = g[curg].v.xv2;
	}
    break;

  case 759:
/* Line 1792 of yacc.c  */
#line 2873 "gram.y"
    {
	    (yyval.val) = g[curg].v.yv1;
	}
    break;

  case 760:
/* Line 1792 of yacc.c  */
#line 2876 "gram.y"
    {
	    (yyval.val) = g[curg].v.yv2;
	}
    break;

  case 761:
/* Line 1792 of yacc.c  */
#line 2879 "gram.y"
    {
	    (yyval.val) = g[curg].w.xg1;
	}
    break;

  case 762:
/* Line 1792 of yacc.c  */
#line 2882 "gram.y"
    {
	    (yyval.val) = g[curg].w.xg2;
	}
    break;

  case 763:
/* Line 1792 of yacc.c  */
#line 2885 "gram.y"
    {
	    (yyval.val) = g[curg].w.yg1;
	}
    break;

  case 764:
/* Line 1792 of yacc.c  */
#line 2888 "gram.y"
    {
	    (yyval.val) = g[curg].w.yg2;
	}
    break;

  case 765:
/* Line 1792 of yacc.c  */
#line 2891 "gram.y"
    {
	    (yyval.val) = setindex;
	}
    break;

  case 766:
/* Line 1792 of yacc.c  */
#line 2894 "gram.y"
    {
	    (yyval.val) = setsetno;
	}
    break;

  case 767:
/* Line 1792 of yacc.c  */
#line 2897 "gram.y"
    {
	    (yyval.val) = (long) (yyvsp[(3) - (4)].val);
	}
    break;

  case 768:
/* Line 1792 of yacc.c  */
#line 2900 "gram.y"
    {
	    (yyval.val) = lrand48() % (long) ((yyvsp[(3) - (4)].val));
	}
    break;

  case 769:
/* Line 1792 of yacc.c  */
#line 2903 "gram.y"
    {
	    (yyval.val) = lgamma((yyvsp[(3) - (4)].val));
	}
    break;

  case 770:
/* Line 1792 of yacc.c  */
#line 2906 "gram.y"
    {
	    (yyval.val) = log((yyvsp[(3) - (4)].val));
	}
    break;

  case 771:
/* Line 1792 of yacc.c  */
#line 2909 "gram.y"
    {
	    (yyval.val) = log10((yyvsp[(3) - (4)].val));
	}
    break;

  case 772:
/* Line 1792 of yacc.c  */
#line 2913 "gram.y"
    {
	    (yyval.val) = 1.0 / (1.0 + exp(-((yyvsp[(3) - (8)].val) - (yyvsp[(5) - (8)].val))/ (yyvsp[(7) - (8)].val)));
	}
    break;

  case 773:
/* Line 1792 of yacc.c  */
#line 2916 "gram.y"
    {
	    (yyval.val) = (yyvsp[(3) - (6)].val) >= (yyvsp[(5) - (6)].val) ? (yyvsp[(3) - (6)].val) : (yyvsp[(5) - (6)].val);
	}
    break;

  case 774:
/* Line 1792 of yacc.c  */
#line 2919 "gram.y"
    {
	    (yyval.val) = (yyvsp[(3) - (6)].val) <= (yyvsp[(5) - (6)].val) ? (yyvsp[(3) - (6)].val) : (yyvsp[(5) - (6)].val);
	}
    break;

  case 775:
/* Line 1792 of yacc.c  */
#line 2922 "gram.y"
    {
	    (yyval.val) = fmod((yyvsp[(3) - (6)].val), (yyvsp[(5) - (6)].val));
	}
    break;

  case 776:
/* Line 1792 of yacc.c  */
#line 2925 "gram.y"
    {
	    (yyval.val) = fx((yyvsp[(3) - (4)].val));
	}
    break;

  case 777:
/* Line 1792 of yacc.c  */
#line 2928 "gram.y"
    {
	    double tmp;
	    (yyval.val) = normp((yyvsp[(3) - (4)].val), &tmp);
	}
    break;

  case 778:
/* Line 1792 of yacc.c  */
#line 2932 "gram.y"
    {
	    (yyval.val) = M_PI;
	}
    break;

  case 779:
/* Line 1792 of yacc.c  */
#line 2935 "gram.y"
    {
	    (yyval.val) = M_PI / 180.0;
	}
    break;

  case 780:
/* Line 1792 of yacc.c  */
#line 2938 "gram.y"
    {
	    (yyval.val) = (double) drand48();
	}
    break;

  case 781:
/* Line 1792 of yacc.c  */
#line 2941 "gram.y"
    {
	    (yyval.val) = sin((yyvsp[(3) - (4)].val));
	}
    break;

  case 782:
/* Line 1792 of yacc.c  */
#line 2944 "gram.y"
    {
	    (yyval.val) = pow((yyvsp[(3) - (4)].val), 2.0);
	}
    break;

  case 783:
/* Line 1792 of yacc.c  */
#line 2947 "gram.y"
    {
	    (yyval.val) = sqrt((yyvsp[(3) - (4)].val));
	}
    break;

  case 784:
/* Line 1792 of yacc.c  */
#line 2950 "gram.y"
    {
	    (yyval.val) = tan((yyvsp[(3) - (4)].val));
	}
    break;

  case 785:
/* Line 1792 of yacc.c  */
#line 2953 "gram.y"
    {
	    if ((int) (yyvsp[(3) - (5)].val))
		(yyval.val) = (yyvsp[(5) - (5)].val);
	}
    break;

  case 786:
/* Line 1792 of yacc.c  */
#line 2957 "gram.y"
    {
	    if ((int) (yyvsp[(3) - (7)].val)) {
		(yyval.val) = (yyvsp[(5) - (7)].val);
	    } else {
		(yyval.val) = (yyvsp[(7) - (7)].val);
	    }
	}
    break;

  case 787:
/* Line 1792 of yacc.c  */
#line 2964 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) > (yyvsp[(3) - (3)].val);
	}
    break;

  case 788:
/* Line 1792 of yacc.c  */
#line 2967 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) < (yyvsp[(3) - (3)].val);
	}
    break;

  case 789:
/* Line 1792 of yacc.c  */
#line 2970 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) <= (yyvsp[(3) - (3)].val);
	}
    break;

  case 790:
/* Line 1792 of yacc.c  */
#line 2973 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) >= (yyvsp[(3) - (3)].val);
	}
    break;

  case 791:
/* Line 1792 of yacc.c  */
#line 2976 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) == (yyvsp[(3) - (3)].val);
	}
    break;

  case 792:
/* Line 1792 of yacc.c  */
#line 2979 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) != (yyvsp[(3) - (3)].val);
	}
    break;

  case 793:
/* Line 1792 of yacc.c  */
#line 2982 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) && (yyvsp[(3) - (3)].val);
	}
    break;

  case 794:
/* Line 1792 of yacc.c  */
#line 2985 "gram.y"
    {
	    (yyval.val) = (yyvsp[(1) - (3)].val) || (yyvsp[(3) - (3)].val);
	}
    break;

  case 795:
/* Line 1792 of yacc.c  */
#line 2988 "gram.y"
    {
	    (yyval.val) = !((yyvsp[(2) - (2)].val));
	}
    break;

  case 796:
/* Line 1792 of yacc.c  */
#line 2991 "gram.y"
    {
	    (yyval.val) = (yyvsp[(2) - (3)].val);
	}
    break;

  case 797:
/* Line 1792 of yacc.c  */
#line 2994 "gram.y"
    {
	    (yyval.val) = -(yyvsp[(2) - (2)].val);
	}
    break;


/* Line 1792 of yacc.c  */
#line 10922 "y.tab.c"
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
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;

  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYEMPTY : YYTRANSLATE (yychar);

  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
# define YYSYNTAX_ERROR yysyntax_error (&yymsg_alloc, &yymsg, \
                                        yyssp, yytoken)
      {
        char const *yymsgp = YY_("syntax error");
        int yysyntax_error_status;
        yysyntax_error_status = YYSYNTAX_ERROR;
        if (yysyntax_error_status == 0)
          yymsgp = yymsg;
        else if (yysyntax_error_status == 1)
          {
            if (yymsg != yymsgbuf)
              YYSTACK_FREE (yymsg);
            yymsg = (char *) YYSTACK_ALLOC (yymsg_alloc);
            if (!yymsg)
              {
                yymsg = yymsgbuf;
                yymsg_alloc = sizeof yymsgbuf;
                yysyntax_error_status = 2;
              }
            else
              {
                yysyntax_error_status = YYSYNTAX_ERROR;
                yymsgp = yymsg;
              }
          }
        yyerror (yymsgp);
        if (yysyntax_error_status == 2)
          goto yyexhaustedlab;
      }
# undef YYSYNTAX_ERROR
#endif
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

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
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
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
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
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

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

#if !defined yyoverflow || YYERROR_VERBOSE
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


/* Line 2055 of yacc.c  */
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
