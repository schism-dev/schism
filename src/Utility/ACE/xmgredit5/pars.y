/*
 * 
 * evaluate an expression
 * 
 */

%{

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include <string.h>
#include <signal.h>
#include <setjmp.h>

#include "defines.h"

char *gettxt (const char *msgid, const char *dflt_str);
void yyerror(const char *s);

double drand48();
long lrand48();

#ifndef lint
static char RCSid[] = "$Id: pars.y,v 1.3 2007/02/21 00:21:21 pturner Exp $";
#endif

double result; /* return value if expression */
static int interr;
static double tmp;

jmp_buf begin;

static char *progname;
static int lineno = 1;
char f_string[256];
static int pos = 0;

static double *aa, *bb, *cc, *dd, *xx, *yy;
static int setindex;
static int setsetno;

%}
 
%union {
	double val;
	double *ptr;
	int func;
}
 
%token	<val> NUMBER
%token	<ptr> VAR 
%token	<func> CEIL FLOOR MOD TAN PI ABS SQR LGAMMA LOG LN
%token	<func> ERF ERFC EXP SIN COS ACOS ASIN ATAN2 ATAN SQRT RAND HYPOT
%token	<func> DEG DX DY RAD MAXP MINP INDEX SETNO INT INVN INVT IRAND NORM NORMP RNORM
%type	<val> expr
%type	<ptr> asgn
%right	'='
%left	OR
%left	AND
%nonassoc GT LT LE GE EQ NE
%left	'+' '-'
%left	'*' '/' '%'
%right '^'
%right	UMINUS NOT
%%
list:
	| list '\n'
	| list asgn '\n'
	| list expr '\n' { result=$2; }
	| error '\n' {yyerror("error:"); yyerrok; }
	;
asgn:	VAR '=' expr { *($$)=$3; result = $3; }
	;
expr: NUMBER
	| VAR	{ $$ = *($1); }
	| expr '+' expr 		{ $$ = $1 + $3; }
	| expr '-' expr 		{ $$ = $1 - $3; }
	| expr '*' expr 		{ $$ = $1 * $3; }
	| expr '/' expr 		{ if ($3 !=0.0) {
						$$ = $1 / $3; 
					  }
					  else {
    						yyerror("Divide by zero");
					  }
					}
	| expr '%' expr 		{ $$ = fmod($1,$3); }
	| expr '^' expr 		{ $$ = pow($1,$3); }
	| ABS '(' expr ')'		{ $$ = fabs($3); }
	| ACOS '(' expr ')'		{ $$ = acos($3); }
	| ASIN '(' expr ')'		{ $$ = asin($3); }
	| ATAN '(' expr ')'		{ $$ = atan($3); }
	| ATAN2 '(' expr ',' expr ')'	{ $$ = atan2($3, $5); }
	| CEIL '(' expr ')'		{ $$ = ceil($3); }
	| COS '(' expr ')'		{ $$ = cos($3); }
	| DEG 				{ $$ = 180.0/M_PI; }
	| DX 				{ $$ = *xx; }
	| DY 				{ $$ = *yy; }
	| ERF '(' expr ')'		{ $$ = erf($3); }
	| ERFC '(' expr ')'		{ $$ = erfc($3); }
	| EXP '(' expr ')'		{ $$ = exp($3); }
	| FLOOR '(' expr ')'		{ $$ = floor($3); }
	| HYPOT '(' expr ',' expr ')'	{ $$ = hypot($3, $5); }
	| INDEX				{ $$ = setindex; }
	| SETNO				{ $$ = setsetno; }
	| INT '(' expr ')'		{ $$ = (long) $3; }
	| IRAND	'(' expr ')'		{ $$ = lrand48(); }
	| LGAMMA '(' expr ')'		{ $$ = lgamma($3); }
	| LN '(' expr ')'		{ $$ = log($3); }
	| LOG '(' expr ')'		{ $$ = log10($3); }
	| MAXP '(' expr ',' expr ')'	{ $$ = $3 >= $5 ? $3 : $5; }
	| MINP '(' expr ',' expr ')'	{ $$ = $3 <= $5 ? $3 : $5; }
	| MOD '(' expr ',' expr ')'	{ $$ = fmod($3,$5); }
	| PI 				{ $$ = M_PI; }
	| RAD 				{ $$ = M_PI/180.0; }
	| RAND				{ $$ = drand48(); }
	| SIN '(' expr ')'		{ $$ = sin($3); }
	| SQR '(' expr ')'		{ $$ = pow($3,2.0); }
	| SQRT '(' expr ')'		{ $$ = sqrt($3); }
	| TAN '(' expr ')'		{ $$ = tan($3); }
	| expr '?' expr ':' expr		{ $$ = $1 ? $3 : $5; }
	| expr GT expr		{ $$ = $1 > $3; }
	| expr LT expr		{ $$ = $1 < $3; }
	| expr LE expr		{ $$ = $1 <= $3; }
	| expr GE expr		{ $$ = $1 >= $3; }
	| expr EQ expr		{ $$ = $1 == $3; }
	| expr NE expr		{ $$ = $1 != $3; }
	| expr AND expr		{ $$ = $1 && $3; }
	| expr OR expr		{ $$ = $1 || $3; }
	| NOT expr		{ $$ = !($2); }
	| '(' expr ')' { $$ = $2; }
	| '-' expr  %prec UMINUS { $$ = -$2; }
	;
%%

void fixupstr(val)
    char val[];

{
    lowtoupper(val);
    val[strlen(val) + 1] = 0;
    val[strlen(val)] = '\n';
}

void scanner(s, x, y, a, b, c, d, i, setno, errpos)
    char s[];
double *x, *y, *a, *b, *c, *d;
int i, setno, *errpos;

{
    interr = 0;
    if (setjmp(begin)) {
	*errpos = interr;
	return;
    }
    pos = 0;
    aa = a;
    bb = b;
    cc = c;
    dd = d;
    xx = x;
    yy = y;
    setindex = i + 1;
    setsetno = setno;
    strcpy(f_string, s);
    yyparse();
    *errpos = interr;
}

struct funcs {
    char *s;
    int type;
} key[] = {
    "A", VAR,
    "ABS", ABS,
    "ACOS", ACOS,
    "ASIN", ASIN,
    "ATAN", ATAN,
    "ATAN2", ATAN2,
    "B", VAR,
    "C", VAR,
    "CEIL", CEIL,
    "COS", COS,
    "D", VAR,
    "DEG", DEG,
    "DEPTH", VAR,
    "DX", DX,
    "DY", DY,
    "ERF", ERF,
    "ERFC", ERFC,
    "EXP", EXP,
    "FLOOR", FLOOR,
    "HYPOT", HYPOT,
    "INDEX", INDEX,
    "INT", INT,
    "IRAND", IRAND,
    "LGAMMA", LGAMMA,
    "LN", LN,
    "LOG", LOG,
    "MAX", MAXP,
    "MIN", MINP,
    "MOD", MOD,
    "PI", PI,
    "RAD", RAD,
    "RAND", RAND,
    "RNORM", RNORM,
    "SETNO", SETNO,
    "SIN", SIN,
    "SQR", SQR,
    "SQRT", SQRT,
    "TAN", TAN,
    "X", VAR,
    "Y", VAR
};

static int maxfunc = sizeof(key) / sizeof(symtab_entry);

int findf(key, s, tlen)
    struct funcs key[];
char *s;
int tlen;

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

int getcharstr()
{
    if (pos >= strlen(f_string))
	return EOF;
    return (f_string[pos++]);
}

void ungetchstr()
{
    if (pos > 0)
	pos--;
}

int yylex()
{
    int c;
    int found;

    while ((c = getcharstr()) == ' ' || c == '\t');
    if (c == EOF)
	return (0);
    if (c == '.' || isdigit(c)) {
	char stmp[80];
	double d;
	int i;

	i = 0;
	while (c == '.' || isdigit(c)) {
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
	    while (c == '.' || isdigit(c)) {
		stmp[i++] = c;
		c = getcharstr();
	    }
	}
	stmp[i] = '\0';
	ungetchstr();
	sscanf(stmp, "%lf", &d);
	yylval.val = d;
	return NUMBER;
    }
    if (isalpha(c)) {
	char sbuf[100], *p = sbuf;

	do {
	    *p++ = c;
	}
	while ((c = getcharstr()) != EOF && isalnum(c));
	ungetchstr();
	*p = '\0';
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
		    if (strcmp(sbuf,"DEPTH")) {
		        yylval.ptr = dd;
		    }
		    else {
		        yylval.ptr = aa;
		    }
		    return VAR;
		case 'X':
		    yylval.ptr = xx;
		    return VAR;
		case 'Y':
		    yylval.ptr = yy;
		    return VAR;
		}
	    }
	    yylval.func = key[found].type;
	    return key[found].type;
	} else {
	    strcat(sbuf, ": No such function or variable");
	    yyerror(sbuf);
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

    if (c == expect)
	return ifyes;
    ungetchstr();
    return ifno;
}

void yyerror(const char *s)
{
    interr = 1;
    errwin(s);
    longjmp(begin, 1);
}
