/* +-------------------------------------------------------------------+ */
/* | Copyright 1990, David Koblas.                                     | */
/* |   Permission to use, copy, modify, and distribute this software   | */
/* |   and its documentation for any purpose and without fee is hereby | */
/* |   granted, provided that the above copyright notice appear in all | */
/* |   copies and that both that copyright notice and this permission  | */
/* |   notice appear in supporting documentation.  This software is    | */
/* |   provided "as is" without express or implied warranty.           | */
/* +-------------------------------------------------------------------+ */


#include <X11/Intrinsic.h>
#include <stdio.h>


#define	MAXCOLORMAPSIZE		256

#define	TRUE	1
#define	FALSE	0

#define CM_RED		0
#define CM_GREEN	1
#define CM_BLUE		2

#define	MAX_LWZ_BITS		12

#define INTERLACE		0x40
#define LOCALCOLORMAP	0x80
#define BitSet(byte, bit)	(((byte) & (bit)) == (bit))

#define	ReadOK(file,buffer,len)	(fread(buffer, len, 1, file) != 0)

#define LM_to_uint(a,b)			(((b)<<8)|(a))

struct {
    unsigned int Width;
    unsigned int Height;
    unsigned char ColorMap[3][MAXCOLORMAPSIZE];
    unsigned int BitPixel;
    unsigned int ColorResolution;
    unsigned int Background;
    unsigned int AspectRatio;
} GifScreen;

struct {
    int transparent;
    int delayTime;
    int inputFlag;
    int disposal;
} Gif89 = {

-1, -1, -1, 0};

int verbose;
int showComment;


/*
unsigned char *ReadGIF ARGS(( FILE	*fd, int imageNumber ));
static int ReadColorMap ARGS(( FILE *fd, int number, unsigned char buffer[3][MAXCOLORMAPSIZE] ));
static int DoExtension ARGS(( FILE *fd, int label ));
static int GetDataBlock ARGS(( FILE *fd, unsigned char  *buf ));
static int GetCode ARGS(( FILE *fd, int code_size, int flag ));
static int LWZReadByte ARGS(( FILE *fd, int flag, int input_code_size ));
static unsigned char *ReadImage ARGS(( FILE *fd, int len, int height, unsigned char cmap[3][MAXCOLORMAPSIZE], int interlace, int ignore ));
*/

static int ReadColorMap();
static int DoExtension();
static int GetDataBlock();
static int GetCode();
static int LWZReadByte();
static unsigned char *ReadImage();

unsigned char *ReadGIF(fd, w, h, colrs)
FILE *fd;
int *w, *h;
XColor *colrs;
{
    unsigned char buf[16];
    unsigned char *data;
    unsigned char c;
    unsigned char localColorMap[3][MAXCOLORMAPSIZE];
    int useGlobalColormap;
    int bitPixel;
    int imageCount = 0;
    char version[4];
    int imageNumber = 1;
    int i;

    verbose = FALSE;
    showComment = FALSE;

    if (!ReadOK(fd, buf, 6)) {
#if 0
	fprintf(stderr, "error reading magic number\n");
#endif
	return (NULL);
    }
    if (strncmp((char *) buf, "GIF", 3) != 0) {
#if 0
	if (verbose)
	    fprintf(stderr, "not a GIF file\n");
#endif
	return (NULL);
    }
    strncpy(version, (char *) buf + 3, 3);
    version[3] = '\0';

    if ((strcmp(version, "87a") != 0) && (strcmp(version, "89a") != 0)) {
#if 0
	fprintf(stderr, "bad version number, not '87a' or '89a'\n");
#endif
	return (NULL);
    }
    if (!ReadOK(fd, buf, 7)) {
#if 0
	fprintf(stderr, "failed to read screen descriptor\n");
#endif
	return (NULL);
    }
    GifScreen.Width = LM_to_uint(buf[0], buf[1]);
    GifScreen.Height = LM_to_uint(buf[2], buf[3]);
    GifScreen.BitPixel = 2 << (buf[4] & 0x07);
    GifScreen.ColorResolution = (((buf[4] & 0x70) >> 3) + 1);
    GifScreen.Background = buf[5];
    GifScreen.AspectRatio = buf[6];

    if (BitSet(buf[4], LOCALCOLORMAP)) {	/* Global Colormap */
	if (ReadColorMap(fd, GifScreen.BitPixel, GifScreen.ColorMap)) {
#if 0
	    fprintf(stderr, "error reading global colormap\n");
#endif
	    return (NULL);
	}
	for (i = 0; i < GifScreen.BitPixel; i++) {
	    int scale = 65536 / MAXCOLORMAPSIZE;

	    colrs[i].red = GifScreen.ColorMap[0][i] * scale;
	    colrs[i].green = GifScreen.ColorMap[1][i] * scale;
	    colrs[i].blue = GifScreen.ColorMap[2][i] * scale;
	    colrs[i].pixel = i;
	    colrs[i].flags = DoRed | DoGreen | DoBlue;
	}
	for (i = GifScreen.BitPixel; i < MAXCOLORMAPSIZE; i++) {
	    colrs[i].red = 0;
	    colrs[i].green = 0;
	    colrs[i].blue = 0;
	    colrs[i].pixel = i;
	    colrs[i].flags = DoRed | DoGreen | DoBlue;
	}

    }
    if (GifScreen.AspectRatio != 0 && GifScreen.AspectRatio != 49) {
	float r;
	r = ((float) GifScreen.AspectRatio + 15.0) / 64.0;
#if 0
	fprintf(stderr, "Warning:  non-square pixels!\n");
#endif
    }
    for (;;) {
	if (!ReadOK(fd, &c, 1)) {
#if 0
	    fprintf(stderr, "EOF / read error on image data\n");
#endif
	    return (NULL);
	}
	if (c == ';') {		/* GIF terminator */
	    if (imageCount < imageNumber) {
#if 0
		fprintf(stderr, "No images found in file\n");
#endif
		return (NULL);
	    }
	    break;
	}
	if (c == '!') {		/* Extension */
	    if (!ReadOK(fd, &c, 1)) {
#if 0
		fprintf(stderr, "EOF / read error on extention function code\n");
#endif
		return (NULL);
	    }
	    DoExtension(fd, c);
	    continue;
	}
	if (c != ',') {		/* Not a valid start character */
#if 0
	    fprintf(stderr, "bogus character 0x%02x, ignoring\n", (int) c);
#endif
	    continue;
	}
	++imageCount;

	if (!ReadOK(fd, buf, 9)) {
#if 0
	    fprintf(stderr, "couldn't read left/top/width/height\n");
#endif
	    return (NULL);
	}
	useGlobalColormap = !BitSet(buf[8], LOCALCOLORMAP);

	bitPixel = 1 << ((buf[8] & 0x07) + 1);

	*w = LM_to_uint(buf[4], buf[5]);
	*h = LM_to_uint(buf[6], buf[7]);
	if (!useGlobalColormap) {
	    if (ReadColorMap(fd, bitPixel, localColorMap)) {
#if 0
		fprintf(stderr, "error reading local colormap\n");
#endif
		return (NULL);
	    }
	    for (i = 0; i < bitPixel; i++) {
		int scale = 65536 / MAXCOLORMAPSIZE;

		colrs[i].red = localColorMap[0][i] * scale;
		colrs[i].green = localColorMap[1][i] * scale;
		colrs[i].blue = localColorMap[2][i] * scale;
		colrs[i].pixel = i;
		colrs[i].flags = DoRed | DoGreen | DoBlue;
	    }
	    for (i = bitPixel; i < MAXCOLORMAPSIZE; i++) {
		colrs[i].red = 0;
		colrs[i].green = 0;
		colrs[i].blue = 0;
		colrs[i].pixel = i;
		colrs[i].flags = DoRed | DoGreen | DoBlue;
	    }
	    data = ReadImage(fd, LM_to_uint(buf[4], buf[5]), LM_to_uint(buf[6], buf[7]), localColorMap, BitSet(buf[8], INTERLACE), imageCount != imageNumber);
	} else {
	    data = ReadImage(fd, LM_to_uint(buf[4], buf[5]), LM_to_uint(buf[6], buf[7]), GifScreen.ColorMap, BitSet(buf[8], INTERLACE), imageCount != imageNumber);
	}

    }
    return (data);
}

static int ReadColorMap(fd, number, buffer)
FILE *fd;
int number;
unsigned char buffer[3][MAXCOLORMAPSIZE];
{
    int i;
    unsigned char rgb[3];

    for (i = 0; i < number; ++i) {
	if (!ReadOK(fd, rgb, sizeof(rgb))) {
#if 0
	    fprintf(stderr, "bad colormap\n");
#endif
	    return (TRUE);
	}
	buffer[CM_RED][i] = rgb[0];
	buffer[CM_GREEN][i] = rgb[1];
	buffer[CM_BLUE][i] = rgb[2];
    }
    return FALSE;
}

static int DoExtension(fd, label)
FILE *fd;
int label;
{
    static char buf[256];
    char *str;

    switch (label) {
    case 0x01:			/* Plain Text Extension */
	str = "Plain Text Extension";
#ifdef notdef
	if (GetDataBlock(fd, (unsigned char *) buf) == 0);

	lpos = LM_to_uint(buf[0], buf[1]);
	tpos = LM_to_uint(buf[2], buf[3]);
	width = LM_to_uint(buf[4], buf[5]);
	height = LM_to_uint(buf[6], buf[7]);
	cellw = buf[8];
	cellh = buf[9];
	foreground = buf[10];
	background = buf[11];

	while (GetDataBlock(fd, (unsigned char *) buf) != 0) {
	    PPM_ASSIGN(image[ypos][xpos], cmap[CM_RED][v], cmap[CM_GREEN][v], cmap[CM_BLUE][v]);
	    ++index;
	}

	return FALSE;
#else
	break;
#endif
    case 0xff:			/* Application Extension */
	str = "Application Extension";
	break;
    case 0xfe:			/* Comment Extension */
	str = "Comment Extension";
	while (GetDataBlock(fd, (unsigned char *) buf) != 0) {
	    if (showComment) {
#if 0
		fprintf(stderr, "gif comment: %s\n", buf);
#endif
	    }
	}
	return FALSE;
    case 0xf9:			/* Graphic Control Extension */
	str = "Graphic Control Extension";
	(void) GetDataBlock(fd, (unsigned char *) buf);
	Gif89.disposal = (buf[0] >> 2) & 0x7;
	Gif89.inputFlag = (buf[0] >> 1) & 0x1;
	Gif89.delayTime = LM_to_uint(buf[1], buf[2]);
	if ((buf[0] & 0x1) != 0)
	    Gif89.transparent = buf[3];

	while (GetDataBlock(fd, (unsigned char *) buf) != 0);
	return FALSE;
    default:
	str = buf;
	sprintf(buf, "UNKNOWN (0x%02x)", label);
	break;
    }

    /* fprintf(stderr, "got a '%s' extension\n", str); */

    while (GetDataBlock(fd, (unsigned char *) buf) != 0);

    return FALSE;
}

static int ZeroDataBlock = FALSE;

static int GetDataBlock(fd, buf)
FILE *fd;
unsigned char *buf;
{
    unsigned char count;

    if (!ReadOK(fd, &count, 1)) {
#if 0
	fprintf(stderr, "error in getting DataBlock size\n");
#endif
	return -1;
    }
    ZeroDataBlock = count == 0;

    if ((count != 0) && (!ReadOK(fd, buf, count))) {
#if 0
	fprintf(stderr, "error in reading DataBlock\n");
#endif
	return -1;
    }
    return count;
}

static int GetCode(fd, code_size, flag)
FILE *fd;
int code_size;
int flag;
{
    static unsigned char buf[280];
    static int curbit, lastbit, done, last_byte;
    int i, j, ret;
    unsigned char count;

    if (flag) {
	curbit = 0;
	lastbit = 0;
	done = FALSE;
	return 0;
    }
    if ((curbit + code_size) >= lastbit) {
	if (done) {
	    if (curbit >= lastbit) {
#if 0
		fprintf(stderr, "ran off the end of my bits\n");
#endif
	    }
	    return -1;
	}
	buf[0] = buf[last_byte - 2];
	buf[1] = buf[last_byte - 1];

	if ((count = GetDataBlock(fd, &buf[2])) == 0)
	    done = TRUE;

	last_byte = 2 + count;
	curbit = (curbit - lastbit) + 16;
	lastbit = (2 + count) * 8;
    }
    ret = 0;
    for (i = curbit, j = 0; j < code_size; ++i, ++j)
	ret |= ((buf[i >> 3] & (1 << (i % 8))) != 0) << j;

    curbit += code_size;

    return ret;
}

static int LWZReadByte(fd, flag, input_code_size)
FILE *fd;
int flag;
int input_code_size;
{
    static int fresh = FALSE;
    int code, incode;
    static int code_size, set_code_size;
    static int max_code, max_code_size;
    static int firstcode, oldcode;
    static int clear_code, end_code;
    static int table[2][(1 << MAX_LWZ_BITS)];
    static int stack[(1 << (MAX_LWZ_BITS)) * 2], *sp;
    register int i;

    if (flag) {
	set_code_size = input_code_size;
	code_size = set_code_size + 1;
	clear_code = 1 << set_code_size;
	end_code = clear_code + 1;
	max_code_size = 2 * clear_code;
	max_code = clear_code + 2;

	GetCode(fd, 0, TRUE);

	fresh = TRUE;

	for (i = 0; i < clear_code; ++i) {
	    table[0][i] = 0;
	    table[1][i] = i;
	}
	for (; i < (1 << MAX_LWZ_BITS); ++i)
	    table[0][i] = table[1][0] = 0;

	sp = stack;

	return 0;
    } else if (fresh) {
	fresh = FALSE;
	do {
	    firstcode = oldcode = GetCode(fd, code_size, FALSE);
	} while (firstcode == clear_code);
	return firstcode;
    }
    if (sp > stack)
	return *--sp;

    while ((code = GetCode(fd, code_size, FALSE)) >= 0) {
	if (code == clear_code) {
	    for (i = 0; i < clear_code; ++i) {
		table[0][i] = 0;
		table[1][i] = i;
	    }
	    for (; i < (1 << MAX_LWZ_BITS); ++i)
		table[0][i] = table[1][i] = 0;
	    code_size = set_code_size + 1;
	    max_code_size = 2 * clear_code;
	    max_code = clear_code + 2;
	    sp = stack;
	    firstcode = oldcode = GetCode(fd, code_size, FALSE);
	    return firstcode;
	} else if (code == end_code) {
	    int count;
	    unsigned char buf[260];

	    if (ZeroDataBlock)
		return -2;

	    while ((count = GetDataBlock(fd, buf)) > 0);

#if 0
	    if (count != 0)
		fprintf(stderr, "missing EOD in data stream (common occurence)\n");
#endif
	    return -2;
	}
	incode = code;

	if (code >= max_code) {
	    *sp++ = firstcode;
	    code = oldcode;
	}
	while (code >= clear_code) {
	    if ((sp - stack) >= ((1 << (MAX_LWZ_BITS)) * 2))
		return -2;	/* stop a code dump */
	    *sp++ = table[1][code];
	    if (code == table[0][code]) {
#if 0
		fprintf(stderr, "circular table entry BIG ERROR\n");
#endif
		return (code);
	    }
	    code = table[0][code];
	}

	*sp++ = firstcode = table[1][code];

	if ((code = max_code) < (1 << MAX_LWZ_BITS)) {
	    table[0][code] = oldcode;
	    table[1][code] = firstcode;
	    ++max_code;
	    if ((max_code >= max_code_size) && (max_code_size < (1 << MAX_LWZ_BITS))) {
		max_code_size *= 2;
		++code_size;
	    }
	}
	oldcode = incode;

	if (sp > stack)
	    return *--sp;
    }
    return code;
}

static unsigned char *ReadImage(fd, len, height, cmap, interlace, ignore)
FILE *fd;
int len, height;
unsigned char cmap[3][MAXCOLORMAPSIZE];
int interlace, ignore;
{
    unsigned char c;
    int v;
    int xpos = 0, ypos = 0, pass = 0;
/*
	pixel		**image;
*/
    unsigned char *data;
    unsigned char *dptr;

    /*
     * *  Initialize the Compression routines
     */
    if (!ReadOK(fd, &c, 1)) {
#if 0
	fprintf(stderr, "EOF / read error on image data\n");
#endif
	return (NULL);
    }
    if (LWZReadByte(fd, TRUE, c) < 0) {
#if 0
	fprintf(stderr, "error reading image\n");
#endif
	return (NULL);
    }
    /*
     * *  If this is an "uninteresting picture" ignore it.
     */
    if (ignore) {
#if 0
	if (verbose)
	    fprintf(stderr, "skipping image...\n");
#endif

	while (LWZReadByte(fd, FALSE, c) >= 0);
	return;
    }
/*
	if ((image = ppm_allocarray(len, height)) == NULL)
	{
		fprintf(stderr, "couldn't alloc space for image\n");
		return(NULL);
	}
*/

    data = (unsigned char *) malloc(len * height);
    if (data == NULL) {
#if 0
	fprintf(stderr, "Cannot allocate space for image data\n");
#endif
	return (NULL);
    }
#if 0
    if (verbose)
	fprintf(stderr, "reading %d by %d%s GIF image\n", len, height, interlace ? " interlaced" : "");
#endif

    while ((v = LWZReadByte(fd, FALSE, c)) >= 0) {
	dptr = (unsigned char *) (data + (ypos * len) + xpos);
	*dptr = (unsigned char) v;
/*
		PPM_ASSIGN(image[ypos][xpos], cmap[CM_RED][v],
					cmap[CM_GREEN][v], cmap[CM_BLUE][v]);
*/

	++xpos;
	if (xpos == len) {
	    xpos = 0;
	    if (interlace) {
		switch (pass) {
		case 0:
		case 1:
		    ypos += 8;
		    break;
		case 2:
		    ypos += 4;
		    break;
		case 3:
		    ypos += 2;
		    break;
		}

		if (ypos >= height) {
		    ++pass;
		    switch (pass) {
		    case 1:
			ypos = 4;
			break;
		    case 2:
			ypos = 2;
			break;
		    case 3:
			ypos = 1;
			break;
		    default:
			goto fini;
		    }
		}
	    } else {
		++ypos;
	    }
	}
	if (ypos >= height)
	    break;
    }

  fini:
    if (LWZReadByte(fd, FALSE, c) >= 0)
	/* fprintf(stderr, "too much input data, ignoring extra...\n") */ ;

/*
	if (verbose)
		pm_message("writing output");
	ppm_writeppm(stdout, image, len, height, (pixval) 255, 0 );
*/
    return (data);
}
