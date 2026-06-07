/****************************************************************
 * fppfile.c
 ****************************************************************/

/******
  Copyright (C) 1993 by Klaus Ehrenfried. 
  
  Permission to use, copy, modify, and distribute this software
  is hereby granted, provided that the above copyright notice appears 
  in all copies and that the software is available to all free of charge. 
  The author disclaims all warranties with regard to this software, 
  including all implied warranties of merchant-ability and fitness. 
  The code is simply distributed as it is.
*******/

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "fpfli.h"

#ifdef WIN32
#ifdef LONG
#undef LONG
#endif
#endif

#define ERRMSG "Image has wrong format !\n"

#define TRUE		1
#define FALSE		0
#define BAD_CHOICE		-1
#define NO_DUPLICATE   		-1
#define UNALLOCATED    		-1

#include <X11/Xos.h>
#include <X11/XWDFile.h>
#include <X11/Xlib.h>

XImage *read_image_from_disk();

/*
save_image_on_disk(disp, cwin, displaybuff, 0, 0, win_w, win_h, "tmp.xwd", NULL);
 */
extern Window cwin;
extern Pixmap displaybuff;
extern int win_w, win_h;
static XImage *ImagePix;

/* An FBM bitmap header in memory */
typedef struct fbm_hdr_struct {
    int cols;			/* Width in pixels */
    int rows;			/* Height in pixels */
    int planes;			/* Depth (1 for B+W, 3 for RGB) */
    int bits;			/* Bits per pixel */
    int physbits;		/* Bits to store each pixel */
    int rowlen;			/* Length of a row in bytes */
    int plnlen;			/* Length of a plane in bytes */
    int clrlen;			/* Length of color map */
    double aspect;		/* ratio of Y to X of one pixel */
} FBMHDR;

/* An FBM bitmap in memory */
typedef struct fbm_struct {
    FBMHDR hdr;			/* Bitmap header */
    unsigned char *cm;		/* Pointer to colormap */
    unsigned char *bm;		/* Pointer to raw bits */
} FBM;

static unsigned char red[256], green[256], blue[256];
static int ccols[256];

int
 get_image(UBYTE * data,
	    LONG color[],
	    int without_data)
{
    FBM image;			/* Image */
    int ncolor, n2color;
    UBYTE *fbm_bm, *pdest, *psource;
    LONG rgb_value;
    int i, j, len, unass, nhelp, image_width;
    int idstart, idend, jdstart, jdend, isstart, jsstart;
    int x_origin, y_origin;
    int histogram[FLI_MAX_COLORS];
    static int ccnt;
    image.bm = (unsigned char *) NULL;
    if (!get_xwd_image(&image)) {
	exitialise(1);
	exit(1);
    }
    ncolor = (image.hdr.clrlen) / 3;
    n2color = ncolor + ncolor;

    for (j = 0; j < FLI_MAX_COLORS; j++) {
	histogram[j] = 0;

	if (j < ncolor) {
	    rgb_value = (long int) blue[j];
	    rgb_value = 256L * rgb_value + (long int) green[j];
	    rgb_value = 256L * rgb_value + (long int) red[j];
	    color[j] = rgb_value;
	} else {
	    color[j] = -1;
	}
    }


    fbm_bm = image.bm;
    image_width = image.hdr.rowlen;

    for (i = 0; i < image.hdr.plnlen; i++) {	/* compute histogram */
	histogram[*(fbm_bm++)]++;
    }

    unass = 0;
    for (j = 0; j < FLI_MAX_COLORS; j++) {
	if ((histogram[j] != 0) && (color[j] == -1)) {
	    color[j] = 0;
	    unass++;
	}
    }

    if (unass != 0) {
	fprintf(stderr, "Warning: %d unassigned color(s) referenced\n", unass);
    }
    if (color[0] == -1)
	color[0] = 0;
    if (color[border_color] == -1)
	color[border_color] = 0;

    memset(data, border_color, fli_size);

    if (Xorigin_flag == 1) {
	x_origin = Xorigin;
    } else {
	nhelp = fli_width - image.hdr.cols;
	x_origin = nhelp / 2;
    }

    if (x_origin >= 0) {
	idstart = x_origin;
	isstart = 0;
    } else {
	idstart = 0;
	isstart = -x_origin;
    }
    nhelp = x_origin + image.hdr.cols;
    idend = (nhelp < fli_width) ? nhelp : fli_width;

    if (Yorigin_flag == 1) {
	y_origin = Yorigin;
    } else {
	nhelp = fli_height - image.hdr.rows;
	y_origin = nhelp / 2;
    }

    if (y_origin >= 0) {
	jdstart = y_origin;
	jsstart = 0;
    } else {
	jdstart = 0;
	jsstart = -y_origin;
    }
    nhelp = y_origin + image.hdr.rows;
    jdend = (nhelp < fli_height) ? nhelp : fli_height;

    psource = image.bm + (jsstart * image_width + isstart);
    pdest = data + (jdstart * fli_width + idstart);

    len = idend - idstart;

    if (len > 0) {
	for (j = jdstart; j < jdend; j++) {
	    memcpy(pdest, psource, len);
	    psource += image_width;
	    pdest += fli_width;
	}
    }
    free(image.bm);
    return (1);
}

int get_xwd_image(image)
    FBM *image;
{
    unsigned long swaptest = TRUE;
    XRectangle box2;
    XRectangle *box;
    FILE *out_file_ptr;
    XColor *colors;
    unsigned buffer_size;
    int win_name_size;
    int header_size;
    int format = ZPixmap;
    int ncolors, i;
    XWindowAttributes win_info;
    XWDFileHeader header;
    int x = 0, y = 0, width = win_w, height = win_h;
    Colormap the_colormap = NULL;
    extern Display *disp;
    Window info_win_id = cwin;
    Pixmap win_id = displaybuff;

    /*
     * convert to form this code understands (I got code from Mark Cook) 
     */
    box = &box2;
    box->x = x;
    box->y = y;
    box->width = width;
    box->height = height;

    the_colormap =
	    XDefaultColormapOfScreen(XDefaultScreenOfDisplay(disp));

    if (!XGetWindowAttributes(disp, info_win_id, &win_info)) {
	printf("Can't get window attributes.\n");
	return (0);
    }
    ImagePix = XGetImage(disp, win_id, box->x, box->y, box->width,
			 box->height, AllPlanes, format);
    XFlush(disp);
    if (ImagePix == NULL) {
	printf("GetImage failed.\n");
	return (0);
    }
    buffer_size = Image_Size(ImagePix, format);	/* determines size of pixmap */

    image->hdr.rowlen = ImagePix->bytes_per_line;
    image->hdr.cols = ImagePix->width;
    image->hdr.rows = ImagePix->height;
    image->hdr.plnlen = ImagePix->width * ImagePix->height;
    image->bm = malloc(buffer_size * sizeof(char));
    memcpy(image->bm, ImagePix->data, buffer_size);

    /*
     * Get the RGB values for the current color cells 
     */
    if ((ncolors = Get_Colors(&colors, disp, the_colormap)) == 0) {
	printf("Cannot alloc memory for color structs.\n");
	return (0);
    }
    image->hdr.clrlen = 3 * ncolors;
    XFlush(disp);

    for (i = 0; i < ncolors; i++) {
	ccols[i] = colors[i].pixel;
	red[i] = colors[i].red >> 8;
	green[i] = colors[i].green >> 8;
	blue[i] = colors[i].blue >> 8;
    }
    if (ncolors > 0)
	free(colors);		/* free the color buffer */

    XFlush(disp);
    XDestroyImage(ImagePix);
}
