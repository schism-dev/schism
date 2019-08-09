/*
 * This file: image_f_io.c (part of the WSCRAWL program)
 *
 * This file contains the "Image File I/O" package for wscrawl (or anything
 * for that matter.)  The format used is the standard X-Window Dump form
 * that the MIT client "xwd" uses.  This File I/O was made possible by the 
 * help and extensive source code of Mark Cook of Hewlett-Packard, who I 
 * bothered endlessly to get this working.  
 * 
 * I tried to make this file of routines as self-contained and portable as
 * possible.  Please feel free to use this file as is, or with modifications,
 * for any and all applications.  -- Brian Wilson
 *
 * This file was last modified: 10/7/91 (initial port to wscrawl 2.0)
 */

#define TRUE		1	
#define FALSE		0	
#define BAD_CHOICE		-1
#define NO_DUPLICATE   		-1
#define UNALLOCATED    		-1

#include <X11/Xos.h>
#include <X11/XWDFile.h>
#include <X11/Xlib.h>
#include <stdio.h>

extern char *malloc();
XImage *read_image_from_disk();


/*
 */
get_xwd_image(disp, info_win_id, win_id, x, y, width, height, name_of_file,
		   the_colormap)
Display *disp;
Window info_win_id;
Pixmap win_id;
int x, y, width, height;
char *name_of_file;
Colormap the_colormap;
{
    unsigned long swaptest = TRUE;
    XRectangle box2;
    XRectangle *box;
    FILE *out_file_ptr;
    XColor *colors;
    unsigned buffer_size;
    int win_name_size;
    int header_size;
    int format=ZPixmap;
    int ncolors, i;
    XWindowAttributes win_info;
    XImage *ImagePix;
    XWDFileHeader header;

    /*
     * convert to form this code understands (I got code from Mark Cook)
     */
    box = &box2;
    box->x = x;
    box->y = y;
    box->width = width;
    box->height = height;

    /*
     * Open the file in which the image is to be stored
     */
    if ((out_file_ptr = fopen(name_of_file, "w")) == NULL) 
    {
        printf("ERROR: Could not open file %s.\n", name_of_file);
        return(0);
    }
    else
    {         /* Dump the image to the specified file */
	if (the_colormap == NULL)
	{
	    the_colormap =
		       XDefaultColormapOfScreen(XDefaultScreenOfDisplay(disp));
	}

        if (!XGetWindowAttributes(disp, info_win_id, &win_info)) 
        {
            printf("Can't get window attributes.\n");
            return(0);
        }
    
        /*
         * sizeof(char) is included for the null string terminator. 
         */
        win_name_size = strlen(name_of_file) + sizeof(char);
    
        ImagePix = XGetImage(disp, win_id, box->x, box->y, box->width,
                               box->height, AllPlanes, format); 
        XFlush(disp);
    
        if (ImagePix == NULL) 
        {
            printf("GetImage failed.\n");
            return(0);
        }
    
        buffer_size = Image_Size(ImagePix,format);/*determines size of pixmap*/
    
        /*
         * Get the RGB values for the current color cells
         */
        if ((ncolors = Get_Colors(&colors, disp, the_colormap)) == 0) 
        {
            printf("Cannot alloc memory for color structs.\n");
            return(0);
        }
        XFlush(disp);
    
        header_size = sizeof(header) +win_name_size; /*Calculates header size*/
    
        /*
         * Assemble the file header information
         */
        header.header_size = (xwdval) header_size;
        header.file_version = (xwdval) XWD_FILE_VERSION;
        header.pixmap_format = (xwdval) format;
        header.pixmap_depth = (xwdval) ImagePix->depth;

        header.pixmap_width = (xwdval) ImagePix->width;
        header.pixmap_height = (xwdval) ImagePix->height;
        header.xoffset = (xwdval) ImagePix->xoffset;
        header.byte_order = (xwdval) ImagePix->byte_order;
        header.bitmap_unit = (xwdval) ImagePix->bitmap_unit;
        header.bitmap_bit_order = (xwdval) ImagePix->bitmap_bit_order;
        header.bitmap_pad = (xwdval) ImagePix->bitmap_pad;
        header.bits_per_pixel = (xwdval) ImagePix->bits_per_pixel;
        header.bytes_per_line = (xwdval) ImagePix->bytes_per_line;
        header.visual_class = (xwdval) win_info.visual->class;
        header.red_mask = (xwdval) win_info.visual->red_mask;
        header.green_mask = (xwdval) win_info.visual->green_mask;
        header.blue_mask = (xwdval) win_info.visual->blue_mask;
        header.bits_per_rgb = (xwdval) win_info.visual->bits_per_rgb;
        header.colormap_entries = (xwdval) win_info.visual->map_entries;
        header.ncolors = ncolors;
        header.window_width = (xwdval) ImagePix->width;
        header.window_height = (xwdval) ImagePix->height;
        header.window_x = (xwdval) 0;
        header.window_y = (xwdval) 0;
        header.window_bdrwidth = (xwdval) 0;
      
        if (*(char *) &swaptest) 
        {
            _swaplong((char *) &header, sizeof(header));
            for (i = 0; i < ncolors; i++) 
	    {
                _swaplong((char *) &colors[i].pixel, sizeof(long));
                _swapshort((char *) &colors[i].red, 3 * sizeof(short));
            }
        }
    
        /*
         * Write out the file header information
         */
        (void) fwrite((char *)&header, sizeof(header), 1, out_file_ptr);
        (void) fwrite(name_of_file, win_name_size, 1, out_file_ptr);
    
        /*
         * Write out the color cell RGB values
         */
        (void) fwrite((char *) colors, sizeof(XColor), ncolors, out_file_ptr);
    
        /*
         * Write out the buffer
         */
        (void) fwrite(ImagePix->data, (int) buffer_size, 1, out_file_ptr);
    
        if(ncolors > 0) 
            free(colors);    /*free the color buffer*/
    
        fclose(out_file_ptr);
        XFlush(disp);
	XDestroyImage(ImagePix);
    }
}
