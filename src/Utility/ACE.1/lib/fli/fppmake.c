/****************************************************************
 * fppmake.c
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
#include <sys/types.h>
#include <sys/stat.h>
#include <memory.h>
#include "fpfli.h"

#define MAXLEN 512

#define SEEK_SET 0

unsigned char *pixel_chunk;

static unsigned char file_header[FLI_FILE_HEADER_SIZE];
static LONG color1[FLI_MAX_COLORS];
static LONG color2[FLI_MAX_COLORS];
static LONG color3[FLI_MAX_COLORS];

int nframes;

/****************************************************************
 * make_fli
 ****************************************************************/

static int help, irun;
static int file_size, frame_size, mem_size_one;
static int read_flag, first_flag;
static UBYTE *pb_run[3];
static UBYTE *first_pixel, *last_pixel, *curr_pixel, *pre_pixel, *prepre_pixel;
static LONG *curr_color, *first_color, *last_color;

int initmake_fli()
{
    int i, j;
    /* ....... prepare pointers ....... */
    first_pixel = big_buffer;
    last_pixel = first_pixel + fli_size;
    pb_run[0] = last_pixel + fli_size;
    pb_run[1] = pb_run[0] + fli_size;
    pb_run[2] = pb_run[1] + fli_size;

    first_color = color1;
    last_color = color2;
    if (!get_image(first_pixel, first_color, 0))
	goto exit_error;
    if (!get_image(last_pixel, last_color, 0))
	goto exit_error;

    memset(file_header, 0, FLI_FILE_HEADER_SIZE);
    if (fwrite(file_header, FLI_FILE_HEADER_SIZE, 1, output) != 1) {
	fprintf(stderr, " write error\n");
	goto exit_error;
    }
    file_size = FLI_FILE_HEADER_SIZE;
    curr_pixel = last_pixel;
    pre_pixel = last_pixel;
    mem_size_one = 0;
    return (1);
exit_error:
    exitialise(1);
    return (0);
}

add_image(i)
    int i;
{
    irun = i % 3;
printf("Adding image %d\n", i);
    prepre_pixel = pre_pixel;
    pre_pixel = curr_pixel;

    if (i == 1) {
	read_flag = 0;
	first_flag = 1;
	curr_pixel = first_pixel;
	curr_color = first_color;
    } else if (i == nframes) {
	read_flag = 0;
	first_flag = 0;
	curr_pixel = last_pixel;
	curr_color = last_color;
    } else {
	read_flag = 1;
	first_flag = 0;
	curr_pixel = pb_run[irun];
	curr_color = color3;
    }

    if (read_flag == 1) {
	if (!get_image(curr_pixel, curr_color, 0)) {
	    goto exit_error;
	}
    }
    frame_size = fli_write_frame(prepre_pixel, pre_pixel, curr_pixel,
				 curr_color, first_flag);

    if (frame_size == 0) {
	fprintf(stderr, " Error writing frame %d\n", i);
	goto exit_error;
    }
    file_size += frame_size;

    if (i == 1)
	mem_size_one = frame_size;
    return (1);
exit_error:
    exitialise(1);
    return (0);
}

close_fli()
{
    prepre_pixel = pre_pixel;
    pre_pixel = curr_pixel;

    frame_size = fli_write_frame(prepre_pixel, pre_pixel, first_pixel,
				 first_color, 0);

    if (frame_size == 0) {
	fprintf(stderr, " Error writing frame\n");
	goto exit_error;
    }
    file_size += frame_size;

    /* ..................... write actual header ......... */

    help = 0;

    add_bytes(file_header, &help, file_size, IOM_LONG);
    add_bytes(file_header, &help, FLI_FILE_MAGIC, IOM_UWORD);
    add_bytes(file_header, &help, nframes, IOM_UWORD);
    add_bytes(file_header, &help, fli_width, IOM_UWORD);
    add_bytes(file_header, &help, fli_height, IOM_UWORD);
    add_bytes(file_header, &help, 0x0008, IOM_UWORD);
    add_bytes(file_header, &help, 0x0000, IOM_UWORD);
    add_bytes(file_header, &help, fli_speed, IOM_UWORD);

    help = 0x0050;
    add_bytes(file_header, &help, 0x0080, IOM_LONG);
    add_bytes(file_header, &help, (0x0080 + mem_size_one), IOM_LONG);

    if (fseek(output, 0L, SEEK_SET) != 0) {
	fprintf(stderr, " fseek error\n");
	goto exit_error;
    }
    if (fwrite(file_header, help, 1, output) != 1) {
	fprintf(stderr, " write error\n");
	goto exit_error;
    }
    exitialise(0);
    return (nframes);

exit_error:
    exitialise(1);
    return (0);
}
