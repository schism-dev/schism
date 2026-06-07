/****************************************************************
 * fppcolor.c
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
#include "fpfli.h"

static int work[FLI_MAX_COLORS];
static long int mem_color[FLI_MAX_COLORS];
static int change[FLI_MAX_COLORS];

/****************************************************************
 * make_256_color_chunk
 ****************************************************************/

int
 make_256_color_chunk
 (
   LONG color[],		/* array of colors */
   int first_flag
) {
    LONG rgb_value;
    int red, green, blue, used_count;
    int skip_count, size_count, help, packets, change_count;
    int i, j, ipos;

    used_count = 0;
    change_count = 0;

    if (first_flag == 1) {	/* 1st frame */
	if (map_color_flag == 1) {
	    for (j = 0; j < FLI_MAX_COLORS; j++) {
		change[j] = 1;
		change_count++;
		color[j] = map_color[j];
		if (color[j] >= 0)
		    used_count++;
	    }
	} else {
	    for (j = 0; j < FLI_MAX_COLORS; j++) {
		change[j] = 1;
		change_count++;
		mem_color[j] = color[j];
		if (color[j] >= 0)
		    used_count++;
	    }
	}
    } else {
	if (map_color_flag == 1)
	    return (0);

	for (j = 0; j < FLI_MAX_COLORS; j++) {
	    if (color[j] >= 0) {/* used */
		used_count++;
		if (mem_color[j] != color[j]) {
		    change[j] = 1;
		    change_count++;
		    mem_color[j] = color[j];
		}
	    }
	}
    }

    fprintf(stdout, " Used colors: %d   new: %d\n", used_count, change_count);

    if (change_count == 0)
	return (0);		/* no new colors  -> no color_chunk */

    if (change[FLI_MAX_COLORS - 1] == 1)
	work[FLI_MAX_COLORS - 1] = 1;
    else
	work[FLI_MAX_COLORS - 1] = -1;

    for (i = (FLI_MAX_COLORS - 2); i >= 0; i--) {
	if (change[i] == 1)
	    work[i] = (work[i + 1] > 0) ? (work[i + 1] + 1) : 1;
	else
	    work[i] = (work[i + 1] < 0) ? (work[i + 1] - 1) : -1;
	/* fprintf(stdout," color: %d  %d\n",i,work[i]); */
    }

    ipos = 8;			/* 4 bytes for size of chunk */
    /* 2 bytes for type of chunk */
    /* 2 bytes for number of packets */
    i = 0;
    skip_count = 0;
    packets = 0;

    while (i < FLI_MAX_COLORS) {
	/* fprintf(stdout," color: %d  %d\n",i,work[i]); */
	if (work[i] < 0) {
	    skip_count = -work[i];
	    i += skip_count;
	} else {
	    size_count = work[i];
	    help = i + size_count;
	    if (size_count == 256)
		size_count = 0;

	    /*
	     * fprintf(stdout," skip: %d  size: %d\n",skip_count,
	     * size_count); 
	     */
	    add_bytes(color_chunk_buffer, &ipos, skip_count, IOM_UBYTE);
	    add_bytes(color_chunk_buffer, &ipos, size_count, IOM_UBYTE);

	    while (i < help) {
		rgb_value = color[i];
		if (rgb_value != -1) {
		    red = rgb_value % 256;
		    rgb_value = (rgb_value - red) / 256;
		    green = rgb_value % 256;
		    rgb_value = (rgb_value - green) / 256;
		    blue = rgb_value % 256;

		    /* printf(" %d   %d %d %d\n",i,red,green,blue); */
		} else {
		    red = 0;
		    green = 0;
		    blue = 0;
		}

		add_bytes(color_chunk_buffer, &ipos, red, IOM_UBYTE);
		add_bytes(color_chunk_buffer, &ipos, green, IOM_UBYTE);
		add_bytes(color_chunk_buffer, &ipos, blue, IOM_UBYTE);
		i++;
	    }
	    packets++;
	}
    }

    if (packets == 0)
	return (0);		/* no packets -> no color_chunk */

    if ((ipos % 2) == 1)	/* add single pad byte to even size */
	add_bytes(color_chunk_buffer, &ipos, 0x0000, IOM_UBYTE);

    help = 0;
    add_bytes(color_chunk_buffer, &help, ipos, IOM_LONG);
    add_bytes(color_chunk_buffer, &help, FLI_256_COLOR, IOM_UWORD);
    add_bytes(color_chunk_buffer, &help, packets, IOM_UWORD);

    return (ipos);
}
