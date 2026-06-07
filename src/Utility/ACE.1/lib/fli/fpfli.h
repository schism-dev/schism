/* fpfli.h */

/******
  Copyright (C) 1993 by Klaus Ehrenfried. 

  Permission to use, copy, modify, and distribute this software
  is hereby granted, provided that the above copyright notice appears 
  in all copies and that the software is available to all free of charge. 
  The author disclaims all warranties with regard to this software, 
  including all implied warranties of merchant-ability and fitness. 
  The code is simply distributed as it is.
*******/

/*
#define UBYTE unsigned char
#define LONG long
##define SHORT short
##define USHORT unsigned short
##define ULONG unsigned long
*/

typedef unsigned char UBYTE;
typedef short SHORT;
typedef unsigned short USHORT;
typedef long LONG;
typedef unsigned long ULONG;

#define FLI_MAX_X  1280
#define FLI_MAX_Y  1024
#define FLI_MAX_COLORS 256
#define FLI_MAX_FRAMES 4000
#define FLI_FILE_MAGIC 0xaf12		/* File header Magic */
#define FLI_FRAME_MAGIC 0xf1fa		/* Frame Magic */

#define FLI_FILE_HEADER_SIZE 128
#define FLI_FRAME_HEADER_SIZE 16

/* types of chunk in a fli_frame */
#define FLI_256_COLOR 4
#define FLI_DELTA 7
#define FLI_COLOR 11
#define FLI_LC	12
#define FLI_BLACK 13
#define FLI_BRUN 15
#define FLI_COPY 16

#define IOM_SBYTE  1
#define IOM_UBYTE  2
#define IOM_SWORD  3
#define IOM_UWORD  4
#define IOM_LONG  5

#define MAP_FIRST_FRAME 1
#define MAP_NEXT_FRAME 2
#define MAP_CLOSE_LOOP 3

#ifdef MAIN
#define EXT
#else
#define EXT extern
#endif

/* #define BORDER_COLOR 0xFF */

/* external variables */

EXT UBYTE *big_buffer;
EXT UBYTE *pixel_chunk_buffer;
EXT UBYTE color_chunk_buffer[3 * FLI_MAX_COLORS + 10];
EXT int fli_width, fli_height, fli_size, fli_speed;
EXT int border_color, double_buffer;
EXT int Xorigin, Yorigin, Xorigin_flag, Yorigin_flag;
EXT LONG map_color[FLI_MAX_COLORS], map_color_flag;
EXT FILE *input, *output;

/* prototypes */

int exitialise(int);

int get_image(UBYTE *data, LONG color[], int without_data);

int make_fli();

int fli_write_frame(UBYTE *prepre_pixel,
	UBYTE *pre_pixel,
	UBYTE *curr_pixel,
	LONG curr_color[],
	int first_flag);

void add_bytes(UBYTE record[], int *ipos, int value, int mode);

int make_256_color_chunk(LONG color[], int first_flag);

int make_brun_chunk(UBYTE *image);

int make_delta_chunk(unsigned char *preprevious,
	unsigned char *previous,
	unsigned char *current);
