/*
 * ACE/vis - Visualization of Flow and Transport
 *
 * Paul J. Turner and Antonio M. Baptista
 *
 * Copyright 1990-2003 Oregon Health and Science University
 * All Rights Reserved
 *
 */

/*
 * stubs
 */

#ifdef WIN32

void drand48()
{
}

void lrand48()
{
}
void srand48()
{
}
void gamma()
{
}
void lgamma()
{
}
void erf()
{
}
void erfc()
{
}
void sleep()
{
}
void gettxt()
{
}
void pclose()
{
}
void popen()
{
}
#endif

void activeset(void)
{
}

void update_draw(void)
{
}

void getdata(void)
{
}

void setdrogscolor(void)
{
}

void set_left_footer(void)
{
}

void readboundary(void)
{
}

void display_textfile(void)
{
}

void init_drog_colors(void)
{
}

void readpar(void)
{
}

void set_right_footer(void)
{
}

void draw_timeinfo(void)
{
}

void data_isactive(void)
{
}

void update_fuzz_items(void)
{
}

void do_drawgrid(void)
{
}

int celevmin;
int cboundmin;

void box(void)
{
}

void my_circle(double x, double y)
{
    drawpolysym(&x, &y, 1, 2, 0, 0, 0.8);
}

void put_text(void)
{
}
