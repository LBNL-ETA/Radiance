#ifndef lint
static const char	RCSid[] = "$Id: imPfuncs.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  imPfuncs - functional interface to imPress, support for imPfuncs.h
 *
 *  Written by William LeFebvre, LCSE, Rice University
 *
 *  This program can be freely redistributed to anyone, with the following
 *  provisions:  that this comment remain intact, and that no money is
 *  charged or collected for that redistribution.
 */

#include <stdio.h>
#include <varargs.h>
#include "imPcodes.h"
#include "imPfuncs.h"

/* commands with byte operands */

/* im_b(code, b) is defined in imPfuncs.h */

/* commands with word operands: */


void
im_w(
	int code,
	unsigned int w
)

{
    putc(code, imout);
    im_putword(w);
}


void
im_www(
	int code,
	unsigned int wa,
	unsigned int wb,
	unsigned int wc
)

{
    putc(code, imout);
    im_putword(wa);
    im_putword(wb);
    im_putword(wc);
}


void
im_wwwww(
	int code,
	unsigned int wa,
	unsigned int wb,
	unsigned int wc,
	unsigned int wd,
	unsigned int we
)

{
    putc(code, imout);
    im_putword(wa);
    im_putword(wb);
    im_putword(wc);
    im_putword(wd);
    im_putword(we);
}

/* commands with bit operands: */

void
im_21(		/* pads 5 bits on left */
	int code,
	unsigned int b2,
	unsigned int b1
)

{
    putc(code, imout);
    im_putbyte((b2 << 1) | b1);
}


void
im_223(		/* pads 1 bit on left */
	int code,
	unsigned int b2a,
	unsigned int b2b,
	unsigned int b3
)

{
    putc(code, imout);
    im_putbyte((b2a << 5) | (b2b << 3) | b3);
}


void
im_77(		/* pads 2 bits on left */
	int code,
	unsigned int b7a,
	unsigned int b7b
)

{
    putc(code, imout);
    im_putword((b7a << 7) | b7b);
}


void
im_putstring(
	char *string
)

{
    fputs(string, imout);
    im_putbyte(0);
}


void
imCreateFamilyTable(va_alist)

va_dcl

{
    va_list vars;
    unsigned int family;
    unsigned int pairs;
    unsigned int map_name;
    char	 *font_name;

    va_start(vars);

    /* first arguments are:  family, pairs */
    family = va_arg(vars, unsigned int);
    pairs  = va_arg(vars, unsigned int);

    /* write the first part of the command */
    putc(imP_CREATE_FAMILY_TABLE, imout);
    im_putbyte(family);
    im_putbyte(pairs);

    /* write each map_name, font_name pair */
    while (pairs-- > 0)
    {
	map_name = va_arg(vars, unsigned int);
	im_putbyte(map_name);
	font_name = va_arg(vars, char *);
	im_putstring(font_name);
    }

    va_end(vars);
}


void
imCreatePath(va_alist)

va_dcl

{
    va_list vars;
    int count;
    int h;
    int v;

    va_start(vars);

    /* first argument is vertex-count */
    count = va_arg(vars, int);

    /* write the first part of the command */
    putc(imP_CREATE_PATH, imout);
    im_putword(count);

    /* write each vertex pair */
    while (count-- > 0)
    {
	h = va_arg(vars, int);
	im_putword(h);
	v = va_arg(vars, int);
	im_putword(v);
    }
}


void
imCreatePathV(
	unsigned int count,
	unsigned int *vec
)

{
    /* write the first part of the command */
    putc(imP_CREATE_PATH, imout);
    im_putword(count);

    /* write each vertex pair */
    while (count-- > 0)
    {
	im_putword(*vec);	/* these are macros... */
	vec++;			/* so im_putword(*vec++) won't always work */
	im_putword(*vec);
	vec++;
    }    
}

#ifdef notyet
/* stuff we still have to do: */
imBitmap(...)
#endif
