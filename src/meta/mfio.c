#ifndef lint
static const char	RCSid[] = "$Id: mfio.c,v 1.2 2003/11/15 02:13:37 schorsch Exp $";
#endif
/*
 *  Meta-file input and output routines
 */

#include  "meta.h"


static PRIMITIVE  peof = {PEOF, 0200, {-1, -1, -1, -1}, NULL, NULL};



#define  INTARG(n)  (((unsigned)inbuf[n] << 7) | (unsigned)inbuf[n+1])


int
readp(		/* read in primitive from file */
	PRIMITIVE  *p,
	FILE  *fp
)

{
 char  inbuf[MAXARGS];
 register int  c, nargs;

 if (fp == NULL) fp = stdin;

 if ((c = getc(fp)) == EOF) {		/* used to be fatal */
    mcopy((char *)p, (char *)&peof, sizeof(PRIMITIVE));
    return(0);
 }
 if (!(c & 0200))
    error(USER, "readp not started on command");

 p->com = c & 0177;
 if (!iscom(p->com))
    error(USER, "bad command in readp");

 for (nargs = 0; (c = getc(fp)) != EOF && !(c & 0200) &&
		 nargs < MAXARGS-1; nargs++)
    inbuf[nargs] = c;

 if (c & 0200)
    ungetc(c, fp);
 else if (nargs >= MAXARGS-1)
    error(USER, "command w/ too many arguments in readp");

 c = 0;

 if (nargs)  {
    p->arg0 = inbuf[c++];
    nargs--;
    }
 else
    p->arg0 = 0200;

 if (isprim(p->com))  {
    if (nargs < 8)
       error(USER, "missing extent in readp");
    p->xy[XMN] = INTARG(c); c += 2;
    p->xy[YMN] = INTARG(c); c += 2;
    p->xy[XMX] = INTARG(c); c += 2;
    p->xy[YMX] = INTARG(c); c += 2;
    nargs -= 8;
    }
 else
    p->xy[XMN] = p->xy[YMN] = p->xy[XMX] = p->xy[YMX] = -1;

 if (nargs)  {
    inbuf[c+nargs] = '\0';
    p->args = savestr(inbuf+c);
    }
 else
    p->args = NULL;

 return(p->com != PEOF);
 }

#undef  INTARG







#define  HI7(i)  ((i >> 7) & 0177)

#define  LO7(i)  (i & 0177)


void
writep(			/* write primitive to file */
	register PRIMITIVE  *p,
	FILE  *fp
)
{
 if (fp == NULL) fp = stdout;

 if (!iscom(p->com))
    error(USER, "bad command in writep");

 putc(p->com | 0200, fp);

 if (isprim(p->com))  {
    putc(p->arg0 & 0177, fp);
    putc(HI7(p->xy[XMN]), fp); putc(LO7(p->xy[XMN]), fp);
    putc(HI7(p->xy[YMN]), fp); putc(LO7(p->xy[YMN]), fp);
    putc(HI7(p->xy[XMX]), fp); putc(LO7(p->xy[XMX]), fp);
    putc(HI7(p->xy[YMX]), fp); putc(LO7(p->xy[YMX]), fp);
    }
 else if (!(p->arg0 & 0200))
    putc(p->arg0, fp);

 if (p->args != NULL)
    fputs(p->args, fp);

 if (p->com == PDRAW || p->com == PEOF)
    if (fflush(fp) == -1)
       error(SYSTEM, "error detected writing file in writep");

 }

#undef  HI7

#undef  LO7


void
writeof(		/* write end of file command to fp */
	FILE  *fp
)
{
 writep(&peof, fp);
}

