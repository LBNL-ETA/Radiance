#ifndef lint
static const char	RCSid[] = "$Id: hfio.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *    Human-readable file I/O
 */


#include  "meta.h"



#if  CPM || MAC
#define  getc  agetc
#define  putc  aputc
#endif


static PRIMITIVE  peof = {PEOF, 0200, -1, -1, -1, -1, NULL};


readp(p, fp)		/* get human-readable primitive */

PRIMITIVE  *p;
FILE  *fp;

{
 char  inbuf[MAXARGS];
 register int  c, nargs;
 int  tmp;

 if (fp == NULL) fp = stdin;

 restart:
 
 if ((c = getc(fp)) == EOF) {		/* used to be fatal */
    mcopy((char *)p, (char *)&peof, sizeof(PRIMITIVE));
    return(0);
 }

 if (c == CDELIM) {			/* skip user comment */
    fgets(inbuf, MAXARGS, fp);
    goto restart;
 } else if (c == '\n')			/* skip empty line */
    goto restart;
    
 if (!iscom(c))
    error(USER, "bad command in readp");

 p->com = c;

 fscanf(fp, "%o", &tmp);
 p->arg0 = tmp & 0377;

 if (isglob(c))
    p->xy[XMN] = p->xy[YMN] = p->xy[XMX] = p->xy[YMX] = -1;
 else if (fscanf(fp, "%d %d %d %d", &p->xy[XMN], &p->xy[YMN],
				&p->xy[XMX], &p->xy[YMX]) != 4)
    error(USER, "missing extent in readp");

 while ((c = getc(fp)) != EOF && c != '\n' && c != ADELIM);

 nargs = 0;

 if (c == ADELIM)
    while ((c = getc(fp)) != EOF && c != '\n' && nargs < MAXARGS-1)
	inbuf[nargs++] = c;

 if (nargs >= MAXARGS)
     error(USER, "too many arguments in readp");

 if (nargs)  {
    inbuf[nargs] = '\0';
    p->args = savestr(inbuf);
    }
 else
    p->args = NULL;

 return(p->com != PEOF);
 }





writep(p, fp)		/* print primitive in human-readable form */

register PRIMITIVE  *p;
FILE  *fp;

{

 if (fp == NULL) fp = stdout;

 if (!iscom(p->com))
     error(USER, "bad command in writep");

 fprintf(fp, "%c ", p->com);
 if (isprim(p->com)) {
    fprintf(fp, "%3o ", p->arg0 & 0177);
    fprintf(fp, "%5d %5d %5d %5d ", p->xy[XMN], p->xy[YMN], p->xy[XMX], p->xy[YMX]);
 } else
    fprintf(fp, "%3o ", p->arg0);

 if (p->args != NULL)  {
    putc(ADELIM, fp);
    fputs(p->args, fp);
    }

 putc('\n', fp);

 if (p->com == PDRAW || p->com == PEOF)
    if (fflush(fp) == -1)
       error(SYSTEM, "error detected writing file in writep");

 }




writeof(fp)		/* write end of file command to fp */

FILE  *fp;

{

 writep(&peof, fp);

 }
