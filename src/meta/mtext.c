#ifndef lint
static const char	RCSid[] = "$Id: mtext.c,v 1.7 2011/08/16 18:09:53 greg Exp $";
#endif
/*
 *  Program to convert ascii file to metafile
 *
 *    6/4/85
 *
 *  cc mtext.c mfio.o syscalls.o misc.o
 */

#include <string.h>

#include  "meta.h"
#include  "plot.h"

#define  MAXLINE  1024

#define  CWIDTH  250

#define  CTHICK  0

#define  CDIR  0

#define  CCOLOR  0

#define  BYASPECT(w)  ((w)*3/2)

static int  cwidth = CWIDTH,
	    cheight = BYASPECT(CWIDTH),
	    cthick = CTHICK,
	    cdir = CDIR,
	    ccolor = CCOLOR;

char  *progname;

static void execute(FILE  *fp);
static void sectout(char  **sect, int  nlines, int  maxlen);
static void plotstr(int  lino, char  *s);



int
main(
	int  argc,
	char  **argv
)
{
    FILE  *fp;
    
    progname = *argv++;
    argc--;

    while (argc && **argv == '-')  {
        switch (*(*argv+1))  {
            case 'w':
            case 'W':
                cwidth = atoi(*argv+2);
		cheight = BYASPECT(cwidth);
		if (cheight < 0 || cheight > XYSIZE)
		    error(USER, "illegal character width");
                break;
	    case 't':
	    case 'T':
		cthick = atoi(*argv+2);
		if (cthick < 0 || cthick > 3)
		    error(USER, "thickness values between 0 and 3 only");
		break;
	    case 'r':
	    case 'R':
	        cdir = 0;
	        break;
	    case 'u':
	    case 'U':
	        cdir = 1;
	        break;
	    case 'l':
	    case 'L':
	        cdir = 2;
	        break;
	    case 'd':
	    case 'D':
	        cdir = 3;
	        break;
            case 'c':
            case 'C':
                ccolor = atoi(*argv+2);
                if (ccolor < 0 || ccolor > 3)
                    error(USER, "color values between 0 and 3 only");
                break;
            default:
                sprintf(errmsg, "unknown option '%s'", *argv);
		error(WARNING, errmsg);
		break;
	}
        argv++;
        argc--;
    }

    if (argc)
        while (argc)  {
	    fp = efopen(*argv, "r");
	    execute(fp);
	    fclose(fp);
	    argv++;
	    argc--;
	}
    else
        execute(stdin);

    pglob(PEOF, 0200, NULL);

return(0);
}


void
execute(		/* execute a file */
	FILE  *fp
)
{
    static char  linbuf[MAXLINE];
    int  nlines;
    char  **section;
    int  maxlen;
    int  done;
    int  j, k;
    
    nlines = XYSIZE/cheight;
    done = FALSE;
    
    if ((section = (char **)calloc(nlines, sizeof(char *))) == NULL)
        error(SYSTEM, "out of memory in execute");
    
    while (!done) {
        maxlen = 0;
        for (j = 0; j < nlines; j++) {
            if ((done = (fgets(linbuf, MAXLINE, fp) == NULL)))
                break;
            k = strlen(linbuf);
            if (linbuf[k-1] == '\n')
                linbuf[--k] = '\0';		/* get rid of newline */
            if (k > maxlen)
                maxlen = k;
            if ((section[j] = malloc(k+1)) == NULL)
                error(SYSTEM, "out of memory in execute");
            strcpy(section[j], linbuf);
        }
	if (maxlen > 0)
	    sectout(section, j, maxlen);
        for (k = 0; k < j; k++) {
            free(section[k]);
            section[k] = NULL;
        }
    }

    free((char *)section);
    
}


void
sectout(		/* write out a section */
	char  **sect,
	int  nlines,
	int  maxlen
)
{
    int  linwidt;
    char  *slin;
    int  i, j;
    
    linwidt = XYSIZE/cwidth;
    
    if ((slin = malloc(linwidt + 1)) == NULL)
        error(SYSTEM, "out of memory in sectout");
        
    for (i = 0; i < maxlen; i += linwidt) {
    
        if (i > 0)
            pglob(PEOP, cdir, NULL);

        for (j = 0; j < nlines; j++)
            if (i < strlen(sect[j])) {
                strncpy(slin, sect[j] + i, linwidt);
                slin[linwidt] = '\0';
                plotstr(j, slin);
	    }

    }

    pglob(PEOP, 0200, NULL);
    free(slin);
    
}


void
plotstr(		/* plot string on line lino */
	int  lino,
	char  *s
)
{
    int  a0;
    register int  bottom, right;
    
    a0 = (cdir<<4) | (cthick<<2) | ccolor;
    bottom = XYSIZE-(lino+1)*cheight;
    right = strlen(s)*cwidth;
    
    switch (cdir) {
        case 0:		/* right */
            pprim(PVSTR, a0, 0, bottom, right, bottom+cheight-1, s);
            break;
        case 1:		/* up */
            pprim(PVSTR, a0, XYSIZE-bottom-cheight+1, 0,
			XYSIZE-bottom, right, s);
            break;
        case 2:		/* left */
            pprim(PVSTR, a0, XYSIZE-right, XYSIZE-bottom-cheight+1,
			XYSIZE-1, XYSIZE-bottom, s);
            break;
        case 3:		/* down */
            pprim(PVSTR, a0, bottom, XYSIZE-right,
			bottom+cheight-1, XYSIZE-1, s);
            break;
    }

}
