#ifndef lint
static const char	RCSid[] = "$Id: tgraph.c,v 1.2 2003/06/30 14:59:12 schorsch Exp $";
#endif
/*
 *   Routines for tel-a-graph file plotting
 *
 *   1/20/85
 */
 
#include  <ctype.h>
#include  <string.h>

#include  "tgraph.h"


#define  isfloat(a)  (isdigit(a) || (a) == '-' || (a) == '.' ||  \
                     (a) == 'E' || (a) == '+' || (a) == 'e')


extern double  log10();





initialize()

{
 int  i;

 for (i = 0; i < NCUR; i++)
    usecurve[i] = 1;

 }






option(s)		/* record option */

char  *s;

{
 double  atof();
 char  *sp;
 short  userror = 0;
 int  i;

 if (s[0] == '-')

    switch(s[1])  {
    
       case 'x':
	  xmnset = atof(s+2);
          break;

       case 'y':
	  ymnset = atof(s+2);
          break;

       case 's':	/* set symbol size */
          symrad = atoi(s+2);
          break;
       
       case 'l':	/* log plot */
	  if (s[2] == 'x')
	     logx++;
	  else if (s[2] == 'y')
	     logy++;
	  else
	     userror++;
	  break;

       case 'g':	/* grid on */
          grid = TRUE;
          break;

       case 'p':	/* polar coordinates */
	  if (s[2] == 'd')
	     polar = DEGREES;
	  else if (s[2] == 'r')
	     polar = RADIANS;
	  else
	     userror++;
	  break;
          
       default:
          for (sp = s+1; *sp && *sp > '@' && *sp < '@'+NCUR; sp++)
             usecurve[*sp-'@'] = 0;
          userror = *sp;
          break;
       }
 else

    switch (s[1])  {

       case 'x':
	  xmxset = atof(s+2);
          break;

       case 'y':
	  ymxset = atof(s+2);
          break;

       default:
          for (i = 0; i < 27; i++)
             usecurve[i] = 0;
          for (sp = s+1; *sp && *sp > '@' && *sp < '@'+NCUR; sp++)
             usecurve[*sp-'@'] = 1;
          userror = *sp;
       }

 if (userror)
    error(USER, "options are [-sSYMRAD][-g][-lx][-ly][-p{dr}][-xXMIN][+xXMAX][-yYMIN][+yYMAX][-C..|+C..]");

 }




normalize(fp, fout)			/* get extrema from file */

FILE  *fp, *fout;

{
 char  line[255];
 double  x, y;
 long  npoints = 0;

 ncurves = 0;
 xmin = ymin = FHUGE;
 xmax = ymax = -FHUGE;

 while (fgets(line, sizeof line, fp) != NULL)  {

    if (fout != NULL)
       fputs(line, fout);

    if (islabel(line))

       ncurves++;

    else if (ncurves < NCUR && usecurve[ncurves] && isdata(line))  {

       if (getdata(line, &x, &y) < 0)
	  continue;

       if (x < xmin)
          xmin = x;
       if (x > xmax)
          xmax = x;
       if (y < ymin)
          ymin = y;
       if (y > ymax)
          ymax = y;
          
       npoints++;

       }

    }

 if (npoints <= 1)
    error(USER, "insufficient data in file");

 }





makeaxis(flag)		/* make and print x and y axis */

int  flag;

{
 double  xstep, ystep, step(), pos;
 int  xorg, yorg;
 char	*format, stemp[32];

					/* set limits */
 if (polar)
    if (xmax-xmin < ymax-ymin)			/* null aspect for polar */
	xmax = xmin + ymax-ymin;
    else
	ymax = ymin + xmax-xmin;
 else {
    if (xmnset > -FHUGE)
       if (logx) {
	  if (xmnset > FTINY)
	     xmin = log10(xmnset);
	  }
       else
	  xmin = xmnset;
    if (xmxset < FHUGE)
       if (logx) {
	  if (xmxset > FTINY)
	     xmax = log10(xmxset);
	  }
       else
	  xmax = xmxset;
    if (ymnset > -FHUGE)
       if (logy) {
	  if (ymnset > FTINY)
	     ymin = log10(ymnset);
	  }
       else
	  ymin = ymnset;
    if (ymxset < FHUGE)
       if (logy) {
	  if (ymxset > FTINY)
	     ymax = log10(ymxset);
	  }
       else
	  ymax = ymxset;
    }
					/* set step */
 if (logx) {
    xmin = floor(xmin);
    xmax = ceil(xmax);
    xstep = 1.0;
 } else
    xstep = step(&xmin, &xmax);

 if (logy) {
    ymin = floor(ymin);
    ymax = ceil(ymax);
    ystep = 1.0;
 } else
    ystep = step(&ymin, &ymax);

 xsize = xmax - xmin;
 ysize = ymax - ymin;

 if (polar)				/* null aspect again */
     if (xsize < ysize)
	xmax = xmin + (xsize = ysize);
     else
	ymax = ymin + (ysize = xsize);

 if (xmin > 0.0)
    xorg = XBEG;
 else if (xmax < 0.0)
    xorg = XBEG+XSIZ-1;
 else
    xorg = XCONV(0.0);
 if (ymin > 0.0)
    yorg = YBEG;
 else if (ymax < 0.0)
    yorg = YBEG+YSIZ-1;
 else
    yorg = YCONV(0.0);
					/* make x-axis */
 if (flag & BOX) {
    plseg(010, XBEG, YBEG, XBEG+XSIZ-1, YBEG);
    plseg(010, XBEG, YBEG+YSIZ-1, XBEG+XSIZ-1, YBEG+YSIZ-1);
    }
 if (flag & ORIGIN)
    plseg(010, XBEG, yorg, XBEG+XSIZ-1, yorg);

#if  CPM || MAC
 if (xstep > 1.0)
    format = "%5.0f";
 else
    format = "%5.3f";
#else
 format = "%5g";
#endif
 if (logx)
    format = "1e%-3.0f";

 for (pos = xmin; pos < xmax+xstep/2; pos += xstep) {
    if (flag & XTICS)
       if (polar)
	  plseg(010, XCONV(pos), yorg-TSIZ/2, XCONV(pos), yorg+TSIZ/2);
       else {
	  plseg(010, XCONV(pos), YBEG, XCONV(pos), YBEG-TSIZ);
	  plseg(010, XCONV(pos), YBEG+YSIZ-1, XCONV(pos), YBEG+YSIZ-1+TSIZ);
	  }
    if (flag & XNUMS) {
       sprintf(stemp, format, pos);
       if (polar)
          pprim(PMSTR, 020, XCONV(pos)-600, yorg-TSIZ-150,
			XCONV(pos)-600, yorg-TSIZ-150, stemp);
       else
          pprim(PMSTR, 020, XCONV(pos)-600, YBEG-2*TSIZ-150,
			XCONV(pos)-600, YBEG-2*TSIZ-150, stemp);
       }
    if (flag & XGRID)
       plseg(040, XCONV(pos), YBEG, XCONV(pos), YBEG+YSIZ-1);
    }
					/* make y-axis */
 if (flag & BOX) {
    plseg(010, XBEG, YBEG, XBEG, YBEG+YSIZ-1);
    plseg(010, XBEG+XSIZ-1, YBEG, XBEG+XSIZ-1, YBEG+YSIZ-1);
    }
 if (flag & ORIGIN)
    plseg(010, xorg, YBEG, xorg, YBEG+YSIZ-1);

#if  CPM || MAC
 if (ystep > 1.0)
    format = "%5.0f";
 else
    format = "%5.3f";
#else
 format = "%5g";
#endif
 if (logy)
    format = "1e%-3.0f";

 for (pos = ymin; pos < ymax+ystep/2; pos += ystep) {
    if (flag & YTICS)
       if (polar)
	  plseg(010, xorg-TSIZ/2, YCONV(pos), xorg+TSIZ/2, YCONV(pos));
       else {
          plseg(010, XBEG, YCONV(pos), XBEG-TSIZ, YCONV(pos));
          plseg(010, XBEG+XSIZ-1, YCONV(pos), XBEG+XSIZ-1+TSIZ, YCONV(pos));
          }
    if (flag & YNUMS) {
       sprintf(stemp, format, pos);
       if (polar)
          pprim(PMSTR, 020, xorg-TSIZ-900, YCONV(pos)+32,
			xorg-TSIZ-900, YCONV(pos)+32, stemp);
       else
          pprim(PMSTR, 020, XBEG-2*TSIZ-900, YCONV(pos)+32,
			XBEG-2*TSIZ-900, YCONV(pos)+32, stemp);
       }
    if (flag & YGRID)
       plseg(040, XBEG, YCONV(pos), XBEG+XSIZ-1, YCONV(pos));
    }

}



isdata(s)

register char  *s;

{
 int  commas = 0;

 while (*s)
    if (isfloat(*s) || isspace(*s))
       s++;
    else if (*s == ',')  {
       commas++;
       s++;
       }
    else
       return 0;

 return commas <= 1;
 }




islabel(s)

char  s[];

{
 int  i;

 i = strlen(s) - 2;

 return(i > 0 && s[0] == '"' && s[i] == '"');
 }




double
step(mn, mx)			/* compute step size for axis */

double	*mn, *mx;

{
    static int	steps[] = {100, 50, 20, 10, 5, 2, 1};
    int		i;
    double	fact, stp;
    double	pown(), floor(), ceil();

    if (*mx-*mn <= FTINY) {
       stp = 1.0;
       }
    else {
	fact = pown(10.0, (int)log10(*mx-*mn)-2);
	stp = (*mx-*mn)/fact/MINDIVS;

	for (i = 0; stp < steps[i]; i++);
	stp = steps[i]*fact;
	}

    *mn = floor(*mn/stp) * stp;
    *mx = ceil(*mx/stp) * stp;

    return(stp);
}



double
pown(x, n)		/* raise x to an integer power */

double	x;
int	n;

{
    register int	i;
    double		p = 1.0;

    if (n > 0)
	for (i = 0; i < n; i++)
	    p *= x;
    else
	for (i = 0; i > n; i--)
	    p /= x;
    
    return(p);
}




istitle(s)

char  *s;

{
 char  word[32];

 if (sscanf(s, "%10s", word) == 1)
    return strcmp(word, "title") == 0;
 else
    return 0;
 }




isdivlab(s)			/* return TRUE if division label(s) */

register char  *s;

{

 return(instr(s, "division") != NULL);
}




isxlabel(s)

register char  *s;

{
 register char  *xindex = instr(s, "x ");

 return(xindex != NULL && instr(xindex, "label ") != NULL);
 }




isylabel(s)

register char  *s;

{
 register char  *yindex = instr(s, "y ");

 return(yindex != NULL && instr(yindex, "label ") != NULL);
 }



char *
instr(s, t)		/* return pointer to first occurrence of t in s */

char  *s,
      *t;

{
 register char  *pt, *ps;

 do  {

    ps = s;
    pt = t;

    while (*pt && *pt == *ps++)
       pt++;

    if (*pt == '\0')
       return(s);

    }  while (*s++);

 return(NULL);
 }




char *
snagquo(s)		/* find and return quoted string within s */

register char  *s;

{
    register char  *rval = NULL;

    for ( ; *s; s++)
	if (*s == '"')
	    if (rval == NULL)
		rval = s+1;
	    else {
		*s = '\0';
		return(rval);
	    }
    
    return(NULL);
}




getdata(s, xp, yp)			/* get data from line */

char  *s;
double  *xp, *yp;

{
    double  sin(), cos();
    int  oobounds = 0;
    double  a;
    register char  *cp;
    
    if ((cp = strchr(s, ',')) != NULL)
	*cp = ' ';

    if (sscanf(s, "%lf %lf", xp, yp) != 2)
	return(-1);

    if (*xp < xmnset || *xp > xmxset)
	oobounds++;
    if (*yp < ymnset || *yp > ymxset)
	oobounds++;

    if (logx)
	if (*xp < FTINY)
	    oobounds++;
	else
	    *xp = log10(*xp);

    if (logy)
	if (*yp < FTINY)
	    oobounds++;
	else
	    *yp = log10(*yp); 

    if (polar) {
	a = *xp;
	if (polar == DEGREES)
	    a *= PI/180.0;
	*xp = *yp * cos(a);
	*yp *= sin(a);
    }

    return(oobounds ? -1 : 0);
}




symout(a0, x, y, sname)			/* output a symbol */

int  a0;
int  x, y;
char  *sname;

{

    pprim(PSEG, a0, x-symrad, y-symrad, x+symrad, y+symrad, sname);
    
}



boxstring(a0, xmn, ymn, xmx, ymx, s)	/* output a string within a box */

int  a0;
int  xmn, ymn, xmx, ymx;
char  *s;

{
    int  start;
    long  size;

    if (a0 & 020) {			/* up or down */

	size = (long)strlen(s)*(xmx-xmn)*2/3;		/* aspect ratio is 1.5 */

	if (size > ymx-ymn) {
	    size = ymx-ymn;				/* squash it */
	    start = ymn;
	} else
	    start = (ymx-ymn-(int)size)/2 + ymn;	/* center it */
	
	pprim(PVSTR, a0, xmn, start, xmx, start+(int)size, s);

    } else {				/* right or left */

	size = (long)strlen(s)*(ymx-ymn)*2/3;		/* aspect ratio is 1.5 */

	if (size > xmx-xmn) {
	    size = xmx-xmn;				/* squash it */
	    start = xmn;
	} else
	    start = (xmx-xmn-(int)size)/2 + xmn;	/* center it */
	
	pprim(PVSTR, a0, start, ymn, start+(int)size, ymx, s);

    }

}
