#ifndef lint
static const char	RCSid[] = "$Id: plotin.c,v 1.1 2003/02/22 02:07:26 greg Exp $";
#endif
/*
 *  Program to convert plot(5) files to metafiles
 *
 *  cc -o plotin plotin.c primout.o mfio.o syscalls.o misc.o -lm
 *
 *     5/20/85
 */

#include  "meta.h"


#define  getsi(f, m)  ICONV(geti(f)-(m), size)



char  *progname;

static int  xmin = 0, ymin = 0,		/* current space settings */
	    size =  XYSIZE;

static int  curx = -1,		/* current position */
	    cury = -1;

static short  curmod = 0;	/* current line drawing mode */



main(argc, argv)

int  argc;
char  **argv;

{
 FILE  *fp;

 progname = *argv++;
 argc--;

 while (argc && **argv == '-')  {
    switch (*(*argv+1))  {
       default:
	  error(WARNING, "unknown option");
	  break;
       }
    argv++;
    argc--;
    }

 if (argc)
    while (argc)  {
	  fp = efopen(*argv, "r");
	  convert(fp);
	  fclose(fp);
	  argv++;
	  argc--;
	  }
 else
    convert(stdin);

 pglob(PEOF, 0200, NULL);

 return(0);
 }





convert(infp)		/* convert to meta-file */

FILE  *infp;

{
    char  *fgets(), sbuf[BUFSIZ];
    int  command;
    int  a1, a2, a3, a4, a5, a6;

    while ((command = getc(infp)) != EOF)
	switch (command) {
	    case 'm':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		move(a1, a2);
		break;
	    case 'n':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		cont(a1, a2);
		break;
	    case 'p':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		point(a1, a2);
		break;
	    case 'l':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		a3 = getsi(infp, xmin);
		a4 = getsi(infp, ymin);
		line(a1, a2, a3, a4);
		break;
	    case 't':
		fgets(sbuf, sizeof(sbuf), infp);
		sbuf[strlen(sbuf)-1] = '\0';
		label(sbuf);
		break;
	    case 'a':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		a3 = getsi(infp, xmin);
		a4 = getsi(infp, ymin);
		a5 = getsi(infp, xmin);
		a6 = getsi(infp, ymin);
		arc(a1, a2, a3, a4, a5, a6);
		break;
	    case 'c':
		a1 = getsi(infp, xmin);
		a2 = getsi(infp, ymin);
		a3 = getsi(infp, 0);
		circle(a1, a2, a3);
		break;
	    case 'e':
		erase();
		break;
	    case 'f':
		fgets(sbuf, sizeof(sbuf), infp);
		linemod(sbuf);
		break;
	    case 's':
		a1 = geti(infp);
		a2 = geti(infp);
		a3 = geti(infp);
		a4 = geti(infp);
		space(a1, a2, a3, a4);
		break;
	    default:
		error(USER, "unknown command in convert");
		break;
	}
}




int
geti(fp)		/* get two-byte integer from file */

register FILE  *fp;

{
    register int  i;

    i = getc(fp);
    i |= getc(fp) << 8;
    if (i & 0100000)
	i |= ~0177777;		/* manual sign extend */
    if (feof(fp))
	error(USER, "unexpected end of file in geti");
    
    return(i);
}




move(x, y)			/* move to new position */

register int  x, y;

{

    curx = x;
    cury = y;

}




cont(x, y)			/* draw line to point */

register int  x, y;

{

    plseg(curmod, curx, cury, x, y);
    curx = x;
    cury = y;

}



point(x, y)			/* draw a point */

register int  x, y;

{

    plseg(0, x, y, x, y);
    curx = x;
    cury = y;

}




line(x0, y0, x1, y1)		/* draw a line segment */

int  x0, y0, x1, y1;

{

    move(x0, y0);
    cont(x1, y1);

}



label(s)			/* print a label */

register char  *s;

{

    pprim(PMSTR, 0, curx, cury, curx, cury, s);

}




static int del = 100;
static step(d){
	del = d;
}
arc(x,y,x0,y0,x1,y1){
	double pc;
	double sqrt();
	int flg,m,xc,yc,xs,ys,qs,qf;
	float dx,dy,r;
	char use;
	dx = x-x0;
	dy = y-y0;
	r = dx*dx+dy*dy;
	pc = r;
	pc = sqrt(pc);
	flg = pc/4;
	if(flg == 0)step(1);
	else if(flg < del)step(flg);
	xc = xs = x0;
	yc = ys = y0;
	move(xs,ys);
	if(x0 == x1 && y0 == y1)flg=0;
	else flg=1;
	qs = quadr(x,y,x0,y0);
	qf = quadr(x,y,x1,y1);
	if(abs(x-x1) < abs(y-y1)){
		use = 'x';
		if(qs == 2 || qs ==3)m = -1;
		else m=1;
	}
	else {
		use = 'y';
		if(qs > 2)m= -1;
		else m= 1;
	}
	while(1){
		switch(use){
		case 'x':	
			if(qs == 2 || qs == 3)yc -= del;
			else yc += del;
			dy = yc-y;
			pc = r-dy*dy;
			xc = m*sqrt(pc)+x;
			if((x < xs && x >= xc) || ( x > xs && x <= xc) ||
			    (y < ys && y >= yc) || ( y > ys && y <=  yc) )
			{
				if(++qs > 4)qs=1;
				if(qs == 2 || qs == 3)m= -1;
				else m=1;
				flg=1;
			}
			cont(xc,yc);
			xs = xc; 
			ys = yc;
			if(qs == qf && flg == 1)
				switch(qf){
				case 3:
				case 4:	
					if(xs >= x1)return;
					continue;
				case 1:
				case 2:
					if(xs <= x1)return;
				}
			continue;
		case 'y':	
			if(qs > 2)xc += del;
			else xc -= del;
			dx = xc-x;
			pc = r-dx*dx;
			yc = m*sqrt(pc)+y;
			if((x < xs && x >= xc) || ( x > xs && x <= xc ) ||
			    (y < ys && y >= yc) || (y > ys && y <= yc) )
			{
				if(++qs > 4)qs=1;
				if(qs > 2)m = -1;
				else m = 1;
				flg=1;
			}
			cont(xc,yc);
			xs = xc; 
			ys = yc;
			if(qs == qf && flg == 1)
				switch(qs){
				case 1:
				case 4:
					if(ys >= y1)return;
					continue;
				case 2:
				case 3:
					if(ys <= y1)return;
				}
		}
	}
}
quadr(x,y,xp,yp){
	if(x < xp)
		if(y <= yp)return(1);
		else return(4);
	else if(x > xp)
		if(y < yp)return(2);
		else return(3);
	else if(y < yp)return(2);
	else return(4);
}



circle(x,y,r){
	arc(x,y,x+r,y,x+r,y);
}




erase()				/* erase plot */

{

    pglob(PEOP, 0200, NULL);

}






linemod(s)		/* set line mode according to s */

char  s[];

{

    switch (s[2]) {
	case 'l':		/* solid */
	    curmod = 0;
	    break;
	case 't':		/* dotted */
	    curmod = 1;
	    break;
	case 'n':		/* long dashed */
	    curmod = 2;
	    break;
	case 'o':		/* short dashed (dot dashed) */
	    curmod = 3;
	    break;
	default:
	    error(WARNING, "unknown line mode in linemod");
	    break;
    }

}




space(xmn, ymn, xmx, ymx)	/* change space */

int  xmn, ymn, xmx, ymx;

{

    if (xmn >= xmx || ymn >= ymx)
	error(USER, "illegal space specification in space");
    
    xmin = xmn;
    ymin = ymn;
    size = min(xmx-xmn, ymx-ymn);

}
