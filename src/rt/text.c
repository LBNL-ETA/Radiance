#ifndef lint
static const char	RCSid[] = "$Id: text.c,v 2.23 2003/07/27 22:12:03 schorsch Exp $";
#endif
/*
 *  text.c - functions for text patterns and mixtures.
 */

#include "copyright.h"

#include  "ray.h"

#include  "otypes.h"

#include  "font.h"

/*
 *	A text pattern is specified as the text (a file or line),
 *  the upper left anchor point, the right motion vector, the down
 *  motion vector, and the foreground and background brightness.
 *  For a file, the description is as follows:
 *
 *	modifier brighttext id
 *	2 fontfile textfile
 *	0
 *	11+
 *		Ax Ay Az
 *		Rx Ry Rz
 *		Dx Dy Dz
 *		foreground background
 *		[spacing]
 *
 *  For a single line, we use:
 *
 *	modifier brighttext id
 *	N+2 fontfile . This is a line with N words...
 *	0
 *	11+
 *		Ax Ay Az
 *		Rx Ry Rz
 *		Dx Dy Dz
 *		foreground background
 *		[spacing]
 *
 *  Colortext is identical, except colors are given rather than
 *  brightnesses.
 *
 *  Mixtext has foreground and background modifiers:
 *
 *	modifier mixtext id
 *	4+ foremod backmod fontfile text..
 *	0
 *	9+
 *		Ax Ay Az
 *		Rx Ry Rz
 *		Dx Dy Dz
 *		[spacing]
 */

#define  fndx(m)	((m)->otype==MIX_TEXT ? 2 : 0)
#define  tndx(m)	((m)->otype==MIX_TEXT ? 3 : 1)
#define  sndx(m)	((m)->otype==PAT_BTEXT ? 11 : \
				(m)->otype==PAT_CTEXT ? 15 : 9)

typedef struct tline {
	struct tline  *next;		/* pointer to next line */
	short  *spc;			/* character spacing */
	int  width;			/* total line width */
					/* followed by the string */
}  TLINE;

#define  TLSTR(l)	((char *)((l)+1))

typedef struct {
	FVECT  right, down;		/* right and down unit vectors */
	FONT  *f;			/* our font */
	TLINE  tl;			/* line list */
}  TEXT;

TEXT  *gettext();

TLINE  *tlalloc();


do_text(m, r)
register OBJREC  *m;
RAY  *r;
{
	FVECT  v;
	int  foreground;
				/* get transformed position */
	if (r->rox != NULL)
		multp3(v, r->rop, r->rox->b.xfm);
	else
		VCOPY(v, r->rop);
				/* check if we are within a text glyph */
	foreground = intext(v, m);
				/* modify */
	if (m->otype == MIX_TEXT) {
		OBJECT  omod;
		char  *modname = m->oargs.sarg[foreground ? 0 : 1];
		if (!strcmp(modname, VOIDID))
			omod = OVOID;
		else if ((omod = lastmod(objndx(m), modname)) == OVOID) {
			sprintf(errmsg, "undefined modifier \"%s\"", modname);
			objerror(m, USER, errmsg);
		}
		if (rayshade(r, omod)) {
			if (m->omod != OVOID)
				objerror(m, USER, "inappropriate modifier");
			return(1);
		}
	} else if (m->otype == PAT_BTEXT) {
		if (foreground)
			scalecolor(r->pcol, m->oargs.farg[9]);
		else
			scalecolor(r->pcol, m->oargs.farg[10]);
	} else { /* PAT_CTEXT */
		COLOR  cval;
		if (foreground)
			setcolor(cval, m->oargs.farg[9],
					m->oargs.farg[10],
					m->oargs.farg[11]);
		else
			setcolor(cval, m->oargs.farg[12],
					m->oargs.farg[13],
					m->oargs.farg[14]);
		multcolor(r->pcol, cval);
	}
	return(0);
}


TLINE *
tlalloc(s)			/* allocate and assign text line */
char  *s;
{
	register int  siz;
	register TLINE  *tl;

	siz = strlen(s) + 1;
	if ((tl=(TLINE *)malloc(sizeof(TLINE)+siz)) == NULL ||
			(tl->spc=(short *)malloc(siz*sizeof(short))) == NULL)
		error(SYSTEM, "out of memory in tlalloc");
	tl->next = NULL;
	strcpy(TLSTR(tl), s);
	return(tl);
}


TEXT *
gettext(tm)			/* get text structure for material */
register OBJREC  *tm;
{
#define  R	(tm->oargs.farg+3)
#define  D	(tm->oargs.farg+6)
	FVECT  DxR;
	double  d;
	FILE  *fp;
	char  linbuf[512];
	TEXT  *t;
	register int  i;
	register TLINE  *tlp;
	register char  *s;

	if ((t = (TEXT *)tm->os) != NULL)
		return(t);
						/* check arguments */
	if (tm->oargs.nsargs - tndx(tm) < 1 || tm->oargs.nfargs < sndx(tm))
		objerror(tm, USER, "bad # arguments");
	if ((t = (TEXT *)malloc(sizeof(TEXT))) == NULL)
		error(SYSTEM, "out of memory in gettext");
						/* compute vectors */
	fcross(DxR, D, R);
	fcross(t->right, DxR, D);
	d = DOT(t->right,t->right);
	if (d <= FTINY*FTINY*FTINY*FTINY)
		objerror(tm, USER, "illegal motion vector");
	d = DOT(D,D)/d;
	for (i = 0; i < 3; i++)
		t->right[i] *= d;
	fcross(t->down, R, DxR);
	d = DOT(R,R)/DOT(t->down,t->down);
	for (i = 0; i < 3; i++)
		t->down[i] *= d;
						/* get text */
	tlp = &t->tl;
	if (tm->oargs.nsargs - tndx(tm) > 1) {	/* single line */
		s = linbuf;
		for (i = tndx(tm)+1; i < tm->oargs.nsargs; i++) {
			strcpy(s, tm->oargs.sarg[i]);
			s += strlen(s);
			*s++ = ' ';
		}
		*--s = '\0';
		tlp->next = tlalloc(linbuf);
		tlp = tlp->next;
	} else {				/* text file */
		if ((s = getpath(tm->oargs.sarg[tndx(tm)],
				getrlibpath(), R_OK)) == NULL) {
			sprintf(errmsg, "cannot find text file \"%s\"",
					tm->oargs.sarg[tndx(tm)]);
			error(USER, errmsg);
		}
		if ((fp = fopen(s, "r")) == NULL) {
			sprintf(errmsg, "cannot open text file \"%s\"", s);
			error(SYSTEM, errmsg);
		}
		while (fgets(linbuf, sizeof(linbuf), fp) != NULL) {
			s = linbuf + strlen(linbuf) - 1;
			if (*s == '\n')
				*s = '\0';
			tlp->next = tlalloc(linbuf);
			tlp = tlp->next;
		}
		fclose(fp);
	}
	tlp->next = NULL;
						/* get the font */
	t->f = getfont(tm->oargs.sarg[fndx(tm)]);
						/* compute character spacing */
	i = sndx(tm);
	d = i < tm->oargs.nfargs ? tm->oargs.farg[i] : 0.0;
	i = d * 255.0;
	t->tl.width = 0;
	for (tlp = t->tl.next; tlp != NULL; tlp = tlp->next) {
		if (i < 0)
			tlp->width = squeeztext(tlp->spc, TLSTR(tlp), t->f, -i);
		else if (i > 0)
			tlp->width = proptext(tlp->spc, TLSTR(tlp), t->f, i, 3);
		else
			tlp->width = uniftext(tlp->spc, TLSTR(tlp), t->f);
		if (tlp->width > t->tl.width)
			t->tl.width = tlp->width;
	}
						/* we're done */
	tm->os = (char *)t;
	return(t);
#undef  R
#undef  D
}


freetext(m)			/* free text structures associated with m */
OBJREC  *m;
{
	register TEXT  *tp;
	register TLINE  *tlp;

	tp = (TEXT *)m->os;
	if (tp == NULL)
		return;
	while ((tlp = tp->tl.next) != NULL) {
		tp->tl.next = tlp->next;
		free((void *)tlp->spc);
		free((void *)tlp);
	}
	freefont(tp->f);	/* release font reference */
	free((void *)tp);
	m->os = NULL;
}


intext(p, m)			/* check to see if p is in text glyph */
FVECT  p;
OBJREC  *m;
{
	register TEXT  *tp;
	register TLINE  *tlp;
	FVECT  v;
	double  y, x;
	register int  i, h;
				/* first, compute position in text */
	tp = gettext(m);
	v[0] = p[0] - m->oargs.farg[0];
	v[1] = p[1] - m->oargs.farg[1];
	v[2] = p[2] - m->oargs.farg[2];
	x = DOT(v, tp->right);
	i = sndx(m);
	if (i < m->oargs.nfargs)
		x *= tp->f->mwidth + 255.*fabs(m->oargs.farg[i]);
	else
		x *= 255.;
	h = x;
	i = y = DOT(v, tp->down);
	if (x < 0.0 || y < 0.0)
		return(0);
	x -= (double)h;
	y = ((i+1) - y)*255.;
				/* find the line position */
	for (tlp = tp->tl.next; tlp != NULL; tlp = tlp->next)
		if (--i < 0)
			break;
	if (tlp == NULL || h >= tlp->width)
		return(0);
	for (i = 0; (h -= tlp->spc[i]) >= 0; i++)
		if (h < 255 && inglyph(h+x, y,
				tp->f->fg[TLSTR(tlp)[i]&0xff]))
			return(1);
	return(0);
}


inglyph(x, y, gl)		/* (x,y) within font glyph gl? */
double  x, y;		/* real coordinates in range [0,255) */
register GLYPH  *gl;
{
	int  n, ncross;
	int  xlb, ylb;
	int  tv;
	register GORD  *p0, *p1;

	if (gl == NULL)
		return(0);
	xlb = x;
	ylb = y;
	if (gl->left > xlb || gl->right <= xlb ||	/* check extent */
			gl->bottom > ylb || gl->top <= ylb)
		return(0);
	xlb = xlb<<1 | 1;		/* add 1/2 to test points... */
	ylb = ylb<<1 | 1;		/* ...so no equal comparisons */
	n = gl->nverts;			/* get # of vertices */
	p0 = gvlist(gl) + 2*(n-1);	/* connect last to first */
	p1 = gvlist(gl);
	ncross = 0;
					/* positive x axis cross test */
	while (n--) {
		if ((p0[1]<<1 > ylb) ^ (p1[1]<<1 > ylb)) {
			tv = (p0[0]<<1 > xlb) | ((p1[0]<<1 > xlb) << 1);
			if (tv == 03)
				ncross++;
			else if (tv)
				ncross += (p1[1] > p0[1]) ^
						((p0[1]-y)*(p1[0]-x) >
						(p0[0]-x)*(p1[1]-y));
		}
		p0 = p1;
		p1 += 2;
	}
	return(ncross & 01);
}
