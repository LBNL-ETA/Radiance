#ifndef lint
static const char	RCSid[] = "$Id: font.c,v 2.12 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Polygonal font handling routines
 */

/* ====================================================================
 * The Radiance Software License, Version 1.0
 *
 * Copyright (c) 1990 - 2002 The Regents of the University of California,
 * through Lawrence Berkeley National Laboratory.   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *         notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in
 *       the documentation and/or other materials provided with the
 *       distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *           if any, must include the following acknowledgment:
 *             "This product includes Radiance software
 *                 (http://radsite.lbl.gov/)
 *                 developed by the Lawrence Berkeley National Laboratory
 *               (http://www.lbl.gov/)."
 *       Alternately, this acknowledgment may appear in the software itself,
 *       if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names "Radiance," "Lawrence Berkeley National Laboratory"
 *       and "The Regents of the University of California" must
 *       not be used to endorse or promote products derived from this
 *       software without prior written permission. For written
 *       permission, please contact radiance@radsite.lbl.gov.
 *
 * 5. Products derived from this software may not be called "Radiance",
 *       nor may "Radiance" appear in their name, without prior written
 *       permission of Lawrence Berkeley National Laboratory.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.   IN NO EVENT SHALL Lawrence Berkeley National Laboratory OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 * ====================================================================
 *
 * This software consists of voluntary contributions made by many
 * individuals on behalf of Lawrence Berkeley National Laboratory.   For more
 * information on Lawrence Berkeley National Laboratory, please see
 * <http://www.lbl.gov/>.
 */

#include "standard.h"

#include "font.h"

#define galloc(nv)	(GLYPH *)malloc(sizeof(GLYPH)+2*sizeof(GORD)*(nv))


int	retainfonts = 0;		/* retain loaded fonts? */

static FONT	*fontlist = NULL;	/* list of loaded fonts */


FONT *
getfont(fname)				/* return font fname */
char  *fname;
{
	FILE  *fp;
	char  *pathname, *err;
	unsigned  wsum, hsum, ngly;
	int  gn, ngv, gv;
	register GLYPH	*g;
	GORD  *gp;
	register FONT  *f;

	for (f = fontlist; f != NULL; f = f->next)
		if (!strcmp(f->name, fname)) {
			f->nref++;
			return(f);
		}
						/* load the font file */
	if ((pathname = getpath(fname, getlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find font file \"%s\"", fname);
		error(USER, errmsg);
	}
	f = (FONT *)calloc(1, sizeof(FONT));
	if (f == NULL)
		goto memerr;
	f->name = savestr(fname);
	f->nref = 1;
	if ((fp = fopen(pathname, "r")) == NULL) {
		sprintf(errmsg, "cannot open font file \"%s\"",
				pathname);
		error(SYSTEM, errmsg);
	}
	wsum = hsum = ngly = 0;			/* get each glyph */
	while ((ngv = fgetval(fp, 'i', (char *)&gn)) != EOF) {
		if (ngv == 0)
			goto nonint;
		if (gn < 1 || gn > 255) {
			err = "illegal";
			goto fonterr;
		}
		if (f->fg[gn] != NULL) {
			err = "duplicate";
			goto fonterr;
		}
		if (fgetval(fp, 'i', (char *)&ngv) <= 0 ||
				ngv < 0 || ngv > 32000) {
			err = "bad # vertices for";
			goto fonterr;
		}
		g = galloc(ngv);
		if (g == NULL)
			goto memerr;
		g->nverts = ngv;
		g->left = g->right = g->top = g->bottom = 128;
		ngv *= 2;
		gp = gvlist(g);
		while (ngv--) {
			if (fgetval(fp, 'i', (char *)&gv) <= 0 ||
					gv < 0 || gv > 255) {
				err = "bad vertex for";
				goto fonterr;
			}
			*gp++ = gv;
			if (ngv & 1) {		/* follow x limits */
				if (gv < g->left)
					g->left = gv;
				else if (gv > g->right)
					g->right = gv;
			} else {		/* follow y limits */
				if (gv < g->bottom)
					g->bottom = gv;
				else if (gv > g->top)
					g->top = gv;
			}
		}
		if (g->right - g->left && g->top - g->bottom) {
			ngly++;
			wsum += g->right - g->left;
			hsum += g->top - g->bottom;
		}
		f->fg[gn] = g;
	}
	fclose(fp);
	if (ngly) {
		f->mwidth = wsum / ngly;
		f->mheight = hsum / ngly;
	}
	f->next = fontlist;
	return(fontlist = f);
nonint:
	sprintf(errmsg, "non-integer in font file \"%s\"", pathname);
	error(USER, errmsg);
fonterr:
	sprintf(errmsg, "%s character (%d) in font file \"%s\"",
			err, gn, pathname);
	error(USER, errmsg);
memerr:
	error(SYSTEM, "out of memory in fontglyph");
}


void
freefont(fnt)			/* release a font (free all if NULL) */
FONT *fnt;
{
	FONT  head;
	register FONT  *fl, *f;
	register int  i;
					/* check reference count */
	if (fnt != NULL && (fnt->nref-- > 1 | retainfonts))
		return;
	head.next = fontlist;
	fl = &head;
	while ((f = fl->next) != NULL)
		if ((fnt == NULL | fnt == f)) {
			fl->next = f->next;
			for (i = 0; i < 256; i++)
				if (f->fg[i] != NULL)
					free((void *)f->fg[i]);
			freestr(f->name);
			free((void *)f);
		} else
			fl = f;
	fontlist = head.next;
}


int
uniftext(sp, tp, f)			/* uniformly space text line */
register short	*sp;		/* returned character spacing */
register char  *tp;		/* text line */
register FONT  *f;		/* font */
{
	int  linelen;

	linelen = *sp++ = 0;
	while (*tp)
		if (f->fg[*tp++&0xff] == NULL)
			*sp++ = 0;
		else
			linelen += *sp++ = 255;
	return(linelen);
}


int
squeeztext(sp, tp, f, cis)		/* squeeze text line */
short  *sp;			/* returned character spacing */
char  *tp;			/* text line */
FONT  *f;			/* font */
int  cis;			/* intercharacter spacing */
{
	int  linelen;
	register GLYPH	*gp;

	linelen = 0;
	gp = NULL;
	while (*tp && (gp = f->fg[*tp++&0xff]) == NULL)
		*sp++ = 0;
	cis /= 2;
	*sp = cis;
	while (gp != NULL) {
		if (gp->nverts) {		/* regular character */
			linelen += *sp++ += cis - gp->left;
			*sp = gp->right + cis;
		} else {			/* space */
			linelen += *sp++;
			*sp = f->mwidth;
		}
		gp = NULL;
		while (*tp && (gp = f->fg[*tp++&0xff]) == NULL) {
			linelen += *sp++;
			*sp = 0;
		}
	}
	linelen += *sp += cis;
	return(linelen);
}


int
proptext(sp, tp, f, cis, nsi)		/* space line proportionally */
short  *sp;			/* returned character spacing */
char  *tp;			/* text line */
FONT  *f;			/* font */
int  cis;			/* target intercharacter spacing */
int  nsi;			/* minimum number of spaces for indent */
{
	register char  *end, *tab;
	GLYPH  *gp;
	short  *nsp;
	int  alen, len, width;
					/* start by squeezing it */
	squeeztext(sp, tp, f, cis);
					/* now, realign spacing */
	width = *sp++;
	while (*tp) {
		len = alen = 0;
		nsp = sp;
		for (end = tp; *end; end = tab) {
			tab = end + 1;
			alen += *nsp++;
			if (f->fg[*end&0xff]) {
				while ((gp = f->fg[*tab&0xff]) != NULL &&
						gp->nverts == 0) { /* tab in */
					alen += *nsp++;
					tab++;
				}
				len += tab - end;
			}
			if (nsi && tab - end > nsi)
				break;
		}
		len *= f->mwidth + cis;		/* compute target length */
		width += len;
		len -= alen;			/* necessary adjustment */
		while (sp < nsp) {
			alen = len/(nsp-sp);
			*sp++ += alen;
			len -= alen;
		}
		tp = tab;
	}
	return(width);
}
