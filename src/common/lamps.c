#ifndef lint
static const char	RCSid[] = "$Id: lamps.c,v 2.6 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Load lamp data.
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

#include  <stdio.h>
#include  <stdlib.h>
#include  <ctype.h>
#include  "color.h"

extern char	*eindex(), *expsave();
extern FILE	*frlibopen();

typedef struct lamp {
	char	*pattern;			/* search pattern */
	float	*color;				/* pointer to lamp value */
	struct lamp	*next;			/* next lamp in list */
} LAMP;					/* a lamp entry */

static LAMP	*lamps = NULL;		/* lamp list */


float *
matchlamp(s)				/* see if string matches any lamp */
char	*s;
{
	register LAMP	*lp;

	for (lp = lamps; lp != NULL; lp = lp->next) {
		expset(lp->pattern);
		if (eindex(s) != NULL)
			return(lp->color);
	}
	return(NULL);
}


loadlamps(file)					/* load lamp type file */
char	*file;
{
	LAMP	*lastp;
	register LAMP	*lp;
	FILE	*fp;
	float	xyz[3];
	char	buf[128], str[128];
	register char	*cp1, *cp2;

	if ((fp = frlibopen(file)) == NULL)
		return(0);
	lastp = NULL;
	while (fgets(buf, sizeof(buf), fp) != NULL) {
					/* work on a copy of buffer */
		strcpy(str, buf);
					/* get pattern for this entry */
		for (cp1 = str; isspace(*cp1); cp1++)
			;
		if (!*cp1 || *cp1 == '#')
			continue;
		for (cp2 = cp1+1; *cp2; cp2++)
			if (*cp2 == *cp1 && cp2[-1] != '\\')
				break;
		if (!*cp2) {
			cp1 = "unclosed pattern";
			goto fmterr;
		}
		cp1++; *cp2++ = '\0';
		if (ecompile(cp1, 1, 0) < 0) {
			cp1 = "bad regular expression";
			goto fmterr;
		}
		if ((lp = (LAMP *)malloc(sizeof(LAMP))) == NULL)
			goto memerr;
		if ((lp->pattern = expsave()) == NULL)
			goto memerr;
						/* get specification */
		for (cp1 = cp2; isspace(*cp1); cp1++)
			;
		if (!isdigit(*cp1) && *cp1 != '.' && *cp1 != '(') {
			cp1 = "missing lamp specification";
			goto fmterr;
		}
		if (*cp1 == '(') {		/* find alias */
			for (cp2 = ++cp1; *cp2; cp2++)
				if (*cp2 == ')' && cp2[-1] != '\\')
					break;
			*cp2 = '\0';
			if ((lp->color = matchlamp(cp1)) == NULL) {
				cp1 = "unmatched alias";
				goto fmterr;
			}
		} else {			/* or read specificaion */
			if ((lp->color=(float *)malloc(6*sizeof(float)))==NULL)
				goto memerr;
			if (sscanf(cp1, "%f %f %f", &lp->color[3], &
					lp->color[4], &lp->color[5]) != 3) {
				cp1 = "bad lamp data";
				goto fmterr;
			}
						/* convert xyY to XYZ */
			xyz[1] = lp->color[5];
			xyz[0] = lp->color[3]/lp->color[4] * xyz[1];
			xyz[2] = xyz[1]*(1./lp->color[4] - 1.) - xyz[0];
						/* XYZ to RGB */
			cie_rgb(lp->color, xyz);
		}
		if (lastp == NULL)
			lamps = lp;
		else
			lastp->next = lp;
		lp->next = NULL;
		lastp = lp;
	}
	fclose(fp);
	return(lastp != NULL);
memerr:
	fputs("Out of memory in loadlamps\n", stderr);
	return(-1);
fmterr:
	fputs(buf, stderr);
	fprintf(stderr, "%s: %s\n", file, cp1);
	return(-1);
}


freelamps()			/* free our lamps list */
{
	register LAMP	*lp1, *lp2;
	
	for (lp1 = lamps; lp1 != NULL; lp1 = lp2) {
		free(lp1->pattern);
		if (lp1->color != NULL) {
			for (lp2 = lp1->next; lp2 != NULL; lp2 = lp2->next)
				if (lp2->color == lp1->color)
					lp2->color = NULL;
			free((void *)lp1->color);
		}
		lp2 = lp1->next;
		free((void *)lp1);
	}
	lamps = NULL;
}
