#ifndef lint
static const char	RCSid[] = "$Id: readfargs.c,v 2.7 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 * Allocate, read and free object arguments
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

#include "object.h"



int
readfargs(fa, fp)		/* read function arguments from stream */
register FUNARGS  *fa;
FILE  *fp;
{
#define getstr(s)	(fgetword(s,sizeof(s),fp)!=NULL)
#define getint(s)	(getstr(s) && isint(s))
#define getflt(s)	(getstr(s) && isflt(s))
	char  sbuf[MAXSTR];
	register int  n, i;

	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
	if (fa->nsargs = n) {
		fa->sarg = (char **)malloc(n*sizeof(char *));
		if (fa->sarg == NULL)
			return(-1);
		for (i = 0; i < fa->nsargs; i++) {
			if (!getstr(sbuf))
				return(0);
			fa->sarg[i] = savestr(sbuf);
		}
	} else
		fa->sarg = NULL;
	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
#ifdef  IARGS
	if (fa->niargs = n) {
		fa->iarg = (long *)malloc(n*sizeof(long));
		if (fa->iarg == NULL)
			return(-1);
		for (i = 0; i < n; i++) {
			if (!getint(sbuf))
				return(0);
			fa->iarg[i] = atol(sbuf);
		}
	} else
		fa->iarg = NULL;
#else
	if (n != 0)
		return(0);
#endif
	if (!getint(sbuf) || (n = atoi(sbuf)) < 0)
		return(0);
	if (fa->nfargs = n) {
		fa->farg = (FLOAT *)malloc(n*sizeof(FLOAT));
		if (fa->farg == NULL)
			return(-1);
		for (i = 0; i < n; i++) {
			if (!getflt(sbuf))
				return(0);
			fa->farg[i] = atof(sbuf);
		}
	} else
		fa->farg = NULL;
	return(1);
#undef getflt
#undef getint
#undef getstr
}


void
freefargs(fa)			/* free object arguments */
register FUNARGS  *fa;
{
	register int  i;

	if (fa->nsargs) {
		for (i = 0; i < fa->nsargs; i++)
			freestr(fa->sarg[i]);
		free((void *)fa->sarg);
		fa->sarg = NULL;
		fa->nsargs = 0;
	}
#ifdef  IARGS
	if (fa->niargs) {
		free((void *)fa->iarg);
		fa->iarg = NULL;
		fa->niargs = 0;
	}
#endif
	if (fa->nfargs) {
		free((void *)fa->farg);
		fa->farg = NULL;
		fa->nfargs = 0;
	}
}
