#ifndef lint
static const char	RCSid[] = "$Id: savestr.c,v 2.5 2003/02/22 02:07:22 greg Exp $";
#endif
/*
 *  savestr.c - routines for efficient string storage.
 *
 *	Savestr(s) stores a shared read-only string.
 *  Freestr(s) indicates a client is finished with a string.
 *	All strings must be null-terminated.  There is
 *  no imposed length limit.
 *	Strings stored with savestr(s) can be equated
 *  reliably using their pointer values.
 *	Calls to savestr(s) and freestr(s) should be
 *  balanced (obviously).  The last call to freestr(s)
 *  frees memory associated with the string; it should
 *  never be referenced again.
 *
 *  External symbols declared in standard.h
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

#ifndef  NHASH
#define  NHASH		509		/* hash table size (prime!) */
#endif

typedef struct s_head {
	struct s_head  *next;		/* next in hash list */
	int  nl;			/* links count */
}  S_HEAD;				/* followed by the string itself */

static S_HEAD  *stab[NHASH];

#define  hash(s)	(shash(s)%NHASH)

extern char  *savestr(), *strcpy(), *malloc();

#define  NULL		0

#define  string(sp)	((char *)((sp)+1))

#define  salloc(str)	(S_HEAD *)malloc(sizeof(S_HEAD)+1+strlen(str))

#define  sfree(sp)	free((void *)(sp))


char *
savestr(str)				/* save a string */
char  *str;
{
	register int  hval;
	register S_HEAD  *sp;

	if (str == NULL)
		return(NULL);
	hval = hash(str);
	for (sp = stab[hval]; sp != NULL; sp = sp->next)
		if (!strcmp(str, string(sp))) {
			sp->nl++;
			return(string(sp));
		}
	if ((sp = salloc(str)) == NULL) {
		eputs("Out of memory in savestr\n");
		quit(1);
	}
	strcpy(string(sp), str);
	sp->nl = 1;
	sp->next = stab[hval];
	stab[hval] = sp;
	return(string(sp));
}


void
freestr(s)				/* free a string */
char  *s;
{
	int  hval;
	register S_HEAD  *spl, *sp;

	if (s == NULL)
		return;
	hval = hash(s);
	for (spl = NULL, sp = stab[hval]; sp != NULL; spl = sp, sp = sp->next)
		if (s == string(sp)) {
			if (--sp->nl > 0)
				return;
			if (spl != NULL)
				spl->next = sp->next;
			else
				stab[hval] = sp->next;
			sfree(sp);
			return;
		}
}


int
shash(s)
register char  *s;
{
	register int  h = 0;

	while (*s)
		h = (h<<1 & 0x7fff) ^ (*s++ & 0xff);
	return(h);
}
