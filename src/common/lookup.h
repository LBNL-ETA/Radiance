/* RCSid: $Id: lookup.h,v 2.6 2003/02/22 02:07:22 greg Exp $ */
/*
 * Header file for general associative table lookup routines
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

typedef struct {
	char	*key;			/* key name */
	unsigned long	hval;		/* key hash value (for efficiency) */
	char	*data;			/* pointer to client data */
} LUENT;

typedef struct {
	unsigned long	(*hashf)();	/* key hash function */
	int	(*keycmp)();		/* key comparison function */
	void	(*freek)();		/* free a key */
	void	(*freed)();		/* free the data */
	int	tsiz;			/* current table size */
	LUENT	*tabl;			/* table, if allocated */
	int	ndel;			/* number of deleted entries */
} LUTAB;

#undef strcmp
#define LU_SINIT(fk,fd)		{lu_shash,strcmp,(void (*)())(fk),\
				(void (*)())(fd),0,NULL,0}

/*
 * The lu_init routine is called to initialize a table.  The number of
 * elements passed is not a limiting factor, as a table can grow to
 * any size permitted by memory.  However, access will be more efficient
 * if this number strikes a reasonable balance between default memory use
 * and the expected (minimum) table size.  The value returned is the
 * actual allocated table size (or zero if there was insufficient memory).
 *
 * The hashf, keycmp, freek and freed member functions must be assigned
 * separately.  If the hash value is sufficient to guarantee equality
 * between keys, then the keycmp pointer may be NULL.  Otherwise, it
 * should return 0 if the two passed keys match.  If it is not necessary
 * (or possible) to free the key and/or data values, then the freek and/or
 * freed member functions may be NULL.
 *
 * It isn't fully necessary to call lu_init to initialize the LUTAB structure.
 * If tsiz is 0, then the first call to lu_find will allocate a minimal table.
 * The LU_SINIT macro provides a convenient static declaration for character
 * string keys.
 *
 * The lu_find routine returns the entry corresponding to the given
 * key.  If the entry does not exist, the corresponding key field will
 * be NULL.  If the entry has been previously deleted but not yet freed,
 * then only the data field will be NULL.  It is the caller's
 * responsibility to (allocate and) assign the key and data fields when
 * creating a new entry.  The only case where lu_find returns NULL is when
 * the system has run out of memory.
 *
 * The lu_delete routine frees an entry's data (if any) by calling
 * the freed member function, but does not free the key field.  This
 * will be freed later during (or instead of) table reallocation.
 * It is therefore an error to reuse or do anything with the key
 * field after calling lu_delete.
 *
 * The lu_doall routine loops through every filled table entry, calling
 * the given function once on each entry.  If a NULL pointer is passed
 * for this function, then lu_doall simply returns the total number of
 * active entries.  Otherwise, it returns the sum of all the function
 * evaluations.
 *
 * The lu_done routine calls the given free function once for each
 * assigned table entry (i.e. each entry with an assigned key value).
 * The user must define these routines to free the key and the data
 * in the LU_TAB structure.  The final action of lu_done is to free the
 * allocated table itself.
 */

extern int	strcmp();

#ifdef NOPROTO

extern int	lu_init();
extern LUENT	*lu_find();
extern void	lu_delete();
extern int	lu_doall();
extern void	lu_done();
extern unsigned long	lu_shash();

#else

extern int	lu_init(LUTAB *tbl, int nel);
extern unsigned long	lu_shash(char *s);
extern LUENT	*lu_find(LUTAB *tbl, char *key);
extern void	lu_delete(LUTAB *tbl, char *key);
extern int	lu_doall(LUTAB *tbl, int (*f)());
extern void	lu_done(LUTAB *tbl);

#endif
