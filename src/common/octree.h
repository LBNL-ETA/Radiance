/* RCSid: $Id: octree.h,v 2.5 2003/02/22 02:07:22 greg Exp $ */
/*
 *  octree.h - header file for routines using octrees.
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

/*
 *	An octree is expressed as an integer which is either
 *	an index to eight other nodes, the empty tree, or an index
 *	to a set of objects.  If the octree has a value:
 *
 *		> -1:	it is an index to eight other nodes.
 *
 *		-1:	it is empty
 *
 *		< -1:	it is an index to a set of objects
 */

#ifndef  OCTREE
#define  OCTREE		int
#endif

#define  EMPTY		(-1)

#define  isempty(ot)	((ot) == EMPTY)
#define  isfull(ot)	((ot) < EMPTY)
#define  istree(ot)	((ot) > EMPTY)

#define  oseti(ot)	(-(ot)-2)	/* object set index */
#define  OCTBLKSIZ	04000		/* octree block size */
#define  octbi(ot)	((ot)>>11)	/* octree block index */
#define  octti(ot)	(((ot)&03777)<<3)/* octree index in block */

#ifndef  MAXOBLK
#ifdef  BIGMEM
#define  MAXOBLK	32767		/* maximum octree block */
#else
#define  MAXOBLK	4095		/* maximum octree block */
#endif
#endif

extern OCTREE  *octblock[MAXOBLK];	/* octree blocks */

#define  octkid(ot,br)	(octblock[octbi(ot)][octti(ot)+br])

/*
 *	The cube structure is used to hold an octree and its cubic
 *	boundaries.
 */

typedef struct {
	OCTREE  cutree;			/* the octree for this cube */
	FVECT  cuorg;			/* the cube origin */
	double  cusize;			/* the cube size */
}  CUBE;

extern CUBE  thescene;			/* the main scene */

				/* flags for reading and writing octrees */
#define  IO_CHECK	0		/* verify file type */
#define  IO_INFO	01		/* information header */
#define  IO_SCENE	02		/* objects */
#define  IO_TREE	04		/* octree */
#define  IO_FILES	010		/* object file names */
#define  IO_BOUNDS	020		/* octree boundary */
#define  IO_ALL		(~0)		/* everything */
				/* octree format identifier */
#define  OCTFMT		"Radiance_octree"
				/* magic number for octree files */
#define  MAXOBJSIZ	8		/* maximum sizeof(OBJECT) */
#define  OCTMAGIC	( 4 *MAXOBJSIZ+251)	/* increment first value */
				/* octree node types */
#define  OT_EMPTY	0
#define  OT_FULL	1
#define  OT_TREE	2
				/* return values for surface functions */
#define  O_MISS		0		/* no intersection */
#define  O_HIT		1		/* intersection */
#define  O_IN		2		/* cube contained entirely */

#ifdef NOPROTO

extern OCTREE	octalloc();
extern void	octfree();
extern void	octdone();
extern OCTREE	combine();
extern void	culocate();
extern void	cucopy();
extern int	incube();

extern int	readoct();

#else

extern OCTREE	octalloc(void);
extern void	octfree(OCTREE ot);
extern void	octdone(void);
extern OCTREE	combine(OCTREE ot);
extern void	culocate(CUBE *cu, FVECT pt);
extern void	cucopy(CUBE *cu1, CUBE *cu2);
extern int	incube(CUBE *cu, FVECT pt);

extern int	readoct(char *fname, int load, CUBE *scene, char *ofn[]);

#endif
