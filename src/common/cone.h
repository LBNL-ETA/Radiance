/* RCSid: $Id: cone.h,v 2.2 2003/02/22 02:07:22 greg Exp $ */
/*
 *  cone.h - header file for cones (cones, cylinders, rings, cups, tubes).
 *
 *	Storage of arguments in the cone structure is a little strange.
 *  To save space, we use an index into the real arguments of the
 *  object structure through ca.  The indices are for the axis
 *  endpoints and radii:  p0, p1, r0 and r1.
 *
 *     2/12/86
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

typedef struct cone {
	FLOAT  *ca;		/* cone arguments (o->oargs.farg) */
	char  p0, p1;		/* indices for endpoints */
	char  r0, r1;		/* indices for radii */
	FVECT  ad;		/* axis direction vector */
	FLOAT  al;		/* axis length */
	FLOAT  sl;		/* side length */
	FLOAT  (*tm)[4];	/* pointer to transformation matrix */
}  CONE;

#define  CO_R0(co)	((co)->ca[(co)->r0])
#define  CO_R1(co)	((co)->ca[(co)->r1])
#define  CO_P0(co)	((co)->ca+(co)->p0)
#define  CO_P1(co)	((co)->ca+(co)->p1)

#ifdef NOPROTO

extern CONE  *getcone();
extern void  freecone();
extern void  conexform();

#else

extern CONE  *getcone(OBJREC *o, int getxf);
extern void  freecone(OBJREC *o);
extern void  conexform(CONE *co);

#endif
