#ifndef lint
static const char	RCSid[] = "$Id: preload.c,v 2.6 2003/02/22 02:07:29 greg Exp $";
#endif
/*
 * Preload associated object structures to maximize memory sharing.
 *
 *  External symbols declared in ray.h
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
#include "octree.h"
#include "object.h"
#include "otypes.h"
#include "face.h"
#include "cone.h"
#include "instance.h"
#include "color.h"
#include "data.h"


/* KEEP THIS ROUTINE CONSISTENT WITH THE DIFFERENT OBJECT FUNCTIONS! */


int
load_os(op)			/* load associated data for object */
register OBJREC	*op;
{
	DATARRAY  *dp;

	switch (op->otype) {
	case OBJ_FACE:		/* polygon */
		getface(op);
		return(1);
	case OBJ_CONE:		/* cone */
	case OBJ_RING:		/* disk */
	case OBJ_CYLINDER:	/* cylinder */
	case OBJ_CUP:		/* inverted cone */
	case OBJ_TUBE:		/* inverted cylinder */
		getcone(op, 1);
		return(1);
	case OBJ_INSTANCE:	/* octree instance */
		getinstance(op, IO_ALL);
		return(1);
	case PAT_CPICT:		/* color picture */
		if (op->oargs.nsargs < 4)
			goto sargerr;
		getpict(op->oargs.sarg[3]);
		getfunc(op, 4, 0x3<<5, 0);
		return(1);
	case PAT_CDATA:		/* color data */
		dp = getdata(op->oargs.sarg[3]);
		getdata(op->oargs.sarg[4]);
		getdata(op->oargs.sarg[5]);
		getfunc(op, 6, ((1<<dp->nd)-1)<<7, 0);
		return(1);
	case PAT_BDATA:		/* brightness data */
		if (op->oargs.nsargs < 2)
			goto sargerr;
		dp = getdata(op->oargs.sarg[1]);
		getfunc(op, 2, ((1<<dp->nd)-1)<<3, 0);
		return(1);
	case PAT_BFUNC:		/* brightness function */
		getfunc(op, 1, 0x1, 0);
		return(1);
	case PAT_CFUNC:		/* color function */
		getfunc(op, 3, 0x7, 0);
		return(1);
	case TEX_DATA:		/* texture data */
		if (op->oargs.nsargs < 6)
			goto sargerr;
		dp = getdata(op->oargs.sarg[3]);
		getdata(op->oargs.sarg[4]);
		getdata(op->oargs.sarg[5]);
		getfunc(op, 6, ((1<<dp->nd)-1)<<7, 1);
		return(1);
	case TEX_FUNC:		/* texture function */
		getfunc(op, 3, 0x7, 1);
		return(1);
	case MIX_DATA:		/* mixture data */
		dp = getdata(op->oargs.sarg[3]);
		getfunc(op, 4, ((1<<dp->nd)-1)<<5, 0);
		return(1);
	case MIX_PICT:		/* mixture picture */
		getpict(op->oargs.sarg[3]);
		getfunc(op, 4, 0x3<<5, 0);
		return(1);
	case MIX_FUNC:		/* mixture function */
		getfunc(op, 3, 0x4, 0);
		return(1);
	case MAT_PLASTIC2:	/* anisotropic plastic */
	case MAT_METAL2:	/* anisotropic metal */
		getfunc(op, 3, 0x7, 1);
		return(1);
	case MAT_BRTDF:		/* BRDTfunc material */
		getfunc(op, 9, 0x3f, 0);
		return(1);
	case MAT_PDATA:		/* plastic BRDF data */
	case MAT_MDATA:		/* metal BRDF data */
	case MAT_TDATA:		/* trans BRDF data */
		if (op->oargs.nsargs < 2)
			goto sargerr;
		getdata(op->oargs.sarg[1]);
		getfunc(op, 2, 0, 0);
		return(1);
	case MAT_PFUNC:		/* plastic BRDF func */
	case MAT_MFUNC:		/* metal BRDF func */
	case MAT_TFUNC:		/* trans BRDF func */
		getfunc(op, 1, 0, 0);
		return(1);
	case MAT_DIRECT1:	/* prism1 material */
		getfunc(op, 4, 0xf, 1);
		return(1);
	case MAT_DIRECT2:	/* prism2 material */
		getfunc(op, 8, 0xff, 1);
		return(1);
	}
			/* nothing to load for the remaining types */
	return(0);
sargerr:
	objerror(op, USER, "too few string arguments");
}


void
preload_objs()		/* preload object data structures */
{
	register OBJECT on;
				/* note that nobjects may change during loop */
	for (on = 0; on < nobjects; on++)
		load_os(objptr(on));
}
