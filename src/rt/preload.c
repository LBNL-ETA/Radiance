#ifndef lint
static const char	RCSid[] = "$Id: preload.c,v 2.9 2005/09/15 18:06:23 greg Exp $";
#endif
/*
 * Preload associated object structures to maximize memory sharing.
 *
 *  External symbols declared in ray.h
 */

#include "copyright.h"

#include "standard.h"
#include "octree.h"
#include "object.h"
#include "otypes.h"
#include "face.h"
#include "cone.h"
#include "instance.h"
#include "mesh.h"
#include "color.h"
#include "data.h"
#include "func.h"
#include "ray.h"


/* KEEP THIS ROUTINE CONSISTENT WITH THE DIFFERENT OBJECT FUNCTIONS! */


extern int
load_os(			/* load associated data for object */
	register OBJREC	*op
)
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
	case OBJ_MESH:		/* mesh instance */
		getmeshinst(op, IO_ALL);
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
	return 0; /* pro forma return */
}


extern void
preload_objs(void)		/* preload object data structures */
{
	register OBJECT on;
				/* note that nobjects may change during loop */
	for (on = 0; on < nobjects; on++)
		load_os(objptr(on));
}
