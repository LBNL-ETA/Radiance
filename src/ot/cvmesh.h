/* RCSid $Id: cvmesh.h,v 2.2 2003/03/14 21:27:46 greg Exp $ */
/*
 *  Header for Radiance triangle mesh conversion
 *
 *  Include after standard.h
 */

#include "octree.h"
#include "object.h"
#include "mesh.h"

extern MESH	*ourmesh;		/* global mesh pointer */

extern FVECT	meshbounds[2];		/* mesh bounding box */

#ifdef NOPROTO

extern MESH	*cvinit();
extern int	cvpoly();
extern int	cvtri();
extern void	cvmeshbounds();
extern MESH	*cvmesh();
void		wfreadobj();

#else

extern MESH	*cvinit(char *nm);
extern int	cvpoly(OBJECT mo, int n, FVECT *vp,
				FVECT *vn, FLOAT (*vc)[2]);
extern int	cvtri(OBJECT mo, FVECT vp1, FVECT vp2, FVECT vp3,
			FVECT vn1, FVECT vn2, FVECT vn3,
			FLOAT vc1[2], FLOAT vc2[2], FLOAT vc3[2]);
extern void	cvmeshbounds(void);
extern MESH	*cvmesh(void);
					/* defined in wfconv.c */
void		wfreadobj(char *objfn);


#endif /* NOPROTO */
