#ifndef lint
static const char RCSid[] = "$Id: initotypes.c,v 2.14 2004/10/25 17:07:35 greg Exp $";
#endif
/*
 * Initialize ofun[] list for renderers
 */

#include  "copyright.h"

#include  "standard.h"
#include  "otypes.h"
#include  "rtotypes.h"
#include  "otspecial.h"


FUN  ofun[NUMOTYPE] = INIT_OTYPE;


extern void
initotypes(void)			/* initialize ofun array */
{
	ofun[OBJ_SPHERE].funp =
	ofun[OBJ_BUBBLE].funp = o_sphere;
	ofun[OBJ_FACE].funp = o_face;
	ofun[OBJ_CONE].funp =
	ofun[OBJ_CUP].funp =
	ofun[OBJ_CYLINDER].funp =
	ofun[OBJ_TUBE].funp =
	ofun[OBJ_RING].funp = o_cone;
	ofun[OBJ_INSTANCE].funp = o_instance;
	ofun[OBJ_MESH].funp = o_mesh;
	ofun[MOD_ALIAS].funp = m_alias;
	ofun[MAT_LIGHT].funp =
	ofun[MAT_ILLUM].funp =
	ofun[MAT_GLOW].funp =
	ofun[MAT_SPOT].funp = m_light;
	ofun[MAT_LIGHT].flags |= T_OPAQUE;
	ofun[MAT_SPOT].flags |= T_OPAQUE;
	ofun[MAT_PLASTIC].funp =
	ofun[MAT_METAL].funp =
	ofun[MAT_TRANS].funp = m_normal;
	ofun[MAT_PLASTIC].flags |= T_OPAQUE;
	ofun[MAT_METAL].flags |= T_OPAQUE;
	ofun[MAT_TRANS].flags |= T_IRR_IGN;
	ofun[MAT_PLASTIC2].funp =
	ofun[MAT_METAL2].funp =
	ofun[MAT_TRANS2].funp = m_aniso;
	ofun[MAT_PLASTIC2].flags |= T_OPAQUE;
	ofun[MAT_METAL2].flags |= T_OPAQUE;
	ofun[MAT_TRANS2].flags |= T_IRR_IGN;
	ofun[MAT_DIELECTRIC].funp =
	ofun[MAT_INTERFACE].funp = m_dielectric;
	ofun[MAT_DIELECTRIC].flags |= T_IRR_IGN;
	ofun[MAT_INTERFACE].flags |= T_IRR_IGN;
	ofun[MAT_MIST].funp = m_mist;
	ofun[MAT_MIST].flags |= T_IRR_IGN;
	ofun[MAT_GLASS].funp = m_glass;
	ofun[MAT_GLASS].flags |= T_IRR_IGN;
	ofun[MAT_MIRROR].funp = m_mirror;
	ofun[MAT_DIRECT1].funp =
	ofun[MAT_DIRECT2].funp = m_direct;
	ofun[MAT_CLIP].funp = m_clip;
	ofun[MAT_BRTDF].funp = m_brdf;
	ofun[MAT_PFUNC].funp =
	ofun[MAT_MFUNC].funp =
	ofun[MAT_PDATA].funp =
	ofun[MAT_MDATA].funp = 
	ofun[MAT_TFUNC].funp =
	ofun[MAT_TDATA].funp = m_brdf2;
	ofun[MAT_PFUNC].flags |= T_OPAQUE;
	ofun[MAT_MFUNC].flags |= T_OPAQUE;
	ofun[MAT_PDATA].flags |= T_OPAQUE;
	ofun[MAT_MDATA].flags |= T_OPAQUE;
	ofun[TEX_FUNC].funp = t_func;
	ofun[TEX_DATA].funp = t_data;
	ofun[PAT_CFUNC].funp = p_cfunc;
	ofun[PAT_BFUNC].funp = p_bfunc;
	ofun[PAT_CPICT].funp = p_pdata;
	ofun[PAT_CDATA].funp = p_cdata;
	ofun[PAT_BDATA].funp = p_bdata;
	ofun[PAT_CTEXT].funp =
	ofun[PAT_BTEXT].funp =
	ofun[MIX_TEXT].funp = do_text;
	ofun[MIX_FUNC].funp = mx_func;
	ofun[MIX_DATA].funp = mx_data;
	ofun[MIX_PICT].funp = mx_pdata;
}


extern int
o_default(OBJREC *o, RAY *r)			/* default action is error */
{
	objerror(o, CONSISTENCY, "unexpected object call");
				/* call to pull in freeobjmem.o */ /* XXX ? */
	free_objs(0, 0);
}
