/* RCSid $Id: rtotypes.h,v 1.5 2019/04/18 22:35:01 greg Exp $ */
/*
 * External functions implementing Radiance object types
 */

#ifndef _RAD_RTOTYPES_H_
#define _RAD_RTOTYPES_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef int otype_implf(OBJREC *o, RAY *r);

extern otype_implf o_sphere, o_face, o_cone, o_instance, o_mesh;

extern otype_implf m_aniso, m_dielectric, m_glass, m_alias, m_light,
	m_normal, m_mist, m_mirror, m_direct, m_clip, m_brdf, m_brdf2,
	m_bsdf, m_ashikhmin;

extern otype_implf t_func, t_data;

extern otype_implf p_cfunc, p_bfunc, p_pdata, p_cdata, p_bdata;

extern otype_implf mx_func, mx_data, mx_pdata;

extern otype_implf do_text;

	/* text.c */
extern void freetext(OBJREC *m);

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RTOTYPES_H_ */
