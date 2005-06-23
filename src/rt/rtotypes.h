/* RCSid $Id: rtotypes.h,v 1.2 2005/06/23 11:51:47 greg Exp $ */
/*
 * External functions implementing Radiance object types
 */

typedef int otype_tracef(OBJREC *o, RAY *r);

extern otype_tracef o_sphere, o_face, o_cone, o_instance, o_mesh;

extern otype_tracef m_aniso, m_dielectric, m_glass, m_alias, m_light,
	m_normal, m_mist, m_mirror, m_direct, m_clip, m_brdf, m_brdf2;

extern otype_tracef t_func, t_data;

extern otype_tracef p_cfunc, p_bfunc, p_pdata, p_cdata, p_bdata;

extern otype_tracef mx_func, mx_data, mx_pdata;

extern otype_tracef do_text;

	/* text.c */
extern void freetext(OBJREC *m);
