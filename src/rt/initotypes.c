#ifndef lint
static const char	RCSid[] = "$Id: initotypes.c,v 2.7 2003/02/22 02:07:28 greg Exp $";
#endif
/*
 * Initialize ofun[] list for renderers
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

#include  "standard.h"

#include  "otypes.h"

#include  "otspecial.h"

extern int  o_sphere();
extern int  o_face();
extern int  o_cone();
extern int  o_instance();
extern int  m_light();
extern int  m_normal();
extern int  m_aniso();
extern int  m_dielectric();
extern int  m_mist();
extern int  m_glass();
extern int  m_clip();
extern int  m_mirror();
extern int  m_direct();
extern int  m_brdf();
extern int  m_brdf2();
extern int  t_func(), t_data();
extern int  p_cfunc(), p_bfunc();
extern int  p_pdata(), p_cdata(), p_bdata();
extern int  mx_func(), mx_data(), mx_pdata();
extern int  do_text();

FUN  ofun[NUMOTYPE] = INIT_OTYPE;


initotypes()			/* initialize ofun array */
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
	ofun[MAT_LIGHT].funp =
	ofun[MAT_ILLUM].funp =
	ofun[MAT_GLOW].funp =
	ofun[MAT_SPOT].funp = m_light;
	ofun[MAT_PLASTIC].funp =
	ofun[MAT_METAL].funp =
	ofun[MAT_TRANS].funp = m_normal;
	ofun[MAT_TRANS].flags |= T_IRR_IGN;
	ofun[MAT_PLASTIC2].funp =
	ofun[MAT_METAL2].funp =
	ofun[MAT_TRANS2].funp = m_aniso;
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


o_default()			/* default action is error */
{
	error(INTERNAL, "unexpected object call");
				/* call to pull in freeobjmem.o */
	free_objs(0, 0);
}
