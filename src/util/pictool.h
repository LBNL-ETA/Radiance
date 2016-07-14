/* RCSid $Id: pictool.h,v 2.3 2016/07/14 17:32:12 greg Exp $ */
#ifndef	__PICTOOL_H
#define __PICTOOL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>

#include "g3vector.h"
#include "g3affine.h"
#include "g3flist.h"

#include "view.h"
#include "color.h"

#define PICT_GLARE

#define	PICT_VIEW	1
#define	PICT_LUT	2

#define	PICT_CH1	1
#define	PICT_CH2	2
#define	PICT_CH3	4

#define	PICT_ALLCH	7

#ifdef PICT_GLARE
typedef	struct	pixinfo {
	int		gsn;
	int		pgs;
	double	omega;
	FVECT	dir;
} pixinfo;
#endif

typedef struct	pict {
	int			use_lut;
	char		comment[4096];
	int			valid_view;
	VIEW		view;
	COLOR*		lut;
	RESOLU		resol;
	COLOR*		pic;
	double		pjitter;
#ifdef PICT_GLARE
	pixinfo*	pinfo;
	g3FList*	glareinfo;
#endif
} pict;

#define	pict_get_origin(p)		(p->view.vp)
#define	pict_get_viewdir(p)		(p->view.vdir)
#define	pict_get_viewup(p)		(p->view.vup)

#define	pict_lut_used(p)		(p->use_lut)

#define	pict_get_comp(p,x,y,i)	((p->pic + x*p->resol.yr)[y][i])
#define	pict_colbuf(p,buf,x,y)	((buf + y*p->resol.xr)[x])
#define	pict_get_color(p,x,y)	pict_colbuf(p,p->pic,x,y)
#define	pict_get_xsize(p)		(p->resol.xr)
#define	pict_get_ysize(p)		(p->resol.yr)



#ifdef PICT_GLARE

#define	pict_get_omega(p,x,y)	((p->pinfo[y + x*p->resol.yr])).omega
#define	pict_get_cached_dir(p,x,y)	(((p->pinfo[y + x*p->resol.yr])).dir)
#define	pict_get_gsn(p,x,y)		((p->pinfo[y + x*p->resol.yr])).gsn
#define	pict_get_pgs(p,x,y)		((p->pinfo[y + x*p->resol.yr])).pgs

enum {
	PICT_GSN = 0,
	PICT_Z1_GSN ,
	PICT_Z2_GSN ,
	PICT_NPIX,
	PICT_AVPOSX,
	PICT_AVPOSY,
	PICT_AVLUM,
	PICT_AVOMEGA,
	PICT_DXMAX,
	PICT_DYMAX,
	PICT_LMIN,
	PICT_LMAX,
	PICT_EGLARE,
	PICT_DGLARE,
	PICT_GLSIZE
};

#define	pict_get_num(p,i)		(g3fl_get(p->glareinfo,i)[PICT_GSN])
#define	pict_get_z1_gsn(p,i)	(g3fl_get(p->glareinfo,i)[PICT_Z1_GSN])
#define	pict_get_z2_gsn(p,i)	(g3fl_get(p->glareinfo,i)[PICT_Z2_GSN])
#define	pict_get_npix(p,i)		(g3fl_get(p->glareinfo,i)[PICT_NPIX])
#define	pict_get_av_posx(p,i)	(g3fl_get(p->glareinfo,i)[PICT_AVPOSX])
#define	pict_get_av_posy(p,i)	(g3fl_get(p->glareinfo,i)[PICT_AVPOSY])
#define	pict_get_av_lum(p,i)	(g3fl_get(p->glareinfo,i)[PICT_AVLUM])
#define	pict_get_av_omega(p,i)	(g3fl_get(p->glareinfo,i)[PICT_AVOMEGA])
#define	pict_get_dx_max(p,i)	(g3fl_get(p->glareinfo,i)[PICT_DXMAX])
#define	pict_get_dy_max(p,i)	(g3fl_get(p->glareinfo,i)[PICT_DYMAX])
#define	pict_get_lum_min(p,i)	(g3fl_get(p->glareinfo,i)[PICT_LMIN])
#define	pict_get_lum_max(p,i)	(g3fl_get(p->glareinfo,i)[PICT_LMAX])
#define	pict_get_Eglare(p,i)	(g3fl_get(p->glareinfo,i)[PICT_EGLARE])
#define	pict_get_Dglare(p,i)	(g3fl_get(p->glareinfo,i)[PICT_DGLARE])

#define pict_get_gli_size(p)	(p->glareinfo->size)

#ifdef PICT_GLARE
#define pict_is_validpixel(p,x,y)	(((p->view.type != VT_ANG) ||  \
									(DOT(pict_get_cached_dir(p,x,y), p->view.vdir) >= 0.0))\
									&& pict_get_omega(p,x,y) > 0.0) 
#else
#define pict_is_validpixel(p,x,y)	(((p->view.type != VT_ANG) ||  \
									(DOT(pict_get_cached_dir(p,x,y), p->view.vdir) >= 0.0))\
									&& pict_get_sangle(p,x,y) > 0)

#endif

void	pict_new_gli(pict* p);			

void pict_update_evalglare_caches(pict* p);

#endif

int		pict_init(pict* p);
pict*	pict_create();
void	pict_free(pict* p);

pict*	pict_get_copy(pict* p);

void	pict_set_comment(pict* p,const char* c);

int		pict_resize(pict* p,int x,int y);

int		pict_zero(pict* p);

int		pict_read(pict* p,char* fn);
int     pict_read_fp(pict* p,FILE* fp);
int		pict_read_tt_txt(pict* p,char* fn);
int     pict_read_from_list(pict* p,FILE* fp);
int		pict_write(pict* p,char* fn);
int     pict_write_fp(pict* p,FILE* fp);
int     pict_write_dir(pict* p,FILE* fp);

int		pict_setup_lut(pict* p,char* theta_fn,char* phi_fn,char* omega_fn);
int		pict_update_lut(pict* p);

int		pict_update_view(pict* p);

int		pict_get_dir(pict* p,int x,int y,FVECT dir);
int		pict_get_dir_ic(pict* p,double x,double y,FVECT dir);
int		pict_get_ray(pict* p,int x,int y,FVECT orig,FVECT dir);
int		pict_get_ray_ic(pict* p,double x,double y,FVECT orig,FVECT dir);
int		pict_locate(pict* p,FVECT pt,int* x,int* y);
double	pict_get_sangle(pict* p,int x,int y);

int pict_get_vangle(pict* p,double x,double y,FVECT vdir,FVECT vup,double* a);
int pict_get_hangle(pict* p,double x,double y,FVECT vdir,FVECT vup,double* a);
int	pict_get_tau(pict* p,double x,double y,FVECT vdir,FVECT vup,double* t);
int	pict_get_sigma(pict* p,double x,double y,FVECT vdir,FVECT vup,double* s);

int		pict_set_pj(pict* p,double pj);
int		pict_setorigin(pict* p,FVECT orig);
int		pict_setdir(pict* p,FVECT dir);
int		pict_setvtype(pict* p,char vt);
int		pict_setvup(pict* p,FVECT vup);
int		pict_sethview(pict* p,double ha);
int		pict_setvview(pict* p,double va);

int		pict_trans_cam(pict* p,g3ATrans t);




#ifdef __cplusplus
}
#endif
#endif
