/* RCSid: $Id: sm_geom.h,v 3.12 2003/02/22 02:07:25 greg Exp $ */
/*
 *  sm_geom.h
 */

/* Assumes included after standard.h  */

#include <values.h>

#ifdef SMLFLT
#define EQUALITY_EPS 1e-6
#else
#define EQUALITY_EPS 1e-10
#endif

#define F_TINY 1e-10
#define FZERO(x) ((x) < F_TINY && (x) > -F_TINY)
#define FEQUAL(a,b) FZERO((a) - (b))

#define ZERO(x) ((x) < FTINY && (x) > -FTINY)
#define EQUAL(a,b) ZERO((a) - (b))

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

typedef struct _FPEQ {
    FVECT n;
    double d;
    char x,y,z;
}FPEQ;

#define FP_N(f) ((f).n)
#define FP_D(f) ((f).d)
#define FP_X(f) ((f).x) 
#define FP_Y(f) ((f).y) 
#define FP_Z(f) ((f).z) 
  
typedef long BCOORD;
typedef long BDIR;

#define BITS_BCOORD     (BITS(long)-2)
#define SHIFT_MAXBCOORD (BITS_BCOORD-1)   
#define MAXBCOORD      ((1L << BITS_BCOORD)-1)
#define MAXBCOORD2      (MAXBCOORD>>1)
#define MAXBCOORD4      (MAXBCOORD2>>1)

#define M_2_3_PI PI*2/3

#ifndef INVALID
#define INVALID -1
#endif
#define IADDV3(v,a)  ((v)[0] += (a)[0],(v)[1] += (a)[1],(v)[2] += (a)[2])
#define ISUBV3(v,a)  ((v)[0] -= (a)[0],(v)[1] -= (a)[1],(v)[2] -= (a)[2])
#define ISCALEV3(v,a)  ((v)[0] *= (a),(v)[1] *= (a),(v)[2] *= (a))
#define IDIVV3(v,a)  ((v)[0] /= (a),(v)[1] /= (a),(v)[2] /= (a))


#define ADDV3(v,a,b) ((v)[0] = (a)[0]+(b)[0],(v)[1] = (a)[1]+(b)[1],\
		      (v)[2] = (a)[2]+(b)[2])
#define SUBV3(v,a,b) ((v)[0] = (a)[0]-(b)[0],(v)[1] = (a)[1]-(b)[1],\
		      (v)[2] = (a)[2]-(b)[2])
#define SUMV3(v,a,b,s) ((v)[0] = (a)[0]+(s)*(b)[0],(v)[1]=(a)[1]+(s)*(b)[1],\
		      (v)[2] = (a)[2]+(s)*(b)[2])
#define SCALEV3(v,a,s)  ((v)[0]=(a)[0]*(s),(v)[1]=(a)[1]*(s),(v)[2]=(a)[2]*(s))
#define ZERO_VEC3(v)     (ZERO(v[0]) && ZERO(v[1]) && ZERO(v[2]) )
#define EQUAL_VEC3(a,b)  (FEQUAL(a[0],b[0])&&FEQUAL(a[1],b[1])&&FEQUAL(a[2],b[2]))
#define OPP_EQUAL_VEC3(a,b)  (EQUAL(a[0],-b[0])&&EQUAL(a[1],-b[1])&&EQUAL(a[2],-b[2]))
#define FZERO_VEC3(v)     (FZERO(v[0]) && FZERO(v[1]) && FZERO(v[2]) )
#define FEQUAL_VEC3(a,b) (FEQUAL(a[0],b[0])&&FEQUAL(a[1],b[1])&&FEQUAL(a[2],b[2]))
#define NEGATE_VEC3(v)   ((v)[0] *= -1.0,(v)[1] *= -1.0,(v)[2] *= -1.0)
#define COPY_VEC2(v1,v2) ((v1)[0]=(v2)[0],(v1)[1]=(v2)[1])
#define DIST(a,b)         (sqrt(((a)[0]-(b)[0])*((a)[0]-(b)[0]) + \
			        ((a)[1]-(b)[1])*((a)[1]-(b)[1]) + \
				((a)[2]-(b)[2])*((a)[2]-(b)[2])))
#define DIST_SQ(a,b)      (((a)[0]-(b)[0])*((a)[0]-(b)[0]) + \
			        ((a)[1]-(b)[1])*((a)[1]-(b)[1]) + \
				((a)[2]-(b)[2])*((a)[2]-(b)[2]))

#define SIGN(x) ((x<0)?-1:1)
#define CROSS_VEC2(v1,v2) (((v1)[0]*(v2)[1]) - ((v1)[1]*(v2)[0])) 
#define DOT_VEC2(v1,v2) ((v1)[0]*(v2)[0] + (v1)[1]*(v2)[1])

#define EDGE_MIDPOINT(a,v1,v2) ((a)[0]=((v1)[0]+(v2)[0])*0.5, \
(a)[1]=((v1)[1]+(v2)[1])*0.5,(a)[2] = ((v1)[2]+(v2)[2])*0.5)

#define MIN_VEC3(v) ((v)[0]<(v)[1]?((v)[0]<(v)[2]?(v)[0]:v[2]): \
		     (v)[1]<(v)[2]?(v)[1]:(v)[2])
#define MAX3(a,b,c) (((b)>(a))?((b) > (c))?(b):(c):((a)>(c))?(a):(c))   
#define MIN3(a,b,c) (((b)<(a))?((b) < (c))?(b):(c):((a)<(c))?(a):(c))

#define MAX(a,b)    (((b)>(a))?(b):(a))   
#define MIN(a,b)  (((b)<(a))?(b):(a))   

#define SUM_3VEC3(r,a,b,c) ((r)[0]=(a)[0]+(b)[0]+(c)[0], \
   (r)[1]=(a)[1]+(b)[1]+(c)[1],(r)[2]=(a)[2]+(b)[2]+(c)[2])

#define NTH_BIT(n,i)         ((n) & (1<<(i)))
#define SET_NTH_BIT(n,i)     ((n) |= (1<<(i)))   

#define PT_ON_PLANE(p,peq) (DOT(FP_N(peq),p)+FP_D(peq))

/* FUNCTIONS:
   int point_in_cone(FVECT p,a,b,c) 
   void triangle_centroid(FVECT v0,v1,v2,c)
   double tri_normal(FVECT v0,v1,v2,n,int norm)
   void tri_plane_equation(FVECT v0,v1,v2,FPEQ *peqptr,int norm)
   int intersect_ray_plane(FVECT orig,dir,FPEQ peq,double *pd,FVECT r)
   double point_on_sphere(FVECT ps,p,c)
   int point_in_stri(FVECT v0,v1,v2,p)
   int ray_intersect_tri(FVECT orig,dir,v0,v1,v2,pt)
   void calculate_view_frustum(FVECT vp,hv,vv,double horiz,vert,near,far,
                               FVECT fnear[4],ffar[4])
   void bary2d(double x1,y1,x2,y2,x3,y3,px,py,coord)
*/
double tri_normal();
double point_on_sphere();



