#ifndef lint
static const char	RCSid[] = "$Id: sm_geom.c,v 3.13 2003/02/22 02:07:25 greg Exp $";
#endif
/*
 *  sm_geom.c
 *    some geometric utility routines
 */

#include "standard.h"
#include "sm_geom.h"

/*
 * int
 * pt_in_cone(p,a,b,c)
 *             : test if point p lies in cone defined by a,b,c and origin
 * double p[3];              : point to test
 * double a[3],b[3],c[3];    : points forming triangle
 * 
 *  Assumes apex at origin, a,b,c are unit vectors defining the
 *  triangle which the cone circumscribes. Assumes p is also normalized
 *  Test is implemented as:
 *             r =  (b-a)X(c-a) 
 *             in = (p.r) > (a.r) 
 *  The center of the cone is r, and corresponds to the triangle normal.
 *  p.r is the proportional to the cosine of the angle between p and the
 *  the cone center ray, and a.r to the radius of the cone. If the cosine
 *  of the angle for p is greater than that for a, the angle between p
 *  and r is smaller, and p lies in the cone.
 */
int
pt_in_cone(p,a,b,c)
double p[3],a[3],b[3],c[3];
{
  double r[3];
  double pr,ar;
  double ab[3],ac[3];
 
#ifdef DEBUG
#if DEBUG > 1
{
  double l;
  VSUB(ab,b,a);
  normalize(ab);
  VSUB(ac,c,a);
  normalize(ac);
  VCROSS(r,ab,ac);
  l = normalize(r);
  /* l = sin@ between ab,ac - if 0 vectors are colinear */
  if( l <= COLINEAR_EPS)
  {
    eputs("pt in cone: null triangle:returning FALSE\n");
    return(FALSE);
  }
}
#endif
#endif

  VSUB(ab,b,a);
  VSUB(ac,c,a);
  VCROSS(r,ab,ac);

  pr = DOT(p,r);	
  ar = DOT(a,r);
  /* Need to check for equality for degeneracy of 4 points on circle */
  if( pr > ar *( 1.0 + EQUALITY_EPS))
    return(TRUE); 
  else
    return(FALSE);
}

/*
 * tri_centroid(v0,v1,v2,c)
 *             : Average triangle vertices to give centroid: return in c
 *FVECT v0,v1,v2,c; : triangle vertices(v0,v1,v2) and vector to hold result(c)
 */                   
tri_centroid(v0,v1,v2,c)
FVECT v0,v1,v2,c;
{
  c[0] = (v0[0] + v1[0] + v2[0])/3.0;
  c[1] = (v0[1] + v1[1] + v2[1])/3.0;
  c[2] = (v0[2] + v1[2] + v2[2])/3.0;
}


/*
 * double  
 * tri_normal(v0,v1,v2,n,norm) : Calculates the normal of a face contour using
 *                            Newell's formula.
 * FVECT v0,v1,v2,n;  : Triangle vertices(v0,v1,v2) and vector for result(n)
 * int norm;          : If true result is normalized
 *
 * Triangle normal is calculated using the following:
 *    A =  SUMi (yi - yi+1)(zi + zi+1);
 *    B =  SUMi (zi - zi+1)(xi + xi+1)
 *    C =  SUMi (xi - xi+1)(yi + yi+1)
*/
double
tri_normal(v0,v1,v2,n,norm)
FVECT v0,v1,v2,n;
int norm;
{
  double mag;

  n[0] = (v0[2] + v1[2]) * (v0[1] - v1[1]) + 
         (v1[2] + v2[2]) * (v1[1] - v2[1])  +
         (v2[2] + v0[2]) * (v2[1] - v0[1]);
  n[1] = (v0[2] - v1[2]) * (v0[0] + v1[0]) +
	   (v1[2] - v2[2]) * (v1[0] + v2[0]) +
	    (v2[2] - v0[2]) * (v2[0] + v0[0]);
  n[2] = (v0[1] + v1[1]) * (v0[0] - v1[0]) +
         (v1[1] + v2[1]) * (v1[0] - v2[0]) +
         (v2[1] + v0[1]) * (v2[0] - v0[0]);
  if(!norm)
     return(0);
  mag = normalize(n);
  return(mag);
}

/*
 * tri_plane_equation(v0,v1,v2,peqptr,norm)
 *             : Calculates the plane equation (A,B,C,D) for triangle
 *               v0,v1,v2 ( Ax + By + Cz = D)
 * FVECT v0,v1,v2;   : Triangle vertices 
 * FPEQ *peqptr;     : ptr to structure to hold result
 *  int norm;        : if TRUE, return unit normal
 */
tri_plane_equation(v0,v1,v2,peqptr,norm)
   FVECT v0,v1,v2; 
   FPEQ *peqptr;
   int norm;
{
    tri_normal(v0,v1,v2,FP_N(*peqptr),norm);
    FP_D(*peqptr) = -(DOT(FP_N(*peqptr),v0));
}

/*
 * int 
 * intersect_ray_plane(orig,dir,peq,pd,r)
 *             : Intersects ray (orig,dir) with plane (peq). Returns TRUE
 *             if intersection occurs. If r!=NULL, sets with resulting i
 *             intersection point, and pd is set with parametric value of the
 *             intersection.
 * FVECT orig,dir;      : vectors defining the ray
 * FPEQ peq;            : plane equation
 * double *pd;          : holds resulting parametric intersection point
 * FVECT r;             : holds resulting intersection point
 *
 *    Plane is Ax + By + Cz +D = 0:
 *    A(orig[0] + dxt) + B(orig[1] + dyt) + C(orig[2] + dzt) + pd = 0
 *     t = -(DOT(plane_n,orig)+ plane_d)/(DOT(plane_n,d))
 *      line is  l = p1 + (p2-p1)t 
 *    Solve for t
 */
int
intersect_ray_plane(orig,dir,peq,pd,r)
   FVECT orig,dir;
   FPEQ peq;
   double *pd;
   FVECT r;
{
  double t,d;
  int hit;

  d = DOT(FP_N(peq),dir);
  if(ZERO(d))
     return(0);
  t =  -(DOT(FP_N(peq),orig) + FP_D(peq))/d; 

  if(t < 0)
       hit = 0;
    else
       hit = 1;
  if(r)
     VSUM(r,orig,dir,t);

    if(pd)
       *pd = t;
  return(hit);
}

/*
 * double
 * point_on_sphere(ps,p,c) : normalize p relative to sphere with center c
 * FVECT ps,p,c;       : ps Holds result vector,p is the original point,
 *                       and c is the sphere center
 */
double
point_on_sphere(ps,p,c)
FVECT ps,p,c;
{
  double d;
  
  VSUB(ps,p,c);    
  d = normalize(ps);
  return(d);
}

/*
 * int
 * point_in_stri(v0,v1,v2,p) : Return TRUE if p is in pyramid defined by
 *                            tri v0,v1,v2 and origin
 * FVECT v0,v1,v2,p;     :Triangle vertices(v0,v1,v2) and point in question(p)
 *
 *   Tests orientation of p relative to each edge (v0v1,v1v2,v2v0), if it is
 *   inside of all 3 edges, returns TRUE, else FALSE.
 */
int
point_in_stri(v0,v1,v2,p)
FVECT v0,v1,v2,p;
{
    double d;
    FVECT n;  

    VCROSS(n,v0,v1);
    /* Test the point for sidedness */
    d  = DOT(n,p);
    if(d > 0.0)
      return(FALSE);
    /* Test next edge */
    VCROSS(n,v1,v2);
    /* Test the point for sidedness */
    d  = DOT(n,p);
    if(d > 0.0)
       return(FALSE);
    /* Test next edge */
    VCROSS(n,v2,v0);
    /* Test the point for sidedness */
    d  = DOT(n,p);
    if(d > 0.0)
       return(FALSE);
    /* Must be interior to the pyramid */
    return(TRUE);
}

/*
 * int
 * ray_intersect_tri(orig,dir,v0,v1,v2,pt)
 *    : test if ray orig-dir intersects triangle v0v1v2, result in pt
 * FVECT orig,dir;   : Vectors defining ray origin and direction
 * FVECT v0,v1,v2;   : Triangle vertices
 * FVECT pt;         : Intersection point (if any)
 */
int
ray_intersect_tri(orig,dir,v0,v1,v2,pt)
FVECT orig,dir;
FVECT v0,v1,v2;
FVECT pt;
{
  FVECT p0,p1,p2,p;
  FPEQ peq;
  int type;

  VSUB(p0,v0,orig);
  VSUB(p1,v1,orig);
  VSUB(p2,v2,orig);

  if(point_in_stri(p0,p1,p2,dir))
  {
      /* Intersect the ray with the triangle plane */
      tri_plane_equation(v0,v1,v2,&peq,FALSE);
      return(intersect_ray_plane(orig,dir,peq,NULL,pt));
  }
  return(FALSE);
}

/*
 * calculate_view_frustum(vp,hv,vv,horiz,vert,near,far,fnear,ffar)
 *    : Calculate vertices defining front and rear clip rectangles of
 *      view frustum defined by vp,hv,vv,horiz,vert,near, and far and
 *      return in fnear and ffar.
 * FVECT vp,hv,vv;        : Viewpoint(vp),hv and vv are the horizontal and
 *                          vertical vectors in the view frame-magnitude is
 *                          the dimension of the front frustum face at z =1
 * double horiz,vert,near,far; : View angle horizontal and vertical(horiz,vert)
 *                              and distance to the near,far clipping planes
 * FVECT fnear[4],ffar[4];     : holds results
 * 
 */
calculate_view_frustum(vp,hv,vv,horiz,vert,near,far,fnear,ffar)
FVECT vp,hv,vv;
double horiz,vert,near,far;
FVECT fnear[4],ffar[4];
{
    double height,width;
    FVECT t,nhv,nvv,ndv;
    double w2,h2;
    /* Calculate the x and y dimensions of the near face */
    VCOPY(nhv,hv);
    VCOPY(nvv,vv);
    w2 = normalize(nhv);
    h2 = normalize(nvv);
    /* Use similar triangles to calculate the dimensions at z=near */
    width  = near*0.5*w2;
    height = near*0.5*h2;

    VCROSS(ndv,nvv,nhv);
    /* Calculate the world space points corresponding to the 4 corners
       of the front face of the view frustum
     */
    fnear[0][0] =  width*nhv[0] + height*nvv[0] + near*ndv[0] + vp[0] ;
    fnear[0][1] =  width*nhv[1] + height*nvv[1] + near*ndv[1] + vp[1];
    fnear[0][2] =  width*nhv[2] + height*nvv[2] + near*ndv[2] + vp[2];
    fnear[1][0] = -width*nhv[0] + height*nvv[0] + near*ndv[0] + vp[0];
    fnear[1][1] = -width*nhv[1] + height*nvv[1] + near*ndv[1] + vp[1];
    fnear[1][2] = -width*nhv[2] + height*nvv[2] + near*ndv[2] + vp[2];
    fnear[2][0] = -width*nhv[0] - height*nvv[0] + near*ndv[0] + vp[0];
    fnear[2][1] = -width*nhv[1] - height*nvv[1] + near*ndv[1] + vp[1];
    fnear[2][2] = -width*nhv[2] - height*nvv[2] + near*ndv[2] + vp[2];
    fnear[3][0] =  width*nhv[0] - height*nvv[0] + near*ndv[0] + vp[0];
    fnear[3][1] =  width*nhv[1] - height*nvv[1] + near*ndv[1] + vp[1];
    fnear[3][2] =  width*nhv[2] - height*nvv[2] + near*ndv[2] + vp[2];

    /* Now do the far face */
    width  = far*0.5*w2;
    height = far*0.5*h2;
    ffar[0][0] =  width*nhv[0] + height*nvv[0] + far*ndv[0] + vp[0] ;
    ffar[0][1] =  width*nhv[1] + height*nvv[1] + far*ndv[1] + vp[1] ;
    ffar[0][2] =  width*nhv[2] + height*nvv[2] + far*ndv[2] + vp[2] ;
    ffar[1][0] = -width*nhv[0] + height*nvv[0] + far*ndv[0] + vp[0] ;
    ffar[1][1] = -width*nhv[1] + height*nvv[1] + far*ndv[1] + vp[1] ;
    ffar[1][2] = -width*nhv[2] + height*nvv[2] + far*ndv[2] + vp[2] ;
    ffar[2][0] = -width*nhv[0] - height*nvv[0] + far*ndv[0] + vp[0] ;
    ffar[2][1] = -width*nhv[1] - height*nvv[1] + far*ndv[1] + vp[1] ;
    ffar[2][2] = -width*nhv[2] - height*nvv[2] + far*ndv[2] + vp[2] ;
    ffar[3][0] =  width*nhv[0] - height*nvv[0] + far*ndv[0] + vp[0] ;
    ffar[3][1] =  width*nhv[1] - height*nvv[1] + far*ndv[1] + vp[1] ;
    ffar[3][2] =  width*nhv[2] - height*nvv[2] + far*ndv[2] + vp[2] ;
}


/*
 * bary2d(x1,y1,x2,y2,x3,y3,px,py,coord)
 *  : Find the normalized barycentric coordinates of p relative to
 *    triangle v0,v1,v2. Return result in coord
 * double x1,y1,x2,y2,x3,y3; : defines triangle vertices 1,2,3
 * double px,py;             : coordinates of pt
 * double coord[3];          : result
 */
bary2d(x1,y1,x2,y2,x3,y3,px,py,coord)
double x1,y1,x2,y2,x3,y3;
double px,py;
double coord[3];
{
  double a;

  a =  (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
  coord[0] = ((x2 - px) * (y3 - py) - (x3 - px) * (y2 - py)) / a; 
  coord[1] = ((x3 - px) * (y1 - py) - (x1 - px) * (y3 - py)) / a;
  coord[2] = ((x1 - px) * (y2 - py) - (x2 - px) * (y1 - py)) / a;
 
}















