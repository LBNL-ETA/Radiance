#ifndef lint
static const char	RCSid[] = "$Id: sm_ogl.c,v 3.18 2003/07/21 22:30:18 schorsch Exp $";
#endif
/*
 * sm_ogl.c
 * 
 *  Rendering routines for triangle mesh representation utilizing OpenGL 
 *
 * smClean(tmflag)	: display has been wiped clean
 *     int tmflag;
 * Called after display has been effectively cleared, meaning that all
 * geometry must be resent down the pipeline in the next call to smUpdate().
 * If tmflag is set, tone-mapping should be performed
 *
 * smUpdate(vp, qua)	: update OpenGL output geometry for view vp
 * VIEW	*vp;		: desired view
 * int	qua;		: quality level (percentage on linear time scale)
 *
 * Draw new geometric representation using OpenGL calls.  Assume that the
 * view has already been set up and the correct frame buffer has been
 * selected for drawing.  The quality level is on a linear scale, where 100%
 * is full (final) quality.  It is not necessary to redraw geometry that has
 * been output since the last call to smClean().  (The last view drawn will
 * be vp==&odev.v each time.)
 */
#include "standard.h"

#include <GL/gl.h>

#include "sm_flag.h"
#include "sm_list.h"
#include "sm_geom.h"
#include "sm.h"

int smClean_notify = TRUE;           /*If TRUE:Do full redraw on next update*/
static int smCompute_mapping = TRUE; /*If TRUE:re-tonemap on next update */ 
static int smIncremental = FALSE;    /*If TRUE: there has been incremental 
				        rendering since last full draw */
#define NEW_TRI_CNT 1000             /* Default number of tris to allocate 
					space to sort for incremental update*/
#define SM_RENDER_FG 0               /* Render foreground tris only*/
#define SM_RENDER_BG 1               /* Render background tris only */
#define SM_RENDER_CULL 8             /* Perform view frustum culling */
#define BASE 1                       /* Indicates base triangle */
#define DIR 2                        /* Indicates triangle w/directional pts*/
/* FOR DISPLAY LIST RENDERING: **********************************************/
#define  SM_DL_LEVELS 2     /* # of levels down to create display lists */
#define SM_DL_LISTS   42    /* # of qtree nodes in tree at above level: 
			       should be 2*(4^(SM_DL_LEVELS+1)-1)/(4-1) */
static GLuint Display_lists[SM_DL_LISTS][2] = {0}; 
/****************************************************************************/

/* FOR APPROXIMATION RENDERING **********************************************/
typedef struct {
	float	dist;		/* average distance */
	BYTE	rgb[3];		/* average color */
} QTRAVG;		/* average quadtree value */

typedef struct {
	QUADTREE	qt;	/* quadtree node (key & hash value) */
	QTRAVG		av;	/* node average */
} QT_LUENT;		/* lookup table entry */

static QT_LUENT	*qt_htbl = NULL;	/* quadtree cache */
static int	qt_hsiz = 0;		/* quadtree cache size */
/****************************************************************************/

/* For DEPTH SORTING ********************************************************/
typedef struct _T_DEPTH {
  int tri;
  double depth;
}T_DEPTH;
/**********************************************************************/

 /*
  * smClean(tmflag)	: display has been wiped clean
  *     int tmflag;
  * Called after display has been effectively cleared, meaning that all
  * geometry must be resent down the pipeline in the next call to smUpdate().
  * If tmflag is set, tone-mapping should be performed
  */
smClean(tmflag)
   int tmflag;
{
    smClean_notify = TRUE;
    smCompute_mapping = tmflag;
}

int
qtCache_init(nel)		/* initialize for at least nel elements */
int	nel;
{
    static int  hsiztab[] = {
	8191, 16381, 32749, 65521, 131071, 262139, 524287, 1048573, 0
    };
    register int  i;

    if (nel <= 0) {			/* call to free table */
	if (qt_hsiz) {
	    free((void *)qt_htbl);
	    qt_htbl = NULL;
	    qt_hsiz = 0;
	}
	return(0);
    }
    nel += nel>>1;			/* 66% occupancy */
    for (i = 0; hsiztab[i]; i++)
       if (hsiztab[i] > nel)
	  break;
    if (!(qt_hsiz = hsiztab[i]))
       qt_hsiz = nel*2 + 1;		/* not always prime */
    qt_htbl = (QT_LUENT *)calloc(qt_hsiz, sizeof(QT_LUENT));
    if (qt_htbl == NULL)
       qt_hsiz = 0;
    for (i = qt_hsiz; i--; )
       qt_htbl[i].qt = EMPTY;
    return(qt_hsiz);
}

QT_LUENT *
qtCache_find(qt)		/* find a quadtree table entry */
QUADTREE qt;
{
    int	i, n;
    register int	ndx;
    register QT_LUENT	*le;

    if (qt_hsiz == 0 && !qtCache_init(1))
       return(NULL);
tryagain:				/* hash table lookup */
    ndx = (unsigned long)qt % qt_hsiz;
    for (i = 0, n = 1; i < qt_hsiz; i++, n += 2) {
	le = &qt_htbl[ndx];
	if (QT_IS_EMPTY(le->qt) || le->qt == qt)
	   return(le);
	if ((ndx += n) >= qt_hsiz)	/* this happens rarely */
	   ndx = ndx % qt_hsiz;
    }
    /* table is full, reallocate */
    le = qt_htbl;
    ndx = qt_hsiz;
    if (!qtCache_init(ndx+1)) {	/* no more memory! */
	qt_htbl = le;
	qt_hsiz = ndx;
	return(NULL);
    }
					/* copy old table to new and free */
    while (ndx--)
       if (!QT_IS_EMPTY(le[ndx].qt))
	  *qtCache_find(le[ndx].qt) = le[ndx];
    free((void *)le);
    goto tryagain;			/* should happen only once! */
}

stCount_level_leaves(lcnt, qt)	/* count quadtree leaf nodes at each level */
int lcnt[];
register QUADTREE qt;
{
  if (QT_IS_EMPTY(qt))
    return;
  if (QT_IS_TREE(qt)) {
    if (!QT_IS_FLAG(qt))	/* not in our frustum */
      return;
    stCount_level_leaves(lcnt+1, QT_NTH_CHILD(qt,0));
    stCount_level_leaves(lcnt+1, QT_NTH_CHILD(qt,1));
    stCount_level_leaves(lcnt+1, QT_NTH_CHILD(qt,2));
    stCount_level_leaves(lcnt+1, QT_NTH_CHILD(qt,3));
  }
  else
    if(QT_LEAF_IS_FLAG(qt))
      lcnt[0]++;
}

QTRAVG *
qtRender_level(qt,v0,v1,v2,sm,lvl)
QUADTREE qt;
FVECT v0,v1,v2;
SM *sm;
int lvl;
{
  FVECT a,b,c;
  register QT_LUENT *le;
  QTRAVG *rc[4];
  TRI *tri;
  
  if (QT_IS_EMPTY(qt))				/* empty leaf node */
    return(NULL);
  if (QT_IS_TREE(qt) && !QT_IS_FLAG(qt))	/* not in our frustum */
    return(NULL);
  if(QT_IS_LEAF(qt)  && !QT_LEAF_IS_FLAG(qt))	/* not in our frustum */
    return(NULL);
					/* else look up node */
  if ((le = qtCache_find(qt)) == NULL)
    error(SYSTEM, "out of memory in qtRender_level");
  if (QT_IS_TREE(qt) && (QT_IS_EMPTY(le->qt) || lvl > 0))
  {					/* compute children */
    qtSubdivide_tri(v0,v1,v2,a,b,c);
    rc[0] = qtRender_level(QT_NTH_CHILD(qt,0),v0,c,b,sm,lvl-1);
    rc[1] = qtRender_level(QT_NTH_CHILD(qt,1),c,v1,a,sm,lvl-1);
    rc[2] = qtRender_level(QT_NTH_CHILD(qt,2),b,a,v2,sm,lvl-1);
    rc[3] = qtRender_level(QT_NTH_CHILD(qt,3),a,b,c,sm,lvl-1);
  }
  if (QT_IS_EMPTY(le->qt))
  {					/* let's make some data! */
    int rgbs[3];
    double distsum;
    register int i, n;
					/* average our triangle vertices */
    rgbs[0] = rgbs[1] = rgbs[2] = 0;
    distsum = 0.; n = 0;
    if(QT_IS_TREE(qt))
    {					/* from subtree */
      for (i = 4; i--; )
	if (rc[i] != NULL)
	{
	  rgbs[0] += rc[i]->rgb[0]; rgbs[1] += rc[i]->rgb[1];
	  rgbs[2] += rc[i]->rgb[2]; distsum += rc[i]->dist; n++;
	}
    }
    else
    {					/* from triangle set */
      S_ID *os;
      S_ID s0, s1, s2,s_id;
      int t_id;
      TRI *tri,*t;

      os = qtqueryset(qt);
      for (i = os[0]; i; i--)
      {
	s_id = os[i];
	t_id = SM_NTH_VERT(smMesh,s_id);
	tri = t = SM_NTH_TRI(smMesh,t_id);
	do
	{
	  if(!SM_IS_NTH_T_BASE(sm,t_id))
	  {
	    n++;
	    s0 = T_NTH_V(t,0);
	    s1 = T_NTH_V(t,1);
	    s2 = T_NTH_V(t,2);
	    VCOPY(a,SM_NTH_WV(sm,s0));
	    VCOPY(b,SM_NTH_WV(sm,s1));
	    VCOPY(c,SM_NTH_WV(sm,s2));	      
	    distsum += SM_BG_SAMPLE(sm,s0) ? dev_zmax
				: sqrt(dist2(a,SM_VIEW_CENTER(sm)));
	    distsum += SM_BG_SAMPLE(sm,s1) ? dev_zmax
	                        : sqrt(dist2(b,SM_VIEW_CENTER(sm)));
	    distsum += SM_BG_SAMPLE(sm,s2) ? dev_zmax
				: sqrt(dist2(c,SM_VIEW_CENTER(sm)));
	    rgbs[0] += SM_NTH_RGB(sm,s0)[0] + SM_NTH_RGB(sm,s1)[0]
	      + SM_NTH_RGB(sm,s2)[0];
	    rgbs[1] += SM_NTH_RGB(sm,s0)[1] + SM_NTH_RGB(sm,s1)[1]
	      + SM_NTH_RGB(sm,s2)[1];
	    rgbs[2] += SM_NTH_RGB(sm,s0)[2] + SM_NTH_RGB(sm,s1)[2]
	      + SM_NTH_RGB(sm,s2)[2];
	  }

	  t_id = smTri_next_ccw_nbr(smMesh,t,s_id);
	  t = SM_NTH_TRI(smMesh,t_id);

	}while(t != tri);
      }
      n *= 3;
    }
    if (!n)
      return(NULL);
    le->qt = qt;
    le->av.rgb[0] = rgbs[0]/n; le->av.rgb[1] = rgbs[1]/n;
    le->av.rgb[2] = rgbs[2]/n; le->av.dist = distsum/(double)n;
  }
  if (lvl == 0 || (lvl > 0 && QT_IS_LEAF(qt)))
  {				/* render this node */
					/* compute pseudo vertices */
    VCOPY(a,v0); VCOPY(b,v1); VCOPY(c,v2);
    normalize(a); normalize(b); normalize(c);
    VSUM(a,SM_VIEW_CENTER(sm),a,le->av.dist);
    VSUM(b,SM_VIEW_CENTER(sm),b,le->av.dist);
    VSUM(c,SM_VIEW_CENTER(sm),c,le->av.dist);
					/* draw triangle */
    glColor3ub(le->av.rgb[0],le->av.rgb[1],le->av.rgb[2]);
    glVertex3d(a[0],a[1],a[2]);
    glVertex3d(b[0],b[1],b[2]);
    glVertex3d(c[0],c[1],c[2]);

  }
  return(&le->av);
}


smRender_approx_stree_level(sm,lvl)
SM *sm;
int lvl;
{
  QUADTREE qt;
  int i;
  FVECT t0,t1,t2;
  STREE *st;
  
  if (lvl < 0)
    return;
  st = SM_LOCATOR(sm);
  glPushAttrib(GL_LIGHTING_BIT);
  glShadeModel(GL_FLAT);
  glBegin(GL_TRIANGLES);
  for(i=0; i < ST_NUM_ROOT_NODES; i++)
  {
    qt = ST_ROOT_QT(st,i);
    qtRender_level(qt,ST_NTH_V(st,i,0),ST_NTH_V(st,i,1),ST_NTH_V(st,i,2),
		   sm,lvl);
  }
  glEnd();
  glPopAttrib();
}

/*
 * smRender_approx(sm,qual,view)
 * SM *sm;                : mesh
 * int qual;              : quality level
 * VIEW *view;            : current view
 *
 *  Renders an approximation to the current mesh based on the quadtree 
 *  subdivision. The quadtree is traversed to a level (based upon the quality:
 *  the lower the quality, the fewer levels visited, and the coarser, and
 *  faster, the approximation). The quadtree triangles are drawn relative to
 *  the current viewpoint, with a depth and color averaged from all of the 
 *  triangles that lie beneath the node. 
 */
smRender_approx(sm, qual,view)	
SM *sm;
int qual;
VIEW *view;
{
  int i, n,ntarget;
  int lvlcnt[QT_MAX_LEVELS];
  STREE *st;
  int32 *active_flag;

  if (qual <= 0)
    return;
  smCull(sm,view,SM_ALL_LEVELS);
				/* compute rendering target */
  ntarget = 0;

  active_flag = SM_NTH_FLAGS(sm,T_ACTIVE_FLAG);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(active_flag[n])
      for(i=0; i < 32; i++)
	if(active_flag[n] & (1L << i))
	  ntarget++;

  ntarget = ntarget*qual/MAXQUALITY;
  if (!ntarget)
    return;
  for (i = QT_MAX_LEVELS; i--; )
    lvlcnt[i] = 0;
  
  st = SM_LOCATOR(sm);
  for(i=0; i < ST_NUM_ROOT_NODES;i++)
    stCount_level_leaves(lvlcnt, ST_ROOT_QT(st,i));
  
  for (ntarget -= lvlcnt[i=0]; i < QT_MAX_LEVELS-1; ntarget -= lvlcnt[++i])
    if (ntarget < lvlcnt[i+1])
      break;
				/* compute and render target level */
  smRender_approx_stree_level(sm,i);
}


#define GLVERTEX3V(v) glVertex3fv(v)


#define render_tri(v0,v1,v2,rgb0,rgb1,rgb2) \
  {glColor3ub(rgb0[0],rgb0[1],rgb0[2]); GLVERTEX3V(v0); \
  glColor3ub(rgb1[0],rgb1[1],rgb1[2]);  GLVERTEX3V(v1); \
  glColor3ub(rgb2[0],rgb2[1],rgb2[2]);  GLVERTEX3V(v2);} 

/*
 * render_mixed_tri(v0,v1,v2,rgb0,rgb1,rgb2,b0,b1,b2)
 *  SFLOAT v0[3],v1[3],v2[3];      : triangle vertex coordinates
 *  BYTE rgb0[3],rgb1[3],rgb2[3]; : vertex RGBs
 *  int b0,b1,b2;                 : background or base vertex flag
 *  
 *  For triangles with one or more base or directional vertices.
 *  render base vertex color as average of the background  and foreground
 *  vertex RGBs. The coordinates for a fg vertex are calculated by
 *  subtracting off the current view,normalizing, then scaling to fit
 * into the current frustum.
 */
render_mixed_tri(v0,v1,v2,rgb0,rgb1,rgb2,vp,vc,bg0,bg1,bg2)
SFLOAT v0[3],v1[3],v2[3];
BYTE rgb0[3],rgb1[3],rgb2[3];
FVECT vp,vc;
int bg0,bg1,bg2;
{
  double d,p[3];
  int j,cnt,rgb[3];
  
  /* Average color from bg vertices */
  cnt = 0;
  if(bg0 == BASE || bg1==BASE || bg2 == BASE)
  {
    rgb[0] = rgb[1] = rgb[2] = 0;
    if(bg0 != BASE)
      {
	IADDV3(rgb,rgb0);
	cnt++;
      }
    if(bg1 != BASE)
      {
	IADDV3(rgb,rgb1);
	cnt++;
      }
    if(bg2 != BASE)
      {
	IADDV3(rgb,rgb2);
	cnt++;
      }
    IDIVV3(rgb,cnt);
  }
  if(bg0 == BASE)
    glColor3i(rgb[0],rgb[1],rgb[2]);
  else
    glColor3ub(rgb0[0],rgb0[1],rgb0[2]);

  if(!bg0)
  {
    VSUB(p,v0,vp);
    normalize(p);
    IADDV3(p,vc);
    glVertex3dv(p); 
  }
  else
    GLVERTEX3V(v0); 

  if(bg1 == BASE)
    glColor3i(rgb[0],rgb[1],rgb[2]);
  else
    glColor3ub(rgb1[0],rgb1[1],rgb1[2]);

  if(!bg1)
  {
    VSUB(p,v1,vp);
    normalize(p);
    IADDV3(p,vc);
    glVertex3dv(p); 
  }
  else
    GLVERTEX3V(v1); 

  if(bg2 == BASE)
    glColor3i(rgb[0],rgb[1],rgb[2]);
  else
    glColor3ub(rgb2[0],rgb2[1],rgb2[2]);

  if(!bg2)
  {
    VSUB(p,v2,vp);
    normalize(p);
    IADDV3(p,vc);
    glVertex3dv(p); 
  }
  else
    GLVERTEX3V(v2); 
}

/*
 * smRender_bg_tris(sm,vp,t_flag,bg_flag,wp,rgb)
 * SM *sm;                         : mesh
 * FVECT vp;                       : current viewpoint
 * int32  *t_flag,*bg_flag;         : triangle flags: t_flag is generic, 
 *                                   and bg_flag indicates if background tri;
 * SFLOAT (*wp)[3];BYTE  (*rgb)[3]; : arrays of sample points and RGB colors
 *
 * Sequentially traverses triangle list and renders all valid tris who 
 * have t_flag set, and bg_flag set.
 */

smRender_bg_tris(sm,vp,t_flag,bg_flag,wp,rgb)
SM *sm;
FVECT vp;
int32 *t_flag,*bg_flag;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
{
  double d;
  int v0_id,v1_id,v2_id;
  int i,n,bg0,bg1,bg2;
  TRI *tri;

  glMatrixMode(GL_MODELVIEW);

  glPushMatrix();
  glTranslated(vp[0],vp[1],vp[2]);
  /* The points are a distance of 1 away from the origin: if necessary scale
     so that they fit in frustum and are therefore not clipped away
   */
  if(dev_zmin >= 0.99)
  {
    d = (dev_zmin+dev_zmax)/2.0;
    glScaled(d,d,d);
  }
  /* move relative to the new view */
  /* move points to unit sphere at origin */
  glTranslated(-SM_VIEW_CENTER(sm)[0],-SM_VIEW_CENTER(sm)[1],
	       -SM_VIEW_CENTER(sm)[2]);
  glBegin(GL_TRIANGLES);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(t_flag[n] & bg_flag[n])
      for(i=0; i < 32; i++)
	if(t_flag[n] & bg_flag[n] & (1L << i))
	 {
	   tri = SM_NTH_TRI(sm,(n<<5)+i);
	   v0_id = T_NTH_V(tri,0);
	   v1_id = T_NTH_V(tri,1);
	   v2_id = T_NTH_V(tri,2);
	   bg0 = SM_DIR_ID(sm,v0_id)?DIR:SM_BASE_ID(sm,v0_id)?BASE:0;
	   bg1 = SM_DIR_ID(sm,v1_id)?DIR:SM_BASE_ID(sm,v1_id)?BASE:0;
	   bg2 = SM_DIR_ID(sm,v2_id)?DIR:SM_BASE_ID(sm,v2_id)?BASE:0;
	   if(bg0==DIR && bg1==DIR && bg2==DIR)
	     render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
		rgb[v2_id])
	   else
	     render_mixed_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
	      rgb[v1_id],rgb[v2_id],vp,SM_VIEW_CENTER(sm),bg0,bg1,bg2);
	 }
  glEnd();

  glPopMatrix();

}
/*
 * smRender_new_bg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb)
 * SM *sm;                         : mesh
 * FVECT vp;                       : current viewpoint
 * int32  *new_flag,*active,*bg_flag; : triangle flags: idicates if tri is, 
 *                                   new,active,bg_flag 
 * SFLOAT (*wp)[3];BYTE  (*rgb)[3]; : arrays of sample points and RGB colors
 *
 * Sequentially traverses triangle list and renders all valid tris who 
 * have new_flag set and active_flag and bg_flag set.
 */
smRender_new_bg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb)
SM *sm;
FVECT vp;
int32 *new_flag,*active_flag,*bg_flag;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
{
  double d;
  int v0_id,v1_id,v2_id;
  int i,n,bg0,bg1,bg2;
  TRI *tri;

  glMatrixMode(GL_MODELVIEW);

  glPushMatrix();
  glTranslated(vp[0],vp[1],vp[2]);
  /* The points are a distance of 1 away from the origin: if necessary scale
     so that they fit in frustum and are therefore not clipped away
   */
  if(dev_zmin >= 0.99)
  {
    d = (dev_zmin+dev_zmax)/2.0;
    glScaled(d,d,d);
  }
  /* move relative to the new view */
  /* move points to unit sphere at origin */
  glTranslated(-SM_VIEW_CENTER(sm)[0],-SM_VIEW_CENTER(sm)[1],
	       -SM_VIEW_CENTER(sm)[2]);
  glBegin(GL_TRIANGLES);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(new_flag[n] & active_flag[n] & bg_flag[n])
      for(i=0; i < 32; i++)
	if(new_flag[n] & active_flag[n] & bg_flag[n] & (1L << i))
	 {
	   tri = SM_NTH_TRI(sm,(n<<5)+i);
	   v0_id = T_NTH_V(tri,0);
	   v1_id = T_NTH_V(tri,1);
	   v2_id = T_NTH_V(tri,2);
	   bg0 = SM_DIR_ID(sm,v0_id)?DIR:SM_BASE_ID(sm,v0_id)?BASE:0;
	   bg1 = SM_DIR_ID(sm,v1_id)?DIR:SM_BASE_ID(sm,v1_id)?BASE:0;
	   bg2 = SM_DIR_ID(sm,v2_id)?DIR:SM_BASE_ID(sm,v2_id)?BASE:0;
	   if(bg0==DIR && bg1==DIR && bg2==DIR)
	     render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
		rgb[v2_id])
	   else
	     render_mixed_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
	      rgb[v1_id],rgb[v2_id],vp,SM_VIEW_CENTER(sm),bg0,bg1,bg2);
	 }
  glEnd();

  glPopMatrix();

}

/*
 * render_base_tri(v0,v1,v2,rgb0,rgb1,rgb2,vp,b0,b1,b2)
 *  SFLOAT v0[3],v1[3],v2[3];       : triangle vertex coordinates
 *  BYTE rgb0[3],rgb1[3],rgb2[3];  : vertex RGBs
 *  FVECT vp;                      : current viewpoint
 *  int b0,b1,b2;                  : vertex base flag
 *  
 *  render base vertex color as average of the non-base vertex RGBs. The
 *  base vertex coordinate is taken as the stored vector, scaled out by
 *  the average distance to the non-base vertices
 */
render_base_tri(v0,v1,v2,rgb0,rgb1,rgb2,vp,b0,b1,b2)
SFLOAT v0[3],v1[3],v2[3];
BYTE rgb0[3],rgb1[3],rgb2[3];
FVECT vp;
int b0,b1,b2;
{
  int cnt;
  int rgb[3];
  double d;
  double p[3];

  cnt = 0;
  rgb[0] = rgb[1] = rgb[2] = 0;
  d = 0.0;
  /* If all vertices are base: don't render */
  if(b0&&b1&&b2)
    return;
  /* First calculate color and coordinates 
     for base vertices based on world space vertices*/
  if(!b0)
  {
    IADDV3(rgb,rgb0);
    d += DIST(v0,vp);
    cnt++;
  }
  if(!b1)
  {
    IADDV3(rgb,rgb1);
    d += DIST(v1,vp);
    cnt++;
  }
  if(!b2)
  {
    IADDV3(rgb,rgb2);
    d += DIST(v2,vp);
    cnt++;
  }
  IDIVV3(rgb,cnt);
  d /= (double)cnt;
  
  /* Now render triangle */
  if(b0)
  {
    glColor3ub(rgb[0],rgb[1],rgb[2]);
    SUBV3(p,v0,vp);
    ISCALEV3(p,d);
    IADDV3(p,vp);
    glVertex3dv(p);    
  }
  else
  {
    glColor3ub(rgb0[0],rgb0[1],rgb0[2]);
    GLVERTEX3V(v0);
  }
  if(b1)
  {
    glColor3ub(rgb[0],rgb[1],rgb[2]);
    SUBV3(p,v1,vp);
    ISCALEV3(p,d);
    IADDV3(p,vp);
    glVertex3dv(p);
  }
  else
  {
    glColor3ub(rgb1[0],rgb1[1],rgb1[2]);
    GLVERTEX3V(v1);
  }
  if(b2)
  {
    glColor3ub(rgb[0],rgb[1],rgb[2]);
    SUBV3(p,v2,vp);
    ISCALEV3(p,d);
    IADDV3(p,vp);
    glVertex3dv(p);
  }
  else
  {
    glColor3ub(rgb2[0],rgb2[1],rgb2[2]);
    GLVERTEX3V(v2);
  }
}
/*
 * smRender_fg_tris(sm,vp,t_flag,bg_flag,wp,rgb)
 * SM *sm;                        : mesh
 * FVECT vp;                      : current viewpoint
 * int32  *t_flag,*bg_flag;        : triangle flags: t_flag is generic,bg_flag 
 *                                  indicates if background tri;
 * SFLOAT (*wp)[3];BYTE (*rgb)[3]; : arrays of sample points and RGB colors
 *
 * Sequentially gos through triangle list and renders all valid tris who 
 * have t_flag set, and NOT bg_flag set.
 */
smRender_fg_tris(sm,vp,t_flag,bg_flag,wp,rgb)
SM *sm;
FVECT vp;
int32  *t_flag,*bg_flag;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
{
  TRI *tri;
  int i,n,b0,b1,b2;
  S_ID v0_id,v1_id,v2_id;
  
  glBegin(GL_TRIANGLES);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(t_flag[n])
      for(i=0; i < 32; i++)
	if(t_flag[n] & (1L << i) & ~bg_flag[n])
	 {
	   tri = SM_NTH_TRI(sm,(n<<5)+i);
	   v0_id = T_NTH_V(tri,0);
	   v1_id = T_NTH_V(tri,1);
	   v2_id = T_NTH_V(tri,2);
	   b0 = SM_BASE_ID(sm,v0_id);
	   b1 = SM_BASE_ID(sm,v1_id);
	   b2 = SM_BASE_ID(sm,v2_id);
	   if(b0 || b1 || b2)
	     render_base_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
	     rgb[v1_id],rgb[v2_id],SM_VIEW_CENTER(sm),b0,b1,b2);
	   else
	     render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
			rgb[v2_id])

	 }
  glEnd();

}

/*
 * smRender_new_fg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb)
 * SM *sm;                        : mesh
 * FVECT vp;                      : current viewpoint
 * int32  *new_flag,*active_flag,*bg_flag; : triangle flags: indicate if
 *                                tri is new,active,background
 * SFLOAT (*wp)[3];BYTE (*rgb)[3]; : arrays of sample points and RGB colors
 *
 * Sequentially gos through triangle list and renders all valid tris who 
 * have new_flag and active_flag set, and NOT bg_flag set.
 */
smRender_new_fg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb)
SM *sm;
FVECT vp;
int32  *new_flag,*active_flag,*bg_flag;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
{
  TRI *tri;
  int i,n,b0,b1,b2;
  S_ID v0_id,v1_id,v2_id;
  
  glBegin(GL_TRIANGLES);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(new_flag[n] & active_flag[n])
      for(i=0; i < 32; i++)
	if(new_flag[n] & active_flag[n] & (1L << i) & ~bg_flag[n])
	 {
	   tri = SM_NTH_TRI(sm,(n<<5)+i);
	   v0_id = T_NTH_V(tri,0);
	   v1_id = T_NTH_V(tri,1);
	   v2_id = T_NTH_V(tri,2);
	   b0 = SM_BASE_ID(sm,v0_id);
	   b1 = SM_BASE_ID(sm,v1_id);
	   b2 = SM_BASE_ID(sm,v2_id);
	   if(b0 || b1 || b2)
	     render_base_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
	     rgb[v1_id],rgb[v2_id],SM_VIEW_CENTER(sm),b0,b1,b2);
	   else
	     render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
			rgb[v2_id])
	 }
  glEnd();

}

/* Test for qsort to depth sort triangles */
int
compare_tri_depths(T_DEPTH *td1,T_DEPTH *td2)
{
  double d;

  d = td2->depth-td1->depth;
  if(d > 0.0)
    return(1);
  if(d < 0.0)
    return(-1);
  return(0);
}


/* 
 * smOrder_new_tris(sm,vp,td)
 * SM *sm;      : mesh 
 * FVECT vp;    : current viewpoint
 * T_DEPTH *td; : holds returned list of depth sorted tris
 *
 * Creates list of all new tris, with their distance from the current 
 * viewpoint, and sorts the list based on this distance
 */
T_DEPTH
*smOrder_new_tris(sm,vp)
SM *sm;
FVECT vp;
{
  T_DEPTH *td;
  int n,i,j,tcnt,v,size;
  TRI *tri;
  double d,min_d;
  FVECT diff;
  int32 *new_flag,*bg_flag,*active_flag;
  
 td = (T_DEPTH *)tempbuf(NEW_TRI_CNT*sizeof(T_DEPTH),FALSE);
 size = NEW_TRI_CNT;

  tcnt=0;
  new_flag = SM_NTH_FLAGS(sm,T_NEW_FLAG);
  bg_flag = SM_NTH_FLAGS(sm,T_BG_FLAG);
  active_flag = SM_NTH_FLAGS(sm,T_ACTIVE_FLAG);
  for(n=((SM_NUM_TRI(sm)+31)>>5) +1; --n;)
    if(active_flag[n] & new_flag[n] & ~bg_flag[n])
      for(i=0; i < 32; i++)
	if(active_flag[n] & new_flag[n] & (1L << i) & ~bg_flag[n])
	 {
	   tri = SM_NTH_TRI(sm,(n<<5)+i);
	   if(tcnt+1 >= size)
	   {
	     size += 100;
	     td = (T_DEPTH *)tempbuf(size*sizeof(T_DEPTH),TRUE);
	   }
	   td[tcnt].tri = (n << 5)+i;
	   min_d = -1;
	   for(j=0;j < 3;j++)
	     {
	       v = T_NTH_V(tri,j);
	       VSUB(diff,SM_NTH_WV(sm,v),vp);
	       d = DOT(diff,diff);
	       if(min_d == -1 || d < min_d)
		 min_d = d;
	     }
	   td[tcnt++].depth = min_d;
  }
  td[tcnt].tri = -1;
  if(tcnt)
      qsort((void *)td,tcnt,sizeof(T_DEPTH),compare_tri_depths);
  return(td);
}

/* 
 * smUpdate_tm(sm)   : Update the tone-mapping
 * SM *sm;           : mesh
 * 
 */
smUpdate_tm(sm)
SM *sm;
{
  int t = SM_TONE_MAP(sm);

  if(t==0 || smCompute_mapping)
  {
    tmClearHisto();
    tmAddHisto(SM_BRT(sm),SM_NUM_SAMP(sm),1);
    if(tmComputeMapping(0.,0.,0.) != TM_E_OK)
      return;
    t = 0;
    smCompute_mapping = FALSE;
  }
  tmMapPixels(SM_NTH_RGB(sm,t),&SM_NTH_BRT(sm,t),SM_NTH_CHR(sm,t), 
	      SM_NUM_SAMP(sm)-t);
  SM_TONE_MAP(sm) = SM_NUM_SAMP(sm);
}

/* 
 *   smRender_inc(sm,vp)      : Incremental update of mesh
 *    SM * sm;                : mesh
 *    FVECT vp;               : current view point
 *
 *  If a relatively small number of new triangles have been created,
 *  do an incremental update. Render new triangles with depth buffering 
 *  turned off, if the current viewpoint is not the same as canonical view- 
 *  point, must use painter's approach to resolve visibility:first depth sort
 *  triangles, then render back-to-front.  
*/
smRender_inc(sm,vp)
SM *sm;
FVECT vp;
{
  S_ID v0_id,v1_id,v2_id;
  int i,n,b0,b1,b2;
  TRI *tri;
  SFLOAT (*wp)[3];
  BYTE  (*rgb)[3];
  int32  *new_flag,*bg_flag,*active_flag;
  T_DEPTH *td = NULL;


  /* For all of the NEW and ACTIVE triangles (since last update): 
     Go through and sort on depth value (from vp). Turn
     Depth Buffer test off and render back-front
   */

  /* Must depth sort if current view point is not same as canonical */
  if(!EQUAL_VEC3(SM_VIEW_CENTER(sm),vp))
    td =  smOrder_new_tris(sm,vp);
  wp = SM_WP(sm);
  rgb =SM_RGB(sm);
  new_flag = SM_NTH_FLAGS(sm,T_NEW_FLAG);
  active_flag = SM_NTH_FLAGS(sm,T_ACTIVE_FLAG);
  bg_flag = SM_NTH_FLAGS(sm,T_BG_FLAG);
  /* Turn Depth Test off -- using Painter's algorithm */
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
  glDepthFunc(GL_ALWAYS);

  smRender_new_bg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb);
  if(!td)
    smRender_new_fg_tris(sm,vp,new_flag,active_flag,bg_flag,wp,rgb);
  else
  {
    glBegin(GL_TRIANGLES);
    for(i=0; td[i].tri != -1;i++)
    {
      tri = SM_NTH_TRI(sm,td[i].tri);
      /* Dont need to check for valid tri because flags are
	 cleared on delete
	 */
      v0_id = T_NTH_V(tri,0);
      v1_id = T_NTH_V(tri,1);
      v2_id = T_NTH_V(tri,2);
      b0 = SM_BASE_ID(sm,v0_id);
      b1 = SM_BASE_ID(sm,v1_id);
      b2 = SM_BASE_ID(sm,v2_id);
      if(b0 || b1 || b2)
	render_base_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
			rgb[v1_id],rgb[v2_id],SM_VIEW_CENTER(sm),b0,b1,b2);
      else
	render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
		   rgb[v2_id])
    }
    glEnd();
    freebuf(td);
  }
  /* Restore Depth Test */
  glPopAttrib();
}

/* 
 * smRender_qtree_dl(sm,qt,vp,wp,rgb,i,level_i,max_level,leaf_cnt,which)
 *  SM *sm;              : mesh
 *  QUADTREE qt;         : quadtree base node
 *  FVECT vp;            : current viewpoint
 *  SFLOAT (*wp)[3];      : array of sample points
 *  BYTE (*rgb)[3];      : array of RGB values for samples
 *  int i,level_i,level,max_level,leaf_cnt; 
 *          : variables to keep track of where
 *         we are in the quadtree traversal in order to map nodes to 
 *         corresponding array locations, where nodes are stored in breadth-
 *         first order. i is the index of the current node,level_i is the 
 *         index of the first node on the current quadtree level, max_level is
 *         the maximum number of levels to traverse, and leaf_cnt is the number
 *         of leaves on the current level
 *  int which;          flag indicates whether to render fg or bg tris
 *
 * 
 *    Render the tris stored in qtree using display lists. For each node at
 *   the leaf or max_level, call the display_list if it exists, else traverse
 *   down the subtree and render the nodes into a new display list which is
 *   stored for future use.
 */
smRender_qtree_dl(sm,qt,vp,wp,rgb,i,level_i,level,max_level,leaf_cnt,which)
SM *sm;
QUADTREE qt;
FVECT vp;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
int i,level_i,level,max_level,leaf_cnt;
int which;
{
  int j;
  
  if(QT_IS_EMPTY(qt))
    return;

  if(QT_IS_LEAF(qt) || level == max_level)
  {
    if(QT_IS_LEAF(qt))
    {
      if(!QT_LEAF_IS_FLAG(qt))
	return;
    }
    else
      if(!QT_IS_FLAG(qt))
	return;

    if(!Display_lists[i][which])
    {
      Display_lists[i][which] = i+1 + which*SM_DL_LISTS;
      glNewList(Display_lists[i][which],GL_COMPILE_AND_EXECUTE);
      smClear_flags(sm,T_NEW_FLAG);
      glBegin(GL_TRIANGLES);
      smRender_qtree(sm,qt,vp,wp,rgb,which,FALSE);
      glEnd();
      glEndList();
    }
    else
    {
      glCallList(Display_lists[i][which]);
    }
  }
  else
    if(QT_IS_FLAG(qt))
    {
      i = ((i - level_i)<< 2) + level_i + leaf_cnt;
      level_i += leaf_cnt;
      leaf_cnt <<= 2;
      for(j=0; j < 4; j++)
	smRender_qtree_dl(sm,QT_NTH_CHILD(qt,j),vp,wp,rgb,
			i+j,level_i,level+1,max_level,leaf_cnt,which);
    }

}

/* 
 * smRender_qtree(sm,qt,vp,wp,rgb,which,cull) : Render the tris stored in qtree
 *  SM *sm;             : mesh
 *  QUADTREE qt;        : quadtree base node
 *  FVECT vp;           : current viewpoint
 *  SFLOAT (*wp)[3]      : array of sample points
 *  BYTE (*rgb)[3]      : array of RGB values for samples
 *  int which;          : flag indicates whether to render fg or bg tris
 *  int cull;           : if true, only traverse active (flagged) nodes
 * 
 */
smRender_qtree(sm,qt,vp,wp,rgb,which,cull)
SM *sm;
QUADTREE qt;
FVECT vp;
SFLOAT (*wp)[3];
BYTE  (*rgb)[3];
int which,cull;
{
  int i;
  
  if(QT_IS_EMPTY(qt))
    return;

  if(QT_IS_LEAF(qt))
  {
    TRI *t,*tri;
    S_ID *optr,s_id,v0_id,v1_id,v2_id;
    int bg0,bg1,bg2,t_id;

    if(cull && !QT_LEAF_IS_FLAG(qt))
      return;

    optr = qtqueryset(qt);
    for (i = QT_SET_CNT(optr);i > 0; i--)
    {
      s_id = QT_SET_NEXT_ELEM(optr);
      t_id = SM_NTH_VERT(smMesh,s_id);
      tri = t = SM_NTH_TRI(smMesh,t_id);
      do
      {
       if((!cull || SM_IS_NTH_T_ACTIVE(sm,t_id)) && !SM_IS_NTH_T_NEW(sm,t_id))
       {
	  bg0 = SM_IS_NTH_T_BG(sm,t_id);
	  if((which == SM_RENDER_FG && !bg0) || (which== SM_RENDER_BG && bg0))
	  {
	    v0_id = T_NTH_V(t,0);
	    v1_id = T_NTH_V(t,1);
	    v2_id = T_NTH_V(t,2);
	    if(bg0)
	    {
	      bg0 = SM_DIR_ID(sm,v0_id)?DIR:SM_BASE_ID(sm,v0_id)?BASE:0;
	      bg1 = SM_DIR_ID(sm,v1_id)?DIR:SM_BASE_ID(sm,v1_id)?BASE:0;
	      bg2 = SM_DIR_ID(sm,v2_id)?DIR:SM_BASE_ID(sm,v2_id)?BASE:0;
	      SM_SET_NTH_T_NEW(sm,t_id);  
	      if(bg0==DIR && bg1==DIR && bg2==DIR)
		render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
			 rgb[v2_id])
	      else
		render_mixed_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
		    rgb[v1_id],rgb[v2_id],vp,SM_VIEW_CENTER(sm),bg0,bg1,bg2);
	    }
	    else
	    {
	      SM_SET_NTH_T_NEW(sm,t_id);  
	      bg0 = SM_BASE_ID(sm,v0_id);
	      bg1 = SM_BASE_ID(sm,v1_id);
	      bg2 = SM_BASE_ID(sm,v2_id);
	      if(bg0 || bg1 || bg2)
		render_base_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],
			rgb[v1_id],rgb[v2_id],SM_VIEW_CENTER(sm),bg0,bg1,bg2);
	      else
		render_tri(wp[v0_id],wp[v1_id],wp[v2_id],rgb[v0_id],rgb[v1_id],
			   rgb[v2_id])
	    }
	  }
       }
       t_id = smTri_next_ccw_nbr(smMesh,t,s_id);
       t = SM_NTH_TRI(smMesh,t_id);
      }while(t!=tri);
    }
  }
  else
    if(!cull  || QT_IS_FLAG(qt))
      for(i=0; i < 4; i++)
	smRender_qtree(sm,QT_NTH_CHILD(qt,i),vp,wp,rgb,which,cull);
}


/* 
 * smRender_mesh(sm,view,cull) : Render mesh Triangles
 *   SM *sm;                   : mesh
 *   VIEW *view;               : current view
 *   int cull;                 : cull Flag 
 *
 *   If cull is TRUE, first mark tris in current 
 *   frustum and only render them. Normally, cull will be FALSE only if 
 *   it is known that all tris lie in frustum, e.g. after a rebuild
 *
 */
smRender_mesh(sm,view,cull)
SM *sm;
VIEW *view;
int cull;
{
  SFLOAT (*wp)[3];
  BYTE  (*rgb)[3];
  int i;
  STREE *st= SM_LOCATOR(sm);

  wp = SM_WP(sm);
  rgb =SM_RGB(sm);

  smClear_flags(sm,T_NEW_FLAG);

  if(cull)
    smCull(sm,view,SM_ALL_LEVELS);

  glPushAttrib(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_DEPTH_TEST);
  
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  /* move relative to the new view */
  glTranslated(view->vp[0],view->vp[1],view->vp[2]);

  /* The points are a distance of 1 away from the origin: if necessary 
     scale so that they fit in frustum and are therefore not clipped away
     */
  if(dev_zmin >= 0.99)
  {
    double d;

    d = (dev_zmin+dev_zmax)/2.0;
    glScaled(d,d,d);
  }
  /* move points to unit sphere at origin */
  glTranslated(-SM_VIEW_CENTER(sm)[0],-SM_VIEW_CENTER(sm)[1],
	       -SM_VIEW_CENTER(sm)[2]);

  glBegin(GL_TRIANGLES);
  for(i=0; i < ST_NUM_ROOT_NODES; i++)
      smRender_qtree(sm,ST_ROOT_QT(st,i),view->vp,wp,rgb,SM_RENDER_BG,cull);
  glEnd();

  glPopMatrix();

  glEnable(GL_DEPTH_TEST);

  glBegin(GL_TRIANGLES);
  for(i=0; i < ST_NUM_ROOT_NODES; i++)
    smRender_qtree(sm,ST_ROOT_QT(st,i),view->vp,wp,rgb,SM_RENDER_FG,cull);
  glEnd();

  glPopAttrib();
}

/* 
 * smRender_mesh_dl(sm,view) : Render stree utilizing display lists 
 * SM *sm;                   : mesh
 * VIEW *view;               : current view
 */
smRender_mesh_dl(sm,view)
SM *sm;
VIEW *view;
{
  SFLOAT (*wp)[3];
  BYTE  (*rgb)[3];
  STREE *st;
  int i;

  if(SM_DL_LEVELS == 0)
  {
    if(!Display_lists[0][0])
    {
      Display_lists[0][0] = 1;
      glNewList(Display_lists[0][0],GL_COMPILE_AND_EXECUTE);
      smRender_mesh(sm,view,FALSE);
      glEndList();
    }
    else
      glCallList(Display_lists[0][0]);

    return;
  }    
  smClear_flags(sm,T_NEW_FLAG);

  smCull(sm,view,SM_DL_LEVELS);

  st = SM_LOCATOR(sm);

  wp = SM_WP(sm);
  rgb =SM_RGB(sm);

    /* For all active quadtree nodes- first render bg tris, then fg */
    /* If display list exists, use otherwise create/display list */
   glPushAttrib(GL_DEPTH_BUFFER_BIT);
   glDisable(GL_DEPTH_TEST);

   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();

   /* move relative to the new view */
   glTranslated(view->vp[0],view->vp[1],view->vp[2]);

   /* The points are a distance of 1 away from the origin: if necessary 
      scale so that they fit in frustum and are therefore not clipped away
      */
   if(dev_zmin >= 0.99)
   {
     double d;
     d = (dev_zmin+dev_zmax)/2.0;
     glScaled(d,d,d);
   }
    /* move points to unit sphere at origin */
   glTranslated(-SM_VIEW_CENTER(sm)[0],-SM_VIEW_CENTER(sm)[1],
		   -SM_VIEW_CENTER(sm)[2]);
   for(i=0; i < ST_NUM_ROOT_NODES; i++)
     smRender_qtree_dl(sm,ST_ROOT_QT(st,i),view->vp,wp,rgb,i,0,1,
		       SM_DL_LEVELS,8,SM_RENDER_BG);
   glPopMatrix();

   glEnable(GL_DEPTH_TEST);
   for(i=0; i < ST_NUM_ROOT_NODES; i++)
     smRender_qtree_dl(sm,ST_ROOT_QT(st,i),view->vp,wp,rgb,i,0,1,
			SM_DL_LEVELS,8,SM_RENDER_FG);
   glPopAttrib();
}



/*
 * smRender_tris(sm,view,render_flag)  : Render all of the mesh triangles
 *  SM *sm             : current geometry
 *  VIEW *view         : current view
 *  int render_flag    : if render_flag & SM_RENDER_CULL: do culling first
 * 
 * Renders mesh by traversing triangle list and drawing all active tris-
 * background tris first, then foreground and mixed tris
 */
smRender_tris(sm,view,render_flag)
SM *sm;
VIEW *view;
int render_flag;
{
  int32  *active_flag,*bg_flag;
  SFLOAT (*wp)[3];
  BYTE  (*rgb)[3];

  wp = SM_WP(sm);
  rgb = SM_RGB(sm);
  active_flag = SM_NTH_FLAGS(sm,T_ACTIVE_FLAG);
  bg_flag = SM_NTH_FLAGS(sm,T_BG_FLAG);

  if(render_flag & SM_RENDER_CULL)
    smCull(sm,view,SM_ALL_LEVELS);
   
  /* Render triangles made up of points at infinity by turning off
     depth-buffering and projecting the points onto a sphere around the view*/
  glPushAttrib(GL_DEPTH_BUFFER_BIT);
  glDisable(GL_DEPTH_TEST);
  smRender_bg_tris(sm,view->vp,active_flag,bg_flag,wp,rgb);

  /* Render triangles containing world-space points */
  glEnable(GL_DEPTH_TEST);
  smRender_fg_tris(sm,view->vp,active_flag,bg_flag,wp,rgb);

  glPopAttrib();

}

/* Clear all of the display lists */
clear_display_lists()
{
  int i;
  for(i=0; i< SM_DL_LISTS; i++)
  {
    if(Display_lists[i][0])
    { /* Clear the foreground display list */
      glDeleteLists(Display_lists[i][0],1);
      Display_lists[i][0] = 0;
    }
    if(Display_lists[i][1])
    { /* Clear the background display list */
      glDeleteLists(Display_lists[i][1],1);
      Display_lists[i][1] = 0;
    }
  }
}

/*
 * qtClear_dl(qt,i,level_i,level,max_level,leaf_cnt) :clear display lists
 * QUADTREE *qt;               : Quadtree node   
 * int i;                      : index into list of display lists for this node
 * int level_i;                : index for first node at this level
 * int level,max_level;        : current level, maximum level to descend
 * int leaf_cnt;               : number of leaves at this level
 *
 *  For each node under this node that has its flag set: delete all 
 *  existing display lists. Display lists are stored in an array indexed as 
 *  if the quadtree was traversed in a breadth first order (indices 0-7 are
 *  the 8 quadtree roots, indices 8-11 the first level children of root 0, 
 *  indices 12-15 the children of root 1, etc). It is assumes that the display
 *  lists will only be stored for a small number of levels: if this is not
 *  true, a hashing scheme would work better for storing/retrieving the
 *  display lists
 */
qtClear_dl(qt,i,level_i,level,max_level,leaf_cnt)
QUADTREE qt;
int i,level_i,level,max_level,leaf_cnt;
{
  int j;

  if(QT_IS_EMPTY(qt))
    return;
  if(QT_IS_LEAF(qt) || level== max_level)
  {
    if(QT_IS_LEAF(qt))
    {
      if(!QT_LEAF_IS_FLAG(qt))
	return;
    }
    else
      if(!QT_IS_FLAG(qt))
	return;
    if(Display_lists[i][0])
    {
      glDeleteLists(Display_lists[i][0],1);
      Display_lists[i][0] = 0;
    }
    if(Display_lists[i][1])
    {
      glDeleteLists(Display_lists[i][1],1);
      Display_lists[i][1] = 0;
    }
  }
  else
    if(QT_IS_FLAG(qt))
    {
      /* Calculate the index for the first child given the values
	 of the parent at  the current level
       */
      i = ((i - level_i)<< 2) + level_i + leaf_cnt;
      level_i += leaf_cnt;
      leaf_cnt <<= 2;
      for(j=0; j < 4; j++)
	qtClear_dl(QT_NTH_CHILD(qt,j),i+j,level_i,level+1,max_level,
			    leaf_cnt);
    }
}

/* 
 * smInvalidate_view(sm,view) : Invalidate rendering representation for view
 * SM *sm;                    : mesh
 * VIEW *view;                : current view
 *
 * Delete the existing display lists for geometry in the current
 * view frustum: Called when the geometry in the frustum has been changed
 */
smInvalidate_view(sm,view)
SM *sm;
VIEW *view;
{
  int i;

  if(SM_DL_LEVELS == 0)
  {
    if(Display_lists[0][0])
    {
      glDeleteLists(Display_lists[0][0],1);
      Display_lists[0][0] = 0;
    }
    return;
  }	
  /* Mark qtree nodes/tris in frustum */
  smCull(sm,view,SM_DL_LEVELS);

  /* Invalidate display_lists in marked qtree nodes */
   for(i=0; i < ST_NUM_ROOT_NODES; i++)
     qtClear_dl(ST_ROOT_QT(SM_LOCATOR(sm),i),i,0,1,SM_DL_LEVELS,8);

}


/*
 * smRender(sm,view, qual): render OpenGL output geometry 
 * SM *sm;                : current mesh representation
 * VIEW	*view;		  : desired view
 * int	qual;             : quality level (percentage on linear time scale)
 *
 * Render the current mesh: 
 * recompute tone mapping if full redraw and specified: 
 * if moving (i.e. qual < MAXQUALITY)
 *     render the cached display lists, if quality drops
 *     below threshold, render approximation instead
 *  if stationary
 *     render mesh geometry without display lists, unless up-to-date
 *     display lists already exist.
 */
smRender(sm,view,qual)
SM *sm;
VIEW *view;
int qual;
{
  /* Unless quality > MAXQUALITY, render using display lists */
  if(qual <= MAXQUALITY)
  {
    /* If quality above threshold: render mesh*/
    if(qual > (MAXQUALITY*2/4) )
      /* render stree using display lists */
      smRender_mesh_dl(sm,view);
    else
    {
      /* If quality below threshold, use approximate rendering */
      smRender_approx(sm,qual,view);
    }
  }
  else
  {
    /* render stree without display lists */
      smRender_mesh(sm,view,TRUE);
  }
}


/*
 * smUpdate(view, qual)	: update OpenGL output geometry 
 * VIEW	*view;		: desired view
 * int	qual;           : quality level (percentage on linear time scale)
 *
 * Draw new geometric representation using OpenGL calls.  Assume that the
 * view has already been set up and the correct frame buffer has been
 * selected for drawing.  The quality level is on a linear scale, where 100%
 * is full (final) quality.  It is not necessary to redraw geometry that has
 * been output since the last call to smClean().(The last view drawn will
 * be view==&odev.v each time.)
 */
smUpdate(view,qual)
   VIEW *view;
   int qual;
{
  
  /* Is there anything to render? */
  if(!smMesh || SM_NUM_TRI(smMesh)<=0)
    return;
  
  /* Is viewer MOVING?*/
  if(qual < MAXQUALITY)
  {
    if(smIncremental)
      smUpdate_tm(smMesh);

    /* Render mesh using display lists */
    smRender(smMesh,view,qual);
    return;
  }
  /* Viewer is STATIONARY */
  /* Has view moved epsilon from canonical view? (epsilon= percentage
     (SM_VIEW_FRAC) of running average of the distance of the sample points 
     from the canonical view */
  if(DIST(view->vp,SM_VIEW_CENTER(smMesh)) > SM_ALLOWED_VIEW_CHANGE(smMesh))
  {
    /* Must rebuild mesh with current view as new canonical view  */
    smRebuild(smMesh,view);
    /* Existing display lists and tonemapping are no longer valid */
    clear_display_lists();
    smCompute_mapping = FALSE;
    smUpdate_tm(smMesh);
    /* Render all the triangles in the new mesh */
    smRender(smMesh,view,qual+1);
  }
  else 
    /* Has a complete redraw been requested ?*/
    if(smClean_notify)
    {
      if(smIncremental)
	smUpdate_tm(smMesh);
      smIncremental = FALSE;
      smRender(smMesh,view,qual);
    }
    else
    {
      smUpdate_tm(smMesh);
      /* If number of new triangles relatively small: do incremental update */
      /* Mark Existing display lists in frustum invalid */
      if(!smIncremental)
      {
	smInvalidate_view(smMesh,view);
	smIncremental = TRUE;
      }
      smRender_inc(smMesh,view->vp);
  }
  /* This is our final update iff qual==MAXQUALITY and view==&odev.v */
  if( (qual >= MAXQUALITY) && (view == &(odev.v)))
  {
    /* reset rendering flags */
    smClean_notify = FALSE;
    if(smIncremental)
      smClear_flags(smMesh,T_NEW_FLAG);
    qtCache_init(0);
  }

}







