#ifndef lint
static const char	RCSid[] = "$Id: sm_stree.c,v 3.13 2003/02/22 02:07:25 greg Exp $";
#endif
/*
 * sm_stree.c
 *  An stree (spherical quadtree) is defined by an octahedron in 
 *  canonical form,and a world center point. Each face of the 
 *  octahedron is adaptively subdivided as a planar triangular quadtree.
 *  World space geometry is projected onto the quadtree faces from the
 *  sphere center.
 */
#include "standard.h"
#include "sm_list.h"
#include "sm_flag.h"
#include "sm_geom.h"
#include "object.h"
#include "sm_qtree.h"
#include "sm_stree.h"


#ifdef TEST_DRIVER
extern FVECT Pick_point[500],Pick_v0[500],Pick_v1[500],Pick_v2[500];
extern int Pick_cnt;
#endif
/* octahedron coordinates */
FVECT stDefault_base[6] = {  {1.,0.,0.},{0.,1.,0.}, {0.,0.,1.},
			    {-1.,0.,0.},{0.,-1.,0.},{0.,0.,-1.}};
/* octahedron triangle vertices */
int stBase_verts[8][3] = { {0,1,2},{3,1,2},{0,4,2},{3,4,2},
			   {0,1,5},{3,1,5},{0,4,5},{3,4,5}};
/* octahedron triangle nbrs ; nbr i is the face opposite vertex i*/
int stBase_nbrs[8][3] =  { {1,2,4},{0,3,5},{3,0,6},{2,1,7},
			   {5,6,0},{4,7,1},{7,4,2},{6,5,3}};
int stRoot_indices[8][3] = {{1,1,1},{-1,1,1},{1,-1,1},{-1,-1,1},
			    {1,1,-1},{-1,1,-1},{1,-1,-1},{-1,-1,-1}};
/*
 +z   y                -z   y
      |                     |
 1    |   0             5   |   4
______|______ x      _______|______ x
 3    |   2             7   |   6   
      |                     |

Nbrs
  +z   y                -z   y
     /0|1\                 /1|0\
 5  /  |  \ 4             /  |  \
   /(1)|(0)\           1 /(5)|(4)\ 0
  /    |    \           /    |    \
 /2   1|0   2\         /2   0|1   2\
/------|------\x      /------|------\x
\0    1|2    0/       \0    2|2    1/ 
 \     |     /         \     |     /
 7\ (3)|(2) / 6       3 \ (7)|(6) / 2
   \   |   /             \   |   /
    \ 2|1 /               \ 1|0 /
*/


stInit(st)
STREE *st;
{
  int i,j;

  qtDone();

  ST_TOP_QT(st) = qtAlloc();
  ST_BOTTOM_QT(st) = qtAlloc();
  /* Clear the children */

   QT_CLEAR_CHILDREN(ST_TOP_QT(st));
   QT_CLEAR_CHILDREN(ST_BOTTOM_QT(st));
}

stFree(st)
STREE *st;
{
  qtDone();
  free(st);
}

/* Allocates a stree structure  and creates octahedron base */
STREE
*stAlloc(st)
STREE *st;
{
  int i,m;
  FVECT v0,v1,v2;
  FVECT n;
  
  if(!st)
    if(!(st = (STREE *)malloc(sizeof(STREE))))
       error(SYSTEM,"stAlloc(): Unable to allocate memory\n");

  /* Allocate the top and bottom quadtree root nodes */
  stInit(st);
  
  return(st);
}

#define BARY_INT(v,b)  if((v)>2.0) (b) = MAXBCOORD;else \
  if((v)<-2.0) (b)=-MAXBCOORD;else (b)=(BCOORD)((v)*MAXBCOORD2);

vert_to_qt_frame(root,v,b)
int root;
FVECT v;
BCOORD b[3];
{
  int i;
  double scale;
  double d0,d1,d2;

  if(STR_NTH_INDEX(root,0)==-1)
    d0 = -v[0];
  else
    d0 = v[0];
  if(STR_NTH_INDEX(root,1)==-1)
    d1 = -v[1];
  else
    d1 = v[1];
  if(STR_NTH_INDEX(root,2)==-1)
    d2 = -v[2];
  else
    d2 = v[2];

  /* Plane is now x+y+z = 1 - intersection of pt ray is qtv/den */
  scale = 1.0/ (d0 + d1 + d2);
  d0 *= scale;
  d1 *= scale;
  d2 *= scale;

  BARY_INT(d0,b[0])
  BARY_INT(d1,b[1])
  BARY_INT(d2,b[2])
}




ray_to_qt_frame(root,v,dir,b,db)
int root;
FVECT v,dir;
BCOORD b[3],db[3];
{
  int i;
  double scale;
  double d0,d1,d2;
  double dir0,dir1,dir2;

  if(STR_NTH_INDEX(root,0)==-1)
  {
    d0 = -v[0];
    dir0 = -dir[0];
  }
  else
  {
    d0 = v[0];
    dir0 = dir[0];
  }
  if(STR_NTH_INDEX(root,1)==-1)
  {
    d1 = -v[1];
    dir1 = -dir[1];
  }
  else
  {
    d1 = v[1];
    dir1 = dir[1];
  }
  if(STR_NTH_INDEX(root,2)==-1)
  {
    d2 = -v[2];
    dir2 = -dir[2];
  }
  else
  {
    d2 = v[2];
    dir2 = dir[2];
  }
  /* Plane is now x+y+z = 1 - intersection of pt ray is qtv/den */
  scale = 1.0/ (d0 + d1 + d2);
  d0 *= scale;
  d1 *= scale;
  d2 *= scale;

  /* Calculate intersection point of orig+dir: This calculation is done
     after the origin is projected into the plane in order to constrain
     the projection( i.e. the size of the projection of the unit direction
     vector translated to the origin depends on how close
     the origin is to the view center
     */
  /* Must divide by at least root2 to insure that projection will fit 
     int [-2,2] bounds: assumed length is 1: therefore greatest projection
     from endpoint of triangle is at 45 degrees or projected length of root2
  */
  dir0 = d0 + dir0*0.5;
  dir1 = d1 + dir1*0.5;
  dir2 = d2 + dir2*0.5;

  scale = 1.0/ (dir0 + dir1 + dir2);
  dir0 *= scale;
  dir1 *= scale;
  dir2 *= scale;

  BARY_INT(d0,b[0])
  BARY_INT(d1,b[1])
  BARY_INT(d2,b[2])
  BARY_INT(dir0,db[0])
  BARY_INT(dir1,db[1])
  BARY_INT(dir2,db[2])

  db[0]  -= b[0];
  db[1]  -= b[1];
  db[2]  -= b[2];
}

qt_frame_to_vert(root,b,v)
int root;
BCOORD b[3];
FVECT v;
{
  int i;
  double d0,d1,d2;

  d0 = b[0]/(double)MAXBCOORD2;
  d1 = b[1]/(double)MAXBCOORD2;
  d2 = b[2]/(double)MAXBCOORD2;
  
  if(STR_NTH_INDEX(root,0)==-1)
    v[0] = -d0;
  else
    v[0] = d0;
  if(STR_NTH_INDEX(root,1)==-1)
    v[1] = -d1;
  else
    v[1] = d1;
  if(STR_NTH_INDEX(root,2)==-1)
    v[2] = -d2;
  else
    v[2] = d2;
}


/* Return quadtree leaf node containing point 'p'*/
QUADTREE
stPoint_locate(st,p)
    STREE *st;
    FVECT p;
{
    QUADTREE qt;
    BCOORD bcoordi[3];
    int i;

    /* Find root quadtree that contains p */
    i = stLocate_root(p);
    qt = ST_ROOT_QT(st,i);
    
     /* Will return lowest level triangle containing point: It the
       point is on an edge or vertex: will return first associated
       triangle encountered in the child traversal- the others can
       be derived using triangle adjacency information
    */
    if(QT_IS_TREE(qt))
    {  
      vert_to_qt_frame(i,p,bcoordi);
      i = bary_child(bcoordi);
      return(qtLocate(QT_NTH_CHILD(qt,i),bcoordi));
    }
    else
      return(qt);
}

static unsigned int nbr_b[8][3] ={{2,4,16},{1,8,32},{8,1,64},{4,2,128},
			   {32,64,1},{16,128,2},{128,16,4},{64,32,8}};
unsigned int
stTri_cells(st,v)
     STREE *st;
     FVECT v[3];
{
  unsigned int cells,cross;
  unsigned int vcell[3];
  double t0,t1;
  int i,inext;

  /* First find base cells that tri vertices are in (0-7)*/
  vcell[0] = stLocate_root(v[0]);
  vcell[1] = stLocate_root(v[1]);
  vcell[2] = stLocate_root(v[2]);

  /* If all in same cell- return that bit only */
  if(vcell[0] == vcell[1] && vcell[1] == vcell[2])
    return( 1 << vcell[0]);

  cells = 0;
  for(i=0;i<3; i++)
  {
    if(i==2)
      inext = 0;
    else
      inext = i+1;
    /* Mark cell containing initial vertex */
    cells |= 1 << vcell[i];

    /* Take the exclusive or: will have bits set where edge crosses axis=0*/
    cross = vcell[i] ^ vcell[inext];
    /* If crosses 2 planes: then have 2 options for edge crossing-pick closest
     otherwise just hits two*/
    /* Neighbors are zyx */
    switch(cross){
    case 3: /* crosses x=0 and y=0 */
      t0 = -v[i][0]/(v[inext][0]-v[i][0]);
      t1 = -v[i][1]/(v[inext][1]-v[i][1]);
      if(t0==t1)
	break;
      else if(t0 < t1)
	cells |= nbr_b[vcell[i]][0];
	  else
	    cells |= nbr_b[vcell[i]][1];
      break;
    case 5: /* crosses x=0 and z=0 */
      t0 = -v[i][0]/(v[inext][0]-v[i][0]);
      t1 = -v[i][2]/(v[inext][2]-v[i][2]);
      if(t0==t1)
	break;
      else if(t0 < t1)
	cells |= nbr_b[vcell[i]][0];
	  else
	    cells |=nbr_b[vcell[i]][2];

      break;
    case 6:/* crosses  z=0 and y=0 */
      t0 = -v[i][2]/(v[inext][2]-v[i][2]);
      t1 = -v[i][1]/(v[inext][1]-v[i][1]);
      if(t0==t1)
	break;
      else if(t0 < t1)
      {
	cells |= nbr_b[vcell[i]][2];
      }
      else
      {
	cells |=nbr_b[vcell[i]][1];
      }
      break;
    case 7:
      error(CONSISTENCY," Insert:Edge shouldnt be able to span 3 cells");
      break;
    }
  }
  return(cells);
}


stRoot_trace_ray(qt,root,orig,dir,nextptr,func,f)
   QUADTREE qt;
   int root;
   FVECT orig,dir;
   int *nextptr;
   FUNC func;
   int *f;
{
  double br[3];
  BCOORD bi[3],dbi[3];
  
  /* Project the origin onto the root node plane */
  /* Find the intersection point of the origin */
  ray_to_qt_frame(root,orig,dir,bi,dbi);

  /* trace the ray starting with this node */
  qtTrace_ray(qt,bi,dbi[0],dbi[1],dbi[2],nextptr,0,0,func,f);
  if(!QT_FLAG_IS_DONE(*f))
    qt_frame_to_vert(root,bi,orig);

}

/* Trace ray 'orig-dir' through stree and apply 'func(qtptr,f,argptr)' at each
   node that it intersects
*/
int
stTrace_ray(st,orig,dir,func)
   STREE *st;
   FVECT orig,dir;
   FUNC func;
{
    int next,last,i,f=0;
    QUADTREE qt;
    FVECT o,n,v;
    double pd,t,d;

    VCOPY(o,orig);
#ifdef TEST_DRIVER
       Pick_cnt=0;
#endif;
    /* Find the qt node that o falls in */
    i = stLocate_root(o);
    qt = ST_ROOT_QT(st,i);
    
    stRoot_trace_ray(qt,i,o,dir,&next,func,&f);

    if(QT_FLAG_IS_DONE(f))
      return(TRUE);
    /*
    d = DOT(orig,dir)/sqrt(DOT(orig,orig));
    VSUM(v,orig,dir,-d);
    */
    /* Crossed over to next cell: id = nbr */
    while(1) 
      {
	/* test if ray crosses plane between this quadtree triangle and
	   its neighbor- if it does then find intersection point with 
	   ray and plane- this is the new origin
	   */
	if(next == INVALID)
	  return(FALSE);
	/*
	if(DOT(o,v) < 0.0)
	  return(FALSE);
	  */
	i = stBase_nbrs[i][next];
	qt = ST_ROOT_QT(st,i);
	stRoot_trace_ray(qt,i,o,dir,&next,func,&f);
	if(QT_FLAG_IS_DONE(f))
	  return(TRUE);
      }
}


stVisit_poly(st,verts,l,root,func,n)
STREE *st;
FVECT *verts;
LIST *l;
unsigned int root;
FUNC func;
int n;
{
  int id0,id1,id2;
  FVECT tri[3];

  id0 = pop_list(&l);
  id1 = pop_list(&l);
  while(l)
  {
    id2 = pop_list(&l);
    VCOPY(tri[0],verts[id0]);
    VCOPY(tri[1],verts[id1]);
    VCOPY(tri[2],verts[id2]);
    stRoot_visit_tri(st,root,tri,func,n);
    id1 = id2;
  }
}
/* Assumption: know crosses plane:dont need to check for 'on' case */
intersect_edge_coord_plane(v0,v1,w,r)
FVECT v0,v1;
int w;
FVECT r;
{
  FVECT dv;
  int wnext;
  double t;

  VSUB(dv,v1,v0);
  t = -v0[w]/dv[w];
  r[w] = 0.0;
  wnext = (w+1)%3;
  r[wnext] = v0[wnext] + dv[wnext]*t;
  wnext = (w+2)%3;
  r[wnext] = v0[wnext] + dv[wnext]*t;
}


stVisit_clip(st,i,verts,vcnt,l,cell,func,n)
     STREE *st;
     int i;
     FVECT *verts;
     int *vcnt;
     LIST *l;
     unsigned int cell;
     FUNC func;
     int n;
{

  LIST *labove,*lbelow,*endb,*enda;
  int last = -1;
  int id,last_id;
  int first,first_id;
  unsigned int cellb;

  labove = lbelow = NULL;
  enda = endb = NULL;
  while(l)
  {
    id = pop_list(&l);
    if(ZERO(verts[id][i]))
    {
      if(last==-1)
      {/* add below and above */
	first = 2;
	first_id= id;
      }
      lbelow=add_data(lbelow,id,&endb);
      labove=add_data(labove,id,&enda);
      last_id = id;
      last = 2;
      continue;
    }
    if(verts[id][i] < 0)
    {
      if(last != 1)
      {
	lbelow=add_data(lbelow,id,&endb);
	if(last==-1)
	{
	  first = 0;
	  first_id = id;
	}
	last_id = id;
	last = 0;
	continue;
      }
      /* intersect_edges */
      intersect_edge_coord_plane(verts[last_id],verts[id],i,verts[*vcnt]);
      /*newpoint goes to above and below*/
      lbelow=add_data(lbelow,*vcnt,&endb);
      lbelow=add_data(lbelow,id,&endb);
      labove=add_data(labove,*vcnt,&enda);
      last = 0;
      last_id = id;
      (*vcnt)++;
    }
    else
    {
      if(last != 0)
      {
	labove=add_data(labove,id,&enda);
	if(last==-1)
	{
	  first = 1;
	  first_id = id;
	}
	last_id = id;
	last = 1;
	continue;
      }
      /* intersect_edges */
      /*newpoint goes to above and below*/
      intersect_edge_coord_plane(verts[last_id],verts[id],i,verts[*vcnt]);
      lbelow=add_data(lbelow,*vcnt,&endb);
      labove=add_data(labove,*vcnt,&enda);
      labove=add_data(labove,id,&enda);
      last_id = id;
      (*vcnt)++;
      last = 1;
    }
  }
  if(first != 2 && first != last)
  {
    intersect_edge_coord_plane(verts[id],verts[first_id],i,verts[*vcnt]);
    /*newpoint goes to above and below*/
    lbelow=add_data(lbelow,*vcnt,&endb);
    labove=add_data(labove,*vcnt,&enda);
    (*vcnt)++;

  }
  if(i==2)
  {
    if(lbelow)
    {
      if(LIST_NEXT(lbelow) && LIST_NEXT(LIST_NEXT(lbelow)))
      {
	cellb = cell | (1 << i);
	stVisit_poly(st,verts,lbelow,cellb,func,n);
      }
      else
	free_list(lbelow);
    }
    if(labove)
     {
      if(LIST_NEXT(labove) && LIST_NEXT(LIST_NEXT(labove)))
	stVisit_poly(st,verts,labove,cell,func,n);
      else
	free_list(labove);
     }
  }
  else
  {
    if(lbelow)
    {
      if(LIST_NEXT(lbelow) && LIST_NEXT(LIST_NEXT(lbelow)))
	{
	  cellb = cell | (1 << i);
	  stVisit_clip(st,i+1,verts,vcnt,lbelow,cellb,func,n);
	}
      else
	free_list(lbelow);
    }
    if(labove)
     {
       if(LIST_NEXT(labove) && LIST_NEXT(LIST_NEXT(labove)))
	 stVisit_clip(st,i+1,verts,vcnt,labove,cell,func,n);
       else
	 free_list(labove);
     }
  }

}

stVisit(st,tri,func,n)
   STREE *st;
   FVECT tri[3];
   FUNC func;
   int n;
{
    int r0,r1,r2;
    LIST *l;

    r0 = stLocate_root(tri[0]);
    r1 = stLocate_root(tri[1]);
    r2 = stLocate_root(tri[2]);
    if(r0 == r1 && r1==r2)
      stRoot_visit_tri(st,r0,tri,func,n);
    else
      {
	FVECT verts[ST_CLIP_VERTS];
	int cnt;

	VCOPY(verts[0],tri[0]);
	VCOPY(verts[1],tri[1]);
	VCOPY(verts[2],tri[2]); 
	
	l = add_data(NULL,0,NULL);
	l = add_data(l,1,NULL);
	l = add_data(l,2,NULL);
	cnt = 3;
	stVisit_clip(st,0,verts,&cnt,l,0,func,n);
      }
}


BCOORD qtRoot[3][3] = { {MAXBCOORD2,0,0},{0,MAXBCOORD2,0},{0,0,MAXBCOORD2}};


convert_tri_to_frame(root,tri,b0,b1,b2,db10,db21,db02)
int root;
FVECT tri[3];
BCOORD b0[3],b1[3],b2[3];
BCOORD db10[3],db21[3],db02[3];
{
  /* Project the vertex into the qtree plane */
  vert_to_qt_frame(root,tri[0],b0);
  vert_to_qt_frame(root,tri[1],b1);
  vert_to_qt_frame(root,tri[2],b2);

  /* calculate triangle edge differences in new frame */
  db10[0] = b1[0] - b0[0]; db10[1] = b1[1] - b0[1]; db10[2] = b1[2] - b0[2];
  db21[0] = b2[0] - b1[0]; db21[1] = b2[1] - b1[1]; db21[2] = b2[2] - b1[2];
  db02[0] = b0[0] - b2[0]; db02[1] = b0[1] - b2[1]; db02[2] = b0[2] - b2[2];
}


QUADTREE
stRoot_insert_tri(st,root,tri,f)
   STREE *st;
   int root;
   FVECT tri[3];
   FUNC f;
{
  BCOORD b0[3],b1[3],b2[3];
  BCOORD db10[3],db21[3],db02[3];
  unsigned int s0,s1,s2,sq0,sq1,sq2;
  QUADTREE qt;

  /* Map the triangle vertices into the canonical barycentric frame */
  convert_tri_to_frame(root,tri,b0,b1,b2,db10,db21,db02);

  /* Calculate initial sidedness info */
  SIDES_GTR(b0,b1,b2,s0,s1,s2,qtRoot[1][0],qtRoot[0][1],qtRoot[0][2]);
  SIDES_GTR(b0,b1,b2,sq0,sq1,sq2,qtRoot[0][0],qtRoot[1][1],qtRoot[2][2]);

  qt = ST_ROOT_QT(st,root);
  /* Visit cells that triangle intersects */
  qt = qtInsert_tri(root,qt,qtRoot[0],qtRoot[1],qtRoot[2],
       b0,b1,b2,db10,db21,db02,MAXBCOORD2 >> 1,s0,s1,s2, sq0,sq1,sq2,f,0);

  return(qt);
}

stRoot_visit_tri(st,root,tri,f,n)
   STREE *st;
   int root;
   FVECT tri[3];
   FUNC f;
   int n;
{
  BCOORD b0[3],b1[3],b2[3];
  BCOORD db10[3],db21[3],db02[3];
  unsigned int s0,s1,s2,sq0,sq1,sq2;
  QUADTREE qt;

  /* Map the triangle vertices into the canonical barycentric frame */
  convert_tri_to_frame(root,tri,b0,b1,b2,db10,db21,db02);

  /* Calculate initial sidedness info */
  SIDES_GTR(b0,b1,b2,s0,s1,s2,qtRoot[1][0],qtRoot[0][1],qtRoot[0][2]);
  SIDES_GTR(b0,b1,b2,sq0,sq1,sq2,qtRoot[0][0],qtRoot[1][1],qtRoot[2][2]);

  qt = ST_ROOT_QT(st,root);
  QT_SET_FLAG(ST_QT(st,root));
  /* Visit cells that triangle intersects */
  qtVisit_tri(root,qt,qtRoot[0],qtRoot[1],qtRoot[2],
       b0,b1,b2,db10,db21,db02,MAXBCOORD2 >> 1,s0,s1,s2, sq0,sq1,sq2,f,n);

}

stInsert_tri(st,tri,f)
   STREE *st;
   FVECT tri[3];
   FUNC f;
{
  unsigned int cells,which;
  int root;
  

  /* calculate entry/exit points of edges through the cells */
  cells = stTri_cells(st,tri);

  /* For each cell that quadtree intersects: Map the triangle vertices into
     the canonical barycentric frame of (1,0,0), (0,1,0),(0,0,1). Insert 
     by first doing a trivial reject on the interior nodes, and then a 
     tri/tri intersection at the leaf nodes.
  */
  for(root=0,which=1; root < ST_NUM_ROOT_NODES; root++,which <<= 1)
  {
    /* For each of the quadtree roots: check if marked as intersecting tri*/
    if(cells & which)
      /* Visit tri cells */
      ST_ROOT_QT(st,root) = stRoot_insert_tri(st,root,tri,f);
  }
}

stInsert_samp(st,p,f)
   STREE *st;
   FVECT p;
   FUNC f;
{

    QUADTREE qt;
    BCOORD bcoordi[3];
    int i,done;

    /* Find root quadtree that contains p */
    i = stLocate_root(p);
    qt = ST_ROOT_QT(st,i);
    
    vert_to_qt_frame(i,p,bcoordi);
    ST_ROOT_QT(st,i) =  qtInsert_point(i,qt,EMPTY,qtRoot[0],qtRoot[1],
			  qtRoot[2],bcoordi,MAXBCOORD2>>1,f,0,&done);

}










