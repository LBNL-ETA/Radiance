#ifndef lint
static const char	RCSid[] = "$Id: sm_del.c,v 3.14 2003/04/23 00:52:34 greg Exp $";
#endif
/*
 *  sm_del.c
 */
#include "standard.h"
#include "sm_flag.h"
#include "sm_list.h"
#include "sm_geom.h"
#include "sm.h"

#define MAX_EDGE_INIT 100
static int Max_edges= MAX_EDGE_INIT;
static EDGE *Edges=NULL;
static int Ecnt=0;

smFree_tri(sm,id,t)
SM *sm;
int id;
TRI *t;
{
  /* Add to the free_list */
  T_NEXT_FREE(t) = SM_FREE_TRIS(sm);
  SM_FREE_TRIS(sm) = id;
#ifdef DEBUG
  T_VALID_FLAG(t) = INVALID;
#endif
}

/* Assumes mesh pointers have been cleaned up appropriately: just deletes from
   Point location and triangle data structure
*/
smDelete_tri(sm,t_id,t)
SM *sm;
int t_id;
TRI *t;
{

  /* NOTE: Assumes that a new triangle adjacent to each vertex
     has been added- before the deletion: replacing
     the old tri- and therefore dont need to dereference any pointers
     to id because the vertices can no longer
     point to tri id as being the first triangle pointer
  */
  SM_CLR_NTH_T_ACTIVE(sm,t_id);
  smFree_tri(sm,t_id,t);

}

int
eNew_edge()
{
  if(!Edges)
    if(!(Edges = (EDGE *)realloc(NULL,(Max_edges+1)*sizeof(EDGE))))
      goto memerr;

  if(Ecnt >= Max_edges)
    {
#ifdef DEBUG
      if(Max_edges > 10000)
	error(CONSISTENCY,"Too many edges in vertex loop\n");
#endif
      Max_edges += 100;
      if(!(Edges = (EDGE *)realloc((void *)Edges,(Max_edges+1)*sizeof(EDGE))))
	goto memerr;
    }
#ifdef DEBUG
  SET_E_NTH_TRI(Ecnt+1,0,INVALID);
  SET_E_NTH_TRI(Ecnt+1,1,INVALID);
#endif
  return(++Ecnt);

memerr:
  error(SYSTEM,"eNew_edge): Unable to allocate memory");
}

/* Return list of edges defining polygon formed by boundary of triangles 
adjacent to id. Return set of triangles adjacent to id to delete in delptr
*/
LIST 
*smVertexStar(sm,id)
SM *sm;
S_ID id;
{
    TRI *tri,*t_next;
    LIST *elist,*end;
    int e,t_id,t_next_id,b_id,v_id,t_last_id;
    S_ID v_next; 

    elist = end =  NULL;

    /* Get the first triangle adjacent to vertex id */
    t_id = SM_NTH_VERT(sm,id);
    tri = SM_NTH_TRI(sm,t_id);

    e = eNew_edge();
    /* Get the  next vertex on the polygon boundary */
    v_id = T_WHICH_V(tri,id);
    b_id = (v_id + 1)%3;
    /* Create an edge */
    SET_E_NTH_VERT(e,0,T_NTH_V(tri,b_id));
    SET_E_NTH_TRI(e,0,INVALID);
    SET_E_NTH_TRI(e,1,T_NTH_NBR(tri,v_id));
     v_next = T_NTH_V(tri,(b_id+1)%3);
    SET_E_NTH_VERT(e,1,v_next);
  
    elist = add_data_to_circular_list(elist,&end,e);
    t_next_id = t_id;
    t_next = tri;
    t_last_id = -1;
 
    /* Create a set to hold all of the triangles for deletion later */

    while((t_next_id = T_NTH_NBR(t_next,b_id)) != t_id)
    {
      e = eNew_edge();
      if(t_last_id != -1)
      {
	smDelete_tri(sm,t_last_id,t_next);
      }
      t_next = SM_NTH_TRI(sm,t_next_id);
      t_last_id = t_next_id;
      SET_E_NTH_VERT(e,0,v_next);
      SET_E_NTH_TRI(e,0,INVALID);
      v_id = T_WHICH_V(t_next,id);
      b_id = (v_id + 1)%3;
      SET_E_NTH_TRI(e,1,T_NTH_NBR(t_next,v_id));
      v_next = T_NTH_V(t_next,(b_id+1)%3);
      SET_E_NTH_VERT(e,1,v_next);
      elist = add_data_to_circular_list(elist,&end,e);

    }
    smDelete_tri(sm,t_last_id,t_next);
    smDelete_tri(sm,t_id,tri); 
   return(elist);
}

int
smTriangulate_add_tri(sm,id0,id1,id2,e0,e1,e2ptr)
SM *sm;
S_ID id0,id1,id2;
int e0,e1,*e2ptr;
{
  int t_id,e2;
  TRI *t;

  t_id = smAdd_tri(sm,id0,id1,id2,&t);
  if(*e2ptr == 0)
  {
    e2 = eNew_edge();
    SET_E_NTH_VERT(e2,0,id2);
    SET_E_NTH_VERT(e2,1,id0);
  }
  else
    e2 = *e2ptr;
  /* set appropriate tri for each edge*/ 
  SET_E_NTH_TRI(e0,0,t_id);
  SET_E_NTH_TRI(e1,0,t_id);
  SET_E_NTH_TRI(e2,0,t_id);
#ifdef DEBUG
#if DEBUG > 1 
  if(E_NTH_TRI(e0,1) == E_NTH_TRI(e0,0) ||
     E_NTH_TRI(e1,1) == E_NTH_TRI(e1,0) ||
     E_NTH_TRI(e2,1) == E_NTH_TRI(e2,0))
    error(CONSISTENCY,"Non manifold\n");
#endif
#endif
  *e2ptr = e2;
  return(t_id);

}

int 
smTriangulate_quad(sm,l)
     SM *sm;
     LIST *l;
{
  int e1,e2,e3,e4,e,t_id;
  S_ID id0,id1,id2,id3;
  FVECT v0,v1,v2,v3,n;
  double dp;
  TRI *tri;
  LIST *lptr;

#ifdef DEBUG
  if(LIST_NEXT(LIST_NEXT(LIST_NEXT(LIST_NEXT(l)))) != l)
  {
    eputs("smTriangulate_quad(): not a quadrilateral\n");
    return(FALSE);
  }
    eputs("smTriangulate_quad():\n");
#endif
  lptr=l;
  e1 = (int)LIST_DATA(l);
  id0 = E_NTH_VERT(e1,0);
  id1 = E_NTH_VERT(e1,1);
  VSUB(v0,SM_NTH_WV(sm,id0),SM_VIEW_CENTER(sm));
  VSUB(v1,SM_NTH_WV(sm,id1),SM_VIEW_CENTER(sm));   
  /* Get v2 */
  l = LIST_NEXT(l);
  e2 = (int)LIST_DATA(l);
  id2 = E_NTH_VERT(e2,1);
  VSUB(v2,SM_NTH_WV(sm,id2),SM_VIEW_CENTER(sm));
  /* Get v3 */
  l = LIST_NEXT(l);
  e3 = (int)LIST_DATA(l);
  id3 = E_NTH_VERT(e3,1);
  VSUB(v3,SM_NTH_WV(sm,id3),SM_VIEW_CENTER(sm));
  l = LIST_NEXT(l);
  e4 = (int)LIST_DATA(l);
  free_list(lptr);

  VCROSS(n,v0,v2);
  normalize(n);
  normalize(v1);
  dp = DOT(n,v1);
  e = 0;
  if(dp > 0) 
  {
    if(dp  >= EV_EPS)
    {
      normalize(v3);
      if(DOT(n,v3) < 0)	
      {
	t_id = smTriangulate_add_tri(sm,id0,id1,id2,e1,e2,&e);
	e *= -1;
	t_id = smTriangulate_add_tri(sm,id2,id3,id0,e3,e4,&e);
	return(TRUE);
      }
    }

  }
#ifdef DEBUG
#if DEBUG > 1
  VCROSS(n,v1,v3);
  normalize(n);
  normalize(v0);
  normalize(v2);
  dp = DOT(n,v2);
  if((dp < 0) || (dp < EV_EPS) || (DOT(n,v0) >= 0.0))	
    error(CONSISTENCY,"smTriangulate_quad: cannot triangulate\n");
#endif
#endif
  t_id = smTriangulate_add_tri(sm,id1,id2,id3,e2,e3,&e);
  e *= -1;
  t_id = smTriangulate_add_tri(sm,id3,id0,id1,e4,e1,&e);
  return(TRUE);
}

eIn_tri(e,t)
int e;
TRI *t;
{

  if(T_NTH_V(t,0)==E_NTH_VERT(e,0)) 
    return(T_NTH_V(t,1)==E_NTH_VERT(e,1)||T_NTH_V(t,2)==E_NTH_VERT(e,1));
  else
    if(T_NTH_V(t,1)==E_NTH_VERT(e,0))
      return(T_NTH_V(t,0)==E_NTH_VERT(e,1)||T_NTH_V(t,2)==E_NTH_VERT(e,1));
    else if(T_NTH_V(t,2)==E_NTH_VERT(e,0))
      return(T_NTH_V(t,0)==E_NTH_VERT(e,1)||T_NTH_V(t,1)==E_NTH_VERT(e,1));

  return(FALSE);
}


/* Triangulate the polygon defined by plist, and generating vertex p_id.
   Return list of added triangles in list add_ptr. Returns TRUE if 
   successful, FALSE otherwise. This is NOT a general triangulation routine,
   assumes polygon star relative to id
*/

int
smTriangulate(sm,id,plist)
SM *sm;
S_ID id;
LIST *plist;
{
    LIST *l,*prev,*t;
    FVECT v0,v1,v2,n,p;
    int is_tri,cut,t_id,e2,e1,enew;
    S_ID id0,id1,id2;
    double dp1,dp2;

    VSUB(p,SM_NTH_WV(sm,id),SM_VIEW_CENTER(sm));
    normalize(p);
    l = plist;

    enew = 0;
    cut = is_tri=  FALSE;
    prev = l;
    /* get v0,v1 */
    e1 = (int)LIST_DATA(l);
    id0 = E_NTH_VERT(e1,0);
    id1 = E_NTH_VERT(e1,1);
    VSUB(v0,SM_NTH_WV(sm,id0),SM_VIEW_CENTER(sm));
    normalize(v0);
    VSUB(v1,SM_NTH_WV(sm,id1),SM_VIEW_CENTER(sm));   
    normalize(v1);
    while(l)
    {
      l = LIST_NEXT(l);
      /* Get v2 */
      e2 = (int)LIST_DATA(l);
      id2 = E_NTH_VERT(e2,1);
      VSUB(v2,SM_NTH_WV(sm,id2),SM_VIEW_CENTER(sm));
      normalize(v2);
      if(LIST_NEXT(LIST_NEXT(l)) == prev)
      {/* Check if have a triangle */
	   is_tri = TRUE;	
	   break;
      }
      /* determine if v0-v1-v2 is convex:defined clockwise on the sphere-
       so switch orientation
       */
      VCROSS(n,v0,v2);
      normalize(n);
      dp1 = DOT(n,p);      
      if(dp1  <= 0.0)
      {
	  /* test if safe to cut off v0-v1-v2 by testing if p lies outside of
         triangle v0-v1-v2: if so, because plist is the star polygon around p,
          the new edge v2-v0 cannot intersect any existing edges
        */
	dp1 = DOT(n,v1);
	/* Distance from point s to plane is d = |N.p0s|/||N|| */
        /* sp1 = s-p0 for point on plane p0, a.(b+c) = a.b + a.c */
	if(dp1 >= EV_EPS)
        {
	    /* remove edges e1,e2 and add triangle id0,id1,id2 */
	    enew = 0;
	    t_id = smTriangulate_add_tri(sm,id0,id1,id2,e1,e2,&enew);
	    cut = TRUE;
	    /* Insert edge enew into the list, reuse prev list element */
	    LIST_NEXT(prev) = LIST_NEXT(l);
	    LIST_DATA(prev) = e1 = -enew;
	    /* If removing head of list- reset plist pointer */
	    if(l== plist)
	      plist = prev;
	    /* free list element for e2 */
	    LIST_NEXT(l)=NULL;
	    free_list(l);
	    l = prev;
	    VCOPY(v1,v2);
	    id1 = id2;
	    continue;
	}
    }
      VCOPY(v0,v1);
      VCOPY(v1,v2);
      id0 = id1;
      id1 = id2;
      e1 = e2;	
      /* check if gone around circular list without adding any
	 triangles: prevent infinite loop */
      if(l == plist)
      {
	if(LIST_NEXT(LIST_NEXT(l)) == prev)
	  {/* Check if have a triangle */
	    is_tri = TRUE;	
	    break;
	  }
	if(!cut)
	  break;
        cut = FALSE;
    }
      prev = l;
  }
  if(is_tri)
  {
      l = LIST_NEXT(l);
      enew = (int)LIST_DATA(l);
      t_id = smTriangulate_add_tri(sm,id0,id1,id2,e1,e2,&enew);
      free_list(l);
  }
  else
     if(!smTriangulate_quad(sm,l))
	goto Terror;
    /* Set triangle adjacencies based on edge adjacencies */ 
    FOR_ALL_EDGES(enew)
    {
      id0 = E_NTH_TRI(enew,0);
      id1 = E_NTH_TRI(enew,1);
	
      e1 = (T_WHICH_V(SM_NTH_TRI(sm,id0),E_NTH_VERT(enew,0))+2)%3;
      T_NTH_NBR(SM_NTH_TRI(sm,id0),e1) = id1;
      
      e2 = (T_WHICH_V(SM_NTH_TRI(sm,id1),E_NTH_VERT(enew,1))+2)%3;
      T_NTH_NBR(SM_NTH_TRI(sm,id1),e2) = id0;
  }

#ifdef DEBUG
#if DEBUG > 1
    {
      TRI *t0,*t1,*n;

    FOR_ALL_EDGES(enew)
    {
      id0 = E_NTH_TRI(enew,0);
      id1 = E_NTH_TRI(enew,1);
      t0 = SM_NTH_TRI(sm,id0);
      t1 = SM_NTH_TRI(sm,id1);
      if(T_NTH_NBR(t0,0) == T_NTH_NBR(t0,1)  || 
       T_NTH_NBR(t0,1) == T_NTH_NBR(t0,2)  ||
       T_NTH_NBR(t0,0) == T_NTH_NBR(t0,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t0,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t0,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t0,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t0,0));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t0,1));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t0,2));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");

      if(T_NTH_NBR(t1,0) == T_NTH_NBR(t1,1)  || 
       T_NTH_NBR(t1,1) == T_NTH_NBR(t1,2)  ||
       T_NTH_NBR(t1,0) == T_NTH_NBR(t1,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t1,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t1,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(t1,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t1,0));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t1,1));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(t1,2));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
  }

  }
#endif
#endif    

    return(TRUE);

Terror:	  

	  error(CONSISTENCY,"smTriangulate():Unable to triangulate\n");

}


void
smTris_swap_edge(sm,t_id,t1_id,e,e1,tn_id,tn1_id)
   SM *sm;
   int t_id,t1_id;
   int e,e1;
   int *tn_id,*tn1_id;
{
    S_ID verts[3];
    int enext,eprev,e1next,e1prev;
    TRI *n,*ta,*tb,*t,*t1;
    FVECT p1,p2,p3;
    int ta_id,tb_id;
    /* form new diagonal (e relative to t, and e1 relative to t1)
      defined by quadrilateral formed by t,t1- swap for the opposite diagonal
   */
    t = SM_NTH_TRI(sm,t_id);
    t1 = SM_NTH_TRI(sm,t1_id);
    enext = (e+1)%3;
    eprev = (e+2)%3;
    e1next = (e1+1)%3;
    e1prev = (e1+2)%3;
    verts[e] = T_NTH_V(t,e);
    verts[enext] = T_NTH_V(t,enext);
    verts[eprev] = T_NTH_V(t1,e1);
    ta_id = smAdd_tri(sm,verts[0],verts[1],verts[2],&ta);
#if 0
  fprintf(stderr,"Added tri %d %d %d %d\n",ta_id,T_NTH_V(ta,0),
	  T_NTH_V(ta,1), T_NTH_V(ta,2));
#endif
    verts[e1] = T_NTH_V(t1,e1);
    verts[e1next] = T_NTH_V(t1,e1next);
    verts[e1prev] = T_NTH_V(t,e);
    tb_id = smAdd_tri(sm,verts[0],verts[1],verts[2],&tb);
#if 0
  fprintf(stderr,"Added tri %d %d %d %d\n",tb_id,T_NTH_V(tb,0),
	  T_NTH_V(tb,1), T_NTH_V(tb,2));
#endif
    /* set the neighbors */
    T_NTH_NBR(ta,e) = T_NTH_NBR(t1,e1next);
    T_NTH_NBR(tb,e1) = T_NTH_NBR(t,enext);
    T_NTH_NBR(ta,enext)= tb_id;
    T_NTH_NBR(tb,e1next)= ta_id;
    T_NTH_NBR(ta,eprev)=T_NTH_NBR(t,eprev);
    T_NTH_NBR(tb,e1prev)=T_NTH_NBR(t1,e1prev);    

    /* Reset neighbor pointers of original neighbors */
    n = SM_NTH_TRI(sm,T_NTH_NBR(t,enext));
    T_NTH_NBR(n,T_NTH_NBR_PTR(t_id,n)) = tb_id;
    n = SM_NTH_TRI(sm,T_NTH_NBR(t,eprev));
    T_NTH_NBR(n,T_NTH_NBR_PTR(t_id,n)) = ta_id;
    n = SM_NTH_TRI(sm,T_NTH_NBR(t1,e1next));
    T_NTH_NBR(n,T_NTH_NBR_PTR(t1_id,n)) = ta_id;
    n = SM_NTH_TRI(sm,T_NTH_NBR(t1,e1prev));
    T_NTH_NBR(n,T_NTH_NBR_PTR(t1_id,n)) = tb_id;

    smDelete_tri(sm,t_id,t); 
    smDelete_tri(sm,t1_id,t1); 

#ifdef DEBUG
#if DEBUG > 1
    if(T_NTH_NBR(ta,0) == T_NTH_NBR(ta,1)  || 
       T_NTH_NBR(ta,1) == T_NTH_NBR(ta,2)  ||
       T_NTH_NBR(ta,0) == T_NTH_NBR(ta,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(ta,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(ta,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(ta,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(ta,0));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(ta,1));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(ta,2));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    if(T_NTH_NBR(ta,0) == T_NTH_NBR(ta,1)  || 
       T_NTH_NBR(ta,1) == T_NTH_NBR(ta,2)  ||
       T_NTH_NBR(ta,0) == T_NTH_NBR(ta,2)) 
      error(CONSISTENCY,"Invalid topology\n");

    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(tb,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(tb,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(tb,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(tb,0));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(tb,1));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
    n = SM_NTH_TRI(sm,T_NTH_NBR(tb,2));
    if(T_NTH_NBR(n,0) == T_NTH_NBR(n,1)  || 
       T_NTH_NBR(n,1) == T_NTH_NBR(n,2)  ||
       T_NTH_NBR(n,0) == T_NTH_NBR(n,2)) 
      error(CONSISTENCY,"Invalid topology\n");
    if(!T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,0))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,1))) ||
       !T_IS_VALID(SM_NTH_TRI(sm,T_NTH_NBR(n,2))))
      error(CONSISTENCY,"Invalid topology\n");
#endif
#endif
    *tn_id = ta_id;
    *tn1_id = tb_id;
    
    return;
}

/* Test the new set of triangles for Delaunay condition. 'Edges' contains
   all of the new edges added. The CCW triangle assoc with each edge is
   tested against the opposite vertex of the CW triangle. If the vertex
   lies inside the circle defined by the CCW triangle- the edge is swapped 
   for the opposite diagonal
*/
smFixEdges(sm)
   SM *sm;
{
    int e,t0_id,t1_id,e_new,e0,e1,e0_next,e1_next;
    S_ID v0_id,v1_id,v2_id,p_id;
    int  i,t0_nid,t1_nid;
    FVECT v0,v1,v2,p,np,v;
    TRI *t0,*t1;

    FOR_ALL_EDGES(e)
    {
	t0_id = E_NTH_TRI(e,0);
	t1_id = E_NTH_TRI(e,1);
#ifdef DEBUG
	if((t0_id==INVALID) || (t1_id==INVALID))
	  error(CONSISTENCY,"smFix_edges: Unassigned edge nbr\n");
#endif
	t0 = SM_NTH_TRI(sm,t0_id);
	t1 = SM_NTH_TRI(sm,t1_id);
	e0 = T_NTH_NBR_PTR(t1_id,t0);
	e1 = T_NTH_NBR_PTR(t0_id,t1);

	v0_id = E_NTH_VERT(e,0);
	v1_id = E_NTH_VERT(e,1);
	v2_id = T_NTH_V(t0,e0);
	p_id = T_NTH_V(t1,e1);

	smDir(sm,v0,v0_id);
	smDir(sm,v1,v1_id);
	smDir(sm,v2,v2_id);
	
	VSUB(p,SM_NTH_WV(sm,p_id),SM_VIEW_CENTER(sm));
	normalize(p);
	if(pt_in_cone(p,v2,v1,v0))
	{
	   smTris_swap_edge(sm,t0_id,t1_id,e0,e1,&t0_nid,&t1_nid);
	   /* Adjust the triangle pointers of the remaining edges to be
	      processed
            */
	    FOR_ALL_EDGES_FROM(e,i)
	    {
	      if(E_NTH_TRI(i,0)==t0_id || E_NTH_TRI(i,0)==t1_id)
	      {
		if(eIn_tri(i,SM_NTH_TRI(sm,t0_nid)))
		  SET_E_NTH_TRI(i,0,t0_nid);
		else
		  SET_E_NTH_TRI(i,0,t1_nid);
	      }

	      if(E_NTH_TRI(i,1)==t0_id || E_NTH_TRI(i,1)==t1_id)
	      {
		if(eIn_tri(i,SM_NTH_TRI(sm,t0_nid)))
		  SET_E_NTH_TRI(i,1,t0_nid);
		else
		  SET_E_NTH_TRI(i,1,t1_nid);
	      }
#ifdef DEBUG 
	      if(E_NTH_TRI(i,1) == E_NTH_TRI(i,0) )
		error(CONSISTENCY,"invalid edge\n");
#endif
	    }
	    t0_id = t0_nid;
	    t1_id = t1_nid;
	    e_new = eNew_edge();
	    SET_E_NTH_VERT(e_new,0,p_id);
	    SET_E_NTH_VERT(e_new,1,v2_id);
	    SET_E_NTH_TRI(e_new,0,t0_id);
	    SET_E_NTH_TRI(e_new,1,t1_id);
#ifdef DEBUG 
	      if(E_NTH_TRI(i,1) == E_NTH_TRI(i,0) )
		error(CONSISTENCY,"invalid edge\n");
#endif
	}
    }
}


smDelete_samp(sm,s_id)
SM *sm;
S_ID s_id;
{
  QUADTREE qt;
  S_ID *os;

  /* Mark as free */
  smUnalloc_samp(sm,s_id);

#ifdef DEBUG  
  SM_NTH_VERT(sm,s_id) = INVALID;
  /* fprintf(stderr,"deleting sample %d\n",s_id); */
#endif
  /* remove from its set */
  qt = SM_S_NTH_QT(sm,s_id);
  os = qtqueryset(qt);
  deletuelem(os, s_id);	/* delete obj from unsorted os, no questions */
}
/* Remove vertex "id" from the mesh- and retriangulate the resulting
   hole: Returns TRUE if successful, FALSE otherwise.
 */
int
smRemoveVertex(sm,id)
   SM *sm;
   S_ID id;
{
    LIST *b_list;
    /* generate list of edges that form the boundary of the
       polygon formed by the triangles adjacent to vertex 'id'*/
    b_list = smVertexStar(sm,id);
#if 0
 {int i;
    eputs("\n\n");
    for(i=1;i<=Ecnt;i++)
      fprintf(stderr,"%d verts %d %d tris  %d %d\n",
	      i,Edges[i].verts[0],Edges[i].verts[1],
	      Edges[i].tris[0],Edges[i].tris[1]);
 }
#endif

    /* Triangulate polygonal hole  */
    smTriangulate(sm,id,b_list);
    
    /* Fix up new triangles to be Delaunay*/

    smFixEdges(sm);
    smDelete_samp(sm,id);
    eClear_edges();
    return(TRUE);
}
   


 











