/* RCSid: $Id: sm_qtree.h,v 3.14 2003/05/15 05:13:35 greg Exp $ */
/*
 *  sm_qtree.h - header file for routines using spherical quadtrees.
 *
 *    adapted from octree.h 
 */

/*
 *	An quadtree is expressed as an integer which is either
 *	an index to 4 other nodes, the empty tree, or an index
 *	to a set of objects.  If the quadtree has a value:
 *
 *		> -1:	it is an index to four other nodes.
 *
 *		-1:	it is empty
 *
 *		< -1:	it is an index to a set of objects
 */

typedef struct _FUNC {
    int (*func)();
    int (*func_after)();
    int *argptr;
}FUNC;

#define F_FUNC(f) (f.func)
#define F_FUNC2(f) (f.func_after)
#define F_ARGS(f) (f.argptr)
#define  QUADTREE		int

#define  EMPTY		(-1)

#define  QT_IS_EMPTY(qt)	((qt) == EMPTY)
#define  QT_IS_LEAF(qt)		((qt) < EMPTY)
#define  QT_IS_TREE(qt)		((qt) > EMPTY)

#define  QT_INDEX(qt)		(-(qt)-2)	/* quadtree node from set */
#define  QT_SET_INDEX(i)	(-((i)+2))	/* object set from node */
#define  QT_BLOCK_SIZE		1024		/* quadtrees in block */
#define  QT_BLOCK(qt)		((qt)>>10)	/* quadtree block index */
#define  QT_BLOCK_INDEX(qt)   (((qt)&0x3ff)<<2)	/* quadtree index in block */


#ifndef  QT_MAX_BLK
#ifdef  SMLMEM
#define  QT_MAX_BLK	2047		/* maximum quadtree block */
#else
#define  QT_MAX_BLK	16383		/* maximum quadtree block */
#endif
#endif


#define  QT_NTH_CHILD(qt,br) (quad_block[QT_BLOCK(qt)][QT_BLOCK_INDEX(qt)+br])
#define  QT_NTH_CHILD_PTR(qt,br) (&QT_NTH_CHILD(qt,br))
#define  QT_CLEAR_CHILDREN(qt) (QT_NTH_CHILD(qt,0)=QT_NTH_CHILD(qt,1)= \
			QT_NTH_CHILD(qt,2)=QT_NTH_CHILD(qt,3)=EMPTY)


/* QUADTREE NODE FLAGS */
#define QT_IS_FLAG(qt)		IS_FLAG(quad_flag,qt)
#define QT_SET_FLAG(qt)		SET_FLAG(quad_flag,qt)
#define QT_CLR_FLAG(qt)		CLR_FLAG(quad_flag,qt)
#define QT_LEAF_IS_FLAG(qt)	IS_FLAG(qtsetflag,QT_INDEX(qt))
#define QT_LEAF_SET_FLAG(qt)	SET_FLAG(qtsetflag,QT_INDEX(qt))
#define QT_LEAF_CLR_FLAG(qt)	CLR_FLAG(qtsetflag,QT_INDEX(qt))

/* OBJECT SET CODE */
#define QT_SET_CNT(s)          ((s)[0])
#define QT_SET_NTH_ELEM(s,n)   ((s)[(n)])   

#define QT_CLEAR_SET(s)        ((s)[0] = 0)
#define QT_SET_NEXT_ELEM(p)    (*(++p))
#define QT_SET_PTR(s)          (&((s)[1]))


#define QT_MAXSET        511
#define MAXCSET          2*QT_MAXSET
#define QT_MAXCSET       MAXCSET
#ifndef QT_SET_THRESHOLD
#define QT_SET_THRESHOLD 32  
#endif

#ifndef QT_MAX_LEVELS
#define QT_MAX_LEVELS     10
#endif

#define QT_FILL_THRESHOLD 2
#define QT_EXPAND   8
#define QT_COMPRESS 16
#define QT_DONE 32
#define QT_MODIFIED 64

#define QT_FLAG_FILL_TRI(f)  (((f)&0x7) == QT_FILL_THRESHOLD)
#define QT_FLAG_UPDATE(f)    ((f)& (QT_EXPAND | QT_COMPRESS))
#define QT_FLAG_IS_DONE(f)   ((f)& QT_DONE)
#define QT_FLAG_SET_DONE(f)   ((f) |= QT_DONE)
#define QT_FLAG_IS_MODIFIED(f)   ((f)& QT_MODIFIED)
#define QT_FLAG_SET_MODIFIED(f)   ((f) |= QT_MODIFIED)

#define qtSubdivide(qt) (qt = qtAlloc(),QT_CLEAR_CHILDREN(qt))
#define qtSubdivide_tri(v0,v1,v2,a,b,c) (EDGE_MIDPOINT(a,v1,v2), \
					 EDGE_MIDPOINT(b,v2,v0), \
					 EDGE_MIDPOINT(c,v0,v1))
 
extern QUADTREE  qtnewleaf(), qtaddelem(), qtdelelem();

extern QUADTREE  *quad_block[QT_MAX_BLK];	/* quadtree blocks */
extern int4  *quad_flag;			/* zeroeth quadtree flag */

extern OBJECT	**qtsettab;		/* quadtree leaf node table */
extern QUADTREE  qtnumsets;		/* number of used set indices */
extern int4   *qtsetflag;
#ifdef DEBUG
extern OBJECT	*qtqueryset();
#else
#define qtqueryset(qt)	(qtsettab[QT_SET_INDEX(qt)])
#endif
#if 0
#define qtremovelast(qt) ((*(qtqueryset(qt)))--)
#endif
#define qtinset(qt,id)	inset(qtqueryset(qt),id)
#define qtgetset(os,qt)	setcopy(os,qtqueryset(qt))


#define SIDES_GTR(b0,b1,b2,s0,s1,s2,a,b,c) \
    (s0 = ((b0[0] > a?4:0) | (b1[0] > a?2:0) | (b2[0] > a?1:0)), \
     s1 = ((b0[1] > b?4:0) | (b1[1] > b?2:0) | (b2[1] > b?1:0)), \
     s2 = ((b0[2] > c?4:0) | (b1[2] > c?2:0) | (b2[2] > c?1:0)))

#define SIDES_LESS(b0,b1,b2,s0,s1,s2,a,b,c) \
    (s0 = ((b0[0] < a?4:0) | (b1[0] < a?2:0) | (b2[0] < a?1:0)), \
     s1 = ((b0[1] < b?4:0) | (b1[1] < b?2:0) | (b2[1] < b?1:0)), \
     s2 = ((b0[2] < c?4:0) | (b1[2] < c?2:0) | (b2[2] < c?1:0)))

/* 
QUADTREE qtRoot_point_locate(qt,q0,q1,q2,peq,pt,r0,r1,r2)
   QUADTREE qt;
   FVECT q0,q1,q2;
   FPEQ peq;
   FVECT pt;
   FVECT r0,r1,r2;

   Return the quadtree node containing pt. It is assumed that pt is in
   the root node qt with ws vertices q0,q1,q2 and plane equation peq.
   If r0 != NULL will return coordinates of node in (r0,r1,r2). 
*/

extern QUADTREE qtRoot_point_locate();
extern QUADTREE qtRoot_add_tri();
extern QUADTREE qtRoot_remove_tri();
extern QUADTREE qtAdd_tri();
extern QUADTREE qtRoot_visit_tri_edges();
extern QUADTREE qtRoot_trace_ray();






