#ifndef lint
static const char	RCSid[] = "$Id: sm_qtree.c,v 3.16 2003/06/20 00:25:49 greg Exp $";
#endif
/*
*  sm_qtree.c: adapted from octree.c from radiance code
*/
/*
 *  octree.c - routines dealing with octrees and cubes.
 *
 *     7/28/85
 */

#include "standard.h"
#include "sm_flag.h"
#include "sm_geom.h"
#include "object.h"
#include "sm_qtree.h"

QUADTREE  *quad_block[QT_MAX_BLK];	  /* our quadtree */
static QUADTREE  quad_free_list = EMPTY;  /* freed octree nodes */
static QUADTREE  treetop = 0;		  /* next free node */
int32 *quad_flag= NULL;                    /* quadtree flags */


qtremovelast(qt,id)
QUADTREE qt;
OBJECT id;
{
#ifdef DEBUG
  if(qtqueryset(qt)[(*(qtqueryset(qt)))] != id)
    error(CONSISTENCY,"not removed\n");
#endif
  ((*(qtqueryset(qt)))--);
}
QUADTREE
qtAlloc()			/* allocate a quadtree */
{
    register QUADTREE  freet;

    if ((freet = quad_free_list) != EMPTY)
    {
	quad_free_list = QT_NTH_CHILD(freet, 0);
	return(freet);
    }
    freet = treetop;
    if (QT_BLOCK_INDEX(freet) == 0)
    {
	if (QT_BLOCK(freet) >= QT_MAX_BLK)
	   return(EMPTY);
	if ((quad_block[QT_BLOCK(freet)] = (QUADTREE *)malloc(
			QT_BLOCK_SIZE*4*sizeof(QUADTREE))) == NULL)
	   error(SYSTEM,"qtAlloc(): Unable to allocate memory\n");

	/* Realloc the per/node flags */
	quad_flag = (int32 *)realloc((void *)quad_flag,
			(QT_BLOCK(freet)+1)*((QT_BLOCK_SIZE+7)>>3));
	if (quad_flag == NULL)
	   error(SYSTEM,"qtAlloc(): Unable to allocate memory\n");
    }
    treetop += 4;
    return(freet);
}


qtClearAllFlags()		/* clear all quadtree branch flags */
{
  if (!treetop) 
    return;
  
  /* Clear the node flags*/
  bzero((char *)quad_flag, (QT_BLOCK(treetop-4)+1)*((QT_BLOCK_SIZE+7)>>3));
  /* Clear set flags */
  qtclearsetflags();
}

qtFree(qt)			/* free a quadtree */
   register QUADTREE  qt;
{
	register int  i;

	if (!QT_IS_TREE(qt))
	{
	  qtfreeleaf(qt);
	  return;
	}
	for (i = 0; i < 4; i++)
	   qtFree(QT_NTH_CHILD(qt, i));
	QT_NTH_CHILD(qt, 0) = quad_free_list;
	quad_free_list = qt;
}


qtDone()			/* free EVERYTHING */
{
	register int	i;
	
	qtfreeleaves();
	for (i = 0; i < QT_MAX_BLK; i++)
	{
	    if (quad_block[i] == NULL)
		break;
	    free((void *)quad_block[i]);
	    quad_block[i] = NULL;
	}
	/* Free the flags */
	if (i) free((void *)quad_flag);
	quad_flag = NULL;
	quad_free_list = EMPTY;
	treetop = 0;
}

/*
 * bary_parent(coord,i)  : Set coord to be value of pt at one level up in
 *                         quadtree, given it was the ith child
 * BCOORD coord[3];      :barycentric coordinates of point, also will hold
 *                        result as side effect
 * int i;                :designates which child coord was
 */
bary_parent(coord,i)
BCOORD coord[3];
int i;
{
  switch(i) {
  case 0:
    /* update bary for child */
    coord[0] = (coord[0] >> 1) + MAXBCOORD4;
    coord[1] >>= 1;
    coord[2] >>= 1;
    break;
  case 1:
    coord[0] >>= 1;
    coord[1]  = (coord[1] >> 1) + MAXBCOORD4;
    coord[2] >>= 1;
    break;
    
  case 2:
    coord[0] >>= 1;
    coord[1] >>= 1;
    coord[2] = (coord[2] >> 1) + MAXBCOORD4;
    break;
    
  case 3:
    coord[0] = MAXBCOORD4 - (coord[0] >> 1);
    coord[1] = MAXBCOORD4 - (coord[1] >> 1);
    coord[2] = MAXBCOORD4 - (coord[2] >> 1);
    break;
#ifdef DEBUG
  default:
    eputs("bary_parent():Invalid child\n");
    break;
#endif
  }
}

/*
 * bary_from_child(coord,child,next): Given that coord was the (child) child
 *                      of the quadtree node, calculate what the (next) child
 *                      is at the same level.
 * BCOORD coord[3];   :At current level (child)th child of node,will also hold
 *                     result as side effect
 * int child,next;    :Which child coord was (child), and which one
 *                     is now desired(next)
 */
bary_from_child(coord,child,next)
BCOORD coord[3];
int child,next;
{
#ifdef DEBUG
  if(child <0 || child > 3)
  {
    eputs("bary_from_child():Invalid child\n");
    return;
  }
  if(next <0 || next > 3)
  {
    eputs("bary_from_child():Invalid next\n");
    return;
  }
#endif
  if(next == child)
    return;

  switch(child){
  case 0:
      coord[0] = 0;
      coord[1] = MAXBCOORD2 - coord[1];
      coord[2] = MAXBCOORD2 - coord[2];
      break;
  case 1:
      coord[0] = MAXBCOORD2 - coord[0];
      coord[1] = 0;
      coord[2] = MAXBCOORD2 - coord[2];
      break;
  case 2:
      coord[0] = MAXBCOORD2 - coord[0];
      coord[1] = MAXBCOORD2 - coord[1];
      coord[2] = 0;
    break;
  case 3:
    switch(next){
    case 0:
      coord[0] = 0;
      coord[1] = MAXBCOORD2 - coord[1];
      coord[2] = MAXBCOORD2 - coord[2];
      break;
    case 1:
      coord[0] = MAXBCOORD2 - coord[0];
      coord[1] = 0;
      coord[2] = MAXBCOORD2 - coord[2];
      break;
    case 2:
      coord[0] = MAXBCOORD2 - coord[0];
      coord[1] = MAXBCOORD2 - coord[1];
      coord[2] = 0;
      break;
    }
    break;
  }
}

/*
 * int
 * bary_child(coord):  Return which child coord is of node at current level
 *                     As side effect recalculate coord at next level
 * BCOORD coord[3];   :At current level coordinates of point, will also set
 *                     to coords at next level as side effect
 */
int
bary_child(coord)
BCOORD coord[3];
{
  int i;

  if(coord[0] > MAXBCOORD4)
  { 
      /* update bary for child */
      coord[0] = (coord[0]<< 1) - MAXBCOORD2;
      coord[1] <<= 1;
      coord[2] <<= 1;
      return(0);
  }
  else
    if(coord[1] > MAXBCOORD4)
    {
      coord[0] <<= 1;
      coord[1] = (coord[1] << 1) - MAXBCOORD2;
      coord[2] <<= 1;
      return(1);
    }
    else
      if(coord[2] > MAXBCOORD4)
      {
	coord[0] <<= 1;
	coord[1] <<= 1;
	coord[2] = (coord[2] << 1) - MAXBCOORD2;
	return(2);
      }
      else
	 {
	   coord[0] = MAXBCOORD2 - (coord[0] << 1);
	   coord[1] = MAXBCOORD2 - (coord[1] << 1);
	   coord[2] = MAXBCOORD2 - (coord[2] << 1);
	   return(3);
	 }
}



QUADTREE
qtLocate(qt,bcoord)
QUADTREE qt;
BCOORD bcoord[3];
{
  int i;

  if(QT_IS_TREE(qt))
  {  
    i = bary_child(bcoord);

    return(qtLocate(QT_NTH_CHILD(qt,i),bcoord));
  }
  else
    return(qt);
}

int
move_to_nbr(b,db0,db1,db2,tptr)
BCOORD b[3];
BDIR db0,db1,db2;
double *tptr;
{
  double t,dt;
  BCOORD bc;
  int nbr;
  
  nbr = -1;
  *tptr = 0.0;
  /* Advance to next node */
  if(b[0]==0 && db0 < 0)
    return(0);
  if(b[1]==0 && db1 < 0)
    return(1);
  if(b[2]==0 && db2 < 0)
    return(2);
  if(db0 < 0)
   {
     t = ((double)b[0])/-db0;
     nbr = 0;
   }
  else
    t = MAXFLOAT;
  if(db1 < 0)
  {
     dt = ((double)b[1])/-db1;
    if( dt < t)
    {
      t = dt;
      nbr = 1;
    }
  }
  if(db2 < 0)
  {
     dt = ((double)b[2])/-db2;
       if( dt < t)
      {
	t = dt;
	nbr = 2;
      }
    }
  *tptr = t;
  return(nbr);
}

qtTrace_ray(qt,b,db0,db1,db2,nextptr,sign,sfactor,func,f)
   QUADTREE qt;
   BCOORD b[3];
   BDIR db0,db1,db2;
   int *nextptr;
   int sign,sfactor;
   FUNC func;
   int *f;
{
    int i,found;
    QUADTREE child;
    int nbr,next,w;
    double t;

    if(QT_IS_TREE(qt))
    {
      /* Find the appropriate child and reset the coord */
      i = bary_child(b);

      for(;;)
      {
	child = QT_NTH_CHILD(qt,i);
	if(i != 3)
	  qtTrace_ray(child,b,db0,db1,db2,nextptr,sign,sfactor+1,func,f);
	else
	  /* If the center cell- must flip direction signs */
	  qtTrace_ray(child,b,-db0,-db1,-db2,nextptr,1-sign,sfactor+1,func,f);

	if(QT_FLAG_IS_DONE(*f))
	  return;
	/* If in same block: traverse */
	if(i==3)
	  next = *nextptr;
	else
	  if(*nextptr == i)
	    next = 3;
	  else
	 {
	   /* reset the barycentric coordinates in the parents*/
	   bary_parent(b,i);
	   /* Else pop up to parent and traverse from there */
	   return(qt);
	 }
	bary_from_child(b,i,next);
	i = next;
      }
    }
    else
   {
#ifdef TEST_DRIVER
       if(Pick_cnt < 500)
	  Pick_q[Pick_cnt++]=qt;
#endif;
       F_FUNC(func)(qt,F_ARGS(func),f);
     if(QT_FLAG_IS_DONE(*f))
       return;
     /* Advance to next node */
     nbr = move_to_nbr(b,db0,db1,db2,&t);

     if(nbr==-1)
     {
       QT_FLAG_SET_DONE(*f); 
       return;
     }
     b[0] += (BCOORD)(t *db0);
     b[1] += (BCOORD)(t *db1);
     b[2] += (BCOORD)(t *db2);
     *nextptr = nbr;
     return;
     
   }
}    

#define TEST_INT(tl,th,d,q0,q1,h,l) \
                  tl=d>q0;th=d<q1; if(tl&&th) goto Lfunc_call; \
		 if(tl) if(l) goto Lfunc_call; else h=TRUE; \
		 if(th) if(h) goto Lfunc_call; else l = TRUE; \
		 if(th) if(h) goto Lfunc_call; else l = TRUE;

/* Leaf node: Do definitive test */
QUADTREE
qtLeaf_insert_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
		 scale, s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
unsigned int scale,s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{ 
  double db;
  unsigned int e0,e1,e2;
  char al,ah,bl,bh,cl,ch,testl,testh;
  char test_t0t1,test_t1t2,test_t2t0;
  BCOORD a,b,c;

  /* First check if can trivial accept: if vertex in cell */
  if(s0 & s1 & s2)
  {
    goto Lfunc_call;
  }
  /* Assumption: Could not trivial reject*/
  /* IF any coords lie on edge- must be in if couldnt reject */
  a = q1[0];b= q0[1]; c= q0[2];
  if(t0[0]== a || t0[1] == b || t0[2] == c)
  {
    goto Lfunc_call;
  }
  if(t1[0]== a || t1[1] == b || t1[2] == c)
  {
    goto Lfunc_call;
  }
  if(t2[0]== a || t2[1] == b || t2[2] == c)
  {
    goto Lfunc_call;
  }
  /* Test for edge crossings */
  /* Test t0t1,t1t2,t2t0 against a */
  /* Calculate edge crossings */
  e0  = (s0 ^ ((s0 >> 1) | (s0 << 2 & 4)));
  /* See if edge can be trivially rejected from intersetion testing */
  test_t0t1 = !(((s0 & 6)==0) || ((s1 & 6)==0)|| ((s2 & 6)==0) || 
       ((sq0 & 6)==6) ||((sq1 & 6)==6)|| ((sq2 & 6)==6));
  bl=bh=0;
  if(test_t0t1 && (e0 & 2) )
  {
    /* Must do double calculation to avoid overflow */
    db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
    TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  test_t1t2= !(((s0 & 3)==0) || ((s1 & 3)==0)|| ((s2 & 3)==0) || 
       ((sq0 & 3)==3) ||((sq1 & 3)==3)|| ((sq2 & 3)==3));
  if(test_t1t2 && (e0 & 1))
  {    /* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  test_t2t0 = !(((s0 & 5)==0) || ((s1 & 5)==0)|| ((s2 & 5)==0) || 
       ((sq0 & 5)==5) ||((sq1 & 5)==5)|| ((sq2 & 5)==5));
  if(test_t2t0 && (e0 & 4))
  {
      db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  e1  = (s1 ^ ((s1 >> 1) | (s1 << 2 & 4)));
  cl = ch = 0;
  if(test_t0t1 && (e1 & 2))
  {/* test t0t1 against b */
      db = t0[2] + dt10[2]*((double)(b-t0[1]))/dt10[1];
      TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
  if(test_t1t2 && (e1 & 1))
  {/* test t1t2 against b */
    db = t1[2] + dt21[2]*((double)(q0[1] - t1[1]))/dt21[1];
    TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
  if(test_t2t0 && (e1 & 4))
  {/* test t2t0 against b */
    db = t2[2] + dt02[2]*((double)(q0[1] - t2[1]))/dt02[1];
    TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
 
  /* test edges against c */
  e2  = (s2 ^ ((s2 >> 1) | (s2 << 2 & 4)));
  al = ah = 0;
  if(test_t0t1 && (e2 & 2))
  { /* test t0t1 against c */
    db = t0[0] + dt10[0]*((double)(c-t0[2]))/dt10[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
   }
  if(test_t1t2 && (e2 & 1))
  {
    db = t1[0] + dt21[0]*((double)(c - t1[2]))/dt21[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
  }
  if(test_t2t0 && (e2 & 4))
  { /* test t2t0 against c */
    db = t2[0] + dt02[0]*((double)(c - t2[2]))/dt02[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
  }
  /* Only remaining case is if some edge trivially rejected */
  if(!e0 || !e1 || !e2)
    return(qt);

  /* Only one/none got tested - pick either of other two/any two */
  /* Only need to test one edge */
  if(!test_t0t1 && (e0 & 2))
  {
     db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
     TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  if(!test_t1t2 && (e0 & 1))
    {/* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
   }
  if(!test_t2t0 && (e0 & 4))
  {
    db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
    TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }

  return(qt);

Lfunc_call:
  qt = f.func(f.argptr,root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	      s0,s1,s2,sq0,sq1,sq2,0,f,n);
  return(qt);

}



/* Leaf node: Do definitive test */
QUADTREE
qtLeaf_insert_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
		 scale, s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
unsigned int scale,s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{ 
  double db;
  unsigned int e0,e1,e2;
  BCOORD a,b,c;
  double p0[2], p1[2],cp;
  char al,ah,bl,bh,cl,ch;
  char testl,testh,test_t0t1,test_t1t2,test_t2t0;
  unsigned int ls0,ls1,ls2;
  
  
  /* May have gotten passed trivial reject if one/two vertices are ON */
  a = q1[0];b= q0[1]; c= q0[2];
  SIDES_LESS(t0,t1,t2,ls0,ls1,ls2,a,b,c);
  
  /* First check if can trivial accept: if vertex in cell */
  if(ls0 & ls1 & ls2)
  {
    goto Lfunc_call;
  }
  if(ls0==0 || ls1 == 0 || ls2 ==0)
    return(qt);
  /* Assumption: Could not trivial reject*/
  /* IF any coords lie on edge- must be in if couldnt reject */

  if(t0[0]== a || t0[1] == b || t0[2] == c)
  {
    goto Lfunc_call;
  }
  if(t1[0]== a || t1[1] == b || t1[2] == c)
  {
    goto Lfunc_call;
  }
  if(t2[0]== a || t2[1] == b || t2[2] == c)
  {
    goto Lfunc_call;
  }
  /* Test for edge crossings */
  /* Test t0t1 against a,b,c */
  /* First test if t0t1 can be trivially rejected */
  /* If both edges are outside bounds- dont need to test */
  
  /* Test t0t1,t1t2,t2t0 against a */
  test_t0t1 = !(((ls0 & 6)==0) || ((ls1 & 6)==0)|| ((ls2 & 6)==0) || 
       ((sq0 & 6)==0) ||((sq1 & 6)==0)|| ((sq2 & 6)==0));
  e0  = (ls0 ^ ((ls0 >> 1) | (ls0 << 2 & 4)));
  bl=bh=0;
  /* Test t0t1,t1t2,t2t0 against a */
  if(test_t0t1 && (e0 & 2) )
  {
      db = t0[1] + dt10[1]*((double)(a-t0[0])/dt10[0]);
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  test_t1t2= !(((ls0 & 3)==0) || ((ls1 & 3)==0)|| ((ls2 & 3)==0) || 
       ((sq0 & 3)==0) ||((sq1 & 3)==0)|| ((sq2 & 3)==0));
  if(test_t1t2 && (e0 & 1))
  {    /* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
 test_t2t0 = !(((ls0 & 5)==0) || ((ls1 & 5)==0)|| ((ls2 & 5)==0) || 
       ((sq0 & 5)==0) ||((sq1 & 5)==0)|| ((sq2 & 5)==0));
  if(test_t2t0 && (e0 & 4))
  {
      db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  cl = ch = 0;
  e1  = (ls1 ^ ((ls1 >> 1) | (ls1 << 2 & 4)));
  if(test_t0t1 && (e1 & 2))
  {/* test t0t1 against b */
      db = t0[2] + dt10[2]*(((double)(b-t0[1]))/dt10[1]);
      TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  if(test_t1t2 && (e1 & 1))
  {/* test t1t2 against b */
    db = t1[2] + dt21[2]*((double)(b - t1[1]))/dt21[1];
    TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  if(test_t2t0 && (e1 & 4))
  {/* test t2t0 against b */
    db = t2[2] + dt02[2]*((double)(b - t2[1]))/dt02[1];
    TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  al = ah = 0;
  e2  = (ls2 ^ ((ls2 >> 1) | (ls2 << 2 & 4)));
  if(test_t0t1 && (e2 & 2))
  { /* test t0t1 against c */
    db = t0[0] + dt10[0]*(((double)(c-t0[2]))/dt10[2]);
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
   }
  if(test_t1t2 && (e2 & 1))
  {
    db = t1[0] + dt21[0]*((double)(c - t1[2]))/dt21[2];
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
  }
  if(test_t2t0 && (e2 & 4))
  { /* test t2t0 against c */
    db = t2[0] + dt02[0]*((double)(c - t2[2]))/dt02[2];
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
  }
  /* Only remaining case is if some edge trivially rejected */
  if(!e0 || !e1 || !e2)
    return(qt);

  /* Only one/none got tested - pick either of other two/any two */
  /* Only need to test one edge */
  if(!test_t0t1 && (e0 & 2))
  {
     db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
     TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  if(!test_t1t2 && (e0 & 1))
    {/* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
   }
  if(!test_t2t0 && (e0 & 4))
  {
    db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
    TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }

  return(qt);
Lfunc_call:
  qt = f.func(f.argptr,root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	      s0,s1,s2,sq0,sq1,sq2,1,f,n);
  return(qt);
  }



/* ASSUMPTION: that triangle t0,t1,t2 intersects quad cell q0,q1,q2 */

/*-------q2--------- sq2
	/ \
s1     /sc \ s0
     qb_____qa
     / \   / \
\sq0/sa \ /sb \   /sq1
 \ q0____qc____q1/
  \             /
   \     s2    /
*/

int Find_id=0;

QUADTREE
qtInsert_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	    s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
BCOORD scale;
unsigned int s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];
  unsigned int sa,sb,sc;

  /* If a tree: trivial reject and recurse on potential children */
  if(QT_IS_TREE(qt))
  {
    /* Test against new a edge: if entirely to left: only need
       to visit child 0
    */
    a = q1[0] + scale;
    b = q0[1] + scale;
    c = q0[2] + scale;

    SIDES_GTR(t0,t1,t2,sa,sb,sc,a,b,c);

    if(sa==7)
    {
      vb[1] = q0[1];
      vc[2] = q0[2];
      vc[1] = b; 
      vb[2] = c;
      vb[0] = vc[0] = a;
      QT_NTH_CHILD(qt,0) = qtInsert_tri(root,QT_NTH_CHILD(qt,0),q0,vc,
	  vb,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,s1,s2,sq0,sb,sc,f,n+1);
      return(qt);
    }
   if(sb==7)
   {
     va[0] = q1[0]; 
     vc[2] = q0[2];
     va[1] = vc[1] = b; 
     va[2] = c;
     vc[0] = a;
     QT_NTH_CHILD(qt,1) = qtInsert_tri(root,QT_NTH_CHILD(qt,1),vc,q1,va,
	     t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n+1);
     return(qt);
   }
   if(sc==7)
   {
     va[0] = q1[0]; 
     vb[1] = q0[1];
     va[1] = b; 
     va[2] = vb[2] = c;
     vb[0] = a;
     QT_NTH_CHILD(qt,2) = qtInsert_tri(root,QT_NTH_CHILD(qt,2),vb,va,
       q2,t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n+1);
     return(qt);
   }

   va[0] = q1[0]; 
   vb[1] = q0[1];
   vc[2] = q0[2];
   va[1] = vc[1] = b; 
   va[2] = vb[2] = c;
   vb[0] = vc[0] = a;
   /* If outside of all edges: only need to Visit 3  */
   if(sa==0 && sb==0 && sc==0)
   {
      QT_NTH_CHILD(qt,3) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,3),va,vb,
       vc,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,sb,sc,s0,s1,s2,f,n+1);
      return(qt);
   }

   if(sa)
     QT_NTH_CHILD(qt,0) = qtInsert_tri(root,QT_NTH_CHILD(qt,0),q0,vc,vb,t0,
	  t1,t2,dt10,dt21,dt02,scale>>1,sa,s1,s2,sq0,sb,sc,f,n+1);
   if(sb)
     QT_NTH_CHILD(qt,1) = qtInsert_tri(root,QT_NTH_CHILD(qt,1),vc,q1,va,t0,
	      t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n+1);
   if(sc)
	QT_NTH_CHILD(qt,2) = qtInsert_tri(root,QT_NTH_CHILD(qt,2),vb,va,q2,t0,
	      t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n+1);
   QT_NTH_CHILD(qt,3) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,3),va,vb,vc,t0,
	     t1,t2,dt10,dt21,dt02,scale>>1,sa,sb,sc,s0,s1,s2,f,n+1);
   return(qt);
  }
  /* Leaf node: Do definitive test */
  else
    return(qt = qtLeaf_insert_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
			 scale,s0,s1,s2,sq0,sq1,sq2,f,n));
}


QUADTREE
qtInsert_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	    s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
BCOORD scale;
unsigned int s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];
  unsigned int sa,sb,sc;

  /* If a tree: trivial reject and recurse on potential children */
  if(QT_IS_TREE(qt))
  {
    /* Test against new a edge: if entirely to left: only need
       to visit child 0
    */
    a = q1[0] - scale;
    b = q0[1] - scale;
    c = q0[2] - scale;

    SIDES_GTR(t0,t1,t2,sa,sb,sc,a,b,c);

    if(sa==0)
    {
      vb[1] = q0[1];
      vc[2] = q0[2];
      vc[1] = b; 
      vb[2] = c;
      vb[0] = vc[0] = a;

      QT_NTH_CHILD(qt,0) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,0),q0,vc,
	  vb,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,s1,s2,sq0,sb,sc,f,n+1);
      return(qt);
    }
    if(sb==0)
    {
      va[0] = q1[0]; 
      vc[2] = q0[2];
      va[1] = vc[1] = b; 
      va[2] = c;
      vc[0] = a;
      QT_NTH_CHILD(qt,1) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,1),vc,q1,va,
         t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n+1);
      return(qt);
    }
    if(sc==0)
    {
      va[0] = q1[0]; 
      vb[1] = q0[1];
      va[1] = b; 
      va[2] = vb[2] = c;
      vb[0] = a;
      QT_NTH_CHILD(qt,2) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,2),vb,va,
	 q2,t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n+1);
      return(qt);
    }
    va[0] = q1[0]; 
    vb[1] = q0[1];
    vc[2] = q0[2];
    va[1] = vc[1] = b; 
    va[2] = vb[2] = c;
    vb[0] = vc[0] = a;
    /* If outside of all edges: only need to Visit 3  */
    if(sa==7 && sb==7 && sc==7)
    {
      QT_NTH_CHILD(qt,3) = qtInsert_tri(root,QT_NTH_CHILD(qt,3),va,vb,
	   vc,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,sb,sc,s0,s1,s2,f,n+1);
      return(qt);
    }
    if(sa!=7)
      QT_NTH_CHILD(qt,0) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,0),q0,vc,vb,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,sa,s1,s2,sq0,sb,sc,f,n+1);
    if(sb!=7)
      QT_NTH_CHILD(qt,1) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,1),vc,q1,va,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n+1);
    if(sc!=7)
      QT_NTH_CHILD(qt,2) = qtInsert_tri_rev(root,QT_NTH_CHILD(qt,2),vb,va,q2,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n+1);

    QT_NTH_CHILD(qt,3) = qtInsert_tri(root,QT_NTH_CHILD(qt,3),va,vb,vc,t0,t1,
	      t2,dt10,dt21,dt02,scale>>1,sa,sb,sc,s0,s1,s2,f,n+1);
    return(qt);
  }
  /* Leaf node: Do definitive test */
  else
    return(qt = qtLeaf_insert_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
			 scale,s0,s1,s2,sq0,sq1,sq2,f,n));
}




/*************************************************************************/
/* Visit code for applying functions: NOTE THIS IS THE SAME AS INSERT CODE:
  except sets flags: wanted insert to be as efficient as possible: so 
  broke into two sets of routines. Also traverses only to n levels.
*/

qtVisit_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	    s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
BCOORD scale;
unsigned int s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];
  unsigned int sa,sb,sc;

  if(n == 0)
    return;
  /* If a tree: trivial reject and recurse on potential children */
  if(QT_IS_TREE(qt))
  {
    QT_SET_FLAG(qt);

    /* Test against new a edge: if entirely to left: only need
       to visit child 0
    */
    a = q1[0] + scale;
    b = q0[1] + scale;
    c = q0[2] + scale;

    SIDES_GTR(t0,t1,t2,sa,sb,sc,a,b,c);

    if(sa==7)
    {
      vb[1] = q0[1];
      vc[2] = q0[2];
      vc[1] = b; 
      vb[2] = c;
      vb[0] = vc[0] = a;
      qtVisit_tri(root,QT_NTH_CHILD(qt,0),q0,vc,
	  vb,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,s1,s2,sq0,sb,sc,f,n-1);
      return;
    }
   if(sb==7)
   {
     va[0] = q1[0]; 
     vc[2] = q0[2];
     va[1] = vc[1] = b; 
     va[2] = c;
     vc[0] = a;
     qtVisit_tri(root,QT_NTH_CHILD(qt,1),vc,q1,va,
	     t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n-1);
     return;
   }
   if(sc==7)
   {
     va[0] = q1[0]; 
     vb[1] = q0[1];
     va[1] = b; 
     va[2] = vb[2] = c;
     vb[0] = a;
     qtVisit_tri(root,QT_NTH_CHILD(qt,2),vb,va,
       q2,t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n-1);
     return;
   }

   va[0] = q1[0]; 
   vb[1] = q0[1];
   vc[2] = q0[2];
   va[1] = vc[1] = b; 
   va[2] = vb[2] = c;
   vb[0] = vc[0] = a;
   /* If outside of all edges: only need to Visit 3  */
   if(sa==0 && sb==0 && sc==0)
   {
     qtVisit_tri_rev(root,QT_NTH_CHILD(qt,3),va,vb,
       vc,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,sb,sc,s0,s1,s2,f,n-1);
      return;
   }

   if(sa)
     qtVisit_tri(root,QT_NTH_CHILD(qt,0),q0,vc,vb,t0,
	  t1,t2,dt10,dt21,dt02,scale>>1,sa,s1,s2,sq0,sb,sc,f,n-1);
   if(sb)
     qtVisit_tri(root,QT_NTH_CHILD(qt,1),vc,q1,va,t0,
	      t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n-1);
   if(sc)
     qtVisit_tri(root,QT_NTH_CHILD(qt,2),vb,va,q2,t0,
	      t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n-1);
   qtVisit_tri_rev(root,QT_NTH_CHILD(qt,3),va,vb,vc,t0,
	     t1,t2,dt10,dt21,dt02,scale>>1,sa,sb,sc,s0,s1,s2,f,n-1);
  }
  /* Leaf node: Do definitive test */
  else
    qtLeaf_visit_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
			 scale,s0,s1,s2,sq0,sq1,sq2,f,n);

}


qtVisit_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,scale,
	    s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
BCOORD scale;
unsigned int s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];
  unsigned int sa,sb,sc;

  if(n==0)
    return;
  /* If a tree: trivial reject and recurse on potential children */
  if(QT_IS_TREE(qt))
  {
    QT_SET_FLAG(qt);
    /* Test against new a edge: if entirely to left: only need
       to visit child 0
    */
    a = q1[0] - scale;
    b = q0[1] - scale;
    c = q0[2] - scale;

    SIDES_GTR(t0,t1,t2,sa,sb,sc,a,b,c);

    if(sa==0)
    {
      vb[1] = q0[1];
      vc[2] = q0[2];
      vc[1] = b; 
      vb[2] = c;
      vb[0] = vc[0] = a;
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,0),q0,vc,
	  vb,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,s1,s2,sq0,sb,sc,f,n-1);
      return;
    }
    if(sb==0)
    {
      va[0] = q1[0]; 
      vc[2] = q0[2];
      va[1] = vc[1] = b; 
      va[2] = c;
      vc[0] = a;
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,1),vc,q1,va,
         t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n-1);
      return;
    }
    if(sc==0)
    {
      va[0] = q1[0]; 
      vb[1] = q0[1];
      va[1] = b; 
      va[2] = vb[2] = c;
      vb[0] = a;
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,2),vb,va,
	 q2,t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n-1);
      return;
    }
    va[0] = q1[0]; 
    vb[1] = q0[1];
    vc[2] = q0[2];
    va[1] = vc[1] = b; 
    va[2] = vb[2] = c;
    vb[0] = vc[0] = a;
    /* If outside of all edges: only need to Visit 3  */
    if(sa==7 && sb==7 && sc==7)
    {
     qtVisit_tri(root,QT_NTH_CHILD(qt,3),va,vb,
	   vc,t0,t1,t2,dt10,dt21,dt02, scale>>1,sa,sb,sc,s0,s1,s2,f,n-1);
      return;
    }
    if(sa!=7)
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,0),q0,vc,vb,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,sa,s1,s2,sq0,sb,sc,f,n-1);
    if(sb!=7)
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,1),vc,q1,va,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,sb,s2,sa,sq1,sc,f,n-1);
    if(sc!=7)
      qtVisit_tri_rev(root,QT_NTH_CHILD(qt,2),vb,va,q2,
	   t0,t1,t2,dt10,dt21,dt02,scale>>1,s0,s1,sc,sa,sb,sq2,f,n-1);

    qtVisit_tri(root,QT_NTH_CHILD(qt,3),va,vb,vc,t0,t1,
	      t2,dt10,dt21,dt02,scale>>1,sa,sb,sc,s0,s1,s2,f,n-1);
    return;
  }
  /* Leaf node: Do definitive test */
  else
    qtLeaf_visit_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
			 scale,s0,s1,s2,sq0,sq1,sq2,f,n);
}



qtLeaf_visit_tri(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
		 scale, s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
unsigned int scale,s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{ 
  double db;
  unsigned int e0,e1,e2;
  char al,ah,bl,bh,cl,ch,testl,testh;
  char test_t0t1,test_t1t2,test_t2t0;
  BCOORD a,b,c;
  
  /* First check if can trivial accept: if vertex in cell */
  if(s0 & s1 & s2)
    goto Lfunc_call;

  /* Assumption: Could not trivial reject*/
  /* IF any coords lie on edge- must be in if couldnt reject */
  a = q1[0];b= q0[1]; c= q0[2];
  if(t0[0]== a || t0[1] == b || t0[2] == c)
    goto Lfunc_call;
  if(t1[0]== a || t1[1] == b || t1[2] == c)
    goto Lfunc_call;
  if(t2[0]== a || t2[1] == b || t2[2] == c)
    goto Lfunc_call;
  
  /* Test for edge crossings */
  /* Test t0t1,t1t2,t2t0 against a */
  /* Calculate edge crossings */
  e0  = (s0 ^ ((s0 >> 1) | (s0 << 2 & 4)));
  /* See if edge can be trivially rejected from intersetion testing */
  test_t0t1 = !(((s0 & 6)==0) || ((s1 & 6)==0)|| ((s2 & 6)==0) || 
       ((sq0 & 6)==6) ||((sq1 & 6)==6)|| ((sq2 & 6)==6));
  bl=bh=0;
  if(test_t0t1 && (e0 & 2) )
  {
    /* Must do double calculation to avoid overflow */
    db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
    TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  test_t1t2= !(((s0 & 3)==0) || ((s1 & 3)==0)|| ((s2 & 3)==0) || 
       ((sq0 & 3)==3) ||((sq1 & 3)==3)|| ((sq2 & 3)==3));
  if(test_t1t2 && (e0 & 1))
  {    /* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  test_t2t0 = !(((s0 & 5)==0) || ((s1 & 5)==0)|| ((s2 & 5)==0) || 
       ((sq0 & 5)==5) ||((sq1 & 5)==5)|| ((sq2 & 5)==5));
  if(test_t2t0 && (e0 & 4))
  {
      db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  e1  = (s1 ^ ((s1 >> 1) | (s1 << 2 & 4)));
  cl = ch = 0;
  if(test_t0t1 && (e1 & 2))
  {/* test t0t1 against b */
      db = t0[2] + dt10[2]*((double)(b-t0[1]))/dt10[1];
      TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
  if(test_t1t2 && (e1 & 1))
  {/* test t1t2 against b */
    db = t1[2] + dt21[2]*((double)(q0[1] - t1[1]))/dt21[1];
    TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
  if(test_t2t0 && (e1 & 4))
  {/* test t2t0 against b */
    db = t2[2] + dt02[2]*((double)(q0[1] - t2[1]))/dt02[1];
    TEST_INT(testl,testh,db,c,q2[2],cl,ch)
  }
 
  /* test edges against c */
  e2  = (s2 ^ ((s2 >> 1) | (s2 << 2 & 4)));
  al = ah = 0;
  if(test_t0t1 && (e2 & 2))
  { /* test t0t1 against c */
    db = t0[0] + dt10[0]*((double)(c-t0[2]))/dt10[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
   }
  if(test_t1t2 && (e2 & 1))
  {
    db = t1[0] + dt21[0]*((double)(c - t1[2]))/dt21[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
  }
  if(test_t2t0 && (e2 & 4))
  { /* test t2t0 against c */
    db = t2[0] + dt02[0]*((double)(c - t2[2]))/dt02[2];
    TEST_INT(testl,testh,db,a,q0[0],al,ah)
  }
  /* Only remaining case is if some edge trivially rejected */
  if(!e0 || !e1 || !e2)
    return;

  /* Only one/none got tested - pick either of other two/any two */
  /* Only need to test one edge */
  if(!test_t0t1 && (e0 & 2))
  {
     db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
     TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }
  if(!test_t1t2 && (e0 & 1))
    {/* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,b,q1[1],bl,bh)
   }
  if(!test_t2t0 && (e0 & 4))
  {
    db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
    TEST_INT(testl,testh,db,b,q1[1],bl,bh)
  }

  return;

Lfunc_call:

  f.func(f.argptr,root,qt,n);
  if(!QT_IS_EMPTY(qt))
    QT_LEAF_SET_FLAG(qt);

}



/* Leaf node: Do definitive test */
QUADTREE
qtLeaf_visit_tri_rev(root,qt,q0,q1,q2,t0,t1,t2,dt10,dt21,dt02,
		 scale, s0,s1,s2,sq0,sq1,sq2,f,n)
int root;
QUADTREE qt;
BCOORD q0[3],q1[3],q2[3];
BCOORD t0[3],t1[3],t2[3];
BCOORD dt10[3],dt21[3],dt02[3];
unsigned int scale,s0,s1,s2,sq0,sq1,sq2;
FUNC f;
int n;
{ 
  double db;
  unsigned int e0,e1,e2;
  BCOORD a,b,c;
  double p0[2], p1[2],cp;
  char al,ah,bl,bh,cl,ch;
  char testl,testh,test_t0t1,test_t1t2,test_t2t0;
  unsigned int ls0,ls1,ls2;
  
  /* May have gotten passed trivial reject if one/two vertices are ON */
  a = q1[0];b= q0[1]; c= q0[2];
  SIDES_LESS(t0,t1,t2,ls0,ls1,ls2,a,b,c);
  
  /* First check if can trivial accept: if vertex in cell */
  if(ls0 & ls1 & ls2)
    goto Lfunc_call;

  if(ls0==0 || ls1 == 0 || ls2 ==0)
    return;
  /* Assumption: Could not trivial reject*/
  /* IF any coords lie on edge- must be in if couldnt reject */

  if(t0[0]== a || t0[1] == b || t0[2] == c)
    goto Lfunc_call;
  if(t1[0]== a || t1[1] == b || t1[2] == c)
    goto Lfunc_call;
  if(t2[0]== a || t2[1] == b || t2[2] == c)
    goto Lfunc_call;

  /* Test for edge crossings */
  /* Test t0t1 against a,b,c */
  /* First test if t0t1 can be trivially rejected */
  /* If both edges are outside bounds- dont need to test */
  
  /* Test t0t1,t1t2,t2t0 against a */
  test_t0t1 = !(((ls0 & 6)==0) || ((ls1 & 6)==0)|| ((ls2 & 6)==0) || 
       ((sq0 & 6)==0) ||((sq1 & 6)==0)|| ((sq2 & 6)==0));
  e0  = (ls0 ^ ((ls0 >> 1) | (ls0 << 2 & 4)));
  bl=bh=0;
  /* Test t0t1,t1t2,t2t0 against a */
  if(test_t0t1 && (e0 & 2) )
  {
      db = t0[1] + dt10[1]*((double)(a-t0[0])/dt10[0]);
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  test_t1t2= !(((ls0 & 3)==0) || ((ls1 & 3)==0)|| ((ls2 & 3)==0) || 
       ((sq0 & 3)==0) ||((sq1 & 3)==0)|| ((sq2 & 3)==0));
  if(test_t1t2 && (e0 & 1))
  {    /* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  test_t2t0 = !(((ls0 & 5)==0) || ((ls1 & 5)==0)|| ((ls2 & 5)==0) || 
       ((sq0 & 5)==0) ||((sq1 & 5)==0)|| ((sq2 & 5)==0));
  if(test_t2t0 && (e0 & 4))
  {
      db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  cl = ch = 0;
  e1  = (ls1 ^ ((ls1 >> 1) | (ls1 << 2 & 4)));
  if(test_t0t1 && (e1 & 2))
  {/* test t0t1 against b */
      db = t0[2] + dt10[2]*(((double)(b-t0[1]))/dt10[1]);
      TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  if(test_t1t2 && (e1 & 1))
  {/* test t1t2 against b */
    db = t1[2] + dt21[2]*((double)(b - t1[1]))/dt21[1];
    TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  if(test_t2t0 && (e1 & 4))
  {/* test t2t0 against b */
    db = t2[2] + dt02[2]*((double)(b - t2[1]))/dt02[1];
    TEST_INT(testl,testh,db,q1[2],c,cl,ch)
  }
  al = ah = 0;
  e2  = (ls2 ^ ((ls2 >> 1) | (ls2 << 2 & 4)));
  if(test_t0t1 && (e2 & 2))
  { /* test t0t1 against c */
    db = t0[0] + dt10[0]*(((double)(c-t0[2]))/dt10[2]);
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
   }
  if(test_t1t2 && (e2 & 1))
  {
    db = t1[0] + dt21[0]*((double)(c - t1[2]))/dt21[2];
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
  }
  if(test_t2t0 && (e2 & 4))
  { /* test t2t0 against c */
    db = t2[0] + dt02[0]*((double)(c - t2[2]))/dt02[2];
    TEST_INT(testl,testh,db,q0[0],a,al,ah)
  }
  /* Only remaining case is if some edge trivially rejected */
  if(!e0 || !e1 || !e2)
    return;

  /* Only one/none got tested - pick either of other two/any two */
  /* Only need to test one edge */
  if(!test_t0t1 && (e0 & 2))
  {
     db = t0[1] + dt10[1]*((double)(a-t0[0]))/dt10[0];
     TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }
  if(!test_t1t2 && (e0 & 1))
    {/* test t1t2 against a */
      db = t1[1] + dt21[1]*((double)(a - t1[0]))/dt21[0];
      TEST_INT(testl,testh,db,q1[1],b,bl,bh)
   }
  if(!test_t2t0 && (e0 & 4))
  {
    db = t2[1] + dt02[1]*((double)(a - t2[0]))/dt02[0];
    TEST_INT(testl,testh,db,q1[1],b,bl,bh)
  }

  return;

Lfunc_call:
  f.func(f.argptr,root,qt,n);
  if(!QT_IS_EMPTY(qt))
    QT_LEAF_SET_FLAG(qt);
}


QUADTREE
qtInsert_point(root,qt,parent,q0,q1,q2,bc,scale,f,n,doneptr)
int root;
QUADTREE qt,parent;
BCOORD q0[3],q1[3],q2[3],bc[3],scale; 
FUNC f;
int n,*doneptr;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];

  if(QT_IS_TREE(qt))
  {  
    a = q1[0] + scale;
    b = q0[1] + scale;
    c = q0[2] + scale;
    if(bc[0] > a)
    {
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      QT_NTH_CHILD(qt,0) = qtInsert_point(root,QT_NTH_CHILD(qt,0),
				  qt,q0,vc,vb,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	      
      return(qt);
    }
    if(bc[1] > b)
    {
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      va[0] = q1[0]; va[1] = b; va[2] = c;
      QT_NTH_CHILD(qt,1) =qtInsert_point(root,QT_NTH_CHILD(qt,1),
				 qt,vc,q1,va,bc,scale >>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	      
      return(qt);
    }
    if(bc[2] > c)
    {
      va[0] = q1[0]; va[1] = b; va[2] = c;
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      QT_NTH_CHILD(qt,2) = qtInsert_point(root,QT_NTH_CHILD(qt,2),
				  qt,vb,va,q2,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	
      return(qt);
    }
    else
    {
      va[0] = q1[0]; va[1] = b; va[2] = c;
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      QT_NTH_CHILD(qt,3) = qtInsert_point_rev(root,QT_NTH_CHILD(qt,3),
				      qt,va,vb,vc,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
      }	
      return(qt);
    }
  }
  else
  {
    qt = f.func(f.argptr,root,qt,parent,q0,q1,q2,bc,scale,0,f,n,doneptr);
    return(qt);
  }
}


QUADTREE
qtInsert_point_rev(root,qt,parent,q0,q1,q2,bc,scale,f,n,doneptr)
int root;
QUADTREE qt,parent;
BCOORD q0[3],q1[3],q2[3],bc[3]; 
BCOORD scale;
FUNC f;
int n,*doneptr;
{
  BCOORD a,b,c;
  BCOORD va[3],vb[3],vc[3];

  if(QT_IS_TREE(qt))
  {  
    a = q1[0] - scale;
    b = q0[1] - scale;
    c = q0[2] - scale;
    if(bc[0] < a)
    {
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      QT_NTH_CHILD(qt,0) = qtInsert_point_rev(root,QT_NTH_CHILD(qt,0),
				  qt,q0,vc,vb,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	
      return(qt);
    }
    if(bc[1] < b)
    {
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      va[0] = q1[0]; va[1] = b; va[2] = c;
      QT_NTH_CHILD(qt,1) = qtInsert_point_rev(root,QT_NTH_CHILD(qt,1),
				   qt,vc,q1,va,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	
      return(qt);
    }
    if(bc[2] < c)
    {
      va[0] = q1[0]; va[1] = b; va[2] = c;
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      QT_NTH_CHILD(qt,2) = qtInsert_point_rev(root,QT_NTH_CHILD(qt,2),
				qt,vb,va,q2,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,3),f,doneptr);
      }	
      return(qt);
    }
    else
    {
      va[0] = q1[0]; va[1] = b; va[2] = c;
      vb[0] = a; vb[1] = q0[1]; vb[2] = c;
      vc[0] = a; vc[1] = b; vc[2] = q0[2];
      QT_NTH_CHILD(qt,3) = qtInsert_point(root,QT_NTH_CHILD(qt,3),
				   qt,va,vb,vc,bc,scale>>1,f,n+1,doneptr);
      if(!*doneptr)
      {
	f.func_after(f.argptr,QT_NTH_CHILD(qt,0),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,1),f,doneptr);
	if(!*doneptr)
	  f.func_after(f.argptr,QT_NTH_CHILD(qt,2),f,doneptr);
      }	
      return(qt);
    }
  }
  else
  {
    qt = f.func(f.argptr,root,qt,parent,q0,q1,q2,bc,scale,1,f,n,doneptr);
    return(qt);
  }
}











