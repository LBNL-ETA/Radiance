#ifndef lint
static const char RCSid[] = "$Id: triangulate.c,v 2.7 2021/04/19 20:37:04 greg Exp $";
#endif
/*
 *  triangulate.c
 *  
 *  Adapted by Greg Ward on 1/23/14.
 *  Fixes for polygons with seams/holes and co-linear vertices added
 *  by Nathaniel Jones on 12/21/16.
 *  Copyright 2016 Anyhere Software. All rights reserved.
 *
 */

/* COTD Entry submitted by John W. Ratcliff [jratcliff@verant.com]

// ** THIS IS A CODE SNIPPET WHICH WILL EFFICIEINTLY TRIANGULATE ANY
// ** POLYGON/CONTOUR (without holes) AS A STATIC CLASS.  THIS SNIPPET
// ** IS COMPRISED OF 3 FILES, TRIANGULATE.H, THE HEADER FILE FOR THE
// ** TRIANGULATE BASE CLASS, TRIANGULATE.CPP, THE IMPLEMENTATION OF
// ** THE TRIANGULATE BASE CLASS, AND TEST.CPP, A SMALL TEST PROGRAM
// ** DEMONSTRATING THE USAGE OF THE TRIANGULATOR.  THE TRIANGULATE
// ** BASE CLASS ALSO PROVIDES TWO USEFUL HELPER METHODS, ONE WHICH
// ** COMPUTES THE AREA OF A POLYGON, AND ANOTHER WHICH DOES AN EFFICENT
// ** POINT IN A TRIANGLE TEST.
// ** SUBMITTED BY JOHN W. RATCLIFF (jratcliff@verant.com) July 22, 2000
*/

#include <stdio.h>
#include <stdlib.h>
#include "triangulate.h"

#ifndef true
#define true	1
#define false	0
#endif

#define EPSILON		0.0000000001

static int
polySnip(const Vert2_list *contour, int u, int v, int w, int n, int *V)
{
  int p;
  double Ax, Ay, Bx, By, Cx, Cy, Px, Py, cross;

  Ax = contour->v[V[u]].mX;
  Ay = contour->v[V[u]].mY;

  Bx = contour->v[V[v]].mX;
  By = contour->v[V[v]].mY;

  Cx = contour->v[V[w]].mX;
  Cy = contour->v[V[w]].mY;

  cross = ((Bx - Ax)*(Cy - Ay)) - ((By - Ay)*(Cx - Ax));
  if (cross < EPSILON)
  	return cross > -EPSILON ? -1 : false; /* Negative if colinear points */

  for (p=0;p<n;p++)
  {
    if( (p == u) | (p == v) | (p == w) ) continue;
    Px = contour->v[V[p]].mX;
    Py = contour->v[V[p]].mY;
    if ((Px == Ax) & (Py == Ay) || (Px == Bx) & (Py == By) ||
		(Px == Cx) & (Py == Cy)) continue; /* Handle donuts */
    if (insideTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
  }

  return true;
}

Vert2_list *
polyAlloc(int nv)
{
	Vert2_list	*pnew;

	if (nv < 3) return NULL;
	
	pnew = (Vert2_list *)malloc(sizeof(Vert2_list) + sizeof(Vert2)*(nv-3));
	if (pnew == NULL) return NULL;
	pnew->nv = nv;
	pnew->p = NULL;
	
	return pnew;
}

/*
  Area is positive if vertices listed counter-clockwise, negative if clockwise
 */
double
polyArea(const Vert2_list *contour)
{
  double A=0.0;
  int	p, q;

  for(p = contour->nv-1, q = 0; q < contour->nv; p=q++)
  {
    A += contour->v[p].mX*contour->v[q].mY - contour->v[q].mX*contour->v[p].mY;
  }
  return A*0.5;
}

/*
  InsideTriangle decides if a point P is Inside of the triangle
  defined by A, B, C.
*/
int
insideTriangle(double Ax, double Ay,
                      double Bx, double By,
                      double Cx, double Cy,
                      double Px, double Py)
{
  double ax, ay, bx, by, cx, cy, apx, apy, bpx, bpy, cpx, cpy;
  double cCROSSap, bCROSScp, aCROSSbp;

  ax = Cx - Bx;  ay = Cy - By;
  bx = Ax - Cx;  by = Ay - Cy;
  cx = Bx - Ax;  cy = By - Ay;
  apx= Px - Ax;  apy= Py - Ay;
  bpx= Px - Bx;  bpy= Py - By;
  cpx= Px - Cx;  cpy= Py - Cy;

  aCROSSbp = ax*bpy - ay*bpx;
  cCROSSap = cx*apy - cy*apx;
  bCROSScp = bx*cpy - by*cpx;

  return ((aCROSSbp >= 0.0) & (bCROSScp >= 0.0) & (cCROSSap >= 0.0));
};

int
polyTriangulate(const Vert2_list *contour, tri_out_t *cb)
{
  /* allocate and initialize list of Vertices in polygon */

  int nv, u, v, w, count, result;
  int *V;

  if ( contour->nv < 3 ) return false;

  V = (int *)malloc(sizeof(int)*contour->nv);
  if (V == NULL) return false;

  /* we want a counter-clockwise polygon in V */

  if ( polyArea(contour) > 0.0 )
    for (v=0; v<contour->nv; v++) V[v] = v;
  else
    for (v=0; v<contour->nv; v++) V[v] = (contour->nv-1)-v;

  nv = contour->nv;

  /*  remove nv-2 Vertices, creating 1 triangle every time */
  count = 2*nv;   /* error detection */

  v = nv-1;
  while (nv > 2)
  {
    /* if we loop, it is probably a non-simple polygon */
    if (count-- <= 0)
    {
      /* Triangulate: ERROR - probable bad polygon */
      free(V);
      return false;
    }

    /* three consecutive vertices in current polygon, <u,v,w> */
    u = v  ; u *= (nv > u);     /* previous */
    v = u+1; v *= (nv > v);     /* new v    */
    w = v+1; w *= (nv > w);     /* next     */

    result = polySnip(contour, u, v, w, nv, V);
    if (result > 0) /* successfully found a triangle */
    {
      /* output Triangle */
      if (!(*cb)(contour, V[u], V[v], V[w])) {
      	free(V);
      	return false;
      }
    }
    if (result) /* successfully found a triangle or three consecutive colinear points */
    {
      int s,t;

      /* remove v from remaining polygon */
      for(s=v,t=v+1;t<nv;s++,t++) V[s] = V[t]; nv--;

      /* reset error detection counter */
      count = 2*nv;
    }
  }

  free(V);

  return true;
}
