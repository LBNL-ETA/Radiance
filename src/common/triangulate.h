/* RCSid $Id: triangulate.h,v 2.3 2014/02/07 18:58:40 greg Exp $ */
/*
 *  triangulate.h
 *  
 *  Adapted by Greg Ward on 1/23/14.
 *  Copyright 2014 Anyhere Software. All rights reserved.
 *
 */

/* COTD Entry submitted by John W. Ratcliff [jratcliff@verant.com] */

/**** THIS IS A CODE SNIPPET WHICH WILL EFFICIEINTLY TRIANGULATE ANY
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

/**********************************************************************/
/************ HEADER FILE FOR TRIANGULATE.C ***************************/
/**********************************************************************/

#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	double	mX, mY;
} Vert2;

typedef struct {
	int	nv;		/* number of vertices */
	void	*p;		/* client data pointer */
	Vert2	v[3];		/* extends struct */
} Vert2_list;

/* allocate a polygon list */
extern Vert2_list *polyAlloc(int nv);

#define polyFree	free

/* callback for output triangle */
typedef int tri_out_t(const Vert2_list *tp, int a, int b, int c);

/* triangulate a polygon, send results to the given callback function */
extern int polyTriangulate(const Vert2_list *contour, tri_out_t *cb);

/* compute area of a contour/polygon */
extern double polyArea(const Vert2_list *contour);

/* decide if point P is inside triangle defined by ABC */
extern int insideTriangle(double Ax, double Ay,
                      double Bx, double By,
                      double Cx, double Cy,
                      double Px, double Py);

#ifdef __cplusplus
}
#endif

#endif
