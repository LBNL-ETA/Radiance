#ifndef lint
static const char	RCSid[] = "$Id: mat4.c,v 2.4 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 *  mat4.c - routines dealing with 4 X 4 homogeneous transformation matrices.
 */

#include "copyright.h"

#include  "mat4.h"

MAT4  m4ident = MAT4IDENT;

static MAT4  m4tmp;		/* for efficiency */


void
multmat4(m4a, m4b, m4c)		/* multiply m4b X m4c and put into m4a */
MAT4  m4a;
register MAT4  m4b, m4c;
{
	register int  i, j;
	
	for (i = 4; i--; )
		for (j = 4; j--; )
			m4tmp[i][j] = m4b[i][0]*m4c[0][j] +
				      m4b[i][1]*m4c[1][j] +
				      m4b[i][2]*m4c[2][j] +
				      m4b[i][3]*m4c[3][j];
	
	copymat4(m4a, m4tmp);
}


void
multv3(v3a, v3b, m4)	/* transform vector v3b by m4 and put into v3a */
FVECT  v3a;
register FVECT  v3b;
register MAT4  m4;
{
	m4tmp[0][0] = v3b[0]*m4[0][0] + v3b[1]*m4[1][0] + v3b[2]*m4[2][0];
	m4tmp[0][1] = v3b[0]*m4[0][1] + v3b[1]*m4[1][1] + v3b[2]*m4[2][1];
	m4tmp[0][2] = v3b[0]*m4[0][2] + v3b[1]*m4[1][2] + v3b[2]*m4[2][2];
	
	v3a[0] = m4tmp[0][0];
	v3a[1] = m4tmp[0][1];
	v3a[2] = m4tmp[0][2];
}


void
multp3(p3a, p3b, m4)		/* transform p3b by m4 and put into p3a */
register FVECT  p3a;
FVECT  p3b;
register MAT4  m4;
{
	multv3(p3a, p3b, m4);	/* transform as vector */
	p3a[0] += m4[3][0];	/* translate */
	p3a[1] += m4[3][1];
	p3a[2] += m4[3][2];
}
