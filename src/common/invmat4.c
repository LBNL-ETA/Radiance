#ifndef lint
static const char	RCSid[] = "$Id: invmat4.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 * invmat4 - computes the inverse of mat into inverse.  Returns 1
 * if there exists an inverse, 0 otherwise.  It uses Gaussian Elimination
 * method.
 */

#include "copyright.h"

#include "mat4.h"


#define  ABS(x)		((x)>=0 ? (x) : -(x))

#define SWAP(a,b,t) (t=a,a=b,b=t)


int
invmat4(inverse, mat)
MAT4  inverse, mat;
{
	MAT4  m4tmp;
	register int i,j,k;
	register double temp;

	copymat4(m4tmp, mat);
					/* set inverse to identity */
	setident4(inverse);

	for(i = 0; i < 4; i++) {
		/* Look for row with largest pivot and swap rows */
		temp = FTINY; j = -1;
		for(k = i; k < 4; k++)
			if(ABS(m4tmp[k][i]) > temp) {
				temp = ABS(m4tmp[k][i]);
				j = k;
				}
		if(j == -1)	/* No replacing row -> no inverse */
			return(0);
		if (j != i)
			for(k = 0; k < 4; k++) {
				SWAP(m4tmp[i][k],m4tmp[j][k],temp);
				SWAP(inverse[i][k],inverse[j][k],temp);
				}

		temp = m4tmp[i][i];
		for(k = 0; k < 4; k++) {
			m4tmp[i][k] /= temp;
			inverse[i][k] /= temp;
			}
		for(j = 0; j < 4; j++) {
			if(j != i) {
				temp = m4tmp[j][i];
				for(k = 0; k < 4; k++) {
					m4tmp[j][k] -= m4tmp[i][k]*temp;
					inverse[j][k] -= inverse[i][k]*temp;
					}
				}
			}
		}
	return(1);
}
