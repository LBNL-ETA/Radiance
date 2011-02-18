#ifndef lint
static const char	RCSid[] = "$Id: vect3ds.c,v 2.1 2011/02/18 00:40:25 greg Exp $";
#endif
#include <math.h>
#include <string.h>
#include "vect3ds.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif
#define PI	((double)M_PI)

#define EPSILON 1e-6

void   adjoint (Matrix mat);
double det4x4 (Matrix mat);
double det3x3 (double a1, double a2, double a3, double b1, double b2,
	       double b3, double c1, double c2, double c3);
double det2x2 (double a, double b, double c, double d);


void vect_init (Vector v, float  x, float  y, float  z)
{
    v[X] = x;
    v[Y] = y;
    v[Z] = z;
}


void vect_copy (Vector v1, Vector v2)
{
    v1[X] = v2[X];
    v1[Y] = v2[Y];
    v1[Z] = v2[Z];
}


int vect_equal (Vector v1, Vector v2)
{
    if (v1[X] == v2[X] && v1[Y] == v2[Y] && v1[Z] == v2[Z])
	return 1;
    else
	return 0;
}


void vect_add (Vector v1, Vector v2, Vector v3)
{
    v1[X] = v2[X] + v3[X];
    v1[Y] = v2[Y] + v3[Y];
    v1[Z] = v2[Z] + v3[Z];
}


void vect_sub (Vector v1, Vector v2, Vector v3)
{
    v1[X] = v2[X] - v3[X];
    v1[Y] = v2[Y] - v3[Y];
    v1[Z] = v2[Z] - v3[Z];
}


void vect_scale (Vector v1, Vector v2, float  k)
{
    v1[X] = k * v2[X];
    v1[Y] = k * v2[Y];
    v1[Z] = k * v2[Z];
}


float vect_mag (Vector v)
{
    float mag = sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);

    return mag;
}


void vect_normalize (Vector v)
{
    float mag = vect_mag (v);

    if (mag > 0.0)
	vect_scale (v, v, 1.0/mag);
}


float vect_dot (Vector v1, Vector v2)
{
    return (v1[X]*v2[X] + v1[Y]*v2[Y] + v1[Z]*v2[Z]);
}


void vect_cross (Vector v1, Vector v2, Vector v3)
{
    v1[X] = (v2[Y] * v3[Z]) - (v2[Z] * v3[Y]);
    v1[Y] = (v2[Z] * v3[X]) - (v2[X] * v3[Z]);
    v1[Z] = (v2[X] * v3[Y]) - (v2[Y] * v3[X]);
}

void vect_min (Vector v1, Vector v2, Vector v3)
{
    v1[X] = (v2[X] < v3[X]) ? v2[X] : v3[X];
    v1[Y] = (v2[Y] < v3[Y]) ? v2[Y] : v3[Y];
    v1[Z] = (v2[Z] < v3[Z]) ? v2[Z] : v3[Z];
}


void vect_max (Vector v1, Vector v2, Vector v3)
{
    v1[X] = (v2[X] > v3[X]) ? v2[X] : v3[X];
    v1[Y] = (v2[Y] > v3[Y]) ? v2[Y] : v3[Y];
    v1[Z] = (v2[Z] > v3[Z]) ? v2[Z] : v3[Z];
}


/* Return the angle between two vectors */
float vect_angle (Vector v1, Vector v2)
{
    float  mag1, mag2, angle, cos_theta;

    mag1 = vect_mag(v1);
    mag2 = vect_mag(v2);

    if (mag1 * mag2 == 0.0)
	angle = 0.0;
    else {
	cos_theta = vect_dot(v1,v2) / (mag1 * mag2);

	if (cos_theta <= -1.0)
	    angle = 180.0;
	else if (cos_theta >= +1.0)
	    angle = 0.0;
	else
	    angle = (180.0/PI) * acos(cos_theta);
    }

    return angle;
}


void vect_print (FILE *f, Vector v, int dec, char sep)
{
    char fstr[] = "%.4f, %.4f, %.4f";

    if (dec < 0) dec = 0;
    if (dec > 9) dec = 9;

    fstr[2]  = '0' + dec;
    fstr[8]  = '0' + dec;
    fstr[14] = '0' + dec;

    fstr[4]  = sep;
    fstr[10] = sep;

    fprintf (f, fstr, v[X], v[Y], v[Z]);
}


/* Rotate a vector about the X, Y or Z axis */
void vect_rotate (Vector v1, Vector v2, int axis, float angle)
{
    float  cosa, sina;

    cosa = cos ((PI/180.0) * angle);
    sina = sin ((PI/180.0) * angle);

    switch (axis) {
	case X:
	    v1[X] =  v2[X];
	    v1[Y] =  v2[Y] * cosa + v2[Z] * sina;
	    v1[Z] =  v2[Z] * cosa - v2[Y] * sina;
	    break;

	case Y:
	    v1[X] = v2[X] * cosa - v2[Z] * sina;
	    v1[Y] = v2[Y];
	    v1[Z] = v2[Z] * cosa + v2[X] * sina;
	    break;

	case Z:
	    v1[X] = v2[X] * cosa + v2[Y] * sina;
	    v1[Y] = v2[Y] * cosa - v2[X] * sina;
	    v1[Z] = v2[Z];
	    break;
    }
}


/* Rotate a vector about a specific axis */
void vect_axis_rotate (Vector v1, Vector v2, Vector axis, float angle)
{
    float  cosa, sina;
    Matrix mat;

    cosa = cos ((PI/180.0) * angle);
    sina = sin ((PI/180.0) * angle);

    mat[0][0] = (axis[X] * axis[X]) + ((1.0 - (axis[X] * axis[X]))*cosa);
    mat[0][1] = (axis[X] * axis[Y] * (1.0 - cosa)) - (axis[Z] * sina);
    mat[0][2] = (axis[X] * axis[Z] * (1.0 - cosa)) + (axis[Y] * sina);
    mat[0][3] = 0.0;

    mat[1][0] = (axis[X] * axis[Y] * (1.0 - cosa)) + (axis[Z] * sina);
    mat[1][1] = (axis[Y] * axis[Y]) + ((1.0 - (axis[Y] * axis[Y])) * cosa);
    mat[1][2] = (axis[Y] * axis[Z] * (1.0 - cosa)) - (axis[X] * sina);
    mat[1][3] = 0.0;

    mat[2][0] = (axis[X] * axis[Z] * (1.0 - cosa)) - (axis[Y] * sina);
    mat[2][1] = (axis[Y] * axis[Z] * (1.0 - cosa)) + (axis[X] * sina);
    mat[2][2] = (axis[Z] * axis[Z]) + ((1.0 - (axis[Z] * axis[Z])) * cosa);
    mat[2][3] = 0.0;

    mat[3][0] = mat[3][1] = mat[3][2] = mat[3][3] = 0.0;

    vect_transform (v1, v2, mat);
}


/* Transform the given vector */
void vect_transform (Vector v1, Vector v2, Matrix mat)
{
    Vector tmp;

    tmp[X] = (v2[X] * mat[0][0]) + (v2[Y] * mat[1][0]) + (v2[Z] * mat[2][0]) + mat[3][0];
    tmp[Y] = (v2[X] * mat[0][1]) + (v2[Y] * mat[1][1]) + (v2[Z] * mat[2][1]) + mat[3][1];
    tmp[Z] = (v2[X] * mat[0][2]) + (v2[Y] * mat[1][2]) + (v2[Z] * mat[2][2]) + mat[3][2];

    vect_copy (v1, tmp);
}


/* Create an identity matrix */
void mat_identity (Matrix mat)
{
    int i, j;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    mat[i][j] = 0.0;

    for (i = 0; i < 4; i++)
	mat[i][i] = 1.0;
}


void mat_copy (Matrix mat1, Matrix mat2)
{
    int i, j;

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    mat1[i][j] = mat2[i][j];
}


/* Rotate a matrix about the X, Y or Z axis */
void mat_rotate (Matrix mat1, Matrix mat2, int axis, float angle)
{
    Matrix mat;
    float  cosa, sina;

    cosa = cos ((PI/180.0) * angle);
    sina = sin ((PI/180.0) * angle);

    mat_identity (mat);

    switch (axis) {
	case X:
	    mat[1][1] = cosa;
	    mat[1][2] = sina;
	    mat[2][1] = -sina;
	    mat[2][2] = cosa;
	    break;

	case Y:
	    mat[0][0] = cosa;
	    mat[0][2] = -sina;
	    mat[2][0] = sina;
	    mat[2][2] = cosa;
	    break;

	case Z:
	    mat[0][0] = cosa;
	    mat[0][1] = sina;
	    mat[1][0] = -sina;
	    mat[1][1] = cosa;
	    break;
    }

    mat_mult (mat1, mat2, mat);
}


void mat_axis_rotate (Matrix mat1, Matrix mat2, Vector axis, float angle)
{
    float  cosa, sina;
    Matrix mat;

    cosa = cos ((PI/180.0) * angle);
    sina = sin ((PI/180.0) * angle);

    mat[0][0] = (axis[X] * axis[X]) + ((1.0 - (axis[X] * axis[X]))*cosa);
    mat[0][1] = (axis[X] * axis[Y] * (1.0 - cosa)) - (axis[Z] * sina);
    mat[0][2] = (axis[X] * axis[Z] * (1.0 - cosa)) + (axis[Y] * sina);
    mat[0][3] = 0.0;

    mat[1][0] = (axis[X] * axis[Y] * (1.0 - cosa)) + (axis[Z] * sina);
    mat[1][1] = (axis[Y] * axis[Y]) + ((1.0 - (axis[Y] * axis[Y])) * cosa);
    mat[1][2] = (axis[Y] * axis[Z] * (1.0 - cosa)) - (axis[X] * sina);
    mat[1][3] = 0.0;

    mat[2][0] = (axis[X] * axis[Z] * (1.0 - cosa)) - (axis[Y] * sina);
    mat[2][1] = (axis[Y] * axis[Z] * (1.0 - cosa)) + (axis[X] * sina);
    mat[2][2] = (axis[Z] * axis[Z]) + ((1.0 - (axis[Z] * axis[Z])) * cosa);
    mat[2][3] = 0.0;

    mat[3][0] = mat[3][1] = mat[3][2] = mat[3][3] = 0.0;

    mat_mult (mat1, mat2, mat);
}


/*  mat1 <-- mat2 * mat3 */
void mat_mult (Matrix mat1, Matrix mat2, Matrix mat3)
{
    float sum;
    int   i, j, k;
    Matrix result;

    for (i = 0; i < 4; i++) {
	for (j = 0; j < 4; j++) {
	    sum = 0.0;

	    for (k = 0; k < 4; k++)
		sum = sum + mat2[i][k] * mat3[k][j];

	    result[i][j] = sum;
	}
    }

    for (i = 0; i < 4; i++)
	for (j = 0; j < 4; j++)
	    mat1[i][j] = result[i][j];
}


/*
   Decodes a 3x4 transformation matrix into separate scale, rotation,
   translation, and shear vectors. Based on a program by Spencer W.
   Thomas (Graphics Gems II)
*/
void mat_decode (Matrix mat, Vector scale,  Vector shear, Vector rotate,
		Vector transl)
{
    int i;
    Vector row[3], temp;

    for (i = 0; i < 3; i++)
	transl[i] = mat[3][i];

    for (i = 0; i < 3; i++) {
	row[i][X] = mat[i][0];
	row[i][Y] = mat[i][1];
	row[i][Z] = mat[i][2];
    }

    scale[X] = vect_mag (row[0]);
    vect_normalize (row[0]);

    shear[X] = vect_dot (row[0], row[1]);
    row[1][X] = row[1][X] - shear[X]*row[0][X];
    row[1][Y] = row[1][Y] - shear[X]*row[0][Y];
    row[1][Z] = row[1][Z] - shear[X]*row[0][Z];

    scale[Y] = vect_mag (row[1]);
    vect_normalize (row[1]);

    if (scale[Y] != 0.0)
	shear[X] /= scale[Y];

    shear[Y] = vect_dot (row[0], row[2]);
    row[2][X] = row[2][X] - shear[Y]*row[0][X];
    row[2][Y] = row[2][Y] - shear[Y]*row[0][Y];
    row[2][Z] = row[2][Z] - shear[Y]*row[0][Z];

    shear[Z] = vect_dot (row[1], row[2]);
    row[2][X] = row[2][X] - shear[Z]*row[1][X];
    row[2][Y] = row[2][Y] - shear[Z]*row[1][Y];
    row[2][Z] = row[2][Z] - shear[Z]*row[1][Z];

    scale[Z] = vect_mag (row[2]);
    vect_normalize (row[2]);

    if (scale[Z] != 0.0) {
	shear[Y] /= scale[Z];
	shear[Z] /= scale[Z];
    }

    vect_cross (temp, row[1], row[2]);
    if (vect_dot (row[0], temp) < 0.0) {
	for (i = 0; i < 3; i++) {
	    scale[i]  *= -1.0;
	    row[i][X] *= -1.0;
	    row[i][Y] *= -1.0;
	    row[i][Z] *= -1.0;
	}
    }

    if (row[0][Z] < -1.0) row[0][Z] = -1.0;
    if (row[0][Z] > +1.0) row[0][Z] = +1.0;

    rotate[Y] = asin(-row[0][Z]);

    if (fabs(cos(rotate[Y])) > EPSILON) {
	rotate[X] = atan2 (row[1][Z], row[2][Z]);
	rotate[Z] = atan2 (row[0][Y], row[0][X]);
    }
    else {
	rotate[X] = atan2 (row[1][X], row[1][Y]);
	rotate[Z] = 0.0;
    }

    /* Convert the rotations to degrees */
    rotate[X] = (180.0/PI)*rotate[X];
    rotate[Y] = (180.0/PI)*rotate[Y];
    rotate[Z] = (180.0/PI)*rotate[Z];
}


/* Matrix inversion code from Graphics Gems */

/* mat1 <-- mat2^-1 */
float mat_inv (Matrix mat1, Matrix mat2)
{
    int i, j;
    float det;

    if (mat1 != mat2) {
	for (i = 0; i < 4; i++)
	    for (j = 0; j < 4; j++)
		mat1[i][j] = mat2[i][j];
    }

    det = det4x4 (mat1);

    if (fabs (det) < EPSILON)
	return 0.0;

    adjoint (mat1);

    for (i = 0; i < 4; i++)
	for(j = 0; j < 4; j++)
	    mat1[i][j] = mat1[i][j] / det;

    return det;
}


void adjoint (Matrix mat)
{
    double a1, a2, a3, a4, b1, b2, b3, b4;
    double c1, c2, c3, c4, d1, d2, d3, d4;

    a1 = mat[0][0]; b1 = mat[0][1];
    c1 = mat[0][2]; d1 = mat[0][3];

    a2 = mat[1][0]; b2 = mat[1][1];
    c2 = mat[1][2]; d2 = mat[1][3];

    a3 = mat[2][0]; b3 = mat[2][1];
    c3 = mat[2][2]; d3 = mat[2][3];

    a4 = mat[3][0]; b4 = mat[3][1];
    c4 = mat[3][2]; d4 = mat[3][3];

    /* row column labeling reversed since we transpose rows & columns */
    mat[0][0]  =  det3x3 (b2, b3, b4, c2, c3, c4, d2, d3, d4);
    mat[1][0]  = -det3x3 (a2, a3, a4, c2, c3, c4, d2, d3, d4);
    mat[2][0]  =  det3x3 (a2, a3, a4, b2, b3, b4, d2, d3, d4);
    mat[3][0]  = -det3x3 (a2, a3, a4, b2, b3, b4, c2, c3, c4);

    mat[0][1]  = -det3x3 (b1, b3, b4, c1, c3, c4, d1, d3, d4);
    mat[1][1]  =  det3x3 (a1, a3, a4, c1, c3, c4, d1, d3, d4);
    mat[2][1]  = -det3x3 (a1, a3, a4, b1, b3, b4, d1, d3, d4);
    mat[3][1]  =  det3x3 (a1, a3, a4, b1, b3, b4, c1, c3, c4);

    mat[0][2]  =  det3x3 (b1, b2, b4, c1, c2, c4, d1, d2, d4);
    mat[1][2]  = -det3x3 (a1, a2, a4, c1, c2, c4, d1, d2, d4);
    mat[2][2]  =  det3x3 (a1, a2, a4, b1, b2, b4, d1, d2, d4);
    mat[3][2]  = -det3x3 (a1, a2, a4, b1, b2, b4, c1, c2, c4);

    mat[0][3]  = -det3x3 (b1, b2, b3, c1, c2, c3, d1, d2, d3);
    mat[1][3]  =  det3x3 (a1, a2, a3, c1, c2, c3, d1, d2, d3);
    mat[2][3]  = -det3x3 (a1, a2, a3, b1, b2, b3, d1, d2, d3);
    mat[3][3]  =  det3x3 (a1, a2, a3, b1, b2, b3, c1, c2, c3);
}


double det4x4 (Matrix mat)
{
    double ans;
    double a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, d1, d2, d3, 			d4;

    a1 = mat[0][0]; b1 = mat[0][1];
    c1 = mat[0][2]; d1 = mat[0][3];

    a2 = mat[1][0]; b2 = mat[1][1];
    c2 = mat[1][2]; d2 = mat[1][3];

    a3 = mat[2][0]; b3 = mat[2][1];
    c3 = mat[2][2]; d3 = mat[2][3];

    a4 = mat[3][0]; b4 = mat[3][1];
    c4 = mat[3][2]; d4 = mat[3][3];

    ans = a1 * det3x3 (b2, b3, b4, c2, c3, c4, d2, d3, d4) -
	  b1 * det3x3 (a2, a3, a4, c2, c3, c4, d2, d3, d4) +
	  c1 * det3x3 (a2, a3, a4, b2, b3, b4, d2, d3, d4) -
	  d1 * det3x3 (a2, a3, a4, b2, b3, b4, c2, c3, c4);

    return ans;
}


double det3x3 (double a1, double a2, double a3, double b1, double b2,
	       double b3, double c1, double c2, double c3)
{
    double ans;

    ans = a1 * det2x2 (b2, b3, c2, c3)
	- b1 * det2x2 (a2, a3, c2, c3)
	+ c1 * det2x2 (a2, a3, b2, b3);

    return ans;
}


double det2x2 (double a, double b, double c, double d)
{
    double ans;
    ans = a * d - b * c;
    return ans;
}

