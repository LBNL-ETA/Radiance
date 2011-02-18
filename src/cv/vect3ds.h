/* RCSid: $Id: vect3ds.h,v 2.1 2011/02/18 00:40:25 greg Exp $ */
#ifndef _MGF_VECT_H_
#define _MGF_VECT_H_
#ifdef __cplusplus
extern "C" {
#endif


#include <stdio.h>

typedef float Vector[3];

#define X 0
#define Y 1
#define Z 2

/* Transformation matrix */
typedef float Matrix[4][4];

void vect_init (Vector v, float x, float y, float z);
void vect_copy (Vector v1, Vector v2);
int vect_equal (Vector v1, Vector v2);
void vect_add (Vector v1, Vector v2, Vector v3);
void vect_sub (Vector v1, Vector v2, Vector v3);
void vect_scale (Vector v1, Vector v2, float k);
float vect_mag (Vector v);
void vect_normalize (Vector v);
float vect_dot (Vector v1, Vector v2);
void vect_cross (Vector v1, Vector v2, Vector v3);
void vect_min (Vector v1, Vector v2, Vector v3);
void vect_max (Vector v1, Vector v2, Vector v3);
float vect_angle (Vector v1, Vector v2);
void vect_print (FILE *f, Vector v, int dec, char sep);
void vect_rotate (Vector v1, Vector v2, int axis, float angle);
void vect_axis_rotate (Vector v1, Vector v2, Vector axis, float angle);
void mat_identity (Matrix mat);
void mat_copy (Matrix mat1, Matrix mat2);
void mat_rotate (Matrix mat1, Matrix mat2, int axis, float angle);
void mat_axis_rotate (Matrix mat1, Matrix mat2, Vector axis, float angle);
void mat_mult (Matrix mat1, Matrix mat2, Matrix mat3);
void vect_transform (Vector v1, Vector v2, Matrix mat);
void mat_decode (Matrix mat, Vector scale, Vector shear, Vector rotate, Vector
	 transl);
float mat_inv (Matrix mat1, Matrix mat2);


#ifdef __cplusplus
}
#endif
#endif /* _MGF_VECT_H_ */
