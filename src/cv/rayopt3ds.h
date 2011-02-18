/* RCSid: $Id: rayopt3ds.h,v 2.1 2011/02/18 00:40:25 greg Exp $ */
/*-------------------------------------------------------------------------

		 Triangle Bounder/Smoother for POV-Ray
		    by Steve Anger 1993

    A number of C routines that can be used to generate POV-Ray ray tracer
 files from triangle data.  Supports generation of smooth triangles and an
 optimal set of bounding shapes for much faster traces.

--------------------------------------------------------------------------*/

#ifndef _MGF_RAYOPT_H
#define _MGF_RAYOPT_H
#ifdef __cplusplus
extern "C" {
#endif


#include "vect3ds.h"

void opt_set_format (int format);
void opt_set_fname (char *pov_name, char *inc_name);
void opt_set_quiet (int quiet);
void opt_set_bound (int bound);
void opt_set_smooth (float smooth);
void opt_set_vert (unsigned vert);
void opt_set_dec (int dec);

void opt_set_color (float red, float green, float blue);
void opt_set_texture (char *texture_name);
void opt_set_transform (Matrix mat);
void opt_clear_transform();
int  opt_add_tri (float ax, float ay, float az,
		  float bx, float by, float bz,
		  float cx, float cy, float cz);

void opt_write_pov (char *obj_name);
void opt_write_file (char *obj_name);
void opt_write_box (char *obj_name);
void opt_finish (void);

void opt_get_limits (float *min_x, float *min_y, float *min_z,
		     float *max_x, float *max_y, float *max_z);
void opt_get_glimits (float *min_x, float *min_y, float *min_z,
		      float *max_x, float *max_y, float *max_z);
unsigned opt_get_vert_cnt (void);
unsigned opt_get_tri_cnt (void);
float    opt_get_index (void);
unsigned opt_get_bounds (void);

void abortmsg (char *msg, int exit_code);
void add_ext (char *fname, char *ext, int force);
void cleanup_name (char *name);


#ifdef __cplusplus
}
#endif
#endif /* _MGF_RAYOPT_H */

