#ifndef lint
static const char	RCSid[] = "$Id: initotypes.c,v 2.4 2004/03/27 12:41:45 schorsch Exp $";
#endif
/*
 * Initialize ofun[] list for octree generator
 */

#include  "standard.h"
#include  "octree.h"
#include  "otypes.h"
#include  "oconv.h"

extern int  o_sphere(); /* XXX way too much linker magic involved here */ 
extern int  o_face();
extern int  o_cone();
extern int  o_instance();
extern int  o_mesh();

FUN  ofun[NUMOTYPE] = INIT_OTYPE;

void
ot_initotypes(void)			/* initialize ofun array */
{
	ofun[OBJ_SPHERE].funp =
	ofun[OBJ_BUBBLE].funp = o_sphere;
	ofun[OBJ_FACE].funp = o_face;
	ofun[OBJ_CONE].funp =
	ofun[OBJ_CUP].funp =
	ofun[OBJ_CYLINDER].funp =
	ofun[OBJ_TUBE].funp =
	ofun[OBJ_RING].funp = o_cone;
	ofun[OBJ_INSTANCE].funp = o_instance;
	ofun[OBJ_MESH].funp = o_mesh;
}


int
o_default(void)			/* default action is no intersection */
{
	return(O_MISS);
}
