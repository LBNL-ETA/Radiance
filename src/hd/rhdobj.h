/* RCSid: $Id: rhdobj.h,v 3.7 2003/07/14 22:24:00 schorsch Exp $ */
/*
 * Header file for object display routines for rholo drivers.
 */
#ifndef _RAD_RHDOBJ_H_
#define _RAD_RHDOBJ_H_

#ifdef __cplusplus
extern "C" {
#endif

				/* additional user commands */
#define DO_LOAD		0		/* load octree object */
#define DO_UNLOAD	1		/* unload (free) octree object */
#define DO_XFORM	2		/* set/print object transform */
#define DO_MOVE		3		/* add to object transform */
#define DO_UNMOVE	4		/* undo last transform change */
#define DO_DUP		5		/* duplicate object */
#define DO_SHOW		6		/* show object (low quality) */
#define DO_LIGHT	7		/* illuminate object (high quality) */
#define DO_HIDE		8		/* hide object */
#define DO_OBJECT	9		/* print object statistics */

#define DO_NCMDS	10		/* number of object display commands */

					/* commands entered on stdin only */
#define DO_INIT		{"load","clear","xform","move","unmove","dup",\
				"show","light","hide","object"}

/*******************************************************************
 * The following routines are available for calling from the driver:

int
(*dobj_lightsamp)(c, p, v)	: apply next sample for local light sources
COLR	c;			: pixel color (RGBE)
FVECT	p;			: world intersection point
FVECT	v;			: ray direction vector

If the pointer to the function dobj_lightsamp is non-NULL, then the
dev_value() driver routine should call this rather than use a
given sample in its own data structures.  This pointer is set
in the dobj_lighting() function described below, which may be
called indirectly by the dobj_command() function.


int
dobj_command(cmd, args)		: check/run object display command
char	*cmd, *args;		: command name and argument string

Check to see if this is an object display command, and return -1 if
it isn't.  If it is a valid command that results in some visible
change, return non-zero, otherwise return 0.  Error messages should
be printed with error(COMMAND,err).


double
dobj_trace(nm, rorg, rdir)	: check for ray intersection with objects
char	nm[];			: object name (modified)
FVECT	rorg, rdir;		: ray origin and direction

Check to see if the given ray intersects the given object,
returning the distance the ray traveled if it did, or FHUGE if it didn't.
If nm contains "*", then all visible objects are checked and the name
of the intersected object is returned, or "" if none.  If nm is NULL,
then all visible objects are checked, but the name is not returned.


int
dobj_render()			: render visible objects to OpenGL

Renders all display objects using OpenGL, assuming desired view has
been set.  This routine also assumes that the tone-mapping library
is being used to set exposure, and it queries this information to
adjust light sources as necessary for illuminated objects.  Returns
the number of objects rendered.


void
dobj_cleanup()			: clean up data structures

Frees all allocated data structures and objects.


++++ The following routines are usually called through dobj_command() ++++


int
dobj_load(oct, nam)		: create/load an octree object
char	*oct, *nam;		: octree and object name

Load an octree for rendering as a local object, giving it the indicated name.
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_unload(nam)		: free the named object
char	*nam;			: object name

Free all memory associated with the named object.
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_xform(nam, add, ac, av)	: set/add transform for nam
char	*nam;			: object name
int	add;			: add to transform or create new one?
int	ac;			: transform argument count
char	*av[];			: transform arguments

Set or add to the transform for the named object.
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_unmove()			: undo last transform change

Undo the last transform change.
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_dup(oldnm, nam)		: duplicate object oldnm as nam
char	*oldnm, *nam;

Duplicate the named object and give the copy a new name.
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_lighting(nam, cn)		: set up lighting for display object
char	*nam;			: object name
int	cn;			: new drawing code

Set the lighting of the named object to DO_SHOW, DO_HIDE or DO_LIGHT.
Change lighting of all objects if nam is "*".
Returns 1 on success, 0 on failure (with COMMAND error message).


int
dobj_putstats(nam, fp)		: print object statistics
char	*nam;			: object name
FILE	*fp;			: output file

Print current position, size and transform for the named object,
or all objects if nam is "*".
Returns 1 on success, 0 on failure (with COMMAND error message).


 ******************************************************************/


extern double	dobj_trace();

extern char	rhdcmd[DO_NCMDS][8];

extern int	(*dobj_lightsamp)();	/* pointer to function to get lights */

#ifdef __cplusplus
}
#endif
#endif /* _RAD_RHDOBJ_H_ */

