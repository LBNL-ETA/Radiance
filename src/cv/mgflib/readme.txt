		MGF PACKAGE DESCRIPTION
		RCSid "$Id: readme.txt,v 1.8 2003/02/28 20:19:26 greg Exp $"

This package includes a description and parser for a new scene
description standard, called for the lack of a better name, MGF
for Materials and Geometry Format.  It was developed by Greg
Ward of the Lawrence Berkeley Laboratory <GJWard@lbl.gov> with
help and advice from Rob Shakespeare of Indiana University
<TCVC@ucs.indiana.edu>, Ian Ashdown of Ledalite Corporation
<72060.2420@CompuServe.COM> and Holly Rushmeier of the National
Institute for Standards and Technology <holly@cam.nist.gov>.

The language itself is described in the file "spec.txt", and
the included Makefile should make building the parser library
fairly straightforward on most systems.  What's left then, is
explaining the why and how of using this package.

The initial purpose of developing a scene description standard
was for inclusion in the Illumination Engineering Society's (IES)
standard data representation for luminaires.  It occurred to us
early on that such a standard might have broader applications,
so an effort was made to create a fairly general description
language, while keeping it as simple as possible for the people
who have to create descriptions with it as well as the programmers
who have to support it.

Why create a new standard rather than exploiting an existing one?
Some of the rationale for our decision is explained at the end of
the specification document, but it mostly boils down to materials.
As easy as it is to describe physically valid materials, most
scene description languages cannot do it.  The material specification
included in the MGF standard may not be perfect, but at least it's
physically plausible.  Furthermore, we are committed to making any
future modifications to the standard backwards-compatible -- a rather
tricky proposition.

This takes us to the how of supporting this new standard.  The basic
approach is to use the standard parser, which does a lot of the work
in supporting the language itself.  The programmer tells the parser
which entities it will support, and the parser does the rest.
That way, it isn't necessary to modify the program when a new version
of the standard comes out; all one has to do is link to the new
standard's parser.  (The include file will change as well, so it's
not QUITE that simple, but close.)

There are two ways to support the language, by linking the parser to
the program itself, or by linking the parser to a translator program
that expresses MGF entities in the native scene description format.
The differences in the two approaches are slight, and we will mention
them following a general explanation of the parser and support library.

The Parser
==========
The MGF parser is written in ANSI-C (though the -DNOPROTO flag may be
used to get back K&R compatibility).  All of the declarations and
definitions needed are in the single include file "parser.h".  This
file is a good place to look for details on using the various support
routines as well.  The parser itself is parser.c, though it relies for
some translations on other C modules.  These same support routines will
no doubt be useful for applications programmers, and we will explain
some of them in the following sections.

Initializing the parser is the most important part of writing an MGF
program, and it is done through the mg_ehand array and a call to mg_init.
The global mg_ehand variable is an array of pointers to entity handler
functions.  The arguments to these functions are always the same, an
argument count and an array of argument pointers (ala main).  The return
value for these integer functions is one of the error codes defined in
parser.h, or MG_OK if the entity was handled correctly.  You must
set the appropriate entries for the entities you can support, then call
mg_init to fill in the rest.  Most of the entities you cannot support
will be translated into (approximately) equivalent ones you can.
Entities that have no equivalent (such as color), will be safely
ignored on the input.  If you have specified support for some entities
without offering support to their prerequisites, mg_init will report an
error and exit.

Once the parser has been properly initialized, MGF input files may be
loaded at will with the mg_load call.  This function takes a single
argument, which is the name of the MGF file.  (The NULL pointer may be
used to specify standard input.)  The behavior of the parser in part
depends on input history, so the mg_clear call should be used after
each file if starting fresh is important.  This also frees any data
structures used by the parser, which may be desirable if the program
is going to do something after loading besides exit.

Support Functions
=================
In translating unsupported entities, the parser makes use of a number
of support functions, contained in associated C modules.  The most
important of these are in context.c, which includes three handler
functions that can support all color, material and vertex entities.
To understand what these functions do, it is necessary to know a
little about the MGF language itself, so please familiarize yourself
with it now if you haven't already.  (See the file "spec.txt".)

Context Support
===============
The MGF language defines three named contexts, the current vertex,
the current color and the current material.  (The current color is
used mostly for setting parameters in the current material.)  There
are three handler routines defined in context.c, and they can handle
all entities related to these three contexts.  The simplest way to
support materials, for example, is to initialize the mg_ehand array
such that the MG_E_MATERIAL, MG_E_RD, MG_E_RS, etc. entries all point
to c_hmaterial.  Then, whenever a material is needed, the global
c_cmaterial variable will be pointing to a structure with all the
current settings.  (Note that you would have to also set the color
mg_ehand entries to c_hcolor if you intended to support color
materials.)  A list of related mg_ehand assignments is given below:

	mg_ehand[MG_E_COLOR] = c_hcolor;
	mg_ehand[MG_E_CCT] = c_hcolor;
	mg_ehand[MG_E_CMIX] = c_hcolor;
	mg_ehand[MG_E_CSPEC] = c_hcolor;
	mg_ehand[MG_E_CXY] = c_hcolor;
	mg_ehand[MG_E_ED] = c_hmaterial;
	mg_ehand[MG_E_IR] = c_hmaterial;
	mg_ehand[MG_E_MATERIAL] = c_hmaterial;
	mg_ehand[MG_E_NORMAL] = c_hvertex;
	mg_ehand[MG_E_POINT] = c_hvertex;
	mg_ehand[MG_E_RD] = c_hmaterial;
	mg_ehand[MG_E_RS] = c_hmaterial;
	mg_ehand[MG_E_SIDES] = c_hmaterial;
	mg_ehand[MG_E_TD] = c_hmaterial;
	mg_ehand[MG_E_TS] = c_hmaterial;
	mg_ehand[MG_E_VERTEX] = c_hvertex;

In addition to the three handler functions, context.c contains a
few support routines that make life simpler.  For vertices, there
is the c_getvertex call, which returns a pointer to a named vertex
structure (or NULL if there is no corresponding definition for the
given name).  This function is needed for support of most surface
entities.  For color support, there is the analogous c_getcolor call,
and the c_ccvt routine, which is used to convert from one color
representation to another (e.g.  spectral color to xy chromaticity
coordinates).  Also, there is a function called c_isgrey, which
simply returns 1 or 0 based on whether the passed color structure
is close to grey or not.  Finally, there is the c_clearall routine,
which clears and frees all context data structures, and is the
principal action of the parser's mg_clear function.

Transform Support
=================
If your program is supporting any geometry at all (and what would be
the point if it wasn't?) you will need to support the transform
entity (MG_E_XF).  This would be tricky, if it weren't for the support
routines provided, which make the task fairly painless.  First, there
is the transform handler itself, xf_handler.  Just set the MG_E_XF
entry of the mg_ehand array to this function.  Then, anytime you want
to transform something, call one of the associated functions, xf_xfmpoint,
xf_xfmvect, xf_rotvect or xf_scale.  These functions transform a 3-D
point, 3-D vector (without translation), rotate a 3-D vector (without
scaling) and scale a floating-point value, respectively.

Object Support
==============
The MGF language includes a single entity for naming objects, MG_E_OBJECT.
It is mostly provided as a convenience for the user, so that individual
geometric parts may be easily identified.  Although supporting this entity
directly is possible, it's hierarchical nature requires maintaining a stack
of object names.  The object handler in object.c provides this functionality.
Simply set the MG_E_OBJECT entry of the mg_ehand array to obj_handler,
and the current object name list will be kept in the global array obj_name.
The number of names is stored in the global obj_nnames variable.  To clear
this array (freeing any memory used in the process), call obj_clear.

Loading vs. Translating
=======================
As mentioned in the introduction, the parser may be used either to load
data into a rendering program directly, or to get MGF input for translation
to another file format.  In either case, the procedure is nearly identical.
The only important difference is what you do with the parser data structures
after loading.  For a translator, this is not an issue, but rendering
programs usually need all the memory they can get.  Therefore, once the
input process is complete, you should call the mg_clear function to free
the parser data structures and return to an initialized state (i.e. it
is never necessary to recall the mg_init routine).

Also, if you use some of the support functions, you should call their
specific clearing functions.  For the transform module, the call is
xf_clear.  For the object support module, the call is obj_clear.  The
context routines use the c_clearall function, but this is actually
called by mg_clear, so calling it again is unnecessary.

Linking Vertices
================
Although the MGF language was designed with linking vertices in mind,
there are certain aspects which make this goal more challenging.
Specifically, the ability to redefine values for a previously named
vertex is troublesome for the programmer, since the same vertex can
have different values at different points in the input.  Likewise, the
effect of the transform entity on surfaces rather than vertices means
that the same named vertex can appear in many positions.

It is not possible to use the parser data structures directly for
linking vertices, but we've taken a couple of steps in the support
routines to make the task of organizing your own data structures a
little easier.  First, there is a clock member in the C_VERTEX
structure that is incremented on each change.  (The same member is
contained in the C_COLOR and C_MATERIAL structures.)  Second, the
current transform (pointed to by xf_context) contains a unique
identifier, xf_context->xid.  This is a long integer that will be
different for each unique transform.  (It is actually a hash key on the
transformation matrix, and there is about 1 chance in 2 billion that
two different matrices will hash to the same value.  Is this a bug?
I guess it depends on how long the programmer lives -- or vice versa.)

There are two ways to use of this additional information.  One
is to record the vertex clock value along with it's id and the
current xf_context->xid value.  If another vertex comes along with
the same name, but one of these two additional values fails to match,
then it (probably) is a different vertex.  Alternatively, one can reset
the clock member every time a new vertex is stored.  That way, it is
only necessary to check the clock against zero rather than storing this
value along with the vertex name and transform id.  If the name and
transform are the same and the clock is zero, then it's the same vertex
as last time.

Yet another approach is to ignore the parser structures entirely and
focus on the actual vertex values.  After all, the user is not compelled
to reuse the same vertex names for the same points.  It is just as likely
that the same vertices will appear under different names, so that none
of the above would help to merge them.  The most sure-fire approach to
linking identical vertices is therefore to hash the point and normal
values directly and use the functions in lookup.c to associate them.
You will have to write your own hash function, and we recommend making
one that allows a little slop so that nearly identical points hash to
the same value.

Examples
========
Two example translator programs are included with this package.

The simplest is a translator from MGF to MGF called mgfilt.c, which
produces on the standard output only those entities from the standard
input that are supported according to the first command line argument.
For example, one could remove everything but the raw, flat polygonal
geometry with the following command:

	mgfilt v,p,f,xf any.mgf > faces.mgf

Note that the xf entity must also be included, for its support is
required by all geometric entities.

The second translator converts from MGF to the Radiance scene description
language, and is a more practical example of parser use.  Unfortunately,
we did not include all of the support functions required by this translator,
so it serves as a source code example only.  If you wish to get the rest
of it because you intend to run it, contact Greg Ward <GJWard@lbl.gov>
and he'll be happy to provide you with the missing pieces.

Copyright
=========
At this point, the legal issues related to this parser have not been
worked out.  The intent is to offer it free of charge to all those who
wish to use it (with no guarantees, of course).  However, we may decide
that copyright protections are necessary to prevent unauthorized versions
of the parser, which do not properly support the MGF standard, from
getting spread around.  Since this is a pre-release, we trust that you
will not share it with anyone without getting our permission first.

Questions
=========
Questions should be directed to Greg Ward <GJWard@lbl.gov>, who will be
happy to offer any reasonable assistance in using this standard.  (Greg's
telephone is 1-510-486-4757, fax 1-510-486-4089.)
