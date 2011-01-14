/* RCSid: $Id: parser.h,v 1.39 2011/01/14 05:46:12 greg Exp $ */
/*
 * Header file for MGF interpreter
 */
#ifndef _MGF_PARSER_H_
#define _MGF_PARSER_H_

/* must include stdio.h and stdlib.h before us */

#include "ccolor.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MG_VMAJOR	2		/* major version number */
#define MG_VMINOR	0		/* minor version number */

			/* Entities (list is only appended, never modified) */
#define MG_E_COMMENT	0		/* #		*/
#define MG_E_COLOR	1		/* c		*/
#define MG_E_CCT	2		/* cct		*/
#define MG_E_CONE	3		/* cone		*/
#define MG_E_CMIX	4		/* cmix		*/
#define MG_E_CSPEC	5		/* cspec	*/
#define MG_E_CXY	6		/* cxy		*/
#define MG_E_CYL	7		/* cyl		*/
#define MG_E_ED		8		/* ed		*/
#define MG_E_FACE	9		/* f		*/
#define MG_E_INCLUDE	10		/* i		*/
#define MG_E_IES	11		/* ies		*/
#define MG_E_IR		12		/* ir		*/
#define MG_E_MATERIAL	13		/* m		*/
#define MG_E_NORMAL	14		/* n		*/
#define MG_E_OBJECT	15		/* o		*/
#define MG_E_POINT	16		/* p		*/
#define MG_E_PRISM	17		/* prism	*/
#define MG_E_RD		18		/* rd		*/
#define MG_E_RING	19		/* ring		*/
#define MG_E_RS		20		/* rs		*/
#define MG_E_SIDES	21		/* sides	*/
#define MG_E_SPH	22		/* sph		*/
#define MG_E_TD		23		/* td		*/
#define MG_E_TORUS	24		/* torus	*/
#define MG_E_TS		25		/* ts		*/
#define MG_E_VERTEX	26		/* v		*/
#define MG_E_XF		27		/* xf		*/
			/* end of Version 1 entities */
#define MG_E_FACEH	28		/* fh		*/
			/* end of Version 2 entities */

#define MG_NENTITIES	29		/* total # entities */

#define MG_NELIST	{28,29}		/* entity count for version 1 and up */

#define MG_NAMELIST	{"#","c","cct","cone","cmix","cspec","cxy","cyl","ed",\
			"f","i","ies","ir","m","n","o","p","prism","rd",\
			"ring","rs","sides","sph","td","torus","ts","v","xf",\
			"fh"}

#define MG_MAXELEN	6

extern char	mg_ename[MG_NENTITIES][MG_MAXELEN];

			/* Handler routines for each entity and unknown ones */
extern int	(*mg_ehand[MG_NENTITIES])(int argc, char **argv);
extern int	(*mg_uhand)(int argc, char **argv);
extern int	mg_defuhand(int, char **);

extern unsigned	mg_nunknown;		/* count of unknown entities */

			/* Error codes */
#define MG_OK		0		/* normal return value */
#define MG_EUNK		1		/* unknown entity */
#define MG_EARGC	2		/* wrong number of arguments */
#define MG_ETYPE	3		/* argument type error */
#define MG_EILL		4		/* illegal argument value */
#define MG_EUNDEF	5		/* undefined reference */
#define MG_ENOFILE	6		/* cannot open input file */
#define MG_EINCL	7		/* error in included file */
#define MG_EMEM		8		/* out of memory */
#define MG_ESEEK	9		/* file seek error */
#define MG_EBADMAT	10		/* bad material specification */
#define MG_ELINE	11		/* input line too long */
#define MG_ECNTXT	12		/* unmatched context close */

#define MG_NERRS	13

extern char	*mg_err[MG_NERRS];	/* list of error messages */

/*
 * The general process for running the parser is to fill in the mg_ehand
 * array with handlers for each entity you know how to handle.
 * Then, call mg_init to fill in the rest.  This function will report
 * an error and quit if you try to support an inconsistent set of entities.
 * For each file you want to parse, call mg_load with the file name.
 * To read from standard input, use NULL as the file name.
 * For additional control over error reporting and file management,
 * use mg_open, mg_read, mg_parse and mg_close instead of mg_load.
 * To pass an entity of your own construction to the parser, use
 * the mg_handle function rather than the mg_ehand routines directly.
 * (The first argument to mg_handle is the entity #, or -1.)
 * To free any data structures and clear the parser, use mg_clear.
 * If there is an error, mg_load, mg_open, mg_parse, mg_handle and
 * mg_fgoto will return an error from the list above.  In addition,
 * mg_load will report the error to stderr.  The mg_read routine
 * returns 0 when the end of file has been reached.
 */

#define MG_MAXLINE	4096		/* maximum input line length */
#define MG_MAXARGC	(MG_MAXLINE/4)	/* maximum argument count */

typedef struct mg_fctxt {
	char	fname[96];			/* file name */
	FILE	*fp;				/* stream pointer */
	int	fid;				/* unique file context id */
	char	inpline[MG_MAXLINE];		/* input line */
	int	lineno;				/* line number */
	struct mg_fctxt	*prev;			/* previous context */
} MG_FCTXT;

typedef struct {
	int	fid;				/* file this position is for */
	int	lineno;				/* line number in file */
	long	offset;				/* offset from beginning */
} MG_FPOS;

extern MG_FCTXT	*mg_file;		/* current file context */

extern void	mg_init(void);		/* fill in mg_ehand array */
extern int	mg_load(char *);	/* parse a file */
extern int	mg_open(MG_FCTXT *, char *);	/* open new input file */
extern int	mg_read(void);		/* read next line */
extern int	mg_parse(void);		/* parse current line */
extern void	mg_fgetpos(MG_FPOS *);	/* get position on input file */
extern int	mg_fgoto(MG_FPOS *);	/* go to position on input file */
extern void	mg_close(void);		/* close input file */
extern void	mg_clear(void);		/* clear parser */
extern int	mg_handle(int, int, char **);	/* handle an entity */

#ifndef MG_NQCD
#define MG_NQCD		5		/* default number of divisions */
#endif

extern int	mg_nqcdivs;		/* divisions per quarter circle */

/*
 * The following library routines are included for your convenience:
 */

extern int mg_entity(char *);		/* get entity number from its name */
extern int isint(char *);		/* non-zero if integer format */
extern int isintd(char *, char *);	/* same with delimiter set */
extern int isflt(char *);		/* non-zero if floating point format */
extern int isfltd(char *, char *);	/* same with delimiter set */
extern int isname(char *);		/* non-zero if legal identifier name */
extern int badarg(int, char **, char *);/* check argument format */
extern int e_include(int, char **);	/* expand include entity */
extern int e_pipe(int, char **);	/* expand piped command */
extern int e_sph(int, char **);		/* expand sphere as other entities */
extern int e_torus(int, char **);	/* expand torus as other entities */
extern int e_cyl(int, char **);		/* expand cylinder as other entities */
extern int e_ring(int, char **);	/* expand ring as other entities */
extern int e_cone(int, char **);	/* expand cone as other entities */
extern int e_prism(int, char **);	/* expand prism as other entities */
extern int e_faceh(int, char **);	/* expand face w/ holes as face */

/************************************************************************
 *	Definitions for 3-d vector manipulation functions
 */

#ifdef  SMLFLT
#define  RREAL		float
#define  FTINY		(1e-3)
#else
#define  RREAL		double
#define  FTINY		(1e-6)
#endif
#define  FHUGE		(1e10)

typedef RREAL  FVECT[3];

#define  VCOPY(v1,v2)	((v1)[0]=(v2)[0],(v1)[1]=(v2)[1],(v1)[2]=(v2)[2])
#define  DOT(v1,v2)	((v1)[0]*(v2)[0]+(v1)[1]*(v2)[1]+(v1)[2]*(v2)[2])
#define  VSUM(vr,v1,v2,f)	((vr)[0]=(v1)[0]+(f)*(v2)[0], \
				(vr)[1]=(v1)[1]+(f)*(v2)[1], \
				(vr)[2]=(v1)[2]+(f)*(v2)[2])

#define is0vect(v)	(DOT(v,v) <= FTINY*FTINY)

#define round0(x)	if (x <= FTINY && x >= -FTINY) x = 0

extern double	normalize(FVECT);	/* normalize a vector */
extern void	fcross(FVECT,FVECT,FVECT);/* cross product of two vectors */

/************************************************************************
 *	Definitions for context handling routines
 *	(materials, colors, vectors)
 */

#define C_1SIDEDTHICK	0.005		/* assumed thickness of 1-sided mat. */

typedef struct {
	int	clock;		/* incremented each change -- resettable */
	void	*client_data;	/* pointer to private client-owned data */
	int	sided;		/* 1 if surface is 1-sided, 0 for 2-sided */
	float	nr, ni;		/* index of refraction, real and imaginary */
	float	rd;		/* diffuse reflectance */
	C_COLOR	rd_c;		/* diffuse reflectance color */
	float	td;		/* diffuse transmittance */
	C_COLOR	td_c;		/* diffuse transmittance color */
	float	ed;		/* diffuse emittance */
	C_COLOR	ed_c;		/* diffuse emittance color */
	float	rs;		/* specular reflectance */
	C_COLOR	rs_c;		/* specular reflectance color */
	float	rs_a;		/* specular reflectance roughness */
	float	ts;		/* specular transmittance */
	C_COLOR	ts_c;		/* specular transmittance color */
	float	ts_a;		/* specular transmittance roughness */
} C_MATERIAL;		/* material context */

typedef struct {
	int	clock;		/* incremented each change -- resettable */
	void	*client_data;	/* pointer to private client-owned data */
	FVECT	p, n;		/* point and normal */
} C_VERTEX;		/* vertex context */

#define C_DEFMATERIAL	{1,NULL,0,1.,0.,0.,C_DEFCOLOR,0.,C_DEFCOLOR,0.,\
				C_DEFCOLOR,0.,C_DEFCOLOR,0.,0.,C_DEFCOLOR,0.}
#define C_DEFVERTEX	{1,NULL,{0.,0.,0.},{0.,0.,0.}}

extern C_COLOR		*c_ccolor;	/* the current color */
extern char		*c_ccname;	/* current color name */
extern C_MATERIAL	*c_cmaterial;	/* the current material */
extern char		*c_cmname;	/* current material name */
extern C_VERTEX		*c_cvertex;	/* the current vertex */
extern char		*c_cvname;	/* current vertex name */

extern int	c_hcolor(int, char **);		/* handle color entity */
extern int	c_hmaterial(int, char **);	/* handle material entity */
extern int	c_hvertex(int, char **);	/* handle vertex entity */
extern void	c_clearall(void);		/* clear context tables */
extern C_MATERIAL	*c_getmaterial(char *);	/* get a named material */
extern C_VERTEX	*c_getvert(char *);		/* get a named vertex */
extern C_COLOR	*c_getcolor(char *);		/* get a named color */

/*************************************************************************
 *	Definitions for hierarchical object name handler
 */

extern int	obj_nnames;		/* depth of name hierarchy */
extern char	**obj_name;		/* names in hierarchy */

extern int	obj_handler(int, char **);	/* handle an object entity */
extern void	obj_clear(void);		/* clear object stack */

/**************************************************************************
 *	Definitions for hierarchical transformation handler
 */

typedef RREAL  MAT4[4][4];

#define  copymat4(m4a,m4b)	(void)memcpy((char *)m4a,(char *)m4b,sizeof(MAT4))

#define  MAT4IDENT		{ {1.,0.,0.,0.}, {0.,1.,0.,0.}, \
				{0.,0.,1.,0.}, {0.,0.,0.,1.} }

extern MAT4  m4ident;

#define  setident4(m4)		copymat4(m4, m4ident)

				/* regular transformation */
typedef struct {
	MAT4	xfm;				/* transform matrix */
	RREAL	sca;				/* scalefactor */
}  XF;

#define identxf(xp)		(void)(setident4((xp)->xfm),(xp)->sca=1.0)

#define XF_MAXDIM	8		/* maximum array dimensions */

struct xf_array {
	MG_FPOS	spos;			/* starting position on input */
	int	ndim;			/* number of array dimensions */
	struct {
		short	i, n;		/* current count and maximum */
		char	arg[8];		/* string argument value */
	} aarg[XF_MAXDIM];
};

typedef struct xf_spec {
	long	xid;			/* unique transform id */
	short	xac;			/* context argument count */
	short	rev;			/* boolean true if vertices reversed */
	XF	xf;			/* cumulative transformation */
	struct xf_array	*xarr;		/* transformation array pointer */
	struct xf_spec	*prev;		/* previous transformation context */
} XF_SPEC;			/* followed by argument buffer */

extern XF_SPEC	*xf_context;			/* current transform context */
extern char	**xf_argend;			/* last transform argument */

#define xf_ac(xf)	((xf)==NULL ? 0 : (xf)->xac)
#define xf_av(xf)	(xf_argend - (xf)->xac)

#define xf_argc		xf_ac(xf_context)
#define xf_argv		xf_av(xf_context)

/*
 * The transformation handler should do most of the work that needs
 * doing.  Just pass it any xf entities, then use the associated
 * functions to transform and translate points, transform vectors
 * (without translation), rotate vectors (without scaling) and scale
 * values appropriately.
 *
 * The routines xf_xfmpoint, xf_xfmvect and xf_rotvect take two
 * 3-D vectors (which may be identical), transforms the second and
 * puts the result into the first.
 */


extern int	xf_handler(int, char **);	/* handle xf entity */
extern void	xf_xfmpoint(FVECT, FVECT);	/* transform point */
extern void	xf_xfmvect(FVECT, FVECT);	/* transform vector */
extern void	xf_rotvect(FVECT, FVECT);	/* rotate vector */
extern double	xf_scale(double);		/* scale a value */
extern void	xf_clear(void);			/* clear xf stack */

/* The following are support routines you probably won't call directly */

XF_SPEC		*new_xf(int, char **);		/* allocate new transform */
void		free_xf(XF_SPEC *);		/* free a transform */
int		xf_aname(struct xf_array *);	/* name this instance */
long		comp_xfid(MAT4);		/* compute unique ID */
extern void	multmat4(MAT4, MAT4, MAT4);	/* m4a = m4b X m4c */
extern void	multv3(FVECT, FVECT, MAT4);	/* v3a = v3b X m4 (vectors) */
extern void	multp3(FVECT, FVECT, MAT4);	/* p3a = p3b X m4 (points) */
extern int	xf(XF *, int, char **);		/* interpret transform spec. */

	/* cvrgb.c */
extern void mgf2rgb(C_COLOR *cin, double intensity, float cout[3]);

/************************************************************************
 *	Miscellaneous definitions
 */

#ifndef  PI
#ifdef	M_PI
#define	 PI		((double)M_PI)
#else
#define	 PI		3.14159265358979323846
#endif
#endif

#ifdef __cplusplus
}
#endif
#endif /* _MGF_PARSER_H_ */

