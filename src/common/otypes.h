/* RCSid $Id: otypes.h,v 2.17 2012/07/29 19:01:39 greg Exp $ */
/*
 *  otypes.h - defines for object types.
 */
#ifndef _RAD_OTYPES_H_
#define _RAD_OTYPES_H_
#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	char  *funame;			/* function name */
	int  flags;			/* type flags */
#ifdef FUN_ARGLIST
	int  (*funp)(FUN_ARGLIST);	/* pointer to function */
#else
	int  (*funp)();			/* pointer to function */
#endif
}  FUN;

#ifdef FUN_ARGLIST
extern int  o_default(FUN_ARGLIST);
#else
extern int  o_default(); /* XXX conflict with radogl.h */
#endif

				/* object types in decreasing frequency */
#define  OBJ_FACE	0		/* polygon */
#define  OBJ_CONE	1		/* cone */
#define  OBJ_SPHERE	2		/* sphere */
#define  TEX_FUNC	3		/* surface texture function */
#define  OBJ_RING	4		/* disk */
#define  OBJ_CYLINDER	5		/* cylinder */
#define  OBJ_INSTANCE	6		/* octree instance */
#define  OBJ_CUP	7		/* inverted cone */
#define  OBJ_BUBBLE	8		/* inverted sphere */
#define  OBJ_TUBE	9		/* inverted cylinder */
#define  OBJ_MESH	10		/* mesh instance */
#define  MOD_ALIAS	11		/* modifier alias */
#define  MAT_PLASTIC	12		/* plastic surface */
#define  MAT_METAL	13		/* metal surface */
#define  MAT_GLASS	14		/* thin glass surface */
#define  MAT_TRANS	15		/* translucent material */
#define  MAT_DIELECTRIC	16		/* dielectric material */
#define  MAT_PLASTIC2	17		/* anisotropic plastic */
#define  MAT_METAL2	18		/* anisotropic metal */
#define  MAT_TRANS2	19		/* anisotropic translucent material */
#define  MAT_INTERFACE	20		/* dielectric interface */
#define  MAT_PFUNC	21		/* plastic brdf function */
#define  MAT_MFUNC	22		/* metal brdf function */
#define  PAT_BFUNC	23		/* brightness function */
#define  PAT_BDATA	24		/* brightness data */
#define  PAT_BTEXT	25		/* monochromatic text */
#define  PAT_CPICT	26		/* color picture */
#define  MAT_GLOW	27		/* proximity light source */
#define  OBJ_SOURCE	28		/* distant source */
#define  MAT_LIGHT	29		/* primary light source */
#define  MAT_ILLUM	30		/* secondary light source */
#define  MAT_SPOT	31		/* spot light source */
#define  MAT_MIST	32		/* mist medium */
#define  MAT_MIRROR	33		/* mirror (secondary source) */
#define  MAT_TFUNC	34		/* trans brdf function */
#define  MAT_BRTDF	35		/* BRTD function */
#define  MAT_BSDF	36		/* BSDF data file */
#define  MAT_PDATA	37		/* plastic brdf data */
#define  MAT_MDATA	38		/* metal brdf data */
#define  MAT_TDATA	39		/* trans brdf data */
#define  PAT_CFUNC	40		/* color function */
#define  MAT_CLIP	41		/* clipping surface */
#define  PAT_CDATA	42		/* color data */
#define  PAT_CTEXT	43		/* colored text */
#define  TEX_DATA	44		/* surface texture data */
#define  MIX_FUNC	45		/* mixing function */
#define  MIX_DATA	46		/* mixing data */
#define  MIX_TEXT	47		/* mixing text */
#define  MIX_PICT	48		/* mixing picture */
#define  MAT_DIRECT1	49		/* unidirecting material */
#define  MAT_DIRECT2	50		/* bidirecting material */
#define  MAT_ASHIKHMIN	51		/* Ashikhmin-Shirley BRDF material */
				/* number of object types */
#define  NUMOTYPE	52
				/* type flags */
#define  T_S		01		/* surface (object) */
#define  T_M		02		/* material */
#define  T_P		04		/* pattern */
#define  T_T		010		/* texture */
#define  T_X		020		/* mixture */
#define  T_V		040		/* volume */
#define  T_L		0100		/* light source modifier */
#define  T_LV		0200		/* virtual light source modifier */
#define  T_F		0400		/* function */
#define  T_D		01000		/* data */
#define  T_I		02000		/* picture */
#define  T_E		04000		/* text */
				/* user-defined types */
#define  T_SP1		010000
#define  T_SP2		020000
#define  T_SP3		040000

extern FUN  ofun[];			/* our type list */

#define  issurface(t)	(ofun[t].flags & T_S)
#define  isvolume(t)	(ofun[t].flags & T_V)
#define  ismodifier(t)	(!(ofun[t].flags & (T_S|T_V)))
#define  ismaterial(t)	(ofun[t].flags & T_M)
#define  ispattern(t)	(ofun[t].flags & T_P)
#define  istexture(t)	(ofun[t].flags & T_T)
#define  ismixture(t)	(ofun[t].flags & T_X)
#define  islight(t)	(ofun[t].flags & T_L)
#define  isvlight(t)	(ofun[t].flags & T_LV)
#define  hasdata(t)	(ofun[t].flags & (T_D|T_I))
#define  hasfunc(t)	(ofun[t].flags & (T_F|T_D|T_I))
#define  hastext(t)	(ofun[t].flags & T_E)
#define  isflat(t)	((t)==OBJ_FACE || (t)==OBJ_RING)

#define  ALIASKEY	"alias"			/* alias keyword */
#define  ALIASMOD	"inherit"		/* inherit target modifier */

					/* type list initialization */
#define  INIT_OTYPE	{	{ "polygon",	T_S,		o_default }, \
				{ "cone",	T_S,		o_default }, \
				{ "sphere",	T_S,		o_default }, \
				{ "texfunc",	T_T|T_F,	o_default }, \
				{ "ring",	T_S,		o_default }, \
				{ "cylinder",	T_S,		o_default }, \
				{ "instance",	T_V,		o_default }, \
				{ "cup",	T_S,		o_default }, \
				{ "bubble",	T_S,		o_default }, \
				{ "tube",	T_S,		o_default }, \
				{ "mesh",	T_V,		o_default }, \
				{ ALIASKEY,	0,		o_default }, \
				{ "plastic",	T_M,		o_default }, \
				{ "metal",	T_M,		o_default }, \
				{ "glass",	T_M,		o_default }, \
				{ "trans",	T_M,		o_default }, \
				{ "dielectric",	T_M,		o_default }, \
				{ "plastic2",	T_M|T_F,	o_default }, \
				{ "metal2",	T_M|T_F,	o_default }, \
				{ "trans2",	T_M|T_F,	o_default }, \
				{ "interface",	T_M,		o_default }, \
				{ "plasfunc",	T_M|T_F,	o_default }, \
				{ "metfunc",	T_M|T_F,	o_default }, \
				{ "brightfunc",	T_P|T_F,	o_default }, \
				{ "brightdata",	T_P|T_D,	o_default }, \
				{ "brighttext",	T_P|T_E,	o_default }, \
				{ "colorpict",	T_P|T_I,	o_default }, \
				{ "glow",	T_M|T_L,	o_default }, \
				{ "source",	T_S,		o_default }, \
				{ "light",	T_M|T_L,	o_default }, \
				{ "illum",	T_M|T_L,	o_default }, \
				{ "spotlight",	T_M|T_L,	o_default }, \
				{ "mist",	T_M,		o_default }, \
				{ "mirror",	T_M|T_LV,	o_default }, \
				{ "transfunc",	T_M|T_F,	o_default }, \
				{ "BRTDfunc",	T_M|T_F,	o_default }, \
				{ "BSDF",	T_M|T_D,	o_default }, \
				{ "plasdata",	T_M|T_D,	o_default }, \
				{ "metdata",	T_M|T_D,	o_default }, \
				{ "transdata",	T_M|T_D,	o_default }, \
				{ "colorfunc",	T_P|T_F,	o_default }, \
				{ "antimatter",	T_M,		o_default }, \
				{ "colordata",	T_P|T_D,	o_default }, \
				{ "colortext",	T_P|T_E,	o_default }, \
				{ "texdata",	T_T|T_D,	o_default }, \
				{ "mixfunc",	T_X|T_F,	o_default }, \
				{ "mixdata",	T_X|T_D,	o_default }, \
				{ "mixtext",	T_X|T_E,	o_default }, \
				{ "mixpict",	T_X|T_I,	o_default }, \
				{ "prism1",	T_M|T_F|T_LV,	o_default }, \
				{ "prism2",	T_M|T_F|T_LV,	o_default }, \
				{ "ashik2",	T_M|T_F,	o_default }, \
			}


#ifdef __cplusplus
}
#endif
#endif /* _RAD_OTYPES_H_ */
