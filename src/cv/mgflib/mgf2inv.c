#ifndef lint
static const char	RCSid[] = "$Id: mgf2inv.c,v 1.9 2003/07/27 22:12:02 schorsch Exp $";
#endif
/*
 * Convert MGF to Inventor file.
 *
 *	December 1995	Greg Ward
 */

#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include <ctype.h>

#include <string.h>

#include "parser.h"

#include "lookup.h"

#define O_INV1		1	/* Inventor 1.0 output */
#define O_INV2		2	/* Inventor 2.0 output */
#define O_VRML1		3	/* VRML 1.0 output */

#define MAXID		48	/* maximum identifier length */

#define VERTFMT		"%+16.9e %+16.9e %+16.9e\n%+6.3f %+6.3f %+6.3f"
#define VZVECT		"+0.000 +0.000 +0.000"
#define VFSEPPOS	50	/* position of newline in above */
#define VFLEN		72	/* total vertex string length */
#define MAXVERT		10240	/* maximum cached vertices */

#define setvkey(k,v)	sprintf(k,VERTFMT,(v)->p[0],(v)->p[1],(v)->p[2],\
					(v)->n[0],(v)->n[1],(v)->n[2]);

char	vlist[MAXVERT][VFLEN];	/* our vertex cache */
int	nverts;			/* current cache size */

LUTAB	vert_tab = LU_SINIT(NULL,NULL);

struct face {
	struct face	*next;		/* next face in list */
	short		nv;		/* number of vertices */
	short		vl[3];		/* vertex index list (variable) */
}	*flist, *flast;		/* our face cache */

#define newface(n)	(struct face *)malloc(sizeof(struct face) + \
				((n) > 3 ? (n)-3 : 0)*sizeof(short))
#define freeface(f)	free((MEM_PTR)f)

#define	TABSTOP		8	/* assumed number of characters per tab */
#define	SHIFTW		2	/* nesting shift width */
#define MAXIND		15	/* maximum indent level */

char	tabs[MAXIND*SHIFTW+1];	/* current tab-in string */

#define curmatname	(c_cmname == NULL ? "mat" : to_id(c_cmname))

int	outtype = O_INV2;	/* output format */

extern int	i_comment(), i_object(), i_xf(),
		i_cyl(), i_face(), i_sph();

extern char	*to_id();


main(argc, argv)
int	argc;
char	*argv[];
{
	int	i;
				/* initialize dispatch table */
	mg_ehand[MG_E_COMMENT] = i_comment;	/* we pass comments */
	mg_ehand[MG_E_COLOR] = c_hcolor;	/* they get color */
	mg_ehand[MG_E_CMIX] = c_hcolor;		/* they mix colors */
	mg_ehand[MG_E_CSPEC] = c_hcolor;	/* they get spectra */
	mg_ehand[MG_E_CXY] = c_hcolor;		/* they get chromaticities */
	mg_ehand[MG_E_CCT] = c_hcolor;		/* they get color temp's */
	mg_ehand[MG_E_CYL] = i_cyl;		/* we do cylinders */
	mg_ehand[MG_E_ED] = c_hmaterial;	/* they get emission */
	mg_ehand[MG_E_FACE] = i_face;		/* we do faces */
	mg_ehand[MG_E_MATERIAL] = c_hmaterial;	/* they get materials */
	mg_ehand[MG_E_NORMAL] = c_hvertex;	/* they get normals */
	mg_ehand[MG_E_OBJECT] = i_object;	/* we track object names */
	mg_ehand[MG_E_POINT] = c_hvertex;	/* they get points */
	mg_ehand[MG_E_RD] = c_hmaterial;	/* they get diffuse refl. */
	mg_ehand[MG_E_RS] = c_hmaterial;	/* they get specular refl. */
	mg_ehand[MG_E_SIDES] = c_hmaterial;	/* they get # sides */
	mg_ehand[MG_E_SPH] = i_sph;		/* we do spheres */
	mg_ehand[MG_E_TD] = c_hmaterial;	/* they get diffuse trans. */
	mg_ehand[MG_E_TS] = c_hmaterial;	/* they get specular trans. */
	mg_ehand[MG_E_VERTEX] = c_hvertex;	/* they get vertices */
	mg_ehand[MG_E_XF] = i_xf;		/* we track transforms */
	mg_init();		/* initialize the parser */
				/* get options and print format line */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		if (!strcmp(argv[i], "-vrml"))
			outtype = O_VRML1;
		else if (!strcmp(argv[i], "-1"))
			outtype = O_INV1;
		else if (!strcmp(argv[i], "-2"))
			outtype = O_INV2;
		else
			goto userr;
	switch (outtype) {
	case O_INV1:
		printf("#Inventor V1.0 ascii\n");
		break;
	case O_INV2:
		printf("#Inventor V2.0 ascii\n");
		break;
	case O_VRML1:
		printf("#VRML V1.0 ascii\n");
		break;
	}
	printf("## Translated from MGF Version %d.%d\n", MG_VMAJOR, MG_VMINOR);
	printf("Separator {\n");		/* begin root node */
						/* general properties */
	printf("MaterialBinding { value OVERALL }\n");
	printf("NormalBinding { value PER_VERTEX_INDEXED }\n");
	if (outtype != O_INV1) {
		printf("ShapeHints {\n");
		printf("\tvertexOrdering CLOCKWISE\n");
		printf("\tfaceType UNKNOWN_FACE_TYPE\n");
		printf("}\n");
	}
	if (i == argc) {	/* load standard input */
		if (mg_load(NULL) != MG_OK)
			exit(1);
		if (mg_nunknown)
			printf("## %s: %u unknown entities\n",
					argv[0], mg_nunknown);
	}
				/* load MGF files */
	for ( ; i < argc; i++) {
		printf("## %s %s ##############################\n",
				argv[0], argv[i]);
		mg_nunknown = 0;
		if (mg_load(argv[i]) != MG_OK)
			exit(1);
		if (mg_nunknown)
			printf("## %s %s: %u unknown entities\n",
					argv[0], argv[i], mg_nunknown);
	}
	flush_cache();		/* flush face cache, just in case */
	printf("}\n");		/* close root node */
	exit(0);
userr:
	fprintf(stderr, "%s: [-1|-2|-vrml] [file] ..\n", argv[0]);
	exit(1);
}


indent(deeper)				/* indent in or out */
int	deeper;
{
	static int	level;		/* current nesting level */
	register int	i;
	register char	*cp;

	if (deeper) level++;		/* in or out? */
	else if (level > 0) level--;
					/* compute actual shift */
	if ((i = level) > MAXIND) i = MAXIND;
	cp = tabs;
	for (i *= SHIFTW; i >= TABSTOP; i -= TABSTOP)
		*cp++ = '\t';
	while (i--)
		*cp++ = ' ';
	*cp = '\0';
}


int
i_comment(ac, av)			/* transfer comment as is */
int	ac;
char	**av;
{
	fputs(tabs, stdout);
	putchar('#');			/* Inventor comment character */
	while (--ac > 0) {
		putchar(' ');
		fputs(*++av, stdout);
	}
	putchar('\n');
	return(MG_OK);
}


int
i_object(ac, av)			/* group object name */
int	ac;
char	**av;
{
	static int	objnest;

	flush_cache();			/* flush cached objects */
	if (ac == 2) {				/* start group */
		printf("%sDEF %s Group {\n", tabs, to_id(av[1]));
		indent(1);
		objnest++;
		return(MG_OK);
	}
	if (ac == 1) {				/* end group */
		if (--objnest < 0)
			return(MG_ECNTXT);
		indent(0);
		fputs(tabs, stdout);
		fputs("}\n", stdout);
		return(MG_OK);
	}
	return(MG_EARGC);
}


int
i_xf(ac, av)				/* transform object(s) */
int	ac;
char	**av;
{
	static long	xfid;
	register XF_SPEC	*spec;

	flush_cache();			/* flush cached objects */
	if (ac == 1) {			/* end of transform */
		if ((spec = xf_context) == NULL)
			return(MG_ECNTXT);
		indent(0);			/* close original segment */
		printf("%s}\n", tabs);
		indent(0);
		printf("%s}\n", tabs);
		if (spec->xarr != NULL) {	/* check for iteration */
			register struct xf_array	*ap = spec->xarr;
			register int	n;

			ap->aarg[ap->ndim-1].i = 1;	/* iterate array */
			for ( ; ; ) {
				n = ap->ndim-1;
				while (ap->aarg[n].i < ap->aarg[n].n) {
					sprintf(ap->aarg[n].arg, "%d",
							ap->aarg[n].i);
					printf("%sSeparator {\n", tabs);
					indent(1);
					(void)put_xform(spec);
					printf("%sUSE _xf%ld\n", tabs,
							spec->xid);
					indent(0);
					printf("%s}\n", tabs);
					++ap->aarg[n].i;
				}
				ap->aarg[n].i = 0;
				(void)strcpy(ap->aarg[n].arg, "0");
				while (n-- && ++ap->aarg[n].i >= ap->aarg[n].n) {
					ap->aarg[n].i = 0;
					(void)strcpy(ap->aarg[n].arg, "0");
				}
				if (n < 0)
					break;
				sprintf(ap->aarg[n].arg, "%d", ap->aarg[n].i);
			}
		}
						/* pop transform */
		xf_context = spec->prev;
		free_xf(spec);
		return(MG_OK);
	}
					/* else allocate new transform */
	if ((spec = new_xf(ac-1, av+1)) == NULL)
		return(MG_EMEM);
	spec->xid = ++xfid;		/* assign unique ID */
	spec->prev = xf_context;	/* push onto stack */
	xf_context = spec;
					/* translate xf specification */
	printf("%sSeparator {\n", tabs);
	indent(1);
	if (put_xform(spec) < 0)
		return(MG_ETYPE);
	printf("%sDEF _xf%ld Group {\n", tabs, spec->xid);	/* begin */
	indent(1);
	return(MG_OK);
}


int
put_xform(spec)			/* translate and print transform */
register XF_SPEC	*spec;
{
	register char	**av;
	register int	n;

	n = xf_ac(spec) - xf_ac(spec->prev);
	if (xf(&spec->xf, n, av=xf_av(spec)) != n)
		return(-1);
	printf("%sMatrixTransform {\n", tabs);
	indent(1);
	printf("%s# xf", tabs);		/* put out original as comment */
	while (n--) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	putchar('\n');			/* put out computed matrix */
	printf("%smatrix %13.9g %13.9g %13.9g %13.9g\n", tabs,
			spec->xf.xfm[0][0], spec->xf.xfm[0][1],
			spec->xf.xfm[0][2], spec->xf.xfm[0][3]);
	for (n = 1; n < 4; n++)
		printf("%s       %13.9g %13.9g %13.9g %13.9g\n", tabs,
				spec->xf.xfm[n][0], spec->xf.xfm[n][1],
				spec->xf.xfm[n][2], spec->xf.xfm[n][3]);
	indent(0);
	printf("%s}\n", tabs);
	return(0);
}


int
put_material()			/* put out current material */
{
	char	*mname = curmatname;
	float	rgbval[3];

	if (!c_cmaterial->clock) {	/* current, just use it */
		printf("%sUSE %s\n", tabs, mname);
		return(0);
	}
				/* else update definition */
	printf("%sDEF %s Group {\n", tabs, mname);
	indent(1);
	printf("%sMaterial {\n", tabs);
	indent(1);
	mgf2rgb(&c_cmaterial->rd_c, c_cmaterial->rd, rgbval);
	printf("%sambientColor %.4f %.4f %.4f\n", tabs,
			rgbval[0], rgbval[1], rgbval[2]);
	printf("%sdiffuseColor %.4f %.4f %.4f\n", tabs,
			rgbval[0], rgbval[1], rgbval[2]);
	if (c_cmaterial->rs > FTINY) {
		mgf2rgb(&c_cmaterial->rs_c, c_cmaterial->rs, rgbval);
		printf("%sspecularColor %.4f %.4f %.4f\n", tabs,
				rgbval[0], rgbval[1], rgbval[2]);
		printf("%sshininess %.3f\n", tabs, 1.-sqrt(c_cmaterial->rs_a));
	}
	if (c_cmaterial->ed > FTINY) {
		mgf2rgb(&c_cmaterial->ed_c, 1.0, rgbval);
		printf("%semissiveColor %.4f %.4f %.4f\n", tabs,
				rgbval[0], rgbval[1], rgbval[2]);
	}
	if (c_cmaterial->ts > FTINY)
		printf("%stransparency %.4f\n", tabs,
				c_cmaterial->ts + c_cmaterial->td);
	indent(0);
	printf("%s}\n", tabs);
	if (outtype != O_INV1)
		printf("%sShapeHints { shapeType %s }\n", tabs,
			c_cmaterial->sided ? "SOLID" : "UNKNOWN_SHAPE_TYPE");
	indent(0);
	printf("%s}\n", tabs);
	c_cmaterial->clock = 0;
	return(0);
}


int
i_face(ac, av)			/* translate an N-sided face */
int	ac;
char	**av;
{
	static char	lastmat[MAXID];
	struct face	*newf;
	register C_VERTEX	*vp;
	register LUENT	*lp;
	register int	i;

	if (ac < 4)
		return(MG_EARGC);
	if ( strcmp(lastmat, curmatname) || c_cmaterial->clock ||
			nverts == 0 || nverts+ac-1 >= MAXVERT) {
		flush_cache();			/* new cache */
		lu_init(&vert_tab, MAXVERT);
		printf("%sSeparator {\n", tabs);
		indent(1);
		if (put_material() < 0)		/* put out material */
			return(MG_EBADMAT);
		(void)strcpy(lastmat, curmatname);
	}
				/* allocate new face */
	if ((newf = newface(ac-1)) == NULL)
		return(MG_EMEM);
	newf->nv = ac-1;
				/* get vertex references */
	for (i = 0; i < newf->nv; i++) {
		if ((vp = c_getvert(av[i+1])) == NULL)
			return(MG_EUNDEF);
		setvkey(vlist[nverts], vp);
		lp = lu_find(&vert_tab, vlist[nverts]);
		if (lp == NULL)
			return(MG_EMEM);
		if (lp->key == NULL)
			lp->key = (char *)vlist[nverts++];
		newf->vl[i] = ((char (*)[VFLEN])lp->key - vlist);
	}
				/* add to face list */
	newf->next = NULL;
	if (flist == NULL)
		flist = newf;
	else
		flast->next = newf;
	flast = newf;
	return(MG_OK);		/* we'll actually put it out later */
}


int
i_sph(ac, av)			/* translate sphere description */
int	ac;
char	**av;
{
	register C_VERTEX	*cent;

	if (ac != 3)
		return(MG_EARGC);
	flush_cache();		/* flush vertex cache */
	printf("%sSeparator {\n", tabs);
	indent(1);
				/* put out current material */
	if (put_material() < 0)
		return(MG_EBADMAT);
				/* get center */
	if ((cent = c_getvert(av[1])) == NULL)
		return(MG_EUNDEF);
				/* get radius */
	if (!isflt(av[2]))
		return(MG_ETYPE);
	printf("%sTranslation { translation %13.9g %13.9g %13.9g }\n", tabs,
			cent->p[0], cent->p[1], cent->p[2]);
	printf("%sSphere { radius %s }\n", tabs, av[2]);
	indent(0);
	printf("%s}\n", tabs);
	return(MG_OK);
}


int
i_cyl(ac, av)			/* translate a cylinder description */
int	ac;
char	**av;
{
	register C_VERTEX	*v1, *v2;
	FVECT	va;
	double	length, angle;

	if (ac != 4)
		return(MG_EARGC);
	flush_cache();		/* flush vertex cache */
	printf("%sSeparator {\n", tabs);
	indent(1);
				/* put out current material */
	if (put_material() < 0)
		return(MG_EBADMAT);
				/* get endpoints */
	if (((v1 = c_getvert(av[1])) == NULL) | ((v2 = c_getvert(av[3])) == NULL))
		return(MG_EUNDEF);
				/* get radius */
	if (!isflt(av[2]))
		return(MG_ETYPE);
				/* compute transform */
	va[0] = v2->p[0] - v1->p[0];
	va[1] = v2->p[1] - v1->p[1];
	va[2] = v2->p[2] - v1->p[2];
	length = sqrt(DOT(va,va));
	if (va[1] >= length)
		angle = 0.;
	else if (va[1] <= -length)
		angle = PI;
	else
		angle = acos(va[1]/length);
	printf("%sTranslation { translation %13.9g %13.9g %13.9g }\n", tabs,
			.5*(v1->p[0]+v2->p[0]), .5*(v1->p[1]+v2->p[1]),
			.5*(v1->p[2]+v2->p[2]));
	printf("%sRotation { rotation %.9g %.9g %.9g %.9g }\n", tabs,
			va[2], 0., -va[0], angle);
				/* open-ended */
	printf("%sCylinder { parts SIDES height %13.9g radius %s }\n", tabs,
			length, av[2]);
	indent(0);
	printf("%s}\n", tabs);
	return(MG_OK);
}


char *
to_id(name)			/* make sure a name is a valid Inventor ID */
register char	*name;
{
	static char	id[MAXID];
	register char	*cp;

	for (cp = id; *name && cp < MAXID-1+id; name++)
		if (isalnum(*name) || *name == '_')
			*cp++ = *name;
		else
			*cp++ = '_';
	*cp = '\0';
	return(id);
}


flush_cache()			/* put out cached faces */
{
	int	donorms = 0;
	register struct face	*f;
	register int	i;

	if (nverts == 0)
		return;
					/* put out coordinates */
	printf("%sCoordinate3 {\n", tabs);
	indent(1);
	vlist[0][VFSEPPOS] = '\0';
	printf("%spoint [ %s", tabs, vlist[0]);
	for (i = 1; i < nverts; i++) {
		vlist[i][VFSEPPOS] = '\0';
		printf(",\n%s        %s", tabs, vlist[i]);
		if (strcmp(VFSEPPOS+1+vlist[i], VZVECT))
			donorms++;
	}
	indent(0);
	printf(" ]\n%s}\n", tabs);
	if (donorms) {			/* put out normals */
		printf("%sNormal {\n", tabs);
		indent(1);
		printf("%svector [ %s", tabs, VFSEPPOS+1+vlist[0]);
		for (i = 1; i < nverts; i++)
			printf(",\n%s         %s", tabs, VFSEPPOS+1+vlist[i]);
		indent(0);
		printf(" ]\n%s}\n", tabs);
	}
					/* put out faces */
	printf("%sIndexedFaceSet {\n", tabs);
	indent(1);
	f = flist;			/* coordinate indices */
	printf("%scoordIndex [ %d", tabs, f->vl[0]);
	for (i = 1; i < f->nv; i++)
		printf(", %d", f->vl[i]);
	for (f = f->next; f != NULL; f = f->next) {
		printf(", -1,\n%s             %d", tabs, f->vl[0]);
		for (i = 1; i < f->nv; i++)
			printf(", %d", f->vl[i]);
	}
	printf(" ]\n");
	if (donorms) {
		f = flist;			/* normal indices */
		printf("%snormalIndex [ %d", tabs, f->vl[0]);
		for (i = 1; i < f->nv; i++)
			printf(", %d", f->vl[i]);
		for (f = f->next; f != NULL; f = f->next) {
			printf(", -1,\n%s              %d", tabs, f->vl[0]);
			for (i = 1; i < f->nv; i++)
				printf(", %d", f->vl[i]);
		}
		printf(" ]\n");
	}
	indent(0);			/* close IndexedFaceSet */
	printf("%s}\n", tabs);
	indent(0);			/* close face group */
	printf("%s}\n", tabs);
	while ((f = flist) != NULL) {	/* free face list */
		flist = f->next;
		freeface(f);
	}
	lu_done(&vert_tab);		/* clear lookup table */
	nverts = 0;
}
