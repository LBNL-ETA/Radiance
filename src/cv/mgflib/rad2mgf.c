#ifndef lint
static const char	RCSid[] = "$Id: rad2mgf.c,v 1.1 2003/02/22 02:07:23 greg Exp $";
#endif
/*
 * Convert Radiance scene description to MGF
 */

#include "standard.h"
#include <ctype.h>
#include <string.h>
#include "object.h"
#include "color.h"
#include "lookup.h"

#define C_1SIDEDTHICK	0.005

int	o_face(), o_cone(), o_sphere(), o_ring(), o_cylinder();
int	o_instance(), o_illum();
int	o_plastic(), o_metal(), o_glass(), o_dielectric(),
	o_mirror(), o_trans(), o_light();

LUTAB	rmats = LU_SINIT(free,NULL);		/* defined material table */

LUTAB	rdispatch = LU_SINIT(NULL,NULL);	/* function dispatch table */

char	curmat[80];				/* current material */
char	curobj[128] = "Untitled";		/* current object name */

double	unit_mult = 1.;				/* units multiplier */

#define hasmult		(unit_mult < .999 || unit_mult > 1.001)

/*
 * Stuff for tracking and reusing vertices:
 */

char	VKFMT[] = "%+16.9e %+16.9e %+16.9e";
#define VKLEN		64

#define mkvkey(k,v)	sprintf(k, VKFMT, (v)[0], (v)[1], (v)[2])

#define NVERTS		256

long	vclock;		/* incremented at each vertex request */

struct vert {
	long	lused;		/* when last used (0 if unassigned) */
	FVECT	p;		/* track point position only */
} vert[NVERTS];		/* our vertex cache */

LUTAB	vertab = LU_SINIT(free,NULL);	/* our vertex lookup table */


main(argc, argv)
int	argc;
char	**argv;
{
	int	i;

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'd':			/* units */
			switch (argv[i][2]) {
			case 'm':			/* meters */
				unit_mult = 1.;
				break;
			case 'c':			/* centimeters */
				unit_mult = .01;
				break;
			case 'f':			/* feet */
				unit_mult = 12.*.0254;
				break;
			case 'i':			/* inches */
				unit_mult = .0254;
				break;
			default:
				goto unkopt;
			}
			break;
		default:
			goto unkopt;
		}
	init();
	if (i >= argc)
		rad2mgf(NULL);
	else
		for ( ; i < argc; i++)
			rad2mgf(argv[i]);
	uninit();
	exit(0);
unkopt:
	fprintf(stderr, "Usage: %s [-d{m|c|f|i}] file ..\n", argv[0]);
	exit(1);
}


rad2mgf(inp)		/* convert a Radiance file to MGF */
char	*inp;
{
#define mod	buf
#define typ	(buf+128)
#define id	(buf+256)
#define alias	(buf+384)
	char	buf[512];
	FUNARGS	fa;
	register FILE	*fp;
	register int	c;

	if (inp == NULL) {
		inp = "standard input";
		fp = stdin;
	} else if (inp[0] == '!') {
		if ((fp = popen(inp+1, "r")) == NULL) {
			fputs(inp, stderr);
			fputs(": cannot execute\n", stderr);
			exit(1);
		}
	} else if ((fp = fopen(inp, "r")) == NULL) {
		fputs(inp, stderr);
		fputs(": cannot open\n", stderr);
		exit(1);
	}
	printf("# Begin conversion from: %s\n", inp);
	while ((c = getc(fp)) != EOF)
		switch (c) {
		case ' ':		/* white space */
		case '\t':
		case '\n':
		case '\r':
		case '\f':
			break;
		case '#':		/* comment */
			if (fgets(buf, sizeof(buf), fp) != NULL)
				printf("# %s", buf);
			break;
		case '!':		/* inline command */
			ungetc(c, fp);
			fgetline(buf, sizeof(buf), fp);
			rad2mgf(buf);
			break;
		default:		/* Radiance primitive */
			ungetc(c, fp);
			if (fscanf(fp, "%s %s %s", mod, typ, id) != 3) {
				fputs(inp, stderr);
				fputs(": unexpected EOF\n", stderr);
				exit(1);
			}
			if (!strcmp(typ, "alias")) {
				strcpy(alias, "EOF");
				fscanf(fp, "%s", alias);
				newmat(id, alias);
			} else {
				if (!readfargs(&fa, fp)) {
					fprintf(stderr,
				"%s: bad argument syntax for %s \"%s\"\n",
							inp, typ, id);
					exit(1);
				}
				cvtprim(inp, mod, typ, id, &fa);
				freefargs(&fa);
			}
			break;
		}
	printf("# End conversion from: %s\n", inp);
	if (inp[0] == '!')
		pclose(fp);
	else
		fclose(fp);
#undef mod
#undef typ
#undef id
#undef alias
}


cvtprim(inp, mod, typ, id, fa)	/* process Radiance primitive */
char	*inp, *mod, *typ, *id;
FUNARGS	*fa;
{
	int	(*df)();

	df = (int (*)())lu_find(&rdispatch, typ)->data;
	if (df != NULL) {				/* convert */
		if ((*df)(mod, typ, id, fa) < 0) {
			fprintf(stderr, "%s: bad %s \"%s\"\n", typ, id);
			exit(1);
		}
	} else {					/* unsupported */
		o_unsupported(mod, typ, id, fa);
		if (lu_find(&rmats, mod)->data != NULL)	/* make alias */
			newmat(id, mod);
	}
}


newmat(id, alias)		/* add a modifier to the alias list */
char	*id;
char	*alias;
{
	register LUENT	*lp, *lpa;

	if (alias != NULL) {			/* look up alias */
		if ((lpa = lu_find(&rmats, alias)) == NULL)
			goto memerr;
		if (lpa->data == NULL)
			alias = NULL;		/* doesn't exist! */
	}
	if ((lp = lu_find(&rmats, id)) == NULL)	/* look up material */
		goto memerr;
	if (alias != NULL && lp->data == lpa->key)
		return;			/* alias set already */
	if (lp->data == NULL) {			/* allocate material */
		if ((lp->key = (char *)malloc(strlen(id)+1)) == NULL)
			goto memerr;
		strcpy(lp->key, id);
	}
	if (alias == NULL) {			/* set this material */
		lp->data = lp->key;
		printf("m %s =\n", id);
	} else {				/* set this alias */
		lp->data = lpa->key;
		printf("m %s = %s\n", id, alias);
	}
	strcpy(curmat, id);
	return;
memerr:
	fputs("Out of memory in newmat!\n", stderr);
	exit(1);
}


setmat(id)			/* set material to this one */
char	*id;
{
	if (!strcmp(id, curmat))	/* already set? */
		return;
	if (!strcmp(id, VOIDID))	/* cannot set */
		return;
	printf("m %s\n", id);
	strcpy(curmat, id);
}


setobj(id)			/* set object name to this one */
char	*id;
{
	register char	*cp, *cp2;
	char	*end = NULL;
	int	diff = 0;
				/* use all but final suffix */
	for (cp = id; *cp; cp++)
		if (*cp == '.')
			end = cp;
	if (end == NULL)
		end = cp;
				/* copy to current object */
	cp2 = curobj;
	if (!isalpha(*id)) {	/* start with letter */
		diff = *cp2 != 'O';
		*cp2++ = 'O';
	}
	for (cp = id; cp < end; *cp2++ = *cp++) {
		if (*cp < '!' | *cp > '~')	/* limit to visible chars */
			*cp = '?';
		diff += *cp != *cp2;
	}
	if (!diff && !*cp2)
		return;
	*cp2 = '\0';
	fputs("o\no ", stdout);
	puts(curobj);
}


init()			/* initialize dispatch table and output */
{
	lu_init(&vertab, NVERTS);
	lu_init(&rdispatch, 22);
	add2dispatch("polygon", o_face);
	add2dispatch("cone", o_cone);
	add2dispatch("cup", o_cone);
	add2dispatch("sphere", o_sphere);
	add2dispatch("bubble", o_sphere);
	add2dispatch("cylinder", o_cylinder);
	add2dispatch("tube", o_cylinder);
	add2dispatch("ring", o_ring);
	add2dispatch("instance", o_instance);
	add2dispatch("plastic", o_plastic);
	add2dispatch("plastic2", o_plastic);
	add2dispatch("metal", o_metal);
	add2dispatch("metal2", o_metal);
	add2dispatch("glass", o_glass);
	add2dispatch("dielectric", o_dielectric);
	add2dispatch("trans", o_trans);
	add2dispatch("trans2", o_trans);
	add2dispatch("mirror", o_mirror);
	add2dispatch("light", o_light);
	add2dispatch("spotlight", o_light);
	add2dispatch("glow", o_light);
	add2dispatch("illum", o_illum);
	puts("# The following was converted from RADIANCE scene input");
	if (hasmult)
		printf("xf -s %.4e\n", unit_mult);
	printf("o %s\n", curobj);
}


uninit()			/* mark end of MGF file */
{
	puts("o");
	if (hasmult)
		puts("xf");
	puts("# End of data converted from RADIANCE scene input");
	lu_done(&rdispatch);
	lu_done(&rmats);
	lu_done(&vertab);
}


clrverts()			/* clear vertex table */
{
	register int	i;

	lu_done(&vertab);
	for (i = 0; i < NVERTS; i++)
		vert[i].lused = 0;
	lu_init(&vertab, NVERTS);
}


add2dispatch(name, func)	/* add function to dispatch table */
char	*name;
int	(*func)();
{
	register LUENT	*lp;

	lp = lu_find(&rdispatch, name);
	if (lp->key != NULL) {
		fputs(name, stderr);
		fputs(": duplicate dispatch entry!\n", stderr);
		exit(1);
	}
	lp->key = name;
	lp->data = (char *)func;
}


char *
getvertid(vname, vp)		/* get/set vertex ID for this point */
char	*vname;
FVECT	vp;
{
	static char	vkey[VKLEN];
	register LUENT	*lp;
	register int	i, vndx;

	vclock++;			/* increment counter */
	mkvkey(vkey, vp);
	if ((lp = lu_find(&vertab, vkey)) == NULL)
		goto memerr;
	if (lp->data == NULL) {		/* allocate new vertex entry */
		if (lp->key != NULL)		/* reclaim deleted entry */
			vertab.ndel--;
		else {
			if ((lp->key = (char *)malloc(VKLEN)) == NULL)
				goto memerr;
			strcpy(lp->key, vkey);
		}
		vndx = 0;			/* find oldest vertex */
		for (i = 1; i < NVERTS; i++)
			if (vert[i].lused < vert[vndx].lused)
				vndx = i;
		if (vert[vndx].lused) {		/* free old entry first */
			mkvkey(vkey, vert[vndx].p);
			lu_delete(&vertab, vkey);
		}
		VCOPY(vert[vndx].p, vp);			/* assign it */
		printf("v v%d =\n\tp %.15g %.15g %.15g\n",	/* print it */
				vndx, vp[0], vp[1], vp[2]);
		lp->data = (char *)&vert[vndx];			/* set it */
	} else
		vndx = (struct vert *)lp->data - vert;
	vert[vndx].lused = vclock;		/* record this use */
	sprintf(vname, "v%d", vndx);
	return(vname);
memerr:
	fputs("Out of memory in getvertid!\n", stderr);
	exit(1);
}


int
o_unsupported(mod, typ, id, fa)		/* mark unsupported primitive */
char	*mod, *typ, *id;
FUNARGS	*fa;
{
	register int	i;

	fputs("\n# Unsupported RADIANCE primitive:\n", stdout);
	printf("# %s %s %s", mod, typ, id);
	printf("\n# %d", fa->nsargs);
	for (i = 0; i < fa->nsargs; i++)
		printf(" %s", fa->sarg[i]);
#ifdef IARGS
	printf("\n# %d", fa->niargs);
	for (i = 0; i < fa->niargs; i++)
		printf(" %ld", fa->iarg[i]);
#else
	fputs("\n# 0", stdout);
#endif
	printf("\n# %d", fa->nfargs);
	for (i = 0; i < fa->nfargs; i++)
		printf(" %g", fa->farg[i]);
	fputs("\n\n", stdout);
	return(0);
}


int
o_face(mod, typ, id, fa)		/* print out a polygon */
char	*mod, *typ, *id;
FUNARGS	*fa;
{
	char	entbuf[2048], *linestart;
	register char	*cp;
	register int	i;

	if (fa->nfargs < 9 | fa->nfargs % 3)
		return(-1);
	setmat(mod);
	setobj(id);
	cp = linestart = entbuf;
	*cp++ = 'f';
	for (i = 0; i < fa->nfargs; i += 3) {
		*cp++ = ' ';
		if (cp - linestart > 72) {
			*cp++ = '\\'; *cp++ = '\n';
			linestart = cp;
			*cp++ = ' '; *cp++ = ' ';
		}
		getvertid(cp, fa->farg + i);
		while (*cp)
			cp++;
	}
	puts(entbuf);
	return(0);
}


int
o_cone(mod, typ, id, fa)	/* print out a cone */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	char	v1[6], v2[6];

	if (fa->nfargs != 8)
		return(-1);
	setmat(mod);
	setobj(id);
	getvertid(v1, fa->farg);
	getvertid(v2, fa->farg + 3);
	if (typ[1] == 'u')			/* cup -> inverted cone */
		printf("cone %s %.12g %s %.12g\n",
				v1, -fa->farg[6], v2, -fa->farg[7]);
	else
		printf("cone %s %.12g %s %.12g\n",
				v1, fa->farg[6], v2, fa->farg[7]);
	return(0);
}


int
o_sphere(mod, typ, id, fa)	/* print out a sphere */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	char	cent[6];

	if (fa->nfargs != 4)
		return(-1);
	setmat(mod);
	setobj(id);
	printf("sph %s %.12g\n", getvertid(cent, fa->farg),
			typ[0]=='b' ? -fa->farg[3] : fa->farg[3]);
	return(0);
}


int
o_cylinder(mod, typ, id, fa)	/* print out a cylinder */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	char	v1[6], v2[6];

	if (fa->nfargs != 7)
		return(-1);
	setmat(mod);
	setobj(id);
	getvertid(v1, fa->farg);
	getvertid(v2, fa->farg + 3);
	printf("cyl %s %.12g %s\n", v1,
			typ[0]=='t' ? -fa->farg[6] : fa->farg[6], v2);
	return(0);
}


int
o_ring(mod, typ, id, fa)	/* print out a ring */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	if (fa->nfargs != 8)
		return(-1);
	setmat(mod);
	setobj(id);
	printf("v cent =\n\tp %.12g %.12g %.12g\n",
			fa->farg[0], fa->farg[1], fa->farg[2]);
	printf("\tn %.12g %.12g %.12g\n",
			fa->farg[3], fa->farg[4], fa->farg[5]);
	if (fa->farg[6] < fa->farg[7])
		printf("ring cent %.12g %.12g\n",
				fa->farg[6], fa->farg[7]);
	else
		printf("ring cent %.12g %.12g\n",
				fa->farg[7], fa->farg[6]);
	return(0);
}


int
o_instance(mod, typ, id, fa)	/* convert an instance */
char	*mod, *typ, *id;
FUNARGS	*fa;
{
	register int	i;
	register char	*cp;
	char	*start = NULL, *end = NULL;
	/*
	 * We don't really know how to do this, so we just create
	 * a reference to an undefined MGF file and it's the user's
	 * responsibility to create this file and put the appropriate
	 * stuff into it.
	 */
	if (fa->nsargs < 1)
		return(-1);
	setmat(mod);			/* only works if surfaces are void */
	setobj(id);
	for (cp = fa->sarg[0]; *cp; cp++)	/* construct MGF file name */
		if (*cp == '/')
			start = cp+1;
		else if (*cp == '.')
			end = cp;
	if (start == NULL)
		start = fa->sarg[0];
	if (end == NULL || start >= end)
		end = cp;
	fputs("i ", stdout);			/* print include entity */
	for (cp = start; cp < end; cp++)
		putchar(*cp);
	fputs(".mgf", stdout);			/* add MGF suffix */
	for (i = 1; i < fa->nsargs; i++) {	/* add transform */
		putchar(' ');
		fputs(fa->sarg[i], stdout);
	}
	putchar('\n');
	clrverts();			/* vertex id's no longer reliable */
	return(0);
}


int
o_illum(mod, typ, id, fa)	/* convert an illum material */
char	*mod, *typ, *id;
FUNARGS	*fa;
{
	if (fa->nsargs == 1 && strcmp(fa->sarg[0], VOIDID)) {
		newmat(id, fa->sarg[0]);	/* just create alias */
		return(0);
	}
					/* else create invisible material */
	newmat(id, NULL);
	puts("\tts 1 0");
	return(0);
}


int
o_plastic(mod, typ, id, fa)	/* convert a plastic material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb;
	double	d;

	if (fa->nfargs != (typ[7]=='2' ? 6 : 5))
		return(-1);
	newmat(id, NULL);
	rrgb[0] = fa->farg[0]; rrgb[1] = fa->farg[1]; rrgb[2] = fa->farg[2];
	rgb_cie(cxyz, rrgb);
	puts("\tc");				/* put diffuse component */
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\trd %.4f\n", cxyz[1]*(1. - fa->farg[3]));
	if (fa->farg[3] > FTINY) {		/* put specular component */
		puts("\tc");
		printf("\trs %.4f %.4f\n", fa->farg[3],
				typ[7]=='2' ? .5*(fa->farg[4] + fa->farg[5]) :
						fa->farg[4]);
	}
	return(0);
}


int
o_metal(mod, typ, id, fa)	/* convert a metal material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb;
	double	d;

	if (fa->nfargs != (typ[5]=='2' ? 6 : 5))
		return(-1);
	newmat(id, NULL);
	rrgb[0] = fa->farg[0]; rrgb[1] = fa->farg[1]; rrgb[2] = fa->farg[2];
	rgb_cie(cxyz, rrgb);
	puts("\tc");				/* put diffuse component */
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\trd %.4f\n", cxyz[1]*(1. - fa->farg[3]));
						/* put specular component */
	printf("\trs %.4f %.4f\n", cxyz[1]*fa->farg[3],
			typ[5]=='2' ? .5*(fa->farg[4] + fa->farg[5]) :
					fa->farg[4]);
	return(0);
}


int
o_glass(mod, typ, id, fa)	/* convert a glass material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb, trgb;
	double	nrfr = 1.52, F, d;
	register int	i;

	if (fa->nfargs != 3 && fa->nfargs != 4)
		return(-1);
	newmat(id, NULL);
	if (fa->nfargs == 4)
		nrfr = fa->farg[3];
	printf("\tir %f 0\n", nrfr);
	F = (1. - nrfr)/(1. + nrfr);		/* use normal incidence */
	F *= F;
	for (i = 0; i < 3; i++) {
		trgb[i] = fa->farg[i] * (1. - F)*(1. - F) /
				(1. - F*F*fa->farg[i]*fa->farg[i]);
		rrgb[i] = F * (1. + (1. - 2.*F)*fa->farg[i]) /
				(1. - F*F*fa->farg[i]*fa->farg[i]);
	}
	rgb_cie(cxyz, rrgb);			/* put reflected component */
	puts("\tc");
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\trs %.4f 0\n", cxyz[1]);
	rgb_cie(cxyz, trgb);			/* put transmitted component */
	puts("\tc");
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\tts %.4f 0\n", cxyz[1]);
	return(0);
}


int
o_dielectric(mod, typ, id, fa)	/* convert a dielectric material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, trgb;
	double	F, d;
	register int	i;

	if (fa->nfargs != 5)
		return(-1);
	newmat(id, NULL);
	F = (1. - fa->farg[3])/(1. + fa->farg[3]);	/* normal incidence */
	F *= F;
	for (i = 0; i < 3; i++)
		trgb[i] = (1. - F)*pow(fa->farg[i], C_1SIDEDTHICK/unit_mult);
	printf("\tir %f 0\n", fa->farg[3]);	/* put index of refraction */
	printf("\tsides 1\n");
	puts("\tc");				/* put reflected component */
	printf("\trs %.4f 0\n", F);
	rgb_cie(cxyz, trgb);			/* put transmitted component */
	puts("\tc");
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\tts %.4f 0\n", cxyz[1]);
	return(0);
}


int
o_mirror(mod, typ, id, fa)	/* convert a mirror material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb;
	double	d;

	if (fa->nsargs == 1) {			/* use alternate material */
		newmat(id, fa->sarg[0]);
		return(0);
	}
	if (fa->nfargs != 3)
		return(-1);
	newmat(id, NULL);
	rrgb[0] = fa->farg[0]; rrgb[1] = fa->farg[1]; rrgb[2] = fa->farg[2];
	rgb_cie(cxyz, rrgb);
	puts("\tc");				/* put specular component */
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\trs %.4f 0\n", cxyz[1]);
	return(0);
}


int
o_trans(mod, typ, id, fa)	/* convert a trans material */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb;
	double	rough, trans, tspec, d;

	if (typ[4] == '2') {		/* trans2 */
		if (fa->nfargs != 8)
			return(-1);
		rough = .5*(fa->farg[4] + fa->farg[5]);
		trans = fa->farg[6];
		tspec = fa->farg[7];
	} else {			/* trans */
		if (fa->nfargs != 7)
			return(-1);
		rough = fa->farg[4];
		trans = fa->farg[5];
		tspec = fa->farg[6];
	}
	newmat(id, NULL);
	rrgb[0] = fa->farg[0]; rrgb[1] = fa->farg[1]; rrgb[2] = fa->farg[2];
	rgb_cie(cxyz, rrgb);
	puts("\tc");				/* put transmitted diffuse */
	d = cxyz[0] + cxyz[1] + cxyz[2];
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\ttd %.4f\n", cxyz[1]*trans*(1. - fa->farg[3])*(1. - tspec));
						/* put transmitted specular */
	printf("\tts %.4f %.4f\n", cxyz[1]*trans*tspec*(1. - fa->farg[3]), rough);
						/* put reflected diffuse */
	printf("\trd %.4f\n", cxyz[1]*(1. - fa->farg[3])*(1. - trans));
	puts("\tc");				/* put reflected specular */
	printf("\trs %.4f %.4f\n", fa->farg[3], rough);
	return(0);
}


int
o_light(mod, typ, id, fa)		/* convert a light type */
char	*mod, *typ, *id;
register FUNARGS	*fa;
{
	COLOR	cxyz, rrgb;
	double	d;

	if (fa->nfargs < 3)
		return(-1);
	newmat(id, NULL);
	rrgb[0] = fa->farg[0]; rrgb[1] = fa->farg[1]; rrgb[2] = fa->farg[2];
	rgb_cie(cxyz, rrgb);
	d = cxyz[0] + cxyz[1] + cxyz[2];
	puts("\tc");
	if (d > FTINY)
		printf("\t\tcxy %.4f %.4f\n", cxyz[0]/d, cxyz[1]/d);
	printf("\ted %.4g\n", cxyz[1]*(PI*WHTEFFICACY));
	return(0);
}
