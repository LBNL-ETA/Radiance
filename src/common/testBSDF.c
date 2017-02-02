#ifndef lint
static const char RCSid[] = "$Id: testBSDF.c,v 1.7 2017/02/02 00:30:16 greg Exp $";
#endif
/*
 * Simple test program to demonstrate BSDF operation.
 *
 *	G. Ward		June 2015
 */

#define	_USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "rtio.h"
#include "bsdf.h"

static void
Usage(const char *prog)
{
	printf("Usage: %s [bsdf_directory]\n", prog);
	printf("Input commands:\n");
	printf("  L bsdf.xml\t\t\t Load (make active) given BSDF input file\n");
	printf("  i\t\t\t\t Report general information (metadata)\n");
	printf("  c\t\t\t\t Report diffuse and specular components\n");
	printf("  q theta_i phi_i theta_o phi_o\t Query BSDF for given path (CIE-XYZ)\n");
	printf("  s N theta phi\t\t\t Generate N ray directions at given incidence\n");
	printf("  h theta phi\t\t\t Report hemispherical total at given incidence\n");
	printf("  r theta phi\t\t\t Report hemispherical reflection at given incidence\n");
	printf("  t theta phi\t\t\t Report hemispherical transmission at given incidence\n");
	printf("  a theta phi [t2 p2]\t\t Report resolution (in proj. steradians) for given direction(s)\n");
	printf("  ^D\t\t\t\t Quit program\n");
}

static void
vec_from_deg(FVECT v, double theta, double phi)
{
	const double	DEG = M_PI/180.;

	theta *= DEG; phi *= DEG;
	v[0] = v[1] = sin(theta);
	v[0] *= cos(phi);
	v[1] *= sin(phi);
	v[2] = cos(theta);
}

static void
printXYZ(const char *intro, const SDValue *vp)
{
	printf("%s%.3e %.3e %.3e\n", intro,
			vp->spec.cx/vp->spec.cy*vp->cieY,
			vp->cieY,
			(1.-vp->spec.cx-vp->spec.cy)/
			vp->spec.cy*vp->cieY);
}

int
main(int argc, char *argv[])
{
	const char	*directory = NULL;
	char		inp[512], path[512];
	const SDData	*bsdf = NULL;

	if (argc > 2 || (argc == 2 && argv[1][0] == '-')) {
		Usage(argv[0]);
		return 1;
	}
	if (argc == 2)
		directory = argv[1];
	
	SDretainSet = SDretainBSDFs;		/* keep BSDFs in memory */

						/* loop on command */
	while (fgets(inp, sizeof(inp), stdin) != NULL) {
		int	sflags = SDsampAll;
		char	*cp = inp;
		char	*cp2;
		FVECT	vin, vout;
		double	proja[2];
		int	n;
		SDValue	val;
		
		while (isspace(*cp)) cp++;

		switch (toupper(*cp)) {
		case 'L':			/* load/activate BSDF input */
			cp2 = cp = sskip2(cp, 1);
			if (!*cp)
				break;
			while (*cp) cp++;
			while (isspace(*--cp)) *cp = '\0';
			if (directory)
				sprintf(path, "%s/%s", directory, cp2);
			else
				strcpy(path, cp2);
			if (bsdf)
				SDfreeCache(bsdf);
			bsdf = SDcacheFile(path);
			continue;
		case 'I':			/* report general info. */
			if (bsdf == NULL)
				goto noBSDFerr;
			printf("Material: '%s'\n", bsdf->matn);
			printf("Manufacturer: '%s'\n", bsdf->makr);
			printf("Width, Height, Thickness (m): %.4e, %.4e, %.4e\n",
					bsdf->dim[0], bsdf->dim[1], bsdf->dim[2]);
			printf("Has geometry: %s\n", bsdf->mgf!=NULL ? "yes" : "no");
			continue;
		case 'C':			/* report constant values */
			if (bsdf == NULL)
				goto noBSDFerr;
			if (bsdf->rf != NULL)
				printf("Peak front hemispherical reflectance: %.3e\n",
						bsdf->rLambFront.cieY +
						bsdf->rf->maxHemi);
			if (bsdf->rb != NULL)
				printf("Peak back hemispherical reflectance: %.3e\n",
						bsdf->rLambBack.cieY +
						bsdf->rb->maxHemi);
			if (bsdf->tf != NULL)
				printf("Peak front hemispherical transmittance: %.3e\n",
						bsdf->tLamb.cieY + bsdf->tf->maxHemi);
			if (bsdf->tb != NULL)
				printf("Peak back hemispherical transmittance: %.3e\n",
						bsdf->tLamb.cieY + bsdf->tb->maxHemi);
			printXYZ("Diffuse Front Reflectance: ", &bsdf->rLambFront);
			printXYZ("Diffuse Back Reflectance: ", &bsdf->rLambBack);
			printXYZ("Diffuse Transmittance: ", &bsdf->tLamb);
			continue;
		case 'Q':			/* query BSDF value */
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,4))
				break;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			vec_from_deg(vout, atof(sskip2(cp,3)), atof(sskip2(cp,4)));
			if (!SDreportError(SDevalBSDF(&val, vout, vin, bsdf), stderr))
				printXYZ("", &val);
			continue;
		case 'S':			/* sample BSDF */
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,3))
				break;
			n = atoi(sskip2(cp,1));
			vec_from_deg(vin, atof(sskip2(cp,2)), atof(sskip2(cp,3)));
			while (n-- > 0) {
				if (SDreportError(SDsampBSDF(&val, vin,
						rand()*(1./(RAND_MAX+.5)),
						sflags, bsdf), stderr))
					break;
				printf("%.8f %.8f %.8f\n", vin[0], vin[1], vin[2]);
			}
			continue;
		case 'H':			/* hemispherical totals */
		case 'R':
		case 'T':
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,2))
				break;
			if (tolower(*cp) == 'r')
				sflags &= ~SDsampT;
			else if (tolower(*cp) == 't')
				sflags &= ~SDsampR;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			printf("%.4e\n", SDdirectHemi(vin, sflags, bsdf));
			continue;
		case 'A':			/* resolution in proj. steradians */
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,2))
				break;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			if (*sskip2(cp,4)) {
				vec_from_deg(vout, atof(sskip2(cp,3)), atof(sskip2(cp,4)));
				if (SDreportError(SDsizeBSDF(proja, vin, vout,
						SDqueryMin+SDqueryMax, bsdf), stderr))
					continue;
			} else if (SDreportError(SDsizeBSDF(proja, vin, NULL,
						SDqueryMin+SDqueryMax, bsdf), stderr))
					continue;
			printf("%.4e %.4e\n", proja[0], proja[1]);
			continue;
		}
		Usage(argv[0]);
		continue;
noBSDFerr:
		fprintf(stderr, "%s: First, use 'L' command to load BSDF\n", argv[0]);
	}
	return 0;
}
