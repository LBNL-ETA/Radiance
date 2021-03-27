#ifndef lint
static const char RCSid[] = "$Id: testBSDF.c,v 1.16 2021/03/27 17:50:18 greg Exp $";
#endif
/*
 * Simple test program to demonstrate BSDF operation.
 *
 *	G. Ward		June 2015
 */

#define	_USE_MATH_DEFINES
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
	printf("  s[r|t][s|d] N theta phi\t Generate N ray directions & colors at given incidence\n");
	printf("  h[s|d] theta phi\t\t Report hemispherical scattering at given incidence\n");
	printf("  r[s|d] theta phi\t\t Report hemispherical reflection at given incidence\n");
	printf("  t[s|d] theta phi\t\t Report hemispherical transmission at given incidence\n");
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
	if (vp->cieY <= 1e-9) {
		printf("%s0 0 0\n", intro);
		return;
	}
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
	while (fgets(inp, sizeof(inp), stdin)) {
		int	sflags = SDsampAll;
		char	*cp = inp;
		char	*cp2;
		FVECT	vin, vout;
		double	proja[2];
		int	n, i;
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
			if (!bsdf)
				goto noBSDFerr;
			printf("Material: '%s'\n", bsdf->matn);
			printf("Manufacturer: '%s'\n", bsdf->makr);
			printf("Width, Height, Thickness (m): %.4e, %.4e, %.4e\n",
					bsdf->dim[0], bsdf->dim[1], bsdf->dim[2]);
			if (bsdf->mgf)
				printf("Has geometry: %lu bytes\n",
						(unsigned long)strlen(bsdf->mgf));
			else
				printf("Has geometry: no\n");
			break;
		case 'C':			/* report constant values */
			if (!bsdf)
				goto noBSDFerr;
			if (bsdf->rf)
				printf("Peak front hemispherical reflectance: %.3e\n",
						bsdf->rLambFront.cieY +
						bsdf->rf->maxHemi);
			if (bsdf->rb)
				printf("Peak back hemispherical reflectance: %.3e\n",
						bsdf->rLambBack.cieY +
						bsdf->rb->maxHemi);
			if (bsdf->tf)
				printf("Peak front hemispherical transmittance: %.3e\n",
						bsdf->tLambFront.cieY + bsdf->tf->maxHemi);
			if (bsdf->tb)
				printf("Peak back hemispherical transmittance: %.3e\n",
						bsdf->tLambBack.cieY + bsdf->tb->maxHemi);
			printXYZ("Diffuse Front Reflectance: ", &bsdf->rLambFront);
			printXYZ("Diffuse Back Reflectance: ", &bsdf->rLambBack);
			printXYZ("Diffuse Front Transmittance: ", &bsdf->tLambFront);
			printXYZ("Diffuse Back Transmittance: ", &bsdf->tLambBack);
			break;
		case 'Q':			/* query BSDF value */
			if (!bsdf)
				goto noBSDFerr;
			if (!*sskip2(cp,4))
				break;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			vec_from_deg(vout, atof(sskip2(cp,3)), atof(sskip2(cp,4)));
			if (!SDreportError(SDevalBSDF(&val, vout, vin, bsdf), stderr))
				printXYZ("", &val);
			break;
		case 'S':			/* sample BSDF */
			if (!bsdf)
				goto noBSDFerr;
			if (!*sskip2(cp,3))
				break;
			if (toupper(cp[1]) == 'R') {
				sflags &= ~SDsampT;
				++cp;
			} else if (toupper(cp[1]) == 'T') {
				sflags &= ~SDsampR;
				++cp;
			}
			if (toupper(cp[1]) == 'S')
				sflags &= ~SDsampDf;
			else if (toupper(cp[1]) == 'D')
				sflags &= ~SDsampSp;
			i = n = atoi(sskip2(cp,1));
			vec_from_deg(vin, atof(sskip2(cp,2)), atof(sskip2(cp,3)));
			while (i-- > 0) {
				VCOPY(vout, vin);
				if (SDreportError(SDsampBSDF(&val, vout,
						(i+rand()*(1./(RAND_MAX+.5)))/(double)n,
						sflags, bsdf), stderr))
					break;
				printf("%.8f %.8f %.8f ", vout[0], vout[1], vout[2]);
				printXYZ("", &val);
			}
			break;
		case 'H':			/* hemispherical values */
		case 'R':
		case 'T':
			if (!bsdf)
				goto noBSDFerr;
			if (!*sskip2(cp,2))
				break;
			if (toupper(cp[0]) == 'R')
				sflags &= ~SDsampT;
			else if (toupper(cp[0]) == 'T')
				sflags &= ~SDsampR;
			if (toupper(cp[1]) == 'S')
				sflags &= ~SDsampDf;
			else if (toupper(cp[1]) == 'D')
				sflags &= ~SDsampSp;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			printf("%.4e\n", SDdirectHemi(vin, sflags, bsdf));
			break;
		case 'A':			/* resolution in proj. steradians */
			if (!bsdf)
				goto noBSDFerr;
			if (!*sskip2(cp,2))
				break;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			if (*sskip2(cp,4)) {
				vec_from_deg(vout, atof(sskip2(cp,3)), atof(sskip2(cp,4)));
				if (SDreportError(SDsizeBSDF(proja, vout, vin,
						SDqueryMin+SDqueryMax, bsdf), stderr))
					continue;
			} else if (SDreportError(SDsizeBSDF(proja, vin, NULL,
						SDqueryMin+SDqueryMax, bsdf), stderr))
					continue;
			printf("%.4e %.4e\n", proja[0], proja[1]);
			break;
		default:
			Usage(argv[0]);
			break;
		}
		fflush(stdout);			/* in case we're on remote */
		continue;
noBSDFerr:
		fprintf(stderr, "%s: First, use 'L' command to load BSDF\n", argv[0]);
	}
	return 0;
}
