#ifndef lint
static const char RCSid[] = "$Id: testBSDF.c,v 1.1 2015/06/05 17:55:18 greg Exp $";
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
	printf("  l bsdf.xml\t\t\t Load (make active) given BSDF input file\n");
	printf("  q theta_i phi_i theta_o phi_o\t Query BSDF for given path\n");
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

int
main(int argc, char *argv[])
{
	const char	*directory = "";
	char		inp[512], path[512];
	const SDData	*bsdf = NULL;

	if (argc > 2 || (argc == 2 && argv[1][0] == '-')) {
		Usage(argv[0]);
		return 1;
	}
	if (argc == 2)
		directory = argv[1];
	
	SDretainSet = SDretainBSDFs;		/* keep BSDFs loaded */

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

		switch (*cp) {
		case 'l':
		case 'L':			/* load/activate BSDF input */
			cp2 = cp = sskip2(cp, 1);
			if (!*cp)
				break;
			while (*cp) cp++;
			while (isspace(*--cp)) *cp = '\0';
			if (directory[0])
				sprintf(path, "%s/%s", directory, cp2);
			else
				strcpy(path, cp2);
			bsdf = SDcacheFile(path);
			continue;
		case 'q':
		case 'Q':			/* query BSDF value */
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,4))
				break;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			vec_from_deg(vout, atof(sskip2(cp,3)), atof(sskip2(cp,4)));
			if (!SDreportError(SDevalBSDF(&val, vout, vin, bsdf), stderr))
				printf("%.3e\n", val.cieY);
			continue;
		case 's':
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
		case 'h':
		case 'H':			/* hemispherical totals */
		case 'r':
		case 'R':
		case 't':
		case 'T':
			if (bsdf == NULL)
				goto noBSDFerr;
			if (!*sskip2(cp,2))
				break;
			if (tolower(*cp) == 'r')
				sflags ^= SDsampT;
			else if (tolower(*cp) == 't')
				sflags ^= SDsampR;
			vec_from_deg(vin, atof(sskip2(cp,1)), atof(sskip2(cp,2)));
			printf("%.4e\n", SDdirectHemi(vin, sflags, bsdf));
			continue;
		case 'a':
		case 'A':			/* resolution in degrees */
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
		fprintf(stderr, "%s: need to use 'l' command to load BSDF\n", argv[0]);
	}
	return 0;
}
