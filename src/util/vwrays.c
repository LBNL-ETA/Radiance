#ifndef lint
static const char	RCSid[] = "$Id: vwrays.c,v 3.5 2003/06/26 00:58:11 schorsch Exp $";
#endif
/*
 * Compute rays corresponding to a given picture or view.
 */


#include "standard.h"

#include "view.h"

extern int	putf(), putd(), puta();

int	(*putr)() = puta;

VIEW	vw = STDVIEW;

RESOLU	rs = {PIXSTANDARD, 512, 512};

double	pa = 1.;

int	zfd = -1;

int	fromstdin = 0;

char	*progname;


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*err;
	int	rval, getdim = 0;
	register int	i;

	progname = argv[0];
	if (argc < 2)
		goto userr;
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'f':			/* output format */
			switch (argv[i][2]) {
			case 'a':			/* ASCII */
				putr = puta;
				break;
			case 'f':			/* float */
				putr = putf;
				break;
			case 'd':			/* double */
				putr = putd;
				break;
			default:
				goto userr;
			}
			break;
		case 'v':			/* view file or option */
			if (argv[i][2] == 'f') {
				rval = viewfile(argv[++i], &vw, NULL);
				if (rval <= 0) {
					fprintf(stderr,
						"%s: no view in file\n",
							argv[i]);
					exit(1);
				}
				break;
			}
			rval = getviewopt(&vw, argc-i, argv+i);
			if (rval < 0)
				goto userr;
			i += rval;
			break;
		case 'd':			/* report dimensions only */
			getdim++;
			break;
		case 'x':			/* x resolution */
			rs.xr = atoi(argv[++i]);
			if (rs.xr <= 0) {
				fprintf(stderr, "%s: bad x resolution\n",
						progname);
				exit(1);
			}
			break;
		case 'y':			/* y resolution */
			rs.yr = atoi(argv[++i]);
			if (rs.yr <= 0) {
				fprintf(stderr, "%s: bad y resolution\n",
						progname);
				exit(1);
			}
			break;
		case 'p':			/* pixel aspect ratio */
			pa = atof(argv[++i]);
			break;
		case 'i':			/* get pixels from stdin */
			fromstdin = 1;
			break;
		default:
			goto userr;
		}
	if (i > argc | i+2 < argc)
		goto userr;
	if (i < argc) {
		rval = viewfile(argv[i], &vw, &rs);
		if (rval <= 0) {
			fprintf(stderr, "%s: no view in picture\n", argv[i]);
			exit(1);
		}
		if (i+1 < argc) {
			zfd = open(argv[i+1], O_RDONLY);
			if (zfd < 0) {
				fprintf(stderr,
					"%s: cannot open depth buffer\n",
						argv[i+1]);
				exit(1);
			}
		}
	}
	if ((err = setview(&vw)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	if (i == argc)
		normaspect(viewaspect(&vw), &pa, &rs.xr, &rs.yr);
	if (getdim) {
		printf("-x %d -y %d -ld%c\n", rs.xr, rs.yr,
				vw.vaft > FTINY ? '+' : '-');
		exit(0);
	}
	if (fromstdin)
		pix2rays(stdin);
	else
		putrays();
	exit(0);
userr:
	fprintf(stderr,
	"Usage: %s [ -i -f{a|f|d} | -d ] { view opts .. | picture [zbuf] }\n",
			progname);
	exit(1);
}


pix2rays(FILE *fp)
{
	static FVECT	rorg, rdir;
	float	zval;
	double	px, py;
	int	pp[2];
	double	d;
	register int	i;

	while (fscanf(fp, "%lf %lf", &px, &py) == 2) {
		if (px < 0 || px >= rs.xr ||
				py < 0 || py >= rs.yr) {
			fprintf(stderr,
				"%s: (x,y) pair (%.0f,%.0f) out of range\n",
					progname, px, py);
			exit(1);
		}
		if (zfd >= 0) {
			loc2pix(pp, &rs, px/rs.xr, py/rs.yr);
			if (lseek(zfd,
				(pp[1]*scanlen(&rs)+pp[0])*sizeof(float), 0)
					< 0 ||
					read(zfd, &zval, sizeof(float))
					< sizeof(float)) {
				fprintf(stderr, "%s: depth buffer read error\n",
						progname);
				exit(1);
			}
		}
		d = viewray(rorg, rdir, &vw, px/rs.xr, py/rs.yr);
		if (d < -FTINY)
			rorg[0] = rorg[1] = rorg[2] =
			rdir[0] = rdir[1] = rdir[2] = 0.;
		else if (zfd >= 0)
			for (i = 0; i < 3; i++) {
				rorg[i] += rdir[i]*zval;
				rdir[i] = -rdir[i];
			}
		else if (d > FTINY) {
			rdir[0] *= d; rdir[1] *= d; rdir[2] *= d;
		}
		(*putr)(rorg, rdir);
	}
	if (!feof(fp)) {
		fprintf(stderr, "%s: expected px py on input\n", progname);
		exit(1);
	}
}


putrays()
{
	static RREAL	loc[2];
	static FVECT	rorg, rdir;
	float	*zbuf;
	int	sc;
	double	d;
	register int	si, i;

	if (zfd >= 0) {
		zbuf = (float *)malloc(scanlen(&rs)*sizeof(float));
		if (zbuf == NULL) {
			fprintf(stderr, "%s: not enough memory\n", progname);
			exit(1);
		}
	}
	for (sc = 0; sc < numscans(&rs); sc++) {
		if (zfd >= 0) {
			if (read(zfd, zbuf, scanlen(&rs)*sizeof(float)) <
					scanlen(&rs)*sizeof(float)) {
				fprintf(stderr, "%s: depth buffer read error\n",
						progname);
				exit(1);
			}
		}
		for (si = 0; si < scanlen(&rs); si++) {
			pix2loc(loc, &rs, si, sc);
			d = viewray(rorg, rdir, &vw, loc[0], loc[1]);
			if (d < -FTINY)
				rorg[0] = rorg[1] = rorg[2] =
				rdir[0] = rdir[1] = rdir[2] = 0.;
			else if (zfd >= 0)
				for (i = 0; i < 3; i++) {
					rorg[i] += rdir[i]*zbuf[si];
					rdir[i] = -rdir[i];
				}
			else if (d > FTINY) {
				rdir[0] *= d; rdir[1] *= d; rdir[2] *= d;
			}
			(*putr)(rorg, rdir);
		}
	}
	if (zfd >= 0)
		free((void *)zbuf);
}


puta(ro, rd)		/* put out ray in ASCII format */
FVECT	ro, rd;
{
	printf("%.5e %.5e %.5e %.5e %.5e %.5e\n",
			ro[0], ro[1], ro[2],
			rd[0], rd[1], rd[2]);
}


putf(ro, rd)		/* put out ray in float format */
FVECT	ro, rd;
{
	float v[6];

	v[0] = ro[0]; v[1] = ro[1]; v[2] = ro[2];
	v[3] = rd[0]; v[4] = rd[1]; v[5] = rd[2];
	fwrite(v, sizeof(float), 6, stdout);
}


putd(ro, rd)		/* put out ray in double format */
FVECT	ro, rd;
{
	double v[6];

	v[0] = ro[0]; v[1] = ro[1]; v[2] = ro[2];
	v[3] = rd[0]; v[4] = rd[1]; v[5] = rd[2];
	fwrite(v, sizeof(double), 6, stdout);
}
