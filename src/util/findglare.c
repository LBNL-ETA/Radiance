#ifndef lint
static const char	RCSid[] = "$Id: findglare.c,v 2.11 2004/01/02 12:51:54 schorsch Exp $";
#endif
/*
 * Find glare sources in a scene or image.
 *
 *	Greg Ward	March 1991
 */

#include "glare.h"

#define FEQ(a,b)	((a)-(b)<=FTINY&&(b)-(a)<=FTINY)
#define VEQ(v1,v2)	(FEQ((v1)[0],(v2)[0])&&FEQ((v1)[1],(v2)[1]) \
				&&FEQ((v1)[2],(v2)[2]))

char	*rtargv[64] = {"rtrace", "-h-", "-ov", "-fff", "-ld-", "-i-", "-I-"};
int	rtargc = 7;

VIEW	ourview = STDVIEW;		/* our view */
VIEW	pictview = STDVIEW;		/* picture view */
VIEW	leftview, rightview;		/* leftmost and rightmost views */

char	*picture = NULL;		/* picture file name */
char	*octree = NULL;			/* octree file name */

int	verbose = 0;			/* verbose reporting */
char	*progname;			/* global argv[0] */

double	threshold = 0.;			/* glare threshold */

int	sampdens = SAMPDENS;		/* sample density */
ANGLE	glarang[180] = {AEND};		/* glare calculation angles */
int	nglarangs = 0;
double	maxtheta;			/* maximum angle (in radians) */
int	hsize;				/* horizontal size */

struct illum	*indirect;		/* array of indirect illuminances */

long	npixinvw;			/* number of pixels in view */
long	npixmiss;			/* number of pixels missed */

static int angcmp(const void	*ap1, const void	*ap2);
static void init(void);
static void cleanup(void);
static void printsources(void);
static void printillum(void);



int
main(
	int	argc,
	char	*argv[]
)
{
	int	combine = 1;
	int	gotview = 0;
	int	rval, i;
	char	*err;

	progname = argv[0];
					/* process options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++) {
						/* expand arguments */
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			fprintf(stderr, "%s: cannot expand '%s'",
					argv[0], argv[i]);
			exit(1);
		}
		rval = getviewopt(&ourview, argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			gotview++;
			continue;
		}
		switch (argv[i][1]) {
		case 't':
			threshold = atof(argv[++i]);
			break;
		case 'r':
			sampdens = atoi(argv[++i])/2;
			break;
		case 'v':
			if (argv[i][2] == '\0') {
				verbose++;
				break;
			}
			if (argv[i][2] != 'f')
				goto userr;
			rval = viewfile(argv[++i], &ourview, NULL);
			if (rval < 0) {
				fprintf(stderr,
				"%s: cannot open view file \"%s\"\n",
						progname, argv[i]);
				exit(1);
			} else if (rval == 0) {
				fprintf(stderr,
					"%s: bad view file \"%s\"\n",
						progname, argv[i]);
				exit(1);
			} else
				gotview++;
			break;
		case 'g':
			if (argv[i][2] != 'a')
				goto userr;
			if (setscan(glarang, argv[++i]) < 0) {
				fprintf(stderr, "%s: bad angle spec \"%s\"\n",
					progname, argv[i]);
				exit(1);
			}
			break;
		case 'p':
			picture = argv[++i];
			break;
		case 'c':
			combine = !combine;
			break;
		case 'd':
			rtargv[rtargc++] = argv[i];
			if (argv[i][2] != 'v')
				rtargv[rtargc++] = argv[++i];
			break;
		case 'l':
			if (argv[i][2] == 'd')
				break;
			/* FALL THROUGH */
		case 's':
		case 'P':
			rtargv[rtargc++] = argv[i];
			rtargv[rtargc++] = argv[++i];
			break;
		case 'w':
			rtargv[rtargc++] = argv[i];
			break;
		case 'a':
			rtargv[rtargc++] = argv[i];
			if (argv[i][2] == 'v') {
				rtargv[rtargc++] = argv[++i];
				rtargv[rtargc++] = argv[++i];
			}
			rtargv[rtargc++] = argv[++i];
			break;
		default:
			goto userr;
		}
	}
						/* get octree */
	if (i < argc-1)
		goto userr;
	if (i == argc-1) {
		rtargv[rtargc++] = octree = argv[i];
		rtargv[rtargc] = NULL;
	}
						/* get view */
	if (picture != NULL) {
		rval = viewfile(picture, &pictview, NULL);
		if (rval < 0) {
			fprintf(stderr, "%s: cannot open picture file \"%s\"\n",
					progname, picture);
			exit(1);
		} else if (rval == 0) {
			fprintf(stderr,
				"%s: cannot get view from picture \"%s\"\n",
					progname, picture);
			exit(1);
		}
		if (pictview.type == VT_PAR) {
			fprintf(stderr, "%s: %s: cannot use parallel view\n",
					progname, picture);
			exit(1);
		}
		if ((err = setview(&pictview)) != NULL) {
			fprintf(stderr, "%s: %s\n", picture, err);
			exit(1);
		}
	}
	if (!gotview) {
		if (picture == NULL) {
			fprintf(stderr, "%s: must have view or picture\n",
					progname);
			exit(1);
		}
		ourview = pictview;
	} else if (picture != NULL && !VEQ(ourview.vp, pictview.vp)) {
		fprintf(stderr, "%s: picture must have same viewpoint\n",
				progname);
		exit(1);
	}
	ourview.type = VT_HEM;
	ourview.horiz = ourview.vert = 180.0;
	ourview.hoff = ourview.voff = 0.0;
	fvsum(ourview.vdir, ourview.vdir, ourview.vup,
			-DOT(ourview.vdir,ourview.vup));
	if ((err = setview(&ourview)) != NULL) {
		fprintf(stderr, "%s: %s\n", progname, err);
		exit(1);
	}
	if (octree == NULL && picture == NULL) {
		fprintf(stderr,
		"%s: must specify at least one of picture or octree\n",
				progname);
		exit(1);
	}
	init();					/* initialize program */
	if (threshold <= FTINY)
		comp_thresh();			/* compute glare threshold */
	analyze();				/* analyze view */
	if (combine)
		absorb_specks();		/* eliminate tiny sources */
	cleanup();				/* tidy up */
						/* print header */
	newheader("RADIANCE", stdout);
	printargs(argc, argv, stdout);
	fputs(VIEWSTR, stdout);
	fprintview(&ourview, stdout);
	printf("\n");
	fputformat("ascii", stdout);
	printf("\n");
	printsources();				/* print glare sources */
	printillum();				/* print illuminances */
	exit(0);
userr:
	fprintf(stderr,
"Usage: %s [view options][-ga angles][-p picture][[rtrace options] octree]\n",
			progname);
	exit(1);
}


static int
angcmp(		/* compare two angles */
	const void	*ap1,
	const void	*ap2
)
{
	register int	a1, a2;

	a1 = *(ANGLE *)ap1;
	a2 = *(ANGLE *)ap2;
	if (a1 == a2) {
		fprintf(stderr, "%s: duplicate glare angle (%d)\n",
				progname, a1);
		exit(1);
	}
	return(a1-a2);
}


static void
init(void)				/* initialize global variables */
{
	double	d;
	register int	i;

	if (verbose)
		fprintf(stderr, "%s: initializing data structures...\n",
				progname);
						/* set direction vectors */
	for (i = 0; glarang[i] != AEND; i++)
		;
	qsort(glarang, i, sizeof(ANGLE), angcmp);
	if (i > 0 && (glarang[0] <= 0 || glarang[i-1] > 180)) {
		fprintf(stderr, "%s: glare angles must be between 1 and 180\n",
				progname);
		exit(1);
	}
	nglarangs = i;
	/* nglardirs = 2*nglarangs + 1; */
	/* vsize = sampdens - 1; */
	if (nglarangs > 0)
		maxtheta = (PI/180.)*glarang[nglarangs-1];
	else
		maxtheta = 0.0;
	hsize = hlim(0) + sampdens - 1;
	if (hsize > (int)(PI*sampdens))
		hsize = PI*sampdens;
	indirect = (struct illum *)calloc(nglardirs, sizeof(struct illum));
	if (indirect == NULL)
		memerr("indirect illuminances");
	npixinvw = npixmiss = 0L;
	leftview = ourview;
	rightview = ourview;
	spinvector(leftview.vdir, ourview.vdir, ourview.vup, maxtheta);
	spinvector(rightview.vdir, ourview.vdir, ourview.vup, -maxtheta);
	setview(&leftview);
	setview(&rightview);
	indirect[nglarangs].lcos =
	indirect[nglarangs].rcos = cos(maxtheta);
	indirect[nglarangs].rsin =
	-(indirect[nglarangs].lsin = sin(maxtheta));
	indirect[nglarangs].theta = 0.0;
	for (i = 0; i < nglarangs; i++) {
		d = (glarang[nglarangs-1] - glarang[i])*(PI/180.);
		indirect[nglarangs-i-1].lcos =
		indirect[nglarangs+i+1].rcos = cos(d);
		indirect[nglarangs+i+1].rsin =
		-(indirect[nglarangs-i-1].lsin = sin(d));
		d = (glarang[nglarangs-1] + glarang[i])*(PI/180.);
		indirect[nglarangs-i-1].rcos =
		indirect[nglarangs+i+1].lcos = cos(d);
		indirect[nglarangs-i-1].rsin =
		-(indirect[nglarangs+i+1].lsin = sin(d));
		indirect[nglarangs-i-1].theta = (PI/180.)*glarang[i];
		indirect[nglarangs+i+1].theta = -(PI/180.)*glarang[i];
	}
						/* open picture */
	if (picture != NULL) {
		if (verbose)
			fprintf(stderr, "%s: opening picture file \"%s\"...\n",
					progname, picture);
		open_pict(picture);
	}
						/* start rtrace */
	if (octree != NULL) {
		if (verbose) {
			fprintf(stderr,
				"%s: starting luminance calculation...\n\t",
					progname);
			printargs(rtargc, rtargv, stderr);
		}
		fork_rtrace(rtargv);
	}
}


static void
cleanup(void)				/* close files, wait for children */
{
	if (verbose)
		fprintf(stderr, "%s: cleaning up...        \n", progname);
	if (picture != NULL)
		close_pict();
	if (octree != NULL)
		done_rtrace();
	if (npixinvw < 100*npixmiss)
		fprintf(stderr, "%s: warning -- missing %d%% of samples\n",
				progname, (int)(100L*npixmiss/npixinvw));
}


extern int
compdir(			/* compute direction for x,y */
	FVECT	vd,
	int	x,
	int	y
)
{
	int	hl;
	FVECT	org;			/* dummy variable */

	hl = hlim(y);
	if (x <= -hl) {			/* left region */
		if (x <= -hl-sampdens)
			return(-1);
		return(viewray(org, vd, &leftview,
				(double)(x+hl)/(2*sampdens)+.5,
				(double)y/(2*sampdens)+.5));
	}
	if (x >= hl) {			/* right region */
		if (x >= hl+sampdens)
			return(-1);
		return(viewray(org, vd, &rightview,
				(double)(x-hl)/(2*sampdens)+.5,
				(double)y/(2*sampdens)+.5));
	}
					/* central region */
	if (viewray(org, vd, &ourview, .5, (double)y/(2*sampdens)+.5) < 0)
		return(-1);
	spinvector(vd, vd, ourview.vup, h_theta(x,y));
	return(0);
}


extern double
pixsize(		/* return the solid angle of pixel at (x,y) */
	int	x,
	int	y
)
{
	register int	hl, xo;
	double	disc;

	hl = hlim(y);
	if (x < -hl)
		xo = x+hl;
	else if (x > hl)
		xo = x-hl;
	else
		xo = 0;
	disc = 1. - (double)((long)xo*xo + (long)y*y)/((long)sampdens*sampdens);
	if (disc <= FTINY*FTINY)
		return(0.);
	return(1./(sampdens*sampdens*sqrt(disc)));
}


extern void
memerr(			/* malloc failure */
	char	*s
)
{
	fprintf(stderr, "%s: out of memory for %s\n", progname, s);
	exit(1);
}


static void
printsources(void)			/* print out glare sources */
{
	register struct source	*sp;

	printf("BEGIN glare source\n");
	for (sp = donelist; sp != NULL; sp = sp->next)
		printf("\t%f %f %f\t%f\t%f\n",
				sp->dir[0], sp->dir[1], sp->dir[2],
				sp->dom, sp->brt);
	printf("END glare source\n");
}


static void
printillum(void)			/* print out indirect illuminances */
{
	register int	i;

	printf("BEGIN indirect illuminance\n");
	for (i = 0; i < nglardirs; i++)
		if (indirect[i].n > FTINY)
			printf("\t%.0f\t%f\n", (180.0/PI)*indirect[i].theta,
					PI * indirect[i].sum / indirect[i].n);
	printf("END indirect illuminance\n");
}
