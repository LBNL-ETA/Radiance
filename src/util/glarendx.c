#ifndef lint
static const char	RCSid[] = "$Id: glarendx.c,v 2.10 2004/01/02 12:48:36 schorsch Exp $";
#endif
/*
 * Compute Glare Index given by program name or -t option:
 *
 *	dgi -		Daylight Glare Index
 *	brs_gi -	Building Research Station Glare Index (Petherbridge
 *							       & Hopkinson)
 *	ugr -		Unified Glare Rating System (Fischer)
 *	guth_dgr -	Guth discomfort glare rating
 *	guth_vcp -	Guth visual comfort probability
 *	cie_cgi -	CIE Glare Index (1983, due to Einhorn)
 *	vert_dir -	Direct vertical illuminance
 *	vert_ind -	Indirect vertical illuminance (from input)
 *	vert_ill -	Total vertical illuminance
 *
 *		12 April 1991	Greg Ward	EPFL
 *		19 April 1993   R. Compagnon    EPFL (added dgi, brs_gi, ugr)
 */

#include <string.h>

#include "standard.h"
#include "view.h"


struct glare_src {
	FVECT	dir;		/* source direction */
	double	dom;		/* solid angle */
	double	lum;		/* average luminance */
	struct glare_src	*next;
} *all_srcs = NULL;

struct glare_dir {
	double	ang;		/* angle (in radians) */
	double	indirect;	/* indirect illuminance */
	struct glare_dir	*next;
} *all_dirs = NULL;

typedef double gdfun(struct glare_dir *gd);

static gdfun dgi;
static gdfun brs_gi;
static gdfun ugr;
static gdfun guth_dgr;
static gdfun guth_vcp;
static gdfun cie_cgi;
static gdfun direct;
static gdfun indirect;
static gdfun total;

static gethfunc headline;
static void init(void);
static void read_input(void);
static void print_values(gdfun *func);
static double posindex(FVECT sd, FVECT vd, FVECT vu);


struct named_func {
	char	*name;
	gdfun	*func;
	char	*descrip;
} all_funcs[] = {
	{"dgi", dgi, "Daylight Glare Index"},
	{"brs_gi", brs_gi, "BRS Glare Index"},
	{"ugr", ugr, "Unified Glare Rating"},
	{"guth_vcp", guth_vcp, "Guth Visual Comfort Probability"},
	{"cie_cgi", cie_cgi, "CIE Glare Index (Einhorn)"},
	{"guth_dgr", guth_dgr, "Guth Disability Glare Rating"},
	{"vert_dir", direct, "Direct Vertical Illuminance"},
	{"vert_ill", total, "Total Vertical Illuminance"},
	{"vert_ind", indirect, "Indirect Vertical Illuminance"},
	{NULL}
};

#define newp(type)	(type *)malloc(sizeof(type))

char	*progname;
int	print_header = 1;

VIEW	midview = STDVIEW;

int	wrongformat = 0;


int
main(
	int	argc,
	char	*argv[]
)
{
	struct named_func	*funp;
	char	*progtail;
	int	i;
					/* get program name */
	progname = argv[0];
	progtail = strrchr(progname, '/');	/* final component */
	if (progtail == NULL)
		progtail = progname;
	else
		progtail++;
					/* get options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 't':
			progtail = argv[++i];
			break;
		case 'h':
			print_header = 0;
			break;
		default:
			goto userr;
		}
	if (i < argc-1)
		goto userr;
	if (i == argc-1)		/* open file */
		if (freopen(argv[i], "r", stdin) == NULL) {
			perror(argv[i]);
			exit(1);
		}
					/* find and run calculation */
	for (funp = all_funcs; funp->name != NULL; funp++)
		if (!strcmp(funp->name, progtail)) {
			init();
			read_input();
			if (print_header) {
				printargs(i, argv, stdout);
				putchar('\n');
			}
			print_values(funp->func);
			exit(0);		/* we're done */
		}
					/* invalid function */
userr:
	fprintf(stderr, "Usage: %s -t type [-h] [input]\n", progname);
	fprintf(stderr, "\twhere 'type' is one of the following:\n");
	for (funp = all_funcs; funp->name != NULL; funp++)
		fprintf(stderr, "\t%12s\t%s\n", funp->name, funp->descrip);
	exit(1);
}


static int
headline(			/* get line from header */
	char	*s,
	void	*p
)
{
	char	fmt[32];

	if (print_header)		/* copy to output */
		fputs(s, stdout);
	if (isview(s))
		sscanview(&midview, s);
	else if (isformat(s)) {
		formatval(fmt, s);
		wrongformat = strcmp(fmt, "ascii");
	}
	return(0);
}


static void
init(void)				/* initialize calculation */
{
					/* read header */
	getheader(stdin, headline, NULL);
	if (wrongformat) {
		fprintf(stderr, "%s: bad input format\n", progname);
		exit(1);
	}
					/* set view */
	if (setview(&midview) != NULL) {
		fprintf(stderr, "%s: bad view information in input\n",
				progname);
		exit(1);
	}
}


static void
read_input(void)			/* read glare sources from stdin */
{
#define	S_SEARCH	0
#define S_SOURCE	1
#define S_DIREC		2
	int	state = S_SEARCH;
	char	buf[128];
	register struct glare_src	*gs;
	register struct glare_dir	*gd;

	while (fgets(buf, sizeof(buf), stdin) != NULL)
		switch (state) {
		case S_SEARCH:
			if (!strcmp(buf, "BEGIN glare source\n"))
				state = S_SOURCE;
			else if (!strcmp(buf, "BEGIN indirect illuminance\n"))
				state = S_DIREC;
			break;
		case S_SOURCE:
			if (!strncmp(buf, "END", 3)) {
				state = S_SEARCH;
				break;
			}
			if ((gs = newp(struct glare_src)) == NULL)
				goto memerr;
			if (sscanf(buf, "%lf %lf %lf %lf %lf",
					&gs->dir[0], &gs->dir[1], &gs->dir[2],
					&gs->dom, &gs->lum) != 5)
				goto readerr;
			normalize(gs->dir);
			gs->next = all_srcs;
			all_srcs = gs;
			break;
		case S_DIREC:
			if (!strncmp(buf, "END", 3)) {
				state = S_SEARCH;
				break;
			}
			if ((gd = newp(struct glare_dir)) == NULL)
				goto memerr;
			if (sscanf(buf, "%lf %lf",
					&gd->ang, &gd->indirect) != 2)
				goto readerr;
			gd->ang *= PI/180.0;	/* convert to radians */
			gd->next = all_dirs;
			all_dirs = gd;
			break;
		}
	return;
memerr:
	perror(progname);
	exit(1);
readerr:
	fprintf(stderr, "%s: read error on input\n", progname);
	exit(1);
#undef S_SEARCH
#undef S_SOURCE
#undef S_DIREC
}


static void
print_values(		/* print out calculations */
	gdfun *func
)
{
	register struct glare_dir	*gd;

	for (gd = all_dirs; gd != NULL; gd = gd->next)
		printf("%f\t%f\n", gd->ang*(180.0/PI), (*func)(gd));
}


static double
direct(			/* compute direct vertical illuminance */
	struct glare_dir	*gd
)
{
	FVECT	mydir;
	double	d, dval;
	register struct glare_src	*gs;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	dval = 0.0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {
		d = DOT(mydir,gs->dir);
		if (d > FTINY)
			dval += d * gs->dom * gs->lum;
	}
	return(dval);
}


static double
indirect(			/* return indirect vertical illuminance */
	struct glare_dir	*gd
)
{
	return(gd->indirect);
}


static double
total(			/* return total vertical illuminance */
	struct glare_dir	*gd
)
{
	return(direct(gd)+gd->indirect);
}


/*
 * posindex -	compute glare position index from:
 *
 *	Source Direction
 *	View Direction
 *	View Up Direction
 *
 * All vectors are assumed to be normalized.
 * This function is an implementation of the method proposed by
 * Robert Levin in his 1975 JIES article.
 * This calculation presumes the view direction and up vectors perpendicular.
 * We return a value less than zero for improper positions.
 */

static double
posindex(			/* compute position index */
	FVECT	sd,
	FVECT	vd,
	FVECT	vu
)
{
	double	sigma, tau;
	double	d;

	d = DOT(sd,vd);
	if (d <= 0.0)
		return(-1.0);
	if (d >= 1.0)
		return(1.0);
	sigma = acos(d) * (180./PI);
	d = fabs(DOT(sd,vu)/sqrt(1.0-d*d));
	if (d >= 1.0)
		tau = 0.0;
	else
		tau = acos(d) * (180./PI);
	return( exp( sigma*( (35.2 - tau*.31889 - 1.22*exp(-.22222*tau))*1e-3
			+ sigma*(21. + tau*(.26667 + tau*-.002963))*1e-5 )
		) );
}


static double
dgi(		/* compute Daylight Glare Index */
	struct glare_dir	*gd
)
{
	register struct glare_src	*gs;
	FVECT	mydir,testdir[7],vhor;
	double	r,posn,omega,p[7],sum;
	int	i,n;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	sum = 0.0; n = 0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {

		/* compute 1/p^2 weighted solid angle of the source */
		r = sqrt(1 - pow(1.-gs->dom/2./PI,2.));
		fcross(vhor,gs->dir,midview.vup);
		normalize(vhor);
		VCOPY(testdir[0],gs->dir);
		fvsum(testdir[1],gs->dir,vhor,r);
		fvsum(testdir[2],gs->dir,vhor,0.5*r);
		fvsum(testdir[5],testdir[2],midview.vup,-0.866*r);
		fvsum(testdir[2],testdir[2],midview.vup,0.866*r);
		fvsum(testdir[3],gs->dir,vhor,-r);
		fvsum(testdir[4],gs->dir,vhor,-0.5*r);
		fvsum(testdir[6],testdir[4],midview.vup,0.866*r);
		fvsum(testdir[4],testdir[4],midview.vup,-0.866*r);
		for (i = 0; i < 7; i++) {
			normalize(testdir[i]);
			posn = posindex(testdir[i],mydir,midview.vup);
			if (posn <= FTINY)
				p[i] = 0.0;
			else
				p[i] = 1./(posn*posn);
		}
		r = 1-gs->dom/2./PI;
		omega = gs->dom*p[0];
		omega += (r*PI*(1+1/r/r)-2*PI)*(-p[0]+(p[1]+p[2])*0.5);
		omega += (2*PI-r*PI*(1+1/r/r))*(-p[0]-0.1667*(p[1]+p[3])
			  +0.3334*(p[2]+p[4]+p[5]+p[6]));

		sum += pow(gs->lum,1.6) * pow(omega,0.8) /
		       (gd->indirect/PI + 0.07*sqrt(gs->dom)*gs->lum);
		n++;
	}
	if (n == 0)
		return(0.0);
	return( 10*log10(0.478*sum) );
}


static double
brs_gi(		/* compute BRS Glare Index */
	struct glare_dir	*gd
)
{
	register struct glare_src	*gs;
	FVECT	mydir;
	double	p;
	double	sum;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	sum = 0.0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {
		p = posindex(gs->dir, mydir, midview.vup);
		if (p <= FTINY)
			continue;
		sum += pow(gs->lum/p,1.6) * pow(gs->dom,0.8);
	}
	if (sum <= FTINY)
		return(0.0);
	sum /= gd->indirect/PI;
	return(10*log10(0.478*sum));
}


static double
guth_dgr(		/* compute Guth discomfort glare rating */
	struct glare_dir	*gd
)
{
#define q(w)	(20.4*w+1.52*pow(w,.2)-.075)
	register struct glare_src	*gs;
	FVECT	mydir;
	double	p;
	double	sum;
	double	wtot, brsum;
	int	n;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	sum = wtot = brsum = 0.0; n = 0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {
		p = posindex(gs->dir, mydir, midview.vup);
		if (p <= FTINY)
			continue;
		sum += gs->lum * q(gs->dom) / p;
		brsum += gs->lum * gs->dom;
		wtot += gs->dom;
		n++;
	}
	if (n == 0)
		return(0.0);
	return( pow(.5*sum/pow((brsum+(5.-wtot)*gd->indirect/PI)/5.,.44),
			pow((double)n, -.0914) ) );
#undef q
}


#ifndef M_SQRT2
#define	M_SQRT2	1.41421356237309504880
#endif

#define norm_integral(z)	(1.-.5*erfc((z)/M_SQRT2))


static double
guth_vcp(		/* compute Guth visual comfort probability */
	struct glare_dir	*gd
)
{
	extern double	erfc();
	double	dgr;

	dgr = guth_dgr(gd);
	if (dgr <= FTINY)
		return(100.0);
	return(100.*norm_integral(6.374-1.3227*log(dgr)));
}


static double
cie_cgi(		/* compute CIE Glare Index */
	struct glare_dir	*gd
)
{
	register struct glare_src	*gs;
	FVECT	mydir;
	double	dillum;
	double	p;
	double	sum;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	sum = 0.0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {
		p = posindex(gs->dir, mydir, midview.vup);
		if (p <= FTINY)
			continue;
		sum += gs->lum*gs->lum * gs->dom / (p*p);
	}
	if (sum <= FTINY)
		return(0.0);
	dillum = direct(gd);
	return(8.*log10(2.*sum*(1.+dillum/500.)/(dillum+gd->indirect)));
}


static double
ugr(		/* compute Unified Glare Rating */
	struct glare_dir	*gd
)
{
	register struct glare_src	*gs;
	FVECT	mydir;
	double	p;
	double	sum;

	spinvector(mydir, midview.vdir, midview.vup, gd->ang);
	sum = 0.0;
	for (gs = all_srcs; gs != NULL; gs = gs->next) {
		p = posindex(gs->dir, mydir, midview.vup);
		if (p <= FTINY)
			continue;
		sum += gs->lum*gs->lum * gs->dom / (p*p);
	}
	if (sum <= FTINY)
		return(0.0);
	return(8.*log10(0.25*sum*PI/gd->indirect));
}
