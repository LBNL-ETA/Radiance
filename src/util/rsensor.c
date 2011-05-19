#ifndef lint
static const char RCSid[] = "$Id: rsensor.c,v 2.10 2011/05/19 16:09:38 greg Exp $";
#endif

/*
 * Compute sensor signal based on spatial sensitivity.
 *
 *	Created Feb 2008 for Architectural Energy Corp.
 */

#include "ray.h"
#include "source.h"
#include "view.h"
#include "random.h"

#define DEGREE		(PI/180.)

#define	MAXNT		180	/* maximum number of theta divisions */
#define MAXNP		360	/* maximum number of phi divisions */

extern char	*progname;	/* global argv[0] */
extern int	nowarn;         /* don't report warnings? */

				/* current sensor's perspective */
VIEW		ourview =  {VT_ANG,{0.,0.,0.},{0.,0.,1.},{1.,0.,0.},
                                1.,180.,180.,0.,0.,0.,0.,
                                {0.,0.,0.},{0.,0.,0.},0.,0.};

unsigned long	nsamps = 10000;	/* desired number of initial samples */
unsigned long	nssamps = 9000;	/* number of super-samples */
int		ndsamps = 32;	/* number of direct samples */
int		nprocs = 1;     /* number of rendering processes */

float		*sensor = NULL;	/* current sensor data */
int		sntp[2];	/* number of sensor theta and phi angles */
float		maxtheta;	/* maximum theta value for this sensor */
float		tvals[MAXNT+1];	/* theta prob. values (1-D table of 1-cos(t)) */
float		*pvals = NULL;	/* phi prob. values (2-D table in radians) */
int		ntheta = 0;	/* polar angle divisions */
int		nphi = 0;	/* azimuthal angle divisions */
double		gscale = 1.;	/* global scaling value */

#define	s_theta(t)	sensor[(t+1)*(sntp[1]+1)]
#define s_phi(p)	sensor[(p)+1]
#define s_val(t,p)	sensor[(p)+1+(t+1)*(sntp[1]+1)]

static void	comp_sensor(char *sfile);

static void
over_options()			/* overriding options */
{
	directvis = (ndsamps <= 0);
	do_irrad = 0;
}

static void
print_defaults()		/* print out default parameters */
{
	over_options();
	printf("-n %-9d\t\t\t# number of processes\n", nprocs);
	printf("-rd %-9ld\t\t\t# ray directions\n", nsamps);
	/* printf("-rs %-9ld\t\t\t# ray super-samples\n", nssamps); */
	printf("-dn %-9d\t\t\t# direct number of samples\n", ndsamps);
	printf("-vp %f %f %f\t# view point\n",
			ourview.vp[0], ourview.vp[1], ourview.vp[2]);
	printf("-vd %f %f %f\t# view direction\n",
			ourview.vdir[0], ourview.vdir[1], ourview.vdir[2]);
	printf("-vu %f %f %f\t# view up\n",
			ourview.vup[0], ourview.vup[1], ourview.vup[2]);
	printf("-vo %f\t\t\t# view fore clipping distance\n", ourview.vfore);
	print_rdefaults();
}


void
quit(ec)			/* make sure exit is called */
int	ec;
{
	if (ray_pnprocs > 0)	/* close children if any */
		ray_pclose(0);		
	exit(ec);
}


int
main(
	int	argc,
	char	*argv[]
)
{
	int	doheader = 1;
	int	optwarn = 0;
	int	i, rval;

	progname = argv[0];
				/* set up rendering defaults */
	rand_samp = 1;
	dstrsrc = 0.65;
	srcsizerat = 0.1;
	directrelay = 3;
	ambounce = 1;
	maxdepth = -10;
				/* get options from command line */
	for (i = 1; i < argc; i++) {
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", argv[i]);
			error(SYSTEM, errmsg);
		}
		if (argv[i][0] != '-') {
			if (i >= argc-1)
				break;		/* final octree argument */
			if (!ray_pnprocs) {
				over_options();
				if (doheader) {	/* print header */
					newheader("RADIANCE", stdout);
					printargs(argc, argv, stdout);
					fputformat("ascii", stdout);
					putchar('\n');
				}
						/* start process(es) */
				if (strcmp(argv[argc-1], "."))
					ray_pinit(argv[argc-1], nprocs);
			}
			comp_sensor(argv[i]);	/* process a sensor file */
			continue;
		}
		if (argv[i][1] == 'r') {	/* sampling options */
			if (argv[i][2] == 'd')
				nsamps = atol(argv[++i]);
			else if (argv[i][2] == 's')
				nssamps = atol(argv[++i]);
			else {
				sprintf(errmsg, "bad option at '%s'", argv[i]);
				error(USER, errmsg);
			}
			continue;
		}
						/* direct component samples */
		if (argv[i][1] == 'd' && argv[i][2] == 'n') {
			ndsamps = atoi(argv[++i]);
			continue;
		}
		if (argv[i][1] == 'v') {	/* next sensor view */
			if (argv[i][2] == 'f') {
				rval = viewfile(argv[++i], &ourview, NULL);
				if (rval < 0) {
					sprintf(errmsg,
					"cannot open view file \"%s\"",
							argv[i]);
					error(SYSTEM, errmsg);
				} else if (rval == 0) {
					sprintf(errmsg,
						"bad view file \"%s\"",
							argv[i]);
					error(USER, errmsg);
				}
				continue;
			}
			rval = getviewopt(&ourview, argc-i, argv+i);
			if (rval >= 0) {
				i += rval;
				continue;
			}
			sprintf(errmsg, "bad view option at '%s'", argv[i]);
			error(USER, errmsg);
		}
		if (!strcmp(argv[i], "-w")) {	/* toggle warnings */
			nowarn = !nowarn;
			continue;
		}
		if (ray_pnprocs) {
			if (!optwarn++)
				error(WARNING,
			"rendering options should appear before first sensor");
		} else if (!strcmp(argv[i], "-defaults")) {
			print_defaults();
			return(0);
		}
		if (argv[i][1] == 'h') {	/* header toggle */
			doheader = !doheader;
			continue;
		}
		if (!strcmp(argv[i], "-n")) {	/* number of processes */
			nprocs = atoi(argv[++i]);
			if (nprocs <= 0)
				error(USER, "illegal number of processes");
			continue;
		}
		rval = getrenderopt(argc-i, argv+i);
		if (rval < 0) {
			sprintf(errmsg, "bad render option at '%s'", argv[i]);
			error(USER, errmsg);
		}
		i += rval;
	}
	if (sensor == NULL)
		error(USER, i<argc ? "missing sensor file" : "missing octree");
	quit(0);
}

/* Load sensor sensitivities (first row and column are angles) */
static float *
load_sensor(
	int	ntp[2],
	char	*sfile
)
{
	int	warnedneg;
	char	linebuf[8192];
	int	nelem = 1000;
	float	*sarr = (float *)malloc(sizeof(float)*nelem);
	FILE	*fp;
	char	*cp;
	int	i;
	
	fp = frlibopen(sfile);
	if (fp == NULL) {
		sprintf(errmsg, "cannot open sensor file '%s'", sfile);
		error(SYSTEM, errmsg);
	}
	fgets(linebuf, sizeof(linebuf), fp);
	if (!strncmp(linebuf, "Elevation ", 10))
		fgets(linebuf, sizeof(linebuf), fp);
						/* get phi values */
	sarr[0] = .0f;
	if (strncmp(linebuf, "degrees", 7)) {
		sprintf(errmsg, "Missing 'degrees' in sensor file '%s'", sfile);
		error(USER, errmsg);
	}
	cp = sskip(linebuf);
	ntp[1] = 0;
	for ( ; ; ) {
		sarr[ntp[1]+1] = atof(cp);
		cp = fskip(cp);
		if (cp == NULL)
			break;
		if (ntp[1] > 1 && sarr[ntp[1]+1] <= sarr[ntp[1]]) {
			sprintf(errmsg,
		"Phi values not monotinically increasing in sensor file '%s'",
					sfile);
			error(USER, errmsg);
		}
		++ntp[1];
	}
	warnedneg = 0;
	ntp[0] = 0;				/* get thetas + data */
	while (fgets(linebuf, sizeof(linebuf), fp) != NULL) {
		++ntp[0];
		if ((ntp[0]+1)*(ntp[1]+1) > nelem) {
			nelem += (nelem>>2) + ntp[1];
			sarr = (float *)realloc((void *)sarr,
					sizeof(float)*nelem);
			if (sarr == NULL)
				error(SYSTEM, "out of memory in load_sensor()");
		}
		cp = linebuf;
		i = ntp[0]*(ntp[1]+1);
		for ( ; ; ) {
			sarr[i] = atof(cp);
			cp = fskip(cp);
			if (cp == NULL)
				break;
			if (i && sarr[i] < .0) {
				if (!warnedneg++) {
					sprintf(errmsg,
		"Negative value(s) in sensor file '%s' (ignored)\n", sfile);
					error(WARNING, errmsg);
				}
				sarr[i] = .0;
			}
			++i;
		}
		if (i == ntp[0]*(ntp[1]+1))
			break;
		if (ntp[0] > 1 && sarr[ntp[0]*(ntp[1]+1)] <=
					sarr[(ntp[0]-1)*(ntp[1]+1)]) {
			sprintf(errmsg,
		"Theta values not monotinically increasing in sensor file '%s'",
					sfile);
			error(USER, errmsg);
		}
		if (i != (ntp[0]+1)*(ntp[1]+1)) {
			sprintf(errmsg,
			"bad column count near line %d in sensor file '%s'",
					ntp[0]+1, sfile);
			error(USER, errmsg);
		}
	}
	nelem = i;
	fclose(fp);
	errmsg[0] = '\0';			/* sanity checks */
	if (ntp[0] <= 0)
		sprintf(errmsg, "no data in sensor file '%s'", sfile);
	else if (fabs(sarr[ntp[1]+1]) > FTINY)
		sprintf(errmsg, "minimum theta must be 0 in sensor file '%s'",
				sfile);
	else if (fabs(sarr[1]) > FTINY)
		sprintf(errmsg, "minimum phi must be 0 in sensor file '%s'",
				sfile);
	else if (sarr[ntp[1]] <= FTINY)
		sprintf(errmsg,
			"maximum phi must be positive in sensor file '%s'",
				sfile);
	else if (sarr[ntp[0]*(ntp[1]+1)] <= FTINY)
		sprintf(errmsg,
			"maximum theta must be positive in sensor file '%s'",
				sfile);
	if (errmsg[0])
		error(USER, errmsg);
	return((float *)realloc((void *)sarr, sizeof(float)*nelem));
}

/* Initialize probability table */
static void
init_ptable(
	char	*sfile
)
{
	int	samptot = nsamps;
	float	*rowp, *rowp1;
	double	rowsum[MAXNT], rowomega[MAXNT];
	double	thdiv[MAXNT+1], phdiv[MAXNP+1];
	double	tsize, psize;
	double	prob, frac, frac1;
	int	i, j, t, p;
					/* free old table */
	if (sensor != NULL)
		free((void *)sensor);
	if (pvals != NULL)
		free((void *)pvals);
	if (sfile == NULL || !*sfile) {
		sensor = NULL;
		sntp[0] = sntp[1] = 0;
		pvals = NULL;
		ntheta = nphi = 0;
		return;
	}
					/* load sensor table */
	sensor = load_sensor(sntp, sfile);
	if (sntp[0] > MAXNT) {
		sprintf(errmsg, "Too many theta rows in sensor file '%s'",
				sfile);
		error(INTERNAL, errmsg);
	}
	if (sntp[1] > MAXNP) {
		sprintf(errmsg, "Too many phi columns in sensor file '%s'",
				sfile);
		error(INTERNAL, errmsg);
	}
					/* compute boundary angles */
	maxtheta = DEGREE*(1.5f*s_theta(sntp[0]-1) - 0.5f*s_theta(sntp[0]-2));
	if (maxtheta > PI)
		maxtheta = PI;
	thdiv[0] = .0;
	for (t = 1; t < sntp[0]; t++)
		thdiv[t] = DEGREE/2.*(s_theta(t-1) + s_theta(t));
	thdiv[sntp[0]] = maxtheta;
	phdiv[0] = DEGREE*(1.5f*s_phi(0) - 0.5f*s_phi(1));
	for (p = 1; p < sntp[1]; p++)
		phdiv[p] = DEGREE/2.*(s_phi(p-1) + s_phi(p));
	phdiv[sntp[1]] = DEGREE*(1.5f*s_phi(sntp[1]-1) - 0.5f*s_phi(sntp[1]-2));
					/* size our table */
	tsize = 1. - cos(maxtheta);
	psize = PI*tsize/maxtheta;
	if (sntp[0]*sntp[1] < samptot)	/* don't overdo resolution */
		samptot = sntp[0]*sntp[1];
	ntheta = (int)(sqrt((double)samptot*tsize/psize) + 0.5);
	if (ntheta > MAXNT)
		ntheta = MAXNT;
	nphi = samptot/ntheta;
	pvals = (float *)malloc(sizeof(float)*(ntheta+1)*(nphi+1));
	if (pvals == NULL)
		error(SYSTEM, "out of memory in init_ptable()");
	gscale = .0;			/* compute our inverse table */
	for (i = 0; i < sntp[0]; i++) {
		rowp = &s_val(i,0);
		rowsum[i] = 1e-20;
		for (j = 0; j < sntp[1]; j++)
			rowsum[i] += *rowp++;
		rowomega[i] = cos(thdiv[i]) - cos(thdiv[i+1]);
		rowomega[i] *= 2.*PI / (double)sntp[1];
		gscale += rowsum[i] * rowomega[i];
	}
	if (gscale <= FTINY) {
		sprintf(errmsg, "Sensor values sum to zero in file '%s'", sfile);
		error(USER, errmsg);
	}
	for (i = 0; i < ntheta; i++) {
		prob = (double)i / (double)ntheta;
		for (t = 0; t < sntp[0]; t++)
			if ((prob -= rowsum[t]*rowomega[t]/gscale) <= .0)
				break;
		if (t >= sntp[0])
			error(INTERNAL, "code error 1 in init_ptable()");
		frac = 1. + prob/(rowsum[t]*rowomega[t]/gscale);
		tvals[i] = 1. - ( (1.-frac)*cos(thdiv[t]) +
						frac*cos(thdiv[t+1]) );
				/* offset b/c sensor values are centered */
		if (!t || (t < sntp[0]-1) & (frac >= 0.5))
			frac -= 0.5;
		else {
			frac += 0.5;
			--t;
		}
		pvals[i*(nphi+1)] = phdiv[0];
		for (j = 1; j < nphi; j++) {
			prob = (double)j / (double)nphi;
			rowp = &s_val(t,0);
			rowp1 = &s_val(t+1,0);
			for (p = 0; p < sntp[1]; p++)
				if ((prob -= (1.-frac)*rowp[p]/rowsum[t] +
					    frac*rowp1[p]/rowsum[t+1]) <= .0)
					break;
			if (p >= sntp[1]) {	/* should never happen? */
				p = sntp[1] - 1;
				prob = .5;
			}
			frac1 = 1. + prob/((1.-frac)*rowp[p]/rowsum[t]
					+ frac*rowp1[p]/rowsum[t+1]);
			pvals[i*(nphi+1) + j] = (1.-frac1)*phdiv[p] +
						frac1*phdiv[p+1];
		}
		pvals[i*(nphi+1) + nphi] = phdiv[sntp[1]];
	}
	tvals[0] = .0f;
	tvals[ntheta] = (float)tsize;
}

/* Get normalized direction from random variables in [0,1) range */
static void
get_direc(
	FVECT dvec,
	double	x,
	double	y
)
{
	double	xfrac = x*ntheta;
	int	tndx = (int)xfrac;
	double	yfrac = y*nphi;
	int	pndx = (int)yfrac;
	double	rad, phi;
	FVECT	dv;
	int	i;

	xfrac -= (double)tndx;
	yfrac -= (double)pndx;
	pndx += tndx*(nphi+1);
	
	dv[2] = 1. - ((1.-xfrac)*tvals[tndx] + xfrac*tvals[tndx+1]);
	rad = sqrt(1. - dv[2]*dv[2]);
	phi = (1.-yfrac)*pvals[pndx] + yfrac*pvals[pndx+1];
	dv[0] = -rad*sin(phi);
	dv[1] = rad*cos(phi);
	for (i = 3; i--; )
		dvec[i] = dv[0]*ourview.hvec[i] +
				dv[1]*ourview.vvec[i] +
				dv[2]*ourview.vdir[i] ;
}

/* Get sensor value in the specified direction (normalized) */
static float
sens_val(
	FVECT	dvec
)
{
	FVECT	dv;
	float	theta, phi;
	int	t, p;
	
	dv[2] = DOT(dvec, ourview.vdir);
	theta = (float)((1./DEGREE) * acos(dv[2]));
	if (theta >= maxtheta)
		return(.0f);
	dv[0] = DOT(dvec, ourview.hvec);
	dv[1] = DOT(dvec, ourview.vvec);
	phi = (float)((1./DEGREE) * atan2(-dv[0], dv[1]));
	while (phi < .0f) phi += 360.f;
	t = (int)(theta/maxtheta * sntp[0]);
	p = (int)(phi*(1./360.) * sntp[1]);
			/* hack for non-uniform sensor grid */
	while (t+1 < sntp[0] && theta >= s_theta(t+1))
		++t;
	while (t-1 >= 0 && theta <= s_theta(t-1))
		--t;
	while (p+1 < sntp[1] && phi >= s_phi(p+1))
		++p;
	while (p-1 >= 0 && phi <= s_phi(p-1))
		--p;
	return(s_val(t,p));
}

/* Print origin and direction */
static void
print_ray(
	FVECT rorg,
	FVECT rdir
)
{
	printf("%.6g %.6g %.6g %.8f %.8f %.8f\n",
			rorg[0], rorg[1], rorg[2],
			rdir[0], rdir[1], rdir[2]);
}

/* Compute sensor output */
static void
comp_sensor(
	char *sfile
)
{
	int	ndirs = dstrsrc > FTINY ? ndsamps :
				ndsamps > 0 ? 1 : 0;
	char	*err;
	int	nt, np;
	COLOR	vsum;
	RAY	rr;
	double	sf;
	int	i, j;
						/* set view */
	ourview.type = VT_ANG;
	ourview.horiz = ourview.vert = 180.;
	ourview.hoff = ourview.voff = .0;
	err = setview(&ourview);
	if (err != NULL)
		error(USER, err);
						/* assign probability table */
	init_ptable(sfile);
						/* stratified MC sampling */
	setcolor(vsum, .0f, .0f, .0f);
	nt = (int)(sqrt((double)nsamps*ntheta/nphi) + .5);
	np = nsamps/nt;
	sf = gscale/nsamps;
	for (i = 0; i < nt; i++)
		for (j = 0; j < np; j++) {
			VCOPY(rr.rorg, ourview.vp);
			get_direc(rr.rdir, (i+frandom())/nt, (j+frandom())/np);
			if (ourview.vfore > FTINY)
				VSUM(rr.rorg, rr.rorg, rr.rdir, ourview.vfore);
			if (!ray_pnprocs) {
				print_ray(rr.rorg, rr.rdir);
				continue;
			}
			rr.rmax = .0;
			rayorigin(&rr, PRIMARY, NULL, NULL);
			scalecolor(rr.rcoef, sf);
			if (ray_pqueue(&rr) == 1)
				addcolor(vsum, rr.rcol);
		}
						/* remaining rays pure MC */
	for (i = nsamps - nt*np; i-- > 0; ) {
		VCOPY(rr.rorg, ourview.vp);
		get_direc(rr.rdir, frandom(), frandom());
		if (ourview.vfore > FTINY)
			VSUM(rr.rorg, rr.rorg, rr.rdir, ourview.vfore);
		if (!ray_pnprocs) {
			print_ray(rr.rorg, rr.rdir);
			continue;
		}
		rr.rmax = .0;
		rayorigin(&rr, PRIMARY, NULL, NULL);
		scalecolor(rr.rcoef, sf);
		if (ray_pqueue(&rr) == 1)
			addcolor(vsum, rr.rcol);
	}
	if (!ray_pnprocs)			/* just printing rays */
		return;
						/* scale partial result */
	scalecolor(vsum, sf);
						/* add direct component */
	for (i = ndirs; i-- > 0; ) {
		SRCINDEX	si;
		initsrcindex(&si);
		while (srcray(&rr, NULL, &si)) {
			sf = sens_val(rr.rdir);
			if (sf <= FTINY)
				continue;
			sf *= si.dom/ndirs;
			scalecolor(rr.rcoef, sf);
			if (ray_pqueue(&rr) == 1) {
				multcolor(rr.rcol, rr.rcoef);
				addcolor(vsum, rr.rcol);
			}
		}
	}
						/* finish our calculation */
	while (ray_presult(&rr, 0) > 0) {
		multcolor(rr.rcol, rr.rcoef);
		addcolor(vsum, rr.rcol);
	}
						/* print our result */
	printf("%.4e %.4e %.4e\n", colval(vsum,RED),
				colval(vsum,GRN), colval(vsum,BLU));
}
