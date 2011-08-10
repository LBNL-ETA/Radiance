#ifndef lint
static const char RCSid[] = "$Id: mksource.c,v 2.6 2011/08/10 22:28:45 greg Exp $";
#endif
/*
 * Generate distant sources corresponding to the given environment map
 */

#include "ray.h"
#include "random.h"
#include "resolu.h"

#define NTRUNKBR	4		/* number of branches at trunk */
#define NTRUNKVERT	4		/* number of vertices at trunk */
#define DEF_NSAMPS	262144L		/* default # sphere samples */
#define DEF_MAXANG	15.		/* maximum source angle (deg.) */

/* Data structure for geodesic samples */

typedef struct tritree {
	int32		gdv[3];		/* spherical triangle vertex direc. */
	int32		sd;		/* sample direction if leaf */
	struct tritree	*kid;		/* 4 children if branch node */
	COLR		val;		/* sampled color value */
} TRITREE;

typedef struct lostlight {
	struct lostlight	*next;	/* next in list */
	int32		sd;		/* lost source direction */
	COLOR		intens;		/* output times solid angle */
} LOSTLIGHT;

extern char	*progname;

FVECT	scene_cent;		/* center of octree cube */
RREAL	scene_rad;		/* radius to get outside cube from center */

const COLR	blkclr = BLKCOLR;

#define isleaf(node)	((node)->kid == NULL)

/* Compute signum of signed volume for three vectors */
int
vol_sign(const FVECT v1, const FVECT v2, const FVECT v3)
{
	double	vol;
	
	vol  = v1[0]*(v2[1]*v3[2] - v2[2]*v3[1]);
	vol += v1[1]*(v2[2]*v3[0] - v2[0]*v3[2]);
	vol += v1[2]*(v2[0]*v3[1] - v2[1]*v3[0]);
	if (vol > .0)
		return(1);
	if (vol < .0)
		return(-1);
	return(0);
}

/* Is the given direction contained within the specified spherical triangle? */
int
intriv(const int32 trid[3], const FVECT sdir)
{
	int	sv[3];
	FVECT	tri[3];

	decodedir(tri[0], trid[0]);
	decodedir(tri[1], trid[1]);
	decodedir(tri[2], trid[2]);
	sv[0] = vol_sign(sdir, tri[0], tri[1]);
	sv[1] = vol_sign(sdir, tri[1], tri[2]);
	sv[2] = vol_sign(sdir, tri[2], tri[0]);
	if ((sv[0] == sv[1]) & (sv[1] == sv[2]))
		return(1);
	return(!sv[0] | !sv[1] | !sv[2]);
}

/* Find leaf containing the given sample direction */
TRITREE *
findleaf(TRITREE *node, const FVECT sdir)
{
	int	i;

	if (isleaf(node))
		return(intriv(node->gdv,sdir) ? node : (TRITREE *)NULL);
	for (i = 0; i < 4; i++) {
		TRITREE	*chknode = &node->kid[i];
		if (intriv(chknode->gdv, sdir))
			return(isleaf(chknode) ? chknode :
					findleaf(chknode, sdir));
	}
	return(NULL);
}

/* Initialize leaf with random sample inside the given spherical triangle */
void
leafsample(TRITREE *leaf)
{
	RAY	myray;
	RREAL	wt[3];
	FVECT	sdir;
	int	i, j;
						/* random point on triangle */
	i = random() % 3;
	wt[i] = frandom();
	j = random() & 1;
	wt[(i+2-j)%3] = 1. - wt[i] - 
			(wt[(i+1+j)%3] = (1.-wt[i])*frandom());
	sdir[0] = sdir[1] = sdir[2] = .0;
	for (i = 0; i < 3; i++) {
		FVECT	vt;
		decodedir(vt, leaf->gdv[i]);
		VSUM(sdir, sdir, vt, wt[i]);
	}
	normalize(sdir);			/* record sample direction */
	leaf->sd = encodedir(sdir);
						/* evaluate at inf. */
	VSUM(myray.rorg, scene_cent, sdir, scene_rad);
	VCOPY(myray.rdir, sdir);
	myray.rmax = 0.;
	ray_trace(&myray);
	setcolr(leaf->val, colval(myray.rcol,RED),
			colval(myray.rcol,GRN),
			colval(myray.rcol,BLU));
}

/* Initialize a branch node contained in the given spherical triangle */
void
subdivide(TRITREE *branch, const int32 dv[3])
{
	FVECT	dvv[3], sdv[3];
	int32	sd[3];
	int	i;
	
	for (i = 0; i < 3; i++) {	/* copy spherical triangle */
		branch->gdv[i] = dv[i];
		decodedir(dvv[i], dv[i]);
	}
	for (i = 0; i < 3; i++) {	/* create new vertices */
		int	j = (i+1)%3;
		VADD(sdv[i], dvv[i], dvv[j]);
		normalize(sdv[i]);
		sd[i] = encodedir(sdv[i]);
	}
					/* allocate leaves */
	branch->kid = (TRITREE *)calloc(4, sizeof(TRITREE));
	if (branch->kid == NULL)
		error(SYSTEM, "out of memory in subdivide()");
					/* assign subtriangle directions */
	branch->kid[0].gdv[0] = dv[0];
	branch->kid[0].gdv[1] = sd[0];
	branch->kid[0].gdv[2] = sd[2];
	branch->kid[1].gdv[0] = sd[0];
	branch->kid[1].gdv[1] = dv[1];
	branch->kid[1].gdv[2] = sd[1];
	branch->kid[2].gdv[0] = sd[1];
	branch->kid[2].gdv[1] = dv[2];
	branch->kid[2].gdv[2] = sd[2];
	branch->kid[3].gdv[0] = sd[0];
	branch->kid[3].gdv[1] = sd[1];
	branch->kid[3].gdv[2] = sd[2];
}

/* Recursively subdivide the given node to the specified quadtree depth */
void
branchsample(TRITREE *node, int depth)
{
	int	i;

	if (depth <= 0)
		return;
	if (isleaf(node)) {			/* subdivide leaf node */
		TRITREE	branch, *moved_leaf;
		FVECT	sdir;
		subdivide(&branch, node->gdv);
		decodedir(sdir, node->sd);
		moved_leaf = findleaf(&branch, sdir);
		if (moved_leaf != NULL) {	/* bequeath old sample */
			moved_leaf->sd = node->sd;
			copycolr(moved_leaf->val, node->val);
		}
		for (i = 0; i < 4; i++)		/* compute new samples */
			if (&branch.kid[i] != moved_leaf)
				leafsample(&branch.kid[i]);
		*node = branch;			/* replace leaf with branch */
	}
	for (i = 0; i < 4; i++)			/* subdivide children */
		branchsample(&node->kid[i], depth-1);
}

/* Sample sphere using triangular geodesic mesh */
TRITREE	*
geosample(int nsamps)
{
	int	depth;
	TRITREE	*tree;
	FVECT	trunk[NTRUNKVERT];
	int	i, j;
					/* figure out depth */
	if ((nsamps -= 4) < 0)
		error(USER, "minimum number of samples is 4");
	nsamps = nsamps*3/NTRUNKBR;	/* round up */
	for (depth = 0; nsamps > 1; depth++)
		nsamps >>= 2;
					/* make base tetrahedron */
	tree = (TRITREE *)malloc(sizeof(TRITREE));
	if (tree == NULL) goto memerr;
	trunk[0][0] = trunk[0][1] = 0; trunk[0][2] = 1;
	trunk[1][0] = 0;
	trunk[1][2] = cos(2.*asin(sqrt(2./3.)));
	trunk[1][1] = sqrt(1. - trunk[1][2]*trunk[1][2]);
	spinvector(trunk[2], trunk[1], trunk[0], 2.*PI/3.);
	spinvector(trunk[3], trunk[1], trunk[0], 4.*PI/3.);
	tree->gdv[0] = tree->gdv[1] = tree->gdv[2] = encodedir(trunk[0]);
	tree->kid = (TRITREE *)calloc(NTRUNKBR, sizeof(TRITREE));
	if (tree->kid == NULL) goto memerr;
					/* grow our tree from trunk */
	for (i = 0; i < NTRUNKBR; i++) {
		for (j = 0; j < 3; j++)		/* XXX works for tetra only */
			tree->kid[i].gdv[j] = encodedir(trunk[(i+j)%NTRUNKVERT]);
		leafsample(&tree->kid[i]);
		branchsample(&tree->kid[i], depth);
	}
	return(tree);
memerr:
	error(SYSTEM, "out of memory in geosample()");
	return NULL; /* dummy return */
}

/* Compute leaf exponent histogram */
void
get_ehisto(const TRITREE *node, long exphisto[256])
{
	int	i;

	if (isleaf(node)) {
		++exphisto[node->val[EXP]];
		return;
	}
	for (i = 0; i < 4; i++)
		get_ehisto(&node->kid[i], exphisto);
}

/* Get reasonable source threshold */
double
get_threshold(const TRITREE *tree)
{
	long	samptotal = 0;
	long	exphisto[256];
	int	i;
						/* compute sample histogram */
	memset((void *)exphisto, 0, sizeof(exphisto));
	for (i = 0; i < NTRUNKBR; i++)
		get_ehisto(&tree->kid[i], exphisto);
						/* use 98th percentile */
	for (i = 0; i < 256; i++)
		samptotal += exphisto[i];
	samptotal /= 50; 
	for (i = 256; (--i > 0) & (samptotal > 0); )
		samptotal -= exphisto[i];
	return(ldexp(.75, i-COLXS));
}

/* Find leaf containing the maximum exponent */
TRITREE *
findemax(TRITREE *node, int *expp)
{
	if (!isleaf(node)) {
		TRITREE	*maxleaf;
		TRITREE *rleaf;
		maxleaf = findemax(&node->kid[0], expp);
		rleaf = findemax(&node->kid[1], expp);
		if (rleaf != NULL) maxleaf = rleaf;
		rleaf = findemax(&node->kid[2], expp);
		if (rleaf != NULL) maxleaf = rleaf;
		rleaf = findemax(&node->kid[3], expp);
		if (rleaf != NULL) maxleaf = rleaf;
		return(maxleaf);
	}
	if (node->val[EXP] <= *expp)
		return(NULL);
	*expp = node->val[EXP];
	return(node);
}

/* Compute solid angle of spherical triangle (approx.) */
double
tri_omegav(const int32 vd[3])
{
	FVECT	v[3], e1, e2, vcross;

	decodedir(v[0], vd[0]);
	decodedir(v[1], vd[1]);
	decodedir(v[2], vd[2]);
	VSUB(e1, v[1], v[0]);
	VSUB(e2, v[2], v[1]);
	fcross(vcross, e1, e2);
	return(.5*VLEN(vcross));
}

/* Sum intensity times direction for above-threshold perimiter within radius */
void
vector_sum(FVECT vsum, TRITREE *node,
		FVECT cent, double maxr2, int ethresh)
{
	if (isleaf(node)) {
		double	intens;
		FVECT	sdir;
		if (node->val[EXP] < ethresh)
			return;				/* below threshold */
		if (fdir2diff(node->sd,cent) > maxr2)
			return;				/* too far away */
		intens = colrval(node->val,GRN) * tri_omegav(node->gdv);
		decodedir(sdir, node->sd);
		VSUM(vsum, vsum, sdir, intens);
		return;
	}
	if (dir2diff(node->gdv[0],node->gdv[1]) > maxr2 &&
			fdir2diff(node->gdv[0],cent) < maxr2 &&
			fdir2diff(node->gdv[1],cent) < maxr2 &&
			fdir2diff(node->gdv[2],cent) < maxr2)
		return;					/* containing node */
	vector_sum(vsum, &node->kid[0], cent, maxr2, ethresh);
	vector_sum(vsum, &node->kid[1], cent, maxr2, ethresh);
	vector_sum(vsum, &node->kid[2], cent, maxr2, ethresh);
	vector_sum(vsum, &node->kid[3], cent, maxr2, ethresh);
}

/* Claim source contributions within the given solid angle */
void
claimlight(COLOR intens, TRITREE *node, FVECT cent, double maxr2)
{
	int	remaining;
	int	i;
	if (isleaf(node)) {		/* claim contribution */
		COLOR	contrib;
		if (node->val[EXP] <= 0)
			return;		/* already claimed */
		if (fdir2diff(node->sd,cent) > maxr2)
			return;		/* too far away */
		colr_color(contrib, node->val);
		scalecolor(contrib, tri_omegav(node->gdv));
		addcolor(intens, contrib);
		copycolr(node->val, blkclr);
		return;
	}
	if (dir2diff(node->gdv[0],node->gdv[1]) > maxr2 &&
			fdir2diff(node->gdv[0],cent) < maxr2 &&
			fdir2diff(node->gdv[1],cent) < maxr2 &&
			fdir2diff(node->gdv[2],cent) < maxr2)
		return;			/* previously claimed node */
	remaining = 0;			/* recurse on children */
	for (i = 0; i < 4; i++) {
		claimlight(intens, &node->kid[i], cent, maxr2);
		if (!isleaf(&node->kid[i]) || node->kid[i].val[EXP] != 0)
			++remaining;
	}
	if (remaining)
		return;
					/* consolidate empties */
	free((void *)node->kid); node->kid = NULL;
	copycolr(node->val, blkclr);
	node->sd = node->gdv[0];	/* doesn't really matter */
}

/* Add lost light contribution to the given list */
void
add2lost(LOSTLIGHT **llp, COLOR intens, FVECT cent)
{
	LOSTLIGHT	*newll = (LOSTLIGHT *)malloc(sizeof(LOSTLIGHT));

	if (newll == NULL)
		return;
	copycolor(newll->intens, intens);
	newll->sd = encodedir(cent);
	newll->next = *llp;
	*llp = newll;
}

/* Check lost light list for contributions */
void
getlost(LOSTLIGHT **llp, COLOR intens, FVECT cent, double omega)
{
	const double	maxr2 = omega/PI;
	LOSTLIGHT	lhead, *lastp, *thisp;

	lhead.next = *llp;
	lastp = &lhead;
	while ((thisp = lastp->next) != NULL)
		if (fdir2diff(thisp->sd,cent) <= maxr2) {
			LOSTLIGHT	*mynext = thisp->next;
			addcolor(intens, thisp->intens);
			free((void *)thisp);
			lastp->next = mynext;
		} else
			lastp = thisp;
	*llp = lhead.next;
}

/* Create & print distant sources */
void
mksources(TRITREE *samptree, double thresh, double maxang)
{
	const int	ethresh = (int)(log(thresh)/log(2.) + (COLXS+.5));
	const double	maxomega = 2.*PI*(1. - cos(PI/180./2.*maxang));
	const double	minintens = .05*thresh*maxomega;
	int		nsrcs = 0;
	LOSTLIGHT	*lostlightlist = NULL;
	int		emax;
	TRITREE		*startleaf;
	double		growstep;
	FVECT		curcent;
	double		currad;
	double		curomega;
	COLOR		curintens;
	double		thisthresh;
	int		thisethresh;
	int		i;
	/*
	 * General algorithm:
	 *	1) Start with brightest unclaimed pixel
	 *	2) Grow source toward brightest unclaimed perimeter until:
	 *		a) Source exceeds maximum size, or
	 *		b) Perimeter values all below threshold, or
	 *		c) Source average drops below threshold
	 *	3) Loop until nothing over threshold
	 *
	 * Complexity added to absorb insignificant sources in larger ones.
	 */
	if (thresh <= FTINY)
		return;
	for ( ; ; ) {
		emax = ethresh;		/* find brightest unclaimed */
		startleaf = NULL;
		for (i = 0; i < NTRUNKBR; i++) {
			TRITREE	*bigger = findemax(&samptree->kid[i], &emax);
			if (bigger != NULL)
				startleaf = bigger;
		}
		if (startleaf == NULL)
			break;
					/* claim it */
		decodedir(curcent, startleaf->sd);
		curomega = tri_omegav(startleaf->gdv);
		currad = sqrt(curomega/PI);
		growstep = 3.*currad;
		colr_color(curintens, startleaf->val);
		thisthresh = .15*bright(curintens);
		if (thisthresh < thresh) thisthresh = thresh;
		thisethresh = (int)(log(thisthresh)/log(2.) + (COLXS+.5));
		scalecolor(curintens, curomega);
		copycolr(startleaf->val, blkclr);
		do {			/* grow source */
			FVECT	vsum;
			double	movedist;
			vsum[0] = vsum[1] = vsum[2] = .0;
			for (i = 0; i < NTRUNKBR; i++)
				vector_sum(vsum, &samptree->kid[i],
					curcent, 2.-2.*cos(currad+growstep),
						thisethresh);
			if (normalize(vsum) == .0)
				break;
			movedist = acos(DOT(vsum,curcent));
			if (movedist > growstep) {
				VSUB(vsum, vsum, curcent);
				movedist = growstep/VLEN(vsum); 
				VSUM(curcent, curcent, vsum, movedist);
				normalize(curcent);
			} else
				VCOPY(curcent, vsum);
			currad += growstep;
			curomega = 2.*PI*(1. - cos(currad));
			for (i = 0; i < NTRUNKBR; i++)
				claimlight(curintens, &samptree->kid[i],
						curcent, 2.-2.*cos(currad));
		} while (curomega < maxomega &&
				bright(curintens)/curomega > thisthresh);
		if (bright(curintens) < minintens) {
			add2lost(&lostlightlist, curintens, curcent);
			continue;
		}
					/* absorb lost contributions */
		getlost(&lostlightlist, curintens, curcent, curomega);
		++nsrcs;		/* output source */
		scalecolor(curintens, 1./curomega);
		printf("\nvoid illum IBLout\n");
		printf("0\n0\n3 %f %f %f\n",
				colval(curintens,RED),
				colval(curintens,GRN),
				colval(curintens,BLU));
		printf("\nIBLout source IBLsrc%d\n", nsrcs);
		printf("0\n0\n4 %f %f %f %f\n",
				curcent[0], curcent[1], curcent[2],
				2.*180./PI*currad);
	}
}

int
main(int argc, char *argv[])
{
	long	nsamps = DEF_NSAMPS;
	double	maxang = DEF_MAXANG;
	TRITREE	*samptree;
	double	thresh = 0;
	int	i;
	
	progname = argv[0];
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'd':		/* number of samples */
			if (i >= argc-1) goto userr;
			nsamps = atol(argv[++i]);
			break;
		case 't':		/* manual threshold */
			if (i >= argc-1) goto userr;
			thresh = atof(argv[++i]);
			break;
		case 'a':		/* maximum source angle */
			if (i >= argc-1) goto userr;
			maxang = atof(argv[++i]);
			if (maxang <= FTINY)
				goto userr;
			if (maxang > 180.)
				maxang = 180.;
			break;
		default:
			goto userr;
		}
	if (i < argc-1)
		goto userr;
					/* start our ray calculation */
	directvis = 0;
	ray_init(i == argc-1 ? argv[i] : (char *)NULL);
	VCOPY(scene_cent, thescene.cuorg);
	scene_cent[0] += 0.5*thescene.cusize;
	scene_cent[1] += 0.5*thescene.cusize;
	scene_cent[2] += 0.5*thescene.cusize;
	scene_rad = 0.86603*thescene.cusize;
					/* sample geodesic mesh */
	samptree = geosample(nsamps);
					/* get source threshold */
	if (thresh <= FTINY)
		thresh = get_threshold(samptree);
					/* done with ray samples */
	ray_done(1);
					/* print header */
	printf("# ");
	printargs(argc, argv, stdout);
					/* create & print sources */
	mksources(samptree, thresh, maxang);
					/* all done, no need to clean up */
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-d nsamps][-t thresh][-a maxang] [octree]\n",
			argv[0]);
	exit(1);
}

void
eputs(char  *s)
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	if (s[strlen(s)-1] == '\n') {
		fflush(stderr);
		midline = 0;
	}
}

void
wputs(char *s)
{
	/* no warnings */
}
