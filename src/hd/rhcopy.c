#ifndef lint
static const char	RCSid[] = "$Id: rhcopy.c,v 3.38 2023/12/19 20:49:05 greg Exp $";
#endif
/*
 * Copy data into a holodeck file
 */

#include "platform.h"
#include "rterror.h"
#include "holo.h"
#include "view.h"

#ifndef BKBSIZE
#define BKBSIZE		256		/* beam clump size (kilobytes) */
#endif
				/* possible operations */
#define FROM_HOLO	1		/* copy between holodecks */
#define FROM_PICZ	2		/* copy from HDR + depth to holodeck */
#define FROM_STDIN	3		/* copy from stdin to holodeck */
#define TO_STDOUT	4		/* copy rays from holodeck to stdout */

int	operation = 0;		/* what we are doing */
char	*rspec = "";		/* ray details for i/o */
int	checkdepth = 1;		/* check depth (!-d option)? */
int	checkrepeats = 0;	/* check for repeats (-u option)? */
int	nholosects;		/* number of original holodeck sections */
int	iofmt = 'a';		/* input/output format for rays */

				/* holodeck flags */
#define H_BADF		01		/* bad format */
#define H_OBST		02		/* OBSTRUCTIONS= True */
#define H_OBSF		04		/* OBSTRUCTIONS= False */
#define H_VDST		010		/* VDISTANCE= True */
#define H_SWAP		020		/* byte order is different */

char	*progname;		/* global argv[0] */

struct phead {
	VIEW	vw;
	double	expos;
	short	gotview;
	short	badfmt;
	short	altprims;
};

typedef struct {
	FVECT	ro;
	FVECT	rd;
	RREAL	d;
	COLR	cv;
} RAYPAR;

static int openholo(char *fname, int append);
static int addray(RAYPAR *rp);
static int readval(RREAL *v, int n, FILE *fp);
static void readrays(FILE *fp);
static int writeval(RREAL *v, int n, FILE *fp);
static int write_ray(RAYPAR *rp, FILE *fp);
static void writerays(FILE *fp);
static gethfunc holheadline;
static int bpcmp(const void *b1p, const void *b2p);
static int addclump(HOLO *hp, int *bq, int nb);
static void addholo(char *hdf);
static gethfunc picheadline;
static void addpicz(char *pcf, char *zbf);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	i;

	progname = argv[0];
	for (i = 2; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'u':
			checkrepeats = 1;
			break;
		case 'd':
			checkdepth = 0;
			break;
		case 'f':
			iofmt = argv[i][2];
			if (!strchr("afd", iofmt))
				error(USER, "-f? i/o format must be 'a', 'f', or 'd'");
			break;
		case 'h':
			operation = FROM_HOLO;
			break;
		case 'p':
			operation = FROM_PICZ;
			break;
		case 'i':
			operation = FROM_STDIN;
			rspec = argv[i]+2;
			break;
		case 'o':
			operation = TO_STDOUT;
			rspec = argv[i]+2;
			break;
		default:
			goto userr;
		}
	if (!operation | (i > argc-((operation==FROM_HOLO)|(operation==FROM_PICZ))))
		goto userr;
	if (operation == FROM_PICZ && (argc-i)%2)
		goto userr;
	nholosects = openholo(argv[1], (operation != TO_STDOUT));
					/* check requested i/o is compatible */
	if (strchr(rspec, 'l') && !(*(int *)hdlist[0]->priv & H_VDST))
		error(USER, "i/o parameter 'l' incompatible with VDISTANCE=False");
	if (strchr(rspec, 'L') && *(int *)hdlist[0]->priv & H_VDST)
		error(USER, "i/o parameter 'L' incompatible with VDISTANCE=True");

	switch (operation) {		/* perform requested operation */
	case FROM_PICZ:
		for ( ; i < argc; i += 2)
			addpicz(argv[i], argv[i+1]);
		break;
	case FROM_HOLO:
		if (BKBSIZE*1024*1.5 > hdcachesize)
			hdcachesize = BKBSIZE*1024*1.5;
		for ( ; i < argc; i++)
			addholo(argv[i]);
		break;
	case FROM_STDIN:
		readrays(stdin);
		break;
	case TO_STDOUT:
		writerays(stdout);
		break;
	}
	quit(0);
userr:
	fprintf(stderr, "Usage: %s dest.hdk [-u][-d] -h inp1.hdk ..\n",
			progname);
	fprintf(stderr, "   Or: %s dest.hdk [-u][-d] -p inp1.hdr inp1.zbf ..\n",
			progname);
	fprintf(stderr, "   Or: %s dest.hdk [-f{a|f|d}][-u][-d] -i[odplLv]\n",
			progname);
	fprintf(stderr, "   Or: %s src.hdk [-f{a|f|d}] -o[odplLv] ..\n",
			progname);
	exit(1);
}

static int
holheadline(		/* check holodeck header line */
	char	*s,
	void	*vhf
)
{
	int	be;
	char	fmt[MAXFMTLEN];
	int	*hf = (int *)vhf;

	if (formatval(fmt, s)) {
		if (strcmp(fmt, HOLOFMT))
			*hf |= H_BADF;
		else
			*hf &= ~H_BADF;
		return(0);
	}
	if (!strncmp(s, "OBSTRUCTIONS=", 13)) {
		s += 13;
		while (*s == ' ') s++;
		if ((*s == 't') | (*s == 'T'))
			*hf |= H_OBST;
		else if ((*s == 'f') | (*s == 'F'))
			*hf |= H_OBSF;
		else
			error(WARNING, "bad OBSTRUCTIONS value in holodeck");
		return(0);
	}
	if (!strncmp(s, "VDISTANCE=", 10)) {
		s += 10;
		while (*s == ' ') s++;
		if ((*s == 't') | (*s == 'T'))
			*hf |= H_VDST;
		else if ((*s != 'f') % (*s != 'F'))
			error(WARNING, "bad VDISTANCE value in holodeck");
		return(0);
	}
	if ((be = isbigendian(s)) >= 0) {
		if (be != nativebigendian())
			*hf |= H_SWAP;
		return(0);
	}
	return(0);
}

int
openholo(		/* open existing holodeck file for i/o */
	char	*fname,
	int	append
)
{
	FILE	*fp;
	int	fd;
	int	hflags = 0;
	int	*hfstore;
	off_t	nextloc;
	int	n;
					/* open holodeck file */
	if ((fp = fopen(fname, append ? "rb+" : "rb")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\" for %s", fname,
				append ? "appending" : "reading");
		error(SYSTEM, errmsg);
	}
					/* check header and magic number */
	if (getheader(fp, holheadline, &hflags) < 0 ||
			hflags&(H_BADF|H_SWAP) || getw(fp) != HOLOMAGIC) {
		sprintf(errmsg, "holodeck \"%s\" not in expected format", fname);
		error(USER, errmsg);
	}
	fd = dup(fileno(fp));			/* dup file handle */
	nextloc = ftell(fp);			/* get stdio position */
	fclose(fp);				/* done with stdio */
	hfstore = (int *)malloc(sizeof(int));	/* tiny memory leak but who cares? */
	*hfstore = hflags;
	for (n = 0; nextloc > 0L; n++) {	/* initialize each section */
		lseek(fd, nextloc, SEEK_SET);
		read(fd, (char *)&nextloc, sizeof(nextloc));
		hdinit(fd, NULL)->priv = hfstore;
	}
	return(n);
}

int
addray(		/* add a ray to our output holodeck */
	RAYPAR *rp
)
{
	int	sn, bi, n;
	HOLO	*hp;
	GCOORD	gc[2];
	uby8	rr[2][2];
	BEAM	*bp;
	double	d0, d1;
	unsigned	dc;
	RAYVAL	*rv;
	int	nsects = 0;
				/* check each output section */
	for (sn = nholosects; sn--; ) {
		hp = hdlist[sn];
		d0 = hdinter(gc, rr, &d1, hp, rp->ro, rp->rd);
		if (rp->d <= d0 || d1 < -0.001)
			continue;	/* missed section */
		if (checkdepth) {		/* check depth */
			if (*(int *)hp->priv & H_OBST && d0 < -0.001)
				continue;	/* ray starts too late */
			if (*(int *)hp->priv & H_OBSF && rp->d < 0.999*d1)
				continue;	/* ray ends too soon */
		}
		dc = hdcode(hp, rp->d-d0);
		bi = hdbindex(hp, gc);		/* check for duplicates */
		if (checkrepeats && (bp = hdgetbeam(hp, bi)) != NULL) {
			for (n = bp->nrm, rv = hdbray(bp); n--; rv++)
				if ((rv->d == dc || *(int *)hp->priv & (H_OBST|H_OBSF)) &&
						rv->r[0][0] == rr[0][0] &&
						rv->r[0][1] == rr[0][1] &&
						rv->r[1][0] == rr[1][0] &&
						rv->r[1][1] == rr[1][1])
					break;
			if (n >= 0)
				continue;	/* found a matching ray */
		}
		rv = hdnewrays(hp, bi, 1);
		rv->d = dc;
		rv->r[0][0] = rr[0][0]; rv->r[0][1] = rr[0][1];
		rv->r[1][0] = rr[1][0]; rv->r[1][1] = rr[1][1];
		copycolr(rv->v, rp->cv);
		++nsects;
	}
	return nsects;
}

/* Read n-vector from file stream */
static int
readval(RREAL *v, int n, FILE *fp)
{
	int	i;
#ifdef SMLFLT
	double	vd[3];
	switch (iofmt) {
	case 'f':
		return getbinary(v, sizeof(float), n, fp);
	case 'd':
		n = getbinary(vd, sizeof(double), n, fp);
		for (i = n; i-- > 0; ) v[i] = vd[i];
		return n;
	case 'a':
		for (i = 0; i < n; i++)
			if (fscanf(fp, "%f ", &v[i]) != 1)
				break;
		return i;
	}
#else
	float	vf[3];
	switch (iofmt) {
	case 'd':
		return getbinary(v, sizeof(double), n, fp);
	case 'f':
		n = getbinary(vf, sizeof(float), n, fp);
		for (i = n; i-- > 0; ) v[i] = vf[i];
		return n;
	case 'a':
		for (i = 0; i < n; i++)
			if (fscanf(fp, "%lf ", &v[i]) != 1)
				break;
		return i;
	}
#endif
	return -1;
}

#define	GOT_ORG		0x01
#define GOT_DIR		0x02
#define GOT_LEN		0x04
#define GOT_VAL		0x10
#define ALSO_POS	0x20
#define BAD_DIR		0x40
#define BAD_LEN		0x80

/* Read rays from stream and add to holodeck */
static void
readrays(FILE *fp)
{
	unsigned long	nread=0, ngood=0;

	if (iofmt != 'a')
		SET_FILE_BINARY(fp);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	while (!feof(fp)) {		/* read entirety of input */
		RAYPAR	ryp;
		FVECT	pos;
		FVECT	col;
		int	flags = 0;
		int	i;
		for (i = 0; rspec[i]; i++) {
			switch (rspec[i]) {
			case 'o':		/* ray origin */
				if (readval(ryp.ro, 3, fp) < 3)
					break;
				flags |= GOT_ORG;
				continue;
			case 'd':		/* ray direction */
				if (readval(ryp.rd, 3, fp) < 3)
					break;
				if (normalize(ryp.rd) == 0)
					flags |= BAD_DIR;
				else
					flags |= GOT_DIR;
				continue;
			case 'p':		/* ray intersection */
				if (readval(pos, 3, fp) < 3)
					break;
				flags |= ALSO_POS;
				continue;
			case 'L':		/* ray first length */
			case 'l':		/* ray virtual length */
				if (readval(&ryp.d, 1, fp) < 1)
					break;
				if (ryp.d <= FTINY)
					flags |= BAD_LEN;
				else
					flags |= GOT_LEN;
				continue;
			case 'v':		/* ray value */
				if (readval(col, 3, fp) < 3)
					break;
				setcolr(ryp.cv, col[0], col[1], col[2]);
				flags |= GOT_VAL;
				continue;
			default:
				sprintf(errmsg, "unsupported parameter '%c' in -i%s",
						rspec[i], rspec);
				error(USER, errmsg);
			}
			if (!flags)	/* got nothing, so may be normal EOF */
				return;
		}
		++nread;
		if (flags & (BAD_DIR|BAD_LEN))
			continue;	/* just a bad ray is all -- skip */
		if (!(flags & GOT_VAL))
			goto missingData;
		if ((flags & (GOT_ORG|GOT_DIR|GOT_LEN)) != (GOT_ORG|GOT_DIR|GOT_LEN)) {
			if (!(flags & ALSO_POS))
				goto missingData;
			if (flags & GOT_ORG) {
				VSUB(ryp.rd, pos, ryp.ro);
				ryp.d = normalize(ryp.rd);
				if (ryp.d == 0)
					continue;
			} else if ((flags & (GOT_DIR|GOT_LEN)) == (GOT_DIR|GOT_LEN)) {
				VSUM(ryp.ro, pos, ryp.rd, -ryp.d);
			} else
				goto missingData;
		}
		ngood += (addray(&ryp) > 0);	/* add our ray to holodeck */
	}
	return;
missingData:
	sprintf(errmsg, "insufficient data or read error with -i%s after %lu rays read (%lu used)",
			rspec, nread, ngood);
	error(USER, errmsg);
}

/* Write vector value to file stream */
static int
writeval(RREAL *v, int n, FILE *fp)
{
	int	i;

	if (iofmt == 'a') {
		for (i = 0; i < n; i++)
			if (fprintf(fp, "\t%.4e", v[i]) < 0)
				break;
		return i;
	}
#ifdef SMLFLT
	if (iofmt == 'd') {
		double	vd[3];
		for (i = n; i--; ) vd[i] = v[i];
		return putbinary(vd, sizeof(double), n, fp);
	}
#else
	if (iofmt == 'f') {
		float	vf[3];
		for (i = n; i--; ) vf[i] = v[i];
		return putbinary(vf, sizeof(float), n, fp);
	}
#endif
	return putbinary(v, sizeof(*v), n, fp);
}

/* Write out an individual ray as requested */
static int
write_ray(RAYPAR *rp, FILE *fp)
{
	COLOR	cval;
	FVECT	v3;
	char	*typ = rspec;

	for ( ; ; ) {
		switch (*typ++) {
		case 'o':		/* ray origin */
			if (writeval(rp->ro, 3, fp) < 3)
				break;
			continue;
		case 'd':		/* ray direction */
			if (writeval(rp->rd, 3, fp) < 3)
				break;
			continue;
		case 'p':		/* ray intersection */
			VSUM(v3, rp->ro, rp->rd, rp->d);
			if (writeval(v3, 3, fp) < 3)
				break;
			continue;
		case 'L':		/* ray first length */
		case 'l':		/* ray virtual length */
			if (writeval(&rp->d, 1, fp) < 1)
				break;
			continue;
		case 'v':		/* ray value */
			colr_color(cval, rp->cv);
			VCOPY(v3, cval);
			if (writeval(v3, 3, fp) < 3)
				break;
			continue;
		case '\0':		/* end of spec -- success */
			if (iofmt == 'a')
				fputc('\n', fp);
			return(1);
		default:
			sprintf(errmsg, "unsupported parameter '%c' in -o%s", typ[-1], rspec);
		}
		break;			/* land here on error */
	}
	return 0;			/* write error? */
}

static BEAMI	*beamdir;

static int
bpcmp(			/* compare beam positions on disk */
	const void	*b1p,
	const void	*b2p
)
{
	off_t	pdif = beamdir[*(int *)b1p].fo - beamdir[*(int *)b2p].fo;

	if (pdif > 0L) return(1);
	if (pdif < 0L) return(-1);
	return(0);
}

/* Write all rays from holodeck to stream */
static void
writerays(FILE *fp)
{
	int	sn, bi, k;
	GCOORD	gc[2];
	RAYVAL	*rv;
	RAYPAR	ryp;

	if (!*rspec) {
		error(WARNING, "empty -o* output spec, quitting");
		return;
	}
	if (iofmt != 'a')
		SET_FILE_BINARY(fp);
#ifdef getc_unlocked
	flockfile(fp);
#endif
	for (sn = 0; sn < nholosects; sn++) {	/* write each holodeck section */
		HOLO	*hp = hdlist[sn];
		int	nb = nbeams(hp);	/* sort beams by file location */
		int	*bq = (int *)malloc(nb*sizeof(int));
		if (!bq)
			error(SYSTEM, "out of memory in writerays()");
		for (bi = nb; bi--; ) bq[bi] = bi+1;
		beamdir = hp->bi;
		qsort(bq, nb, sizeof(*bq), bpcmp);
		for (bi = 0; bi < nb; bi++) {
			BEAM	*bp = hdgetbeam(hp, bq[bi]);
			if (!bp)		/* empty beam? */
				continue;
			hdbcoord(gc, hp, bq[bi]);
			rv = hdbray(bp);
			for (k = bp->nrm; k--; rv++) {
				ryp.d = hdray(ryp.ro, ryp.rd, hp, gc, rv->r);
				if (*(int *)hp->priv & H_OBSF)
					VSUM(ryp.ro, ryp.ro, ryp.rd, ryp.d);
				else
					ryp.d = 0.;
				ryp.d = hddepth(hp, rv->d) - ryp.d;
				copycolr(ryp.cv, rv->v);
				if (!write_ray(&ryp, fp)) {
					free(bq);
					goto writError;
				}
			}
			hdfreebeam(hp, bq[bi]);
		}
		free(bq);
	}
	if (fflush(fp) != EOF)
		return;
writError:
	error(SYSTEM, "error writing holodeck rays");
}

static int
addclump(		/* transfer the given clump and free */
	HOLO	*hp,
	int	*bq,
	int	nb
)
{
	GCOORD	gc[2];
	RAYPAR	ryp;
	RAYVAL	*rv;
	int	i;
	int	k;
	BEAM	*bp;
					/* sort based on file position */
	beamdir = hp->bi;
	qsort(bq, nb, sizeof(*bq), bpcmp);
					/* transfer each beam */
	for (i = 0; i < nb; i++) {
		bp = hdgetbeam(hp, bq[i]);
		hdbcoord(gc, hp, bq[i]);
		rv = hdbray(bp);			/* add each ray to output */
		for (k = bp->nrm; k--; rv++) {
			ryp.d = hdray(ryp.ro, ryp.rd, hp, gc, rv->r);
			if (*(int *)hp->priv & H_OBSF)
				VSUM(ryp.ro, ryp.ro, ryp.rd, ryp.d);
			else
				ryp.d = 0.;
			ryp.d = hddepth(hp, rv->d) - ryp.d;
			copycolr(ryp.cv, rv->v);
			addray(&ryp);
		}
		hdfreebeam(hp, bq[i]);		/* free the beam */
	}
	return(0);
}


void
addholo(			/* add a holodeck file */
	char	*hdf
)
{
	int	fd;
					/* open the holodeck for reading */
	openholo(hdf, 0);
	fd = hdlist[nholosects]->fd;	/* remember the file handle */
	while (hdlist[nholosects] != NULL) {	/* load each section */
							/* clump the beams */
		clumpbeams(hdlist[nholosects], 0, BKBSIZE*1024, addclump);
		hddone(hdlist[nholosects]);		/* free the section */
	}
	close(fd);			/* close input file */
	hdflush(NULL);			/* flush output */
}



static int
picheadline(		/* process picture header line */
	char	*s,
	void	*vph
)
{
	char	fmt[32];
	struct phead *ph = (struct phead *)vph;

	if (formatval(fmt, s)) {
		ph->badfmt = strcmp(fmt, COLRFMT);
		return(0);
	}
	if (isprims(s)) {
		ph->altprims++;		/* don't want to deal with this */
		return(0);
	}
	if (isexpos(s)) {
		ph->expos *= exposval(s);
		return(0);
	}
	if (isview(s)) {
		ph->gotview += sscanview(&ph->vw, s);
		return(0);
	}
	return(0);
}


void
addpicz(		/* add a picture + depth-buffer */
	char	*pcf,
	char	*zbf
)
{
	FILE	*pfp;
	int	zfd;
	COLR	*cscn;
	float	*zscn;
	struct phead	phd;
	int	eshft;
	double	emult;
	RESOLU	prs;
	RREAL	vl[2];
	RAYPAR	ryp;
	double	aftd;
	int	j, i;
				/* open picture & get header */
	if ((pfp = fopen(pcf, "rb")) == NULL) {
		sprintf(errmsg, "cannot open picture file \"%s\"", pcf);
		error(SYSTEM, pcf);
	}
	phd.vw = stdview;
	phd.expos = 1.0;
	phd.badfmt = phd.gotview = phd.altprims = 0;
	if (getheader(pfp, picheadline, &phd) < 0 ||
			phd.badfmt || !fgetsresolu(&prs, pfp)) {
		sprintf(errmsg, "bad format for picture file \"%s\"", pcf);
		error(USER, errmsg);
	}
	if (!phd.gotview || setview(&phd.vw) != NULL) {
		sprintf(errmsg, "missing/illegal view in picture \"%s\"",
				pcf);
		error(USER, errmsg);
	}
	if (phd.altprims) {
		sprintf(errmsg, "ignoring color primaries in picture \"%s\"",
				pcf);
		error(WARNING, errmsg);
	}
				/* open depth buffer */
	if ((zfd = open_float_depth(zbf, prs.xr*prs.yr)) < 0)
		quit(1);
				/* figure out what to do about exposure */
	if ((phd.expos < 0.99) | (phd.expos > 1.01)) {
		emult = -log(phd.expos)/log(2.);
		eshft = emult >= 0. ? emult+.5 : emult-.5;
		emult -= (double)eshft;
		if ((emult <= 0.01) & (emult >= -0.01))
			emult = -1.;
		else {
			emult = 1./phd.expos;
			eshft = 0;
		}
	} else {
		emult = -1.;
		eshft = 0;
	}
				/* allocate buffers */
	cscn = (COLR *)malloc(scanlen(&prs)*sizeof(COLR));
	zscn = (float *)malloc(scanlen(&prs)*sizeof(float));
	if ((cscn == NULL) | (zscn == NULL))
		error(SYSTEM, "out of memory in addpicz");
				/* read and process each scanline */
	for (j = 0; j < numscans(&prs); j++) {
		i = scanlen(&prs);			/* read colrs */
		if (freadcolrs(cscn, i, pfp) < 0) {
			sprintf(errmsg, "error reading picture \"%s\"", pcf);
			error(USER, errmsg);
		}
		if (eshft)				/* shift exposure */
			shiftcolrs(cscn, i, eshft);
							/* read depth */
		if (read(zfd, zscn, i*sizeof(float)) != i*sizeof(float)) {
			sprintf(errmsg, "error reading depth file \"%s\"", zbf);
			error(USER, errmsg);
		}
		while (i--) {				/* process each pixel */
			if (zscn[i] <= 0.0)
				continue;		/* illegal depth */
			pix2loc(vl, &prs, i, j);
			aftd = viewray(ryp.ro, ryp.rd, &phd.vw, vl[0], vl[1]);
			if (aftd < -FTINY)
				continue;		/* off view */
			if (aftd > FTINY && zscn[i] > aftd)
				continue;		/* aft clipped */
			ryp.d = (RREAL)zscn[i];
			copycolr(ryp.cv, cscn[i]);
			if (emult > 0.) {		/* whatta pain */
				COLOR	ctmp;
				colr_color(ctmp, ryp.cv);
				scalecolor(ctmp, emult);
				setcolr(ryp.cv, colval(ctmp,RED),
					colval(ctmp,GRN), colval(ctmp,BLU));
			}
			addray(&ryp);
		}
	}
				/* write output and free beams */
	hdflush(NULL);
				/* clean up */
	free((void *)cscn);
	free((void *)zscn);
	fclose(pfp);
	close(zfd);
}


void
eputs(			/* put error message to stderr */
	const char  *s
)
{
	static int  midline = 0;

	if (!*s)
		return;
	if (!midline++) {	/* prepend line with program name */
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
quit(			/* exit the program gracefully */
	int	code
)
{
	hdsync(NULL, 1);	/* write out any buffered data */
	exit(code);
}
