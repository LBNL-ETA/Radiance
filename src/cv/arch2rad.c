#ifndef lint
static const char	RCSid[] = "$Id: arch2rad.c,v 2.2 2003/02/22 02:07:23 greg Exp $";
#endif
/*
 * Convert Architrion file to Radiance
 *
 *	Greg Ward
 */

#include <stdio.h>

#include <stdlib.h>

#include "trans.h"

#define DEFMAPFILE	"/usr/local/lib/ray/lib/arch.map"

				/* qualifiers */
#define Q_COL		0
#define Q_FAC		1
#define Q_LAY		2
#define Q_REF		3
#define NQUALS		4

char	*qname[NQUALS] = {
	"Color",
	"Face",
	"Layer",
	"RefId",
};

QLIST	qlist = {NQUALS, qname};
				/* face ids */
#define F_BOT		0
#define F_END		1
#define F_OPP		2
#define F_REF		3
#define F_SIL		4
#define F_TOP		5
#define NFACES		6
ID	faceId[NFACES] = {
	{"bottom"},{"end"},{"opposite"},{"reference"},{"sill"},{"top"}
};
				/* valid qualifier ids */
IDLIST	qual[NQUALS] = {
	{ 0, NULL },
	{ NFACES, faceId },
	{ 0, NULL },
	{ 0, NULL },
};
				/* mapping rules */
RULEHD	*ourmap = NULL;
				/* file header */
struct header {
	char	*filename;
	char	*layer[9];
	char	length_u[16], area_u[16];
	double	length_f, area_f;
	int	nblocks, nopenings;
} fhead;
				/* block or opening geometry */
typedef struct {
	int	x[4], y[4], z[4], h[4];
} PRISM;
				/* macros for x,y,z */
#define p_x(p,i)	((p)->x[(i)&3])
#define p_y(p,i)	((p)->y[(i)&3])
#define p_z(p,i)	((i)&4 ? (p)->h[(i)&3] : (p)->z[i])
				/* opening */
typedef struct {
	PRISM	p;
	int	corner;
	int	depth;
	ID	frame;
} OPNG;
				/* block */
typedef struct {
	short	layer;
	short	color;
	ID	refid;
	PRISM	p;
	int	nopenings;
	OPNG	*opening;
} BLOCK;
				/* point format */
char	ptfmt[] = "\t%12.9g %12.9g %12.9g\n";
				/* antimatter modifier id for openings */
char	openmod[] = "opening";
				/* return flags for checkface() */
#define T1_OK		1
#define T2_OK		2
#define QL_OK		4

char	*progname;		/* argv[0] */


main(argc, argv)		/* translate Architrion file */
int	argc;
char	*argv[];
{
	int	donames = 0;		/* -n flag, produce namelist */
	int	i;
	
	progname = argv[0];
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'n':		/* just produce name list */
			donames++;
			break;
		case 'm':		/* use custom mapfile */
			ourmap = getmapping(argv[++i], &qlist);
			break;
		default:
			goto userr;
		}
	if (i < argc-1)
		goto userr;
	if (i == argc-1)
		if (freopen(argv[i], "r", stdin) == NULL) {
			fprintf(stderr, "%s: cannot open\n", argv[i]);
			exit(1);
		}
	getfhead(stdin);		/* get Architrion header */
	if (donames) {				/* scan for ids */
		arch2names(stdin);
		printf("filename \"%s\"\n", fhead.filename);
		printf("filetype \"Architrion\"\n");
		write_quals(&qlist, qual, stdout);
	} else {				/* translate file */
		if (ourmap == NULL)
			ourmap = getmapping(DEFMAPFILE, &qlist);
		arch2rad(stdin, stdout);
	}
	return(0);
userr:
	fprintf(stderr, "Usage: %s [-n][-m mapfile] [input]\n", argv[0]);
	exit(1);
}


arch2rad(inp, out)	/* translate Architrion file to Radiance */
FILE	*inp, *out;
{
	int	nbs, nos;
	BLOCK	blk;
	
	puthead(out);
	nbs = nos = 0;
	while (getblock(&blk, inp) != EOF) {
		putblock(&blk, out);
		nbs++;
		nos += blk.nopenings;
		doneblock(&blk);
	}
	if (nbs != fhead.nblocks)
		fprintf(stderr,
		"%s: warning -- block count incorrect (%d != %d)\n",
				progname, nbs, fhead.nblocks);
	if (nos != fhead.nopenings)
		fprintf(stderr,
		"%s: warning -- opening count incorrect (%d != %d)\n",
				progname, nos, fhead.nopenings);
}


arch2names(inp)		/* get name list from an Architrion file */
FILE	*inp;
{
	BLOCK	blk;
	
	while (getblock(&blk, inp) != EOF) {
		if (ourmap == NULL || hasmatch(&blk, ourmap))
			add2quals(&blk);
		doneblock(&blk);
	}
}


hasmatch(bp, mp)		/* check for any match in rule list */
BLOCK	*bp;
RULEHD	*mp;
{
	if (mp == NULL)
		return(0);
	if (hasmatch(bp, mp->next))	/* check in order -- more likely */
		return(1);
	return(matchrule(bp, mp));
}


getfhead(fp)			/* get file header */
FILE	*fp;
{
	char	buf[MAXSTR];
	int	i, n;
	register int	c;
					/* get file name */
	if (fgets(buf, MAXSTR, fp) == NULL)
		goto readerr;
	buf[strlen(buf)-1] = '\0';
	fhead.filename = savestr(buf);
					/* get layers */
	fhead.layer[0] = "Worksheet";
	for (i = 1; i <= 8; i++) {
		if (fscanf(fp, "L%*[^0-8]%d", &n) != 1 || n != i)
			goto readerr;
		while ((c = getc(fp)) != EOF && (c == ' ' || c == '\t'))
			;
		if (c == EOF)
			goto readerr;
		ungetc(c, fp);
		if (fgets(buf, MAXSTR, fp) == NULL)
			goto readerr;
		buf[strlen(buf)-1] = '\0';
		if (buf[0])
			fhead.layer[i] = savestr(buf);
		else
			fhead.layer[i] = NULL;
	}
					/* get units */
	if (fgets(buf, MAXSTR, fp) == NULL)
		goto readerr;
	if (sscanf(buf, "0 %*f %s %*f %s\n",
			fhead.length_u, fhead.area_u) != 2) {
		fhead.length_u[0] = '\0';
		fhead.area_u[0] = '\0';
	}
	if (fgets(buf, MAXSTR, fp) == NULL ||
			sscanf(buf, "0 %lf\n", &fhead.length_f) != 1)
		goto readerr;
	if (fgets(buf, MAXSTR, fp) == NULL ||
			sscanf(buf, "0 %lf\n", &fhead.area_f) != 1)
		goto readerr;
					/* get number of blocks and openings */
	if (fgets(buf, MAXSTR, fp) == NULL ||
			sscanf(buf, "0 %d %d\n", &fhead.nblocks,
			&fhead.nopenings) != 2)
		goto readerr;
	return;
readerr:
	fprintf(stderr, "%s: error reading Architrion header\n", progname);
	exit(1);
}


puthead(fp)			/* put out header information */
FILE	*fp;
{
	register int	i;
	
	fprintf(fp, "# File created by: %s\n", progname);
	fprintf(fp, "# Input file: %s\n", fhead.filename);
	fprintf(fp, "# Input units: %s\n", fhead.length_u);
	fprintf(fp, "# Output units: meters\n");
	fprintf(fp, "# Named layers:\n");
	for (i = 0; i < 8; i++)
		if (fhead.layer[i] != NULL)
			fprintf(fp, "#\tLayer No. %d\t%s\n", i, fhead.layer[i]);
}


getblock(bp, fp)			/* get an Architrion block */
register BLOCK	*bp;
FILE	*fp;
{
	char	word[32];
	int	i;
	
	if (fgets(word, sizeof(word), fp) == NULL)
		return(EOF);
	if (strncmp(word, "Block", 5))
		return(EOF);
	if (fscanf(fp, "1 %hd %hd %d", &bp->layer,
			&bp->color, &bp->nopenings) != 3)
		return(EOF);
	if (fgetid(&bp->refid, "\n", fp) == EOF)
		return(EOF);
	if (fscanf(fp, "1 %d %d %d %d %d %d %d %d\n",
			&bp->p.x[0], &bp->p.y[0], &bp->p.z[0], &bp->p.h[0],
			&bp->p.x[1], &bp->p.y[1], &bp->p.z[1], &bp->p.h[1]) != 8)
		return(EOF);
	if (fscanf(fp, "1 %d %d %d %d %d %d %d %d\n",
			&bp->p.x[2], &bp->p.y[2], &bp->p.z[2], &bp->p.h[2],
			&bp->p.x[3], &bp->p.y[3], &bp->p.z[3], &bp->p.h[3]) != 8)
		return(EOF);
	if (bp->nopenings == 0) {
		bp->opening = NULL;
		return(0);
	}
	bp->opening = (OPNG *)malloc(bp->nopenings*sizeof(OPNG));
	if (bp->opening == NULL)
		goto memerr;
	for (i = 0; i < bp->nopenings; i++)
		if (getopening(&bp->opening[i], fp) < 0)
			return(EOF);
	return(0);
memerr:
	fprintf(stderr, "%s: out of memory in getblock\n", progname);
	exit(1);
}


getopening(op, fp)		/* read in opening from fp */
register OPNG	*op;
FILE	*fp;
{
	register int	c;
	char	word[32];

	if (fgets(word, sizeof(word), fp) == NULL)
		return(EOF);
	if (strncmp(word, "Opening", 7))
		return(EOF);
	if (fscanf(fp, "2 %d %d %d %d %d %d %d %d\n",
			&op->p.x[0], &op->p.y[0],
			&op->p.z[0], &op->p.h[0],
			&op->p.x[1], &op->p.y[1],
			&op->p.z[1], &op->p.h[1]) != 8)
		return(EOF);
	if (fscanf(fp, "2 %d %d %d %d %d %d %d %d\n",
			&op->p.x[2], &op->p.y[2],
			&op->p.z[2], &op->p.h[2],
			&op->p.x[3], &op->p.y[3],
			&op->p.z[3], &op->p.h[3]) != 8)
		return(EOF);
	c = getc(fp);
	if (c == '3') {
		if (fscanf(fp, "%d %d", &op->corner, &op->depth) != 2)
			return(EOF);
		if (fgetid(&op->frame, "\n", fp) == EOF)
			return(EOF);
	} else {
		op->corner = -1;
		op->frame.name = NULL;
		if (c != EOF)
			ungetc(c, fp);
	}
	return(0);
}


doneblock(bp)			/* free data associated with bp */
register BLOCK	*bp;
{
	register int	i;

	if (bp->nopenings > 0) {
		for (i = 0; i < bp->nopenings; i++)
			doneid(&bp->opening[i].frame);
		free((void *)bp->opening);
	}
	doneid(&bp->refid);
}


add2quals(bp)			/* add to qualifier lists */
register BLOCK	*bp;
{
	ID	tmpid;

	if (fhead.layer[bp->layer] == NULL) {
		tmpid.name = NULL;
		tmpid.number = bp->layer;
	} else {
		tmpid.name = fhead.layer[bp->layer];
		tmpid.number = 0;
	}
	findid(&qual[Q_LAY], &tmpid, 1);
	tmpid.name = NULL;
	tmpid.number = bp->color;
	findid(&qual[Q_COL], &tmpid, 1);
	findid(&qual[Q_REF], &bp->refid, 1);
}


putblock(bp, fp)			/* put out a block */
BLOCK	*bp;
FILE	*fp;
{
	RULEHD	*rp;
	char	*mod[NFACES], *sillmod;
	int	nmods, i;
	int	donematch, newmatch;
	
	nmods = 0; sillmod = NULL;
	donematch = 0;
	for (rp = ourmap; rp != NULL; rp = rp->next) {
		newmatch = matchrule(bp, rp);		/* test this rule */
		newmatch &= ~donematch;			/* don't repeat */
		if (newmatch == 0)
			continue;
		mod[nmods++] = rp->mnam;
		if (sillmod == NULL && newmatch&F_SIL)
			sillmod = rp->mnam;
		putfaces(rp->mnam, bp, newmatch, fp);	/* put out new faces */
		donematch |= newmatch;
		if (donematch == (1<<NFACES)-1)
			break;				/* done all faces */
	}
						/* put out openings */
	if (donematch && bp->nopenings > 0) {
		putopenmod(sillmod, mod, nmods, fp);
		for (i = 0; i < bp->nopenings; i++)
			putopening(&bp->opening[i], fp);
	}
	if (ferror(fp)) {
		fprintf(stderr, "%s: write error in putblock\n", progname);
		exit(1);
	}
}


putopenmod(sm, ml, nm, fp)		/* put out opening modifier */
char	*sm, *ml[];
int	nm;
FILE	*fp;
{
	int	rept, nrepts;
	register int	i, j;
	
	if (sm == NULL) {		/* if no sill modifier, use list */
		sm = *ml++;
		nm--;
	}
					/* check for repeats */
	rept = 0; nrepts = 0;
	for (i = 0; i < nm; i++) {
		if (ml[i] == sm || !strcmp(ml[i], VOIDID)) {
			rept |= 1<<i;
			nrepts++;
			continue;
		}
		for (j = 0; j < i; j++)
			if (!(rept & 1<<j) && ml[j] == ml[i]) {
				rept |= 1<<j;
				nrepts++;
			}
	}
					/* print antimatter and modlist */
	fprintf(fp, "\n%s antimatter %s\n", VOIDID, openmod);
	fprintf(fp, "%d %s", 1+nm-nrepts, sm);
	for (i = 0; i < nm; i++)
		if (!(rept & 1<<i))
			fprintf(fp, " %s", ml[i]);
	fprintf(fp, "\n0\n0\n");
}


putopening(op, fp)			/* put out an opening */
OPNG	*op;
FILE	*fp;
{
	static int	nopens = 0;
	char	buf[32];
	register PRISM	*p = &op->p;
	PRISM	newp;
	register int	i;
					/* copy original prism */
	for (i = 0; i < 4; i++) {
		newp.x[i] = p->x[i];
		newp.y[i] = p->y[i];
		newp.z[i] = p->z[i];
		newp.h[i] = p->h[i];
	}
					/* spread reference and opposite */
	if (p->x[2] > p->x[0]) {
		newp.x[0] -= 2;
		newp.x[1] -= 2;
		newp.x[2] += 2;
		newp.x[3] += 2;
	} else if (p->x[0] > p->x[2]) {
		newp.x[0] += 2;
		newp.x[1] += 2;
		newp.x[2] -= 2;
		newp.x[3] -= 2;
	}
	if (p->y[2] > p->y[0]) {
		newp.y[0] -= 2;
		newp.y[1] -= 2;
		newp.y[2] += 2;
		newp.y[3] += 2;
	} else if (p->y[0] > p->y[2]) {
		newp.y[0] += 2;
		newp.y[1] += 2;
		newp.y[2] -= 2;
		newp.y[3] -= 2;
	}
						/* put out faces */
	sprintf(buf, "op%d", ++nopens);
	putface(openmod, buf, "ref", &newp, 4, 5, 1, 0, fp);
	putface(openmod, buf, "opp", &newp, 2, 6, 7, 3, fp);
	putface(openmod, buf, "end1", &newp, 5, 6, 2, 1, fp);
	putface(openmod, buf, "end2", &newp, 3, 7, 4, 0, fp);
	putface(openmod, buf, "bot", &newp, 1, 2, 3, 0, fp);
	putface(openmod, buf, "top", &newp, 7, 6, 5, 4, fp);
}


int
matchrule(bp, rp)			/* see if block matches this rule */
register BLOCK	*bp;
register RULEHD	*rp;
{
	register int	i;
	ID	tmpid;
	
	if (rp->qflg & FL(Q_LAY)) {		/* check layer */
		tmpid.name = fhead.layer[bp->layer];
		tmpid.number = bp->layer;
		if (!matchid(&tmpid, &idm(rp)[Q_LAY]))
			return(0);
	}
	if (rp->qflg & FL(Q_COL)) {		/* check color */
		tmpid.name = NULL;
		tmpid.number = bp->color;
		if (!matchid(&tmpid, &idm(rp)[Q_COL]))
			return(0);
	}
	if (rp->qflg & FL(Q_REF)) {		/* check reference id */
		if (!matchid(&bp->refid, &idm(rp)[Q_REF]))
			return(0);
	}
	if (rp->qflg & FL(Q_FAC)) {		/* check faces */
		for (i = 0; i < NFACES; i++)
			if (matchid(&faceId[i], &idm(rp)[Q_FAC]))
				return(1<<i);	/* match single face */
		return(0);
	}
	return((1<<NFACES)-1);			/* matched all faces */
}		


putfaces(m, bp, ff, fp)			/* put out faces */
char	*m;
BLOCK	*bp;
register int	ff;
FILE	*fp;
{
	char	*blkname(), *bn;

	if (!strcmp(m, VOIDID))
		return;
	bn = blkname(bp);
	if (ff & 1<<F_REF)
		putface(m, bn, "ref", &bp->p, 4, 5, 1, 0, fp);
	if (ff & 1<<F_OPP)
		putface(m, bn, "opp", &bp->p, 2, 6, 7, 3, fp);
	if (ff & 1<<F_END) {
		putface(m, bn, "end1", &bp->p, 5, 6, 2, 1, fp);
		putface(m, bn, "end2", &bp->p, 3, 7, 4, 0, fp);
	}
	if (ff & 1<<F_BOT)
		putface(m, bn, "bot", &bp->p, 1, 2, 3, 0, fp);
	if (ff & 1<<F_TOP)
		putface(m, bn, "top", &bp->p, 7, 6, 5, 4, fp);
}


char *
blkname(bp)		/* think up a good name for this block */
register BLOCK	*bp;
{
	static char	nambuf[32];
	static int	blkcnt = 0;
	register char	*nam;
	register int	i, j;

	sprintf(nambuf, "l%d.", bp->layer);
	i = strlen(nambuf);
	nam = bp->refid.name;
	if (nam == NULL) {
		nam = fhead.layer[bp->layer];
		if (nam != NULL)
			i = 0;
	}
	if (nam != NULL) {
		for (j = 0; j < 12 && nam[j]; j++) {
			if (nam[j] == ' ' || nam[j] == '\t')
				nambuf[i++] = '_';
			else
				nambuf[i++] = nam[j];
		}
		nambuf[i++] = '.';
	}
	if (bp->refid.number != 0) {
		sprintf(nambuf+i, "r%d.", bp->refid.number);
		i = strlen(nambuf);
	}
	sprintf(nambuf+i, "b%d", ++blkcnt);
	return(nambuf);
}


putface(m, bn, fn, p, a, b, c, d, fp)	/* put out a face */
char	*m;
char	*bn, *fn;
PRISM	*p;
int	a, b, c, d;
FILE	*fp;
{
	int	cf;
	
	cf = checkface(p, a, b, c, d);
	if (cf & QL_OK) {
		fprintf(fp, "\n%s polygon %s.%s\n", m, bn, fn);
		fprintf(fp, "0\n0\n12\n");
		putpoint(p, a, fp);
		putpoint(p, b, fp);
		putpoint(p, c, fp);
		putpoint(p, d, fp);
		return;
	}
	if (cf & T1_OK) {
		fprintf(fp, "\n%s polygon %s.%sA\n", m, bn, fn);
		fprintf(fp, "0\n0\n9\n");
		putpoint(p, d, fp);
		putpoint(p, a, fp);
		putpoint(p, b, fp);
	}
	if (cf & T2_OK) {
		fprintf(fp, "\n%s polygon %s.%sB\n", m, bn, fn);
		fprintf(fp, "0\n0\n9\n");
		putpoint(p, b, fp);
		putpoint(p, c, fp);
		putpoint(p, d, fp);
	}
}


#define ldot(v1,v2)	(v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])


long
lcross(vr, v1, v2)			/* compute cross product */
register long	vr[3], v1[3], v2[3];
{
	vr[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vr[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vr[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


#define labs(a)		((a)>0 ? (a) : -(a))


checkface(p, a, b, c, d)		/* check a face for validity */
register PRISM	*p;
int	a, b, c, d;
{
	int	rval = 0;
	long	lt;
	long	vc1[3], vc2[3], vt1[3], vt2[3];
	long	vc1l, vc2l, vc1vc2;
					/* check first triangle */
	vt1[0] = p_x(p,b) - p_x(p,a);
	vt1[1] = p_y(p,b) - p_y(p,a);
	vt1[2] = p_z(p,b) - p_z(p,a);
	vt2[0] = p_x(p,d) - p_x(p,a);
	vt2[1] = p_y(p,d) - p_y(p,a);
	vt2[2] = p_z(p,d) - p_z(p,a);
	lcross(vc1, vt1, vt2);
	lt = labs(vc1[0]) + labs(vc1[1]) + labs(vc1[2]);
	if (lt > 1L<<18) lt = 16;
	else if (lt > 1L<<12) lt = 8;
	else lt = 0;
	if (lt) { vc1[0] >>= lt; vc1[1] >>= lt; vc1[2] >>= lt; }
	vc1l = ldot(vc1,vc1);
	if (vc1l > 4)
		rval |= T1_OK;
					/* check second triangle */
	vt1[0] = p_x(p,d) - p_x(p,c);
	vt1[1] = p_y(p,d) - p_y(p,c);
	vt1[2] = p_z(p,d) - p_z(p,c);
	vt2[0] = p_x(p,b) - p_x(p,c);
	vt2[1] = p_y(p,b) - p_y(p,c);
	vt2[2] = p_z(p,b) - p_z(p,c);
	lcross(vc2, vt1, vt2);
	lt = labs(vc2[0]) + labs(vc2[1]) + labs(vc2[2]);
	if (lt > 1L<<18) lt = 16;
	else if (lt > 1L<<12) lt = 8;
	else lt = 0;
	if (lt) { vc2[0] >>= lt; vc2[1] >>= lt; vc2[2] >>= lt; }
	vc2l = ldot(vc2,vc2);
	if (vc2l > 4)
		rval |= T2_OK;
					/* check quadrilateral */
	if (rval == (T1_OK|T2_OK)) {
		vc1vc2 = ldot(vc1,vc2);
		if (vc1vc2*vc1vc2 >= vc1l*vc2l-8)
			rval |= QL_OK;
	}
	return(rval);
}


putpoint(p, n, fp)			/* put out a point */
register PRISM	*p;
int	n;
FILE	*fp;
{
	register int	i = n&3;
	
	fprintf(fp, ptfmt, p->x[i]*fhead.length_f, p->y[i]*fhead.length_f,
			(n&4 ? p->h[i] : p->z[i])*fhead.length_f);
}


void
eputs(s)
char	*s;
{
	fputs(s, stderr);
}


void
quit(code)
int	code;
{
	exit(code);
}
