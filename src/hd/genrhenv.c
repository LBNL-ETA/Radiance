#ifndef lint
static const char	RCSid[] = "$Id: genrhenv.c,v 3.6 2003/10/20 16:01:55 greg Exp $";
#endif
/*
 * Create a closed environment from a holodeck section
 */

#include "holo.h"
#include "platform.h"

#define	ourhp		(hdlist[0])

char	*progname;		/* global argv[0] */
char	*prefix;		/* file name prefix */

HDBEAMI	*gabmi;			/* beam index array */


main(argc, argv)
int	argc;
char	*argv[];
{
	int	sect;
	register int	n;

	progname = argv[0];
	if (argc < 3 | argc > 4)
		goto userr;
	prefix = argv[1];
	sect = argc==3 ? 0 : atoi(argv[3]);
	openholo(argv[2], sect);
	mkcalfile(argv[2], sect);
	n = ourhp->grid[0] * ourhp->grid[1];
	if (ourhp->grid[0] * ourhp->grid[2] > n)
		n = ourhp->grid[0] * ourhp->grid[2];
	if (ourhp->grid[1] * ourhp->grid[2] > n)
		n = ourhp->grid[1] * ourhp->grid[2];
	gabmi = (HDBEAMI *)malloc(n*sizeof(HDBEAMI));
	if (gabmi == NULL)
		error(SYSTEM, "out of memory in main");
	while (n--)
		gabmi[n].h = ourhp;
	fputs("# ", stdout);
	printargs(argc, argv, stdout);
	hdcachesize = 0;
	for (n = 0; n < 6; n++)
		mkwall(argv[2], sect, n);
	quit(0);
userr:
	fprintf(stderr, "Usage: %s prefix input.hdk [section]\n", progname);
	exit(1);
}


openholo(fname, sect)		/* open holodeck section for input */
char	*fname;
int	sect;
{
	FILE	*fp;
	int	fd;
	int32	nextloc;
	int	n;
					/* open holodeck file */
	if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
					/* check header and magic number */
	if (checkheader(fp, HOLOFMT, NULL) < 0 || getw(fp) != HOLOMAGIC) {
		sprintf(errmsg, "file \"%s\" not in holodeck format", fname);
		error(USER, errmsg);
	}
	fd = dup(fileno(fp));			/* dup file handle */
	nextloc = ftell(fp);			/* get stdio position */
	fclose(fp);				/* done with stdio */
	for (n = 0; nextloc > 0L; n++) {	/* get the indicated section */
		lseek(fd, (off_t)nextloc, SEEK_SET);
		read(fd, (char *)&nextloc, sizeof(nextloc));
		if (n == sect) {
			hdinit(fd, NULL);
			return;
		}
	}
	error(USER, "holodeck section not found");
}


mkwall(hdk, sect, wn)		/* make map for indicated wall */
char	*hdk;
int	sect;
int	wn;
{
	static char	cname[3][6] = {"Red", "Green", "Blue"};
	char	fnbuf[128];
	FILE	*cdat[3];
	int	st = 0;
	int	wo, i, j;
					/* loop through each opposing wall */
	for (wo = 0; wo < 6; wo++) {
		if (wo == wn) continue;
						/* create data files */
		sprintf(fnbuf, "%s_%dr%d.dat", prefix, wn, wo);
		if ((cdat[0] = fopen(fnbuf, "w")) == NULL)
			goto openerr;
		sprintf(fnbuf, "%s_%dg%d.dat", prefix, wn, wo);
		if ((cdat[1] = fopen(fnbuf, "w")) == NULL)
			goto openerr;
		sprintf(fnbuf, "%s_%db%d.dat", prefix, wn, wo);
		if ((cdat[2] = fopen(fnbuf, "w")) == NULL)
			goto openerr;
						/* write dimensions */
		for (i = 0; i < 3; i++) {
			fprintf(cdat[i],
"# %s channel for wall %d -> wall %d of section %d in holodeck \"%s\"\n",
					cname[i], wo, wn, sect, hdk);
			fprintf(cdat[i],
		"4\n0.5 %.1f %d\n0.5 %.1f %d\n0.5 %.1f %d\n0.5 %.1f %d\n",
			ourhp->grid[hdwg0[wn]]-.5, ourhp->grid[hdwg0[wn]],
			ourhp->grid[hdwg1[wn]]-.5, ourhp->grid[hdwg1[wn]],
			ourhp->grid[hdwg0[wo]]-.5, ourhp->grid[hdwg0[wo]],
			ourhp->grid[hdwg1[wo]]-.5, ourhp->grid[hdwg1[wo]]);
		}
						/* run each grid cell */
		for (i = 0; i < ourhp->grid[hdwg0[wn]]; i++)
			for (j = 0; j < ourhp->grid[hdwg1[wn]]; j++)
				dogrid(wn, i, j, wo, cdat);
						/* close files */
		if (fclose(cdat[0]) == EOF | fclose(cdat[1]) == EOF |
				fclose(cdat[2]) == EOF)
			error(SYSTEM, "write error in mkwall");
						/* put out pattern */
		printf("\n%s colordata wallmap\n", st++ ? "wallmap" : "void");
		printf("11 f f f %s_%dr%d.dat %s_%dg%d.dat %s_%db%d.dat\n",
				prefix,wn,wo, prefix,wn,wo, prefix,wn,wo);
		printf("\t%s.cal nu nv ou ov\n0\n2 %d %d\n", prefix, wn, wo);
	}
						/* put out polygon */
	printf("\nwallmap glow wallglow\n0\n0\n4 1 1 1 0\n");
	printf("\nwallglow polygon wall%d\n0\n0\n12\n", wn);
	wallverts(wn);
	return;
openerr:
	sprintf(errmsg, "cannot open \"%d\" for writing", fnbuf);
	error(SYSTEM, errmsg);
}


wallverts(wn)			/* print out vertices for wall wn */
int	wn;
{
	int	rev, i0, i1;
	FVECT	corn, vt;

	if ((rev = wn & 1))
		VSUM(corn, ourhp->orig, ourhp->xv[wn>>1], 1.);
	else
		VCOPY(corn, ourhp->orig);
	fcross(vt, ourhp->xv[i0=hdwg0[wn]], ourhp->xv[i1=hdwg1[wn]]);
	if (DOT(vt, ourhp->wg[wn>>1]) < 0.)
		rev ^= 1;
	if (rev) {
		i0 = hdwg1[wn]; i1 = hdwg0[wn];
	}
	printf("\t%.10e %.10e %.10e\n", corn[0], corn[1], corn[2]);
	printf("\t%.10e %.10e %.10e\n", corn[0] + ourhp->xv[i0][0],
			corn[1] + ourhp->xv[i0][1], corn[2] + ourhp->xv[i0][2]);
	printf("\t%.10e %.10e %.10e\n",
			corn[0] + ourhp->xv[i0][0] + ourhp->xv[i1][0],
			corn[1] + ourhp->xv[i0][1] + ourhp->xv[i1][1],
			corn[2] + ourhp->xv[i0][2] + ourhp->xv[i1][2]);
	printf("\t%.10e %.10e %.10e\n", corn[0] + ourhp->xv[i1][0],
			corn[1] + ourhp->xv[i1][1], corn[2] + ourhp->xv[i1][2]);
}


dogrid(wn, ci, cj, wo, dfp)	/* compute and write grid cell data */
int	wn, ci, cj, wo;
FILE	*dfp[3];
{
	GCOORD	gc[2];
	COLOR	cavg;
	register int	n;

	gc[1].w = wn; gc[1].i[0] = ci; gc[1].i[1] = cj;
	gc[0].w = wo;
	n = 0;				/* load beams in optimal order */
	for (gc[0].i[0] = ourhp->grid[hdwg0[gc[0].w]]; gc[0].i[0]--; )
		for (gc[0].i[1] = ourhp->grid[hdwg1[gc[0].w]]; gc[0].i[1]--; )
			gabmi[n++].b = hdbindex(ourhp, gc);
	hdloadbeams(gabmi, n, NULL);
					/* run beams in regular order */
	fprintf(dfp[0], "\n# Begin grid cell (%d,%d)\n", ci, cj);
	fprintf(dfp[1], "\n# Begin grid cell (%d,%d)\n", ci, cj);
	fprintf(dfp[2], "\n# Begin grid cell (%d,%d)\n", ci, cj);
	for (gc[0].i[0] = 0; gc[0].i[0] < ourhp->grid[hdwg0[gc[0].w]];
			gc[0].i[0]++) {
		for (gc[0].i[1] = 0; gc[0].i[1] < ourhp->grid[hdwg1[gc[0].w]];
				gc[0].i[1]++) {
			avgbeam(cavg, gc);
			fprintf(dfp[0], " %.2e", colval(cavg,RED));
			fprintf(dfp[1], " %.2e", colval(cavg,GRN));
			fprintf(dfp[2], " %.2e", colval(cavg,BLU));
		}
		fputc('\n', dfp[0]);
		fputc('\n', dfp[1]);
		fputc('\n', dfp[2]);
	}
	hdflush(NULL);			/* flush cache */
}


avgbeam(cavg, gc)		/* compute average beam color */
COLOR	cavg;
GCOORD	gc[2];
{
	static COLOR	crun = BLKCOLOR;
	static int	ncrun = 0;
	register BEAM	*bp;
	COLOR	cray;
	unsigned	beamlen;
	int	ls[2][3];
	FVECT	gp[2], wp[2];
	double	d;
	int	n;
	register int	i;

	bp = hdgetbeam(ourhp, hdbindex(ourhp, gc));
	if (bp == NULL || !bp->nrm) {		/* use running average */
		copycolor(cavg, crun);
		if (ncrun > 1) {
			d = 1./ncrun;
			scalecolor(cavg, d);
		}
	} else {				/* use beam average */
		hdlseg(ls, ourhp, gc);			/* get beam length */
		for (i = 3; i--; ) {
			gp[0][i] = ls[0][i] + .5;
			gp[1][i] = ls[1][i] + .5;
		}
		hdworld(wp[0], ourhp, gp[0]);
		hdworld(wp[1], ourhp, gp[1]);
		d = sqrt(dist2(wp[0], wp[1]));
		beamlen = hdcode(ourhp, d);
		setcolor(cavg, 0., 0., 0.);
		n = 0;
		for (i = bp->nrm; i--; ) {
			if (hdbray(bp)[i].d < beamlen)
				continue;		/* inside section */
			colr_color(cray, hdbray(bp)[i].v);
			addcolor(cavg, cray);
			n++;
		}
		if (n > 1) {
			d = 1./n;
			scalecolor(cavg, d);
		}
		if (n) {
			addcolor(crun, cavg);
			ncrun++;
		}
	}
}


mkcalfile(hdk, sect)		/* create .cal file */
char	*hdk;
int	sect;
{
	char	fname[128];
	FILE	*fp;
	register int	i;

	sprintf(fname, "%s.cal", prefix);
	if ((fp = fopen(fname, "w")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\" for writing", fname);
		error(SYSTEM, errmsg);
	}
	fprintf(fp,
	"{\n\tCoordinate calculation for section %d of holodeck \"%s\"\n}\n",
			sect, hdk);

	fprintf(fp, "\n{ Constants }\n");
	fprintf(fp, "hox : %.7e; hoy : %.7e; hoz : %.7e;\n",
			ourhp->orig[0], ourhp->orig[1], ourhp->orig[2]);
	fprintf(fp, "grid(i) : select(i+1, %d, %d, %d);\n",
			ourhp->grid[0], ourhp->grid[1], ourhp->grid[2]);
	fprintf(fp, "hwgx(i) : select(i+1, %.7e, %.7e, %.7e);\n",
			ourhp->wg[0][0], ourhp->wg[1][0], ourhp->wg[2][0]);
	fprintf(fp, "hwgy(i) : select(i+1, %.7e, %.7e, %.7e);\n",
			ourhp->wg[0][1], ourhp->wg[1][1], ourhp->wg[2][1]);
	fprintf(fp, "hwgz(i) : select(i+1, %.7e, %.7e, %.7e);\n",
			ourhp->wg[0][2], ourhp->wg[1][2], ourhp->wg[2][2]);
	fprintf(fp, "hwo(i) : select(i+1, %.7e, %.7e, %.7e,\n",
			ourhp->wo[0], ourhp->wo[1], ourhp->wo[2]);
	fprintf(fp, "\t\t%.7e, %.7e, %.7e);\n",
			ourhp->wo[3], ourhp->wo[4], ourhp->wo[5]);
	fprintf(fp,
"wgp(i,x,y,z) : (x-hox)*hwgx(i) + (y-hoy)*hwgy(i) + (z-hoz)*hwgz(i);\n");
	fprintf(fp, "wg0(w) : select(.5*w+.75, 1, 2, 0);\n");
	fprintf(fp, "wg1(w) : select(.5*w+.75, 2, 0, 1);\n");

	fprintf(fp, "\n{ Variables }\n");
	fprintf(fp, "nw = arg(1); ow = arg(2);\n");
	fprintf(fp, "oax = .5*ow-.25;\n");
	fprintf(fp, "ng0 = wg0(nw); ng1 = wg1(nw);\n");
	fprintf(fp, "og0 = wg0(ow); og1 = wg1(ow);\n");
	fprintf(fp,
	"odenom = Dx*hwgx(oax) + Dy*hwgy(oax) + Dz*hwgz(oax);\n");
	fprintf(fp, "odist = if(and(FTINY-odenom, FTINY+odenom), -1,\n");
	fprintf(fp,
"\t\t(Px*hwgx(oax) + Py*hwgy(oax) + Pz*hwgz(oax) - hwo(ow))/odenom);\n");
	fprintf(fp,
	"opx = Px - odist*Dx; opy = Py - odist*Dy; opz = Pz - odist*Dz;\n");
	fprintf(fp, "nu = wgp(ng0,Px,Py,Pz);\n");
	fprintf(fp, "nv = wgp(ng1,Px,Py,Pz);\n");
	fprintf(fp, "ou = wgp(og0,opx,opy,opz);\n");
	fprintf(fp, "ov = wgp(og1,opx,opy,opz);\n");
	fprintf(fp, "f(v) = if(and(and(ou-FTINY,grid(og0)-ou-FTINY),\n");
	fprintf(fp, "\t\tand(ov-FTINY,grid(og1)-ov-FTINY)), if(v,v,0), 1);\n");
	fclose(fp);
}


void
eputs(s)			/* put error message to stderr */
register char  *s;
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
