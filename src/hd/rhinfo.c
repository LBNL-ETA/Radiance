#ifndef lint
static const char	RCSid[] = "$Id: rhinfo.c,v 3.6 2003/05/29 16:26:22 greg Exp $";
#endif
/*
 * Get general information on holodeck file
 */

#include "holo.h"

#ifndef NHBINS
#define NHBINS		12	/* number of histogram bins to use for stats */
#endif

char	*progname;		/* global argv[0] */

long	beamtot, samptot;	/* total beams and samples */


main(argc, argv)
int	argc;
char	*argv[];
{
	int	sect;

	progname = argv[0];
	if (argc != 2)
		goto userr;
	gethdinfo(argv[1], stdout);
	quit(0);
userr:
	fprintf(stderr, "Usage: %s input.hdk\n", progname);
	exit(1);
}


gethdinfo(fname, fout)		/* get information on holodeck */
char	*fname;
FILE	*fout;
{
	FILE	*fp;
	HOLO	*hdsect;
	int	fd;
	int4	nextloc;
	int	n;
					/* open holodeck file */
	if ((fp = fopen(fname, "r")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\"", fname);
		error(SYSTEM, errmsg);
	}
					/* check header and magic number */
	if (checkheader(fp, HOLOFMT, fout) < 0 || getw(fp) != HOLOMAGIC) {
		sprintf(errmsg, "file \"%s\" not in holodeck format", fname);
		error(USER, errmsg);
	}
	fd = dup(fileno(fp));			/* dup file handle */
	nextloc = ftell(fp);			/* get stdio position */
	fclose(fp);				/* done with stdio */
	for (n = 0; nextloc > 0L; n++) {	/* get the section(s) */
		lseek(fd, (off_t)nextloc, 0);
		read(fd, (char *)&nextloc, sizeof(nextloc));
		fprintf(fout, "Section %d:\n", n);
		hdsect = hdinit(fd, NULL);	/* load section directory */
		psectstats(hdsect, fout);	/* print section statistics */
	}
	nextloc = hdfilen(fd);			/* print global statistics */
	fputs("=====================================================\n", fout);
	fprintf(fout, "Total samples/beams: %ld/%ld (%.2f samples/beam)\n",
			samptot, beamtot, (double)samptot/beamtot);
	fprintf(fout, "%.1f Mbyte file, %.1f%% fragmentation\n",
			nextloc/(1024.*1024.),
			100.*(nextloc-hdfiluse(fd,1))/nextloc);
						/* don't bother with cleanup */
#if 0
	hddone(NULL);				/* free sections */
	close(fd);				/* done with the holodeck */
#endif
}


psectstats(hp, fp)		/* print statistical information for section */
register HOLO	*hp;
FILE	*fp;
{
	int	scount[NHBINS];
	int	minsamp = 10000, maxsamp = 0;
	FVECT	vt;
	double	sqrtmaxp;
	int	bmin, bmax, cnt;
	register int	i;

	fcross(vt, hp->xv[0], hp->xv[1]);
	fprintf(fp, "\tWorld volume:            %g\n", fabs(DOT(vt,hp->xv[2])));
	fprintf(fp, "\tGrid resolution:         %d x %d x %d\n",
			hp->grid[0], hp->grid[1], hp->grid[2]);
	fprintf(fp, "\tNumber of beams:         %ld\n", (long)nbeams(hp));
	beamtot += nbeams(hp);
	fprintf(fp, "\tNumber of ray samples:   %ld\n", (long)biglob(hp)->nrd);
	samptot += biglob(hp)->nrd;
	if (biglob(hp)->nrd <= 0)
		return;				/* no samples to stat! */
	for (i = nbeams(hp); i > 0; i--) {
		if (hp->bi[i].nrd < minsamp)
			minsamp = hp->bi[i].nrd;
		if (hp->bi[i].nrd > maxsamp)
			maxsamp = hp->bi[i].nrd;
	}
	sqrtmaxp = sqrt(maxsamp+1.0);
	for (i = NHBINS; i--; )
		scount[i] = 0;
	for (i = nbeams(hp); i > 0; i--)
		scount[(int)(NHBINS*sqrt((double)hp->bi[i].nrd)/sqrtmaxp)]++;
	for (cnt = 0, i = 0; i < NHBINS && cnt<<1 < nbeams(hp); i++)
		cnt += scount[i];
	fprintf(fp, "\tSamples per beam:        [min,med,max]= [%d,%.0f,%d]\n",
			minsamp,
			(i-.5)*(i-.5)*(maxsamp+1)/(NHBINS*NHBINS),
			maxsamp);
	fprintf(fp, "\tHistogram: [minsamp,maxsamp)= #beams\n");
	bmax = 0;
	for (i = 0; i < NHBINS; i++) {
		bmin = bmax;
		bmax = (i+1)*(i+1)*(maxsamp+1)/(NHBINS*NHBINS);
		fprintf(fp, "\t\t[%d,%d)= %d\n", bmin, bmax, scount[i]);
	}
}
