#ifndef lint
static const char	RCSid[] = "$Id: rhoptimize.c,v 3.13 2003/07/14 20:02:29 schorsch Exp $";
#endif
/*
 * Optimize holodeck for quick access.
 *
 *	11/4/98		Greg Ward Larson
 */

#include <signal.h>
#include <string.h>

#include "rtprocess.h" /* getpid() */
#include "holo.h"

#ifndef BKBSIZE
#define BKBSIZE		256		/* beam clump size (kilobytes) */
#endif

char	*progname;
char	tempfile[128];
int	dupchecking = 0;

extern long	rhinitcopy();


main(argc, argv)
int	argc;
char	*argv[];
{
	char	*inpname, *outname;
	int	hdfd[2];
	long	nextipos, lastopos, thisopos;

	progname = argv[0];
	argv++; argc--;			/* duplicate checking flag? */
	if (argc > 1 && !strcmp(argv[0], "-u")) {
		dupchecking++;
		argv++; argc--;
	}
	if (argc < 1 | argc > 2) {
		fprintf(stderr, "Usage: %s [-u] input.hdk [output.hdk]\n",
				progname);
		exit(1);
	}
	inpname = argv[0];		/* get input file */
	argv++; argc--;
	if (argc == 1)			/* use given output file */
		outname = argv[0];
	else {				/* else use temporary file */
		if (access(inpname, R_OK|W_OK) < 0) {	/* check permissions */
			sprintf(errmsg, "cannot access \"%s\"", inpname);
			error(SYSTEM, errmsg);
		}
		strcpy(tempfile, inpname);
		if ((outname = strrchr(tempfile, '/')) != NULL)
			outname++;
		else
			outname = tempfile;
		sprintf(outname, "rho%d.hdk", getpid());
		outname = tempfile;
	}
					/* copy holodeck file header */
	nextipos = rhinitcopy(hdfd, inpname, outname);
	lastopos = 0L;			/* copy sections one by one */
	while (nextipos != 0L) {
					/* set input position; get next */
		lseek(hdfd[0], (off_t)nextipos, 0);
		read(hdfd[0], (char *)&nextipos, sizeof(nextipos));
					/* get output position; set last */
		thisopos = lseek(hdfd[1], (off_t)0, 2);
		if (lastopos > 0L) {
			lseek(hdfd[1], (off_t)lastopos, 0);
			write(hdfd[1], (char *)&thisopos, sizeof(thisopos));
			lseek(hdfd[1], (off_t)0, 2);
		}
		lastopos = thisopos;
		thisopos = 0L;		/* write place holder */
		write(hdfd[1], (char *)&thisopos, sizeof(thisopos));
					/* copy holodeck section */
		copysect(hdfd[0], hdfd[1]);
	}
					/* clean up */
	close(hdfd[0]);
	close(hdfd[1]);
	if (outname == tempfile && rename(outname, inpname) < 0) {
		sprintf(errmsg, "cannot rename \"%s\" to \"%s\"",
				outname, inpname);
		error(SYSTEM, errmsg);
	}
	exit(0);
}


long
rhinitcopy(hfd, infn, outfn)	/* open files and copy header */
int	hfd[2];			/* returned file descriptors */
char	*infn, *outfn;
{
	FILE	*infp, *outfp;
	long	ifpos;
					/* open files for i/o */
	if ((infp = fopen(infn, "r")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\" for reading", infn);
		error(SYSTEM, errmsg);
	}
	if (access(outfn, F_OK) == 0) {
		sprintf(errmsg, "output file \"%s\" already exists!", outfn);
		error(USER, errmsg);
	}
	if ((outfp = fopen(outfn, "w+")) == NULL) {
		sprintf(errmsg, "cannot open \"%s\" for writing", outfn);
		error(SYSTEM, errmsg);
	}
					/* set up signal handling */
	if (signal(SIGINT, quit) == SIG_IGN) signal(SIGINT, SIG_IGN);
	if (signal(SIGHUP, quit) == SIG_IGN) signal(SIGHUP, SIG_IGN);
	if (signal(SIGTERM, quit) == SIG_IGN) signal(SIGTERM, SIG_IGN);
#ifdef SIGXCPU
	if (signal(SIGXCPU, quit) == SIG_IGN) signal(SIGXCPU, SIG_IGN);
	if (signal(SIGXFSZ, quit) == SIG_IGN) signal(SIGXFSZ, SIG_IGN);
#endif
					/* copy and verify header */
	if (checkheader(infp, HOLOFMT, outfp) < 0 || getw(infp) != HOLOMAGIC)
		error(USER, "input not in holodeck format");
	fputformat(HOLOFMT, outfp);
	fputc('\n', outfp);
	putw(HOLOMAGIC, outfp);
					/* get descriptors and free stdio */
	if ((hfd[0] = dup(fileno(infp))) < 0 ||
			(hfd[1] = dup(fileno(outfp))) < 0)
		error(SYSTEM, "dup call failed in rhinitcopy");
	ifpos = ftell(infp);
	fclose(infp);
	if (fclose(outfp) == EOF)
		error(SYSTEM, "file flushing error in rhinitcopy");
					/* check cache size */
	if (BKBSIZE*1024*1.5 > hdcachesize)
		hdcachesize = BKBSIZE*1024*1.5;
					/* return input position */
	return(ifpos);
}


int
nuniq(rva, n)			/* sort unique rays to front of beam list */
register RAYVAL	*rva;
int	n;
{
	register int	i, j;
	RAYVAL	rtmp;

	for (j = 0; j < n; j++)
		for (i = j+1; i < n; i++)
			if ( rva[i].d == rva[j].d &&
					rva[i].r[0][0]==rva[j].r[0][0] &&
					rva[i].r[0][1]==rva[j].r[0][1] &&
					rva[i].r[1][0]==rva[j].r[1][0] &&
					rva[i].r[1][1]==rva[j].r[1][1] ) {
				n--;		/* swap duplicate with end */
				copystruct(&rtmp, rva+n);
				copystruct(rva+n, rva+i);
				copystruct(rva+i, &rtmp);
				i--;		/* recheck one we swapped */
			}
	return(n);
}


static BEAMI	*beamdir;

static int
bpcmp(b1p, b2p)			/* compare beam positions on disk */
int	*b1p, *b2p;
{
	register off_t	pdif = beamdir[*b1p].fo - beamdir[*b2p].fo;

	if (pdif < 0L) return(-1);
	return(pdif > 0L);
}

static HOLO	*hout;

static int
xferclump(hp, bq, nb)		/* transfer the given clump to hout and free */
HOLO	*hp;
int	*bq, nb;
{
	register int	i;
	register BEAM	*bp;
	int	n;

	beamdir = hp->bi;		/* sort based on file position */
	qsort((void *)bq, nb, sizeof(*bq), bpcmp);
					/* transfer and free each beam */
	for (i = 0; i < nb; i++) {
		bp = hdgetbeam(hp, bq[i]);
		DCHECK(bp==NULL, CONSISTENCY, "empty beam in xferclump");
		n = dupchecking ? nuniq(hdbray(bp),bp->nrm) : bp->nrm;
		memcpy((void *)hdnewrays(hout,bq[i],n),(void *)hdbray(bp), 
				n*sizeof(RAYVAL));
		hdfreebeam(hp, bq[i]);
	}
	hdfreebeam(hout, 0);		/* write & free clump */
	return(0);
}

copysect(ifd, ofd)		/* copy holodeck section from ifd to ofd */
int	ifd, ofd;
{
	HOLO	*hinp;
					/* load input section directory */
	hinp = hdinit(ifd, NULL);
					/* create output section directory */
	hout = hdinit(ofd, (HDGRID *)hinp);
					/* clump the beams */
	clumpbeams(hinp, 0, BKBSIZE*1024, xferclump);
					/* clean up */
	hddone(hinp);
	hddone(hout);
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


void
quit(code)			/* exit the program gracefully */
int	code;
{
	if (tempfile[0])
		unlink(tempfile);
	exit(code);
}
