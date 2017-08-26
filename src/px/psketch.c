#ifndef lint
static const char	RCSid[] = "$Id: psketch.c,v 2.1 2017/08/26 16:07:22 greg Exp $";
#endif
/*
 *  psketch.c - modify picture to sketch objects with named modifiers
 *
 *	9/23/2017
 */

#include  "copyright.h"

#include  <string.h>

#include  "standard.h"
#include  "platform.h"
#include  "resolu.h"
#include  "color.h"
#include  "random.h"
					/* our probabilities */
#define PROB_LEFT1	0.1
#define PROB_RIGHT1	0.1
#define PROB_LEFT2	0.2
#define PROB_RIGHT2	0.2
#define PROB_LEFT3	0.07
#define PROB_RIGHT3	0.07
#define PROB_DOWN1	0.1

#define	MAXMOD		2048		/* maximum number of modifiers */
#define USESORT		12		/* switch to sorted search */

char	*modlist[MAXMOD];		/* (sorted) modifier list */
int	nmods = 0;			/* number of modifiers */

unsigned char  *hasmod, *moved[2];	/* scanline bitmaps */
COLR	*scan[2];			/* i/o scanlines */
RESOLU	pres;				/* input resolution */
int	bmwidth;			/* bytes per bitmap */
					/* conflicting def's in param.h */
#undef	tstbit
#undef	setbit
#undef	clrbit
#undef	tglbit

#define	 bitop(bm,x,op)		((bm)[(x)>>3] op (1<<((x)&7)))
#define	 tstbit(bm,x)		bitop(bm,x,&)
#define	 setbit(bm,x)		bitop(bm,x,|=)
#define	 clrbit(bm,x)		bitop(bm,x,&=~)
#define	 tglbit(bm,x)		bitop(bm,x,^=)

const char	*octree, *infname;	/* input octree and picture names */
FILE		*infp;			/* picture input stream */
FILE		*mafp;			/* pipe from rtrace with modifiers */
int		linesread = 0;		/* scanlines read so far */

static int
find_mod(const char *s)
{
	int     ilower, iupper;
	int	c, i;

	if (nmods < USESORT) {			/* linear search */
		for (i = nmods; i-- > 0; )
			if (!strcmp(s, modlist[i]))
				return(i);
		return(-1);
	}
	ilower = 0; iupper = nmods;
	c = iupper;				/* binary search */
	while ((i = (iupper + ilower) >> 1) != c) {
		c = strcmp(s, modlist[i]);
		if (c > 0)
			ilower = i;
		else if (c < 0)
			iupper = i;
		else
			return(i);
		c = i;
	}
	return(-1);
}

static int
read_scan(void)
{
	int		x, width = scanlen(&pres);
	unsigned char	*tbmp;
	COLR		*tscn;
	char		modbuf[516];
					/* advance buffers */
	tbmp = moved[0];
	moved[0] = moved[1];
	moved[1] = tbmp;
	tscn = scan[0];
	scan[0] = scan[1];
	scan[1] = tscn;
					/* check if we are at the end */
	if (linesread >= numscans(&pres))
		return(0);
					/* clear bitmaps */
	if (linesread)
		memset(hasmod, 0, bmwidth);
	memset(moved[1], 0, bmwidth);
					/* load & check materials */
	for (x = 0; x < width*(linesread>0); x++) {
		int	len;
		if (fgets(modbuf, sizeof(modbuf), mafp) == NULL) {
			fprintf(stderr, "Error reading from rtrace!\n");
			return(-1);
		}
		len = strlen(modbuf);
		if (len < 3 || (modbuf[len-1] != '\n') |
					(modbuf[len-2] != '\t')) {
			fprintf(stderr, "Garbled rtrace output: %s", modbuf);
			return(-1);
		}
		modbuf[len-2] = '\0';
		if (find_mod(modbuf) >= 0)
			setbit(hasmod, x);
	}
					/* read next picture scanline */
	if (freadcolrs(scan[1], width, infp) < 0) {
		fprintf(stderr, "%s: error reading scanline %d\n",
					infname, linesread);
		return(-1);
	}
	return(++linesread);
}

static int
spcmp(const void *p1, const void *p2)
{
	return strcmp(*(const char **)p1, *(const char **)p2);
}

static int
get_started(void)
{
	char	combuf[1024];

	if (nmods >= USESORT)		/* need to sort for search? */
		qsort(modlist, nmods, sizeof(char *), &spcmp);
					/* open pipe from rtrace */
	sprintf(combuf, "vwrays -fd %s | rtrace -h -fda -om %s",
			infname, octree);
	mafp = popen(combuf, "r");
	if (mafp == NULL) {
		perror("popen");
		return(0);
	}
					/* allocate bitmaps & buffers */
	bmwidth = (scanlen(&pres)+7) >> 3;
	hasmod = (unsigned char *)malloc(bmwidth);
	moved[0] = (unsigned char *)malloc(bmwidth);
	moved[1] = (unsigned char *)malloc(bmwidth);
	scan[0] = (COLR *)malloc(scanlen(&pres)*sizeof(COLR));
	scan[1] = (COLR *)malloc(scanlen(&pres)*sizeof(COLR));
	if (!hasmod | !moved[0] | !moved[1] | !scan[0] | !scan[1]) {
		perror("malloc");
		return(0);
	}
	if (!read_scan())		/* read first scanline */
		return(0);
	return(1);
}

static int
advance_scanline(void)
{
	static int	alldone = 0;
	int		width = scanlen(&pres);
	int		height = numscans(&pres);
	int		x, xstart = 0, xstop = width, xstep = 1;
	COLR		tclr;

#define CSWAP(y0,x0,y1,x1) { copycolr(tclr,scan[y0][x0]); \
		copycolr(scan[y1][x1],scan[y0][x0]); \
		copycolr(scan[y0][x0],tclr); \
		setbit(moved[y0],x0); setbit(moved[y1],x1); }

	if (alldone)			/* finished last scanline? */
		return(0);
	alldone = (linesread >= height);
					/* advance scanline */
	if (read_scan() < !alldone)
		return(-1);
	if (linesread & 1) {		/* process in horizontal zig-zag */
		xstart = width-1;
		xstop = -1;
		xstep = -1;
	}				/* process this scanline */
	for (x = xstart; x != xstop; x += xstep) {
		double	runif = frandom();
		if (tstbit(moved[0],x) || !tstbit(hasmod,x))
			continue;
		if (runif < PROB_LEFT1) {
			if (x > 0 && !tstbit(moved[0],x-1))
				CSWAP(0,x,0,x-1);
			continue;
		}
		runif -= PROB_LEFT1;
		if (runif < PROB_RIGHT1) {
			if (x < width-1 && !tstbit(moved[0],x+1))
				CSWAP(0,x,0,x+1);
			continue;
		}
		runif -= PROB_RIGHT1;
		if (runif < PROB_LEFT2) {
			if (x > 1 && !tstbit(moved[0],x-2))
				CSWAP(0,x,0,x-2);
			continue;
		}
		runif -= PROB_LEFT2;
		if (runif < PROB_RIGHT2) {
			if (x < width-2 && !tstbit(moved[0],x+2))
				CSWAP(0,x,0,x+2);
			continue;
		}
		runif -= PROB_RIGHT2;
		if (runif < PROB_LEFT3) {
			if (x > 2 && !tstbit(moved[0],x-3))
				CSWAP(0,x,0,x-3);
			continue;
		}
		runif -= PROB_LEFT3;
		if (runif < PROB_RIGHT3) {
			if (x < width-3 && !tstbit(moved[0],x+3))
				CSWAP(0,x,0,x+3);
			continue;
		}
		runif -= PROB_RIGHT3;
		if (runif < PROB_DOWN1) {
			if (linesread < height-1)
				CSWAP(0,x,1,x);
			continue;
		}
		runif -= PROB_DOWN1;
	}
					/* write it out */
	if (fwritecolrs(scan[0], width, stdout) < 0) {
		perror("write error");
		return(-1);
	}
	return(1);
#undef CSWAP
}

static int
clean_up(void)
{
	char	linebuf[128];
					/* read unused materials */
	while (fgets(linebuf, sizeof(linebuf), mafp) != NULL)
		;
	if (pclose(mafp) != 0) {
		fprintf(stderr, "Error running rtrace!\n");
		return(0);
	}
	fclose(infp);
	free(hasmod);
	free(moved[0]);
	free(moved[1]);
	free(scan[0]);
	free(scan[1]);
	return(1);
}

int
main(int argc, char *argv[])
{
	int		i, rval;
	char		pfmt[LPICFMT+1];
					/* process options */
	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'm':		/* new modifier name */
			if (nmods >= MAXMOD) {
				fprintf(stderr, "%s: too many modifiers\n", argv[0]);
				return(1);
			}
			modlist[nmods++] = argv[++i];
			break;
		case 'M':		/* modifier file */
			rval = wordfile(modlist, MAXMOD-nmods, argv[++i]);
			if (rval < 0) {
				fprintf(stderr, "%s: cannot open modifier file '%s'\n",
							argv[0], argv[i]);
				return(1);
			}
			nmods += rval;
			break;
		default:
			fprintf(stderr, "%s: unknown option '%s'\n",
							argv[0], argv[i]);
			return(1);
		}
	if ((argc-i < 2) | (argc-i > 3)) {
		fprintf(stderr, "Usage: %s [-m modname][-M modfile] octree input.hdr [output.hdr]\n",
						argv[0]);
		return(1);
	}
	if (!nmods) {
		fprintf(stderr, "%s: at least one '-m' or '-M' option needed\n",
					argv[0]);
		return(1);
	}
	octree = argv[i];
	infname = argv[i+1];		/* open input picture */
	infp = fopen(infname, "rb");
	if (infp == NULL) {
		fprintf(stderr, "%s: cannot open input '%s' for reading\n",
				argv[0], argv[i+1]);
		return(1);
	}
	if (i+2 < argc &&		/* open output picture */
		freopen(argv[i+2], "w", stdout) == NULL) {
			fprintf(stderr, "%s: cannot open output '%s' for writing\n",
					argv[0], argv[i+2]);
			return(1);
	}
	SET_FILE_BINARY(stdout);
	strcpy(pfmt, PICFMT);		/* copy format/resolution */
	if (checkheader(infp, pfmt, stdout) < 0 || !fgetsresolu(&pres, infp)) {
		fprintf(stderr, "%s: bad format for input picture\n",
				argv[0]);
		return(1);
	}
	printargs(argc, argv, stdout);	/* complete header */
	fputformat(pfmt, stdout);
	fputc('\n', stdout);
	fputsresolu(&pres, stdout);	/* resolution does not change */
	if (!get_started())
		return(1);
	while ((rval = advance_scanline()))
		if (rval < 0)
			return(1);
	if (!clean_up())
		return(1);
	return(0);
}
