#ifndef lint
static const char	RCSid[] = "$Id: psketch.c,v 2.4 2018/08/02 18:33:46 greg Exp $";
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
#include  "paths.h"
#include  "resolu.h"
#include  "color.h"

#define	MAXMOD		2048		/* maximum number of modifiers */
#define USESORT		12		/* switch to sorted search */

double	smoothing = 0.8;		/* weight for moving average filter */
double	mixfact = 0.3;			/* amount to mix above/below */

char	*modlist[MAXMOD];		/* (sorted) modifier list */
int	nmods = 0;			/* number of modifiers */

unsigned char  *hasmod;			/* scanline bitmap */
unsigned char  *hasmod1;
COLOR	*scan[2];			/* i/o scanlines */
RESOLU	pres;				/* input resolution */
int	bmwidth;			/* bytes per bitmap */

#undef	tstbit				/* conflicting def's in param.h */
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
	const int	spread = (int)(smoothing*10);
	int		x, width = scanlen(&pres);
	unsigned char	*tbmp;
	COLOR		*tscn;
	char		modbuf[516];
					/* advance buffer */
	tscn = scan[0];
	scan[0] = scan[1];
	scan[1] = tscn;
					/* check if we are at the end */
	if (linesread >= numscans(&pres))
		return(0);
					/* set scanline bitmap */
	if (linesread) {
					/* load & check materials */
		memset(hasmod1, 0, bmwidth);
		for (x = 0; x < width; x++) {
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
				setbit(hasmod1, x);
		}
					/* smear to cover around object */
		memset(hasmod, 0, bmwidth);
		for (x = spread; x < width-spread; x++) {
			int	ox;
			if (!tstbit(hasmod1,x)) continue;
			for (ox = -spread; ox <= spread; ox++)
				setbit(hasmod,x+ox);
		}
	}
					/* read next picture scanline */
	if (freadscan(scan[1], width, infp) < 0) {
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
	hasmod1 = (unsigned char *)malloc(bmwidth);
	scan[0] = (COLOR *)malloc(scanlen(&pres)*sizeof(COLOR));
	scan[1] = (COLOR *)malloc(scanlen(&pres)*sizeof(COLOR));
	if (!hasmod | !hasmod1 | !scan[0] | !scan[1]) {
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
	COLOR		cmavg;

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
	}
	setcolor(cmavg, .0f, .0f, .0f);	/* process this scanline */
	for (x = xstart; x != xstop; x += xstep) {
		COLOR	cmix;
		if (!tstbit(hasmod,x)) continue;
					/* apply moving average */
		scalecolor(scan[0][x], 1.-smoothing);
		scalecolor(cmavg, smoothing);
		addcolor(cmavg, scan[0][x]);
		copycolor(scan[0][x], cmavg);
					/* mix pixel into next scanline */
		copycolor(cmix, scan[0][x]);
		scalecolor(cmix, mixfact);
		scalecolor(scan[1][x], 1.-mixfact);
		addcolor(scan[1][x], cmix);
	}
					/* write out result */
	if (fwritescan(scan[0], width, stdout) < 0) {
		perror("write error");
		return(-1);
	}
	return(1);
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
	free(hasmod1);
	free(scan[0]);
	free(scan[1]);
	return(1);
}

int
main(int argc, char *argv[])
{
	int		i, rval;
	char		pfmt[MAXFMTLEN];
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
		case 's':		/* filter smoothing amount */
			smoothing = atof(argv[++i]);
			if ((smoothing <= FTINY) | (smoothing >= 1.-FTINY)) {
				fprintf(stderr, "%s: smoothing factor must be in (0,1) range\n",
						argv[0]);
				return(1);
			}
			break;
		default:
			fprintf(stderr, "%s: unknown option '%s'\n",
							argv[0], argv[i]);
			return(1);
		}
	if ((argc-i < 2) | (argc-i > 3)) {
		fprintf(stderr, "Usage: %s [-m modname][-M modfile][-s smoothing] octree input.hdr [output.hdr]\n",
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
