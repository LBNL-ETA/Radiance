/* RCSid $Id: rcontrib.h,v 2.5 2012/06/16 17:09:49 greg Exp $ */

/*
 * Header file for rcontrib modules
 */

#include "ray.h"
#include "func.h"
#include "lookup.h"

extern int		gargc;		/* global argc */
extern char		**gargv;	/* global argv */
extern char		*octname;	/* global octree name */

extern int		nproc;		/* number of processes requested */
extern int		nchild;		/* number of children (-1 in child) */

extern int		inpfmt;		/* input format */
extern int		outfmt;		/* output format */

extern int		header;		/* output header? */
extern int		force_open;	/* truncate existing output? */
extern int		recover;	/* recover previous output? */
extern int		accumulate;	/* input rays per output record */
extern int		contrib;	/* computing contributions? */

extern int		xres;		/* horizontal (scan) size */
extern int		yres;		/* vertical resolution */

extern int		using_stdout;	/* are we using stdout? */

extern int		imm_irrad;	/* compute immediate irradiance? */
extern int		lim_dist;	/* limit distance? */

extern int		account;	/* current accumulation count */
extern RNUMBER		raysleft;	/* number of rays left to trace */
extern long		waitflush;	/* how long until next flush */

extern RNUMBER		lastray;	/* last ray number sent */
extern RNUMBER		lastdone;	/* last ray processed */

typedef double		DCOLOR[3];	/* double-precision color */

/*
 * The MODCONT structure is used to accumulate ray contributions
 * for a particular modifier, which may be subdivided into bins
 * if binv evaluates > 0.  If outspec contains a %s in it, this will
 * be replaced with the modifier name.  If outspec contains a %d in it,
 * this will be used to create one output file per bin, otherwise all bins
 * will be written to the same file, in order.  If the global outfmt
 * is 'c', then a 4-byte RGBE pixel will be output for each bin value
 * and the file will conform to a RADIANCE image if xres & yres are set.
 */
typedef struct {
	const char	*outspec;	/* output file specification */
	const char	*modname;	/* modifier name */
	EPNODE		*binv;		/* bin value expression */
	int		nbins;		/* number of contribution bins */
	DCOLOR		cbin[1];	/* contribution bins (extends struct) */
} MODCONT;			/* modifier contribution */

extern LUTAB		modconttab;	/* modifier contribution table */

/*
 * The STREAMOUT structure holds an open FILE pointer and a count of
 * the number of RGB triplets per record, or 0 if unknown.
 */
typedef struct {
	FILE		*ofp;		/* output file pointer */
	int		outpipe;	/* output is to a pipe */
	int		reclen;		/* triplets/record */
	int		xr, yr;		/* output resolution for picture */
} STREAMOUT;

extern LUTAB		ofiletab;	/* output stream table */

#ifndef MAXPROCESS
#ifdef _WIN32
#define MAXPROCESS	1
#else
#define MAXPROCESS	128
#endif
#endif

#ifndef	MAXMODLIST
#define	MAXMODLIST	2048		/* maximum modifiers we'll track */
#endif

extern const char	*modname[MAXMODLIST];	/* ordered modifier name list */
extern int		nmods;			/* number of modifiers */

extern char		RCCONTEXT[];		/* special evaluation context */

extern char		*formstr(int f);	/* return format identifier */

extern void		process_rcontrib(void);	/* trace ray contributions */

extern STREAMOUT *	getostream(const char *ospec, const char *mname,
							int bn, int noopen);

extern void		mod_output(MODCONT *mp);
extern void		end_record(void);

extern MODCONT		*addmodifier(char *modn, char *outf,
						char *binv, int bincnt);
extern void		addmodfile(char *fname, char *outf,
						char *binv, int bincnt);

extern void		reload_output(void);
extern void		recover_output(void);

extern int		getvec(FVECT vec);

extern int		in_rchild(void);
extern void		end_children(void);

extern void		put_zero_record(int ndx);

extern void		parental_loop(void);	/* controlling process */

extern void		rcontrib(void);		/* main calculation loop */
