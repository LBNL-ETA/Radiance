#ifndef lint
static const char	RCSid[] = "$Id: rxpiece.cpp,v 2.6 2024/10/31 19:22:36 greg Exp $";
#endif
/*
 *  rxpiece.cpp - main for rxpiece tile rendering program
 */

#include "copyright.h"

#include  <time.h>
#include  <signal.h>
#include  <sys/mman.h>
#include  <unistd.h>

#include  "platform.h"
#include  "RpictSimulManager.h"
#include  "ambient.h"
#include  "pmapray.h"
#include  "random.h"

extern char  *progname;			/* argv[0] */
const char  *sigerr[NSIG];		/* signal error messages */

VIEW  ourview = STDVIEW;		/* global view parameters */
int  hresolu = 1024;			/* horizontal resolution */
int  vresolu = 1024;			/* vertical resolution */
double	pixaspect = 1.0;		/* pixel aspect ratio */
int  hres, vres;			/* current image resolution for srcdraw.c */

int	tileGrid[2] = {5,5};		// tile subdivisions

int  psample = 4;			/* pixel sample size */
double	maxdiff = .05;			/* max. difference for interpolation */
double	dstrpix = 0.67;			/* square pixel distribution */

double  mblur = 0.;			/* motion blur parameter (unused) */

double  dblur = 0.;			/* depth-of-field blur parameter */

int  nproc = 1;				/* number of processes to run */

RpictSimulManager	myRPmanager;	// global simulation manager

static void onsig(int signo);
static void onalrm(int signo);
static void sigdie(int  signo, const char  *msg);
static void printdefaults(void);
static RenderDataType rpiece(char *pout, RenderDataType dt, char *zout);

					/* rxpiece additional features */
#define RXPIECE_FEATURES	"Recovery\nIrradianceCalc\nViewTypes=v,l,a,h,s,c\n" \
		"ParticipatingMedia=Mist\n" \
		"HessianAmbientCache\nAmbientAveraging\nAmbientValueSharing\n" \
		"PixelJitter\nPixelSampling\nPixelDepthOfField\n" \
		"SmallSourceDrawing\n" \
		"AdaptiveShadowTesting\nOutputs=v,l\n" \
		"OutputCS=RGB,XYZ,prims,spec\n"


// We could call myRPmanager.Cleanup() but why waste time
// unwinding data structures when the whole frame is going away?
void
quit(int code)				/* quit program */
{
	ambsync();			// flush ambient cache

	ray_done_pmap();		/* PMAP: free photon maps */

	exit(code);
}


int
main(int  argc, char  *argv[])
{
#define	 check(ol,al)		if (argv[i][ol] || \
				badarg(argc-i-1,argv+i+1,al)) \
				goto badopt
#define	 check_bool(olen,var)		switch (argv[i][olen]) { \
				case '\0': var = !var; break; \
				case 'y': case 'Y': case 't': case 'T': \
				case '+': case '1': var = 1; break; \
				case 'n': case 'N': case 'f': case 'F': \
				case '-': case '0': var = 0; break; \
				default: goto badopt; }
	RGBPRIMS  our_prims;		/* private output color primitives */
	RenderDataType	dtype = RDTrgbe;	// output data flags
	char  *outfile = NULL;
	char  *zfile = NULL;
	int  outfmt = 'c';
	int  rval;
	int  i;
					/* global program name */
	progname = argv[0];
					/* feature check only? */
	strcat(RFeatureList, RXPIECE_FEATURES);
	if (argc > 1 && !strcmp(argv[1], "-features"))
		return feature_status(argc-2, argv+2);
					/* option city */
	for (i = 1; i < argc; i++) {
						/* expand arguments */
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", argv[i]);
			error(SYSTEM, errmsg);
		}
		if (argv[i] == NULL || argv[i][0] != '-')
			break;			/* break from options */
		if (!strcmp(argv[i], "-version")) {
			puts(VersionID);
			quit(0);
		}
		if (!strcmp(argv[i], "-defaults") ||
				!strcmp(argv[i], "-help")) {
			printdefaults();
			quit(0);
		}
		rval = getrenderopt(argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
		rval = getviewopt(&ourview, argc-i, argv+i);
		if (rval >= 0) {
			i += rval;
			continue;
		}
						/* rxpiece options */
		switch (argv[i][1]) {
		case 'v':				/* view file */
			if (argv[i][2] != 'f')
				goto badopt;
			check(3,"s");
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
			break;
		case 'n':				/* number of processes */
			check(2,"i");
			nproc = atoi(argv[++i]);
			if (nproc < 0 && (nproc += RadSimulManager::GetNCores()) <= 0)
				nproc = 1;
			break;
		case 'f':				/* output format */
			if ((argv[i][2] != 'c') & (argv[i][2] != 'f')
					|| argv[i][3])
				goto badopt;
			outfmt = argv[i][2];
			break;
		case 'p':				/* pixel */
			switch (argv[i][2]) {
			case 's':				/* sample */
				check(3,"i");
				psample = atoi(argv[++i]);
				if (psample < 1) psample = 1;
				break;
			case 't':				/* threshold */
				check(3,"f");
				maxdiff = atof(argv[++i]);
				break;
			case 'j':				/* jitter */
				check(3,"f");
				dstrpix = atof(argv[++i]);
				break;
			case 'a':				/* aspect */
				check(3,"f");
				pixaspect = atof(argv[++i]);
				break;
			case 'd':				/* aperture */
				check(3,"f");
				dblur = atof(argv[++i]);
				dblur *= (dblur > 0);
				break;
			case 'R':				/* standard RGB output */
				if (strcmp(argv[i]+2, "RGB"))
					goto badopt;
				myRPmanager.prims = stdprims;
				dtype = RDTnewCT(dtype, RDTrgbe);
				break;
			case 'X':				/* XYZ output */
				if (strcmp(argv[i]+2, "XYZ"))
					goto badopt;
				myRPmanager.prims = xyzprims;
				dtype = RDTnewCT(dtype, RDTxyze);
				break;
			case 'c':				/* chromaticities */
				check(3,"ffffffff");
				rval = 0;
				for (int j = 0; j < 8; j++) {
					our_prims[0][j] = atof(argv[++i]);
					rval |= fabs(our_prims[0][j]-stdprims[0][j]) > .001;
				}
				if (rval) {
					if (!colorprimsOK(our_prims))
						error(USER, "illegal primary chromaticities");
					myRPmanager.prims = our_prims;
				} else
					myRPmanager.prims = stdprims;
				dtype = RDTnewCT(dtype, RDTrgbe);
				break;
			default:
				goto badopt;
			}
			break;
		case 'd':				/* reference depth */
			if (argv[i][2] || !myRPmanager.SetReferenceDepth(argv[++i]))
				goto badopt;
			dtype = RDTnewDT(dtype, RDTdshort);
			break;
		case 'x':				/* x resolution */
			check(2,"i");
			hresolu = atoi(argv[++i]);
			break;
		case 'y':				/* y resolution */
			check(2,"i");
			vresolu = atoi(argv[++i]);
			break;
		case 'X':				/* horizontal tile subdivisions */
			check(2,"i");
			tileGrid[0] = atoi(argv[++i]);
			break;
		case 'Y':				/* vertical tile subdivisions */
			check(2,"i");
			tileGrid[1] = atoi(argv[++i]);
			break;
		case 'o':				/* output file */
			check(2,"s");
			outfile = argv[++i];
			break;
		case 'z':				/* z file */
			check(2,"s");
			zfile = argv[++i];
			break;
#if MAXCSAMP>3
		case 'c':				/* output spectral results */
			if (argv[i][2] != 'o')
				goto badopt;
			rval = (myRPmanager.prims == NULL);
			check_bool(3,rval);
			if (rval)
				myRPmanager.prims = NULL;
			else if (myRPmanager.prims == NULL)
				myRPmanager.prims = stdprims;
			dtype = RDTnewCT(dtype, rval ? RDTscolr : RDTrgbe);
			break;
#endif
		case 'w':				/* warnings */
			rval = erract[WARNING].pf != NULL;
			check_bool(2,rval);
			if (rval) erract[WARNING].pf = wputs;
			else erract[WARNING].pf = NULL;
			break;
		default:
			goto badopt;
		}
	}
	if (maxdiff <= FTINY)		/* check for useless sampling */
		psample = 1;
	if (outfile == NULL)
		error(USER, "missing output file (-o option)");
	if (zfile == NULL)		/* set up depth output */
		dtype = RDTnewDT(dtype, RDTnone);
	else if (!RDTdepthT(dtype))
		dtype = RDTnewDT(dtype, RDTdfloat);
					/* check pixel output type */
	if ((myRPmanager.prims == NULL) & (NCSAMP == 3)) {
		myRPmanager.prims = stdprims;
		dtype = RDTnewCT(dtype, RDTrgbe);
	}
	if (outfmt == 'f')
		switch (RDTcolorT(dtype)) {
		case RDTrgbe:
			dtype = RDTnewCT(dtype, RDTrgb);
			break;
		case RDTxyze:
			dtype = RDTnewCT(dtype, RDTxyz);
			break;
		case RDTscolr:
			dtype = RDTnewCT(dtype, RDTscolor);
			break;
		case RDTrgb:
		case RDTxyz:
		case RDTscolor:
			break;
		default:
			error(INTERNAL, "botched color output type");
		}
					/* set up signal handling */
	sigdie(SIGINT, "Interrupt");
	sigdie(SIGHUP, "Hangup");
	sigdie(SIGTERM, "Terminate");
	sigdie(SIGPIPE, "Broken pipe");
	signal(SIGALRM, onalrm);	// used to gracefully terminate
#ifdef	SIGXCPU
	sigdie(SIGXCPU, "CPU limit exceeded");
	sigdie(SIGXFSZ, "File size exceeded");
#endif
#ifdef	NICE
	nice(NICE);			/* lower priority */
#endif
	if (i < argc-1)
		goto badopt;
					// load octree
	if (!myRPmanager.LoadOctree(argv[i]))
		error(USER, "missing octree argument");
					// add new header info
	myRPmanager.AddHeader(i, argv);
	{
		char	buf[128] = "SOFTWARE= ";
		strcpy(buf+10, VersionID);
		myRPmanager.AddHeader(buf);
	}
					// render tiles
	dtype = rpiece(outfile, dtype, zfile);

	quit(dtype==RDTnone);		// status is 1 on failure

badopt:
	sprintf(errmsg, "command line error at '%s'", argv[i]);
	error(USER, errmsg);
	return 1; /* pro forma return */

#undef	check
#undef	check_bool
}


void
wputs(				/* warning output function */
	const char	*s
)
{
	if (!erract[WARNING].pf)
		return;		// warnings were disabled!
	int  lasterrno = errno;
	eputs(s);
	errno = lasterrno;
}


void
eputs(				/* put string to stderr */
	const char  *s
)
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


static void
onsig(				/* fatal signal */
	int  signo
)
{
	static int  gotsig = 0;

	if (gotsig++)			/* two signals and we're gone! */
		_exit(signo);

	alarm(30);			/* allow 30 seconds to clean up */
	signal(SIGALRM, SIG_DFL);	/* make certain we do die */
	eputs("signal - ");
	eputs(sigerr[signo]);
	eputs("\n");
	quit(3);
}


static bool	gotALRM = false;	// flag for ALRM signal

static void
onalrm(int signo)
{
	gotALRM = true;
}


static void
sigdie(			/* set fatal signal */
	int  signo,
	const char  *msg
)
{
	if (signal(signo, onsig) == SIG_IGN)
		signal(signo, SIG_IGN);
	sigerr[signo] = msg;
}


static void
printdefaults(void)			/* print default values to stdout */
{
	printf("-n %-2d\t\t\t\t# number of rendering processes\n", nproc);
	printf("-vt%c\t\t\t\t# view type %s\n", ourview.type,
			ourview.type==VT_PER ? "perspective" :
			ourview.type==VT_PAR ? "parallel" :
			ourview.type==VT_HEM ? "hemispherical" :
			ourview.type==VT_ANG ? "angular" :
			ourview.type==VT_CYL ? "cylindrical" :
			ourview.type==VT_PLS ? "planisphere" :
			"unknown");
	printf("-vp %f %f %f\t# view point\n",
			ourview.vp[0], ourview.vp[1], ourview.vp[2]);
	printf("-vd %f %f %f\t# view direction\n",
			ourview.vdir[0], ourview.vdir[1], ourview.vdir[2]);
	printf("-vu %f %f %f\t# view up\n",
			ourview.vup[0], ourview.vup[1], ourview.vup[2]);
	printf("-vh %f\t\t\t# view horizontal size\n", ourview.horiz);
	printf("-vv %f\t\t\t# view vertical size\n", ourview.vert);
	printf("-vo %f\t\t\t# view fore clipping plane\n", ourview.vfore);
	printf("-va %f\t\t\t# view aft clipping plane\n", ourview.vaft);
	printf("-vs %f\t\t\t# view shift\n", ourview.hoff);
	printf("-vl %f\t\t\t# view lift\n", ourview.voff);
	printf("-x  %-9d\t\t\t# x resolution\n", hresolu);
	printf("-y  %-9d\t\t\t# y resolution\n", vresolu);
	printf("-X  %-9d\t\t\t# horizontal tile divisions\n", tileGrid[0]);
	printf("-Y  %-9d\t\t\t# vertical tile divisions\n", tileGrid[1]);
	if (myRPmanager.prims == stdprims)
		printf("-pRGB\t\t\t\t# standard RGB color output\n");
	else if (myRPmanager.prims == xyzprims)
		printf("-pXYZ\t\t\t\t# CIE XYZ color output\n");
	else if (myRPmanager.prims != NULL)
		printf("-pc %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\t# output color primaries and white point\n",
				myRPmanager.prims[RED][0], myRPmanager.prims[RED][1],
				myRPmanager.prims[GRN][0], myRPmanager.prims[GRN][1],
				myRPmanager.prims[BLU][0], myRPmanager.prims[BLU][1],
				myRPmanager.prims[WHT][0], myRPmanager.prims[WHT][1]);
	if (NCSAMP > 3)
		printf(myRPmanager.prims != NULL ? "-co-\t\t\t\t# output tristimulus colors\n" :
				"-co+\t\t\t\t# output spectral values\n");
	printf("-pa %f\t\t\t# pixel aspect ratio\n", pixaspect);
	printf("-pj %f\t\t\t# pixel jitter\n", dstrpix);
	printf("-pd %f\t\t\t# pixel depth-of-field\n", dblur);
	printf("-ps %-9d\t\t\t# pixel sample\n", psample);
	printf("-pt %f\t\t\t# pixel threshold\n", maxdiff);
	printf(erract[WARNING].pf != NULL ?
			"-w+\t\t\t\t# warning messages on\n" :
			"-w-\t\t\t\t# warning messages off\n");
	print_rdefaults();
}


// Struct for tracking tiles being rendered in mapped file / shared memory
struct TileProg {
	short		status;		// 0==Unstarted, -1==InProgress, 1==Done
	pid_t		pID;		// process operating on tile
} *tprog = NULL;		// shared tile progress array

#define tile_p(ti)	(tprog + (ti)[1]*tileGrid[0] + (ti)[0])

// Return true if tile is renderable
static bool
renderable_tile(TileProg *tp)
{
	if (tp->status < 0 && kill(tp->pID, 0) < 0)
		tp->status = 0;		// dead process - reset

	return !tp->status;
}


// handle multi-processing if requested, return true if all done
static bool
children_finished()
{
	if (nproc <= 1)			// single process -> run in parent
		return false;
	int	cnt = 0;		// else count ready-to-go tiles
	int	ti[2];
	for (ti[1] = 0; ti[1] < tileGrid[1]; ti[1]++)
		for (ti[0] = 0; ti[0] < tileGrid[0]; ti[0]++)
			cnt += renderable_tile(tile_p(ti));
	if (!cnt)
		return false;		// parent can do nothing
	if (cnt < nproc) {
		sprintf(errmsg, "only %d renderable tiles, reducing process count", cnt);
		error(WARNING, errmsg);
		if ((nproc = cnt) == 1)
			return false;	// back to single process
	}
	cow_memshare();			// else we'll be sharing memory
	fflush(NULL);			// and forking children
	pid_t	cpid;			// create nproc children
	for (cnt = nproc; cnt && (cpid = fork()) != 0; cnt--)
		if (cpid < 0)
			error(SYSTEM, "fork error!");

	if (cpid == 0) {		// children render tiles
		sleep(nproc - cnt);	// avoid race conditions
		return false;
	}
	cow_doneshare();		// parent frees memory and waits
	signal(SIGALRM, SIG_IGN);
	myRPmanager.Cleanup(true);
	int	nfailed = 0;
	int	status;
	for (cnt = nproc; cnt && wait(&status) > 0; cnt--)
		if (status) {
			sprintf(errmsg, "child exited with status %d", status);
			error(WARNING, errmsg);
			if (!nfailed++ & (cnt > 1)) {
				kill(0, SIGALRM);
				error(WARNING, "waiting for other tiles to finish...");
			}
		}
	if (cnt) {
		sprintf(errmsg, "lost track of %d children", cnt);
		error(WARNING, errmsg);
	}
	if (nfailed) {
		sprintf(errmsg, "%d tiles were not completed", nfailed);
		error(USER, errmsg);
	}
	return true;			// all done!
}

// return next renderable tile, false if everything is done
static bool
nexttile(int ti[2])
{
	static pid_t	ourpID = 0;
	static short *	tlist = NULL;
	static int	tlen = 0;
	static int	tnext = 0;

	if (gotALRM) {			// pre-empting new work?
		if (tlist) {
			sprintf(errmsg, "process %d got alarm, exiting", ourpID);
			CHECK(tnext<tlen, WARNING, errmsg);
			free(tlist); tlist = NULL;
		}
		return false;
	}
	if (!tlist) {			// initialize random tile list
		ABitMap2	todoMap(tileGrid[0], tileGrid[1]);
		tlen = 0;
		for (ti[1] = 0; ti[1] < tileGrid[1]; ti[1]++)
			for (ti[0] = 0; ti[0] < tileGrid[0]; ti[0]++)
				if (renderable_tile(tile_p(ti))) {
					todoMap.Set(ti[0], ti[1]);
					tlen++;
				}
		if (!tlen)
			return false;	// nothing to do!
		tlist = (short *)malloc(sizeof(short)*2*tlen);
		CHECK(!tlist, SYSTEM, "out of memory in nexttile()");
		tlen = 0;		// assign entries
		for (ti[0] = ti[1] = 0; todoMap.Find(&ti[0], &ti[1]); ti[0]++) {
			tlist[2*tlen] = ti[0];
			tlist[2*tlen+1] = ti[1];
			tlen++;
		}
				// shuffle order w/ Fisher-Yates
		for (int i = 0; i < tlen-1; i++) {
			const int	ix = irandom(tlen-i) + i;
			ti[0] = tlist[2*i];
			ti[1] = tlist[2*i+1];
			tlist[2*i] = tlist[2*ix];
			tlist[2*i+1] = tlist[2*ix+1];
			tlist[2*ix] = ti[0];
			tlist[2*ix+1] = ti[1];
		}
		ourpID = getpid();	// save time on system calls
	}
	while (tnext < tlen) {		// find first available
		ti[0] = tlist[2*tnext];
		ti[1] = tlist[2*tnext+1];
		tnext++;		// take if still unclaimed
		if (renderable_tile(tile_p(ti))) {
			tile_p(ti)->status = -1;
			tile_p(ti)->pID = ourpID;
			return true;
		}
	}
	free(tlist); tlist = NULL;	// exhausted our list?
	return false;
}


// Principal function for rpiece
static RenderDataType
rpiece(char *pout, RenderDataType dt, char *zout)
{
	if (zout && *zout == '!')
		error(USER, "cannot send depth to a command");

	const bool	newOutput = (access(pout, F_OK) < 0);
	FILE		*pdfp[2];
	if (newOutput) {			// new output file?
		CHECK((tileGrid[0] <= 1) & (tileGrid[1] <= 1),
				 USER, "bad tiling specification");
	} else {
		dt = myRPmanager.ReopenOutput(pdfp, pout, zout);
		if (dt == RDTnone)
			quit(1);
		if (!fscnresolu(&hresolu, &vresolu, pdfp[0]))
			error(USER, "missing picture resolution");
		pixaspect = .0;			// need to leave this as is
		myRPmanager.NewHeader(pout);	// get prev. header info
		const char *	tval = myRPmanager.GetHeadStr("TILED=");
		if (tval) sscanf(tval, "%d %d", &tileGrid[0], &tileGrid[1]);
		CHECK(myRPmanager.GetView()==NULL,
				USER, "missing view in picture file");
		ourview = *myRPmanager.GetView();
	}
	int	hvdim[2] = {hresolu, vresolu};	// set up tiled frame
	if (!myRPmanager.NewFrame(ourview, hvdim, &pixaspect, tileGrid))
		error(USER, "tiling setup error in rpiece");

	if ((hvdim[0] != hresolu) | (hvdim[1] != vresolu)) {
		if (!newOutput)
			error(USER, "unexpected output size adjustment");
		sprintf(errmsg, "resolution adjusted from %dx%d to %dx%d",
				hresolu, vresolu, hvdim[0], hvdim[1]);
		error(WARNING, errmsg);
		hresolu = hvdim[0];
		vresolu = hvdim[1];
	}
	if (newOutput){				// open new output here
		char	buf[64];
		sprintf(buf, "TILED= %d %d\n", tileGrid[0], tileGrid[1]);
		myRPmanager.AddHeader(buf);
		dt = myRPmanager.NewOutput(pdfp, pout, dt, zout);
		if (dt == RDTnone)
			quit(1);
		fprtresolu(hresolu, vresolu, pdfp[0]);
		fflush(pdfp[0]);
		if (RDTdepthT(dt) == RDTdshort) {
			fprtresolu(hresolu, vresolu, pdfp[1]);
			fflush(pdfp[1]);
		}
	} else if (RDTdepthT(dt) == RDTdshort &&
			(!fscnresolu(&hvdim[0], &hvdim[1], pdfp[1]) ||
				(hvdim[0] != hresolu) | (hvdim[1] != vresolu)))
		error(USER, "mismatched depth file resolution");
						// prepare (flat) pixel buffer
	const long	pdata_beg = ftell(pdfp[0]);
	const size_t	pixSiz = (RDTcolorT(dt)==RDTrgbe)|(RDTcolorT(dt)==RDTxyze) ? sizeof(COLR)
			: (RDTcolorT(dt)==RDTrgb)|(RDTcolorT(dt)==RDTxyz) ? sizeof(COLORV)*3
			: RDTcolorT(dt)==RDTscolr ? LSCOLR : sizeof(COLORV)*NCSAMP;
	size_t		pmlen = pdata_beg + pixSiz*hresolu*vresolu;
						// put tile progress array at end
	if (pmlen&7) pmlen += 8 - (pmlen&7);	// 8-byte alignment to be safe
	pmlen += sizeof(TileProg)*tileGrid[0]*tileGrid[1];
						// map picture file to memory
	if (newOutput && ftruncate(fileno(pdfp[0]), pmlen) < 0)
		error(SYSTEM, "cannot extend picture buffer");
	uby8 *		pixMap = (uby8 *)mmap(NULL, pmlen, PROT_READ|PROT_WRITE,
						MAP_SHARED, fileno(pdfp[0]), 0);
	if ((void *)pixMap == MAP_FAILED)
		error(SYSTEM, "cannot map picture file into memory");
						// map depth buffer to memory
	const long	zdata_beg = RDTdepthT(dt) ? ftell(pdfp[1]) : 0L;
	const size_t	zdpSiz = RDTdepthT(dt)==RDTdshort ? sizeof(short) :
				RDTdepthT(dt)==RDTdfloat ? sizeof(float) : 0;
	const size_t	zmlen = zdata_beg + zdpSiz*hresolu*vresolu;
	uby8 *		zdMap = NULL;
	if (RDTdepthT(dt)) {
		if (newOutput && ftruncate(fileno(pdfp[1]), zmlen) < 0)
			error(SYSTEM, "cannot extend depth buffer");
		zdMap = (uby8 *)mmap(NULL, zmlen, PROT_READ|PROT_WRITE,
						MAP_SHARED, fileno(pdfp[1]), 0);
		if ((void *)zdMap == MAP_FAILED)
			error(SYSTEM, "cannot map depth file into memory");
	}
	fclose(pdfp[0]);			// done with file pointers
	if (RDTdepthT(dt)) fclose(pdfp[1]);
						// point to tile progress array
	tprog = (TileProg *)(pixMap + pmlen - sizeof(TileProg)*tileGrid[0]*tileGrid[1]);

	if (children_finished())		// work done in children?
		return dt;

	int	ndone = 0;			// else render tiles
	int	ti[2];
	while (nexttile(ti)) {
		const int	offset = (tileGrid[1]-1-ti[1])*myRPmanager.GetWidth()*myRPmanager.THeight() +
						(myRPmanager.THeight()-1)*myRPmanager.GetWidth() +
						ti[0]*myRPmanager.TWidth();
		uby8 *		pptr = pixMap + pdata_beg + pixSiz*offset;
		uby8 *		zptr = zdMap + zdata_beg + zdpSiz*offset;
		bool		ok = false;
		switch (RDTcommonE(dt)<<1 | (RDTdepthT(dt)==RDTdshort)) {
		case 2:		// common-exponent color, float/no depth
			ok = myRPmanager.RenderTile((COLRV *)pptr, -myRPmanager.GetWidth(),
							(float *)zptr, ti);
			break;
		case 0:		// float color, float/no depth
			ok = myRPmanager.RenderTile((COLORV *)pptr, -myRPmanager.GetWidth(),
							(float *)zptr, ti);
			break;
		case 3:		// common-exponent color, encoded depth
			ok = myRPmanager.RenderTile((COLRV *)pptr, -myRPmanager.GetWidth(),
							(short *)zptr, ti);
			break;
		case 1:		// float color, encoded depth
			ok = myRPmanager.RenderTile((COLORV *)pptr, -myRPmanager.GetWidth(),
							(short *)zptr, ti);
			break;
		}
		if (!ok) {			// got an error
			sprintf(errmsg, "error rendering tile (%d,%d)/(%d,%d)",
					ti[0], ti[1], tileGrid[0], tileGrid[1]);
			error(USER, errmsg);
		}
		tile_p(ti)->status = 1;		// mark tile completed
		ndone++;
	}
	if (!ndone)
		error(WARNING, "no tiles need rendering, exit");
	/*
	munmap(pixMap, pmlen);			// technically unnecessary...
	if (zdMap) munmap(zdMap, zmlen);
	*/
	return dt;				// we're done here
}
