#ifndef lint
static const char	RCSid[] = "$Id: rhdisp.c,v 3.50 2004/01/01 11:21:55 schorsch Exp $";
#endif
/*
 * Holodeck display process.
 */

#include <string.h>
#include <ctype.h>

#include "rterror.h"
#include "rholo.h"
#include "rhdisp.h"
#include "rhdriver.h"
#include "selcall.h"

#ifndef VIEWHISTLEN
#define VIEWHISTLEN	4	/* number of remembered views */
#endif

#ifndef FSIZDEF
#define FSIZDEF		0.125	/* default focus frame size */
#endif

#if defined(freebsd)
#define fbufcnt(f)	((f)->_r)
#elif defined(__GNUC__)
#define fbufcnt(f)	((f)->_IO_read_end - (f)->_IO_read_ptr)
#else
#define fbufcnt(f)	((f)->_cnt)
#endif

HOLO	*hdlist[HDMAX+1];	/* global holodeck list */

char	*hdgfn[HDMAX];		/* holodeck section geometry list */
char	*hdpfn[HDMAX];		/* section portal list */

char	cmdlist[DC_NCMDS][8] = DC_INIT;

int	imm_mode = 0;		/* bundles are being delivered immediately */

int	do_outside = 0;		/* render from outside sections */

double	eyesepdist = 1;		/* eye separation distance */

char	*progname;		/* global argv[0] */

FILE	*sstdin, *sstdout;	/* server's standard input and output */

#ifdef DEBUG
extern time_t	time();
static time_t	tmodesw;
static time_t	timm, tadd;
static long	nimmrays, naddrays;
#endif

#define RDY_SRV		01
#define RDY_DEV		02
#define RDY_SIN		04

static int disp_wait(void);
static void add_holo(HDGRID *hdg, char *gfn, char *pfn);
static void disp_bundle(PACKHEAD *p);
static void new_view(register VIEW *v);
static void set_focus(char *args);
static int usr_input(void);
static void printview(void);


int
main(
	int	argc,
	char	*argv[]
)
{
	int	rdy, inp, res = 0, pause = 0;

	progname = argv[0];
	if (argc < 3)
		error(USER, "bad command line arguments");
					/* open our device */
	dev_open(argv[1]);
					/* open server process i/o */
	sstdout = fdopen(atoi(argv[2]), "w");
	if (argc < 4 || (inp = atoi(argv[3])) < 0)
		sstdin = NULL;
	else
		sstdin = fdopen(inp, "r");
					/* set command error vector */
	erract[COMMAND].pf = eputs;
#ifdef DEBUG
	tmodesw = time(NULL);
#endif
					/* enter main loop */
	do {
		rdy = disp_wait();
		if (rdy & RDY_SRV) {		/* process server result */
			res = serv_result();
			if (pause && res != DS_SHUTDOWN) {
				serv_request(DR_ATTEN, 0, NULL);
				while ((res = serv_result()) != DS_ACKNOW &&
						res != DS_SHUTDOWN)
					;
			}
		}
		if (rdy & RDY_DEV) {		/* user input from driver */
			inp = dev_input();
			if (inp & DFL(DC_SETVIEW))
				new_view(&odev.v);
			else if (inp & DFL(DC_LASTVIEW))
				new_view(NULL);
			if (inp & DFL(DC_REDRAW))
				imm_mode = beam_sync(1) > 0;
			if (inp & DFL(DC_GETVIEW))
				printview();
			if (inp & DFL(DC_FOCUS))
				set_focus(odev_args);
			if (inp & DFL(DC_KILL)) {
				serv_request(DR_KILL, 0, NULL);
				pause = 0;
			}
			if (inp & DFL(DC_CLOBBER))
				serv_request(DR_CLOBBER, 0, NULL);
			if (inp & DFL(DC_RESTART)) {
				serv_request(DR_RESTART, 0, NULL);
				pause = 0;
			}
			if (inp & DFL(DC_RESUME)) {
				serv_request(DR_NOOP, 0, NULL);
				pause = 0;
			}
			if (inp & DFL(DC_PAUSE))
				pause = 1;
			if (inp & DFL(DC_QUIT))
				serv_request(DR_SHUTDOWN, 0, NULL);
		}
		if (rdy & RDY_SIN && !imm_mode)	/* user input from sstdin */
			switch (usr_input()) {
			case DC_PAUSE:
				pause = 1;
				break;
			case DC_RESUME:
				serv_request(DR_NOOP, 0, NULL);
				/* fall through */
			case DC_KILL:
			case DC_RESTART:
				pause = 0;
				break;
			}
	} while (res != DS_SHUTDOWN);
#ifdef DEBUG
	if (timm && nimmrays)
		fprintf(stderr,
			"%s: %.1f rays recalled/second (%ld rays total)\n",
				progname, (double)nimmrays/timm, nimmrays);
	if (tadd && naddrays)
		fprintf(stderr,
			"%s: %.1f rays calculated/second (%ld rays total)\n",
				progname, (double)naddrays/tadd, naddrays);
#endif
					/* all done */
	quit(0);
	return 0; /* pro forma return */
}


static int
disp_wait(void)			/* wait for more input */
{
	fd_set	readset, errset;
	int	flgs;
	int	n;
				/* see if we can avoid select call */
	if (hdlist[0] == NULL)
		return(RDY_SRV);	/* initialize first */
	flgs = 0;		/* flag what's ready already */
	if (imm_mode || fbufcnt(stdin) > 0)
		flgs |= RDY_SRV;
	if (sstdin != NULL && fbufcnt(sstdin) > 0)
		flgs |= RDY_SIN;
	if (odev.inpready)
		flgs |= RDY_DEV;
	if (flgs)		/* got something? */
		return(flgs);
	if (dev_flush())	/* else flush output & check keyboard+mouse */
		return(RDY_DEV);
				/* if nothing, we need to call select */
	FD_ZERO(&readset); FD_ZERO(&errset);
	FD_SET(0, &readset);
	FD_SET(0, &errset);
	FD_SET(odev.ifd, &readset);
	FD_SET(odev.ifd, &errset);
	n = odev.ifd+1;
	if (sstdin != NULL) {
		FD_SET(fileno(sstdin), &readset);
		FD_SET(fileno(sstdin), &errset);
		if (fileno(sstdin) >= n)
			n = fileno(sstdin) + 1;
	}
	n = select(n, &readset, NULL, &errset, NULL);
	if (n < 0) {
		if (errno == EINTR)
			return(0);
		error(SYSTEM, "select call failure in disp_wait");
	}
	if (FD_ISSET(0, &readset) || FD_ISSET(0, &errset))
		flgs |= RDY_SRV;
	if (FD_ISSET(odev.ifd, &readset) || FD_ISSET(odev.ifd, &errset))
		flgs |= RDY_DEV;
	if (sstdin != NULL && (FD_ISSET(fileno(sstdin), &readset) ||
				FD_ISSET(fileno(sstdin), &errset)))
		flgs |= RDY_SIN;
	return(flgs);
}


static void
add_holo(		/* register a new holodeck section */
	HDGRID	*hdg,
	char	*gfn,
	char	*pfn
)
{
	VIEW	nv;
	double	d;
	register int	hd;

	for (hd = 0; hd < HDMAX && hdlist[hd] != NULL; hd++)
		;
	if (hd >= HDMAX)
		error(INTERNAL, "too many holodeck sections in add_holo");
	hdlist[hd] = (HOLO *)malloc(sizeof(HOLO));
	if (hdlist[hd] == NULL)
		error(SYSTEM, "out of memory in add_holo");
	memcpy((void *)hdlist[hd], (void *)hdg, sizeof(HDGRID));
	hdcompgrid(hdlist[hd]);
	hdgfn[hd] = savestr(gfn);
	hdpfn[hd] = pfn && *pfn ? savestr(pfn) : (char *)NULL;
	if (hd)
		return;
					/* set initial viewpoint */
	nv = odev.v;
	VSUM(nv.vp, hdlist[0]->orig, hdlist[0]->xv[0], 0.5);
	VSUM(nv.vp, nv.vp, hdlist[0]->xv[1], 0.5);
	VSUM(nv.vp, nv.vp, hdlist[0]->xv[2], 0.5);
	fcross(nv.vdir, hdlist[0]->xv[1], hdlist[0]->xv[2]);
	VCOPY(nv.vup, hdlist[0]->xv[2]);
	if (do_outside) {
		normalize(nv.vdir);
		d = VLEN(hdlist[0]->xv[1]);
		d += VLEN(hdlist[0]->xv[2]);
		VSUM(nv.vp, nv.vp, nv.vdir, -d);
	}
	new_view(&nv);
}


static void
disp_bundle(			/* display a ray bundle */
	register PACKHEAD	*p
)
{
	GCOORD	gc[2];
	FVECT	ro, rd, wp;
	double	d;
	register int	i;
					/* get beam coordinates */
	if ((p->hd < 0) | (p->hd >= HDMAX) || hdlist[p->hd] == NULL)
		error(INTERNAL, "bad holodeck number in disp_bundle");
	if (!hdbcoord(gc, hdlist[p->hd], p->bi))
		error(INTERNAL, "bad beam index in disp_bundle");
					/* display each ray */
	for (i = p->nr; i--; ) {
		hdray(ro, rd, hdlist[p->hd], gc, packra(p)[i].r);
		d = hddepth(hdlist[p->hd], packra(p)[i].d);
		if (d < .99*FHUGE) {
			VSUM(wp, ro, rd, d);	/* might be behind viewpoint */
			dev_value(packra(p)[i].v, rd, wp);
		} else
			dev_value(packra(p)[i].v, rd, NULL);
	}
#ifdef DEBUG
	if (imm_mode) nimmrays += p->nr;
	else naddrays += p->nr;
#endif
}


static void
new_view(			/* change view parameters */
	register VIEW	*v
)
{
	static VIEW	viewhist[VIEWHISTLEN];
	static unsigned	nhist;
	VIEW	*dv;
	int	i, res[2];
	int16	*slist;
	char	*err;
				/* restore previous view? */
	if (v == NULL) {
		if (nhist > 1)		/* get one before last setting */
			nhist--;
		else			/* else go to end of list */
			while (nhist < VIEWHISTLEN && viewhist[nhist].type)
				nhist++;
		v = viewhist + ((nhist-1)%VIEWHISTLEN);
	} else
again:
	if ((err = setview(v)) != NULL) {
		error(COMMAND, err);
		return;
	}
	if (!dev_view(v))	/* notify display driver */
		goto again;
	if (v->type == VT_PAR) {
		error(COMMAND, "cannot handle parallel views");
		return;
	}
	beam_init(odev.firstuse);	/* compute new beam set */
	for (i = 0; (dv = dev_auxview(i, res)) != NULL; i++) {
		if ((slist = beam_view(dv, res[0], res[1])) == NULL) {
			if (!nhist) {
				error(COMMAND, "invalid starting view");
				return;
			}
			*v = *(viewhist + ((nhist-1)%VIEWHISTLEN));
			goto again;	/* poss. overloading dev_section()? */
		}
		DCHECK(*slist < 0, WARNING, "no visible sections in new_view");
		for ( ; *slist >= 0; slist++)
			dev_section(hdgfn[*slist], hdpfn[*slist]);
	}
	dev_section(NULL,NULL);	/* end section list */
	dev_flush();		/* update display */
				/* update server */
	imm_mode = beam_sync(odev.firstuse) > 0;
				/* record new view */
	if (v < viewhist || v >= viewhist+VIEWHISTLEN) {
		*(viewhist + (nhist%VIEWHISTLEN)) = *v;
		nhist++;
	}
}


static void
set_focus(			/* set focus frame */
	char	*args
)
{
	double	hcent, vcent, hsiz, vsiz;
	VIEW	*dv, vwfocus;
	int	i, res[2];

	i = sscanf(args, "%lf %lf %lf %lf", &hcent, &vcent, &hsiz, &vsiz);
	if (i < 2 || hcent < 0 || hcent > 1 || vcent < 0 || vcent > 1) {
		beam_init(0);				/* restore view */
		for (i = 0; (dv = dev_auxview(i, res)) != NULL; i++)
			beam_view(dv, res[0], res[1]);
		beam_sync(0);				/* update server */
		return;
	}
	if (i < 4 || hsiz <= hcent || hsiz > 1 || vsiz <= vcent || vsiz > 1)
		hsiz = vsiz = FSIZDEF;			/* gave center only */
	else {
		hsiz -= hcent; hcent += 0.5*hsiz;	/* gave min and max */
		vsiz -= vcent; vcent += 0.5*vsiz;
	}
	beam_init(0);					/* add basic views */
	for (i = 0; (dv = dev_auxview(i, res)) != NULL; i++)
		beam_view(dv, res[0]>>4, res[1]>>4);
	vwfocus = odev.v;			/* add focus view */
	switch (odev.v.type) {
	case VT_PER:
		vwfocus.horiz = 2.*180./PI*atan(
				hsiz * tan(PI/180./2.*odev.v.horiz) );
		vwfocus.vert = 2.*180./PI*atan(
				vsiz * tan(PI/180./2.*odev.v.vert) );
		break;
	case VT_PAR:
	case VT_ANG:
		vwfocus.horiz = hsiz * odev.v.horiz;
		vwfocus.vert = vsiz * odev.v.vert;
		break;
	case VT_HEM:
		vwfocus.horiz = 2.*180./PI*asin(
				hsiz * sin(PI/180./2.*odev.v.horiz) );
		vwfocus.vert = 2.*180./PI*asin(
				vsiz * sin(PI/180./2.*odev.v.vert) );
		break;
	case VT_CYL:
		vwfocus.horiz = hsiz * odev.v.horiz;
		vwfocus.vert = 2.*180./PI*atan(
				vsiz * tan(PI/180./2.*odev.v.vert) );
		break;
	default:
		error(INTERNAL, "bad view type in set_focus");
	}
	vwfocus.hoff = (odev.v.hoff + hcent - 0.5)/hsiz;
	vwfocus.voff = (odev.v.voff + vcent - 0.5)/vsiz;
	setview(&vwfocus);
	beam_view(&vwfocus, (int)(3*odev.hres*hsiz)+100,
			(int)(3*odev.vres*vsiz)+100);
	beam_sync(0);					/* update server */
}


static int
usr_input(void)			/* get user input and process it */
{
	VIEW	vparams;
	char	cmd[256];
	register char	*args;
	register int	i;

	if (fgets(cmd, sizeof(cmd), sstdin) == NULL) {
		fclose(sstdin);
		sstdin = NULL;
		return(-1);
	}
	if (*cmd == '\n')
		return(DC_RESUME);
	for (args = cmd; *args && !isspace(*args); args++)
		;
	while (isspace(*args))
		*args++ = '\0';
	if (*args && args[i=strlen(args)-1] == '\n')
		args[i] = '\0';
	for (i = 0; i < DC_NCMDS; i++)
		if (!strcmp(cmd, cmdlist[i]))
			break;
	if (i >= DC_NCMDS) {
		dev_auxcom(cmd, args);
		return(-1);
	}
	switch (i) {
	case DC_SETVIEW:		/* set the view */
		vparams = odev.v;
		if (!sscanview(&vparams, args))
			error(COMMAND, "missing view options");
		else
			new_view(&vparams);
		break;
	case DC_GETVIEW:		/* print the current view */
		printview();
		break;
	case DC_LASTVIEW:		/* restore previous view */
		new_view(NULL);
		break;
	case DC_FOCUS:			/* set focus frame */
		set_focus(args);
		break;
	case DC_PAUSE:			/* pause the current calculation */
	case DC_RESUME:			/* resume the calculation */
		/* handled in main() */
		break;
	case DC_REDRAW:			/* redraw from server */
		imm_mode = beam_sync(1) > 0;
		dev_clear();
		break;
	case DC_KILL:			/* kill rtrace process(es) */
		serv_request(DR_KILL, 0, NULL);
		break;
	case DC_CLOBBER:		/* clobber holodeck */
		serv_request(DR_CLOBBER, 0, NULL);
		break;
	case DC_RESTART:		/* restart rtrace */
		serv_request(DR_RESTART, 0, NULL);
		break;
	case DC_QUIT:			/* quit request */
		serv_request(DR_SHUTDOWN, 0, NULL);
		break;
	default:
		error(CONSISTENCY, "bad command id in usr_input");
	}
	return(i);
}


static void
printview(void)			/* print our current view to server stdout */
{
	fputs(VIEWSTR, sstdout);
	fprintview(&odev.v, sstdout);
	fputc('\n', sstdout);
	fflush(sstdout);
}


extern int
serv_result(void)			/* get next server result and process it */
{
	static char	*buf = NULL;
	static int	bufsiz = 0;
	MSGHEAD	msg;
					/* read message header */
	if (fread((char *)&msg, sizeof(MSGHEAD), 1, stdin) != 1)
		goto readerr;
	if (msg.nbytes > 0) {		/* get the message body */
		if (msg.nbytes > bufsiz) {
			if (buf == NULL)
				buf = (char *)malloc(bufsiz=msg.nbytes);
			else
				buf = (char *)realloc((void *)buf,
						bufsiz=msg.nbytes);
			if (buf == NULL)
				error(SYSTEM, "out of memory in serv_result");
		}
		if (fread(buf, 1, msg.nbytes, stdin) != msg.nbytes)
			goto readerr;
	}
	switch (msg.type) {		/* process results */
	case DS_BUNDLE:
		if (msg.nbytes < sizeof(PACKHEAD) ||
				msg.nbytes != packsiz(((PACKHEAD *)buf)->nr))
			error(INTERNAL, "bad display packet from server");
		disp_bundle((PACKHEAD *)buf);
		break;
	case DS_ADDHOLO:
		if (msg.nbytes < sizeof(HDGRID)+2)
			error(INTERNAL, "bad holodeck record from server");
		add_holo((HDGRID *)buf, buf+sizeof(HDGRID),
			buf+sizeof(HDGRID)+strlen(buf+sizeof(HDGRID))+1);
		break;
	case DS_OUTSECT:
		do_outside = 1;
		goto noargs;
	case DS_EYESEP:
		if (msg.nbytes <= 1 || (eyesepdist = atof(buf)) <= FTINY)
			error(INTERNAL, "bad eye separation from server");
		break;
	case DS_STARTIMM:
	case DS_ENDIMM:
		if (!(imm_mode = msg.type==DS_STARTIMM))
			dev_flush();
#ifdef DEBUG
		{
			time_t	tnow = time(NULL);
			if (msg.type==DS_STARTIMM) tadd += tnow - tmodesw;
			else timm += tnow - tmodesw;
			tmodesw = tnow;
		}
#endif
		goto noargs;
	case DS_ACKNOW:
	case DS_SHUTDOWN:
		goto noargs;
	default:
		error(INTERNAL, "unrecognized result from server process");
	}
	return(msg.type);		/* return message type */
noargs:
	if (msg.nbytes) {
		sprintf(errmsg, "unexpected body with server message %d",
				msg.type);
		error(INTERNAL, errmsg);
	}
	return(msg.type);
readerr:
	if (feof(stdin))
		error(SYSTEM, "server process died");
	error(SYSTEM, "error reading from server process");
	return -1;  
}


extern void
serv_request(	/* send a request to the server process */
	int	type,
	int	nbytes,
	char	*p
)
{
	MSGHEAD	msg;
	int	m;
				/* consistency checks */
	DCHECK(nbytes < 0 || nbytes > 0 & p == NULL,
			CONSISTENCY, "bad buffer handed to serv_request");
				/* get server's attention for big request */
	if (nbytes >= BIGREQSIZ-sizeof(MSGHEAD)) {
		serv_request(DR_ATTEN, 0, NULL);
		while ((m = serv_result()) != DS_ACKNOW)
			if (m == DS_SHUTDOWN)	/* the bugger quit on us */
				quit(0);
	}
	msg.type = type;	/* write and flush the message */
	msg.nbytes = nbytes;
	fwrite((char *)&msg, sizeof(MSGHEAD), 1, stdout);
	if (nbytes > 0)
		fwrite(p, 1, nbytes, stdout);
	if (fflush(stdout) < 0)
		error(SYSTEM, "write error in serv_request");
}


void
eputs(			/* put error message to stderr */
	register char  *s
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
quit(			/* clean up and exit */
	int	code
)
{
	if (code)
		exit(code);
	if (odev.v.type)
		dev_close();
	exit(0);
}
