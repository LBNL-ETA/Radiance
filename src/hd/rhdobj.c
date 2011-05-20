#ifndef lint
static const char	RCSid[] = "$Id: rhdobj.c,v 3.20 2011/05/20 02:06:39 greg Exp $";
#endif
/*
 * Routines for loading and displaying Radiance objects in rholo with GLX.
 */

#include <string.h>
#include <ctype.h>

#include "radogl.h"
#include "tonemap.h"
#include "rhdisp.h"
#include "rhdriver.h"
#include "rhdobj.h"
#include "rtprocess.h"

extern FILE	*sstdout;		/* user standard output */

char	rhdcmd[DO_NCMDS][8] = DO_INIT;	/* user command list */

/* pointer to function to get lights */
void	(*dobj_lightsamp)(COLR clr, FVECT direc, FVECT pos) = NULL;

#define	AVGREFL		0.5		/* assumed average reflectance */

#define	MAXAC		512		/* maximum number of args */

#ifndef MINTHRESH
#define MINTHRESH	5.0		/* source threshold w.r.t. mean */
#endif

#ifndef NALT
#define NALT		11		/* # sampling altitude angles */
#endif
#ifndef NAZI
#define NAZI		17
#endif

typedef struct dlights {
	struct dlights	*next;		/* next in lighting set list */
	FVECT	lcent;			/* computed lighting center */
	double	ravg;			/* harmonic mean free range radius */
	TMbright	larb;		/* average reflected brightness */
	COLOR	lamb;			/* local ambient value */
	short	nl;			/* number of lights in this set */
	struct lsource {
		FVECT	direc;			/* source direction */
		double	omega;			/* source solid angle */
		COLOR	val;			/* source color */
	} li[MAXLIGHTS];		/* light sources */
} DLIGHTS;			/* a light source set */

typedef struct dobject {
	struct dobject	*next;		/* next object in list */
	char	name[64];		/* object name */
	FVECT	center;			/* orig. object center */
	RREAL	radius;			/* orig. object radius */
	int	listid;			/* GL display list identifier */
	int	nlists;			/* number of lists allocated */
	SUBPROC	rtp;			/* associated rtrace process */
	DLIGHTS	*ol;			/* object lights */
	FULLXF	xfb;			/* coordinate transform */
	short	drawcode;		/* drawing code */
	short	xfac;			/* transform argument count */
	char	*xfav[MAXAC+1];		/* transform args */
} DOBJECT;			/* a displayable object */

static DLIGHTS	*dlightsets;		/* lighting sets */
static DOBJECT	*dobjects;		/* display object list */
static DOBJECT	*curobj;		/* current (last referred) object */
static int	lastxfac;		/* last number of transform args */
static char	*lastxfav[MAXAC+1];	/* saved transform arguments */

#define getdcent(c,op)	multp3(c,(op)->center,(op)->xfb.f.xfm)
#define getdrad(op)	((op)->radius*(op)->xfb.f.sca)

#define	RTARGC	8
static char	*rtargv[RTARGC+1] = {"rtrace", "-h-", "-w-", "-fdd",
					"-x", "1", "-oL"};

static struct {
	int	nsamp;			/* number of ray samples */
	COLOR	val;			/* value (sum) */
} ssamp[NALT][NAZI];		/* current sphere samples */

#define	curname		(curobj==NULL ? (char *)NULL : curobj->name)

static DOBJECT *getdobj(char *nm);
static int freedobj(DOBJECT *op);
static int savedxf(DOBJECT *op);
static void ssph_sample(COLR clr, FVECT direc, FVECT pos);
static void ssph_direc(FVECT direc, int alt, int azi);
static int ssph_neigh(int sp[2], int next);
static int ssph_compute(void);
static int getdlights(DOBJECT *op, int force);
static void cmderror(int cn, char *err);


static DOBJECT *
getdobj(			/* get object from list by name */
	char	*nm
)
{
	register DOBJECT	*op;

	if (nm == NULL)
		return(NULL);
	if (nm == curname)
		return(curobj);
	for (op = dobjects; op != NULL; op = op->next)
		if (!strcmp(op->name, nm))
			break;
	return(op);
}


static int
freedobj(			/* free resources and memory assoc. with op */
	register DOBJECT	*op
)
{
	int	foundlink = 0;
	DOBJECT	ohead;
	register DOBJECT	*opl;

	if (op == NULL)
		return(0);
	ohead.next = dobjects;
	for (opl = &ohead; opl->next != NULL; opl = opl->next) {
		if (opl->next == op && (opl->next = op->next) == NULL)
			break;
		foundlink += opl->next->listid == op->listid;
	}
	dobjects = ohead.next;
	if (!foundlink) {
		glDeleteLists(op->listid, op->nlists);
		close_process(&(op->rtp));
	}
	while (op->xfac)
		freestr(op->xfav[--op->xfac]);
	free((void *)op);
	return(1);
}


static int
savedxf(			/* save transform for display object */
	register DOBJECT	*op
)
{
					/* free old */
	while (lastxfac)
		freestr(lastxfav[--lastxfac]);
					/* nothing to save? */
	if (op == NULL) {
		lastxfav[0] = NULL;
		return(0);
	}
					/* else save new */
	for (lastxfac = 0; lastxfac < op->xfac; lastxfac++)
		lastxfav[lastxfac] = savestr(op->xfav[lastxfac]);
	lastxfav[lastxfac] = NULL;
	return(1);
}


static void
ssph_sample(	/* add sample to current source sphere */
	COLR	clr,
	FVECT	direc,
	FVECT	pos
)
{
	COLOR	col;
	double	d;
	register int	alt, azi;

	if (dlightsets == NULL)
		return;
	if (pos == NULL)
		d = FHUGE;		/* sample is at infinity */
	else if ((d = (pos[0] - dlightsets->lcent[0])*direc[0] +
			(pos[1] - dlightsets->lcent[1])*direc[1] +
			(pos[2] - dlightsets->lcent[2])*direc[2]) > FTINY)
		dlightsets->ravg += 1./d;
	else
		return;			/* sample is behind us */
	alt = NALT*(1.-FTINY)*(.5*(1.+FTINY) + direc[2]*.5);
	azi = NAZI*(1.-FTINY)*(.5*(1.+FTINY) + .5/PI*atan2(direc[1],direc[0]));
	colr_color(col, clr);
	addcolor(ssamp[alt][azi].val, col);
	ssamp[alt][azi].nsamp++;
}


static void
ssph_direc(	/* compute sphere sampling direction */
	FVECT	direc,
	int	alt,
	int	azi
)
{
	double	phi, d;

	direc[2] = 2./NALT*(alt+.5) - 1.;
	d = sqrt(1. - direc[2]*direc[2]);
	phi = 2.*PI/NAZI*(azi+.5) - PI;
	direc[0] = d*cos(phi);
	direc[1] = d*sin(phi);
}


static int
ssph_neigh(		/* neighbor counter on sphere */
	register int	sp[2],
	int	next
)
{
	static short	nneigh = 0;		/* neighbor count */
	static short	neighlist[NAZI+6][2];	/* neighbor list (0 is home) */
	register int	i;

	if (next) {
		if (nneigh <= 0)
			return(0);
		sp[0] = neighlist[--nneigh][0];
		sp[1] = neighlist[nneigh][1];
		return(1);
	}
	if ((sp[0] < 0) | (sp[0] >= NALT) | (sp[1] < 0) | (sp[1] >= NAZI))
		return(nneigh=0);
	neighlist[0][0] = sp[0]; neighlist[0][1] = sp[1];
	nneigh = 1;
	if (sp[0] == 0) {
		neighlist[nneigh][0] = 1;
		neighlist[nneigh++][1] = (sp[1]+1)%NAZI;
		neighlist[nneigh][0] = 1;
		neighlist[nneigh++][1] = sp[1];
		neighlist[nneigh][0] = 1;
		neighlist[nneigh++][1] = (sp[1]+(NAZI-1))%NAZI;
		for (i = 0; i < NAZI; i++)
			if (i != sp[1]) {
				neighlist[nneigh][0] = 0;
				neighlist[nneigh++][1] = i;
			}
	} else if (sp[0] == NALT-1) {
		neighlist[nneigh][0] = NALT-2;
		neighlist[nneigh++][1] = (sp[1]+1)%NAZI;
		neighlist[nneigh][0] = NALT-2;
		neighlist[nneigh++][1] = sp[1];
		neighlist[nneigh][0] = NALT-2;
		neighlist[nneigh++][1] = (sp[1]+(NAZI-1))%NAZI;
		for (i = 0; i < NAZI; i++)
			if (i != sp[1]) {
				neighlist[nneigh][0] = NALT-1;
				neighlist[nneigh++][1] = i;
			}
	} else {
		neighlist[nneigh][0] = sp[0]-1;
		neighlist[nneigh++][1] = (sp[1]+1)%NAZI;
		neighlist[nneigh][0] = sp[0]-1;
		neighlist[nneigh++][1] = sp[1];
		neighlist[nneigh][0] = sp[0]-1;
		neighlist[nneigh++][1] = (sp[1]+(NAZI-1))%NAZI;
		neighlist[nneigh][0] = sp[0];
		neighlist[nneigh++][1] = (sp[1]+1)%NAZI;
		neighlist[nneigh][0] = sp[0];
		neighlist[nneigh++][1] = (sp[1]+(NAZI-1))%NAZI;
		neighlist[nneigh][0] = sp[0]+1;
		neighlist[nneigh++][1] = (sp[1]+1)%NAZI;
		neighlist[nneigh][0] = sp[0]+1;
		neighlist[nneigh++][1] = sp[1];
		neighlist[nneigh][0] = sp[0]+1;
		neighlist[nneigh++][1] = (sp[1]+(NAZI-1))%NAZI;
	}
	return(nneigh);
}


static int
ssph_compute(void)			/* compute source set from sphere samples */
{
	int	ncells, nsamps;
	COLOR	csum;
	FVECT	v;
	double	d, thresh, maxbr;
	int	maxalt, maxazi, spos[2];
	register int	alt, azi;
	register struct lsource	*ls;
					/* count & average sampled cells */
	setcolor(csum, 0., 0., 0.);
	ncells = nsamps = 0;
	for (alt = 0; alt < NALT; alt++)
		for (azi = 0; azi < NAZI; azi++)
			if (ssamp[alt][azi].nsamp) {
				if (ssamp[alt][azi].nsamp > 1) {
					d = 1.0/ssamp[alt][azi].nsamp;
					scalecolor(ssamp[alt][azi].val, d);
				}
				addcolor(csum, ssamp[alt][azi].val);
				nsamps += ssamp[alt][azi].nsamp;
				ncells++;
			}
	if ((dlightsets == NULL) | (ncells < NALT*NAZI/4)) {
		ncells = 0;
		goto done;
	}
						/* harmonic mean distance */
	if (dlightsets->ravg > FTINY)
		dlightsets->ravg = nsamps / dlightsets->ravg;
	else
		dlightsets->ravg = FHUGE;
						/* light source threshold */
	thresh = MINTHRESH*bright(csum)/ncells;
	if (thresh <= FTINY) {
		ncells = 0;
		goto done;
	}
						/* avg. reflected brightness */
	d = AVGREFL / (double)ncells;	
	scalecolor(csum, d);
	if (tmCvColors(tmGlobal, &dlightsets->larb,
			TM_NOCHROM, &csum, 1) != TM_E_OK)
		error(CONSISTENCY, "tone mapping problem in ssph_compute");
					/* greedy light source clustering */
	while (dlightsets->nl < MAXLIGHTS) {
		maxbr = 0.;			/* find brightest cell */
		for (alt = 0; alt < NALT; alt++)
			for (azi = 0; azi < NAZI; azi++)
				if ((d = bright(ssamp[alt][azi].val)) > maxbr) {
					maxalt = alt; maxazi = azi;
					maxbr = d;
				}
		if (maxbr < thresh)		/* below threshold? */
			break;
		ls = dlightsets->li + dlightsets->nl++;
		spos[0] = maxalt; spos[1] = maxazi;	/* cluster */
		for (ssph_neigh(spos, 0); ssph_neigh(spos, 1); ) {
			alt = spos[0]; azi = spos[1];
			if ((d = bright(ssamp[alt][azi].val)) < .75*thresh)
				continue;		/* too dim */
			ssph_direc(v, alt, azi);	/* else add it in */
			VSUM(ls->direc, ls->direc, v, d);
			ls->omega += 1.;
			addcolor(ls->val, ssamp[alt][azi].val);
							/* remove from list */
			setcolor(ssamp[alt][azi].val, 0., 0., 0.);
			ssamp[alt][azi].nsamp = 0;
		}
		d = 1./ls->omega;			/* avg. brightness */
		scalecolor(ls->val, d);
		ls->omega *= 4.*PI/(NALT*NAZI);		/* solid angle */
		normalize(ls->direc);			/* direction */
	}
					/* compute ambient remainder */
	for (alt = 0; alt < NALT; alt++)
		for (azi = 0; azi < NAZI; azi++)
			if (ssamp[alt][azi].nsamp)
				addcolor(dlightsets->lamb, ssamp[alt][azi].val);
	d = 1.0/ncells;
	scalecolor(dlightsets->lamb, d);
done:					/* clear sphere sample array */
	memset((void *)ssamp, '\0', sizeof(ssamp));
	return(ncells);
}


static int
getdlights(		/* get lights for display object */
	register DOBJECT	*op,
	int	force
)
{
	double	d2, mind2 = FHUGE*FHUGE;
	FVECT	ocent;
	VIEW	cvw;
	register DLIGHTS	*dl;

	op->ol = NULL;
	if (op->drawcode != DO_LIGHT)
		return(0);
					/* check for usable light set */
	getdcent(ocent, op);
	for (dl = dlightsets; dl != NULL; dl = dl->next)
		if ((d2 = dist2(dl->lcent, ocent)) < mind2) {
			op->ol = dl;
			mind2 = d2;
		}
					/* the following is heuristic */
	d2 = 2.*getdrad(op); d2 *= d2;
	if ((dl = op->ol) != NULL && (mind2 < 0.0625*dl->ravg*dl->ravg ||
			mind2 < 4.*getdrad(op)*getdrad(op)))
		return(1);
	if (!force)
		return(0);
					/* need to compute new light set */
	cvw = stdview;
	cvw.type = VT_PER;
	VCOPY(cvw.vp, ocent);
	cvw.vup[0] = 1.; cvw.vup[1] = cvw.vup[2] = 0.;
	cvw.horiz = 90; cvw.vert = 90.;
	beam_init(1);			/* query beams through center */
	cvw.vdir[0] = cvw.vdir[1] = 0.; cvw.vdir[2] = 1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
	cvw.vdir[0] = cvw.vdir[1] = 0.; cvw.vdir[2] = -1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
	cvw.vup[0] = cvw.vup[1] = 0.; cvw.vup[2] = 1.;
	cvw.vdir[0] = cvw.vdir[2] = 0.; cvw.vdir[1] = 1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
	cvw.vdir[0] = cvw.vdir[2] = 0.; cvw.vdir[1] = -1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
	cvw.vdir[1] = cvw.vdir[2] = 0.; cvw.vdir[0] = 1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
	cvw.vdir[1] = cvw.vdir[2] = 0.; cvw.vdir[0] = -1.;
	setview(&cvw); beam_view(&cvw, 0, 0);
					/* allocate new light set */
	dl = (DLIGHTS *)calloc(1, sizeof(DLIGHTS));
	if (dl == NULL)
		goto memerr;
	VCOPY(dl->lcent, ocent);
					/* push onto our light set list */
	dl->next = dlightsets;
	dlightsets = dl;
	dobj_lightsamp = ssph_sample;	/* get beams from server */
	imm_mode = beam_sync(-1) > 0;
	while (imm_mode)
		if (serv_result() == DS_SHUTDOWN)
			quit(0);
	if (!ssph_compute()) {		/* compute light sources from sphere */
		dlightsets = dl->next;
		free((void *)dl);
		return(0);
	}
	op->ol = dl;
	return(1);
memerr:
	error(SYSTEM, "out of memory in getdlights");
	return 0; /* pro forma return */
}


static void
cmderror(		/* report command error */
	int	cn,
	char	*err
)
{
	sprintf(errmsg, "%s: %s", rhdcmd[cn], err);
	error(COMMAND, errmsg);
}


extern int
dobj_command(		/* run object display command */
	char	*cmd,
	register char	*args
)
{
	int	somechange = 0;
	int	cn, na;
	register int	nn;
	char	*alist[MAXAC+1], *nm;
					/* find command */
	for (cn = 0; cn < DO_NCMDS; cn++)
		if (!strcmp(cmd, rhdcmd[cn]))
			break;
	if (cn >= DO_NCMDS)
		return(-1);		/* not in our list */
					/* make argument list */
	for (na = 0; *args; na++) {
		if (na > MAXAC)
			goto toomany;
		alist[na] = args;
		while (*args && !isspace(*args))
			args++;
		while (isspace(*args))
			*args++ = '\0';
	}
	alist[na] = NULL;
					/* execute command */
	switch (cn) {
	case DO_LOAD:				/* load an octree */
		if (na == 1)
			dobj_load(alist[0], alist[0]);
		else if (na == 2)
			dobj_load(alist[0], alist[1]);
		else {
			cmderror(cn, "need octree [name]");
			return(0);
		}
		break;
	case DO_UNLOAD:				/* clear an object */
		if (na > 1) goto toomany;
		if (na && alist[0][0] == '*')
			somechange += dobj_cleanup();
		else
			somechange += dobj_unload(na ? alist[0] : curname);
		break;
	case DO_XFORM:				/* transform object */
	case DO_MOVE:
		if (na && alist[0][0] != '-') {
			nm = alist[0]; nn = 1;
		} else {
			nm = curname; nn = 0;
		}
		if (cn == DO_MOVE && nn >= na) {
			cmderror(cn, "missing transform");
			return(0);
		}
		somechange += dobj_xform(nm, cn==DO_MOVE, na-nn, alist+nn);
		break;
	case DO_UNMOVE:				/* undo last transform */
		somechange += dobj_unmove();
		break;
	case DO_OBJECT:				/* print object statistics */
		if (dobj_putstats(na ? alist[0] : curname, sstdout))
			if (na && alist[0][0] != '*' && (curobj == NULL ||
					strcmp(alist[0], curobj->name)))
				savedxf(curobj = getdobj(alist[0]));
		break;
	case DO_DUP:				/* duplicate object */
		for (nn = 0; nn < na; nn++)
			if (alist[nn][0] == '-')
				break;
		switch (nn) {
		case 0:
			cmderror(cn, "need new object name");
			return(0);
		case 1:
			nm = curname;
			break;
		case 2:
			nm = alist[0];
			break;
		default:
			goto toomany;
		}
		if (!dobj_dup(nm, alist[nn-1]))
			break;
		if (na > nn)
			somechange += dobj_xform(curname, 1, na-nn, alist+nn);
		else
			curobj->drawcode = DO_HIDE;
		savedxf(curobj);
		break;
	case DO_SHOW:				/* change rendering option */
	case DO_LIGHT:
	case DO_HIDE:
		if (na > 1) goto toomany;
		dobj_lighting(na ? alist[0] : curname, cn);
		somechange++;
		break;
	default:
		error(CONSISTENCY, "bad command id in dobj_command");
	}
	return(somechange);
toomany:
	cmderror(cn, "too many arguments");
	return(-1);
}


extern int
dobj_load(		/* create/load an octree object */
	char	*oct,
	char	*nam
)
{
	char	*fpp, fpath[128];
	register DOBJECT	*op;
					/* check arguments */
	if (oct == NULL) {
		error(COMMAND, "missing octree");
		return(0);
	}
	if (nam == NULL) {
		error(COMMAND, "missing name");
		return(0);
	}
	if ((*nam == '*') | (*nam == '-')) {
		error(COMMAND, "illegal name");
		return(0);
	}
	if (getdobj(nam) != NULL) {
		error(COMMAND, "name already taken (clear first)");
		return(0);
	}
					/* get octree path */
	if ((fpp = getpath(oct, getrlibpath(), R_OK)) == NULL) {
		sprintf(errmsg, "cannot find octree \"%s\"", oct);
		error(COMMAND, errmsg);
		return(0);
	}
	strcpy(fpath, fpp);
	op = (DOBJECT *)malloc(sizeof(DOBJECT));
	if (op == NULL)
		error(SYSTEM, "out of memory in dobj_load");
					/* set struct fields */
	strcpy(op->name, nam);
	op->ol = NULL;
	op->drawcode = DO_HIDE;
	setident4(op->xfb.f.xfm); op->xfb.f.sca = 1.;
	setident4(op->xfb.b.xfm); op->xfb.b.sca = 1.;
	op->xfav[op->xfac=0] = NULL;
					/* load octree into display list */
	dolights = 0;
	domats = 1;
	op->listid = rgl_octlist(fpath, op->center, &op->radius, &op->nlists);
					/* start rtrace */
	rtargv[RTARGC-1] = fpath;
	rtargv[RTARGC] = NULL;
	open_process(&(op->rtp), rtargv);
					/* insert into main list */
	op->next = dobjects;
	curobj = dobjects = op;
	savedxf(NULL);
	return(1);
}


extern int
dobj_unload(			/* free the named object */
	char	*nam
)
{
	register DOBJECT	*op;

	if ((op = getdobj(nam)) == NULL) {
		error(COMMAND, "no object");
		return(0);
	}
	freedobj(op);
	savedxf(curobj = NULL);
	return(1);
}


extern int
dobj_cleanup(void)				/* free all resources */
{
	register DLIGHTS	*lp;

	while (dobjects != NULL)
		freedobj(dobjects);
	savedxf(curobj = NULL);
	while ((lp = dlightsets) != NULL) {
		dlightsets = lp->next;
		free((void *)lp);
	}
	return(1);
}


extern int
dobj_xform(		/* set/add transform for nam */
	char	*nam,
	int	rel,
	int	ac,
	char	**av
)
{
	register DOBJECT	*op;
	FVECT	cent;
	double	rad;
	char	scoord[16];
	int	i;

	if ((op = getdobj(nam)) == NULL) {
		error(COMMAND, "no object");
		return(0);
	}
	if (rel)
		rel = op->xfac + 8;
	if (ac + rel > MAXAC) {
		error(COMMAND, "too many transform arguments");
		return(0);
	}
	savedxf(curobj = op);		/* remember current transform */
	if (rel && ac == 4 && !strcmp(av[0], "-t"))
		rel = -1;			/* don't move for translate */
	else {
		getdcent(cent, op);		/* don't move if near orig. */
		rad = getdrad(op);
		if (DOT(cent,cent) < rad*rad)
			rel = -1;
	}
	if (!rel) {				/* remove old transform */
		while (op->xfac)
			freestr(op->xfav[--op->xfac]);
	} else if (rel > 0) {			/* relative move */
		op->xfav[op->xfac++] = savestr("-t");
		for (i = 0; i < 3; i++) {
			sprintf(scoord, "%.4e", -cent[i]);
			op->xfav[op->xfac++] = savestr(scoord);
		}
	}
	while (ac--)
		op->xfav[op->xfac++] = savestr(*av++);
	if (rel > 0) {				/* move back */
		op->xfav[op->xfac++] = savestr("-t");
		for (i = 0; i < 3; i++) {
			sprintf(scoord, "%.4e", cent[i]);
			op->xfav[op->xfac++] = savestr(scoord);
		}
	}
	op->xfav[op->xfac] = NULL;
	if (fullxf(&op->xfb, op->xfac, op->xfav) != op->xfac) {
		error(COMMAND, "bad transform arguments");
		dobj_unmove();
		savedxf(op);		/* save current transform instead */
		return(0);
	}
					/* don't know local lights anymore */
	getdlights(op, 0);
	return(1);
}


extern int
dobj_putstats(			/* put out statistics for nam */
	char	*nam,
	FILE	*fp
)
{
	FVECT	ocent;
	register DOBJECT	*op;
	register int	i;

	if (nam == NULL) {
		error(COMMAND, "no current object");
		return(0);
	}
	if (nam[0] == '*') {
		i = 0;
		for (op = dobjects; op != NULL; op = op->next)
			i += dobj_putstats(op->name, fp);
		return(i);
	}
	if ((op = getdobj(nam)) == NULL) {
		error(COMMAND, "unknown object");
		return(0);
	}
	getdcent(ocent, op);
	fprintf(fp, "%s: %s, center [%g %g %g], radius %g", op->name,
			op->drawcode==DO_HIDE ? "hidden" :
			op->drawcode==DO_LIGHT && op->ol!=NULL ? "lighted" :
			"shown",
			ocent[0],ocent[1],ocent[2], getdrad(op));
	if (op->xfac)
		fputs(", (xform", fp);
	for (i = 0; i < op->xfac; i++) {
		putc(' ', fp);
		fputs(op->xfav[i], fp);
	}
	if (op->xfac)
		fputc(')', fp);
	fputc('\n', fp);
	return(1);
}


extern int
dobj_unmove(void)				/* undo last transform change */
{
	int	txfac;
	char	*txfav[MAXAC+1];

	if (curobj == NULL) {
		error(COMMAND, "no current object");
		return(0);
	}
					/* hold last transform */
	memcpy((void *)txfav, (void *)lastxfav, 
			(txfac=lastxfac)*sizeof(char *));
					/* save this transform */
	memcpy((void *)lastxfav, (void *)curobj->xfav, 
			(lastxfac=curobj->xfac)*sizeof(char *));
					/* copy back last transform */
	memcpy((void *)curobj->xfav, (void *)txfav, 
			(curobj->xfac=txfac)*sizeof(char *));
					/* set matrices */
	fullxf(&curobj->xfb, curobj->xfac, curobj->xfav);
					/* don't know local lights anymore */
	getdlights(curobj, 0);
	return(1);
}


extern int
dobj_dup(			/* duplicate object oldnm as nam */
	char	*oldnm,
	char	*nam
)
{
	register DOBJECT	*op, *opdup;
					/* check arguments */
	if ((op = getdobj(oldnm)) == NULL) {
		error(COMMAND, "no object");
		return(0);
	}
	if (nam == NULL) {
		error(COMMAND, "missing name");
		return(0);
	}
	if ((*nam == '*') | (*nam == '-')) {
		error(COMMAND, "illegal name");
		return(0);
	}
	if (getdobj(nam) != NULL) {
		error(COMMAND, "name already taken (clear first)");
		return(0);
	}
					/* allocate and copy struct */
	opdup = (DOBJECT *)malloc(sizeof(DOBJECT));
	if (opdup == NULL)
		error(SYSTEM, "out of memory in dobj_dup");
	*opdup = *op;
					/* rename */
	strcpy(opdup->name, nam);
					/* get our own copy of transform */
	for (opdup->xfac = 0; opdup->xfac < op->xfac; opdup->xfac++)
		opdup->xfav[opdup->xfac] = savestr(op->xfav[opdup->xfac]);
	opdup->xfav[opdup->xfac] = NULL;
					/* insert it into our list */
	opdup->next = dobjects;
	curobj = dobjects = opdup;
	return(1);
}


extern int
dobj_lighting(		/* set up lighting for display object */
	char	*nam,
	int	cn
)
{
	int	i, res[2];
	VIEW	*dv;
	register DOBJECT	*op;

	if (nam == NULL) {
		error(COMMAND, "no current object");
		return(0);
	}
	if (nam[0] == '*') {
		for (op = dobjects; op != NULL; op = op->next)
			if ((op->drawcode = cn) == DO_LIGHT)
				getdlights(op, 1);
			else
				op->ol = NULL;
	} else if ((op = getdobj(nam)) == NULL) {
		error(COMMAND, "unknown object");
		return(0);
	} else if ((op->drawcode = cn) == DO_LIGHT) {
		if (!getdlights(op, 1))
			error(COMMAND, "insufficient samples to light object");
	} else
		op->ol = NULL;

	if (dobj_lightsamp != NULL) {		/* restore beam set */
		dobj_lightsamp = NULL;
		beam_init(1);
		for (i = 0; (dv = dev_auxview(i, res)) != NULL; i++)
			beam_view(dv, res[0], res[1]);
		beam_sync(1);			/* update server */
	}
	return 0; /* XXX not sure if this is the right value */
}


extern double
dobj_trace(	/* check for ray intersection with object(s) */
	char	nm[],
	FVECT  rorg,
	FVECT  rdir
)
{
	register DOBJECT	*op;
	FVECT	xorg, xdir;
	double	darr[6];
					/* check each visible object? */
	if (nm == NULL || *nm == '*') {
		double	dist, mindist = 1.01*FHUGE;

		if (nm != NULL) nm[0] = '\0';
		for (op = dobjects; op != NULL; op = op->next) {
			if (op->drawcode == DO_HIDE)
				continue;
			dist = dobj_trace(op->name, rorg, rdir);
			if (dist < mindist) {
				if (nm != NULL) strcpy(nm, op->name);
				mindist = dist;
			}
		}
		return(mindist);
	}
					/* else check particular object */
	if ((op = getdobj(nm)) == NULL) {
		error(COMMAND, "unknown object");
		return(FHUGE);
	}
	if (op->xfac) {		/* put ray in local coordinates */
		multp3(xorg, rorg, op->xfb.b.xfm);
		multv3(xdir, rdir, op->xfb.b.xfm);
		VCOPY(darr, xorg); VCOPY(darr+3, xdir);
	} else {
		VCOPY(darr, rorg); VCOPY(darr+3, rdir);
	}
				/* trace it */
	if (process(&(op->rtp), (char *)darr, (char *)darr, sizeof(double),
			6*sizeof(double)) != sizeof(double))
		error(SYSTEM, "rtrace communication error");
				/* return distance */
	if (darr[0] >= .99*FHUGE)
		return(FHUGE);
	return(darr[0]*op->xfb.f.sca);
}


extern int
dobj_render(void)			/* render our objects in OpenGL */
{
	int	nrendered = 0;
	GLboolean	normalizing;
	GLfloat	vec[4];
	FVECT	v1;
	register DOBJECT	*op;
	register int	i;
					/* anything to render? */
	for (op = dobjects; op != NULL; op = op->next)
		if (op->drawcode != DO_HIDE)
			break;
	if (op == NULL)
		return(0);
					/* set up general rendering params */
	glGetBooleanv(GL_NORMALIZE, &normalizing);
	glPushAttrib(GL_LIGHTING_BIT|GL_TRANSFORM_BIT|GL_ENABLE_BIT|
		GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_POLYGON_BIT);
	glDepthFunc(GL_LESS);
	glEnable(GL_DEPTH_TEST);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glShadeModel(GL_SMOOTH);
	glFrontFace(GL_CCW);
	glDisable(GL_CULL_FACE);
	for (i = MAXLIGHTS; i--; )
		glDisable(glightid[i]);
	glEnable(GL_LIGHTING);
	rgl_checkerr("setting rendering modes in dobj_render");
					/* render each object */
	for (op = dobjects; op != NULL; op = op->next) {
		if (op->drawcode == DO_HIDE)
			continue;
					/* set up lighting */
		if (op->drawcode == DO_LIGHT && op->ol != NULL) {
			uby8	pval;
			double	expval, d;
						/* use computed sources */
			if (tmMapPixels(tmGlobal, &pval, &op->ol->larb,
					TM_NOCHROM, 1) != TM_E_OK)
				error(CONSISTENCY, "dobj_render w/o tone map");
			expval = pval * (WHTEFFICACY/256.) /
					tmLuminance(op->ol->larb);
			vec[0] = expval * colval(op->ol->lamb,RED);
			vec[1] = expval * colval(op->ol->lamb,GRN);
			vec[2] = expval * colval(op->ol->lamb,BLU);
			vec[3] = 1.;
			glLightModelfv(GL_LIGHT_MODEL_AMBIENT, vec);
			for (i = op->ol->nl; i--; ) {
				VCOPY(vec, op->ol->li[i].direc);
				vec[3] = 0.;
				glLightfv(glightid[i], GL_POSITION, vec);
				d = expval * op->ol->li[i].omega;
				vec[0] = d * colval(op->ol->li[i].val,RED);
				vec[1] = d * colval(op->ol->li[i].val,GRN);
				vec[2] = d * colval(op->ol->li[i].val,BLU);
				vec[3] = 1.;
				glLightfv(glightid[i], GL_SPECULAR, vec);
				glLightfv(glightid[i], GL_DIFFUSE, vec);
				vec[0] = vec[1] = vec[2] = 0.; vec[3] = 1.;
				glLightfv(glightid[i], GL_AMBIENT, vec);
				glEnable(glightid[i]);
			}
		} else {			/* fake lighting */
			vec[0] = vec[1] = vec[2] = 0.; vec[3] = 1.;
			glLightModelfv(GL_LIGHT_MODEL_AMBIENT, vec);
			getdcent(v1, op);
			VSUB(v1, odev.v.vp, v1);
			if (normalize(v1) <= getdrad(op)) {
				vec[0] = -odev.v.vdir[0];
				vec[1] = -odev.v.vdir[1];
				vec[2] = -odev.v.vdir[2];
			} else
				VCOPY(vec, v1);
			vec[3] = 0.;
			glLightfv(GL_LIGHT0, GL_POSITION, vec);
			vec[0] = vec[1] = vec[2] = .7; vec[3] = 1.;
			glLightfv(GL_LIGHT0, GL_SPECULAR, vec);
			glLightfv(GL_LIGHT0, GL_DIFFUSE, vec);
			vec[0] = vec[1] = vec[2] = .3; vec[3] = 1.;
			glLightfv(GL_LIGHT0, GL_AMBIENT, vec);
			glEnable(GL_LIGHT0);
		}
					/* set up object transform */
		if (op->xfac) {
			if (!normalizing && (op->xfb.f.sca < 1.-FTINY) |
						(op->xfb.f.sca > 1.+FTINY))
				glEnable(GL_NORMALIZE);
			glMatrixMode(GL_MODELVIEW);
			glPushMatrix();
					/* matrix order works out to same */
#ifdef SMLFLT
			glMultMatrixf((GLfloat *)op->xfb.f.xfm);
#else
			glMultMatrixd((GLdouble *)op->xfb.f.xfm);
#endif
		}
					/* render the display list */
		glCallList(op->listid);
		nrendered++;
					/* restore matrix */
		if (op->xfac) {
			glMatrixMode(GL_MODELVIEW);
			glPopMatrix();
			if (!normalizing)
				glDisable(GL_NORMALIZE);
		}
					/* restore lighting */
		if (op->drawcode == DO_LIGHT && op->ol != NULL)
			for (i = op->ol->nl; i--; )
				glDisable(glightid[i]);
		else
			glDisable(GL_LIGHT0);
					/* check errors */
	}
	glPopAttrib();			/* restore rendering params */
	rgl_checkerr("rendering objects in dobj_render");
	return(nrendered);
}
