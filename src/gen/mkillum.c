#ifndef lint
static const char RCSid[] = "$Id: mkillum.c,v 2.17 2003/09/24 14:55:54 greg Exp $";
#endif
/*
 * Make illum sources for optimizing rendering process
 */

#include  <signal.h>
#include  <ctype.h>

#include  "mkillum.h"
#include  "platform.h"

				/* default parameters */
#define  SAMPDENS	48		/* points per projected steradian */
#define  NSAMPS		32		/* samples per point */
#define  DFLMAT		"illum_mat"	/* material name */
#define  DFLDAT		"illum"		/* data file name */
				/* selection options */
#define  S_NONE		0		/* select none */
#define  S_ELEM		1		/* select specified element */
#define  S_COMPL	2		/* select all but element */
#define  S_ALL		3		/* select all */

				/* rtrace command and defaults */
char  *rtargv[64] = { "rtrace", "-dj", ".25", "-dr", "3", "-dv-",
		"-ab", "2", "-ad", "1024", "-as", "512", "-aa", ".1", };
int  rtargc = 14;
				/* overriding rtrace options */
char  *myrtopts[] = { "-I-", "-i-", "-ld-", "-ov", "-h-",
			"-fff", "-y", "0", NULL };

struct rtproc	rt;		/* our rtrace process */

struct illum_args  thisillum = {	/* our illum and default values */
		0,
		DFLMAT,
		DFLDAT,
		0,
		VOIDID,
		SAMPDENS,
		NSAMPS,
		0.,
	};

char	matcheck[MAXSTR];	/* current material to include or exclude */
int	matselect = S_ALL;	/* selection criterion */

FUN	ofun[NUMOTYPE] = INIT_OTYPE;	/* object types */

int	gargc;			/* global argc */
char	**gargv;		/* global argv */
#define  progname	gargv[0]

int	doneheader = 0;		/* printed header yet? */
#define  checkhead()	if (!doneheader++) printhead(gargc,gargv)

int	warnings = 1;		/* print warnings? */


main(argc, argv)		/* compute illum distributions using rtrace */
int	argc;
char	*argv[];
{
	char	*rtpath;
	FILE	*fp;
	register int	i;
				/* set global arguments */
	gargv = argv;
				/* set up rtrace command */
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '<' && argv[i][1] == '\0')
			break;
		rtargv[rtargc++] = argv[i];
		if (argv[i][0] == '-' && argv[i][1] == 'w')
			switch (argv[i][2]) {
			case '\0':
				warnings = !warnings;
				break;
			case '+':
			case 'T': case 't':
			case 'Y': case 'y':
			case '1':
				warnings = 1;
				break;
			case '-':
			case 'F': case 'f':
			case 'N': case 'n':
			case '0':
				warnings = 0;
				break;
			}
	}
	gargc = i;
	rtargc--;
	for (i = 0; myrtopts[i] != NULL; i++)
		rtargv[rtargc++] = myrtopts[i];
	rtargv[rtargc++] = argv[gargc-1];
	rtargv[rtargc] = NULL;
				/* just asking for defaults? */
	if (!strcmp(argv[gargc-1], "-defaults")) {
		printopts(); fflush(stdout);
		rtpath = getpath(rtargv[0], getenv("PATH"), X_OK);
		if (rtpath == NULL) {
			eputs(rtargv[0]);
			eputs(": command not found\n");
			exit(1);
		}
		execv(rtpath, rtargv);
		perror(rtpath);
		exit(1);
	}
	if (gargc < 2 || argv[gargc-1][0] == '-')
		error(USER, "missing octree argument");
				/* else initialize and run our calculation */
	init();
	if (gargc+1 < argc)
		for (i = gargc+1; i < argc; i++) {
			if ((fp = fopen(argv[i], "r")) == NULL) {
				sprintf(errmsg,
				"cannot open scene file \"%s\"", argv[i]);
				error(SYSTEM, errmsg);
			}
			filter(fp, argv[i]);
			fclose(fp);
		}
	else
		filter(stdin, "standard input");
	quit(0);
}


void
quit(status)			/* exit with status */
int  status;
{
	int	rtstat;

	rtstat = close_process(&(rt.pd));
	if (status == 0) {
		if (rtstat < 0)
			error(WARNING,
			"unknown return status from rtrace process");
		else
			status = rtstat;
	}
	exit(status);
}


init()				/* start rtrace and set up buffers */
{
	extern int  o_face(), o_sphere(), o_ring();
	int	maxbytes;
					/* set up object functions */
	ofun[OBJ_FACE].funp = o_face;
	ofun[OBJ_SPHERE].funp = o_sphere;
	ofun[OBJ_RING].funp = o_ring;
					/* set up signal handling */
#ifdef SIGPIPE /* not present on Windows */
	signal(SIGPIPE, quit);
#endif
					/* start rtrace process */
	errno = 0;
	maxbytes = open_process(&(rt.pd), rtargv);
	if (maxbytes == 0) {
		eputs(rtargv[0]);
		eputs(": command not found\n");
		exit(1);
	}
	if (maxbytes < 0)
		error(SYSTEM, "cannot start rtrace process");
	rt.bsiz = maxbytes/(6*sizeof(float));
	rt.buf = (float *)malloc(6*sizeof(float)*rt.bsiz--);
	rt.dest = (float **)malloc(sizeof(float *)*rt.bsiz);
	if (rt.buf == NULL || rt.dest == NULL)
		error(SYSTEM, "out of memory in init");
	rt.nrays = 0;
					/* set up urand */
	initurand(16384);
}


void
eputs(s)				/* put string to stderr */
register char  *s;
{
	static int  midline = 0;

	if (!*s) return;
	if (!midline) {
		fputs(progname, stderr);
		fputs(": ", stderr);
	}
	fputs(s, stderr);
	midline = s[strlen(s)-1] != '\n';
}


void
wputs(s)			/* print warning if enabled */
char  *s;
{
	if (warnings)
		eputs(s);
}


filter(infp, name)		/* process stream */
register FILE	*infp;
char	*name;
{
	char	buf[512];
	FILE	*pfp;
	register int	c;

	while ((c = getc(infp)) != EOF) {
		if (isspace(c))
			continue;
		if (c == '#') {				/* comment/options */
			buf[0] = c;
			fgets(buf+1, sizeof(buf)-1, infp);
			xoptions(buf, name);
		} else if (c == '!') {			/* command */
			buf[0] = c;
			fgetline(buf+1, sizeof(buf)-1, infp);
			if ((pfp = popen(buf+1, "r")) == NULL) {
				sprintf(errmsg, "cannot execute \"%s\"", buf);
				error(SYSTEM, errmsg);
			}
			filter(pfp, buf);
			pclose(pfp);
		} else {				/* object */
			ungetc(c, infp);
			xobject(infp, name);
		}
	}
}


xoptions(s, nm)			/* process options in string s */
char	*s;
char	*nm;
{
	extern FILE	*freopen();
	char	buf[64];
	int	nerrs = 0;
	register char	*cp;

	if (strncmp(s, "#@mkillum", 9) || !isspace(s[9])) {
		fputs(s, stdout);		/* not for us */
		return;
	}
	cp = s+10;
	while (*cp) {
		switch (*cp) {
		case ' ':
		case '\t':
		case '\n':
		case '\r':
		case '\f':
			cp++;
			continue;
		case 'm':			/* material name */
			if (*++cp != '=')
				break;
			if (!*++cp || isspace(*cp))
				break;
			atos(thisillum.matname, MAXSTR, cp);
			cp = sskip(cp);
			if (!(thisillum.flags & IL_DATCLB)) {
				strcpy(thisillum.datafile, thisillum.matname);
				thisillum.dfnum = 0;
			}
			continue;
		case 'f':			/* data file name */
			if (*++cp != '=')
				break;
			if (!*++cp || isspace(*cp)) {
				strcpy(thisillum.datafile,thisillum.matname);
				thisillum.dfnum = 0;
				thisillum.flags &= ~IL_DATCLB;
				continue;
			}
			atos(thisillum.datafile, MAXSTR, cp);
			cp = sskip(cp);
			thisillum.dfnum = 0;
			thisillum.flags |= IL_DATCLB;
			continue;
		case 'i':			/* include material */
		case 'e':			/* exclude material */
			if (cp[1] != '=')
				break;
			matselect = (*cp == 'i') ? S_ELEM : S_COMPL;
			cp += 2;
			atos(matcheck, MAXSTR, cp);
			cp = sskip(cp);
			continue;
		case 'a':			/* use everything */
			cp = sskip(cp);
			matselect = S_ALL;
			continue;
		case 'n':			/* use nothing (passive) */
			cp = sskip(cp);
			matselect = S_NONE;
			continue;
		case 'c':			/* color calculation */
			if (*++cp != '=')
				break;
			switch (*++cp) {
			case 'a':			/* average */
				thisillum.flags = (thisillum.flags|IL_COLAVG)
							& ~IL_COLDST;
				break;
			case 'd':			/* distribution */
				thisillum.flags |= (IL_COLDST|IL_COLAVG);
				break;
			case 'n':			/* none */
				thisillum.flags &= ~(IL_COLAVG|IL_COLDST);
				break;
			default:
				goto opterr;
			}
			cp = sskip(cp);
			continue;
		case 'd':			/* point sample density */
			if (*++cp != '=')
				break;
			if (!isintd(++cp, " \t\n\r"))
				break;
			thisillum.sampdens = atoi(cp);
			cp = sskip(cp);
			continue;
		case 's':			/* point super-samples */
			if (*++cp != '=')
				break;
			if (!isintd(++cp, " \t\n\r"))
				break;
			thisillum.nsamps = atoi(cp);
			cp = sskip(cp);
			continue;
		case 'l':			/* light sources */
			cp++;
			if (*cp == '+')
				thisillum.flags |= IL_LIGHT;
			else if (*cp == '-')
				thisillum.flags &= ~IL_LIGHT;
			else
				break;
			cp++;
			continue;
		case 'b':			/* brightness */
			if (*++cp != '=')
				break;
			if (!isfltd(++cp, " \t\n\r"))
				break;
			thisillum.minbrt = atof(cp);
			if (thisillum.minbrt < 0.)
				thisillum.minbrt = 0.;
			cp = sskip(cp);
			continue;
		case 'o':			/* output file */
			if (*++cp != '=')
				break;
			if (!*++cp || isspace(*cp))
				break;
			atos(buf, sizeof(buf), cp);
			cp = sskip(cp);
			if (freopen(buf, "w", stdout) == NULL) {
				sprintf(errmsg,
				"cannot open output file \"%s\"", buf);
				error(SYSTEM, errmsg);
			}
			doneheader = 0;
			continue;
		case '!':			/* processed file! */
			sprintf(errmsg, "(%s): already processed!", nm);
			error(WARNING, errmsg);
			matselect = S_NONE;
			return;
		}
	opterr:					/* skip faulty option */
		while (*cp && !isspace(*cp))
			cp++;
		nerrs++;
	}
						/* print header? */
	checkhead();
						/* issue warnings? */
	if (nerrs) {
		sprintf(errmsg, "(%s): %d error(s) in option line:",
				nm, nerrs);
		error(WARNING, errmsg);
		wputs(s);
		printf("# %s: the following option line has %d error(s):\n",
				progname, nerrs);
	}
						/* print pure comment */
	printf("# %s", s+2);
}


printopts()			/* print out option default values */
{
	printf("m=%-15s\t\t# material name\n", thisillum.matname);
	printf("f=%-15s\t\t# data file name\n", thisillum.datafile);
	if (thisillum.flags & IL_COLAVG)
		if (thisillum.flags & IL_COLDST)
			printf("c=d\t\t\t\t# color distribution\n");
		else
			printf("c=a\t\t\t\t# color average\n");
	else
		printf("c=n\t\t\t\t# color none\n");
	if (thisillum.flags & IL_LIGHT)
		printf("l+\t\t\t\t# light type on\n");
	else
		printf("l-\t\t\t\t# light type off\n");
	printf("d=%d\t\t\t\t# density of points\n", thisillum.sampdens);
	printf("s=%d\t\t\t\t# samples per point\n", thisillum.nsamps);
	printf("b=%f\t\t\t# minimum average brightness\n", thisillum.minbrt);
}


printhead(ac, av)			/* print out header */
register int  ac;
register char  **av;
{
	putchar('#');
	while (ac-- > 0) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	fputs("\n#@mkillum !\n", stdout);
}


xobject(fp, nm)				/* translate an object from fp */
FILE  *fp;
char  *nm;
{
	OBJREC  thisobj;
	char  str[MAXSTR];
	int  doit;
					/* read the object */
	if (fgetword(thisillum.altmat, MAXSTR, fp) == NULL)
		goto readerr;
	if (fgetword(str, MAXSTR, fp) == NULL)
		goto readerr;
					/* is it an alias? */
	if (!strcmp(str, ALIASKEY)) {
		if (fgetword(str, MAXSTR, fp) == NULL)
			goto readerr;
		printf("\n%s %s %s", thisillum.altmat, ALIASKEY, str);
		if (fgetword(str, MAXSTR, fp) == NULL)
			goto readerr;
		printf("\t%s\n", str);
		return;
	}
	thisobj.omod = OVOID;		/* unused field */
	if ((thisobj.otype = otype(str)) < 0) {
		sprintf(errmsg, "(%s): unknown type \"%s\"", nm, str);
		error(USER, errmsg);
	}
	if (fgetword(str, MAXSTR, fp) == NULL)
		goto readerr;
	thisobj.oname = str;
	if (readfargs(&thisobj.oargs, fp) != 1)
		goto readerr;
	thisobj.os = NULL;
					/* check for translation */
	switch (matselect) {
	case S_NONE:
		doit = 0;
		break;
	case S_ALL:
		doit = 1;
		break;
	case S_ELEM:
		doit = !strcmp(thisillum.altmat, matcheck);
		break;
	case S_COMPL:
		doit = strcmp(thisillum.altmat, matcheck);
		break;
	}
	doit = doit && issurface(thisobj.otype);
						/* print header? */
	checkhead();
						/* process object */
	if (doit)
		(*ofun[thisobj.otype].funp)(&thisobj, &thisillum, &rt, nm);
	else
		printobj(thisillum.altmat, &thisobj);
						/* free arguments */
	freefargs(&thisobj.oargs);
	return;
readerr:
	sprintf(errmsg, "(%s): error reading scene", nm);
	error(USER, errmsg);
}
