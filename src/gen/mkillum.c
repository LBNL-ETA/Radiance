#ifndef lint
static const char RCSid[] = "$Id: mkillum.c,v 2.37 2011/08/15 19:48:06 greg Exp $";
#endif
/*
 * Make illum sources for optimizing rendering process
 */

#include  <signal.h>
#include  <ctype.h>

#include  "mkillum.h"

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

struct illum_args  thisillum = {	/* our illum and default values */
		0,
		UDzpos,
		0.,
		DFLMAT,
		DFLDAT,
		0,
		VOIDID,
		SAMPDENS,
		NSAMPS,
		NULL,
		0.,
	};

char	matcheck[MAXSTR];	/* current material to include or exclude */
int	matselect = S_ALL;	/* selection criterion */

int	gargc;			/* global argc */
char	**gargv;		/* global argv */

int	doneheader = 0;		/* printed header yet? */
#define  checkhead()	if (!doneheader++) printhead(gargc,gargv)

int	warnings = 1;		/* print warnings? */

void init(char *octnm, int np);
void filter(register FILE	*infp, char	*name);
void xoptions(char	*s, char	*nm);
void printopts(void);
void printhead(register int  ac, register char  **av);
void xobject(FILE  *fp, char  *nm);


int
main(		/* compute illum distributions using rtrace */
	int	argc,
	char	*argv[]
)
{
	int	nprocs = 1;
	FILE	*fp;
	int	rval;
	register int	i;
				/* set global arguments */
	gargv = argv;
	progname = gargv[0];
				/* set up rendering defaults */
	dstrsrc = 0.5;
	directrelay = 3;
	ambounce = 2;
				/* get options from command line */
	for (i = 1; i < argc; i++) {
		while ((rval = expandarg(&argc, &argv, i)) > 0)
			;
		if (rval < 0) {
			sprintf(errmsg, "cannot expand '%s'", argv[i]);
			error(SYSTEM, errmsg);
		}
		if (argv[i][0] != '-')
			break;
		if (!strcmp(argv[i], "-w")) {
			warnings = 0;
			continue;
		}
		if (!strcmp(argv[i], "-n")) {
			nprocs = atoi(argv[++i]);
			if (nprocs <= 0)
				error(USER, "illegal number of processes");
			continue;
		}
		if (!strcmp(argv[i], "-defaults")) {
			printopts();
			print_rdefaults();
			quit(0);
		}
		rval = getrenderopt(argc-i, argv+i);
		if (rval < 0) {
			sprintf(errmsg, "bad render option at '%s'", argv[i]);
			error(USER, errmsg);
		}
		i += rval;
	}
	gargc = ++i;
				/* add "mandatory" render options */
	do_irrad = 0;
	if (gargc > argc || argv[gargc-1][0] == '-')
		error(USER, "missing octree argument");
				/* else initialize and run our calculation */
	init(argv[gargc-1], nprocs);
	if (gargc < argc) {
		if (gargc == argc-1 || argv[gargc][0] != '<' || argv[gargc][1])
			error(USER, "use '< file1 file2 ..' for multiple inputs");
		for (i = gargc+1; i < argc; i++) {
			if ((fp = fopen(argv[i], "r")) == NULL) {
				sprintf(errmsg,
				"cannot open scene file \"%s\"", argv[i]);
				error(SYSTEM, errmsg);
			}
			filter(fp, argv[i]);
			fclose(fp);
		}
	} else
		filter(stdin, "standard input");
	quit(0);
	return 0; /* pro forma return */
}


void
init(char *octnm, int np)		/* start rendering process(es) */
{
					/* set up signal handling */
	signal(SIGINT, quit);
#ifdef SIGHUP
	signal(SIGHUP, quit);
#endif
#ifdef SIGTERM
	signal(SIGTERM, quit);
#endif
#ifdef SIGPIPE
	signal(SIGPIPE, quit);
#endif
					/* start rendering process(es) */
	ray_pinit(octnm, np);
}


void
eputs(				/* put string to stderr */
	register char  *s
)
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


void
quit(ec)			/* make sure exit is called */
int	ec;
{
	if (ray_pnprocs > 0)	/* close children if any */
		ray_pclose(0);		
	exit(ec);
}


void
filter(		/* process stream */
	register FILE	*infp,
	char	*name
)
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


void
xoptions(			/* process options in string s */
	char	*s,
	char	*nm
)
{
	extern FILE	*freopen();
	char	buf[64];
	int	negax;
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
		case 'd':			/* sample density / BSDF data */
			if (*++cp != '=')
				break;
			if (thisillum.sd != NULL) {
				free_BSDF(thisillum.sd);
				thisillum.sd = NULL;
			}
			if (!*++cp || isspace(*cp))
				continue;
			if (isintd(cp, " \t\n\r")) {
				thisillum.sampdens = atoi(cp);
			} else {
				atos(buf, sizeof(buf), cp);
				thisillum.sd = load_BSDF(buf);
				if (thisillum.sd == NULL)
					break;
			}
			cp = sskip(cp);
			continue;
		case 's':			/* surface super-samples */
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
		case 'u':			/* up direction */
			if (*++cp != '=')
				break;
			if (!*++cp || isspace(*cp)) {
				thisillum.udir = UDunknown;
				continue;
			}
			negax = 0;
			if (*cp == '+')
				cp++;
			else if (*cp == '-') {
				negax++;
				cp++;
			}
			switch (*cp++) {
			case 'x':
			case 'X':
				thisillum.udir = negax ? UDxneg : UDxpos;
				break;
			case 'y':
			case 'Y':
				thisillum.udir = negax ? UDyneg : UDypos;
				break;
			case 'z':
			case 'Z':
				thisillum.udir = negax ? UDzneg : UDzpos;
				break;
			default:
				thisillum.udir = UDunknown;
				break;
			}
			if (thisillum.udir == UDunknown || !isspace(*cp))
				break;
			continue;
		case 't':			/* object thickness */
			if (*++cp != '=')
				break;
			if (!isfltd(++cp, " \t\n\r"))
				break;
			thisillum.thick = atof(cp);
			if (thisillum.thick < .0)
				thisillum.thick = .0;
			cp = sskip(cp);
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

void
printopts(void)			/* print out option default values */
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
	printf("d=%d\t\t\t\t# density of directions\n", thisillum.sampdens);
	printf("s=%d\t\t\t\t# samples per direction\n", thisillum.nsamps);
	printf("b=%f\t\t\t# minimum average brightness\n", thisillum.minbrt);
	switch (thisillum.udir) {
	case UDzneg:
		fputs("u=-Z\t\t\t\t# up is negative Z\n", stdout);
		break;
	case UDyneg:
		fputs("u=-Y\t\t\t\t# up is negative Y\n", stdout);
		break;
	case UDxneg:
		fputs("u=-X\t\t\t\t# up is negative X\n", stdout);
		break;
	case UDxpos:
		fputs("u=+X\t\t\t\t# up is positive X\n", stdout);
		break;
	case UDypos:
		fputs("u=+Y\t\t\t\t# up is positive Y\n", stdout);
		break;
	case UDzpos:
		fputs("u=+Z\t\t\t\t# up is positive Z\n", stdout);
		break;
	case UDunknown:
		break;
	}
	printf("t=%f\t\t\t# object thickness\n", thisillum.thick);
}


void
printhead(			/* print out header */
	register int  ac,
	register char  **av
)
{
	putchar('#');
	while (ac-- > 0) {
		putchar(' ');
		fputs(*av++, stdout);
	}
	fputs("\n#@mkillum !\n", stdout);
}


void
xobject(				/* translate an object from fp */
	FILE  *fp,
	char  *nm
)
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
		switch (thisobj.otype) {
		case OBJ_FACE:
			my_face(&thisobj, &thisillum, nm);
			break;
		case OBJ_SPHERE:
			my_sphere(&thisobj, &thisillum, nm);
			break;
		case OBJ_RING:
			my_ring(&thisobj, &thisillum, nm);
			break;
		default:
			my_default(&thisobj, &thisillum, nm);
			break;
		}
	else
		printobj(thisillum.altmat, &thisobj);
						/* free arguments */
	freefargs(&thisobj.oargs);
	return;
readerr:
	sprintf(errmsg, "(%s): error reading input", nm);
	error(USER, errmsg);
}
