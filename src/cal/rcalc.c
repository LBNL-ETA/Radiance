#ifndef lint
static const char RCSid[] = "$Id: rcalc.c,v 1.8 2003/08/01 14:14:23 schorsch Exp $";
#endif
/*
 *  rcalc.c - record calculator program.
 *
 *     9/11/87
 */

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <math.h>
#include  <ctype.h>

#include  "platform.h"
#include  "calcomp.h"
#include  "rterror.h"

#define  isnum(c)       (isdigit(c) || (c)=='-' || (c)=='.' \
				|| (c)=='+' || (c)=='e' || (c)=='E')

#define  isblnk(c)      (igneol ? isspace(c) : (c)==' '||(c)=='\t')

#define  INBSIZ         4096    /* longest record */
#define  MAXCOL         32      /* number of columns recorded */

				/* field type specifications */
#define  F_NUL          0               /* empty */
#define  F_TYP          0x7000          /* mask for type */
#define  F_WID          0x0fff          /* mask for width */
#define  T_LIT          0x1000          /* string literal */
#define  T_STR          0x2000          /* string variable */
#define  T_NUM          0x3000          /* numeric value */

struct strvar {                 /* string variable */
	char  *name;
	char  *val;
	char  *preset;
	struct strvar  *next;
};

struct field {                  /* record format structure */
	int  type;                      /* type of field (& width) */
	union {
		char  *sl;                      /* string literal */
		struct strvar  *sv;             /* string variable */
		char  *nv;                      /* numeric variable */
		EPNODE  *ne;                    /* numeric expression */
	} f;                            /* field contents */
	struct field  *next;            /* next field in record */
};

#define  savqstr(s)     strcpy(emalloc(strlen(s)+1),s)
#define  freqstr(s)     efree(s)

static void scaninp(void), advinp(void), resetinp(void);
static void putrec(void), putout(void), nbsynch(void);
static int getrec(void);
static void execute(char *file);
static void initinp(FILE *fp);
static void svpreset(char *eqn);
static void readfmt(char *spec, int output);
static int readfield(char **pp);
static int getfield(struct field *f);
static void chanset(int n, double v);
static void bchanset(int n, double v);
static struct strvar* getsvar(char *svname);

struct field  *inpfmt = NULL;   /* input record format */
struct field  *outfmt = NULL;   /* output record structure */
struct strvar  *svhead = NULL;  /* string variables */

int  blnkeq = 1;                /* blanks compare equal? */
int  igneol = 0;                /* ignore end of line? */
char  sepchar = '\t';           /* input/output separator */
int  noinput = 0;               /* no input records? */
int  nbicols = 0;		/* number of binary input columns */
int  bocols = 0;		/* produce binary output columns */
char  inpbuf[INBSIZ];           /* input buffer */
double  colval[MAXCOL];         /* input column values */
unsigned long  colflg = 0;      /* column retrieved flags */
int  colpos;                    /* output column position */

int  nowarn = 0;                /* non-fatal diagnostic output */
int  unbuff = 0;		/* unbuffered output (flush each record) */

struct {
	FILE  *fin;                     /* input file */
	int  chr;                       /* next character */
	char  *beg;                     /* home position */
	char  *pos;                     /* scan position */
	char  *end;                     /* read position */
} ipb;                          /* circular lookahead buffer */


int
main(
int  argc,
char  *argv[]
)
{
	int  i;

	esupport |= E_VARIABLE|E_FUNCTION|E_INCHAN|E_OUTCHAN|E_RCONST;
	esupport &= ~(E_REDEFW);

#ifdef  BIGGERLIB
	biggerlib();
#endif
	varset("PI", ':', 3.14159265358979323846);

	for (i = 1; i < argc && argv[i][0] == '-'; i++)
		switch (argv[i][1]) {
		case 'b':
			blnkeq = !blnkeq;
			break;
		case 'l':
			igneol = !igneol;
			break;
		case 't':
			sepchar = argv[i][2];
			break;
		case 's':
			svpreset(argv[++i]);
			break;
		case 'f':
			fcompile(argv[++i]);
			break;
		case 'e':
			scompile(argv[++i], NULL, 0);
			break;
		case 'n':
			noinput = 1;
			break;
		case 'i':
			switch (argv[i][2]) {
			case '\0':
				nbicols = 0;
				readfmt(argv[++i], 0);
				break;
			case 'a':
				nbicols = 0;
				break;
			case 'd':
				if (isdigit(argv[i][3]))
					nbicols = atoi(argv[i]+3);
				else
					nbicols = 1;
				break;
			case 'f':
				if (isdigit(argv[i][3]))
					nbicols = -atoi(argv[i]+3);
				else
					nbicols = -1;
				break;
			default:
				goto userr;
			}
			break;
		case 'o':
			switch (argv[i][2]) {
			case '\0':
				bocols = 0;
				readfmt(argv[++i], 1);
				break;
			case 'a':
				bocols = 0;
				break;
			case 'd':
				bocols = 1;
				break;
			case 'f':
				bocols = -1;
				break;
			}
			break;
		case 'w':
			nowarn = !nowarn;
			break;
		case 'u':
			unbuff = !unbuff;
			break;
		default:;
		userr:
			eputs("Usage: ");
			eputs(argv[0]);
eputs(" [-b][-l][-n][-w][-u][-tS][-s svar=sval][-e expr][-f source][-i infmt][-o outfmt] [file]\n");
			quit(1);
		}

	if (noinput) {          /* produce a single output record */
		eclock++;
		putout();
		quit(0);
	}

	if (blnkeq)             /* for efficiency */
		nbsynch();

	if (i == argc)          /* from stdin */
		execute(NULL);
	else                    /* from one or more files */
		for ( ; i < argc; i++)
			execute(argv[i]);
	
	quit(0);
}


static void
nbsynch(void)               /* non-blank starting synch character */
{
	if (inpfmt == NULL || (inpfmt->type & F_TYP) != T_LIT)
		return;
	while (isblnk(*inpfmt->f.sl))
		inpfmt->f.sl++;
	if (!*inpfmt->f.sl)
		inpfmt = inpfmt->next;
}


int
getinputrec(		/* get next input record */
FILE  *fp
)
{
	if (inpfmt != NULL)
		return(getrec());
	if (nbicols > 0)
		return(fread(inpbuf, sizeof(double),
					nbicols, fp) == nbicols);
	if (nbicols < 0)
		return(fread(inpbuf, sizeof(float),
					-nbicols, fp) == -nbicols);
	return(fgets(inpbuf, INBSIZ, fp) != NULL);
}


static void
execute(           /* process a file */
char  *file
)
{
	int  conditional = vardefined("cond");
	long  nrecs = 0;
	long  nout = 0;
	FILE  *fp;
	
	if (file == NULL)
		fp = stdin;
	else if ((fp = fopen(file, "r")) == NULL) {
		eputs(file);
		eputs(": cannot open\n");
		quit(1);
	}
	if (inpfmt != NULL)
		initinp(fp);
	
	while (getinputrec(fp)) {
		varset("recno", '=', (double)++nrecs);
		colflg = 0;
		eclock++;
		if (!conditional || varvalue("cond") > 0.0) {
			varset("outno", '=', (double)++nout);
			putout();
		}
	}
	fclose(fp);
}


static void
putout(void)                /* produce an output record */
{

	colpos = 0;
	if (outfmt != NULL)
		putrec();
	else if (bocols)
		chanout(bchanset);
	else
		chanout(chanset);
	if (colpos && !bocols)
		putchar('\n');
	if (unbuff)
		fflush(stdout);
}


double
chanvalue(            /* return value for column n */
int  n
)
{
	int  i;
	register char  *cp;

	if (noinput || inpfmt != NULL) {
		eputs("no column input\n");
		quit(1);
	}
	if (n < 1) {
		eputs("illegal channel number\n");
		quit(1);
	}
	if (nbicols > 0) {
		if (n > nbicols)
			return(0.0);
		cp = inpbuf + (n-1)*sizeof(double);
		return(*(double *)cp);
	}
	if (nbicols < 0) {
		if (n > -nbicols)
			return(0.0);
		cp = inpbuf + (n-1)*sizeof(float);
		return(*(float *)cp);
	}
	if (n <= MAXCOL && colflg & 1L<<(n-1))
		return(colval[n-1]);

	cp = inpbuf;
	for (i = 1; i < n; i++)
		if (blnkeq && isspace(sepchar)) {
			while (isspace(*cp))
				cp++;
			while (*cp && !isspace(*cp))
				cp++;
		} else
			while (*cp && *cp++ != sepchar)
				;

	while (isspace(*cp))            /* some atof()'s don't like tabs */
		cp++;

	if (n <= MAXCOL) {
		colflg |= 1L<<(n-1);
		return(colval[n-1] = atof(cp));
	} else
		return(atof(cp));
}


void
chanset(                   /* output column n */
int  n,
double  v
)
{
	if (colpos == 0)                /* no leading separator */
		colpos = 1;
	while (colpos < n) {
		putchar(sepchar);
		colpos++;
	}
	printf("%.9g", v);
}


void
bchanset(                   /* output binary channel n */
int  n,
double  v
)
{
	static char	zerobuf[sizeof(double)];

	while (++colpos < n)
		fwrite(zerobuf,
			bocols>0 ? sizeof(double) : sizeof(float),
			1, stdout);
	if (bocols > 0)
		fwrite(&v, sizeof(double), 1, stdout);
	else {
		float	fval = v;
		fwrite(&fval, sizeof(float), 1, stdout);
	}
}


static void
readfmt(                   /* read record format */
char  *spec,
int  output
)
{
	int  fd;
	char  *inptr;
	struct field  fmt;
	int  res;
	register struct field  *f;
						/* check for inline format */
	for (inptr = spec; *inptr; inptr++)
		if (*inptr == '$')
			break;
	if (*inptr)                             /* inline */
		inptr = spec;
	else {                                  /* from file */
		if ((fd = open(spec, 0)) == -1) {
			eputs(spec);
			eputs(": cannot open\n");
			quit(1);
		}
		res = read(fd, inpbuf+1, INBSIZ-1);
		if (res <= 0 || res >= INBSIZ-1) {
			eputs(spec);
			if (res < 0)
				eputs(": read error\n");
			else if (res == 0)
				eputs(": empty file\n");
			else if (res >= INBSIZ-1)
				eputs(": format too long\n");
			quit(1);
		}
		close(fd);
		(inptr=inpbuf+1)[res] = '\0';
	}
	f = &fmt;                               /* get fields */
	while ((res = readfield(&inptr)) != F_NUL) {
		f->next = (struct field *)emalloc(sizeof(struct field));
		f = f->next;
		f->type = res;
		switch (res & F_TYP) {
		case T_LIT:
			f->f.sl = savqstr(inpbuf);
			break;
		case T_STR:
			f->f.sv = getsvar(inpbuf);
			break;
		case T_NUM:
			if (output)
				f->f.ne = eparse(inpbuf);
			else
				f->f.nv = savestr(inpbuf);
			break;
		}
					/* add final newline if necessary */
		if (!igneol && *inptr == '\0' && inptr[-1] != '\n')
			inptr = "\n";
	}
	f->next = NULL;
	if (output)
		outfmt = fmt.next;
	else
		inpfmt = fmt.next;
}


static int
readfield(                   /* get next field in format */
register char  **pp
)
{
	int  type = F_NUL;
	int  width = 0;
	register char  *cp;
	
	cp = inpbuf;
	while (cp < &inpbuf[INBSIZ-1] && **pp != '\0') {
		width++;
		switch (type) {
		case F_NUL:
			if (**pp == '$') {
				(*pp)++;
				width++;
				if (**pp == '{') {
					type = T_NUM;
					(*pp)++;
					continue;
				} else if (**pp == '(') {
					type = T_STR;
					(*pp)++;
					continue;
				} else if (**pp != '$') {
					eputs("format error\n");
					quit(1);
				}
				width--;
			}
			type = T_LIT;
			*cp++ = *(*pp)++;
			continue;
		case T_LIT:
			if (**pp == '$') {
				width--;
				break;
			}
			*cp++ = *(*pp)++;
			continue;
		case T_NUM:
			if (**pp == '}') {
				(*pp)++;
				break;
			}
			if (!isspace(**pp))
				*cp++ = **pp;
			(*pp)++;
			continue;
		case T_STR:
			if (**pp == ')') {
				(*pp)++;
				break;
			}
			if (!isspace(**pp))
				*cp++ = **pp;
			(*pp)++;
			continue;
		}
		break;
	}
	*cp = '\0';
	return(type | width);
}


struct strvar *
getsvar(                         /* get string variable */
char  *svname
)
{
	register struct strvar  *sv;
	
	for (sv = svhead; sv != NULL; sv = sv->next)
		if (!strcmp(sv->name, svname))
			return(sv);
	sv = (struct strvar *)emalloc(sizeof(struct strvar));
	sv->name = savqstr(svname);
	sv->val = sv->preset = NULL;
	sv->next = svhead;
	svhead = sv;
	return(sv);
}


static void
svpreset(                    /* preset a string variable */
char  *eqn
)
{
	register struct strvar  *sv;
	register char  *val;

	for (val = eqn; *val != '='; val++)
		if (!*val)
			return;
	*val++ = '\0';
	sv = getsvar(eqn);
	if (sv->preset != NULL)
		freqstr(sv->preset);
	if (sv->val != NULL)
		freqstr(sv->val);
	sv->val = sv->preset = savqstr(val);
	*--val = '=';
}


static void
clearrec(void)			/* clear input record variables */
{
	register struct field  *f;

	for (f = inpfmt; f != NULL; f = f->next)
		switch (f->type & F_TYP) {
		case T_NUM:
			dremove(f->f.nv);
			break;
		case T_STR:
			if (f->f.sv->val != f->f.sv->preset) {
				freqstr(f->f.sv->val);
				f->f.sv->val = f->f.sv->preset;
			}
			break;
		}
}


static int
getrec(void)                                /* get next record from file */
{
	int  eatline;
	register struct field  *f;
	
	while (ipb.chr != EOF) {
		eatline = !igneol && ipb.chr != '\n';
		if (blnkeq)             /* beware of nbsynch() */
			while (isblnk(ipb.chr))
				scaninp();
		clearrec();		/* start with fresh record */
		for (f = inpfmt; f != NULL; f = f->next)
			if (getfield(f) == -1)
				break;
		if (f == NULL) {
			advinp();
			return(1);
		}
		resetinp();
		if (eatline) {          /* eat rest of line */
			while (ipb.chr != '\n') {
				if (ipb.chr == EOF)
					return(0);
				scaninp();
			}
			scaninp();
			advinp();
		}
	}
	return(0);
}


static int
getfield(                             /* get next field */
register struct field  *f
)
{
	static char  buf[RMAXWORD+1];            /* no recursion! */
	int  delim, inword;
	double  d;
	char  *np;
	register char  *cp;

	switch (f->type & F_TYP) {
	case T_LIT:
		cp = f->f.sl;
		do {
			if (blnkeq && isblnk(*cp)) {
				if (!isblnk(ipb.chr))
					return(-1);
				do
					cp++;
				while (isblnk(*cp));
				do
					scaninp();
				while (isblnk(ipb.chr));
			} else if (*cp == ipb.chr) {
				cp++;
				scaninp();
			} else
				return(-1);
		} while (*cp);
		return(0);
	case T_STR:
		if (f->next == NULL || (f->next->type & F_TYP) != T_LIT)
			delim = EOF;
		else
			delim = f->next->f.sl[0];
		cp = buf;
		do {
			if (ipb.chr == EOF)
				inword = 0;
			else if (blnkeq && delim != EOF)
				inword = isblnk(delim) ?
						!isblnk(ipb.chr)
						: ipb.chr != delim;
			else
				inword = cp-buf < (f->type & F_WID);
			if (inword) {
				*cp++ = ipb.chr;
				scaninp();
			}
		} while (inword && cp < &buf[RMAXWORD]);
		*cp = '\0';
		if (f->f.sv->val == NULL)
			f->f.sv->val = savqstr(buf);	/* first setting */
		else if (strcmp(f->f.sv->val, buf))
			return(-1);			/* doesn't match! */
		return(0);
	case T_NUM:
		if (f->next == NULL || (f->next->type & F_TYP) != T_LIT)
			delim = EOF;
		else
			delim = f->next->f.sl[0];
		np = NULL;
		cp = buf;
		do {
			if (!((np==NULL&&isblnk(ipb.chr)) || isnum(ipb.chr)))
				inword = 0;
			else if (blnkeq && delim != EOF)
				inword = isblnk(delim) ?
						!isblnk(ipb.chr)
						: ipb.chr != delim;
			else
				inword = cp-buf < (f->type & F_WID);
			if (inword) {
				if (np==NULL && !isblnk(ipb.chr))
					np = cp;
				*cp++ = ipb.chr;
				scaninp();
			}
		} while (inword && cp < &buf[RMAXWORD]);
		*cp = '\0';
		d = np==NULL ? 0. : atof(np);
		if (!vardefined(f->f.nv))
			varset(f->f.nv, '=', d);	/* first setting */
		else if ((d = (varvalue(f->f.nv)-d)/(d==0.?1.:d)) > .001
				|| d < -.001)
			return(-1);			/* doesn't match! */
		return(0);
	}
	return -1; /* pro forma return */
}


static void
putrec(void)                                /* output a record */
{
	char  fmt[32];
	register int  n;
	register struct field  *f;
	int  adlast, adnext;
	
	adlast = 0;
	for (f = outfmt; f != NULL; f = f->next) {
		adnext =        blnkeq &&
				f->next != NULL &&
				!( (f->next->type&F_TYP) == T_LIT &&
					f->next->f.sl[0] == ' ' );
		switch (f->type & F_TYP) {
		case T_LIT:
			fputs(f->f.sl, stdout);
			adlast = f->f.sl[(f->type&F_WID)-1] != ' ';
			break;
		case T_STR:
			if (f->f.sv->val == NULL) {
				eputs(f->f.sv->name);
				eputs(": undefined string\n");
				quit(1);
			}
			n = (int)(f->type & F_WID) - strlen(f->f.sv->val);
			if (adlast)
				fputs(f->f.sv->val, stdout);
			if (!(adlast && adnext))
				while (n-- > 0)
					putchar(' ');
			if (!adlast)
				fputs(f->f.sv->val, stdout);
			adlast = 1;
			break;
		case T_NUM:
			n = f->type & F_WID;
			if (adlast && adnext)
				strcpy(fmt, "%g");
			else if (adlast)
				sprintf(fmt, "%%-%dg", n);
			else
				sprintf(fmt, "%%%dg", n);
			printf(fmt, evalue(f->f.ne));
			adlast = 1;
			break;
		}
	}
}


static void
initinp(FILE  *fp)                     /* prepare lookahead buffer */

{
	ipb.fin = fp;
	ipb.beg = ipb.end = inpbuf;
	ipb.pos = inpbuf-1;             /* position before beginning */
	ipb.chr = '\0';
	scaninp();
}


static void
scaninp(void)                       /* scan next character */
{
	if (ipb.chr == EOF)
		return;
	if (++ipb.pos >= &inpbuf[INBSIZ])
		ipb.pos = inpbuf;
	if (ipb.pos == ipb.end) {               /* new character */
		if ((ipb.chr = getc(ipb.fin)) != EOF) {
			*ipb.end = ipb.chr;
			if (++ipb.end >= &inpbuf[INBSIZ])
				ipb.end = inpbuf;
			if (ipb.end == ipb.beg)
				ipb.beg = NULL;
		}
	} else
		ipb.chr = *ipb.pos;
}


static void
advinp(void)                        /* move home to current position */
{
	ipb.beg = ipb.pos;
}


static void
resetinp(void)                      /* rewind position and advance 1 */
{
	if (ipb.beg == NULL)            /* full */
		ipb.beg = ipb.end;
	ipb.pos = ipb.beg;
	ipb.chr = *ipb.pos;
	if (++ipb.beg >= &inpbuf[INBSIZ])
		ipb.beg = inpbuf;
	scaninp();
}


void
eputs(char  *msg)
{
	fputs(msg, stderr);
}


void
wputs(char  *msg)
{
	if (!nowarn)
		eputs(msg);
}


void
quit(int  code)
{
	exit(code);
}
