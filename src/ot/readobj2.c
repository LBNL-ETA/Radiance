#ifndef lint
static const char	RCSid[] = "$Id: readobj2.c,v 2.5 2003/03/11 19:29:05 greg Exp $";
#endif
/*
 *  readobj2.c - routines for reading in object descriptions.
 */

#include  "standard.h"

#include  "object.h"

#include  "otypes.h"

#include  <ctype.h>

extern char  *fgetword();

readobj2(input, callback)	/* read in an object file or stream */
char  *input;
int  (*callback)();
{
	FILE  *popen();
	char  *fgetline();
	FILE  *infp;
	char  buf[512];
	register int  c;

	if (input == NULL) {
		infp = stdin;
		input = "standard input";
	} else if (input[0] == '!') {
		if ((infp = popen(input+1, "r")) == NULL) {
			sprintf(errmsg, "cannot execute \"%s\"", input);
			error(SYSTEM, errmsg);
		}
	} else if ((infp = fopen(input, "r")) == NULL) {
		sprintf(errmsg, "cannot open scene file \"%s\"", input);
		error(SYSTEM, errmsg);
	}
	while ((c = getc(infp)) != EOF) {
		if (isspace(c))
			continue;
		if (c == '#') {				/* comment */
			fgets(buf, sizeof(buf), infp);
		} else if (c == '!') {			/* command */
			ungetc(c, infp);
			fgetline(buf, sizeof(buf), infp);
			readobj2(buf, callback);
		} else {				/* object */
			ungetc(c, infp);
			getobject2(input, infp, callback);
		}
	}
	if (input[0] == '!')
		pclose(infp);
	else
		fclose(infp);
}


getobject2(name, fp, f)			/* read the next object */
char  *name;
FILE  *fp;
int  (*f)();
{
	char  sbuf[MAXSTR];
	OBJREC  thisobj;
					/* get modifier */
	fgetword(sbuf, MAXSTR, fp);
	thisobj.omod = OVOID;
					/* get type */
	fgetword(sbuf, MAXSTR, fp);
	if ((thisobj.otype = otype(sbuf)) < 0) {
		sprintf(errmsg, "(%s): unknown type \"%s\"", name, sbuf);
		error(USER, errmsg);
	}
					/* get identifier */
	fgetword(sbuf, MAXSTR, fp);
	thisobj.oname = sbuf;
					/* get arguments */
	if (thisobj.otype == MOD_ALIAS) {
		fscanf(fp, "%*s");
		return;
	}
	if (readfargs(&thisobj.oargs, fp) <= 0) {
		sprintf(errmsg, "(%s): bad arguments", name);
		objerror(&thisobj, USER, errmsg);
	}
	thisobj.os = NULL;
					/* call function */
	(*f)(&thisobj);
					/* free memory */
	freefargs(&thisobj.oargs);
	free_os(&thisobj);
}
