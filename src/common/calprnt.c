#ifndef lint
static const char	RCSid[] = "$Id: calprnt.c,v 2.3 2003/02/25 02:47:21 greg Exp $";
#endif
/*
 *  calprint.c - routines for printing calcomp expressions.
 */

#include "copyright.h"

#include  <stdio.h>

#include  "calcomp.h"


void
eprint(ep, fp)			/* print a parse tree */
register EPNODE  *ep;
FILE  *fp;
{
    static EPNODE  *curdef = NULL;
    register EPNODE  *ep1;

    switch (ep->type) {

	case VAR:
	    fputs(ep->v.ln->name, fp);
	    break;

	case SYM:
	    fputs(ep->v.name, fp);
	    break;

	case FUNC:
	    eprint(ep->v.kid, fp);
	    fputc('(', fp);
	    ep1 = ep->v.kid->sibling;
	    while (ep1 != NULL) {
		eprint(ep1, fp);
		if ((ep1 = ep1->sibling) != NULL)
		    fputs(", ", fp);
	    }
	    fputc(')', fp);
	    break;

	case ARG:
	    if (curdef == NULL || curdef->v.kid->type != FUNC ||
			(ep1 = ekid(curdef->v.kid, ep->v.chan)) == NULL) {
		eputs("Bad argument!\n");
		quit(1);
	    }
	    eprint(ep1, fp);
	    break;

	case NUM:
	    fprintf(fp, "%.9g", ep->v.num);
	    break;

	case UMINUS:
	    fputc('-', fp);
	    eprint(ep->v.kid, fp);
	    break;

	case CHAN:
	    fprintf(fp, "$%d", ep->v.chan);
	    break;
	
	case '=':
	case ':':
	    ep1 = curdef;
	    curdef = ep;
	    eprint(ep->v.kid, fp);
	    fputc(' ', fp);
	    fputc(ep->type, fp);
	    fputc(' ', fp);
	    eprint(ep->v.kid->sibling, fp);
	    curdef = ep1;
	    break;
	    
	case '+':
	case '-':
	case '*':
	case '/':
	case '^':
	    fputc('(', fp);
	    eprint(ep->v.kid, fp);
	    fputc(' ', fp);
	    fputc(ep->type, fp);
	    fputc(' ', fp);
	    eprint(ep->v.kid->sibling, fp);
	    fputc(')', fp);
	    break;

	default:
	    eputs("Bad expression!\n");
	    quit(1);

    }

}


void
dprint(name, fp)		/* print a definition (all if no name) */
char  *name;
FILE  *fp;
{
    register EPNODE  *ep;
    
    if (name == NULL)
	for (ep = dfirst(); ep != NULL; ep = dnext()) {
	    eprint(ep, fp);
	    fputs(";\n", fp);
	}
    else if ((ep = dlookup(name)) != NULL) {
	eprint(ep, fp);
	fputs(";\n", fp);
    } else {
	wputs(name);
	wputs(": undefined\n");
    }
}
