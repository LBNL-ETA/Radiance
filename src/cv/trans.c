#ifndef lint
static const char	RCSid[] = "$Id: trans.c,v 2.3 2003/04/23 00:52:33 greg Exp $";
#endif
/*
 * Translator utilities
 *
 *	Greg Ward
 */

#include <stdio.h>

#include <stdlib.h>

#include "trans.h"


fgetid(idp, dls, fp)		/* read an id up to char in dls from fp */
ID	*idp;
char	*dls;
register FILE	*fp;
{
	char	dset[256/8];
	char	buf[MAXSTR];
	register int	c;
	register char	*cp;
					/* create delimiter set */
	for (cp = dset+sizeof(dset); cp-- > dset; )
		*cp = 0;
	for (cp = dls; *cp; cp++)
		dset[*cp>>3] |= 1<<(*cp&7);
					/* get characters up to delimiter */
	cp = buf;
	while ((c = getc(fp)) != EOF && !(dset[c>>3] & 1<<(c&7)))
		*cp++ = c;
					/* check for empty */
	if (cp == buf) {
		if (c == EOF)
			return(EOF);
		idp->number = 0;
		idp->name = NULL;
		return(0);
	}
					/* truncate space at end */
	while (cp[-1] == ' ' || cp[-1] == '\t')
		cp--;
	*cp = '\0';
					/* skip white space at beginning */
	for (cp = buf; *cp && (*cp == ' ' || *cp == '\t'); cp++)
		;
					/* assign value */
	idp->number = atoi(cp);
	idp->name = NULL;
					/* check for ghost string */
	if (!*cp)
		return(0);
					/* check for name */
	if (idp->number == 0 && *cp != '0')
		idp->name = savestr(cp);
	return(0);
}


int
findid(idl, idp, insert)		/* find (or insert) id in list */
register IDLIST	*idl;
ID	*idp;
int	insert;
{
	int  upper, lower;
	register int  cm, i;
					/* binary search */
	lower = 0;
	upper = cm = idl->nids;
	while ((i = (lower + upper) >> 1) != cm) {
		cm = idcmp(idp, &idl->id[i]);
		if (cm > 0)
			lower = i;
		else if (cm < 0)
			upper = i;
		else
			return(i);
		cm = i;
	}
	if (!insert)
    		return(-1);
    	if (idl->nids == 0) {			/* create new list */
    		idl->id = (ID *)malloc(sizeof(ID));
    		if (idl->id == NULL)
    			goto memerr;
    	} else {				/* grow old list */
    		idl->id = (ID *)realloc((void *)idl->id,(idl->nids+1)*sizeof(ID));
    		if (idl->id == NULL)
    			goto memerr;
    		for (i = idl->nids; i > upper; i--) {
    			idl->id[i].number = idl->id[i-1].number;
    			idl->id[i].name = idl->id[i-1].name;
    		}
    	}
    	idl->nids++;				/* insert new element */
	idl->id[i].number = idp->number;
    	if (idp->name == NULL)
    		idl->id[i].name = NULL;
    	else
    		idl->id[i].name = savestr(idp->name);
	return(i);
memerr:
	eputs("Out of memory in findid\n");
	quit(1);
}


int
idcmp(id1, id2)				/* compare two identifiers */
register ID	*id1, *id2;
{
					/* names are greater than numbers */
	if (id1->name == NULL)
		if (id2->name == NULL)
			return(id1->number - id2->number);
		else
			return(-1);
	else
		if (id2->name == NULL)
			return(1);
		else
			return(strcmp(id1->name, id2->name));
}


write_quals(qlp, idl, fp)	/* write out qualifier lists */
QLIST	*qlp;
IDLIST	idl[];
FILE	*fp;
{
	int	i;
	
	for (i = 0; i < qlp->nquals; i++)
		if (idl[i].nids > 0) {
			fprintf(fp, "qualifier %s begin\n", qlp->qual[i]);
			fputidlist(&idl[i], fp);
			fprintf(fp, "end\n");
		}
}


fputidlist(qp, fp)		/* put id list out to fp */
IDLIST	*qp;
FILE	*fp;
{
	int	fi;
	register int	i;
					/* print numbers/ranges */
	fi = 0;
	for (i = 0; i < qp->nids && qp->id[i].name == NULL; i++)
		if (i > 0 && qp->id[i].number > qp->id[i-1].number+1) {
			if (i > fi+1)
				fprintf(fp, "[%d:%d]\n", qp->id[fi].number,
						qp->id[i-1].number);
			else
				fprintf(fp, "%d\n", qp->id[fi].number);
			fi = i;
		}
	if (i-1 > fi)
		fprintf(fp, "[%d:%d]\n", qp->id[fi].number,
				qp->id[i-1].number);
	else if (i > 0)
		fprintf(fp, "%d\n", qp->id[fi].number);
	for ( ; i < qp->nids; i++)
		fprintf(fp, "\"%s\"\n", qp->id[i].name);
}


RULEHD *
getmapping(file, qlp)		/* read in mapping file */
char	*file;
QLIST	*qlp;
{
	char	*err;
	register int	c;
	RULEHD	*mp = NULL;
	int	nrules = 0;
	register RULEHD	*rp;
	register int	qt;
	char	buf[MAXSTR];
	FILE	*fp;
	
	if ((fp = fopen(file, "r")) == NULL) {
		eputs(file);
		eputs(": cannot open\n");
		quit(1);
	}
							/* get each rule */
	while (fscanf(fp, " %[^ 	;\n]", buf) == 1) {
		if (buf[0] == '#') {			/* comment */
			while ((c = getc(fp)) != EOF && c != '\n')
				;
			continue;
		}
		rp = (RULEHD *)calloc(1, rulsiz(qlp->nquals));
		if (rp == NULL)
			goto memerr;
		rp->mnam = savestr(buf);
		for ( ; ; ) {				/* get conditions */
			while ((c = getc(fp)) != '(')
				if (c == ';' || c == EOF)
					goto endloop;
			if (fscanf(fp, " %s ", buf) != 1) {
				err = "missing variable";
				goto fmterr;
			}
			if ((qt = qtype(buf, qlp)) == -1) {
				err = "unknown variable";
				goto fmterr;
			}
			if (rp->qflg & FL(qt)) {
				err = "variable repeated";
				goto fmterr;
			}
			rp->qflg |= FL(qt);
			c = getc(fp);
			switch (c) {
			case '"':			/* id name */
				if (fscanf(fp, "%[^\"]\" )", buf) != 1) {
					err = "bad string value";
					goto fmterr;
				}
				idm(rp)[qt].nam = savestr(buf);
				break;
			case '[':			/* id range */
				if (fscanf(fp, "%d : %d ] )", &idm(rp)[qt].min,
						&idm(rp)[qt].max) != 2) {
					err = "bad range value";
					goto fmterr;
				}
				if (idm(rp)[qt].min > idm(rp)[qt].max) {
					err = "reverse range value";
					goto fmterr;
				}
				break;
			default:			/* id number? */
				if ((c < '0' || c > '9') && c != '-') {
					err = "unrecognizable value";
					goto fmterr;
				}
				ungetc(c, fp);
				if (fscanf(fp, "%d )", &idm(rp)[qt].min) != 1) {
					err = "bad number id";
					goto fmterr;
				}
				idm(rp)[qt].max = idm(rp)[qt].min;
				break;
			}
		}
	endloop:
		rp->next = mp;
		mp = rp;
		nrules++;
	}
	fclose(fp);
	return(mp);
fmterr:
	sprintf(buf, "%s: %s for rule %d\n", file, err, nrules+1);
	eputs(buf);
	quit(1);
memerr:
	eputs("Out of memory in getmapping\n");
	quit(1);
}


int
qtype(qnm, qlp)			/* return number for qualifier name */
char	*qnm;
register QLIST	*qlp;
{
	register int	i;
	
	for (i = 0; i < qlp->nquals; i++)
		if (!strcmp(qnm, qlp->qual[i]))
			return(i);
	return(-1);
}


matchid(it, im)			/* see if we match an id */
register ID	*it;
register IDMATCH	*im;
{
	if (it->name == NULL) {
		if (im->nam != NULL)
			return(0);
		return(it->number >= im->min && it->number <= im->max);
	}
	if (im->nam == NULL)
		return(0);
	return(!strcmp(it->name, im->nam));
}
