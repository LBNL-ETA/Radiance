#ifndef lint
static const char	RCSid[] = "$Id: rexpr.c,v 2.7 2003/06/07 12:50:20 schorsch Exp $";
#endif
/*
 * Regular expression parsing routines.
 *
 * External symbols declared in standard.h
 */

#include "copyright.h"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

/*
 * rexpr.c - regular expression parser (ala grep)
 */

#define CCHR    2
#define CDOT    4
#define CCL     6
#define NCCL    8
#define CDOL    10
#define CEOF    11
#define CBRC    14
#define CLET    15
#define STAR    01

#define ESIZE   255

#define same(a,b) (a==b || (iflag && (a^b)==' ' && isalpha(a)))

#ifdef  BSD
#define  memcpy(to,from,len)	bcopy(from,to,len)
#endif

static int     advance(), cclass();

static char    expbuf[ESIZE];
static int     iflag;
static int     circf;

int	explen;                 /* length of the last expression found */

int
ecompile(sp, iflg, wflag)               /* compile the expression */
register char  *sp;
int iflg, wflag;
{
        register c;
        register char *ep;
        char *lastep;
        int cclcnt;
        
        iflag = iflg;
        ep = expbuf;
	explen = 0;
        if (*sp == '^') {
                circf = 1;
                sp++;
        } else
        	circf = 0;
        if (wflag)
                *ep++ = CBRC;
        for (;;) {
                if (ep >= &expbuf[ESIZE])
                        return(-1);
                if ((c = *sp++) != '*')
                        lastep = ep;
                switch (c) {

                case '\0':
                        if (wflag)
                                *ep++ = CLET;
                        *ep++ = CEOF;
			explen = ep - expbuf;
                        return(0);

                case '.':
                        *ep++ = CDOT;
                        continue;

                case '*':
                        if (lastep==0)
                                goto defchar;
                        *lastep |= STAR;
                        continue;

                case '$':
                        if (*sp != '\0')
                                goto defchar;
                        *ep++ = CDOL;
                        continue;

                case '[':
                        *ep++ = CCL;
                        *ep++ = 0;
                        cclcnt = 1;
                        if ((c = *sp++) == '^') {
                                c = *sp++;
                                ep[-2] = NCCL;
                        }
                        do {
                                *ep++ = c;
                                cclcnt++;
                                if (c=='\0' || ep >= &expbuf[ESIZE])
                                        return(-1);
                        } while ((c = *sp++) != ']');
                        lastep[1] = cclcnt;
                        continue;

                case '\\':
                        if ((c = *sp++) == '\0')
                                return(-1);
                        if (c == '<') {
				if (ep == expbuf || ep[-1] != CBRC)
					*ep++ = CBRC;
                                continue;
                        }
                        if (c == '>') {
                                *ep++ = CLET;
                                continue;
                        }
                defchar:
                default:
                        *ep++ = CCHR;
                        *ep++ = c;
                }
        }
}

char *
expsave()		      /* save compiled string */
{
	register char  *ep;

	if (explen == 0)
		return(NULL);
	if ((ep = (char *)malloc(explen+3)) == NULL)
		return(NULL);
	ep[0] = iflag;
	ep[1] = circf;
	ep[2] = explen;
	(void)memcpy(ep+3, expbuf, explen);
	return(ep);
}

void
expset(ep)			/* install saved string */
register char  *ep;
{
	iflag = ep[0];
	circf = ep[1];
	(void)memcpy(expbuf, ep+3, ep[2]&0xff);
}

char *
eindex(sp)                    /* find the expression in string sp */
register char *sp;
{
	/* check for match at beginning of line, watch CBRC */
	if (advance(sp, expbuf[0]==CBRC ? expbuf+1 : expbuf))
		return(sp);
	if (circf)
                return(NULL);
        /* fast check for first character */
        if (expbuf[0]==CCHR) {
		register c = expbuf[1];
		while (*++sp)
			if (same(*sp, c) && advance(sp, expbuf))
				return(sp);
                return(NULL);
        }
        /* regular algorithm */
	while (*++sp)
                if (advance(sp, expbuf))
                        return(sp);
        return(NULL);
}

static int
advance(alp, ep)
        char *alp;
register char *ep;
{
        register char *lp;
	char *curlp;

        lp = alp;
        for (;;) switch (*ep++) {

        case CCHR:
                if (!same(*ep, *lp))
                        return (0);
                ep++, lp++;
                continue;

        case CDOT:
                if (*lp++)
                        continue;
                return(0);

        case CDOL:
                if (*lp==0)
                        continue;
                return(0);

        case CEOF:
                explen = lp - alp;
                return(1);

        case CCL:
                if (cclass(ep, *lp++, 1)) {
                        ep += *ep;
                        continue;
                }
                return(0);

        case NCCL:
                if (cclass(ep, *lp++, 0)) {
                        ep += *ep;
                        continue;
                }
                return(0);

        case CDOT|STAR:
                curlp = lp;
                while (*lp++);
                goto star;

        case CCHR|STAR:
                curlp = lp;
                while (same(*lp, *ep))
                        lp++;
                lp++;
                ep++;
                goto star;

        case CCL|STAR:
        case NCCL|STAR:
                curlp = lp;
                while (cclass(ep, *lp++, ep[-1]==(CCL|STAR)));
                ep += *ep;
        star:
                do {
                        lp--;
                        if (advance(lp, ep)) {
                                explen += lp - alp;
                                return(1);
                        }
                } while (lp > curlp);
                return(0);

        case CBRC:
                if ((isalnum(*lp) || *lp == '_') && !(isalnum(lp[-1]) || lp[-1] == '_'))
                        continue;
                return (0);

        case CLET:
                if (!isalnum(*lp) && *lp != '_')
                        continue;
                return (0);

        default:
                fprintf(stderr, "RE botch\n");
        }
}

static int
cclass(set, c, af)
register char *set;
register c;
int af;
{
        register n;

        if (c == 0)
                return(0);
        n = *set++;
        while (--n)
                if (n > 2 && set[1] == '-') {
                        if (c >= set[0] && c <= set[2])
                                return (af);
                        set += 3;
                        n -= 2;
                } else
                        if (*set++ == c)
                                return(af);
        return(!af);
}
