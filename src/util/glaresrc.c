#ifndef lint
static const char	RCSid[] = "$Id: glaresrc.c,v 2.4 2003/02/22 02:07:30 greg Exp $";
#endif
/*
 * Gather samples and compute glare sources.
 */

#include "glare.h"
#include "linregr.h"

#define vcont(vd,vu)	((vu)-(vd)<=SEPS)
#define hcont(s1,s2)	((s1)->r-(s2)->l>=-SEPS&&(s2)->r-(s1)->l>=-SEPS)

struct source	*curlist = NULL;	/* current source list */
struct source	*donelist = NULL;	/* finished sources */


struct srcspan *
newspan(l, r, v, sb)		/* allocate a new source span */
int	l, r, v;
float	*sb;
{
	register struct srcspan	*ss;
	register int	i;

	ss = (struct srcspan *)malloc(sizeof(struct srcspan));
	if (ss == NULL)
		memerr("source spans");
	ss->l = l;
	ss->r = r;
	ss->v = v;
	ss->brsum = 0.0;
	for (i = l; i < r; i++)
		ss->brsum += sb[i+hsize];
	return(ss);
}


analyze()			/* analyze our scene */
{
	int	h, v;
	int	left;
	float	*spanbr;

	spanbr = (float *)malloc((2*hsize+1)*sizeof(float));
	if (spanbr == NULL)
		memerr("view span brightness buffer");
	for (v = vsize; v >= -vsize; v--) {
		close_sources(v);
#ifndef DEBUG
		if (verbose) {
			fprintf(stderr, "%s: analyzing... %3ld%%\r",
				progname, 100L*(vsize-v)/(2*vsize));
			fflush(stderr);
		}
#endif
		getviewspan(v, spanbr);
		left = hsize + 1;
		for (h = -hsize; h <= hsize; h++) {
			if (spanbr[h+hsize] < 0.0) {	/* off view */
				if (left < h) {
					addsrcspan(newspan(left,h,v,spanbr));
					left = hsize + 1;
				}
				continue;
			}
			if (spanbr[h+hsize] > threshold) {	/* in source */
				if (left > h)
					left = h;
			} else {			/* out of source */
				if (left < h) {
					addsrcspan(newspan(left,h,v,spanbr));
					left = hsize + 1;
				}
				addindirect(h, v, spanbr[h+hsize]);
			}
		}
		if (left < h)
			addsrcspan(newspan(left,h,v,spanbr));
	}
	free((void *)spanbr);
	close_allsrcs();
}


addindirect(h, v, br)		/* add brightness to indirect illuminances */
int	h, v;
double	br;
{
	double	tanb, d;
	int	hl;
	register int	i;

	hl = hlim(v);
	if (h <= -hl) {			/* left region */
		d = (double)(-h-hl)/sampdens;
		if (d >= 1.0-FTINY)
			return;
		tanb = d/sqrt(1.0-d*d);
		for (i = 0; i < nglardirs; i++) {
			d = indirect[i].lcos - tanb*indirect[i].lsin;
			if (d > 0.0) {
				indirect[i].sum += d * br;
				indirect[i].n += d;
			}
		}
		return;
	}
	if (h >= hl) {			/* right region */
		d = (double)(-h+hl)/sampdens;
		if (d <= -1.0+FTINY)
			return;
		tanb = d/sqrt(1.0-d*d);
		for (i = 0; i < nglardirs; i++) {
			d = indirect[i].rcos - tanb*indirect[i].rsin;
			if (d > 0.0) {
				indirect[i].sum += d * br;
				indirect[i].n += d;
			}
		}
		return;
	}
					/* central region */
	for (i = 0; i < nglardirs; i++) {
		d = cos(h_theta(h,v) - indirect[i].theta);
		if (d > 0.0) {
			indirect[i].sum += d * br;
			indirect[i].n += d;
		}
	}
}


comp_thresh()			/* compute glare threshold */
{
	int	h, v;
	int	nsamps;
	double	brsum, br;

	if (verbose)
		fprintf(stderr, "%s: computing glare threshold...\n",
				progname);
	brsum = 0.0;
	nsamps = 0;
	for (v = vsize; v >= -vsize; v -= TSAMPSTEP) {
		for (h = -hsize; h <= hsize; h += TSAMPSTEP) {
			if ((br = getviewpix(h, v)) < 0.0)
				continue;
			brsum += br;
			nsamps++;
		}
	}
	if (nsamps == 0) {
		fprintf(stderr, "%s: no viewable scene!\n", progname);
		exit(1);
	}
	threshold = GLAREBR * brsum / nsamps;
	if (threshold <= FTINY) {
		fprintf(stderr, "%s: threshold zero!\n", progname);
		exit(1);
	}
	if (verbose) {
#ifdef DEBUG
		pict_stats();
#endif
		fprintf(stderr,
			"%s: threshold set to %f cd/m2 from %d samples\n",
				progname, threshold, nsamps);
	}
}


addsrcspan(nss)			/* add new source span to our list */
struct srcspan	*nss;
{
	struct source	*last, *cs, *this;
	register struct srcspan	*ss;

	cs = NULL;
	for (this = curlist; this != NULL; this = this->next) {
		for (ss = this->first; ss != NULL; ss = ss->next) {
			if (!vcont(nss->v, ss->v))
				break;
			if (hcont(ss, nss)) {
				if (cs == NULL)
					cs = this;
				else {
					last->next = this->next;
					mergesource(cs, this);
					this = last;
				}
				break;
			}
		}
		last = this;
	}
	if (cs == NULL) {
		cs = (struct source *)malloc(sizeof(struct source));
		if (cs == NULL)
			memerr("source records");
		cs->dom = 0.0;
		cs->first = NULL;
		cs->next = curlist;
		curlist = cs;
	}
	nss->next = cs->first;
	cs->first = nss;
}


mergesource(sp, ap)		/* merge source ap into source sp */
struct source	*sp, *ap;
{
	struct srcspan	head;
	register struct srcspan	*alp, *prev, *tp;

	head.next = sp->first;
	prev = &head;
	alp = ap->first;
	while (alp != NULL && prev->next != NULL) {
		if (prev->next->v > alp->v) {
			tp = alp->next;
			alp->next = prev->next;
			prev->next = alp;
			alp = tp;
		}
		prev = prev->next;
	}
	if (prev->next == NULL)
		prev->next = alp;
	sp->first = head.next;
	if (ap->dom > 0.0 && sp->dom > 0.0) {		/* sources are done */
		sp->dir[0] *= sp->dom;
		sp->dir[1] *= sp->dom;
		sp->dir[2] *= sp->dom;
		fvsum(sp->dir, sp->dir, ap->dir, ap->dom);
		normalize(sp->dir);
		sp->brt = (sp->brt*sp->dom + ap->brt*ap->dom)
				/ (sp->dom + ap->dom);
	}
	free((void *)ap);
}


close_sources(v)		/* close sources above v */
int	v;
{
	struct source	head;
	register struct source	*last, *this;

	head.next = curlist;
	last = &head;
	for (this = curlist; this != NULL; this = this->next)
		if (!vcont(v, this->first->v)) {
			last->next = this->next;
			donesource(this);
			this = last;
		} else
			last = this;
	curlist = head.next;
}


close_allsrcs()			/* done with everything */
{
	register struct source	*this, *next;

	this = curlist;
	while (this != NULL) {
		next = this->next;
		donesource(this);
		this = next;
	}
	curlist = NULL;
}


struct srcspan *
splitspan(sso, h, v, m)		/* divide source span at point */
register struct srcspan	*sso;
double	h, v, m;
{
	register struct srcspan	*ssn;
	double	d;
	int	hs;

	d = h - m*(sso->v - v);
	hs = d < 0. ? d-.5 : d+.5;
	if (sso->l >= hs)
		return(NULL);
	if (sso->r <= hs)
		return(sso);
				/* need to split it */
	ssn = (struct srcspan *)malloc(sizeof(struct srcspan));
	if (ssn == NULL)
		memerr("source spans in splitspan");
	ssn->brsum = (double)(hs - sso->l)/(sso->r - sso->l) * sso->brsum;
	sso->brsum -= ssn->brsum;
	ssn->v = sso->v;
	ssn->l = sso->l;
	ssn->r = sso->l = hs;
	return(ssn);
}


struct source *
splitsource(so)			/* divide source in two if it's big and long */
struct source	*so;
{
	LRSUM	lr;
	LRLIN	fit;
	register struct srcspan	*ss, *ssn;
	struct srcspan	*ssl, *ssnl, head;
	int	h;
	double	mh, mv;
	struct source	*sn;

	lrclear(&lr);
	for (ss = so->first; ss != NULL; ss = ss->next)
		for (h = ss->l; h < ss->r; h++)
			lrpoint(h, ss->v, &lr);
	if ((double)lr.n/(sampdens*sampdens) < SABIG)
		return(NULL);			/* too small */
	if (lrfit(&fit, &lr) < 0)
		return(NULL);			/* can't fit a line */
	if (fit.correlation < LCORR && fit.correlation > -LCORR)
		return(NULL);
	if (verbose)
		fprintf(stderr, "%s: splitting large source\n", progname);
	mh = lrxavg(&lr);
	mv = lryavg(&lr);
	sn = (struct source *)malloc(sizeof(struct source));
	if (sn == NULL)
		memerr("source records in splitsource");
	sn->dom = 0.0;
	sn->first = NULL;
	ssnl = NULL;
	head.next = so->first;
	ssl = &head;
	for (ss = so->first; ss != NULL; ssl = ss, ss = ss->next)
		if ((ssn = splitspan(ss, mh, mv, fit.slope)) != NULL) {
			if (ssn == ss) {	/* remove from old */
				ssl->next = ss->next;
				ss = ssl;
			}
			if (ssnl == NULL)	/* add to new */
				sn->first = ssn;
			else
				ssnl->next = ssn;
			ssn->next = NULL;
			ssnl = ssn;
		}
	so->first = head.next;
	return(sn);
}


donesource(sp)			/* finished with this source */
register struct source	*sp;
{
	struct source	*newsrc;
	register struct srcspan	*ss;
	int	h, n;
	double	hsum, vsum, d;

	while ((newsrc = splitsource(sp)) != NULL)	/* split it? */
		donesource(newsrc);
	sp->dom = 0.0;
	hsum = vsum = 0.0;
	sp->brt = 0.0;
	n = 0;
	for (ss = sp->first; ss != NULL; ss = ss->next) {
		sp->brt += ss->brsum;
		n += ss->r - ss->l;
		for (h = ss->l; h < ss->r; h++) {
			d = pixsize(h, ss->v);
			hsum += d*h;
			vsum += d*ss->v;
			sp->dom += d;
		}
	}
	freespans(sp);
	if (sp->dom <= FTINY) {		/* must be right at edge of image */
		free((void *)sp);
		return;
	}
	sp->brt /= (double)n;
	compdir(sp->dir, (int)(hsum/sp->dom), (int)(vsum/sp->dom));
	sp->next = donelist;
	donelist = sp;
	if (verbose)
		fprintf(stderr,
	"%s: source at [%.3f,%.3f,%.3f], dw %.5f, br %.1f (%d samps)\n",
			progname, sp->dir[0], sp->dir[1], sp->dir[2], 
			sp->dom, sp->brt, n);
}


struct source *
findbuddy(s, l)			/* find close enough source to s in l*/
register struct source	*s, *l;
{
	struct source	*bestbuddy = NULL;
	double	d, r, mindist = MAXBUDDY;

	r = sqrt(s->dom/PI);
	for ( ; l != NULL; l = l->next) {
		d = sqrt(dist2(l->dir, s->dir)) - sqrt(l->dom/PI) - r;
		if (d < mindist) {
			bestbuddy = l;
			mindist = d;
		}
	}
	return(bestbuddy);
}


absorb_specks()			/* eliminate too-small sources */
{
	struct source	head, *buddy;
	register struct source	*last, *this;

	if (verbose)
		fprintf(stderr, "%s: absorbing small sources...\n", progname);
	head.next = donelist;
	last = &head;
	for (this = head.next; this != NULL; this = this->next)
		if (TOOSMALL(this)) {
			last->next = this->next;
			buddy = findbuddy(this, head.next);
			if (buddy != NULL)
				mergesource(buddy, this);
			else
				absorb(this);
			this = last;
		} else
			last = this;
	donelist = head.next;
}


absorb(s)			/* absorb a source into indirect */
register struct source	*s;
{
	FVECT	dir;
	double	d;
	register int	i;

	for (i = 0; i < nglardirs; i++) {
		spinvector(dir, ourview.vdir, ourview.vup, indirect[i].theta);
		d = DOT(dir,s->dir)*s->dom*(sampdens*sampdens);
		if (d <= 0.0)
			continue;
		indirect[i].sum += d * s->brt;
		indirect[i].n += d;
	}
	freespans(s);
	free((void *)s);
}


freespans(sp)			/* free spans associated with source */
struct source	*sp;
{
	register struct srcspan	*ss;

	while ((ss = sp->first) != NULL) {
		sp->first = ss->next;
		free((void *)ss);
	}
}
