#ifndef lint
static const char RCSid[] = "$Id: bsdfrep.c,v 2.5 2012/11/07 03:04:23 greg Exp $";
#endif
/*
 * Support BSDF representation as radial basis functions.
 *
 *	G. Ward
 */

#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rtio.h"
#include "resolu.h"
#include "bsdfrep.h"
				/* active grid resolution */
int			grid_res = GRIDRES;

				/* coverage/symmetry using INP_QUAD? flags */
int			inp_coverage = 0;
				/* all incident angles in-plane so far? */
int			single_plane_incident = -1;

				/* input/output orientations */
int			input_orient = 0;
int			output_orient = 0;

				/* processed incident DSF measurements */
RBFNODE			*dsf_list = NULL;

				/* RBF-linking matrices (edges) */
MIGRATION		*mig_list = NULL;

				/* current input direction */
double			theta_in_deg, phi_in_deg;

/* Register new input direction */
int
new_input_direction(double new_theta, double new_phi)
{
	if (!input_orient)		/* check input orientation */
		input_orient = 1 - 2*(new_theta > 90.);
	else if (input_orient > 0 ^ new_theta < 90.) {
		fprintf(stderr,
		"%s: Cannot handle input angles on both sides of surface\n",
				progname);
		return(0);
	}
					/* normalize angle ranges */
	while (new_theta < -180.)
		new_theta += 360.;
	while (new_theta > 180.)
		new_theta -= 360.;
	if (new_theta < 0) {
		new_theta = -new_theta;
		new_phi += 180.;
	}
	while (new_phi < 0)
		new_phi += 360.;
	while (new_phi >= 360.)
		new_phi -= 360.;
	if (single_plane_incident > 0)	/* check input coverage */
		single_plane_incident = (round(new_phi) == round(phi_in_deg));
	else if (single_plane_incident < 0)
		single_plane_incident = 1;
	theta_in_deg = new_theta;	/* assume it's OK */
	phi_in_deg = new_phi;
	if ((1. < new_phi) & (new_phi < 89.))
		inp_coverage |= INP_QUAD1;
	else if ((91. < new_phi) & (new_phi < 179.))
		inp_coverage |= INP_QUAD2;
	else if ((181. < new_phi) & (new_phi < 269.))
		inp_coverage |= INP_QUAD3;
	else if ((271. < new_phi) & (new_phi < 359.))
		inp_coverage |= INP_QUAD4;
	return(1);
}

/* Apply symmetry to the given vector based on distribution */
int
use_symmetry(FVECT vec)
{
	double	phi = get_phi360(vec);

	switch (inp_coverage) {
	case INP_QUAD1|INP_QUAD2|INP_QUAD3|INP_QUAD4:
		break;
	case INP_QUAD1|INP_QUAD2:
		if ((-FTINY > phi) | (phi > 180.+FTINY))
			goto mir_y;
		break;
	case INP_QUAD2|INP_QUAD3:
		if ((90.-FTINY > phi) | (phi > 270.+FTINY))
			goto mir_x;
		break;
	case INP_QUAD3|INP_QUAD4:
		if ((180.-FTINY > phi) | (phi > 360.+FTINY))
			goto mir_y;
		break;
	case INP_QUAD4|INP_QUAD1:
		if ((270.-FTINY > phi) & (phi > 90.+FTINY))
			goto mir_x;
		break;
	case INP_QUAD1:
		if ((-FTINY > phi) | (phi > 90.+FTINY))
			switch ((int)(phi*(1./90.))) {
			case 1: goto mir_x;
			case 2: goto mir_xy;
			case 3: goto mir_y;
			}
		break;
	case INP_QUAD2:
		if ((90.-FTINY > phi) | (phi > 180.+FTINY))
			switch ((int)(phi*(1./90.))) {
			case 0: goto mir_x;
			case 2: goto mir_y;
			case 3: goto mir_xy;
			}
		break;
	case INP_QUAD3:
		if ((180.-FTINY > phi) | (phi > 270.+FTINY))
			switch ((int)(phi*(1./90.))) {
			case 0: goto mir_xy;
			case 1: goto mir_y;
			case 3: goto mir_x;
			}
		break;
	case INP_QUAD4:
		if ((270.-FTINY > phi) | (phi > 360.+FTINY))
			switch ((int)(phi*(1./90.))) {
			case 0: goto mir_y;
			case 1: goto mir_xy;
			case 2: goto mir_x;
			}
		break;
	default:
		fprintf(stderr, "%s: Illegal input coverage (%d)\n",
					progname, inp_coverage);
		exit(1);
	}
	return(0);		/* in range */
mir_x:
	vec[0] = -vec[0];
	return(MIRROR_X);
mir_y:
	vec[1] = -vec[1];
	return(MIRROR_Y);
mir_xy:
	vec[0] = -vec[0];
	vec[1] = -vec[1];
	return(MIRROR_X|MIRROR_Y);
}

/* Reverse symmetry based on what was done before */
void
rev_symmetry(FVECT vec, int sym)
{
	if (sym & MIRROR_X)
		vec[0] = -vec[0];
	if (sym & MIRROR_Y)
		vec[1] = -vec[1];
}

/* Reverse symmetry for an RBF distribution */
void
rev_rbf_symmetry(RBFNODE *rbf, int sym)
{
	int	n;

	rev_symmetry(rbf->invec, sym);
	if (sym & MIRROR_X)
		for (n = rbf->nrbf; n-- > 0; )
			rbf->rbfa[n].gx = grid_res-1 - rbf->rbfa[n].gx;
	if (sym & MIRROR_Y)
		for (n = rbf->nrbf; n-- > 0; )
			rbf->rbfa[n].gy = grid_res-1 - rbf->rbfa[n].gy;
}

/* Compute volume associated with Gaussian lobe */
double
rbf_volume(const RBFVAL *rbfp)
{
	double	rad = R2ANG(rbfp->crad);

	return((2.*M_PI) * rbfp->peak * rad*rad);
}

/* Compute outgoing vector from grid position */
void
ovec_from_pos(FVECT vec, int xpos, int ypos)
{
	double	uv[2];
	double	r2;
	
	SDsquare2disk(uv, (1./grid_res)*(xpos+.5), (1./grid_res)*(ypos+.5));
				/* uniform hemispherical projection */
	r2 = uv[0]*uv[0] + uv[1]*uv[1];
	vec[0] = vec[1] = sqrt(2. - r2);
	vec[0] *= uv[0];
	vec[1] *= uv[1];
	vec[2] = output_orient*(1. - r2);
}

/* Compute grid position from normalized input/output vector */
void
pos_from_vec(int pos[2], const FVECT vec)
{
	double	sq[2];		/* uniform hemispherical projection */
	double	norm = 1./sqrt(1. + fabs(vec[2]));

	SDdisk2square(sq, vec[0]*norm, vec[1]*norm);

	pos[0] = (int)(sq[0]*grid_res);
	pos[1] = (int)(sq[1]*grid_res);
}

/* Evaluate RBF for DSF at the given normalized outgoing direction */
double
eval_rbfrep(const RBFNODE *rp, const FVECT outvec)
{
	double		res = .0;
	const RBFVAL	*rbfp;
	FVECT		odir;
	double		sig2;
	int		n;

	if (rp == NULL)
		return(.0);
	rbfp = rp->rbfa;
	for (n = rp->nrbf; n--; rbfp++) {
		ovec_from_pos(odir, rbfp->gx, rbfp->gy);
		sig2 = R2ANG(rbfp->crad);
		sig2 = (DOT(odir,outvec) - 1.) / (sig2*sig2);
		if (sig2 > -19.)
			res += rbfp->peak * exp(sig2);
	}
	return(res);
}

/* Insert a new directional scattering function in our global list */
int
insert_dsf(RBFNODE *newrbf)
{
	RBFNODE		*rbf, *rbf_last;
	int		pos;
					/* check for redundant meas. */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next)
		if (DOT(rbf->invec, newrbf->invec) >= 1.-FTINY) {
			fprintf(stderr,
				"%s: Duplicate incident measurement (ignored)\n",
					progname);
			free(newrbf);
			return(-1);
		}
					/* keep in ascending theta order */
	for (rbf_last = NULL, rbf = dsf_list; rbf != NULL;
					rbf_last = rbf, rbf = rbf->next)
		if (single_plane_incident && input_orient*rbf->invec[2] <
						input_orient*newrbf->invec[2])
			break;
	if (rbf_last == NULL) {		/* insert new node in list */
		newrbf->ord = 0;
		newrbf->next = dsf_list;
		dsf_list = newrbf;
	} else {
		newrbf->ord = rbf_last->ord + 1;
		newrbf->next = rbf;
		rbf_last->next = newrbf;
	}
	rbf_last = newrbf;
	while (rbf != NULL) {		/* update ordinal positions */
		rbf->ord = rbf_last->ord + 1;
		rbf_last = rbf;
		rbf = rbf->next;
	}
	return(newrbf->ord);
}

/* Get the DSF indicated by its ordinal position */
RBFNODE *
get_dsf(int ord)
{
	RBFNODE		*rbf;

	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next)
		if (rbf->ord == ord)
			return(rbf);
	return(NULL);
}

/* Get triangle surface orientation (unnormalized) */
void
tri_orient(FVECT vres, const FVECT v1, const FVECT v2, const FVECT v3)
{
	FVECT	v2minus1, v3minus2;

	VSUB(v2minus1, v2, v1);
	VSUB(v3minus2, v3, v2);
	VCROSS(vres, v2minus1, v3minus2);
}

/* Determine if vertex order is reversed (inward normal) */
int
is_rev_tri(const FVECT v1, const FVECT v2, const FVECT v3)
{
	FVECT	tor;

	tri_orient(tor, v1, v2, v3);

	return(DOT(tor, v2) < 0.);
}

/* Find vertices completing triangles on either side of the given edge */
int
get_triangles(RBFNODE *rbfv[2], const MIGRATION *mig)
{
	const MIGRATION	*ej1, *ej2;
	RBFNODE		*tv;

	rbfv[0] = rbfv[1] = NULL;
	if (mig == NULL)
		return(0);
	for (ej1 = mig->rbfv[0]->ejl; ej1 != NULL;
				ej1 = nextedge(mig->rbfv[0],ej1)) {
		if (ej1 == mig)
			continue;
		tv = opp_rbf(mig->rbfv[0],ej1);
		for (ej2 = tv->ejl; ej2 != NULL; ej2 = nextedge(tv,ej2))
			if (opp_rbf(tv,ej2) == mig->rbfv[1]) {
				rbfv[is_rev_tri(mig->rbfv[0]->invec,
						mig->rbfv[1]->invec,
						tv->invec)] = tv;
				break;
			}
	}
	return((rbfv[0] != NULL) + (rbfv[1] != NULL));
}

/* Clear our BSDF representation and free memory */
void
clear_bsdf_rep(void)
{
	while (mig_list != NULL) {
		MIGRATION	*mig = mig_list;
		mig_list = mig->next;
		free(mig);
	}
	while (dsf_list != NULL) {
		RBFNODE		*rbf = dsf_list;
		dsf_list = rbf->next;
		free(rbf);
	}
	inp_coverage = 0;
	single_plane_incident = -1;
	input_orient = output_orient = 0;
	grid_res = GRIDRES;
}

/* Write our BSDF mesh interpolant out to the given binary stream */
void
save_bsdf_rep(FILE *ofp)
{
	RBFNODE		*rbf;
	MIGRATION	*mig;
	int		i, n;
					/* finish header */
	fprintf(ofp, "SYMMETRY=%d\n", !single_plane_incident * inp_coverage);
	fprintf(ofp, "IO_SIDES= %d %d\n", input_orient, output_orient);
	fprintf(ofp, "GRIDRES=%d\n", grid_res);
	fputformat(BSDFREP_FMT, ofp);
	fputc('\n', ofp);
					/* write each DSF */
	for (rbf = dsf_list; rbf != NULL; rbf = rbf->next) {
		putint(rbf->ord, 4, ofp);
		putflt(rbf->invec[0], ofp);
		putflt(rbf->invec[1], ofp);
		putflt(rbf->invec[2], ofp);
		putflt(rbf->vtotal, ofp);
		putint(rbf->nrbf, 4, ofp);
		for (i = 0; i < rbf->nrbf; i++) {
			putflt(rbf->rbfa[i].peak, ofp);
			putint(rbf->rbfa[i].crad, 2, ofp);
			putint(rbf->rbfa[i].gx, 1, ofp);
			putint(rbf->rbfa[i].gy, 1, ofp);
		}
	}
	putint(-1, 4, ofp);		/* terminator */
					/* write each migration matrix */
	for (mig = mig_list; mig != NULL; mig = mig->next) {
		int	zerocnt = 0;
		putint(mig->rbfv[0]->ord, 4, ofp);
		putint(mig->rbfv[1]->ord, 4, ofp);
					/* write out as sparse data */
		n = mtx_nrows(mig) * mtx_ncols(mig);
		for (i = 0; i < n; i++) {
			if (zerocnt == 0xff) {
				putint(0xff, 1, ofp); zerocnt = 0;
			}
			if (mig->mtx[i] != 0) {
				putint(zerocnt, 1, ofp); zerocnt = 0;
				putflt(mig->mtx[i], ofp);
			} else
				++zerocnt;
		}
		putint(zerocnt, 1, ofp);
	}
	putint(-1, 4, ofp);		/* terminator */
	putint(-1, 4, ofp);
	if (fflush(ofp) == EOF) {
		fprintf(stderr, "%s: error writing BSDF interpolant\n",
				progname);
		exit(1);
	}
}

/* Check header line for critical information */
static int
headline(char *s, void *p)
{
	char	fmt[32];

	if (!strncmp(s, "SYMMETRY=", 9)) {
		inp_coverage = atoi(s+9);
		single_plane_incident = !inp_coverage;
		return(0);
	}
	if (!strncmp(s, "IO_SIDES=", 9)) {
		sscanf(s+9, "%d %d", &input_orient, &output_orient);
		return(0);
	}
	if (!strncmp(s, "GRIDRES=", 8)) {
		sscanf(s+8, "%d", &grid_res);
		return(0);
	}
	if (formatval(fmt, s) && strcmp(fmt, BSDFREP_FMT))
		return(-1);
	return(0);
}

/* Read a BSDF mesh interpolant from the given binary stream */
int
load_bsdf_rep(FILE *ifp)
{
	RBFNODE		rbfh;
	int		from_ord, to_ord;
	int		i;

	clear_bsdf_rep();
	if (ifp == NULL)
		return(0);
	if (getheader(ifp, headline, NULL) < 0 || single_plane_incident < 0 |
			!input_orient | !output_orient) {
		fprintf(stderr, "%s: missing/bad format for BSDF interpolant\n",
				progname);
		return(0);
	}
	rbfh.next = NULL;		/* read each DSF */
	rbfh.ejl = NULL;
	while ((rbfh.ord = getint(4, ifp)) >= 0) {
		RBFNODE		*newrbf;

		rbfh.invec[0] = getflt(ifp);
		rbfh.invec[1] = getflt(ifp);
		rbfh.invec[2] = getflt(ifp);
		rbfh.vtotal = getflt(ifp);
		rbfh.nrbf = getint(4, ifp);
		newrbf = (RBFNODE *)malloc(sizeof(RBFNODE) +
					sizeof(RBFVAL)*(rbfh.nrbf-1));
		if (newrbf == NULL)
			goto memerr;
		memcpy(newrbf, &rbfh, sizeof(RBFNODE));
		for (i = 0; i < rbfh.nrbf; i++) {
			newrbf->rbfa[i].peak = getflt(ifp);
			newrbf->rbfa[i].crad = getint(2, ifp) & 0xffff;
			newrbf->rbfa[i].gx = getint(1, ifp) & 0xff;
			newrbf->rbfa[i].gy = getint(1, ifp) & 0xff;
		}
		if (feof(ifp))
			goto badEOF;
					/* insert in global list */
		if (insert_dsf(newrbf) != rbfh.ord) {
			fprintf(stderr, "%s: error adding DSF\n", progname);
			return(0);
		}
	}
					/* read each migration matrix */
	while ((from_ord = getint(4, ifp)) >= 0 &&
			(to_ord = getint(4, ifp)) >= 0) {
		RBFNODE		*from_rbf = get_dsf(from_ord);
		RBFNODE		*to_rbf = get_dsf(to_ord);
		MIGRATION	*newmig;
		int		n;

		if ((from_rbf == NULL) | (to_rbf == NULL)) {
			fprintf(stderr,
				"%s: bad DSF reference in migration edge\n",
					progname);
			return(0);
		}
		n = from_rbf->nrbf * to_rbf->nrbf;
		newmig = (MIGRATION *)malloc(sizeof(MIGRATION) +
						sizeof(float)*(n-1));
		if (newmig == NULL)
			goto memerr;
		newmig->rbfv[0] = from_rbf;
		newmig->rbfv[1] = to_rbf;
		memset(newmig->mtx, 0, sizeof(float)*n);
		for (i = 0; ; ) {	/* read sparse data */
			int	zc = getint(1, ifp) & 0xff;
			if ((i += zc) >= n)
				break;
			if (zc == 0xff)
				continue;
			newmig->mtx[i++] = getflt(ifp);
		}
		if (feof(ifp))
			goto badEOF;
					/* insert in edge lists */
		newmig->enxt[0] = from_rbf->ejl;
		from_rbf->ejl = newmig;
		newmig->enxt[1] = to_rbf->ejl;
		to_rbf->ejl = newmig;
					/* push onto global list */
		newmig->next = mig_list;
		mig_list = newmig;
	}
	return(1);			/* success! */
memerr:
	fprintf(stderr, "%s: Out of memory in load_bsdf_rep()\n", progname);
	exit(1);
badEOF:
	fprintf(stderr, "%s: Unexpected EOF in load_bsdf_rep()\n", progname);
	return(0);
}
