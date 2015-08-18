#ifndef lint
static const char RCSid[] = "$Id: g3affine.c,v 2.2 2015/08/18 15:02:53 greg Exp $";
#endif
#include <stdio.h>
#include <stdlib.h>

#include "g3affine.h"
#include "gbasic.h"
#include "g3sphere.h"

static	g3AVec		g3av_zero	= {0.0,0.0,0.0,0.0};

static	g3ATrans	g3at_unit	= {{1.0,0.0,0.0,0.0},
							{0.0,1.0,0.0,0.0},
							{0.0,0.0,1.0,0.0},
							{0.0,0.0,0.0,1.0}};

static	g3ATrans	g3at_zero	= {{0.0,0.0,0.0,0.0},
							{0.0,0.0,0.0,0.0},
							{0.0,0.0,0.0,0.0},
							{0.0,0.0,0.0,0.0}};

void	g3at_tounit(g3ATrans t) 
{
	g3at_copy(t,g3at_unit);
}

static	int 	nr_gaussj(g3Float a[4][4],g3Float b[4])
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	g3Float big,dum,pivinv;

	indxc = (int*) calloc(4,sizeof(int));
	indxr = (int*) calloc(4,sizeof(int));
	ipiv = (int*) calloc(4,sizeof(int));

	icol = irow = 0;
	for (j=0;j<4;j++) 
		ipiv[j]=0;
	for (i=0;i<4;i++) {
		big=0.0;
		for (j=0;j<4;j++)
			if (ipiv[j] != 1)
				for (k=0;k<4;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big = fabs(a[j][k]);
							irow = j;
							icol = k;
						}
					} else if (ipiv[k] > 1) {
						return 0;
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=0;l<4;l++) 
				GB_SWAP(a[irow][l],a[icol][l])
			GB_SWAP(b[irow],b[icol])
		}
		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) 
			return 0;
		pivinv = 1.0/a[icol][icol];
		a[icol][icol] = 1.0;
		for (l=0;l<4;l++) 
			a[icol][l] *= pivinv;
		b[icol] *= pivinv;
		for (ll=0;ll<4;ll++) {
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l=0;l<4;l++) 
					a[ll][l] -= a[icol][l]*dum;
				b[ll] -= b[icol]*dum;
			}
		}
	}
	for (l=3;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<4;k++)
				GB_SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free(ipiv);
	free(indxr);
	free(indxc);
	return 1;
}

void	g3at_translate(g3ATrans t,g3Vec tv)
{
	int i;
	g3at_copy(t,g3at_unit);
	for(i=0;i<3;i++)
		t[i][3] += tv[i];
}

void	g3at_add_trans(g3ATrans t,g3Vec tv)
{
	g3ATrans tt;

	g3at_translate(tt,tv);
	g3at_comb_to(t,tt,t);
}

void	g3at_prepend_trans(g3ATrans t,g3Vec tv)
{
	g3ATrans tt;

	g3at_translate(tt,tv);
	g3at_comb(t,tt);
}

void	g3at_rotate(g3ATrans t,g3Vec rnorm,g3Float ang)
{
	g3Vec axe;
	int aid,i;

	axe = g3v_create();
	g3at_copy(t,g3at_unit);

	for(aid=0;aid<3;aid++) {
		g3v_zero(axe);
		axe[aid] = 1.0;
		g3v_rotate(axe,axe,rnorm,ang);
		for(i=0;i<3;i++)
			t[i][aid] = axe[i];
	}
	g3v_free(axe);
}

void	g3at_add_rot(g3ATrans t,g3Vec rnorm,g3Float ang)
{
	g3ATrans rt;

	g3at_rotate(rt,rnorm,ang);
	g3at_comb_to(t,rt,t);
}

void	g3at_prepend_rot(g3ATrans t,g3Vec rnorm,g3Float ang)
{
	g3ATrans rt;

	g3at_rotate(rt,rnorm,ang);
	g3at_comb(t,rt);
}

void	g3at_btrans(g3ATrans t,g3Vec xv,g3Vec yv,g3Vec zv)
{
	int i;

	g3at_copy(t,g3at_unit);
	for(i=0;i<3;i++) {
		t[i][0] = xv[i];
		t[i][1] = yv[i];
		t[i][2] = zv[i];
	}
}

int		g3at_inverse(g3ATrans t)
{
	g3Float	 v[] = {0,0,0,0};

	if (!nr_gaussj(t,v)) {
		fprintf(stderr,"g3at_inverse: singular matrix\n");
		return 0;
	}
	return 1;
}

void	g3at_apply_h(g3ATrans t,g3AVec v)
{
	int i,j;
	g3AVec vt;

	g3av_copy(vt,v);
	g3av_copy(v,g3av_zero);
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			v[i] += t[i][j]*vt[j];
	g3av_tonorm(v);
}

void	g3av_print(g3AVec vh,FILE* outf)
{
	int i;

	for(i=0;i<4;i++)
		fprintf(outf,"%f ",vh[i]);
}

void	g3at_apply(g3ATrans t,g3Vec v)
{
	g3AVec res;

	g3av_vectohom(res,v);
	g3at_apply_h(t,res);
	g3v_set(v,res[0],res[1],res[2]);
}

void	g3at_iapply(g3Vec v,g3ATrans t)
{
	g3AVec res;

	g3av_copy(res,v);
	g3at_iapply_h(res,t);
	g3v_set(v,res[0],res[1],res[2]);
}

void	g3at_iapply_h(g3Vec v,g3ATrans t)
{
	int i,j;
	g3AVec vt;

	g3av_copy(vt,v);
	g3av_copy(v,g3av_zero);
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			v[i] += t[j][i]*vt[j];
	g3av_tonorm(v);
}

void	g3at_comb(g3ATrans t,g3ATrans tf)
{
	g3ATrans res;

	g3at_comb_to(res,t,tf);
	g3at_copy(t,res);
}

void	g3at_comb_to(g3ATrans res,g3ATrans t,g3ATrans tf)
{
	int i,j,k;
	g3ATrans temp;
	g3at_copy(temp,g3at_zero);
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)
			for(k=0;k<4;k++)
				temp[i][j] += t[i][k]*tf[k][j];
	g3at_copy(res,temp);
}
	

void	g3av_vectohom(g3AVec res,g3Vec v)
{
	int i;

	for(i=0;i<3;i++)
		res[i] = v[i];
	res[3] = 1.0;
}

void	g3av_homtovec(g3Vec res,g3AVec h)
{
	int i;

	g3av_tonorm(h);
	for(i=0;i<3;i++)
		res[i] = h[i];
}

int		g3av_tonorm(g3AVec v)
{
	int i;

	if (gb_epseq(v[3],0.0))
		return 0;
	for(i=0;i<3;i++)
		v[i] /= v[3];
	v[3] = 1.0;
	return 1;
}

void	g3av_copy(g3AVec v,g3AVec vsrc)
{
	int i;
	
	for(i=0;i<4;i++)
		v[i] = vsrc[i];
}

void	g3at_print(g3ATrans t,FILE* outf)
{
	int i,j;
	
	for(i=0;i<4;i++) {
		for(j=0;j<4;j++)	
			fprintf(outf,"%f ",t[i][j]);
		fprintf(outf,"\n");
	}
}

void		g3at_copy(g3ATrans t,g3ATrans tsrc)
{
	int i,j;
	
	for(i=0;i<4;i++)
		for(j=0;j<4;j++)	
			t[i][j] = tsrc[i][j];
}

int		g3at_ctrans(g3ATrans t,g3Vec vdir,g3Vec vup)
{
	g3Vec x,y,z;

	x = g3v_create();
	y = g3v_create();
	z = g3v_create();
	g3v_copy(z,vdir);
	g3v_copy(x,vup);
	g3v_normalize(z);
	g3v_normalize(x);
	g3v_cross(y,z,x);
	if (!g3v_normalize(y)) {
		fprintf(stderr,"g3at_sph_ctrans: parallel\n");
		g3v_free(x);
		g3v_free(y);
		g3v_free(z);
		return 0;
	}
	g3v_cross(x,y,z);
	if (!g3v_normalize(x)) {
		fprintf(stderr,"g3at_sph_ctrans: parallel\n");
		g3v_free(x);
		g3v_free(y);
		g3v_free(z);
		return 0;
	}
	g3at_btrans(t,x,y,z);
	g3v_free(x);
	g3v_free(y);
	g3v_free(z);
	return 1;
}

int		g3at_sph_ctrans(g3ATrans t,g3Vec vdir,g3Vec vup)
{
	g3Vec x,z;
	int res;

	x = g3v_create();
	z = g3v_create();
	g3v_copy(z,vdir);
	g3v_copy(x,vup);
	z[G3S_RAD] = x[G3S_RAD] = 1.0;
	g3s_sphtocc(x,x);
	g3s_sphtocc(z,z);
	res = g3at_ctrans(t,z,x);
	g3v_free(x);
	g3v_free(z);

	return res;
}

#ifdef GL_CONVERSIONS
void	g3at_get_openGL(g3ATrans t,GLdouble* glm)
{
	int r,c;
	
	for(r=0;r<4;r++)
		for(c=0;c<4;c++)
			glm[c+(r*4)] = t[c][r];
}

void	g3at_from_openGL(g3ATrans t,GLdouble* glm)
{
	int r,c;

	for(r=0;r<4;r++)
		for(c=0;c<4;c++)
			t[c][r] = glm[c+(r*4)];
}
#endif

#ifdef	G3AFFINETEST


int	main(int argc,char** argv)
{
	g3Vec v,iv,x,y,z;
	g3ATrans t,it;

	g3at_copy(t,g3at_unit);
	v = g3v_create();
	iv = g3v_create();
	x = g3v_create();
	y = g3v_create();
	z = g3v_create();
	g3v_set(z,1,0,0);
	g3v_set(x,atof(argv[1]),atof(argv[2]),atof(argv[3]));
	//g3v_set(y,atof(argv[4]),atof(argv[5]),atof(argv[6]));
	//g3v_set(z,atof(argv[7]),atof(argv[8]),atof(argv[9]));
	g3at_tounit(t);
	g3at_rotate(t,z,DEG2RAD(-90));
	g3v_set(z,0,0,1);
	g3at_rotate(it,z,DEG2RAD(180));
	g3at_comb_to(t,it,t);
	g3at_apply(t,x);
	g3v_print(x,stdout);
	printf("\n");
	//g3v_print(y,stdout);
	//printf("\n");
	exit(19);
	g3v_set(z,1,DEG2RAD(atof(argv[1])),DEG2RAD(atof(argv[2])));
	g3v_set(x,1,DEG2RAD(atof(argv[3])),DEG2RAD(atof(argv[4])));
	g3at_sph_ctrans(t,z,x);
	g3at_copy(it,t);
	g3at_inverse(it);
	g3at_print(t,stderr);
	fprintf(stderr,"\n");
	g3v_set(v,1,DEG2RAD(atof(argv[5])),DEG2RAD(atof(argv[6])));
	g3s_sphtocc(v,v);
	g3v_copy(iv,v);
	g3at_apply(it,iv);
	g3at_apply(t,v);
	g3s_cctosph(v,v);
	g3s_cctosph(iv,iv);
	g3s_sphwrap(v);
	g3s_sphwrap(iv);
	fprintf(stderr,"tt %f %f\n",RAD2DEG(v[1]),RAD2DEG(v[2]));
	fprintf(stderr,"ise %f %f\n",RAD2DEG(iv[1]),RAD2DEG(iv[2]));
	return EXIT_SUCCESS;
}

#endif
