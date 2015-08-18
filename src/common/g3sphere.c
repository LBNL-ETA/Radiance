#ifndef lint
static const char RCSid[] = "$Id: g3sphere.c,v 2.2 2015/08/18 15:02:53 greg Exp $";
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "g3sphere.h"

static	g3Float	get_equator_rad(g3Vec cc,int c1,int c2)
{
	g3Float res;

	if (gb_epseq(cc[c2],0.0)) {
		if (gb_epseq(cc[c1],0.0)) {
			res = 0;
		} else {
			res = (cc[c1] > 0.0) ? M_PI/2.0 : -M_PI/2.0;
		}
	} else {
		res = (cc[c2] < 0.0) ? M_PI - fabs(atan(cc[c1]/cc[c2])) :
										fabs(atan(cc[c1]/cc[c2]));
		if (cc[c1] < 0.0)
			res *= -1.0;
	}
	return res;
}



g3Vec	g3s_cctomtr(g3Vec res,g3Vec cc)
{
	int copy = 0;
	if (res == cc) {
		res = g3v_create();
		copy = 1;
	}
	res[G3S_RAD] = g3v_length(cc);
	res[G3S_MY] = res[G3S_RAD]*get_equator_rad(cc,1,0);
	if (cc[2] >= res[G3S_RAD])
		res[G3S_MZ] = atanh((1.0 - GB_EPSILON)/res[G3S_RAD]);
	else
		res[G3S_MZ] =  res[G3S_RAD]*atanh(cc[2]/res[G3S_RAD]);
	if (copy) {
		g3v_copy(cc,res);
		g3v_free(res);
		res = cc;
	}
	return res;
}

g3Vec	g3s_mtrtocc(g3Vec res,g3Vec mtr)
{
	g3Float r = mtr[G3S_RAD];

	res[0] = r*cos(mtr[G3S_MY]/r)/cosh(mtr[G3S_MZ]/r);
	res[1] = r*sin(mtr[G3S_MY]/r)/cosh(mtr[G3S_MZ]/r);
	res[2] = r*tanh(mtr[G3S_MZ]/r);
	return res;
}

g3Vec	g3s_cctotr(g3Vec res,g3Vec cc)
{
	g3Float len;
	G3S_RESCOPY(res,cc);
	if (gb_epseq((res[G3S_RAD] = g3v_length(cc)),0)) {
		fprintf(stderr,"g3s_cctotr: zero vector\n");
		G3S_RESFREE(res,cc);
		return NULL;
	}
	g3v_copy(res,cc);
	res[1] = 0;
	len = g3v_length(res);
	if (gb_epseq(len,0.0)) {
		res[G3S_THETA] = M_PI/2.0;
		res[G3S_PHI] = gb_signum(cc[1])*M_PI/2.0;
	} else {
		res[G3S_THETA] = acos(cc[2]/len);
		res[G3S_PHI] = atan(cc[1]/len);
		if (cc[0] < 0.0)
			res[G3S_THETA] *= -1.0;
	}
	res[G3S_RAD] = g3v_length(cc);
	g3s_trwrap(res);
	G3S_RESFREE(res,cc);
	return res;
}
	
g3Vec	g3s_trtocc(g3Vec res,g3Vec tr)
{
	g3Vec v;
	G3S_RESCOPY(res,tr);
	v = g3v_create();
	g3v_set(v,0,-1,0);
	g3v_set(res,0,0,1);
	g3v_rotate(res,res,v,tr[G3S_THETA]);
	g3v_cross(v,res,v);
	g3v_rotate(res,res,v,tr[G3S_PHI]);
	g3v_free(v);
	g3v_scale(res,res,tr[G3S_RAD]);
	G3S_RESFREE(res,tr);
	return res;
}

g3Vec	g3s_sphtotr(g3Vec res,g3Vec sph)
{
	g3s_sphtocc(res,sph);
	return g3s_cctotr(res,res);
}

g3Vec	g3s_trtosph(g3Vec res,g3Vec tr)
{
	g3s_trtocc(res,tr);
	return g3s_cctosph(res,res);
}

g3Vec	g3s_sphtocc(g3Vec res,g3Vec sph)
{
 	g3Float r,s2,t;
    r = sph[G3S_RAD];
    s2 = sin(sph[G3S_THETA]);
    t = sph[G3S_THETA];
    res[0] = r*cos(sph[G3S_PHI])*s2;
    res[1] = r*sin(sph[G3S_PHI])*s2;
    res[2] = r*cos(t);
	return res;
}

g3Vec	g3s_cctosph(g3Vec res,g3Vec cc)
{
	int copy = 0;
	if (res == cc) {
		res = g3v_create();
		copy = 1;
	}
	res[G3S_RAD] = g3v_length(cc);
	res[G3S_THETA] = acos(cc[2]/res[G3S_RAD]);
	res[G3S_PHI] = get_equator_rad(cc,1,0);
	if (copy) {
		g3v_copy(cc,res);
		g3v_free(res);
		res = cc;
	}
	return res;
}

g3Vec	g3s_sphwrap(g3Vec sph)
{
	sph[G3S_THETA] = fmod((fmod(sph[G3S_THETA],2.0*M_PI) + 2.0*M_PI),2.0*M_PI);
	sph[G3S_PHI] = fmod((fmod(sph[G3S_PHI],2.0*M_PI) + 2.0*M_PI),2.0*M_PI);
	return sph;
}

g3Vec	g3s_trwrap(g3Vec tr)
{
	tr[G3S_THETA] = fmod((fmod(tr[G3S_THETA],2.0*M_PI) + 2.0*M_PI),2.0*M_PI);
	tr[G3S_PHI] += M_PI;
	tr[G3S_PHI] = fmod((fmod(tr[G3S_PHI],2.0*M_PI) + 2.0*M_PI),2.0*M_PI);
	tr[G3S_PHI] -= M_PI;
	return tr;
}

g3Float	g3s_dist(const g3Vec cc1,const g3Vec cc2)
{
	return acos(g3v_dot(cc1,cc2));
}

g3Float	g3s_dist_norm(const g3Vec cc1,const g3Vec cc2)
{
	return acos(g3v_dot(g3v_normalize(cc1),g3v_normalize(cc2)));
}




#ifdef G3SPHERE_TEST

int	main(int argc,char** argv)
{
	g3Vec a,b,e,z,ang;
	int i;
	g3Float vo,v,vr;

	if (argc < 4) {
		fprintf(stderr,"usage: %s <s | e> <x1> <y1> <z1>\n",argv[0]);
		return EXIT_FAILURE;
	}
	
	a = g3v_create();
	e = g3v_create();
	if (!strcmp(argv[1],"s")) {
		g3v_set(a,1.0,DEG2RAD(atof(argv[2])),DEG2RAD(atof(argv[3])));
		g3v_print(g3s_sphtocc(a,a),stdout);;
		//g3s_trtosph(a,a);
		//g3v_print(a,stdout);
		//printf("%f %f",RAD2DEG(a[1]),RAD2DEG(a[2]));
		printf("\n");
	} else {
		g3v_set(a,atof(argv[2]),atof(argv[3]),atof(argv[4]));
		g3s_cctosph(a,a);
		printf("%f %f",RAD2DEG(a[1]),RAD2DEG(a[2]));
		printf("\n");
	}
	exit(1);
	for(i=0;i<10000;i++) {
		a[0] = 1.0;	
		a[1] = M_PI/2.0*rand()/RAND_MAX;
		a[2] = 2.0*M_PI*rand()/RAND_MAX;
		g3s_sphtocc(a,a);
		g3s_cctotr(e,a);
		g3s_trtocc(e,e);
		if (!g3v_epseq(e,a)) {
			fprintf(stderr,"aaaa \n");
			g3v_print(a,stderr);
			g3v_print(e,stderr);
		}
	}
	fprintf(stderr,"ok\n");
	exit(1);
	g3v_set(a,atof(argv[1]),DEG2RAD(atof(argv[2])),DEG2RAD(atof(argv[3])));
	b = g3v_create();
	z = g3v_create();
	ang = g3v_create();
	g3v_set(z,0,0,1);
	g3v_normalize(a);
	for(i=0;i<5000;i++) {
		b[0] = 1.0;	
		b[1] = M_PI/2.0*rand()/RAND_MAX;
		b[2] = 2.0*M_PI*rand()/RAND_MAX;
		g3v_copy(ang,b);
		g3s_sphtocc(b,b);
		v = g3s_dist(a,b);
		if (!gb_epseq(sqrt(2.0 - 2.0*cos(v)),g3v_length(g3v_sub(ang,b,a))))
			fprintf(stderr,"oops %f %f\n",sqrt(2.0 - 2.0*cos(v)),g3v_length(g3v_sub(ang,b,a)));
		vo = (v > 0.2) ? 0.1 : (0.4 - v);
		vr = sqrt(2.0 - 2.0*cos(0.2));
		if (fabs(a[0] - b[0]) > vr || fabs(a[1] - b[1]) > vr || fabs(a[2] - b[2]) > vr) {
			vo -= 0.1;
			if (v < 0.2) {
				vo = 3;
				fprintf(stderr,"autsch %f %f\n",v,vr);
			}
		}
		g3s_cctosph(b,b);
		printf("%f %f %f\n",b[G3S_THETA],b[G3S_PHI],vo);
	}
	//printf("%f\n",g3s_dist_norm(a,b));
	//g3v_print(g3s_cctomtr(b,a),stdout);
	//g3v_print(g3s_cctosph(b,a),stdout);
	//printf("\n");
	return EXIT_SUCCESS;
}

#endif
