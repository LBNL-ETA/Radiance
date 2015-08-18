#ifndef lint
static const char RCSid[] = "$Id: g3vector.c,v 2.2 2015/08/18 15:02:53 greg Exp $";
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "g3vector.h"



g3Vec	g3v_create()
{
	return ((g3Vec) calloc(3,sizeof(g3Float)));
}

void	g3v_free(g3Vec v)
{
	if (v) {
		free(v);
		v = NULL;
	}
}

g3Vec	g3v_copy(g3Vec dest,g3Vec src)
{
	if (dest == src)
		return dest;
	memcpy(dest,src,(3*sizeof(g3Float)));
	return dest;
}

g3Vec	g3v_zero(g3Vec v)
{
	memset(v,0,3*sizeof(g3Float));
	return v;
}

g3Vec	g3v_set(g3Vec v,g3Float x,g3Float y,g3Float z)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
	return v;
}

g3Vec	g3v_sane(g3Vec v)
{
	int i;
	for(i=0;i<3;i++) 
		if (gb_epseq(v[i],0.0))
			v[i] = 0.0;
	return v;
}

g3Vec	g3v_fromrad(g3Vec v,FVECT rv)
{
	int i;
	
	for(i=0;i<3;i++) 
		v[i] = rv[i];
	return v;
}

void	g3v_torad(FVECT rv,g3Vec v)
{
	int i;
	
	for(i=0;i<3;i++) 
		rv[i] = v[i];
}

int		g3v_epszero(g3Vec v)
{
	int i;

	for(i=0;i<3;i++) 
		if (!gb_epseq(v[i],0))
			return 0;
	return 1;
}

int		g3v_epsorder(g3Vec v1,g3Vec v2)
{
	int i;
	int res;

	for(i=0;i<3;i++) 
		if ((res = gb_epsorder(v1[i],v2[i]))) 
			return res;
	return 0;
}

int		g3v_epseq(g3Vec v1,g3Vec v2)
{
	return (g3v_epsorder(v1,v2) == 0);
}

g3Float	g3v_dot(g3Vec v1,g3Vec v2)
{
	int i;
	g3Float res = 0.0;

	for(i=0;i<3;i++) 
		res += v1[i]*v2[i];
	return res;
}


g3Float	g3v_length(g3Vec v)
{
	return sqrt(g3v_dot(v,v));
}

g3Vec	g3v_normalize(g3Vec v)
{
	int i;
	g3Float len = g3v_length(v);

	if (gb_epseq(len,0))
		return NULL;
	for(i=0;i<3;i++) 
		v[i] /= len;
	return v;
}

g3Vec	g3v_cross(g3Vec res,g3Vec v1,g3Vec v2)
{
	g3Vec h = NULL;
	if ((res == v1) || (res == v2))
		h = g3v_create();
	else
		h = res;
	h[0] = v1[1]*v2[2] - v1[2]*v2[1];
    h[1] = v1[2]*v2[0] - v1[0]*v2[2];
    h[2] = v1[0]*v2[1] - v1[1]*v2[0];
	if (h != res) {
		g3v_copy(res,h);
		g3v_free(h);
	}
	return res;
}

g3Vec	g3v_sub(g3Vec res,g3Vec v1,g3Vec v2)
{
	int i;

	for(i=0;i<3;i++) 
		res[i] = v1[i] - v2[i];
	return res;
}

g3Vec	g3v_add(g3Vec res,g3Vec v1,g3Vec v2)
{
	int i;

	for(i=0;i<3;i++) 
		res[i] = v1[i] + v2[i];
	return res;
}

g3Vec	g3v_scale(g3Vec res,g3Vec v,g3Float sc)
{
	int i;

	for(i=0;i<3;i++) 
		res[i] = v[i] * sc;
	return res;
}

g3Vec	g3v_rotate(g3Vec res,g3Vec v,g3Vec rnorm,g3Float ang)
{
	g3Vec c1;
	int copy = 0;

	if (res == v) {
		res = g3v_create();
		copy = 1;
	}
	ang = -ang;
	c1 = g3v_create();
	g3v_scale(res,v,cos(ang));
	g3v_add(res,res,g3v_scale(c1,rnorm,g3v_dot(rnorm,v)*(1 - cos(ang))));
	g3v_add(res,res,g3v_scale(c1,g3v_cross(c1,v,rnorm),sin(ang)));
	g3v_free(c1);
	if (copy) {
		g3v_copy(v,res);
		g3v_free(res);
		res = v;
	}
	return res;
}

g3Vec	g3v_print(g3Vec v,FILE* fp)
{
	int i;
	
	for(i=0;i<3;i++) 
		fprintf(fp,"%f ",v[i]);
	return v;
}

int		g3v_read(g3Vec v,FILE* fp)
{
	char st[100];

	if (!fgets(st,99,fp))
		return 0;
#if USE_DOUBLE
		if (sscanf(st,"%lf %lf %lf",&(v[0]),&(v[1]),&(v[2])) != 3)
			return 0;
#else
		if (sscanf(st,"%f %f %f",&(v[0]),&(v[1]),&(v[2])) != 3)
			return 0;
#endif
	return 1;
}

#ifdef VEC_TEST

int 	main(int argc,char** argv)
{
	g3Vec a,b;
	char st[100];

	a = g3v_create();
	b = g3v_create();
	g3v_set(b,atof(argv[1]),atof(argv[2]),atof(argv[3]));
	while (fgets(st,99,stdin)) {
		sscanf(st,"%lf %lf %lf",a,(a+1),(a+2));
		g3v_rotate(a,a,b,DEG2RAD(atof(argv[4])));
		g3v_print(a,stdout);
	}
	return EXIT_SUCCESS;
}

#endif
