#ifndef lint
static const char	RCSid[] = "$Id: fprism.c,v 2.10 2014/07/08 18:25:00 greg Exp $";
#endif
/* Ce programme calcule les directions et les energies des rayons lumineux
   resultant du passage d'un rayon au travers d'un vitrage prismatique
 
   1991, LESO - EPFL, R. Compagnon - F. Di Pasquale  			     */

#include "standard.h"

#include "ray.h"
#include "calcomp.h"
#include "func.h"

#ifdef NOSTRUCTASSIGN
static double err = "No structure assignment!";	/* generate compiler error */
#endif


static double
Sqrt(
	double x
)
{
	if (x < 0.)
		return(0.);
	return(sqrt(x));
}

/* definitions de macros utiles */

#define ALPHA 0
#define BETA 1
#define GAMMA 2
#define DELTA 3
#define AUCUNE 4
#define X(r) r.v[0]
#define Y(r) r.v[1]
#define Z(r) r.v[2]
#define XX(v) v[0]
#define YY(v) v[1]
#define ZZ(v) v[2]
#define alpha_beta(v_alpha,v_beta) tfm(matbt,v_alpha,v_beta)
#define beta_alpha(v_beta,v_alpha) tfm(matb,v_beta,v_alpha)
#define alpha_gamma(v_alpha,v_gamma) tfm(matct,v_alpha,v_gamma)
#define gamma_alpha(v_gamma,v_alpha) tfm(matc,v_gamma,v_alpha)
#define prob_alpha_gamma(r) (1.-prob_alpha_beta(r))
#define prob_beta_gamma(r) (1.-prob_beta_alpha(r))
#define prob_gamma_beta(r) (1.-prob_gamma_alpha(r))
#define prob_delta_gamma(r) (1.-prob_delta_beta(r))
#define prob_beta_delta(r) (prob_beta_alpha(r))
#define prob_gamma_delta(r) (prob_gamma_alpha(r))
#define prob_delta_beta(r) (prob_alpha_beta(r))


/* Definitions des types de donnees */

typedef struct { FVECT v;		/* direction */
		 double ppar1,pper1,
			ppar2,pper2;	/* polarisations */
		 double e;		/* energie */
		 double n;		/* milieu dans lequel on se propage */
		 int orig,dest;		/* origine et destination          */
		} TRAYON;	

typedef struct { double a,b,c,d;	/* longueurs caracteristiques      */
		 double np;		/* indice de refraction            */
		} TPRISM;


/* Definitions des variables globales */

static TPRISM prism;		     	    /* le prisme ! */
static MAT4 matb = MAT4IDENT;		/* matrices de changement de bases */
static MAT4 matbt = MAT4IDENT;
static MAT4 matc = MAT4IDENT;
static MAT4 matct = MAT4IDENT;			
static double seuil;			/* seuil pour l'arret du trace     */
static double sinus,cosinus;		/* sin et cos    		   */
static double rapport;			/* rapport entre les indices des 
					   milieux refracteur et incident  */
static int tot_ref;			/* flag pour les surfaces 
					   reflechissantes		   */
static double fact_ref[4]={1.0,1.0,1.0,1.0};	/* facteurs de reflexion     */ 
static double tolerance;		/* degre de tol. pour les amalgames */
static double tolsource;		/* degre de tol. pour les sources  */
static int bidon;
#define BADVAL	(-10)
static long prismclock = -1;
static int nosource;			/* indique que l'on ne trace pas
					   en direction d'une source */
static int sens;			/* indique le sens de prop. du ray.*/  
static int nbrayons;			/* indice des rayons sortants      */
static TRAYON *ray;			/* tableau des rayons sortants     */
static TRAYON *raytemp;			/* variable temporaire    	   */


static void prepare_matrices(void);
static void tfm(MAT4 mat, FVECT v_old, FVECT v_new);
static double prob_alpha_beta(TRAYON r);
static double prob_beta_alpha(TRAYON r);
static double prob_gamma_alpha(TRAYON r);
static void v_par(FVECT v, FVECT v_out);
static void v_per(FVECT v, FVECT v_out);
static TRAYON transalphabeta(TRAYON r_initial);
static TRAYON transbetaalpha(TRAYON r_initial);
static TRAYON transalphagamma(TRAYON r_initial);
static TRAYON transgammaalpha(TRAYON r_initial);
static int compare(TRAYON r1, TRAYON r2, double marge);
static void sortie(TRAYON r);
static void trigo(TRAYON r);
static TRAYON reflexion(TRAYON r_incident);
static TRAYON transmission(TRAYON r_incident);
static void trace_rayon(TRAYON r_incident);
static void inverser(TRAYON *r1, TRAYON *r2);
static void setprism(void);
static double l_get_val(char *nm);


/* Definition des routines */

#define term(a,b) a/Sqrt(a*a+b*b)
static void
prepare_matrices(void)
{
	/* preparation des matrices de changement de bases */

	matb[0][0] = matbt[0][0] = matb[1][1] = matbt[1][1] = term(prism.a,prism.d);
	matb[1][0] = matbt[0][1] = term(-prism.d,prism.a);
	matb[0][1] = matbt[1][0] = term(prism.d,prism.a);
	matc[0][0] = matct[0][0] = matc[1][1] = matct[1][1] = term(prism.b,prism.d);
	matc[1][0] = matct[0][1] = term(prism.d,prism.b);
	matc[0][1] = matct[1][0] = term(-prism.d,prism.b);
	return;
}
#undef term


static void
tfm(
	MAT4 mat,
	FVECT v_old,
	FVECT v_new
)
{
	/* passage d'un repere old au repere new par la matrice mat */
	FVECT v_temp;

	multv3(v_temp,v_old,mat);
	normalize(v_temp);
	VCOPY(v_new,v_temp);
	return;
}

#define A prism.a
#define B prism.b
#define C prism.c
#define D prism.d


static double
prob_alpha_beta(
	TRAYON r
)
{
	/* calcul de la probabilite de passage de alpha a beta */
	double prob,test;

	if ( X(r) != 0. )
	{
		test = Y(r)/X(r);
		if ( test > B/D ) prob = 1.;
		else if ( test >= -A/D ) prob = (A+test*D)/(A+B);
		else prob = 0.;
	} 
	else prob = 0.;
	return prob;
}


static double
prob_beta_alpha(
	TRAYON r
)
{
	/* calcul de la probabilite de passage de beta a aplha */
	double prob,test;

	if ( X(r) != 0. )
	{
		test = Y(r)/X(r);
		if ( test > B/D ) prob = (A+B)/(A+test*D);
		else if ( test >= -A/D ) prob = 1.;
		else prob = 0.;
	}
	else prob = 0.;
	return prob;
}


static double
prob_gamma_alpha(
	TRAYON r
)
{
	/* calcul de la probabilite de passage de gamma a alpha */
	double prob,test;

	if ( X(r) != 0. )
	{
		test = Y(r)/X(r);
		if ( test > B/D ) prob = 0.;
		else if ( test >= -A/D ) prob = 1.;
		else prob = (A+B)/(B-test*D);
	} 
	else prob = 0.;
	return prob;
}

#undef A
#undef B
#undef C
#undef D


static void
v_par(
	FVECT v,
	FVECT v_out
)
/* calcule le vecteur par au plan d'incidence lie a v */
{
	FVECT v_temp;
	double det;

	det = Sqrt( (YY(v)*YY(v)+ZZ(v)*ZZ(v))*(YY(v)*YY(v)+ZZ(v)*ZZ(v))+
			(XX(v)*XX(v)*YY(v)*YY(v))+(XX(v)*XX(v)*ZZ(v)*ZZ(v)) );
	XX(v_temp) = (YY(v)*YY(v)+ZZ(v)*ZZ(v))/det;
	YY(v_temp) = -( XX(v)*YY(v) )/det;
	ZZ(v_temp) = -( XX(v)*ZZ(v) )/det;
	VCOPY(v_out,v_temp); 
	return;
}


static void
v_per(
	FVECT v,
	FVECT v_out
)
/* calcule le vecteur perp au plan d'incidence lie a v */
{
	FVECT v_temp;
	double det;

	det = Sqrt( (ZZ(v)*ZZ(v)+YY(v)*YY(v)) );
	XX(v_temp) = 0.;
	YY(v_temp) = -ZZ(v)/det;
	ZZ(v_temp) = YY(v)/det;
	VCOPY(v_out,v_temp);
	return;
}


static TRAYON
transalphabeta(
	TRAYON r_initial
)
/* transforme le rayon r_initial de la base associee a alpha dans
   la base associee a beta */ 
{
	TRAYON r_final;
	FVECT vpar_temp1,vpar_temp2,vper_temp1,vper_temp2;

	r_final = r_initial;
	alpha_beta(r_initial.v,r_final.v);
	if ((Y(r_initial) != 0. || Z(r_initial) != 0.)&&(Y(r_final) !=0. || Z(r_final)!= 0.))
	{
		v_par(r_initial.v,vpar_temp1);
		alpha_beta(vpar_temp1,vpar_temp1);
		v_per(r_initial.v,vper_temp1);
		alpha_beta(vper_temp1,vper_temp1);
		v_par(r_final.v,vpar_temp2);
		v_per(r_final.v,vper_temp2); 
		r_final.ppar1 = (r_initial.ppar1*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vpar_temp2));
		r_final.pper1 = (r_initial.ppar1*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vper_temp2));
		r_final.ppar2 = (r_initial.ppar2*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vpar_temp2));
		r_final.pper2 = (r_initial.ppar2*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vper_temp2));
	}
	return r_final;
}


static TRAYON
transbetaalpha(
	TRAYON r_initial
)
{
	/* transforme le rayon r_initial de la base associee a beta dans
	   la base associee a alpha */ 
	TRAYON r_final;
	FVECT vpar_temp1,vpar_temp2,vper_temp1,vper_temp2;

	r_final = r_initial;
	beta_alpha(r_initial.v,r_final.v);
	if ((Y(r_initial) != 0. || Z(r_initial) != 0. )&&(Y(r_final) != 0. || Z(r_final)!= 0.))
	{
		v_par(r_initial.v,vpar_temp1);
		beta_alpha(vpar_temp1,vpar_temp1);
		v_per(r_initial.v,vper_temp1);
		beta_alpha(vper_temp1,vper_temp1);
		v_par(r_final.v,vpar_temp2);
		v_per(r_final.v,vper_temp2);
		r_final.ppar1 = (r_initial.ppar1*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vpar_temp2));
		r_final.pper1 = (r_initial.ppar1*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vper_temp2));
		r_final.ppar2 = (r_initial.ppar2*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vpar_temp2));
		r_final.pper2 = (r_initial.ppar2*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vper_temp2));

	}
	return r_final;
}

 
static TRAYON
transalphagamma(
	TRAYON r_initial
)
/* transforme le rayon r_initial de la base associee a alpha dans
   la base associee a gamma */ 
{
	TRAYON r_final;
	FVECT vpar_temp1,vpar_temp2,vper_temp1,vper_temp2;

	r_final = r_initial;
	alpha_gamma(r_initial.v,r_final.v);
	if (( Y(r_initial) != 0. || Z(r_initial) != 0. )&&(Y(r_final)!= 0. || Z(r_final) !=0.))
	{
		v_par(r_initial.v,vpar_temp1);
		alpha_gamma(vpar_temp1,vpar_temp1);
		v_per(r_initial.v,vper_temp1);
		alpha_gamma(vper_temp1,vper_temp1);
		v_par(r_final.v,vpar_temp2);
		v_per(r_final.v,vper_temp2); 
		r_final.ppar1 = (r_initial.ppar1*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vpar_temp2));
		r_final.pper1 = (r_initial.ppar1*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vper_temp2));
		r_final.ppar2 = (r_initial.ppar2*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vpar_temp2));
		r_final.pper2 = (r_initial.ppar2*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vper_temp2));

	}
	return r_final;
}

 
static TRAYON
transgammaalpha(
	TRAYON r_initial
)
/* transforme le rayon r_initial de la base associee a gamma dans
   la base associee a alpha */ 
{
	TRAYON r_final;
	FVECT vpar_temp1,vpar_temp2,vper_temp1,vper_temp2;

	r_final = r_initial;
	gamma_alpha(r_initial.v,r_final.v);
	if (( Y(r_initial) != 0. || Z(r_initial) != 0. )&&(Y(r_final) !=0. || Z(r_final) != 0.))
	{
		v_par(r_initial.v,vpar_temp1);
		gamma_alpha(vpar_temp1,vpar_temp1);
		v_per(r_initial.v,vper_temp1);
		gamma_alpha(vper_temp1,vper_temp1);
		v_par(r_final.v,vpar_temp2);
		v_per(r_final.v,vper_temp2); 
		r_final.ppar1 = (r_initial.ppar1*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vpar_temp2));
		r_final.pper1 = (r_initial.ppar1*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper1*fdot(vper_temp1,vper_temp2));
		r_final.ppar2 = (r_initial.ppar2*fdot(vpar_temp1,vpar_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vpar_temp2));
		r_final.pper2 = (r_initial.ppar2*fdot(vpar_temp1,vper_temp2))+
			(r_initial.pper2*fdot(vper_temp1,vper_temp2));
	}
	return r_final;
}
 


static int
compare(
	TRAYON r1,
	TRAYON r2,
	double marge
)
{
	double arctg1, arctg2;

	arctg1 = atan2(Y(r1),X(r1));
	arctg2 = atan2(Y(r2),X(r2));
	if ((arctg1 - marge <= arctg2) && (arctg1 + marge >= arctg2)) return 1;
	else return 0;
}




static void
sortie(
	TRAYON r
)
{
	int i = 0;
	int egalite = 0;


	if(r.e > seuil)
	{
		while (i < nbrayons && egalite == 0)
		{
			raytemp = &ray[i];
			egalite = compare(r,*raytemp,tolerance);
			if (egalite) raytemp->e = raytemp->e + r.e;
			else i = i + 1;
		}
		if (egalite == 0)
		{
			if (nbrayons == 0) ray = (TRAYON *)calloc(1,sizeof(TRAYON));
			else ray = (TRAYON *)realloc((void *)ray, (nbrayons+1)*(sizeof(TRAYON)));         
			if (ray == NULL)
				error(SYSTEM, "out of memory in sortie\n");
			raytemp = &ray[nbrayons];
			raytemp->v[0] = X(r);
			raytemp->v[1] = Y(r);
			raytemp->v[2] = Z(r);
			raytemp->e = r.e;
			nbrayons++;
		}
	}
	return;
}


static void
trigo(
	TRAYON r
)
/* calcule les grandeurs trigonometriques relatives au rayon incident
   et le rapport entre les indices du milieu refracteur et incident   */
{
	double det;

	det = Sqrt(X(r)*X(r)+Y(r)*Y(r)+Z(r)*Z(r));
	sinus = Sqrt(Y(r)*Y(r)+Z(r)*Z(r))/det;
	cosinus = Sqrt(X(r)*X(r))/det;
	if (r.n == 1.) rapport = prism.np * prism.np;
	else rapport = 1./(prism.np * prism.np);
	return;
}


static TRAYON
reflexion(
	TRAYON r_incident
)
{
 /* calcul du rayon reflechi par une face */
 TRAYON r_reflechi;

 r_reflechi = r_incident; 
 trigo(r_incident);
 X(r_reflechi) = -X(r_incident);
 Y(r_reflechi) = Y(r_incident);
 Z(r_reflechi) = Z(r_incident);
 if(sinus > Sqrt(rapport) || r_incident.dest == tot_ref)
	{
	 r_reflechi.ppar1 = r_incident.ppar1;
	 r_reflechi.pper1 = r_incident.pper1;
         r_reflechi.ppar2 = r_incident.ppar2;
	 r_reflechi.pper2 = r_incident.pper2;
	 r_reflechi.e = r_incident.e * fact_ref[r_incident.dest];
	}
 else
	{ 
	 r_reflechi.ppar1 = r_incident.ppar1*(rapport*cosinus-Sqrt(rapport-
		 (sinus*sinus)))/(rapport*cosinus+Sqrt(rapport-(sinus*sinus)));
 	 r_reflechi.pper1 = r_incident.pper1*(cosinus-Sqrt
	        (rapport-(sinus*sinus)))/(cosinus+Sqrt(rapport-(sinus*sinus)));
 	 r_reflechi.ppar2 = r_incident.ppar2*(rapport*cosinus-Sqrt(rapport-
		 (sinus*sinus)))/(rapport*cosinus+Sqrt(rapport-(sinus*sinus)));
 	 r_reflechi.pper2 = r_incident.pper2*(cosinus-Sqrt
	        (rapport-(sinus*sinus)))/(cosinus+Sqrt(rapport-(sinus*sinus)));
	 r_reflechi.e = r_incident.e *(((r_reflechi.ppar1*r_reflechi.ppar1+
         r_reflechi.pper1*r_reflechi.pper1)/(r_incident.ppar1*r_incident.ppar1+
	 r_incident.pper1*r_incident.pper1))+((r_reflechi.ppar2*r_reflechi.ppar2
	 +r_reflechi.pper2*r_reflechi.pper2)/(r_incident.ppar2*r_incident.ppar2
	 +r_incident.pper2*r_incident.pper2)))/2;
	}
 
 /* a la sortie de cette routine r_transmis.orig et .dest ne sont pas definis!*/
 return r_reflechi;
}


static TRAYON
transmission(
	TRAYON r_incident
)
{
 /* calcul du rayon refracte par une face */
 TRAYON r_transmis;

 r_transmis = r_incident;
 trigo(r_incident); 
 if (sinus <= Sqrt(rapport) && r_incident.dest != tot_ref) 
 	{
 	 X(r_transmis) = (X(r_incident)/(fabs(X(r_incident))))*
			 (Sqrt(1.-(Y(r_incident)*Y(r_incident)+Z(r_incident)*
	         	           Z(r_incident))/rapport));
 	 Y(r_transmis) = Y(r_incident)/Sqrt(rapport);
	 Z(r_transmis) = Z(r_incident)/Sqrt(rapport);
	 r_transmis.ppar1 = r_incident.ppar1*2.*Sqrt(rapport)*cosinus/
		 	   (Sqrt(rapport-sinus*sinus)+rapport*cosinus);
	 r_transmis.pper1 = r_incident.pper1*2.*cosinus/(cosinus+Sqrt(rapport
			   - sinus*sinus));
	 r_transmis.ppar2 = r_incident.ppar2*2.*Sqrt(rapport)*cosinus/
		 	   (Sqrt(rapport-sinus*sinus)+rapport*cosinus);
	 r_transmis.pper2 = r_incident.pper2*2.*cosinus/(cosinus+Sqrt(rapport
			   - sinus*sinus));
	 r_transmis.e = (r_incident.e/2)*(Sqrt(rapport-sinus*sinus)/cosinus)
	    *(((r_transmis.ppar1*r_transmis.ppar1+r_transmis.pper1*
	        r_transmis.pper1)
	    /(r_incident.ppar1*r_incident.ppar1+r_incident.pper1*
		r_incident.pper1))+
	    ((r_transmis.ppar2*r_transmis.ppar2+r_transmis.pper2*r_transmis.pper2)
	   /(r_incident.ppar2*r_incident.ppar2+r_incident.pper2*r_incident.pper2)));
	 if(r_incident.n == 1.) r_transmis.n = prism.np;
	 else r_transmis.n = 1.;
	}
 else r_transmis.e = 0.;

 /* a la sortie de cette routine r_transmis.orig et .dest ne sont pas definis!*/

 return r_transmis;
}




#define ensuite(rayon,prob_passage,destination) r_suite = rayon; \
				 r_suite.e = prob_passage(rayon)*rayon.e; \
				 r_suite.dest = destination; \
				 if ( r_suite.e > seuil ) trace_rayon(r_suite)


static void
trace_rayon(
	TRAYON r_incident
)
{
 /* trace le rayon donne */
 TRAYON r_reflechi,r_transmis,r_suite;

 switch (r_incident.dest)
	{
	 case ALPHA:
		if ( r_incident.orig == ALPHA )
			{
			 r_reflechi = reflexion(r_incident);
			 sortie(r_reflechi);

			 r_transmis = transmission(r_incident);
			 r_transmis.orig = ALPHA;

			 ensuite(r_transmis,prob_alpha_beta,BETA);
			 ensuite(r_transmis,prob_alpha_gamma,GAMMA);
			} 
		else
			{
			 r_reflechi = reflexion(r_incident);
			 r_reflechi.orig = ALPHA;
			 ensuite(r_reflechi,prob_alpha_beta,BETA);
			 ensuite(r_reflechi,prob_alpha_gamma,GAMMA);

			 r_transmis = transmission(r_incident); 
			 sortie(r_transmis);
			}
		break;
	 case BETA:
	   r_reflechi = transbetaalpha(reflexion(transalphabeta(r_incident)));
	   r_reflechi.orig = BETA;
	   r_transmis = transbetaalpha(transmission(transalphabeta
			(r_incident)));
	   r_transmis.orig = GAMMA;	
		if ( r_incident.n > 1.0 )  /* le rayon vient de l'interieur */
			{
			 ensuite(r_reflechi,prob_beta_alpha,ALPHA);
			 ensuite(r_reflechi,prob_beta_gamma,GAMMA); 

			 ensuite(r_transmis,prob_beta_gamma,GAMMA);
			 ensuite(r_transmis,prob_beta_delta,DELTA);
			}
		 else	/* le rayon vient de l'exterieur */
			{
			 ensuite(r_reflechi,prob_beta_gamma,GAMMA);
			 ensuite(r_reflechi,prob_beta_delta,DELTA);

			 ensuite(r_transmis,prob_beta_alpha,ALPHA);
			 ensuite(r_transmis,prob_beta_gamma,GAMMA);
			}
		 break;
	 case GAMMA:
	   r_reflechi = transgammaalpha(reflexion(transalphagamma(r_incident)));
	   r_reflechi.orig = GAMMA; 
	   r_transmis = transgammaalpha(transmission(transalphagamma
		        (r_incident)));		 
	   r_transmis.orig = GAMMA;	
		if ( r_incident.n > 1.0 )  /* le rayon vient de l'interieur */
			{
			 ensuite(r_reflechi,prob_gamma_alpha,ALPHA);
			 ensuite(r_reflechi,prob_gamma_beta,BETA);

			 ensuite(r_transmis,prob_gamma_beta,BETA);
			 ensuite(r_transmis,prob_gamma_delta,DELTA);
			}
		 else	/* le rayon vient de l'exterieur */
			{
			 ensuite(r_reflechi,prob_gamma_beta,BETA);
			 ensuite(r_reflechi,prob_gamma_delta,DELTA);

			 ensuite(r_transmis,prob_gamma_alpha,ALPHA);
			 ensuite(r_transmis,prob_gamma_beta,BETA);
			}
		 break;
	 case DELTA:
		 if ( r_incident.orig != DELTA ) sortie(r_incident);
		 else
			{
			 ensuite(r_incident,prob_delta_beta,BETA);
			 ensuite(r_incident,prob_delta_gamma,GAMMA);
			}
		 break;
	}
 return;
}

#undef ensuite

static void
inverser(
	TRAYON *r1,
	TRAYON *r2
)
{
	TRAYON temp;
	temp = *r1;
	*r1 = *r2;
	*r2 = temp;
}



static void
setprism(void)
{
 double d;
 TRAYON r_initial,rsource;
 int i,j; 

 prismclock = eclock;
r_initial.ppar1 = r_initial.pper2 = 1.;
r_initial.pper1 = r_initial.ppar2 = 0.;

d = 1; prism.a = funvalue("arg", 1, &d);
 if(prism.a < 0.) goto badopt;
d = 2; prism.b = funvalue("arg", 1, &d);
 if(prism.b < 0.) goto badopt;
d = 3; prism.c = funvalue("arg", 1, &d);
 if(prism.c < 0.) goto badopt;
d = 4; prism.d = funvalue("arg", 1, &d);
 if(prism.d < 0.) goto badopt;
d = 5; prism.np = funvalue("arg", 1, &d);
 if(prism.np < 1.) goto badopt;
d = 6; seuil = funvalue("arg", 1, &d);
 if (seuil < 0. || seuil >=1) goto badopt;
d = 7; tot_ref = (int)(funvalue("arg", 1, &d) + .5);
 if (tot_ref != 1 && tot_ref != 2 && tot_ref != 4) goto badopt;
 if (tot_ref < 4 ) 
	{
	 d = 8; fact_ref[tot_ref] = funvalue("arg", 1, &d);
 	 if (fact_ref[tot_ref] < 0. || fact_ref[tot_ref] > 1.) goto badopt;
	}
d = 9; tolerance = funvalue("arg", 1, &d);
 if (tolerance <= 0.) goto badopt;
d = 10; tolsource = funvalue("arg", 1, &d);
 if (tolsource < 0. ) goto badopt;
X(r_initial) = varvalue("Dx");
Y(r_initial) = varvalue("Dy");
Z(r_initial) = varvalue("Dz");
#ifdef DEBUG
fprintf(stderr,"dx=%lf dy=%lf dz=%lf\n",X(r_initial),Y(r_initial),Z(r_initial));
#endif

 /* initialisation */ 
 prepare_matrices();
 r_initial.e = 1.0;
 r_initial.n = 1.0;

if(ray!=NULL) free(ray);	
nbrayons = 0;
/* determination de l'origine et de la destination du rayon initial */

if ( X(r_initial) != 0.)
 {
  if ( X(r_initial) > 0. )
	{
	 r_initial.orig = r_initial.dest = ALPHA;
	 sens = 1;
	}
  else if ( X(r_initial) < 0. ) 
	{
	 r_initial.orig = r_initial.dest = DELTA;
	 sens = -1;
	}
 
  normalize(r_initial.v); 

  trace_rayon(r_initial);

  X(rsource) = varvalue("DxA");
  Y(rsource) = varvalue("DyA");
  Z(rsource) = varvalue("DzA");
  nosource = ( X(rsource)==0. && Y(rsource)==0. && Z(rsource)==0. );
  if ( !nosource )
	{
	 for (j=0; j<nbrayons; j++)
		{
		 if ( !compare(ray[j],rsource,tolsource) ) ray[j].e =0.;
		}
	}
  for  (j = 0; j < nbrayons; j++)
        {
         for (i = j+1; i < nbrayons; i++)
	  {
	   if (ray[j].e < ray[i].e) inverser(&ray[j],&ray[i]);
	  }
	}

    bidon = 1;
  }
 else bidon = 0;
 return;

 /* message puis sortie si erreur dans la ligne de commande */
 badopt:
 bidon = BADVAL;
 return;
}


static double
l_get_val(
	char *nm
)
{
 int val, dir, i, trouve, curseur;
 int nb;
 double valeur;
 TRAYON *rayt=NULL, raynull;
 
 if (prismclock < 0 || prismclock < eclock) setprism();
 if (bidon == BADVAL) {
	errno = EDOM;
	return(0.0);
 }
 val = (int)(argument(1) + .5);
 dir = (int)(argument(2) + .5);
 nb = (int)(argument(3) + .5);
 X(raynull) = bidon;
 Y(raynull) = Z(raynull) = 0.;
 raynull.e = 0.;
 trouve = curseur = 0;
 if ( !nosource && nb==2 ) nb=1; /* on est en train de tracer la source
				     a partir de sa seconde source virtuelle */ 
#ifdef DEBUG
 fprintf(stderr, " On considere le rayon no: %d\n", nb);
#endif
 for(i=0; i < nbrayons &&!trouve; i++)
  {
   if(ray[i].v[0] * dir * sens >= 0.) curseur ++;
   if(curseur == nb)
   {
    rayt = &ray[i];
    trouve = 1;
   }
  }
 if(!trouve) rayt = &raynull; 
 switch(val) {
	case 0 : valeur = rayt->v[0];
	         break;
	case 1 : valeur = rayt->v[1];
		 break;
	case 2 : valeur = rayt->v[2];
		 break;
	case 3 : valeur = rayt->e;
		 break;
	default : errno = EDOM; return(0.0);
    } 
#ifdef DEBUG
  fprintf(stderr, "get_val( %i, %i, %i) = %lf\n",val,dir,nb,valeur);
#endif
  return valeur;
}


void
setprismfuncs(void)  /* declared in func.h */
{
 funset("fprism_val", 3, '=', l_get_val);
}
