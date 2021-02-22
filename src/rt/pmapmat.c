#ifndef lint
static const char RCSid[] = "$Id: pmapmat.c,v 2.24 2021/02/22 13:27:49 rschregle Exp $";
#endif
/* 

   ======================================================================
   Photon map support routines for scattering by materials. 

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme FARESYS")
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ======================================================================
   
*/



#include "pmapmat.h"
#include "pmapdata.h"
#include "pmaprand.h"
#include "otypes.h"
#include "data.h"
#include "func.h"
#include "bsdf.h"
#include <math.h>



/* Stuff ripped off from material modules */
#define  MAXITER     10
#define  SP_REFL     01
#define  SP_TRAN     02 
#define  SP_PURE     04
#define  SP_FLAT     010
#define  SP_BADU     040
#define  MLAMBDA     500
#define  RINDEX      1.52
#define  FRESNE(ci)  (exp(-5.85*(ci)) - 0.00287989916)



typedef struct {
   OBJREC   *mp;
   RAY      *rp;
   short    specfl;
   COLOR    mcolor, scolor;
   FVECT    vrefl, prdir, pnorm;
   double   alpha2, rdiff, rspec, trans, tdiff, tspec, pdot;
} NORMDAT;

typedef struct {
   OBJREC   *mp;
   RAY      *rp;
   short    specfl;
   COLOR    mcolor, scolor;
   FVECT    vrefl, prdir, u, v, pnorm;
   double   u_alpha, v_alpha, rdiff, rspec, trans, tdiff, tspec, pdot;
} ANISODAT;

typedef struct {
   OBJREC   *mp;
   RAY      *pr;
   DATARRAY *dp;
   COLOR    mcolor;
   COLOR    rdiff;
   COLOR    tdiff;
   double   rspec;
   double   trans;
   double   tspec;
   FVECT    pnorm;
   double   pdot;
} BRDFDAT;

typedef struct {
   OBJREC   *mp;
   RAY      *pr;
   FVECT    pnorm;
   FVECT    vray;
   double   sr_vpsa [2];
   RREAL    toloc [3][3];
   RREAL    fromloc [3][3];
   double   thick;
   SDData   *sd;
   COLOR	   runsamp;
   COLOR	   rdiff;
   COLOR	   tunsamp;
   COLOR	   tdiff;
} BSDFDAT;



extern const SDCDst SDemptyCD;

/* Per-material scattering function dispatch table; return value is usually
 * zero, indicating photon termination */
int (*photonScatter [NUMOTYPE]) (OBJREC*, RAY*);

/* List of antimatter sensor modifier names and associated object set */
char *photonSensorList [MAXSET + 1] = {NULL};
static OBJECT photonSensorSet [MAXSET + 1] = {0};



/* ================ General support routines ================ */


void photonRay (const RAY *rayIn, RAY *rayOut, 
                int rayOutType, COLOR fluxAtten)
/* Spawn a new photon ray from a previous one; this is effectively a
 * customised rayorigin(). 
 * A SPECULAR rayOutType flags this photon as _caustic_ for subsequent hits.
 * It is preserved for transferred rays (of type PMAP_XFER).
 * fluxAtten specifies the RGB attenuation of the photon flux effected by
 * the scattering material. The outgoing flux is then normalised to maintain
 * a uniform average of 1 over RGB. If fluxAtten == NULL, the flux remains
 * unchanged for the outgoing photon. fluxAtten is ignored for transferred
 * rays.
 * The ray direction is preserved for transferred rays, and undefined for
 * scattered rays and must be subsequently set by the caller. */
{
   rayorigin(rayOut, rayOutType, rayIn, NULL);
   
   if (rayIn) {
      /* Transfer flux */
      copycolor(rayOut -> rcol, rayIn -> rcol);
      
      /* Copy caustic flag & direction for transferred rays */
      if (rayOutType == PMAP_XFER) {
         /* rayOut -> rtype |= rayIn -> rtype & SPECULAR; */
         rayOut -> rtype |= rayIn -> rtype;
         VCOPY(rayOut -> rdir, rayIn -> rdir);
      }
      else if (fluxAtten) {
         /* Attenuate and normalise flux for scattered rays */
         multcolor(rayOut -> rcol, fluxAtten);
         colorNorm(rayOut -> rcol);
      }

      /* Propagate index of emitting light source */
      rayOut -> rsrc = rayIn -> rsrc;
      
      /* Update maximum photon path distance */
      rayOut -> rmax = rayIn -> rmax - rayIn -> rot;
   }
}


static void addPhotons (const RAY *r)
/* Insert photon hits, where applicable */
{
   if (!r -> rlvl)
      /* Add direct photon at primary hitpoint */
      newPhoton(directPmap, r);
   else {
      /* Add global or precomputed photon at indirect hitpoint */
      newPhoton(preCompPmap ? preCompPmap : globalPmap, r);

      /* Store caustic photon if specular flag set */
      if (PMAP_CAUSTICRAY(r))
         newPhoton(causticPmap, r);
         
      /* Store in contribution photon map */
      newPhoton(contribPmap, r);
   }
}



void getPhotonSensors (char **sensorList)
/* Find antimatter geometry declared as photon sensors */
{
   OBJECT i;
   OBJREC *obj;
   char **lp;
   
   /* Init sensor set */
   photonSensorSet [0] = 0;
   
   if (!sensorList [0])
      return;
      
   for (i = 0; i < nobjects; i++) {
      obj = objptr(i);
      
      /* Insert object in sensor set if it's in the specified sensor list
       * and of type antimatter */
      for (lp = sensorList; *lp; lp++) {
         if (!strcmp(obj -> oname, *lp)) {
            if (obj -> otype != MAT_CLIP) {
               sprintf(errmsg, "photon sensor modifier %s is not antimatter",
                       obj -> oname);
               error(USER, errmsg);
            }
            
            if (photonSensorSet [0] >= AMBLLEN)
               error(USER, "too many photon sensor modifiers");
            
            insertelem(photonSensorSet, i);
         }
      }
   }
   
   if (!photonSensorSet [0])
      error(USER, "no photon sensors found");
}



/* ================ Material specific scattering routines ================ */   


static int isoSpecPhotonScatter (NORMDAT *nd, RAY *rayOut)
/* Generate direction for isotropically specularly reflected 
   or transmitted ray. Returns 1 if successful. */
{
   FVECT    u, v, h;
   RAY      *rayIn = nd -> rp;
   double   d, d2, sinp, cosp;
   int      niter, i = 0;
   
   /* Set up sample coordinates */  
   getperpendicular(u, nd -> pnorm, 1);
   fcross(v, nd -> pnorm, u);
   
   if (nd -> specfl & SP_REFL) {
      /* Specular reflection; make MAXITER attempts at getting a ray */
      
      for (niter = 0; niter < MAXITER; niter++) {
         d = 2 * PI * pmapRandom(scatterState);
         cosp = cos(d);
         sinp = sin(d);
         d2 = pmapRandom(scatterState);
         d = d2 <= FTINY ? 1 : sqrt(nd -> alpha2 * -log(d2));
         
         for (i = 0; i < 3; i++)
            h [i] = nd -> pnorm [i] + d * (cosp * u [i] + sinp * v [i]);
            
         d = -2 * DOT(h, rayIn -> rdir) / (1 + d * d);
         VSUM(rayOut -> rdir, rayIn -> rdir, h, d);
         
         if (DOT(rayOut -> rdir, rayIn -> ron) > FTINY) 
            return 1;
      }
      
      return 0;
   }
   
   else {
      /* Specular transmission; make MAXITER attempts at getting a ray */
      
      for (niter = 0; niter < MAXITER; niter++) {
         d = 2 * PI * pmapRandom(scatterState);
         cosp = cos(d);
         sinp = sin(d);
         d2 = pmapRandom(scatterState);
         d = d2 <= FTINY ? 1 : sqrt(-log(d2) * nd -> alpha2);
         
         for (i = 0; i < 3; i++)
            rayOut -> rdir [i] = nd -> prdir [i] + 
                                 d * (cosp * u [i] + sinp * v [i]);
                                 
         if (DOT(rayOut -> rdir, rayIn -> ron) < -FTINY) {
            normalize(rayOut -> rdir);
            return 1;
         }
      }
      
      return 0;
   }
}



static void diffPhotonScatter (FVECT normal, RAY* rayOut)
/* Generate cosine-weighted direction for diffuse ray */
{
   const RREAL cosThetaSqr = pmapRandom(scatterState),
               cosTheta = sqrt(cosThetaSqr),
               sinTheta = sqrt(1 - cosThetaSqr),
               phi = 2 * PI * pmapRandom(scatterState),
               du = cos(phi) * sinTheta, dv = sin(phi) * sinTheta;
   FVECT       u, v;
   int         i = 0;

   /* Set up sample coordinates */
   getperpendicular(u, normal, 1);
   fcross(v, normal, u);
   
   /* Convert theta & phi to cartesian */
   for (i = 0; i < 3; i++)
      rayOut -> rdir [i] = du * u [i] + dv * v [i] + cosTheta * normal [i];
      
   normalize(rayOut -> rdir);
}



static int normalPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for isotropic material and recurse */
{
   NORMDAT  nd;
   int      i, hastexture;
   float    xi, albedo, prdiff, ptdiff, prspec, ptspec;
   double   d, fresnel;
   RAY      rayOut;

   if (mat -> oargs.nfargs != (mat -> otype == MAT_TRANS ? 7 : 5))
      objerror(mat, USER, "bad number of arguments");
      
   /* Check for back side; reorient if back is visible */
   if (rayIn -> rod < 0)
      if (!backvis && mat -> otype != MAT_TRANS) 
         return 0;
      else {
         /* Get modifiers */
         raytexture(rayIn, mat -> omod);
         flipsurface(rayIn);
      }
   else raytexture(rayIn, mat -> omod);
   
   nd.mp = mat;
   nd.rp = rayIn;
   
   /* Get material color */
   copycolor(nd.mcolor, mat -> oargs.farg);
   
   /* Get roughness */
   nd.specfl = 0;
   nd.alpha2 = mat -> oargs.farg [4];
   
   if ((nd.alpha2 *= nd.alpha2) <= FTINY) 
      nd.specfl |= SP_PURE;
      
   if (rayIn -> ro != NULL && isflat(rayIn -> ro -> otype)) 
      nd.specfl |= SP_FLAT;  
      
   /* Perturb normal */
   if ((hastexture = (DOT(rayIn -> pert, rayIn -> pert) > sqr(FTINY)) ))
      nd.pdot = raynormal(nd.pnorm, rayIn);
   else {
      VCOPY(nd.pnorm, rayIn -> ron);
      nd.pdot = rayIn -> rod;
   }
   
   nd.pdot = max(nd.pdot, .001);
      
   /* Modify material color */
   multcolor(nd.mcolor, rayIn -> pcol);
   nd.rspec = mat -> oargs.farg [3];
   
   /* Approximate Fresnel term */
   if (nd.specfl & SP_PURE && nd.rspec > FTINY) {
      fresnel = FRESNE(rayIn -> rod);
      nd.rspec += fresnel * (1 - nd.rspec);
   } 
   else fresnel = 0;
   
   /* Transmission params */
   if (mat -> otype == MAT_TRANS) {
      nd.trans = mat -> oargs.farg [5] * (1 - nd.rspec);
      nd.tspec = nd.trans * mat -> oargs.farg [6];
      nd.tdiff = nd.trans - nd.tspec; 
   } 
   else nd.tdiff = nd.tspec = nd.trans = 0;
   
   /* Specular reflection params */
   if (nd.rspec > FTINY) {      
      /* Specular color */
      if (mat -> otype != MAT_METAL) 
         setcolor(nd.scolor, nd.rspec, nd.rspec, nd.rspec);
      else if (fresnel > FTINY) {
         d = nd.rspec * (1 - fresnel);
         for (i = 0; i < 3; i++) 
            nd.scolor [i] = fresnel + nd.mcolor [i] * d;
      }
      else {
         copycolor(nd.scolor, nd.mcolor);
         scalecolor(nd.scolor, nd.rspec);
      }
   }
   else setcolor(nd.scolor, 0, 0, 0);
   
   /* Diffuse reflection params */
   nd.rdiff = 1 - nd.trans - nd.rspec;
   
   /* Set up probabilities */
   prdiff = ptdiff = ptspec = colorAvg(nd.mcolor);
   prdiff *= nd.rdiff;
   ptdiff *= nd.tdiff;
   prspec = colorAvg(nd.scolor);
   ptspec *= nd.tspec;
   albedo = prdiff + ptdiff + prspec + ptspec;
   
   /* Insert direct and indirect photon hits if diffuse component */
   if (prdiff > FTINY || ptdiff > FTINY)
      addPhotons(rayIn);
   
   xi = pmapRandom(rouletteState);

   if (xi > albedo) 
      /* Absorbed */
      return 0;
      
   if (xi > (albedo -= prspec)) {
      /* Specular reflection */
      nd.specfl |= SP_REFL;
      
      if (nd.specfl & SP_PURE) {
         /* Perfect specular reflection */
         for (i = 0; i < 3; i++) {
            /* Reflected ray */
            nd.vrefl [i] = rayIn -> rdir [i] + 2 * nd.pdot * nd.pnorm [i];
         }
         
         /* Penetration? */
         if (hastexture && DOT(nd.vrefl, rayIn -> ron) <= FTINY)
            for (i = 0; i < 3; i++) {
               /* Safety measure */
               nd.vrefl [i] = rayIn -> rdir [i] + 
                              2 * rayIn -> rod * rayIn -> ron [i];
            }
            
         VCOPY(rayOut.rdir, nd.vrefl);
      }
      
      else if (!isoSpecPhotonScatter(&nd, &rayOut)) 
         return 0;
      
      photonRay(rayIn, &rayOut, PMAP_SPECREFL, nd.scolor);
   }
   
   else if (xi > (albedo -= ptspec)) {
      /* Specular transmission */
      nd.specfl |= SP_TRAN;
      
      if (hastexture) {
         /* Perturb */
         for (i = 0; i < 3; i++)
            nd.prdir [i] = rayIn -> rdir [i] - rayIn -> pert [i];
            
         if (DOT(nd.prdir, rayIn -> ron) < -FTINY)
            normalize(nd.prdir);
         else VCOPY(nd.prdir, rayIn -> rdir);
      }
      else VCOPY(nd.prdir, rayIn -> rdir);
      
      if ((nd.specfl & (SP_TRAN | SP_PURE)) == (SP_TRAN | SP_PURE))
         /* Perfect specular transmission */
         VCOPY(rayOut.rdir, nd.prdir);
      else if (!isoSpecPhotonScatter(&nd, &rayOut))
         return 0;
         
      photonRay(rayIn, &rayOut, PMAP_SPECTRANS, nd.mcolor);
   }
   
   else if (xi > (albedo -= prdiff)) {
      /* Diffuse reflection */
      photonRay(rayIn, &rayOut, PMAP_DIFFREFL, nd.mcolor);
      diffPhotonScatter(hastexture ? nd.pnorm : rayIn -> ron, &rayOut);
   }
   
   else {
      /* Diffuse transmission */
      flipsurface(rayIn);
      photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, nd.mcolor);

      if (hastexture) {
         FVECT bnorm;
         bnorm [0] = -nd.pnorm [0];
         bnorm [1] = -nd.pnorm [1];
         bnorm [2] = -nd.pnorm [2];
         diffPhotonScatter(bnorm, &rayOut);
      } 
      else diffPhotonScatter(rayIn -> ron, &rayOut);
   }
   
   tracePhoton(&rayOut);
   return 0;
}



static void getacoords (ANISODAT *nd)
/* Set up coordinate system for anisotropic sampling; cloned from aniso.c */
{
   MFUNC *mf;
   int   i;

   mf = getfunc(nd -> mp, 3, 0x7, 1);
   setfunc(nd -> mp, nd -> rp);
   errno = 0;

   for (i = 0; i < 3; i++)
      nd -> u [i] = evalue(mf -> ep [i]);
   
   if (errno == EDOM || errno == ERANGE)
      nd -> u [0] = nd -> u [1] = nd -> u [2] = 0.0;
      
   if (mf -> fxp != &unitxf)
      multv3(nd -> u, nd -> u, mf -> fxp -> xfm);

   fcross(nd -> v, nd -> pnorm, nd -> u);
   
   if (normalize(nd -> v) == 0.0) {
      if (fabs(nd -> u_alpha - nd -> v_alpha) > 0.001)
         objerror(nd -> mp, WARNING, "illegal orientation vector");
      getperpendicular(nd -> u, nd -> pnorm, 1);
      fcross(nd -> v, nd -> pnorm, nd -> u);
      nd -> u_alpha = nd -> v_alpha = 
         sqrt(0.5 * (sqr(nd -> u_alpha) + sqr(nd -> v_alpha)));
   } 
   else fcross(nd -> u, nd -> v, nd -> pnorm);
}



static int anisoSpecPhotonScatter (ANISODAT *nd, RAY *rayOut)
/* Generate direction for anisotropically specularly reflected 
   or transmitted ray. Returns 1 if successful. */
{
   FVECT    h;
   double   d, d2, sinp, cosp;
   int      niter, i;
   RAY      *rayIn = nd -> rp;

   if (rayIn -> ro != NULL && isflat(rayIn -> ro -> otype)) 
      nd -> specfl |= SP_FLAT;
      
   /* set up coordinates */
   getacoords(nd);
   
   if (rayOut -> rtype & TRANS) {
      /* Specular transmission */

      if (DOT(rayIn -> pert, rayIn -> pert) <= sqr(FTINY)) 
         VCOPY(nd -> prdir, rayIn -> rdir);
      else {
         /* perturb */
         for (i = 0; i < 3; i++) 
            nd -> prdir [i] = rayIn -> rdir [i] - rayIn -> pert [i];
            
         if (DOT(nd -> prdir, rayIn -> ron) < -FTINY) 
            normalize(nd -> prdir);
         else VCOPY(nd -> prdir, rayIn -> rdir);
      }
      
      /* Make MAXITER attempts at getting a ray */
      for (niter = 0; niter < MAXITER; niter++) {
         d = 2 * PI * pmapRandom(scatterState);
         cosp = cos(d) * nd -> u_alpha;
         sinp = sin(d) * nd -> v_alpha;
         d = sqrt(sqr(cosp) + sqr(sinp));
         cosp /= d;
         sinp /= d;
         d2 = pmapRandom(scatterState);
         d = d2 <= FTINY ? 1
                         : sqrt(-log(d2) /
                                (sqr(cosp) / sqr(nd -> u_alpha) + 
                                 sqr(sinp) / (nd -> v_alpha * nd -> u_alpha)));
                                 
         for (i = 0; i < 3; i++)
            rayOut -> rdir [i] = nd -> prdir [i] + d * 
                                 (cosp * nd -> u [i] + sinp * nd -> v [i]);
                                 
         if (DOT(rayOut -> rdir, rayIn -> ron) < -FTINY) {
            normalize(rayOut -> rdir);
            return 1;
         }
      }
      
      return 0;
   }
   
   else {
      /* Specular reflection */
      
      /* Make MAXITER attempts at getting a ray */
      for (niter = 0; niter < MAXITER; niter++) {
         d = 2 * PI * pmapRandom(scatterState);
         cosp = cos(d) * nd -> u_alpha;
         sinp = sin(d) * nd -> v_alpha;
         d = sqrt(sqr(cosp) + sqr(sinp));
         cosp /= d;
         sinp /= d;
         d2 = pmapRandom(scatterState);
         d = d2 <= FTINY ? 1 
                         : sqrt(-log(d2) / 
                                (sqr(cosp) / sqr(nd -> u_alpha) +
                                 sqr(sinp) / (nd->v_alpha * nd->v_alpha)));
                                 
         for (i = 0; i < 3; i++)
            h [i] = nd -> pnorm [i] + 
                    d * (cosp * nd -> u [i] + sinp * nd -> v [i]);
                    
         d = -2 * DOT(h, rayIn -> rdir) / (1 + d * d);
         VSUM(rayOut -> rdir, rayIn -> rdir, h, d);
         
         if (DOT(rayOut -> rdir, rayIn -> ron) > FTINY) 
            return 1;
      }
      
      return 0;
   }
}



static int anisoPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for anisotropic material and recurse */
{
   ANISODAT nd;
   float    xi, albedo, prdiff, ptdiff, prspec, ptspec;
   RAY      rayOut;
   
   if (mat -> oargs.nfargs != (mat -> otype == MAT_TRANS2 ? 8 : 6))
      objerror(mat, USER, "bad number of real arguments");
      
   nd.mp = mat;
   nd.rp = rayIn;
   
   /* get material color */
   copycolor(nd.mcolor, mat -> oargs.farg);
   
   /* get roughness */
   nd.specfl = 0;
   nd.u_alpha = mat -> oargs.farg [4];
   nd.v_alpha = mat -> oargs.farg [5];
   if (nd.u_alpha < FTINY || nd.v_alpha <= FTINY)
      objerror(mat, USER, "roughness too small");
      
   /* check for back side; reorient if back is visible */
   if (rayIn -> rod < 0)
      if (!backvis && mat -> otype != MAT_TRANS2) 
         return 0;
      else {
         /* get modifiers */
         raytexture(rayIn, mat -> omod);
         flipsurface(rayIn);
      }
   else raytexture(rayIn, mat -> omod);
   
   /* perturb normal */
   nd.pdot = max(raynormal(nd.pnorm, rayIn), .001);
   
   /* modify material color */
   multcolor(nd.mcolor, rayIn -> pcol);
   nd.rspec = mat -> oargs.farg [3];
   
   /* transmission params */
   if (mat -> otype == MAT_TRANS2) {
      nd.trans = mat -> oargs.farg [6] * (1 - nd.rspec);
      nd.tspec = nd.trans * mat -> oargs.farg [7];
      nd.tdiff = nd.trans - nd.tspec;
      if (nd.tspec > FTINY) 
         nd.specfl |= SP_TRAN;
   } 
   else nd.tdiff = nd.tspec = nd.trans = 0;
   
   /* specular reflection params */
   if (nd.rspec > FTINY) {
      nd.specfl |= SP_REFL;
      
      /* compute specular color */
      if (mat -> otype == MAT_METAL2) 
         copycolor(nd.scolor, nd.mcolor);
      else setcolor(nd.scolor, 1, 1, 1);
      
      scalecolor(nd.scolor, nd.rspec);
   }
   else setcolor(nd.scolor, 0, 0, 0);
   
   /* diffuse reflection params */
   nd.rdiff = 1 - nd.trans - nd.rspec;
   
   /* Set up probabilities */
   prdiff = ptdiff = ptspec = colorAvg(nd.mcolor);
   prdiff *= nd.rdiff;
   ptdiff *= nd.tdiff;
   prspec = colorAvg(nd.scolor);
   ptspec *= nd.tspec;
   albedo = prdiff + ptdiff + prspec + ptspec;
   
   /* Insert direct and indirect photon hits if diffuse component */
   if (prdiff > FTINY || ptdiff > FTINY)
      addPhotons(rayIn);
   
   xi = pmapRandom(rouletteState);
   
   if (xi > albedo) 
      /* Absorbed */
      return 0;
      
   if (xi > (albedo -= prspec))
      /* Specular reflection */
      if (!(nd.specfl & SP_BADU)) {
         photonRay(rayIn, &rayOut, PMAP_SPECREFL, nd.scolor);

         if (!anisoSpecPhotonScatter(&nd, &rayOut)) 
            return 0;
      }
      else return 0;
      
   else if (xi > (albedo -= ptspec))
      /* Specular transmission */
      
      if (!(nd.specfl & SP_BADU)) {
         /* Specular transmission */
         photonRay(rayIn, &rayOut, PMAP_SPECTRANS, nd.mcolor);

         if (!anisoSpecPhotonScatter(&nd, &rayOut)) 
            return 0;
      }
      else return 0;
      
   else if (xi > (albedo -= prdiff)) {
      /* Diffuse reflection */   
      photonRay(rayIn, &rayOut, PMAP_DIFFREFL, nd.mcolor);
      diffPhotonScatter(nd.pnorm, &rayOut);
   }
   
   else {
      /* Diffuse transmission */
      FVECT bnorm;
      flipsurface(rayIn);
      bnorm [0] = -nd.pnorm [0];
      bnorm [1] = -nd.pnorm [1];
      bnorm [2] = -nd.pnorm [2];
      
      photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, nd.mcolor);
      diffPhotonScatter(bnorm, &rayOut);
   }
   
   tracePhoton(&rayOut);
   return 0;
}


static double mylog (double x)
/* special log for extinction coefficients; cloned from dielectric.c */
{
   if (x < 1e-40)
      return(-100.);
      
   if (x >= 1.)
      return(0.);
      
   return(log(x));
}


static int dielectricPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for dielectric material and recurse */
{
   double   cos1, cos2, nratio, d1, d2, refl;
   COLOR    ctrans, talb;
   FVECT    dnorm;
   int      hastexture, i;
   RAY      rayOut;

   if (mat -> oargs.nfargs != (mat -> otype == MAT_DIELECTRIC ? 5 : 8))
      objerror(mat, USER, "bad arguments");
      
   /* get modifiers */
   raytexture(rayIn, mat -> omod);			
   
   if ((hastexture = (DOT(rayIn -> pert, rayIn -> pert) > sqr(FTINY))))
      /* Perturb normal */
      cos1 = raynormal(dnorm, rayIn);
   else {
      VCOPY(dnorm, rayIn -> ron);
      cos1 = rayIn -> rod;
   }
   
   /* index of refraction */
   nratio = mat -> otype == 
      MAT_DIELECTRIC ? mat->oargs.farg[3] + mat->oargs.farg[4] / MLAMBDA
                     : mat->oargs.farg[3] / mat->oargs.farg[7];
                     
   if (cos1 < 0) {
      /* inside */
      hastexture = -hastexture;
      cos1 = -cos1;
      dnorm [0] = -dnorm [0];
      dnorm [1] = -dnorm [1];
      dnorm [2] = -dnorm [2];            
      setcolor(rayIn -> cext, 
               -mylog(mat -> oargs.farg [0] * rayIn -> pcol [0]),
               -mylog(mat -> oargs.farg [1] * rayIn -> pcol [1]),
               -mylog(mat -> oargs.farg [2] * rayIn -> pcol [2]));
      setcolor(rayIn -> albedo, 0, 0, 0);
      rayIn -> gecc = 0;
      
      if (mat -> otype == MAT_INTERFACE) {
         setcolor(ctrans, 
                  -mylog(mat -> oargs.farg [4] * rayIn -> pcol [0]),
                  -mylog(mat -> oargs.farg [5] * rayIn -> pcol [1]),
                  -mylog(mat -> oargs.farg [6] * rayIn -> pcol [2]));
         setcolor(talb, 0, 0, 0);
      } 
      else {
         copycolor(ctrans, cextinction);
         copycolor(talb, salbedo);
      }
   } 
   
   else {
      /* outside */
      nratio = 1.0 / nratio;
      setcolor(ctrans, 
               -mylog(mat -> oargs.farg [0] * rayIn -> pcol [0]),
               -mylog(mat -> oargs.farg [1] * rayIn -> pcol [1]),
               -mylog(mat -> oargs.farg [2] * rayIn -> pcol [2]));
      setcolor(talb, 0, 0, 0);
      
      if (mat -> otype == MAT_INTERFACE) {
         setcolor(rayIn -> cext,
                  -mylog(mat -> oargs.farg [4] * rayIn -> pcol [0]),
                  -mylog(mat -> oargs.farg [5] * rayIn -> pcol [1]),
                  -mylog(mat -> oargs.farg [6] * rayIn -> pcol [2]));
         setcolor(rayIn -> albedo, 0, 0, 0);
         rayIn -> gecc = 0;
      }
   } 
              
   /* compute cos theta2 */
   d2 = 1 - sqr(nratio) * (1 - sqr(cos1));
   
   if (d2 < FTINY) {
      /* Total reflection */	
      refl = cos2 = 1.0;
   }
   else {
      /* Refraction, compute Fresnel's equations */
      cos2 = sqrt(d2);
      d1 = cos1;
      d2 = nratio * cos2;
      d1 = (d1 - d2) / (d1 + d2);
      refl = sqr(d1);
      d1 = 1 / cos1;
      d2 = nratio / cos2;
      d1 = (d1 - d2) / (d1 + d2);
      refl += sqr(d1);
      refl *= 0.5;
   }
   
   if (pmapRandom(rouletteState) > refl) {
      /* Refraction */
      photonRay(rayIn, &rayOut, PMAP_REFRACT, NULL);
      d1 = nratio * cos1 - cos2;
      
      for (i = 0; i < 3; i++)
         rayOut.rdir [i] = nratio * rayIn -> rdir [i] + d1 * dnorm [i];
         
      if (hastexture && DOT(rayOut.rdir, rayIn->ron)*hastexture >= -FTINY) {
         d1 *= hastexture;
         
         for (i = 0; i < 3; i++)
            rayOut.rdir [i] = nratio * rayIn -> rdir [i] + 
                              d1 * rayIn -> ron [i];
                              
         normalize(rayOut.rdir);
      }
      
      copycolor(rayOut.cext, ctrans);
      copycolor(rayOut.albedo, talb);
   }
   
   else {
      /* Reflection */
      photonRay(rayIn, &rayOut, PMAP_SPECREFL, NULL);
      VSUM(rayOut.rdir, rayIn -> rdir, dnorm, 2 * cos1);
      
      if (hastexture && DOT(rayOut.rdir, rayIn->ron) * hastexture <= FTINY)
         for (i = 0; i < 3; i++)
            rayOut.rdir [i] = rayIn -> rdir [i] + 
                              2 * rayIn -> rod * rayIn -> ron [i];
   }
   
   /* Ray is modified by medium defined by cext and albedo in
    * photonParticipate() */
   tracePhoton(&rayOut);
   
   return 0;
}



static int glassPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for glass material and recurse */
{
   float    albedo, xi, ptrans;
   COLOR    mcolor, refl, trans;
   double   pdot, cos2, d, r1e, r1m, rindex = 0.0;
   FVECT    pnorm, pdir;
   int      hastexture, i;
   RAY      rayOut;

   /* check arguments */
   if (mat -> oargs.nfargs == 3)
	   rindex = RINDEX;
   else if (mat -> oargs.nfargs == 4)
	   rindex = mat -> oargs.farg [3];
   else objerror(mat, USER, "bad arguments");

   copycolor(mcolor, mat -> oargs.farg);   
   
   /* get modifiers */
   raytexture(rayIn, mat -> omod);
   
   /* reorient if necessary */
   if (rayIn -> rod < 0) 
      flipsurface(rayIn);
   if ((hastexture = (DOT(rayIn -> pert, rayIn -> pert) > sqr(FTINY))))
      pdot = raynormal(pnorm, rayIn);
   else {
      VCOPY(pnorm, rayIn -> ron);
      pdot = rayIn -> rod;
   }
   
   /* Modify material color */
   multcolor(mcolor, rayIn -> pcol);
   
   /* angular transmission */
   cos2 = sqrt((1 - 1 / sqr(rindex)) + sqr(pdot / rindex));
   setcolor(mcolor, pow(mcolor [0], 1 / cos2), pow(mcolor [1], 1 / cos2), 
            pow(mcolor [2], 1 / cos2));
            
   /* compute reflection */
   r1e = (pdot - rindex * cos2) / (pdot + rindex * cos2);
   r1e *= r1e;
   r1m = (1 / pdot - rindex / cos2) / (1 / pdot + rindex / cos2);
   r1m *= r1m;
   
   for (i = 0; i < 3; i++) {
      double r1ed2, r1md2, d2;
      
      d = mcolor [i];
      d2 = sqr(d);
      r1ed2 = sqr(r1e) * d2;
      r1md2 = sqr(r1m) * d2;
      
      /* compute transmittance */
      trans [i] = 0.5 * d * 
                  (sqr(1 - r1e) / (1 - r1ed2) + sqr(1 - r1m) / (1 - r1md2));
                             
      /* compute reflectance */
      refl [i] = 0.5 * (r1e * (1 + (1 - 2 * r1e) * d2) / (1 - r1ed2) +
                        r1m * (1 + (1 - 2 * r1m) * d2) / (1 - r1md2));
   }
   
   /* Set up probabilities */
   ptrans = colorAvg(trans);
   albedo = colorAvg(refl) + ptrans;
   xi = pmapRandom(rouletteState);
   

   if (xi > albedo) 
      /* Absorbed */
      return 0;
      
   if (xi > (albedo -= ptrans)) {
      /* Transmitted */

      if (hastexture) {
         /* perturb direction */
         VSUM(pdir, rayIn -> rdir, rayIn -> pert, 2 * (1 - rindex));
         
         if (normalize(pdir) == 0) {
            objerror(mat, WARNING, "bad perturbation");
            VCOPY(pdir, rayIn -> rdir);
         }
      } 
      else VCOPY(pdir, rayIn -> rdir);
      
      VCOPY(rayOut.rdir, pdir);
      photonRay(rayIn, &rayOut, PMAP_SPECTRANS, mcolor);      
   }
   
   else {
      /* reflected ray */
      VSUM(rayOut.rdir, rayIn -> rdir, pnorm, 2 * pdot);
      photonRay(rayIn, &rayOut, PMAP_SPECREFL, mcolor);      
   }

   tracePhoton(&rayOut);
   return 0;
}



static int aliasPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Transfer photon scattering to alias target */
{
   OBJECT   aliasObj;
   OBJREC   aliasRec, *aliasPtr;
   
   /* Straight replacement? */
   if (!mat -> oargs.nsargs) {
      /* Skip void modifier! */
      if (mat -> omod != OVOID) {   
         mat = objptr(mat -> omod);
         photonScatter [mat -> otype] (mat, rayIn);
      }
      
      return 0;
   }
   
   /* Else replace alias */
   if (mat -> oargs.nsargs != 1)
      objerror(mat, INTERNAL, "bad # string arguments");
      
   aliasPtr = mat;
   aliasObj = objndx(aliasPtr);
   
   /* Follow alias trail */
   do {
      aliasObj = aliasPtr -> oargs.nsargs == 1 
                     ? lastmod(aliasObj, aliasPtr -> oargs.sarg [0])
                     : aliasPtr -> omod;
      if (aliasObj < 0)
         objerror(aliasPtr, USER, "bad reference");
         
      aliasPtr = objptr(aliasObj);
   } while (aliasPtr -> otype == MOD_ALIAS);

   /* Copy alias object */
   aliasRec = *aliasPtr;
   
   /* Substitute modifier */
   aliasRec.omod = mat -> omod;
   
   /* Replacement scattering routine */
   photonScatter [aliasRec.otype] (&aliasRec, rayIn);

   /* Avoid potential memory leak? */
   if (aliasRec.os != aliasPtr -> os) {
      if (aliasPtr -> os)
         free_os(aliasPtr);
      aliasPtr -> os = aliasRec.os;
   }

   return 0;
}



static int clipPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for antimatter material and recurse */
{
   OBJECT      obj = objndx(mat), mod, cset [MAXSET + 1], *modset;
   int         entering, inside = 0, i;
   const RAY   *rp;
   RAY         rayOut;

   if ((modset = (OBJECT*)mat -> os) == NULL) {   
      if (mat -> oargs.nsargs < 1 || mat -> oargs.nsargs > MAXSET)
         objerror(mat, USER, "bad # arguments");
         
      modset = (OBJECT*)malloc((mat -> oargs.nsargs + 1) * sizeof(OBJECT));
      
      if (modset == NULL) 
         error(SYSTEM, "out of memory in clipPhotonScatter");
      modset [0] = 0;
      
      for (i = 0; i < mat -> oargs.nsargs; i++) {
         if (!strcmp(mat -> oargs.sarg [i], VOIDID))
            continue;
            
         if ((mod = lastmod(obj, mat -> oargs.sarg [i])) == OVOID) {
            sprintf(errmsg, "unknown modifier \"%s\"", mat->oargs.sarg[i]);
            objerror(mat, WARNING, errmsg);
            continue;
         }
         
         if (inset(modset, mod)) {
            objerror(mat, WARNING, "duplicate modifier");
            continue;
         }
         
         insertelem(modset, mod);
      }
      
      mat -> os = (char*)modset;
   }
   
   if (rayIn -> clipset != NULL) 
      setcopy(cset, rayIn -> clipset);
   else cset [0] = 0;
   
   entering = rayIn -> rod > 0;
   
   /* Store photon incident from front if material defined as sensor */
   if (entering && inset(photonSensorSet, obj))
      addPhotons(rayIn);

   for (i = modset [0]; i > 0; i--) {
      if (entering) {
         if (!inset(cset, modset [i])) {
            if (cset [0] >= MAXSET) 
               error(INTERNAL, "set overflow in clipPhotonScatter");
            insertelem(cset, modset [i]);
         }
      } 
      else if (inset(cset, modset [i])) 
         deletelem(cset, modset [i]);
   }
   
   rayIn -> newcset = cset;
   
   if (strcmp(mat -> oargs.sarg [0], VOIDID)) {   
	   for (rp = rayIn; rp -> parent != NULL; rp = rp -> parent) {	   
		   if ( !(rp -> rtype & RAYREFL) && rp->parent->ro != NULL && 
				inset(modset, rp -> parent -> ro -> omod)) {	

			   if (rp -> parent -> rod > 0)
				   inside++;
			   else inside--;
		   }
	   }
	   
	   if (inside > 0) {
		   flipsurface(rayIn);
		   mat = objptr(lastmod(obj, mat -> oargs.sarg [0]));
		   photonScatter [mat -> otype] (mat, rayIn);
		   return 0;
	   }
   }
   
   /* Else transfer ray */
   photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
   tracePhoton(&rayOut);
   
   return 0;
}



static int mirrorPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for mirror material and recurse */
{
   RAY      rayOut;
   int      rpure = 1, i;
   FVECT    pnorm;
   double   pdot;
   float    albedo;
   COLOR    mcolor;

   /* check arguments */
   if (mat -> oargs.nfargs != 3 || mat -> oargs.nsargs > 1)
      objerror(mat, USER, "bad number of arguments");
      
   /* back is black */
   if (rayIn -> rod < 0) 
      return 0;
      
   /* get modifiers */
   raytexture(rayIn, mat -> omod);
   
   /* assign material color */
   copycolor(mcolor, mat -> oargs.farg); 
   multcolor(mcolor, rayIn -> pcol);
   
   /* Set up probabilities */
   albedo = colorAvg(mcolor);
   
   if (pmapRandom(rouletteState) > albedo) 
      /* Absorbed */
      return 0;
      
   /* compute reflected ray */
   photonRay(rayIn, &rayOut, PMAP_SPECREFL, mcolor);

   if (DOT(rayIn -> pert, rayIn -> pert) > sqr(FTINY)) {
      /* use textures */
      pdot = raynormal(pnorm, rayIn);   
      
      for (i = 0; i < 3; i++)
         rayOut.rdir [i] = rayIn -> rdir [i] + 2 * pdot * pnorm [i];
         
      rpure = 0;
   }
   
   /* Check for penetration */
   if (rpure || DOT(rayOut.rdir, rayIn -> ron) <= FTINY)
      for (i = 0; i < 3; i++)
         rayOut.rdir [i] = rayIn -> rdir [i] + 
                           2 * rayIn -> rod * rayIn -> ron [i];
                           
   tracePhoton(&rayOut);
   return 0;
}



static int mistPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray within mist and recurse */
{
   COLOR    mext;
   RREAL    re, ge, be;
   RAY      rayOut;

   /* check arguments */
   if (mat -> oargs.nfargs > 7) 
      objerror(mat, USER, "bad arguments");
      
   if (mat -> oargs.nfargs > 2) {
      /* compute extinction */
      copycolor(mext, mat -> oargs.farg);
      /* get modifiers */
      raytexture(rayIn, mat -> omod);
      multcolor(mext, rayIn -> pcol);
   } 
   else setcolor(mext, 0, 0, 0);
      
   photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
   
   if (rayIn -> rod > 0) {
      /* entering ray */
      addcolor(rayOut.cext, mext);
      
      if (mat -> oargs.nfargs > 5) 
         copycolor(rayOut.albedo, mat -> oargs.farg + 3);
      if (mat -> oargs.nfargs > 6) 
         rayOut.gecc = mat -> oargs.farg [6];
   } 
   
   else {
      /* leaving ray */
      re = max(rayIn -> cext [0] - mext [0], cextinction [0]);
      ge = max(rayIn -> cext [1] - mext [1], cextinction [1]);
      be = max(rayIn -> cext [2] - mext [2], cextinction [2]);
      setcolor(rayOut.cext, re, ge, be);
      
      if (mat -> oargs.nfargs > 5) 
         copycolor(rayOut.albedo, salbedo);
      if (mat -> oargs.nfargs > 6)
         rayOut.gecc = seccg;
   }
   
   tracePhoton(&rayOut);
   
   return 0;
}



static int mx_dataPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Pass photon on to materials selected by mixture data */
{
   OBJECT   obj;
   double   coef, pt [MAXDIM];
   DATARRAY *dp;
   OBJECT   mod [2];
   MFUNC    *mf;
   int      i;

   if (mat -> oargs.nsargs < 6) 
      objerror(mat, USER, "bad # arguments");
      
   obj = objndx(mat);
   
   for (i = 0; i < 2; i++)
      if (!strcmp(mat -> oargs.sarg [i], VOIDID)) 
         mod [i] = OVOID;
      else if ((mod [i] = lastmod(obj, mat -> oargs.sarg [i])) == OVOID) {
         sprintf(errmsg, "undefined modifier \"%s\"", mat->oargs.sarg[i]);
         objerror(mat, USER, errmsg);
      }
      
   dp = getdata(mat -> oargs.sarg [3]);
   i = (1 << dp -> nd) - 1;
   mf = getfunc(mat, 4, i << 5, 0);
   setfunc(mat, rayIn);
   errno = 0;
   
   for (i = 0; i < dp -> nd; i++) {
      pt [i] = evalue(mf -> ep [i]);
      
      if (errno) {
         objerror(mat, WARNING, "compute error");
         return 0;
      }
   }
   
   coef = datavalue(dp, pt);
   errno = 0;
   coef = funvalue(mat -> oargs.sarg [2], 1, &coef);
   
   if (errno) 
      objerror(mat, WARNING, "compute error");
   else {
      OBJECT mxMod = mod [pmapRandom(rouletteState) < coef ? 0 : 1];
      
      if (mxMod != OVOID) {
         mat = objptr(mxMod);
         photonScatter [mat -> otype] (mat, rayIn);
      }
      else {
         /* Transfer ray if no modifier */
         RAY rayOut;
         
         photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
         tracePhoton(&rayOut);      
      }            
   }
   
   return 0;
}



static int mx_pdataPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Pass photon on to materials selected by mixture picture */
{
   OBJECT   obj;
   double   col [3], coef, pt [MAXDIM];
   DATARRAY *dp;
   OBJECT   mod [2];
   MFUNC    *mf;
   int      i;

   if (mat -> oargs.nsargs < 7) 
      objerror(mat, USER, "bad # arguments");
      
   obj = objndx(mat);
   
   for (i = 0; i < 2; i++)
      if (!strcmp(mat -> oargs.sarg [i], VOIDID)) 
         mod [i] = OVOID;
      else if ((mod [i] = lastmod(obj, mat -> oargs.sarg [i])) == OVOID) {
         sprintf(errmsg, "undefined modifier \"%s\"", mat->oargs.sarg[i]);
         objerror(mat, USER, errmsg);
      }
      
   dp = getpict(mat -> oargs.sarg [3]);
   mf = getfunc(mat, 4, 0x3 << 5, 0);
   setfunc(mat, rayIn);
   errno = 0;
   pt [1] = evalue(mf -> ep [0]);
   pt [0] = evalue(mf -> ep [1]);
   
   if (errno) {
      objerror(mat, WARNING, "compute error");
      return 0;
   }
   
   for (i = 0; i < 3; i++) 
      col [i] = datavalue(dp + i, pt);
      
   errno = 0;
   coef = funvalue(mat -> oargs.sarg [2], 3, col);
   
   if (errno) 
      objerror(mat, WARNING, "compute error");
   else {
      OBJECT mxMod = mod [pmapRandom(rouletteState) < coef ? 0 : 1];
      
      if (mxMod != OVOID) {
         mat = objptr(mxMod);
         photonScatter [mat -> otype] (mat, rayIn);
      }
      else {
         /* Transfer ray if no modifier */
         RAY rayOut;
         
         photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
         tracePhoton(&rayOut);      
      }      
   }   
   
   return 0;
}



static int mx_funcPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Pass photon on to materials selected by mixture function */
{
   OBJECT   obj, mod [2];
   int      i;
   double   coef;
   MFUNC    *mf;

   if (mat -> oargs.nsargs < 4) 
      objerror(mat, USER, "bad # arguments");
      
   obj = objndx(mat);
   
   for (i = 0; i < 2; i++)
      if (!strcmp(mat -> oargs.sarg [i], VOIDID))
         mod [i] = OVOID;
      else if ((mod [i] = lastmod(obj, mat -> oargs.sarg [i])) == OVOID) {
         sprintf(errmsg, "undefined modifier \"%s\"", mat->oargs.sarg[i]);
         objerror(mat, USER, errmsg);
      }
      
   mf = getfunc(mat, 3, 0x4, 0);
   setfunc(mat, rayIn);
   errno = 0;
   
   /* bound coefficient */
   coef = min(1, max(0, evalue(mf -> ep [0])));
   
   if (errno) 
      objerror(mat, WARNING, "compute error");
   else {         
      OBJECT mxMod = mod [pmapRandom(rouletteState) < coef ? 0 : 1];
      
      if (mxMod != OVOID) {
         mat = objptr(mxMod);
         photonScatter [mat -> otype] (mat, rayIn);
      }
      else {
         /* Transfer ray if no modifier */
         RAY rayOut;
         
         photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
         tracePhoton(&rayOut);      
      }
   }
   
   return 0;
}



static int pattexPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for pattern or texture modifier and recurse.
   This code is brought to you by Henkel! :^) */
{
   RAY   rayOut;
   
   /* Get pattern */
   ofun [mat -> otype].funp(mat, rayIn);
   if (mat -> omod != OVOID) {
      /* Scatter using modifier (if any) */
      mat = objptr(mat -> omod);
      photonScatter [mat -> otype] (mat, rayIn);
   }
   else {
      /* Transfer ray if no modifier */
      photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
      tracePhoton(&rayOut);
   }
   
   return 0;
}



static int setbrdfunc(BRDFDAT *bd)
/* Set up brdf function and variables; ripped off from m_brdf.c */
{
   FVECT v;
   
   if (setfunc(bd -> mp, bd -> pr) == 0)
      return 0;

   /* (Re)Assign func variables */
   multv3(v, bd -> pnorm, funcxf.xfm);
   varset("NxP", '=', v [0] / funcxf.sca);
   varset("NyP", '=', v [1] / funcxf.sca);
   varset("NzP", '=', v [2] / funcxf.sca);
   varset("RdotP", '=', 
          bd -> pdot <= -1. ? -1. : bd -> pdot >= 1. ? 1. : bd -> pdot);
   varset("CrP", '=', colval(bd -> mcolor, RED));
   varset("CgP", '=', colval(bd -> mcolor, GRN));
   varset("CbP", '=', colval(bd -> mcolor, BLU));
   
   return 1;
}



static int brdfPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for BRTDfunc material and recurse. Only ideal
   reflection and transmission are sampled for the specular componentent. */
{
   int      hitfront = 1, hastexture, i;
   BRDFDAT  nd;
   RAY      rayOut;
   COLOR    rspecCol, tspecCol;
   double   prDiff, ptDiff, prSpec, ptSpec, albedo, xi;
   MFUNC    *mf;
   FVECT    bnorm;

   /* Check argz */
   if (mat -> oargs.nsargs < 10 || mat -> oargs.nfargs < 9)
      objerror(mat, USER, "bad # arguments");
      
   nd.mp = mat;
   nd.pr = rayIn;
   /* Dummiez */
   nd.rspec = nd.tspec = 1.0;
   nd.trans = 0.5;

   /* Diffuz reflektanz */
   if (rayIn -> rod > 0.0)
      setcolor(nd.rdiff, mat -> oargs.farg[0], mat -> oargs.farg [1],
               mat -> oargs.farg [2]);
   else
      setcolor(nd.rdiff, mat-> oargs.farg [3], mat -> oargs.farg [4],
               mat -> oargs.farg [5]);
   /* Diffuz tranzmittanz */
   setcolor(nd.tdiff, mat -> oargs.farg [6], mat -> oargs.farg [7],
            mat -> oargs.farg [8]);

   /* Get modz */
   raytexture(rayIn, mat -> omod);
   hastexture = (DOT(rayIn -> pert, rayIn -> pert) > sqr(FTINY));
   if (hastexture) {
      /* Perturb normal */
      nd.pdot = raynormal(nd.pnorm, rayIn);
   } 
   else {
      VCOPY(nd.pnorm, rayIn -> ron);
      nd.pdot = rayIn -> rod;
   }

   if (rayIn -> rod < 0.0) {
      /* Orient perturbed valuz */
      nd.pdot = -nd.pdot;
      for (i = 0; i < 3; i++) {
         nd.pnorm [i] = -nd.pnorm [i];
         rayIn -> pert [i] = -rayIn -> pert [i];
      }
      
      hitfront = 0;
   }
   
   /* Get pattern kolour, modify diffuz valuz */
   copycolor(nd.mcolor, rayIn -> pcol);
   multcolor(nd.rdiff, nd.mcolor);
   multcolor(nd.tdiff, nd.mcolor);

   /* Load cal file, evaluate spekula refl/tranz varz */
   nd.dp = NULL;
   mf = getfunc(mat, 9, 0x3f, 0);
   setbrdfunc(&nd);
   errno = 0;
   setcolor(rspecCol, 
            evalue(mf->ep[0]), evalue(mf->ep[1]), evalue(mf->ep[2]));
   setcolor(tspecCol, 
            evalue(mf->ep[3]), evalue(mf->ep[4]), evalue(mf->ep[5]));
   if (errno == EDOM || errno == ERANGE)
      objerror(mat, WARNING, "compute error");
   else {
      /* Set up probz */
      prDiff = colorAvg(nd.rdiff);
      ptDiff = colorAvg(nd.tdiff);
      prSpec = colorAvg(rspecCol);
      ptSpec = colorAvg(tspecCol);
      albedo = prDiff + ptDiff + prSpec + ptSpec;
   }

   /* Insert direct and indirect photon hitz if diffuz komponent */
   if (prDiff > FTINY || ptDiff > FTINY)
      addPhotons(rayIn);

   /* Stochastically sample absorption or scattering evenz */
   if ((xi = pmapRandom(rouletteState)) > albedo)
      /* Absorbed */
      return 0;

   if (xi > (albedo -= prSpec)) {
      /* Ideal spekula reflekzion */
      photonRay(rayIn, &rayOut, PMAP_SPECREFL, rspecCol);
      VSUM(rayOut.rdir, rayIn -> rdir, nd.pnorm, 2 * nd.pdot);
      checknorm(rayOut.rdir);
   }
   else if (xi > (albedo -= ptSpec)) {
      /* Ideal spekula tranzmission */
      photonRay(rayIn, &rayOut, PMAP_SPECTRANS, tspecCol);
      if (hastexture) {
         /* Perturb direkzion */
         VSUB(rayOut.rdir, rayIn -> rdir, rayIn -> pert);
         if (normalize(rayOut.rdir) == 0.0) {
            objerror(mat, WARNING, "illegal perturbation");
            VCOPY(rayOut.rdir, rayIn -> rdir);
         }
      }
      else VCOPY(rayOut.rdir, rayIn -> rdir);
   }
   else if (xi > (albedo -= prDiff)) {
      /* Diffuz reflekzion */
      if (!hitfront)
         flipsurface(rayIn);
      photonRay(rayIn, &rayOut, PMAP_DIFFREFL, nd.mcolor);
      diffPhotonScatter(nd.pnorm, &rayOut);
   }
   else {
      /* Diffuz tranzmission */
      if (hitfront)
         flipsurface(rayIn);
      photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, nd.mcolor);
      bnorm [0] = -nd.pnorm [0];
      bnorm [1] = -nd.pnorm [1];
      bnorm [2] = -nd.pnorm [2];
      diffPhotonScatter(bnorm, &rayOut); 
   }

   tracePhoton(&rayOut);
   return 0;
}



int brdf2PhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for procedural or data driven BRDF material and
   recurse. Only diffuse reflection and transmission are sampled. */
{
   BRDFDAT  nd;
   RAY      rayOut;
   double   dtmp, prDiff, ptDiff, albedo, xi;
   MFUNC    *mf;
   FVECT    bnorm;

   /* Check argz */
   if (mat -> oargs.nsargs < (hasdata(mat -> otype) ? 4 : 2) || 
       mat -> oargs.nfargs < (mat -> otype == MAT_TFUNC || 
                              mat -> otype == MAT_TDATA ? 6 : 4))
      objerror(mat, USER, "bad # arguments");
      
   if (rayIn -> rod < 0.0) {
      /* Hit backside; reorient if visible, else transfer photon */
      if (!backvis) {
         photonRay(rayIn, &rayOut, PMAP_XFER, NULL);
         tracePhoton(&rayOut);
         return 0;
      }
      
      raytexture(rayIn, mat -> omod);
      flipsurface(rayIn);
   } 
   else raytexture(rayIn, mat -> omod);

   nd.mp = mat;
   nd.pr = rayIn;
   
   /* Material kolour */
   setcolor(nd.mcolor, mat -> oargs.farg [0], mat -> oargs.farg [1], 
            mat -> oargs.farg [2]);
   /* Spekula komponent */
   nd.rspec = mat -> oargs.farg [3];
   
   /* Tranzmittanz */
   if (mat -> otype == MAT_TFUNC || mat -> otype == MAT_TDATA) {
      nd.trans = mat -> oargs.farg [4] * (1.0 - nd.rspec);
      nd.tspec = nd.trans * mat -> oargs.farg [5];
      dtmp = nd.trans - nd.tspec;
      setcolor(nd.tdiff, dtmp, dtmp, dtmp);
   } 
   else {
      nd.tspec = nd.trans = 0.0;
      setcolor(nd.tdiff, 0.0, 0.0, 0.0);
   }
   
   /* Reflektanz */
   dtmp = 1.0 - nd.trans - nd.rspec;
   setcolor(nd.rdiff, dtmp, dtmp, dtmp);
   /* Perturb normal */
   nd.pdot = raynormal(nd.pnorm, rayIn);
   /* Modify material kolour */
   multcolor(nd.mcolor, rayIn -> pcol);
   multcolor(nd.rdiff, nd.mcolor);
   multcolor(nd.tdiff, nd.mcolor);
   
   /* Load auxiliary filez */
   if (hasdata(mat -> otype)) {
      nd.dp = getdata(mat -> oargs.sarg [1]);
      getfunc(mat, 2, 0, 0);
   } 
   else {
      nd.dp = NULL;
      getfunc(mat, 1, 0, 0);
   }

   /* Set up probz */
   prDiff = colorAvg(nd.rdiff);
   ptDiff = colorAvg(nd.tdiff);
   albedo = prDiff + ptDiff;

   /* Insert direct and indirect photon hitz if diffuz komponent */
   if (prDiff > FTINY || ptDiff > FTINY)
      addPhotons(rayIn);

   /* Stochastically sample absorption or scattering evenz */
   if ((xi = pmapRandom(rouletteState)) > albedo)
      /* Absorbed */
      return 0;

   if (xi > (albedo -= prDiff)) {
      /* Diffuz reflekzion */
      photonRay(rayIn, &rayOut, PMAP_DIFFREFL, nd.rdiff);
      diffPhotonScatter(nd.pnorm, &rayOut);
   }
   else {
      /* Diffuz tranzmission */
      flipsurface(rayIn);
      photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, nd.tdiff);
      bnorm [0] = -nd.pnorm [0];
      bnorm [1] = -nd.pnorm [1];
      bnorm [2] = -nd.pnorm [2];
      diffPhotonScatter(bnorm, &rayOut); 
   }

   tracePhoton(&rayOut);
   return 0;
}



/* 
   ======================================================================
   The following code is
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ======================================================================
*/

static int bsdfPhotonScatter (OBJREC *mat, RAY *rayIn)
/* Generate new photon ray for BSDF modifier and recurse. */
{
   int            hasthick = (mat->otype == MAT_BSDF);
   int            hitFront;
   SDError        err;
   SDValue        bsdfVal;
   FVECT          upvec;
   MFUNC          *mf;
   BSDFDAT        nd;
   RAY            rayOut;
   COLOR          bsdfRGB;
   int            transmitted;
   double         prDiff, ptDiff, prDiffSD, ptDiffSD, prSpecSD, ptSpecSD, 
                  albedo, xi;
   const double   patAlb = bright(rayIn -> pcol);
   
   /* Following code adapted from m_bsdf() */
   /* Check arguments */
   if (
      mat -> oargs.nsargs < hasthick+5 || 
      mat -> oargs.nfargs > 9 || mat -> oargs.nfargs % 3
   )
      objerror(mat, USER, "bad # arguments");
      
   hitFront = (rayIn -> rod > 0);

   /* Load cal file */
   mf = hasthick ? getfunc(mat, 5, 0x1d, 1) : getfunc(mat, 4, 0xe, 1);

   /* Get thickness */
   nd.thick = 0;
   if (hasthick) {
      nd.thick = evalue(mf -> ep [0]);
      if ((-FTINY <= nd.thick) & (nd.thick <= FTINY))
         nd.thick = .0;
   }

   /* Get BSDF data */
   nd.sd = loadBSDF(mat -> oargs.sarg [hasthick]);
   
   /* Extra diffuse reflectance from material def */
   if (hitFront) {
      if (mat -> oargs.nfargs < 3)
         setcolor(nd.rdiff, .0, .0, .0);
      else setcolor(
         nd.rdiff, 
         mat -> oargs.farg [0], mat -> oargs.farg [1], mat -> oargs.farg [2]
      );
   }    
   else if (mat -> oargs.nfargs < 6) {
      /* Check for absorbing backside */
      if (!backvis && !nd.sd -> rb && !nd.sd -> tf) {
         SDfreeCache(nd.sd);
         return 0;
      }
      
      setcolor(nd.rdiff, .0, .0, .0);
   } 
   else setcolor(
      nd.rdiff, 
      mat -> oargs.farg [3], mat -> oargs.farg [4], mat -> oargs.farg [5]
   );

   /* Extra diffuse transmittance from material def */
   if (mat -> oargs.nfargs < 9)
      setcolor(nd.tdiff, .0, .0, .0);
   else setcolor(
      nd.tdiff, 
      mat -> oargs.farg [6], mat -> oargs.farg [7], mat -> oargs.farg [8]
   );
               
   nd.mp = mat;
   nd.pr = rayIn;

   /* Get modifiers */
   raytexture(rayIn, mat -> omod);
   
   /* Modify diffuse values */
   multcolor(nd.rdiff, rayIn -> pcol);
   multcolor(nd.tdiff, rayIn -> pcol);

   /* Get up vector & xform to world coords */
   upvec [0] = evalue(mf -> ep [hasthick+0]);
   upvec [1] = evalue(mf -> ep [hasthick+1]);
   upvec [2] = evalue(mf -> ep [hasthick+2]);
   
   if (mf -> fxp != &unitxf) {
      multv3(upvec, upvec, mf -> fxp -> xfm);
      nd.thick *= mf -> fxp -> sca;
   }
   
   if (rayIn -> rox) {
      multv3(upvec, upvec, rayIn -> rox -> f.xfm);
      nd.thick *= rayIn -> rox -> f.sca;
   }
   
   /* Perturb normal */
   raynormal(nd.pnorm, rayIn);
   
   /* Xform incident dir to local BSDF coords */
   err = SDcompXform(nd.toloc, nd.pnorm, upvec);
   
   if (!err) {
      nd.vray [0] = -rayIn -> rdir [0];
      nd.vray [1] = -rayIn -> rdir [1];
      nd.vray [2] = -rayIn -> rdir [2];
      err = SDmapDir(nd.vray, nd.toloc, nd.vray);
   }
   
   if (!err)
      err = SDinvXform(nd.fromloc, nd.toloc);
      
   if (err) {
      objerror(mat, WARNING, "Illegal orientation vector");
      return 0;
   }
   
   /* Determine BSDF resolution */
   err = SDsizeBSDF(
      nd.sr_vpsa, nd.vray, NULL, SDqueryMin + SDqueryMax, nd.sd
   );
   
   if (err)
      objerror(mat, USER, transSDError(err));
      
   nd.sr_vpsa [0] = sqrt(nd.sr_vpsa [0]);
   nd.sr_vpsa [1] = sqrt(nd.sr_vpsa [1]);

   /* Orient perturbed normal towards incident side */
   if (!hitFront) {
      nd.pnorm [0] = -nd.pnorm [0];
      nd.pnorm [1] = -nd.pnorm [1];
      nd.pnorm [2] = -nd.pnorm [2];
   }

   /* Get scatter probabilities (weighted by pattern except for spec refl)
    * prDiff, ptDiff:      extra diffuse component in material def
    * prDiffSD, ptDiffSD:  diffuse (constant) component in SDF
    * prSpecSD, ptSpecSD:  non-diffuse ("specular") component in SDF 
    * albedo:              sum of above, inverse absorption probability */
   prDiff   = colorAvg(nd.rdiff);
   ptDiff   = colorAvg(nd.tdiff);
   prDiffSD = patAlb * SDdirectHemi(nd.vray, SDsampDf | SDsampR, nd.sd);
   ptDiffSD = patAlb * SDdirectHemi(nd.vray, SDsampDf | SDsampT, nd.sd);
   prSpecSD = SDdirectHemi(nd.vray, SDsampSp | SDsampR, nd.sd);
   ptSpecSD = patAlb * SDdirectHemi(nd.vray, SDsampSp | SDsampT, nd.sd);
   albedo   = prDiff + ptDiff + prDiffSD + ptDiffSD + prSpecSD + ptSpecSD;

   /*    
   if (albedo > 1)
      objerror(mat, WARNING, "Invalid albedo");
   */
         
   /* Insert direct and indirect photon hits if diffuse component */
   if (prDiff + ptDiff + prDiffSD + ptDiffSD > FTINY)
      addPhotons(rayIn);	

   xi = pmapRandom(rouletteState);
      
   if (xi > albedo)
      /* Absorbtion */
      return 0;
   
   transmitted = 0;

   if ((xi -= prDiff) <= 0) {
      /* Diffuse reflection (extra component in material def) */
      photonRay(rayIn, &rayOut, PMAP_DIFFREFL, nd.rdiff);
      diffPhotonScatter(nd.pnorm, &rayOut);
   }
   
   else if ((xi -= ptDiff) <= 0) {
      /* Diffuse transmission (extra component in material def) */
      photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, nd.tdiff);
      diffPhotonScatter(nd.pnorm, &rayOut);
      transmitted = 1;
   }
  
   else {   /* Sample SDF */
      if ((xi -= prDiffSD) <= 0) {
         /* Diffuse SDF reflection (constant component) */
         if ((err = SDsampBSDF(
            &bsdfVal, nd.vray, pmapRandom(scatterState), 
            SDsampDf | SDsampR, nd.sd
         )))
            objerror(mat, USER, transSDError(err));
         
         /* Apply pattern to spectral component */
         ccy2rgb(&bsdfVal.spec, bsdfVal.cieY, bsdfRGB);
         multcolor(bsdfRGB, rayIn -> pcol);
         photonRay(rayIn, &rayOut, PMAP_DIFFREFL, bsdfRGB);
      }

      else if ((xi -= ptDiffSD) <= 0) {
         /* Diffuse SDF transmission (constant component) */
         if ((err = SDsampBSDF(
            &bsdfVal, nd.vray, pmapRandom(scatterState), 
            SDsampDf | SDsampT, nd.sd
         )))
            objerror(mat, USER, transSDError(err));
         
         /* Apply pattern to spectral component */
         ccy2rgb(&bsdfVal.spec, bsdfVal.cieY, bsdfRGB);
         multcolor(bsdfRGB, rayIn -> pcol);
         addcolor(bsdfRGB, nd.tdiff);      
         photonRay(rayIn, &rayOut, PMAP_DIFFTRANS, bsdfRGB);
         transmitted = 1;
      }

      else if ((xi -= prSpecSD) <= 0) {
         /* Non-diffuse ("specular") SDF reflection */
         if ((err = SDsampBSDF(
            &bsdfVal, nd.vray, pmapRandom(scatterState), 
            SDsampSp | SDsampR, nd.sd
         )))
            objerror(mat, USER, transSDError(err));
         
         ccy2rgb(&bsdfVal.spec, bsdfVal.cieY, bsdfRGB);
         photonRay(rayIn, &rayOut, PMAP_SPECREFL, bsdfRGB);
      }
      
      else {
         /* Non-diffuse ("specular") SDF transmission */
         if ((err = SDsampBSDF(
            &bsdfVal, nd.vray, pmapRandom(scatterState), 
            SDsampSp | SDsampT, nd.sd
         )))
            objerror(mat, USER, transSDError(err));

         /* Apply pattern to spectral component */
         ccy2rgb(&bsdfVal.spec, bsdfVal.cieY, bsdfRGB);
         multcolor(bsdfRGB, rayIn -> pcol);
         photonRay(rayIn, &rayOut, PMAP_SPECTRANS, bsdfRGB);
         transmitted = 1;
      }      
      
      /* Xform outgoing dir to world coords */
      if ((err = SDmapDir(rayOut.rdir, nd.fromloc, nd.vray))) {
         objerror(mat, USER, transSDError(err));
         return 0;
      }
   }
      
   /* Clean up */
   SDfreeCache(nd.sd);

   /* Offset outgoing photon origin by thickness to bypass proxy geometry */
   if (transmitted && nd.thick != 0)
      VSUM(rayOut.rorg, rayOut.rorg, rayIn -> ron, -nd.thick);

   tracePhoton(&rayOut);
   return 0;
}



static int lightPhotonScatter (OBJREC* mat, RAY* ray)
/* Light sources doan' reflect, mang */
{
   return 0;
}



void initPhotonScatterFuncs ()
/* Init photonScatter[] dispatch table */
{
   int i;

   /* Catch-all for inconsistencies */
   for (i = 0; i < NUMOTYPE; i++) 
      photonScatter [i] = o_default;

   photonScatter [MAT_LIGHT] = photonScatter [MAT_ILLUM] =
      photonScatter [MAT_GLOW] = photonScatter [MAT_SPOT] = 
         lightPhotonScatter;

   photonScatter [MAT_PLASTIC] = photonScatter [MAT_METAL] =
      photonScatter [MAT_TRANS] = normalPhotonScatter;
      
   photonScatter [MAT_PLASTIC2] = photonScatter [MAT_METAL2] =
      photonScatter [MAT_TRANS2] = anisoPhotonScatter;
      
   photonScatter [MAT_DIELECTRIC] = photonScatter [MAT_INTERFACE] = 
      dielectricPhotonScatter;

   photonScatter [MAT_MIST] = mistPhotonScatter;
   photonScatter [MAT_GLASS] = glassPhotonScatter;
   photonScatter [MAT_CLIP] = clipPhotonScatter;
   photonScatter [MAT_MIRROR] = mirrorPhotonScatter;
   photonScatter [MIX_FUNC] = mx_funcPhotonScatter;
   photonScatter [MIX_DATA] = mx_dataPhotonScatter;
   photonScatter [MIX_PICT]= mx_pdataPhotonScatter;

   photonScatter [PAT_BDATA] = photonScatter [PAT_CDATA] =
      photonScatter [PAT_BFUNC] = photonScatter [PAT_CFUNC] =
         photonScatter [PAT_CPICT] = photonScatter [TEX_FUNC] = 
            photonScatter [TEX_DATA] = pattexPhotonScatter;

   photonScatter [MOD_ALIAS] = aliasPhotonScatter;
   photonScatter [MAT_BRTDF] = brdfPhotonScatter;
   
   photonScatter [MAT_PFUNC] = photonScatter [MAT_MFUNC] =
      photonScatter [MAT_PDATA] = photonScatter [MAT_MDATA] = 
         photonScatter [MAT_TFUNC] = photonScatter [MAT_TDATA] = 
            brdf2PhotonScatter;

   photonScatter [MAT_BSDF] = photonScatter [MAT_ABSDF] = 
      bsdfPhotonScatter;
}
