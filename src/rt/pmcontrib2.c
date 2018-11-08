#ifndef lint
static const char RCSid[] = "$Id: pmcontrib2.c,v 2.5 2018/11/08 00:54:07 greg Exp $";
#endif

/* 
   ======================================================================
   Photon map support for using light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmcontrib2.c,v 2.5 2018/11/08 00:54:07 greg Exp $
*/


#include "pmapcontrib.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmaprand.h"
#include "pmapio.h"
#include "pmapdiag.h"
#include "rcontrib.h"
#include "otypes.h"
#include "otspecial.h"


static void setPmapContribParams (PhotonMap *pmap, LUTAB *srcContrib)
/* Set parameters for light source contributions */
{
   /* Set light source modifier list and appropriate callback to extract
    * their contributions from the photon map */
   if (pmap) {
      pmap -> srcContrib = srcContrib;
      pmap -> lookup = photonContrib;
      /* Ensure we get all requested photon contribs during lookups */
      pmap -> gatherTolerance = 1.0;
   }
}



static void checkPmapContribs (const PhotonMap *pmap, LUTAB *srcContrib)
/* Check modifiers for light source contributions */
{
   const PhotonPrimary  *primary = pmap -> primaries;
   PhotonPrimaryIdx     i, found = 0;
   OBJREC *srcMod;
   
   /* Make sure at least one of the modifiers is actually in the pmap,
    * otherwise findPhotons() winds up in an infinite loop! */
   for (i = pmap -> numPrimary; i; --i, ++primary) {
      if (primary -> srcIdx < 0 || primary -> srcIdx >= nsources)
         error(INTERNAL, "invalid light source index in photon map");
         
      srcMod = findmaterial(source [primary -> srcIdx].so);
      if ((MODCONT*)lu_find(srcContrib, srcMod -> oname) -> data)
         ++found;
   }
   
   if (!found)
      error(USER, "modifiers not in photon map");
}
 
   

void initPmapContrib (LUTAB *srcContrib, unsigned numSrcContrib)
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (photonMaps [t] && t != PMAP_TYPE_CONTRIB) {
         sprintf(errmsg, "%s photon map does not support contributions",
                 pmapName [t]);
         error(USER, errmsg);
      }
   
   /* Get params */
   setPmapContribParams(contribPmap, srcContrib);
   
   if (contribPhotonMapping) {
      if (contribPmap -> maxGather < numSrcContrib) {
#if 0          
         /* Adjust density estimate bandwidth if lower than modifier
          * count, otherwise contributions are missing */
         error(WARNING, "contrib density estimate bandwidth too low, "
                        "adjusting to modifier count");
         contribPmap -> maxGather = numSrcContrib;
#else
         /* Warn if density estimate bandwidth is lower than modifier
          * count, since some contributions will be missing */
         error(WARNING, "photon density estimate bandwidth too low,"
                        " contributions may be underestimated");
#endif
      }
      
      /* Sanity check */
      checkPmapContribs(contribPmap, srcContrib);
   }
}



void photonContrib (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Sum up light source contributions from photons in pmap->srcContrib */
{
   unsigned                i;
   PhotonSearchQueueNode   *sqn;
   float                   r2, invArea;
   RREAL                   rayCoeff [3];
   Photon                  *photon;
   static char             warnPos = 1, warnDir = 1;
 
   setcolor(irrad, 0, 0, 0);
 
   if (!pmap -> maxGather) 
      return;
      
   /* Ignore sources */
   if (ray -> ro && islight(objptr(ray -> ro -> omod) -> otype)) 
      return;

   /* Get cumulative path coefficient up to photon lookup point */
   raycontrib(rayCoeff, ray, PRIMARY);

   /* Lookup photons */
   pmap -> squeue.tail = 0;
   findPhotons(pmap, ray);
   
   /* Need at least 2 photons */
   if (pmap -> squeue.tail < 2) {
#ifdef PMAP_NONEFOUND
      sprintf(errmsg, "no photons found on %s at (%.3f, %.3f, %.3f)", 
              ray -> ro ? ray -> ro -> oname : "<null>",
              ray -> rop [0], ray -> rop [1], ray -> rop [2]);
      error(WARNING, errmsg);
#endif
      
      return;
   }

   /* Average radius^2 between furthest two photons to improve accuracy and
    * get inverse search area 1 / (PI * r^2), with extra normalisation
    * factor 1 / PI for ambient calculation */
   sqn = pmap -> squeue.node + 1;
   r2 = max(sqn -> dist2, (sqn + 1) -> dist2);
   r2 = 0.25 * (pmap -> maxDist2 + r2 + 2 * sqrt(pmap -> maxDist2 * r2));
   invArea = 1 / (PI * PI * r2);
   
   /* Skip the extra photon */
   for (i = 1 ; i < pmap -> squeue.tail; i++, sqn++) {         
      COLOR flux;
      
      /* Get photon's contribution to density estimate */
      photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
      getPhotonFlux(photon, flux);
      scalecolor(flux, invArea);      
#ifdef PMAP_EPANECHNIKOV
      /* Apply Epanechnikov kernel to photon flux based on photon distance */
      scalecolor(flux, 2 * (1 - sqn -> dist2 / r2));
#endif
      addcolor(irrad, flux);
      
      if (pmap -> srcContrib) {      
         const PhotonPrimary  *primary = pmap -> primaries + 
                                         photon -> primary;
         const SRCREC         *sp = &source [primary -> srcIdx];
         OBJREC               *srcMod = findmaterial(sp -> so);
         MODCONT *srcContrib = (MODCONT*)lu_find(pmap -> srcContrib, 
                                                 srcMod -> oname) -> data;
         double   srcBinReal;
         int      srcBin;
         RAY      srcRay;                                                 
                                                 
         if (!srcContrib)
            continue;

         /* Photon's emitting light source has modifier whose contributions
          * are sought */
         if (srcContrib -> binv -> type != NUM) {
            /* Use intersection function to set shadow ray parameters if
             * it's not simply a constant */
            rayorigin(&srcRay, SHADOW, NULL, NULL);
            srcRay.rsrc = primary -> srcIdx;
#ifdef PMAP_PRIMARYPOS
            VCOPY(srcRay.rorg, primary -> pos);
#else
            /* No primary hitpoints; set dummy ray origin and warn once */
            srcRay.rorg [0] = srcRay.rorg [1] = srcRay.rorg [2] = 0;
            if (warnPos) {
               error(WARNING, 
                     "no photon primary hitpoints for bin evaluation; "
                     "using dummy (0,0,0)! Recompile with -DPMAP_CBDM.");
               warnPos = 0;
            }
#endif           
#ifdef PMAP_PRIMARYDIR
            decodedir(srcRay.rdir, primary -> dir);
#else
            /* No primary incident direction; set dummy and warn once */
            if (warnDir) {
               error(WARNING, 
                     "no photon primary directions for bin evaluation; "
                     "using dummy (0,0,0)! Recompile with -DPMAP_CBDM.");
               warnDir = 0;
            }
            srcRay.rdir [0] = srcRay.rdir [1] = srcRay.rdir [2] = 0;
#endif

            if (!(sp->sflags & SDISTANT 
                  ? sourcehit(&srcRay)
                  : (*ofun[sp -> so -> otype].funp)(sp -> so, &srcRay)))
               continue;		/* XXX shouldn't happen! */

            worldfunc(RCCONTEXT, &srcRay);
            set_eparams((char *)srcContrib -> params);
         }

         if ((srcBinReal = evalue(srcContrib -> binv)) < -.5)
            continue;		/* silently ignore negative bins */
  
         if ((srcBin = srcBinReal + .5) >= srcContrib -> nbins) {
            error(WARNING, "bad bin number (ignored)");
            continue;
         }
            
         if (!contrib) {
            /* Ray coefficient mode; normalise by light source radiance
             * after applying distrib pattern */
            int j;
            
            raytexture(ray, srcMod -> omod);
            setcolor(ray -> rcol, srcMod -> oargs.farg [0], 
                     srcMod -> oargs.farg [1], srcMod -> oargs.farg [2]);
            multcolor(ray -> rcol, ray -> pcol);
            for (j = 0; j < 3; j++)
               flux [j] = ray -> rcol [j] ? flux [j] / ray -> rcol [j] : 0;
         }
                     
         multcolor(flux, rayCoeff);
         addcolor(srcContrib -> cbin [srcBin], flux);
      }
   }
        
   return;
}
