#ifndef lint
static const char RCSid[] = "$Id: pmutil.c,v 2.6 2021/03/23 00:07:13 rschregle Exp $";
#endif

/* 
   ======================================================================
   Photon map utilities

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmutil.c,v 2.6 2021/03/23 00:07:13 rschregle Exp $
*/

#include "pmap.h"
#include "pmapio.h"
#include "pmapbias.h"
#include "otypes.h"
#include <sys/stat.h>


extern char *octname;


/* Photon map lookup functions per type */
void (*pmapLookup [NUM_PMAP_TYPES])(PhotonMap*, RAY*, COLOR) = {
   photonDensity, photonPreCompDensity, photonDensity, volumePhotonDensity,
   photonDensity, photonDensity
};




void colorNorm (COLOR c)
/* Normalise colour channels to average of 1 */
{
   const float avg = colorAvg(c);
   
   if (!avg) 
      return;
      
   c [0] /= avg;
   c [1] /= avg;
   c [2] /= avg;
}




void loadPmaps (PhotonMap **pmaps, const PhotonMapParams *parm)
{
   unsigned t;
   struct stat octstat, pmstat;
   PhotonMap *pm;
   PhotonMapType type;

   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (setPmapParam(&pm, parm + t)) {
         /* Check if photon map newer than octree */
         if (pm -> fileName && octname &&
             !stat(pm -> fileName, &pmstat) && !stat(octname, &octstat) && 
             octstat.st_mtime > pmstat.st_mtime) {
            sprintf(errmsg, "photon map in file %s may be stale", 
                    pm -> fileName);
            error(USER, errmsg);
         }

         /* Load photon map from file and get its type */
         if ((type = loadPhotonMap(pm, pm -> fileName)) == PMAP_TYPE_NONE)
            error(USER, "failed loading photon map");
            
         /* Assign to appropriate photon map type (deleting previously
          * loaded photon map of same type if necessary) */
         if (pmaps [type]) {
            sprintf(errmsg, "multiple %s photon maps, dropping previous",
                    pmapName [type]);
            error(WARNING, errmsg);
            deletePhotons(pmaps [type]);
            free(pmaps [type]);
         }
         pmaps [type] = pm;

         /* Check for valid density estimate bandwidths */
         if ((pm -> minGather > 1 || pm -> maxGather > 1) &&
             (type == PMAP_TYPE_PRECOMP)) {
            /* Force bwidth to 1 for precomputed pmap */
            error(WARNING, "ignoring bandwidth for precomp photon map");
            pm -> minGather = pm -> maxGather = 1;
         }

         if ((pm -> maxGather > pm -> minGather) && 
             (type == PMAP_TYPE_VOLUME)) {            
            /* Biascomp for volume pmaps (see volumePhotonDensity() below) 
               is considered redundant, and there's probably no point in 
               recovering by using the lower bandwidth, since it's probably
               not what the user wants, so bail out. */
            sprintf(errmsg, 
                    "bias compensation is not available with %s photon maps",
                    pmapName [type]);
            error(USER, errmsg);
         }

         if (pm -> maxGather > pm -> numPhotons) {
            /* Clamp lookup bandwidth to total number of photons (minus one,
               since density estimate gets extra photon to obtain averaged 
               radius) */
            sprintf(
               errmsg, "clamping density estimate bandwidth to %ld",
               pm -> numPhotons
            );
            error(WARNING, errmsg);
            pm -> minGather = pm -> maxGather = pm -> numPhotons - 1;
         }
      }
}


   
void cleanUpPmaps (PhotonMap **pmaps)
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      if (pmaps [t]) {
         deletePhotons(pmaps [t]);
         free(pmaps [t]);
      }
   }
}




void photonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon density estimate. Returns irradiance at ray -> rop. */
{
   unsigned                      i;
   float                         r2;
   COLOR                         flux;
   Photon                        *photon;
   const PhotonSearchQueueNode   *sqn;
 
   setcolor(irrad, 0, 0, 0);

   if (!pmap -> maxGather) 
      return;
      
   /* Ignore sources */
   if (ray -> ro && islight(objptr(ray -> ro -> omod) -> otype)) 
      return;
         
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

   if (pmap -> minGather == pmap -> maxGather) {
      /* No bias compensation. Just do a plain vanilla estimate */
      sqn = pmap -> squeue.node + 1;
      
      /* Average radius^2 between furthest two photons to improve accuracy */
      r2 = max(sqn -> dist2, (sqn + 1) -> dist2);
      r2 = 0.25 * (pmap -> maxDist2 + r2 + 2 * sqrt(pmap -> maxDist2 * r2));
      
      /* Skip the extra photon */
      for (i = 1 ; i < pmap -> squeue.tail; i++, sqn++) {
         photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
         getPhotonFlux(photon, flux);         
#ifdef PMAP_EPANECHNIKOV
         /* Apply Epanechnikov kernel to photon flux based on photon dist */
         scalecolor(flux, 2 * (1 - sqn -> dist2 / r2));
#endif   
         addcolor(irrad, flux);
      }
      
      /* Divide by search area PI * r^2, 1 / PI required as ambient 
         normalisation factor */         
      scalecolor(irrad, 1 / (PI * PI * r2)); 
      
      return;
   }
   else 
      /* Apply bias compensation to density estimate */
      biasComp(pmap, irrad);
}




void photonPreCompDensity (PhotonMap *pmap, RAY *r, COLOR irrad)
/* Returns precomputed photon density estimate at ray -> rop. */
{
   Photon p;
   
   setcolor(irrad, 0, 0, 0);

   /* Ignore sources */
   if (r -> ro && islight(objptr(r -> ro -> omod) -> otype)) 
      return;
      
   if (find1Photon(preCompPmap, r, &p))
      /* p contains a found photon, so get its irradiance, otherwise it
       * remains zero under the assumption all photons are too distant
       * to contribute significantly */ 
      getPhotonFlux(&p, irrad);
}




void volumePhotonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon volume density estimate. Returns irradiance at ray -> rop. */
{
   unsigned                      i;
   float                         r2, gecc2, ph;
   COLOR                         flux;
   Photon                        *photon;
   const PhotonSearchQueueNode   *sqn;

   setcolor(irrad, 0, 0, 0);
   
   if (!pmap -> maxGather) 
      return;
      
   findPhotons(pmap, ray);
   
   /* Need at least 2 photons */
   if (pmap -> squeue.tail < 2) 
      return;

#if 0      
   /* Volume biascomp disabled (probably redundant) */
   if (pmap -> minGather == pmap -> maxGather)
#endif   
   {
      /* No bias compensation. Just do a plain vanilla estimate */
      gecc2 = ray -> gecc * ray -> gecc;
      sqn = pmap -> squeue.node + 1;
      
      /* Average radius^2 between furthest two photons to improve accuracy */      
      r2 = max(sqn -> dist2, (sqn + 1) -> dist2);
      r2 = 0.25 * (pmap -> maxDist2 + r2 + 2 * sqrt(pmap -> maxDist2 * r2));
      
      /* Skip the extra photon */
      for (i = 1; i < pmap -> squeue.tail; i++, sqn++) {
         photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
         
         /* Compute phase function for inscattering from photon */
         if (gecc2 <= FTINY) 
            ph = 1;
         else {
            ph = DOT(ray -> rdir, photon -> norm) / 127;
            ph = 1 + gecc2 - 2 * ray -> gecc * ph;
            ph = (1 - gecc2) / (ph * sqrt(ph));
         }
         
         getPhotonFlux(photon, flux);
         scalecolor(flux, ph);
         addcolor(irrad, flux);
      }
      
      /* Divide by search volume 4 / 3 * PI * r^3 and phase function
         normalization factor 1 / (4 * PI) */
      scalecolor(irrad, 3 / (16 * PI * PI * r2 * sqrt(r2)));
      return;
   }
#if 0
   else 
      /* Apply bias compensation to density estimate */
      volumeBiasComp(pmap, ray, irrad);
#endif      
}
