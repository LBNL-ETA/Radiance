#ifndef lint
static const char RCSid[] = "$Id: pmapcontrib.c,v 2.9 2015/08/18 18:45:55 greg Exp $";
#endif
/* 
   ==================================================================
   Photon map support for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapcontrib.c,v 2.9 2015/08/18 18:45:55 greg Exp $
*/


#include "pmapcontrib.h"
#include "pmap.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmaprand.h"
#include "pmapio.h"
#include "pmapdiag.h"
#include "rcontrib.h"
#include "otypes.h"



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
   const PhotonPrimary *primary = pmap -> primary;
   OBJREC *srcMod;
   unsigned long i, found = 0;
   
   /* Make sure at least one of the modifiers is actually in the pmap,
    * otherwise findPhotons() winds up in an infinite loop! */
   for (i = pmap -> primarySize; i; --i, ++primary) {
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
         /* Adjust density estimate bandwidth if lower than modifier
          * count, otherwise contributions are missing */
         error(WARNING, "contrib density estimate bandwidth too low, "
                        "adjusting to modifier count");
         contribPmap -> maxGather = numSrcContrib;
      }
      
      /* Sanity check */
      checkPmapContribs(contribPmap, srcContrib);
   }
}



void photonContrib (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Sum up light source contributions from photons in pmap->srcContrib */
{
   unsigned       i;
   PhotonSQNode   *sq;
   float          r, invArea;
   RREAL          rayCoeff [3];
 
   setcolor(irrad, 0, 0, 0);
 
   if (!pmap -> maxGather) 
      return;
      
   /* Ignore sources */
   if (ray -> ro) 
      if (islight(objptr(ray -> ro -> omod) -> otype)) 
         return;

   /* Get cumulative path
    * coefficient up to photon lookup point */
   raycontrib(rayCoeff, ray, PRIMARY);

   /* Lookup photons */
   pmap -> squeueEnd = 0;
   findPhotons(pmap, ray);
   
   /* Need at least 2 photons */
   if (pmap -> squeueEnd < 2) {
      #ifdef PMAP_NONEFOUND
         sprintf(errmsg, "no photons found on %s at (%.3f, %.3f, %.3f)", 
                 ray -> ro ? ray -> ro -> oname : "<null>",
                 ray -> rop [0], ray -> rop [1], ray -> rop [2]);
         error(WARNING, errmsg);
      #endif
      
      return;
   }

   /* Average (squared) radius between furthest two photons to improve
    * accuracy and get inverse search area 1 / (PI * r^2), with extra
    * normalisation factor 1 / PI for ambient calculation */
   sq = pmap -> squeue + 1;
   r = max(sq -> dist, (sq + 1) -> dist);
   r = 0.25 * (pmap -> maxDist + r + 2 * sqrt(pmap -> maxDist * r));   
   invArea = 1 / (PI * PI * r);
   
   /* Skip the extra photon */
   for (i = 1 ; i < pmap -> squeueEnd; i++, sq++) {         
      COLOR flux;
      
      /* Get photon's contribution to density estimate */
      getPhotonFlux(sq -> photon, flux);
      scalecolor(flux, invArea);
#ifdef PMAP_EPANECHNIKOV
      /* Apply Epanechnikov kernel to photon flux (dists are squared) */
      scalecolor(flux, 2 * (1 - sq -> dist / r));
#endif
      addcolor(irrad, flux);
      
      if (pmap -> srcContrib) {
         const PhotonPrimary *primary = pmap -> primary + 
                                        sq -> photon -> primary;
	 const SRCREC *sp = &source[primary -> srcIdx];
         OBJREC *srcMod = findmaterial(sp -> so);
         MODCONT *srcContrib = (MODCONT*)lu_find(pmap -> srcContrib, 
                                                 srcMod -> oname) -> data;
         if (!srcContrib)
	    continue;

	 /* Photon's emitting light source has modifier whose
          * contributions are sought */
	 double srcBinReal;
         int srcBin;
         RAY srcRay;

         if (srcContrib -> binv -> type != NUM) {
            /* Use intersection function to set shadow ray parameters
	     * if it's not simply a constant
	     */
            rayorigin(&srcRay, SHADOW, NULL, NULL);
	    srcRay.rsrc = primary -> srcIdx;
            VCOPY(srcRay.rorg, primary -> pos);
	    decodedir(srcRay.rdir, primary -> dir);

            if (!(sp->sflags & SDISTANT ? sourcehit(&srcRay)
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
                flux [j] = ray -> rcol [j] ? flux [j] / ray -> rcol [j]
                                             : 0;
	 }
                     
         multcolor(flux, rayCoeff);
         addcolor(srcContrib -> cbin [srcBin], flux);
      }
   }
        
   return;
}



void distribPhotonContrib (PhotonMap* pm)
{
   EmissionMap emap;
   char errmsg2 [128];
   unsigned srcIdx;
   double *srcFlux;                 /* Emitted flux per light source */
   const double srcDistribTarget =  /* Target photon count per source */
      nsources ? (double)pm -> distribTarget / nsources : 0;   

   if (!pm)
      error(USER, "no photon map defined");
      
   if (!nsources)
      error(USER, "no light sources");
   
   /* Allocate photon flux per light source; this differs for every 
    * source as all sources contribute the same number of distributed
    * photons (srcDistribTarget), hence the number of photons emitted per
    * source does not correlate with its emitted flux. The resulting flux
    * per photon is therefore adjusted individually for each source. */
   if (!(srcFlux = calloc(nsources, sizeof(double))))
      error(SYSTEM, "cannot allocate source flux");

   /* ================================================================
    * INITIALISASHUNN - Set up emisshunn and scattering funcs
    * ================================================================ */
   emap.samples = NULL;
   emap.src = NULL;
   emap.maxPartitions = MAXSPART;
   emap.partitions = (unsigned char*)malloc(emap.maxPartitions >> 1);
   if (!emap.partitions)
      error(USER, "can't allocate source partitions");

   initPhotonMap(pm, PMAP_TYPE_CONTRIB);
   initPhotonEmissionFuncs();
   initPhotonScatterFuncs();
   
   /* Get photon ports if specified */
   if (ambincl == 1)
      getPhotonPorts();
      
   /* Get photon sensor modifiers */
   getPhotonSensors(photonSensorList);      
   
   /* Seed RNGs for photon distribution */
   pmapSeed(randSeed, partState);
   pmapSeed(randSeed, emitState);
   pmapSeed(randSeed, cntState);
   pmapSeed(randSeed, mediumState);
   pmapSeed(randSeed, scatterState);
   pmapSeed(randSeed, rouletteState);

   /* Record start time and enable progress report signal handler */
   repStartTime = time(NULL);
   #ifdef SIGCONT   
      signal(SIGCONT, pmapDistribReport);
   #endif   
   
   for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
      unsigned portCnt = 0, passCnt = 0, prePassCnt = 0;
      double srcNumEmit = 0;          /* # photons to emit from source */
      unsigned long srcNumDistrib = pm -> heapEnd; /* # photons stored */

      srcFlux [srcIdx] = repProgress = 0;
      emap.src = source + srcIdx;
      
      if (photonRepTime) 
         eputs("\n");
      
      /* =============================================================
       * FLUX INTEGRATION - Get total flux emitted from light source
       * ============================================================= */
      do {
         emap.port = emap.src -> sflags & SDISTANT 
                     ? photonPorts + portCnt : NULL;
         photonPartition [emap.src -> so -> otype] (&emap);
         
         if (photonRepTime) {
            sprintf(errmsg, "Integrating flux from source %s (mod %s) ", 
                    source [srcIdx].so -> oname, 
                    objptr(source [srcIdx].so -> omod) -> oname);
                    
            if (emap.port) {
               sprintf(errmsg2, "via port %s ", 
                       photonPorts [portCnt].so -> oname);
               strcat(errmsg, errmsg2);
            }
            
            sprintf(errmsg2, "(%lu partitions)...\n", 
                    emap.numPartitions);
            strcat(errmsg, errmsg2);
            eputs(errmsg);
            fflush(stderr);
         }
         
         for (emap.partitionCnt = 0; 
              emap.partitionCnt < emap.numPartitions; 
              emap.partitionCnt++) {
            initPhotonEmission(&emap, pdfSamples);
            srcFlux [srcIdx] += colorAvg(emap.partFlux);
         }
         
         portCnt++;
      } while (portCnt < numPhotonPorts);

      if (srcFlux [srcIdx] < FTINY) {
         sprintf(errmsg, "source %s has zero emission", 
                 source [srcIdx].so -> oname);
         error(WARNING, errmsg);
      }
      else {
         /* ==========================================================
          * 2-PASS PHOTON DISTRIBUTION
          * Pass 1 (pre):  emit fraction of target photon count
          * Pass 2 (main): based on outcome of pass 1, estimate
          *                remaining number of photons to emit to
          *                approximate target count
          * ========================================================== */
         do {
            if (!passCnt) {   
               /* INIT PASS 1 */
               if (++prePassCnt > maxPreDistrib) {
                  /* Warn if no photons contributed after sufficient
                   * iterations */
                  sprintf(errmsg, "too many prepasses, no photons "
                          "from source %s", source [srcIdx].so -> oname);
                  error(WARNING, errmsg);
                  break;
               }
               
               /* Num to emit is fraction of target count */
               srcNumEmit = preDistrib * srcDistribTarget;
            }

            else {            
               /* INIT PASS 2 */                           
               /* Based on the outcome of the predistribution we can now
                * figure out how many more photons we have to emit from
                * the current source to meet the target count,
                * srcDistribTarget. This value is clamped to 0 in case
                * the target has already been exceeded in pass 1.
                * srcNumEmit and srcNumDistrib is the number of photons
                * emitted and distributed (stored) from the current
                * source in pass 1, respectively. */
               srcNumDistrib = pm -> heapEnd - srcNumDistrib;
               srcNumEmit *= srcNumDistrib 
                             ? max(srcDistribTarget/srcNumDistrib, 1) - 1
                             : 0;

               if (!srcNumEmit)
                  /* No photons left to distribute in main pass */
                  break;
            }
         
            /* Set completion count for progress report */
            repComplete = srcNumEmit + repProgress;
            portCnt = 0;
                     
            do {
               emap.port = emap.src -> sflags & SDISTANT 
                           ? photonPorts + portCnt : NULL;
               photonPartition [emap.src -> so -> otype] (&emap);
               
               if (photonRepTime) {
                  if (!passCnt)
                     sprintf(errmsg, "PREPASS %d on source %s (mod %s) ",
                             prePassCnt, source [srcIdx].so -> oname,
                             objptr(source[srcIdx].so->omod) -> oname);
                  else 
                     sprintf(errmsg, "MAIN PASS on source %s (mod %s) ",
                             source [srcIdx].so -> oname,
                             objptr(source[srcIdx].so->omod) -> oname);
                          
                  if (emap.port) {
                     sprintf(errmsg2, "via port %s ", 
                             photonPorts [portCnt].so -> oname);
                     strcat(errmsg, errmsg2);
                  }
                  
                  sprintf(errmsg2, "(%lu partitions)...\n",
                          emap.numPartitions);
                  strcat(errmsg, errmsg2);
                  eputs(errmsg);
                  fflush(stderr);
               }
               
               for (emap.partitionCnt = 0; 
                    emap.partitionCnt < emap.numPartitions; 
                    emap.partitionCnt++) {
                  double partNumEmit;
                  unsigned long partEmitCnt;
                  
                  /* Get photon origin within current source partishunn
                   * and build emission map */
                  photonOrigin [emap.src -> so -> otype] (&emap);
                  initPhotonEmission(&emap, pdfSamples);
                  
                  /* Number of photons to emit from ziss partishunn;
                   * scale according to its normalised contribushunn to
                   * the emitted source flux */
                  partNumEmit = srcNumEmit * colorAvg(emap.partFlux) / 
                                srcFlux [srcIdx];
                  partEmitCnt = (unsigned long)partNumEmit;
                                       
                  /* Probabilistically account for fractional photons */
                  if (pmapRandom(cntState) < partNumEmit - partEmitCnt)
                     partEmitCnt++;
                     
                  /* Integer counter avoids FP rounding errors */
                  while (partEmitCnt--) {                   
                     RAY photonRay;
                  
                     /* Emit photon according to PDF (if any), allocate
                      * associated primary ray, and trace through scene
                      * until absorbed/leaked */
                     emitPhoton(&emap, &photonRay);
                     addPhotonPrimary(pm, &photonRay);
                     tracePhoton(&photonRay);
                     
                     /* Record progress */
                     repProgress++;
                     
                     if (photonRepTime > 0 && 
                         time(NULL) >= repLastTime + photonRepTime)
                        pmapDistribReport();
                     #ifdef SIGCONT
                        else signal(SIGCONT, pmapDistribReport);
                     #endif
                  }
               }
                           
               portCnt++;
            } while (portCnt < numPhotonPorts);

            if (pm -> heapEnd == srcNumDistrib) 
               /* Double preDistrib in case no photons were stored 
                * for this source and redo pass 1 */
               preDistrib *= 2;
            else {
               /* Now do pass 2 */
               passCnt++;
               if (photonRepTime)
                  eputs("\n");
            }
         } while (passCnt < 2);
         
         /* Flux per photon emitted from this source; repProgress is the
          * number of emitted photons after both passes */
         srcFlux [srcIdx] = repProgress ? srcFlux [srcIdx] / repProgress 
                                        : 0;
      }
   }

   /* ================================================================
    * POST-DISTRIBUTION - Set photon flux and build kd-tree, etc.
    * ================================================================ */
   #ifdef SIGCONT    
      signal(SIGCONT, SIG_DFL);
   #endif   
   free(emap.samples);

   if (!pm -> heapEnd)
      error(USER, "empty photon map");

   /* Check for valid primary photon rays */
   if (!pm -> primary)
      error(INTERNAL, "no primary rays in contribution photon map");
      
   if (pm -> primary [pm -> primaryEnd].srcIdx < 0)
      /* Last primary ray is unused, so decrement counter */
      pm -> primaryEnd--;
   
   if (photonRepTime) {
      eputs("\nBuilding contrib photon heap...\n");
      fflush(stderr);
   }
      
   balancePhotons(pm, srcFlux);
}
