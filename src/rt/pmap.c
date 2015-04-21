/* 
   ==================================================================
   Photon map main module

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: pmap.c,v 2.2 2015/04/21 19:16:51 greg Exp $
*/



#include "pmap.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmaprand.h"
#include "pmapio.h"
#include "pmapbias.h"
#include "pmapdiag.h"
#include "otypes.h"
#include <time.h>
#include <sys/stat.h>



extern char *octname;

static char PmapRevision [] = "$Revision: 2.2 $";



/* Photon map lookup functions per type */
void (*pmapLookup [NUM_PMAP_TYPES])(PhotonMap*, RAY*, COLOR) = {
   photonDensity, photonPreCompDensity, photonDensity, volumePhotonDensity,
   photonDensity, NULL
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
         if (!stat(pm -> fileName, &pmstat) && !stat(octname, &octstat) && 
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
            deletePhotons(pmaps [type]);
            free(pmaps [type]);
         }
         pmaps [type] = pm;
         
         /* Check for invalid density estimate bandwidth */                            
         if (pm -> maxGather > pm -> heapSize) {
            error(WARNING, "adjusting density estimate bandwidth");
            pm -> minGather = pm -> maxGather = pm -> heapSize;
         }
      }
}



void savePmaps (const PhotonMap **pmaps, int argc, char **argv)
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      if (pmaps [t])
         savePhotonMap(pmaps [t], pmaps [t] -> fileName, t, argc, argv);
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


     
static int photonParticipate (RAY *ray)
/* Trace photon through participating medium. Returns 1 if passed through,
   or 0 if absorbed and $*%&ed. Analogon to rayparticipate(). */
{
   int i;
   RREAL cosTheta, cosPhi, du, dv;
   const float cext = colorAvg(ray -> cext),
               albedo = colorAvg(ray -> albedo);
   FVECT u, v;
   COLOR cvext;

   /* Mean free distance until interaction with medium */
   ray -> rmax = -log(pmapRandom(mediumState)) / cext;
   
   while (!localhit(ray, &thescene)) {
      setcolor(cvext, exp(-ray -> rmax * ray -> cext [0]),
                      exp(-ray -> rmax * ray -> cext [1]),
                      exp(-ray -> rmax * ray -> cext [2]));
                      
      /* Modify ray color and normalise */
      multcolor(ray -> rcol, cvext);
      colorNorm(ray -> rcol);
      VCOPY(ray -> rorg, ray -> rop);
      
      if (albedo > FTINY)
         /* Add to volume photon map */
         if (ray -> rlvl > 0) addPhoton(volumePmap, ray);
         
      /* Absorbed? */
      if (pmapRandom(rouletteState) > albedo) return 0;
      
      /* Colour bleeding without attenuation (?) */
      multcolor(ray -> rcol, ray -> albedo);
      scalecolor(ray -> rcol, 1 / albedo);    
      
      /* Scatter photon */
      cosTheta = ray -> gecc <= FTINY ? 2 * pmapRandom(scatterState) - 1
                                      : 1 / (2 * ray -> gecc) * 
                                            (1 + ray -> gecc * ray -> gecc - 
                                               (1 - ray -> gecc * ray -> gecc) / 
                                               (1 - ray -> gecc + 2 * ray -> gecc *
                                                  pmapRandom(scatterState)));
                                                  
      cosPhi = cos(2 * PI * pmapRandom(scatterState));
      du = dv = sqrt(1 - cosTheta * cosTheta);   /* sin(theta) */
      du *= cosPhi;
      dv *= sqrt(1 - cosPhi * cosPhi);           /* sin(phi) */
      
      /* Get axes u & v perpendicular to photon direction */
      i = 0;
      do {
         v [0] = v [1] = v [2] = 0;
         v [i++] = 1;
         fcross(u, v, ray -> rdir);
      } while (normalize(u) < FTINY);
      fcross(v, ray -> rdir, u);
      
      for (i = 0; i < 3; i++)
         ray -> rdir [i] = du * u [i] + dv * v [i] + 
                           cosTheta * ray -> rdir [i];
      ray -> rlvl++;
      ray -> rmax = -log(pmapRandom(mediumState)) / cext;
   }  
    
   setcolor(cvext, exp(-ray -> rot * ray -> cext [0]),
                   exp(-ray -> rot * ray -> cext [1]),
                   exp(-ray -> rot * ray -> cext [2]));
                   
   /* Modify ray color and normalise */
   multcolor(ray -> rcol, cvext);
   colorNorm(ray -> rcol);
   
   /* Passed through medium */  
   return 1;
}



void tracePhoton (RAY *ray)
/* Follow photon as it bounces around... */
{
   long mod;
   OBJREC* mat;
 
   if (ray -> rlvl > photonMaxBounce) {
      error(WARNING, "runaway photon!");
      return;
   }
   
   if (colorAvg(ray -> cext) > FTINY && !photonParticipate(ray)) 
      return;
      
   if (localhit(ray, &thescene)) {
      mod = ray -> ro -> omod;
      
      if ((ray -> clipset && inset(ray -> clipset, mod)) || mod == OVOID) {
         /* Transfer ray if modifier is VOID or clipped within antimatta */
         RAY tray;
         photonRay(ray, &tray, PMAP_XFER, NULL);
         tracePhoton(&tray);
      }
      else {
         /* Scatter for modifier material */
         mat = objptr(mod);
         photonScatter [mat -> otype] (mat, ray);
      }
   }
}



static void preComputeGlobal (PhotonMap *pmap)
/* Precompute irradiance from global photons for final gathering using 
   the first finalGather * pmap -> heapSize photons in the heap. Returns 
   new heap with precomputed photons. */
{
   unsigned long i, nuHeapSize;
   unsigned j;
   Photon *nuHeap, *p;
   COLOR irrad;
   RAY ray;
   float nuMinPos [3], nuMaxPos [3];

   repComplete = nuHeapSize = finalGather * pmap -> heapSize;
   
   if (photonRepTime) {
      sprintf(errmsg, 
              "Precomputing irradiance for %ld global photons...\n", 
              nuHeapSize);
      eputs(errmsg);
      fflush(stderr);
   }
   
   p = nuHeap = (Photon*)malloc(nuHeapSize * sizeof(Photon));
   if (!nuHeap) 
      error(USER, "can't allocate photon heap");
      
   for (j = 0; j <= 2; j++) {
      nuMinPos [j] = FHUGE;
      nuMaxPos [j] = -FHUGE;
   }
   
   /* Record start time, baby */
   repStartTime = time(NULL);
   signal(SIGCONT, pmapPreCompReport);
   repProgress = 0;
   memcpy(nuHeap, pmap -> heap, nuHeapSize * sizeof(Photon));
   
   for (i = 0, p = nuHeap; i < nuHeapSize; i++, p++) {
      ray.ro = NULL;
      VCOPY(ray.rop, p -> pos);
      
      /* Update min and max positions & set ray normal */
      for (j = 0; j < 3; j++) {
         if (p -> pos [j] < nuMinPos [j]) nuMinPos [j] = p -> pos [j];
         if (p -> pos [j] > nuMaxPos [j]) nuMaxPos [j] = p -> pos [j];
         ray.ron [j] = p -> norm [j] / 127.0;
      }
      
      photonDensity(pmap, &ray, irrad);
      setPhotonFlux(p, irrad);
      repProgress++;
      
      if (photonRepTime > 0 && time(NULL) >= repLastTime + photonRepTime)
         pmapPreCompReport();
#ifndef BSD
      else signal(SIGCONT, pmapPreCompReport);
#endif
   }
   
   signal(SIGCONT, SIG_DFL);
   
   /* Replace & rebuild heap */
   free(pmap -> heap);
   pmap -> heap = nuHeap;
   pmap -> heapSize = pmap -> heapEnd = nuHeapSize;
   VCOPY(pmap -> minPos, nuMinPos);
   VCOPY(pmap -> maxPos, nuMaxPos);
   
   if (photonRepTime) {
      eputs("Rebuilding global photon heap...\n");
      fflush(stderr);
   }
   
   balancePhotons(pmap, NULL);
}



void distribPhotons (PhotonMap **pmaps)
{
   EmissionMap emap;
   char errmsg2 [128];
   unsigned t, srcIdx, passCnt = 0, prePassCnt = 0;
   double totalFlux = 0;
   PhotonMap *pm;
   
   for (t = 0; t < NUM_PMAP_TYPES && !photonMaps [t]; t++);
   if (t >= NUM_PMAP_TYPES)
      error(USER, "no photon maps defined");
      
   if (!nsources)
      error(USER, "no light sources");

   /* ===================================================================
    * INITIALISATION - Set up emission and scattering funcs
    * =================================================================== */
   emap.samples = NULL;
   emap.maxPartitions = MAXSPART;
   emap.partitions = (unsigned char*)malloc(emap.maxPartitions >> 1);
   if (!emap.partitions)
      error(INTERNAL, "can't allocate source partitions");
      
   /* Initialise all defined photon maps */
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      initPhotonMap(photonMaps [t], t);

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
   
   if (photonRepTime)
      eputs("\n");
   
   /* ===================================================================
    * FLUX INTEGRATION - Get total photon flux from light sources
    * =================================================================== */
   for (srcIdx = 0; srcIdx < nsources; srcIdx++) {         
      unsigned portCnt = 0;
      emap.src = source + srcIdx; 
      
      do {
         emap.port = emap.src -> sflags & SDISTANT ? photonPorts + portCnt 
                                                   : NULL;
         photonPartition [emap.src -> so -> otype] (&emap);
         
         if (photonRepTime) {
            sprintf(errmsg, "Integrating flux from source %s ", 
                    source [srcIdx].so -> oname);
                    
            if (emap.port) {
               sprintf(errmsg2, "via port %s ", 
                       photonPorts [portCnt].so -> oname);
               strcat(errmsg, errmsg2);
            }
            
            sprintf(errmsg2, "(%lu partitions)...\n", emap.numPartitions);
            strcat(errmsg, errmsg2);
            eputs(errmsg);
            fflush(stderr);
         }
         
         for (emap.partitionCnt = 0; emap.partitionCnt < emap.numPartitions;
              emap.partitionCnt++) {
            initPhotonEmission(&emap, pdfSamples);
            totalFlux += colorAvg(emap.partFlux);
         }
         
         portCnt++;
      } while (portCnt < numPhotonPorts);
   }

   if (totalFlux < FTINY)
      error(USER, "zero flux from light sources");

   /* Record start time and enable progress report signal handler */
   repStartTime = time(NULL);
   signal(SIGCONT, pmapDistribReport);
   repProgress = prePassCnt = 0;
   
   if (photonRepTime)
      eputs("\n");
   
   /* ===================================================================
    * 2-PASS PHOTON DISTRIBUTION
    * Pass 1 (pre):  emit fraction of target photon count
    * Pass 2 (main): based on outcome of pass 1, estimate remaining number
    *                of photons to emit to approximate target count
    * =================================================================== */
   do {
      double numEmit;
      
      if (!passCnt) {   
         /* INIT PASS 1 */
         /* Skip if no photons contributed after sufficient iterations; make
          * it clear to user which photon maps are missing so (s)he can
          * check the scene geometry and materials */
         if (++prePassCnt > maxPreDistrib) {
            sprintf(errmsg, "too many prepasses");

            for (t = 0; t < NUM_PMAP_TYPES; t++)
               if (photonMaps [t] && !photonMaps [t] -> heapEnd) {
                  sprintf(errmsg2, ", no %s photons stored", pmapName [t]);
                  strcat(errmsg, errmsg2);
               }
            
            error(USER, errmsg);
            break;
         }

         /* Num to emit is fraction of minimum target count */
         numEmit = FHUGE;
         
         for (t = 0; t < NUM_PMAP_TYPES; t++)
            if (photonMaps [t])
               numEmit = min(photonMaps [t] -> distribTarget, numEmit);
               
         numEmit *= preDistrib;
      }

      else {            
         /* INIT PASS 2 */
         /* Based on the outcome of the predistribution we can now estimate
          * how many more photons we have to emit for each photon map to
          * meet its respective target count. This value is clamped to 0 in
          * case the target has already been exceeded in the pass 1. Note
          * repProgress is the number of photons emitted thus far, while
          * heapEnd is the number of photons stored in each photon map. */
         double maxDistribRatio = 0;

         /* Set the distribution ratio for each map; this indicates how many
          * photons of each respective type are stored per emitted photon,
          * and is used as probability for storing a photon by addPhoton().
          * Since this biases the photon density, addPhoton() promotes the
          * flux of stored photons to compensate. */
         for (t = 0; t < NUM_PMAP_TYPES; t++)
            if ((pm = photonMaps [t])) {
               pm -> distribRatio = (double)pm -> distribTarget / 
                                    pm -> heapEnd - 1;

               /* Check if photon map "overflowed", i.e. exceeded its target
                * count in the prepass; correcting the photon flux via the
                * distribution ratio is no longer possible, as no more
                * photons of this type will be stored, so notify the user
                * rather than deliver incorrect results. 
                * In future we should handle this more intelligently by
                * using the photonFlux in each photon map to individually
                * correct the flux after distribution. */
               if (pm -> distribRatio <= FTINY) {
                  sprintf(errmsg, 
                          "%s photon map overflow in prepass, reduce -apD",
                          pmapName [t]);
                  error(INTERNAL, errmsg);
               }
                  
               maxDistribRatio = max(pm -> distribRatio, maxDistribRatio);
            }
         
         /* Normalise distribution ratios and calculate number of photons to
          * emit in main pass */
         for (t = 0; t < NUM_PMAP_TYPES; t++)
            if ((pm = photonMaps [t]))
               pm -> distribRatio /= maxDistribRatio;
               
         if ((numEmit = repProgress * maxDistribRatio) < FTINY)
            /* No photons left to distribute in main pass */
            break;
      }
      
      /* Set completion count for progress report */
      repComplete = numEmit + repProgress;                             
      
      /* PHOTON DISTRIBUTION LOOP */
      for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
         unsigned portCnt = 0;
         emap.src = source + srcIdx;
                  
         do {
            emap.port = emap.src -> sflags & SDISTANT ? photonPorts + portCnt 
                                                      : NULL;
            photonPartition [emap.src -> so -> otype] (&emap);
            
            if (photonRepTime) {
               if (!passCnt)
                  sprintf(errmsg, "PREPASS %d on source %s ", 
                          prePassCnt, source [srcIdx].so -> oname);
               else 
                  sprintf(errmsg, "MAIN PASS on source %s ",
                          source [srcIdx].so -> oname);
                       
               if (emap.port) {
                  sprintf(errmsg2, "via port %s ", 
                          photonPorts [portCnt].so -> oname);
                  strcat(errmsg, errmsg2);
               }
               
               sprintf(errmsg2, "(%lu partitions)...\n", emap.numPartitions);
               strcat(errmsg, errmsg2);
               eputs(errmsg);
               fflush(stderr);
            }
            
            for (emap.partitionCnt = 0; emap.partitionCnt < emap.numPartitions;
                 emap.partitionCnt++) {             
               double partNumEmit;
               unsigned long partEmitCnt;
               
               /* Get photon origin within current source partishunn and
                * build emission map */
               photonOrigin [emap.src -> so -> otype] (&emap);
               initPhotonEmission(&emap, pdfSamples);
               
               /* Number of photons to emit from ziss partishunn --
                * proportional to flux; photon ray weight and scalar flux
                * are uniform (the latter only varying in RGB). */
               partNumEmit = numEmit * colorAvg(emap.partFlux) / totalFlux;
               partEmitCnt = (unsigned long)partNumEmit;
               
               /* Probabilistically account for fractional photons */
               if (pmapRandom(cntState) < partNumEmit - partEmitCnt)
                  partEmitCnt++;

               /* Integer counter avoids FP rounding errors */
               while (partEmitCnt--) {
                  RAY photonRay;
                  
                  /* Emit photon based on PDF and trace through scene until
                   * absorbed/leaked */
                  emitPhoton(&emap, &photonRay);
                  tracePhoton(&photonRay);
                  
                  /* Record progress */
                  repProgress++;
                  
                  if (photonRepTime > 0 && 
                      time(NULL) >= repLastTime + photonRepTime)
                     pmapDistribReport();
                  #ifndef BSD
                     else signal(SIGCONT, pmapDistribReport);
                  #endif
               }
            }
                        
            portCnt++;
         } while (portCnt < numPhotonPorts);
      }
      
      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (photonMaps [t] && !photonMaps [t] -> heapEnd) {
            /* Double preDistrib in case a photon map is empty and redo
             * pass 1 --> possibility of infinite loop for pathological
             * scenes (e.g. absorbing materials) */
            preDistrib *= 2;
            break;
         }
      
      if (t >= NUM_PMAP_TYPES) {
         /* No empty photon maps found; now do pass 2 */
         passCnt++;
         if (photonRepTime)
            eputs("\n");
      }
   } while (passCnt < 2);

   /* ===================================================================
    * POST-DISTRIBUTION - Set photon flux and build kd-tree, etc.
    * =================================================================== */
   signal(SIGCONT, SIG_DFL);
   free(emap.samples);
   
   /* Set photon flux (repProgress is total num emitted) */
   totalFlux /= repProgress;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (photonMaps [t]) {
         if (photonRepTime) {
            sprintf(errmsg, "\nBuilding %s photon map...\n", pmapName [t]);
            eputs(errmsg);
            fflush(stderr);
         }
      
         balancePhotons(photonMaps [t], &totalFlux);
      }
      
   /* Precompute photon irradiance if necessary */
   if (preCompPmap)
      preComputeGlobal(preCompPmap);
}



void photonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon density estimate. Returns irradiance at ray -> rop. */
{
   unsigned i;
   PhotonSQNode *sq;
   float r;
   COLOR flux;
 
   setcolor(irrad, 0, 0, 0);

   if (!pmap -> maxGather) 
      return;
      
   /* Ignore sources */
   if (ray -> ro) 
      if (islight(objptr(ray -> ro -> omod) -> otype)) 
         return;
         
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
      
   if (pmap -> minGather == pmap -> maxGather) {
      /* No bias compensation. Just do a plain vanilla estimate */
      sq = pmap -> squeue + 1;
      
      /* Average radius between furthest two photons to improve accuracy */      
      r = max(sq -> dist, (sq + 1) -> dist);
      r = 0.25 * (pmap -> maxDist + r + 2 * sqrt(pmap -> maxDist * r));   
      
      /* Skip the extra photon */
      for (i = 1 ; i < pmap -> squeueEnd; i++, sq++) {
         getPhotonFlux(sq -> photon, flux);         
#ifdef PMAP_EPANECHNIKOV
         /* Apply Epanechnikov kernel to photon flux (dists are squared) */
         scalecolor(flux, 2 * (1 - sq -> dist / r));
#endif         
         addcolor(irrad, flux);
      }
      
      /* Divide by search area PI * r^2, 1 / PI required as ambient 
         normalisation factor */         
      scalecolor(irrad, 1 / (PI * PI * r)); 
      
      return;
   }
   else 
      /* Apply bias compensation to density estimate */
      biasComp(pmap, irrad);
}



void photonPreCompDensity (PhotonMap *pmap, RAY *r, COLOR irrad)
/* Returns precomputed photon density estimate at ray -> rop. */
{
   Photon *p;
   
   setcolor(irrad, 0, 0, 0);

   /* Ignore sources */
   if (r -> ro && islight(objptr(r -> ro -> omod) -> otype)) 
      return;
      
   if ((p = find1Photon(preCompPmap, r))) 
      getPhotonFlux(p, irrad);
}



void volumePhotonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon volume density estimate. Returns irradiance at ray -> rop. */
{
   unsigned i;
   PhotonSQNode *sq;
   float gecc2, r, ph;
   COLOR flux;

   setcolor(irrad, 0, 0, 0);
   
   if (!pmap -> maxGather) 
      return;
      
   pmap -> squeueEnd = 0;
   findPhotons(pmap, ray);
   
   /* Need at least 2 photons */
   if (pmap -> squeueEnd < 2) 
      return;
      
   if (pmap -> minGather == pmap -> maxGather) {
      /* No bias compensation. Just do a plain vanilla estimate */
      gecc2 = ray -> gecc * ray -> gecc;
      sq = pmap -> squeue + 1;
      
      /* Average radius between furthest two photons to improve accuracy */      
      r = max(sq -> dist, (sq + 1) -> dist);
      r = 0.25 * (pmap -> maxDist + r + 2 * sqrt(pmap -> maxDist * r));   
      
      /* Skip the extra photon */
      for (i = 1 ; i < pmap -> squeueEnd; i++, sq++) {
         /* Compute phase function for inscattering from photon */
         if (gecc2 <= FTINY) 
            ph = 1;
         else {
            ph = DOT(ray -> rdir, sq -> photon -> norm) / 127;
            ph = 1 + gecc2 - 2 * ray -> gecc * ph;
            ph = (1 - gecc2) / (ph * sqrt(ph));
         }
         
         getPhotonFlux(sq -> photon, flux);
         scalecolor(flux, ph);
         addcolor(irrad, flux);
      }
      
      /* Divide by search volume 4 / 3 * PI * r^3 and phase function
         normalization factor 1 / (4 * PI) */
      scalecolor(irrad, 3 / (16 * PI * PI * r * sqrt(r)));
      
      return;
   }
   
   else 
      /* Apply bias compensation to density estimate */
      volumeBiasComp(pmap, ray, irrad);
}
