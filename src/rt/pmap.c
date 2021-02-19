#ifndef lint
static const char RCSid[] = "$Id: pmap.c,v 2.18 2021/02/19 02:10:35 rschregle Exp $";
#endif


/* 
   ======================================================================
   Photon map main module

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme FARESYS") 
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ======================================================================
   
   $Id: pmap.c,v 2.18 2021/02/19 02:10:35 rschregle Exp $
*/


#include "pmap.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmaprand.h"
#include "pmapio.h"
#include "pmapbias.h"
#include "pmapdiag.h"
#include "otypes.h"
#include "otspecial.h"
#include <time.h>
#if NIX
   #include <sys/stat.h>
   #include <sys/mman.h>
   #include <sys/wait.h>
#endif


void savePmaps (const PhotonMap **pmaps, int argc, char **argv)
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      if (pmaps [t])
         savePhotonMap(pmaps [t], pmaps [t] -> fileName, argc, argv);
   }
}                   


     
static int photonParticipate (RAY *ray)
/* Trace photon through participating medium. Returns 1 if passed through,
   or 0 if absorbed and $*%&ed. Analogon to rayparticipate(). */
{
   int i;
   RREAL xi1, cosTheta, phi, du, dv;
   const float cext = colorAvg(ray -> cext),
               albedo = colorAvg(ray -> albedo),
               gecc = ray -> gecc, gecc2 = sqr(gecc);
   FVECT u, v;
   COLOR cvext;

   /* Mean free distance until interaction with medium */
   ray -> rmax = -log(pmapRandom(mediumState)) / cext;
   
   while (!localhit(ray, &thescene)) {
      if (!incube(&thescene, ray -> rop)) {
         /* Terminate photon if it has leaked from the scene */
#ifdef DEBUG_PMAP
         fprintf(stderr, 
                 "Volume photon leaked from scene at [%.3f %.3f %.3f]\n", 
                 ray -> rop [0], ray -> rop [1], ray -> rop [2]);
#endif                 
         return 0;
      }
         
      setcolor(cvext, exp(-ray -> rmax * ray -> cext [0]),
                      exp(-ray -> rmax * ray -> cext [1]),
                      exp(-ray -> rmax * ray -> cext [2]));
                      
      /* Modify ray color and normalise */
      multcolor(ray -> rcol, cvext);
      colorNorm(ray -> rcol);
      VCOPY(ray -> rorg, ray -> rop);
      
#if 0
      if (albedo > FTINY && ray -> rlvl > 0)
#else
      /* Store volume photons unconditionally in mist to also account for
         direct inscattering from sources */
      if (albedo > FTINY)
#endif
         /* Add to volume photon map */
         newPhoton(volumePmap, ray);
         
      /* Absorbed? */
      if (pmapRandom(rouletteState) > albedo) 
         return 0;
      
      /* Colour bleeding without attenuation (?) */
      multcolor(ray -> rcol, ray -> albedo);
      scalecolor(ray -> rcol, 1 / albedo);    
      
      /* Scatter photon */
      xi1 = pmapRandom(scatterState);
      cosTheta = ray -> gecc <= FTINY 
                    ? 2 * xi1 - 1
                    : 0.5 / gecc * 
                      (1 + gecc2 - sqr((1 - gecc2) / 
                                       (1 + gecc * (2 * xi1 - 1))));

      phi = 2 * PI * pmapRandom(scatterState);
      du = dv = sqrt(1 - sqr(cosTheta));     /* sin(theta) */
      du *= cos(phi);
      dv *= sin(phi);
      
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
   
   /* Passed through medium until intersecting local object */  
   setcolor(cvext, exp(-ray -> rot * ray -> cext [0]),
                   exp(-ray -> rot * ray -> cext [1]),
                   exp(-ray -> rot * ray -> cext [2]));
                   
   /* Modify ray color and normalise */
   multcolor(ray -> rcol, cvext);
   colorNorm(ray -> rcol);   

   return 1;
}



void tracePhoton (RAY *ray)
/* Follow photon as it bounces around... */
{
   long mod;
   OBJREC *mat, *port = NULL;
   
   if (!ray -> parent) {
      /* !!!  PHOTON PORT REJECTION SAMPLING HACK: get photon port for
       * !!!  primary ray from ray -> ro, then reset the latter to NULL so
       * !!!  as not to interfere with localhit() */
      port = ray -> ro;
      ray -> ro = NULL;
   }
 
   if (ray -> rlvl > photonMaxBounce) {
#ifdef PMAP_RUNAWAY_WARN   
      error(WARNING, "runaway photon!");
#endif      
      return;
   }
  
   if (colorAvg(ray -> cext) > FTINY && !photonParticipate(ray)) 
      return;

   if (localhit(ray, &thescene)) {
      mod = ray -> ro -> omod;

      /* XXX: Is port -> omod != mod sufficient here? Probably not... */
      if (
         port && ray -> ro != port && 
         findmaterial(port) != findmaterial(ray -> ro)
      ) {
         /* !!! PHOTON PORT REJECTION SAMPLING HACK !!!
          * Terminate photon if emitted from port without intersecting it or
          * its other associated surfaces or same material. 
          * This can happen when the port's partitions extend beyond its
          * actual geometry, e.g.  with polygons.  Since the total flux
          * relayed by the port is based on the (in this case) larger
          * partition area, it is overestimated; terminating these photons
          * constitutes rejection sampling and thereby compensates any bias
          * incurred by the overestimated flux.  */
#ifdef PMAP_PORTREJECT_WARN
         sprintf(errmsg, "photon outside port %s", ray -> ro -> oname);
         error(WARNING, errmsg);
#endif         
         return;
      }

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
/* Precompute irradiance from global photons for final gathering for   
   a random subset of finalGather * pmap -> numPhotons photons, and builds
   the photon map, discarding the original photons. */
/* !!! NOTE: PRECOMPUTATION WITH OOC CURRENTLY WITHOUT CACHE !!! */   
{
   unsigned long  i, numPreComp;
   unsigned       j;
   PhotonIdx      pIdx;
   Photon         photon;
   RAY            ray;
   PhotonMap      nuPmap;

   repComplete = numPreComp = finalGather * pmap -> numPhotons;
   
   if (verbose) {
      sprintf(errmsg, 
              "\nPrecomputing irradiance for %ld global photons\n", 
              numPreComp);
      eputs(errmsg);
#if NIX      
      fflush(stderr);
#endif      
   }
   
   /* Copy photon map for precomputed photons */
   memcpy(&nuPmap, pmap, sizeof(PhotonMap));

   /* Zero counters, init new heap and extents */   
   nuPmap.numPhotons = 0;   
   initPhotonHeap(&nuPmap);
   
   for (j = 0; j < 3; j++) {   
      nuPmap.minPos [j] = FHUGE;
      nuPmap.maxPos [j] = -FHUGE;
   }

   /* Record start time, baby */
   repStartTime = time(NULL);
#ifdef SIGCONT
   signal(SIGCONT, pmapPreCompReport);
#endif
   repProgress = 0;
   
   photonRay(NULL, &ray, PRIMARY, NULL);
   ray.ro = NULL;
   
   for (i = 0; i < numPreComp; i++) {
      /* Get random photon from stratified distribution in source heap to
       * avoid duplicates and clustering */
      pIdx = firstPhoton(pmap) + 
             (unsigned long)((i + pmapRandom(pmap -> randState)) / 
                             finalGather);
      getPhoton(pmap, pIdx, &photon);
      
      /* Init dummy photon ray with intersection at photon position */
      VCOPY(ray.rop, photon.pos);
      for (j = 0; j < 3; j++)
         ray.ron [j] = photon.norm [j] / 127.0;
      
      /* Get density estimate at photon position */
      photonDensity(pmap, &ray, ray.rcol);
                  
      /* Append photon to new heap from ray */
      newPhoton(&nuPmap, &ray);
      
      /* Update progress */
      repProgress++;
      
      if (photonRepTime > 0 && time(NULL) >= repLastTime + photonRepTime)
         pmapPreCompReport();
#ifdef SIGCONT
      else signal(SIGCONT, pmapPreCompReport);
#endif
   }
   
   /* Flush heap */
   flushPhotonHeap(&nuPmap);
   
#ifdef SIGCONT   
   signal(SIGCONT, SIG_DFL);
#endif
   
   /* Trash original pmap, replace with precomputed one */
   deletePhotons(pmap);
   memcpy(pmap, &nuPmap, sizeof(PhotonMap));
   
   if (verbose) {
      eputs("\nRebuilding precomputed photon map\n");
#if NIX      
      fflush(stderr);
#endif      
   }

   /* Rebuild underlying data structure, destroying heap */   
   buildPhotonMap(pmap, NULL, NULL, 1);
}



typedef struct {
   unsigned long  numPhotons [NUM_PMAP_TYPES],
                  numEmitted, numComplete;
} PhotonCnt;



void distribPhotons (PhotonMap **pmaps, unsigned numProc)
{
   EmissionMap    emap;
   char           errmsg2 [128], shmFname [PMAP_TMPFNLEN];
   unsigned       t, srcIdx, proc;
   double         totalFlux = 0;
   int            shmFile, stat, pid;
   PhotonMap      *pm;
   PhotonCnt      *photonCnt;
   
   for (t = 0; t < NUM_PMAP_TYPES && !pmaps [t]; t++);
   
   if (t >= NUM_PMAP_TYPES)
      error(USER, "no photon maps defined in distribPhotons");
      
   if (!nsources)
      error(USER, "no light sources in distribPhotons");

   /* ===================================================================
    * INITIALISATION - Set up emission and scattering funcs
    * =================================================================== */
   emap.samples = NULL;
   emap.maxPartitions = MAXSPART;
   emap.partitions = (unsigned char*)malloc(emap.maxPartitions >> 1);
   if (!emap.partitions)
      error(INTERNAL, "can't allocate source partitions in distribPhotons");
      
   /* Initialise all defined photon maps */
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (pmaps [t]) {
         initPhotonMap(pmaps [t], t);
         /* Open photon heapfile */
         initPhotonHeap(pmaps [t]);
         /* Per-subprocess target count */
         pmaps [t] -> distribTarget /= numProc;
         
         if (!pmaps [t] -> distribTarget)
            error(INTERNAL, "no photons to distribute in distribPhotons");
      }

   initPhotonEmissionFuncs();
   initPhotonScatterFuncs();
   
   /* Get photon ports from modifier list */
   getPhotonPorts(photonPortList);

   /* Get photon sensor modifiers */
   getPhotonSensors(photonSensorList);
   
#if NIX
   /* Set up shared mem for photon counters (zeroed by ftruncate) */
   strcpy(shmFname, PMAP_TMPFNAME);
   shmFile = mkstemp(shmFname);

   if (shmFile < 0 || ftruncate(shmFile, sizeof(*photonCnt)) < 0)
      error(SYSTEM, "failed shared mem init in distribPhotons");

   photonCnt = mmap(NULL, sizeof(*photonCnt), PROT_READ | PROT_WRITE, 
                    MAP_SHARED, shmFile, 0);
                     
   if (photonCnt == MAP_FAILED)
      error(SYSTEM, "failed mapping shared memory in distribPhotons"); 
#else
   /* Allocate photon counters statically on Windoze */
   if (!(photonCnt = malloc(sizeof(PhotonCnt))))
      error(SYSTEM, "failed trivial malloc in distribPhotons");
   photonCnt -> numEmitted = photonCnt -> numComplete = 0;      
#endif /* NIX */

   if (verbose) {
      sprintf(errmsg, "\nIntegrating flux from %d sources", nsources);
      
      if (photonPorts) {
         sprintf(errmsg2, " via %d ports", numPhotonPorts);
         strcat(errmsg, errmsg2);
      }
      
      strcat(errmsg, "\n");
      eputs(errmsg);
   }   
   
   /* ===================================================================
    * FLUX INTEGRATION - Get total photon flux from light sources
    * =================================================================== */
   for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
      unsigned portCnt = 0;
      emap.src = source + srcIdx; 
      
      do {  /* Need at least one iteration if no ports! */
         emap.port = emap.src -> sflags & SDISTANT ? photonPorts + portCnt 
                                                   : NULL;
         photonPartition [emap.src -> so -> otype] (&emap);
         
         if (verbose) {
            sprintf(errmsg, "\tIntegrating flux from source %s ", 
                    source [srcIdx].so -> oname);

            if (emap.port) {
               sprintf(errmsg2, "via port %s ", 
                       photonPorts [portCnt].so -> oname);
               strcat(errmsg, errmsg2);
            }

            sprintf(errmsg2, "(%lu partitions)\n", emap.numPartitions);
            strcat(errmsg, errmsg2);
            eputs(errmsg);
#if NIX            
            fflush(stderr);
#endif            
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
      
   /* Record start time for progress reports */
   repStartTime = time(NULL);

   if (verbose) {
      sprintf(errmsg, "\nPhoton distribution @ %d procs\n", numProc);
      eputs(errmsg);
   }

   /* MAIN LOOP */   
   for (proc = 0; proc < numProc; proc++) {
#if NIX          
      if (!(pid = fork())) {
         /* SUBPROCESS ENTERS HERE; open and mmapped files inherited */
#else
      if (1) {
         /* No subprocess under Windoze */
#endif
         /* Local photon counters for this subprocess */
         unsigned       passCnt = 0, prePassCnt = 0;
         unsigned long  lastNumPhotons [NUM_PMAP_TYPES];
         unsigned long  localNumEmitted = 0; /* Num photons emitted by this
                                                subprocess alone */
         
         /* Seed RNGs from PID for decorellated photon distribution */
         pmapSeed(randSeed + proc, partState);
         pmapSeed(randSeed + (proc + 1) % numProc, emitState);
         pmapSeed(randSeed + (proc + 2) % numProc, cntState);
         pmapSeed(randSeed + (proc + 3) % numProc, mediumState);
         pmapSeed(randSeed + (proc + 4) % numProc, scatterState);
         pmapSeed(randSeed + (proc + 5) % numProc, rouletteState);
               
#ifdef DEBUG_PMAP          
         /* Output child process PID after random delay to prevent corrupted
          * console output due to race condition */
         usleep(1e6 * pmapRandom(rouletteState));
         fprintf(stderr, "Proc %d: PID = %d "
                 "(waiting 10 sec to attach debugger...)\n", 
                 proc, getpid());
         /* Allow time for debugger to attach to child process */
         sleep(10);
#endif            

         for (t = 0; t < NUM_PMAP_TYPES; t++)
            lastNumPhotons [t] = 0;
            
         /* =============================================================
          * 2-PASS PHOTON DISTRIBUTION
          * Pass 1 (pre):  emit fraction of target photon count
          * Pass 2 (main): based on outcome of pass 1, estimate remaining 
          *                number of photons to emit to approximate target 
          *                count
          * ============================================================= */
         do {
            double numEmit;
            
            if (!passCnt) {   
               /* INIT PASS 1 */               
               /* Skip if no photons contributed after sufficient
                * iterations; make it clear to user which photon maps are
                * missing so (s)he can check geometry and materials */
               if (++prePassCnt > maxPreDistrib) {
                  sprintf(errmsg, "proc %d: too many prepasses", proc);

                  for (t = 0; t < NUM_PMAP_TYPES; t++)
                     if (pmaps [t] && !pmaps [t] -> numPhotons) {
                        sprintf(errmsg2, ", no %s photons stored", 
                                pmapName [t]);
                        strcat(errmsg, errmsg2);
                     }
                  
                  error(USER, errmsg);
                  break;
               }

               /* Num to emit is fraction of minimum target count */
               numEmit = FHUGE;
               
               for (t = 0; t < NUM_PMAP_TYPES; t++)
                  if (pmaps [t])
                     numEmit = min(pmaps [t] -> distribTarget, numEmit);
                     
               numEmit *= preDistrib;
            }
            else {            
               /* INIT PASS 2 */
               /* Based on the outcome of the predistribution we can now
                * estimate how many more photons we have to emit for each
                * photon map to meet its respective target count.  This
                * value is clamped to 0 in case the target has already been
                * exceeded in the pass 1. */
               double maxDistribRatio = 0;

               /* Set the distribution ratio for each map; this indicates
                * how many photons of each respective type are stored per
                * emitted photon, and is used as probability for storing a
                * photon by newPhoton().  Since this biases the photon
                * density, newPhoton() promotes the flux of stored photons
                * to compensate.  */
               for (t = 0; t < NUM_PMAP_TYPES; t++)
                  if ((pm = pmaps [t])) {
                     pm -> distribRatio = (double)pm -> distribTarget / 
                                          pm -> numPhotons - 1;

                     /* Check if photon map "overflowed", i.e. exceeded its
                      * target count in the prepass; correcting the photon
                      * flux via the distribution ratio is no longer
                      * possible, as no more photons of this type will be
                      * stored, so notify the user rather than deliver
                      * incorrect results.  In future we should handle this
                      * more intelligently by using the photonFlux in each
                      * photon map to individually correct the flux after
                      * distribution.  */
                     if (pm -> distribRatio <= FTINY) {
                        sprintf(errmsg, "%s photon map overflow in "
                                "prepass, reduce -apD", pmapName [t]);
                        error(INTERNAL, errmsg);
                     }
                        
                     maxDistribRatio = max(pm -> distribRatio, 
                                           maxDistribRatio);
                  }
               
               /* Normalise distribution ratios and calculate number of
                * photons to emit in main pass */
               for (t = 0; t < NUM_PMAP_TYPES; t++)
                  if ((pm = pmaps [t]))
                     pm -> distribRatio /= maxDistribRatio;
                     
               if ((numEmit = localNumEmitted * maxDistribRatio) < FTINY)
                  /* No photons left to distribute in main pass */
                  break;
            }

            /* Update shared completion counter for progreport by parent */
            photonCnt -> numComplete += numEmit;                             

            /* PHOTON DISTRIBUTION LOOP */
            for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
               unsigned portCnt = 0;
               emap.src = source + srcIdx;

               do {  /* Need at least one iteration if no ports! */
                  emap.port = emap.src -> sflags & SDISTANT 
                              ? photonPorts + portCnt : NULL;
                  photonPartition [emap.src -> so -> otype] (&emap);

                  if (verbose && !proc) {
                     /* Output from subproc 0 only to avoid race condition
                      * on console I/O */
                     if (!passCnt)
                        sprintf(errmsg, "\tPREPASS %d on source %s ",
                                prePassCnt, source [srcIdx].so -> oname);
                     else 
                        sprintf(errmsg, "\tMAIN PASS on source %s ",
                                source [srcIdx].so -> oname);

                     if (emap.port) {
                        sprintf(errmsg2, "via port %s ", 
                                photonPorts [portCnt].so -> oname);
                        strcat(errmsg, errmsg2);
                     }

                     sprintf(errmsg2, "(%lu partitions)\n",
                             emap.numPartitions);
                     strcat(errmsg, errmsg2);
                     eputs(errmsg);
#if NIX                     
                     fflush(stderr);
#endif                     
                  }
                  
                  for (emap.partitionCnt = 0; emap.partitionCnt < emap.numPartitions;
                       emap.partitionCnt++) {
                     double partNumEmit;
                     unsigned long partEmitCnt;
                     
                     /* Get photon origin within current source partishunn
                      * and build emission map */
                     photonOrigin [emap.src -> so -> otype] (&emap);
                     initPhotonEmission(&emap, pdfSamples);
                     
                     /* Number of photons to emit from ziss partishunn --
                      * proportional to flux; photon ray weight and scalar
                      * flux are uniform (latter only varying in RGB). */
                     partNumEmit = numEmit * colorAvg(emap.partFlux) / 
                                   totalFlux;
                     partEmitCnt = (unsigned long)partNumEmit;
                     
                     /* Probabilistically account for fractional photons */
                     if (pmapRandom(cntState) < partNumEmit - partEmitCnt)
                        partEmitCnt++;

                     /* Update local and shared (global) emission counter */
                     photonCnt -> numEmitted += partEmitCnt;
                     localNumEmitted += partEmitCnt;

                     /* Integer counter avoids FP rounding errors during
                      * iteration */
                     while (partEmitCnt--) {
                        RAY photonRay;
                        
                        /* Emit photon based on PDF and trace through scene
                         * until absorbed/leaked */
                        emitPhoton(&emap, &photonRay);
#if 1
                        if (emap.port)
                           /* !!!  PHOTON PORT REJECTION SAMPLING HACK: set
                            * !!!  photon port as fake hit object for
                            * !!!  primary ray to check for intersection in
                            * !!!  tracePhoton() */
                           photonRay.ro = emap.port -> so;
#endif
                        tracePhoton(&photonRay);
                     }                                          

                     /* Update shared global photon count for each pmap */
                     for (t = 0; t < NUM_PMAP_TYPES; t++)
                        if (pmaps [t]) {
                           photonCnt -> numPhotons [t] += 
                              pmaps [t] -> numPhotons - lastNumPhotons [t];
                           lastNumPhotons [t] = pmaps [t] -> numPhotons;
                        }
#if !NIX
                     /* Synchronous progress report on Windoze */
                     if (!proc && photonRepTime > 0 && 
                           time(NULL) >= repLastTime + photonRepTime) {
                        repEmitted = repProgress = photonCnt -> numEmitted;
                        repComplete = photonCnt -> numComplete;                           
                        pmapDistribReport();
                     }
#endif
                  }
                  
                  portCnt++;
               } while (portCnt < numPhotonPorts);
            }
            
            for (t = 0; t < NUM_PMAP_TYPES; t++)
               if (pmaps [t] && !pmaps [t] -> numPhotons) {
                  /* Double preDistrib in case a photon map is empty and
                   * redo pass 1 --> possibility of infinite loop for
                   * pathological scenes (e.g.  absorbing materials) */
                  preDistrib *= 2;
                  break;
               }
            
            if (t >= NUM_PMAP_TYPES)
               /* No empty photon maps found; now do pass 2 */
               passCnt++;
         } while (passCnt < 2);
         
         /* Flush heap buffa for every pmap one final time; 
          * avoids potential data corruption! */
         for (t = 0; t < NUM_PMAP_TYPES; t++)
            if (pmaps [t]) {
               flushPhotonHeap(pmaps [t]);
               /* Heap file closed automatically on exit 
                  fclose(pmaps [t] -> heap); */
#ifdef DEBUG_PMAP               
               sprintf(errmsg, "Proc %d: total %ld photons\n", proc,
                       pmaps [t] -> numPhotons);
               eputs(errmsg);
#endif               
            }
#if NIX
         /* Terminate subprocess */
         exit(0);
#endif
      }
      else if (pid < 0)
         error(SYSTEM, "failed to fork subprocess in distribPhotons");         
   }

#if NIX
   /* PARENT PROCESS CONTINUES HERE */
#ifdef SIGCONT
   /* Enable progress report signal handler */
   signal(SIGCONT, pmapDistribReport);
#endif   
   /* Wait for subprocesses complete while reporting progress */
   proc = numProc;
   while (proc) {
      while (waitpid(-1, &stat, WNOHANG) > 0) {
         /* Subprocess exited; check status */
         if (!WIFEXITED(stat) || WEXITSTATUS(stat))
            error(USER, "failed photon distribution");
      
         --proc;
      }
      
      /* Nod off for a bit and update progress  */
      sleep(1);

      /* Asynchronous progress report from shared subprocess counters */  
      repEmitted = repProgress = photonCnt -> numEmitted;
      repComplete = photonCnt -> numComplete;      

      repProgress = repComplete = 0;
      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if ((pm = pmaps [t])) {
            /* Get global photon count from shmem updated by subprocs */
            repProgress += pm -> numPhotons = photonCnt -> numPhotons [t];
            repComplete += pm -> distribTarget;
         }
      repComplete *= numProc;

      if (photonRepTime > 0 && time(NULL) >= repLastTime + photonRepTime)
         pmapDistribReport();
#ifdef SIGCONT
      else signal(SIGCONT, pmapDistribReport);
#endif
   }
#endif /* NIX */

   /* ===================================================================
    * POST-DISTRIBUTION - Set photon flux and build data struct for photon
    * storage, etc.
    * =================================================================== */
#ifdef SIGCONT    
   /* Reset signal handler */
   signal(SIGCONT, SIG_DFL);
#endif
   free(emap.samples);
   
   /* Set photon flux */
   totalFlux /= photonCnt -> numEmitted;
#if NIX   
   /* Photon counters no longer needed, unmap shared memory */
   munmap(photonCnt, sizeof(*photonCnt));
   close(shmFile);
   unlink(shmFname);
#else
   free(photonCnt);   
#endif      
   if (verbose)
      eputs("\n");
      
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (pmaps [t]) {
         if (verbose) {
            sprintf(errmsg, "Building %s photon map\n", pmapName [t]);
            eputs(errmsg);
#if NIX            
            fflush(stderr);
#endif            
         }
         
         /* Build underlying data structure; heap is destroyed */
         buildPhotonMap(pmaps [t], &totalFlux, NULL, numProc);
      }
      
   /* Precompute photon irradiance if necessary */
   if (preCompPmap) {
      if (verbose)
         eputs("\n");
      preComputeGlobal(preCompPmap);
   }      
   
   if (verbose)
      eputs("\n");
}
