#ifndef lint
static const char RCSid[] = "$Id: pmap.c,v 2.11 2016/05/17 17:39:47 rschregle Exp $";
#endif

/* 
   ======================================================================
   Photon map main module

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmap.c,v 2.11 2016/05/17 17:39:47 rschregle Exp $
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
#include <sys/mman.h>
#include <sys/wait.h>

#define PMAP_REV  "$Revision: 2.11 $"


extern char *octname;



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
            deletePhotons(pmaps [type]);
            free(pmaps [type]);
         }
         pmaps [type] = pm;
         
         /* Check for invalid density estimate bandwidth */                            
         if (pm -> maxGather > pm -> numPhotons) {
            error(WARNING, "adjusting density estimate bandwidth");
            pm -> minGather = pm -> maxGather = pm -> numPhotons;
         }
      }
}



void savePmaps (const PhotonMap **pmaps, int argc, char **argv)
{
   unsigned t;
   
   for (t = 0; t < NUM_PMAP_TYPES; t++) {
      if (pmaps [t])
         savePhotonMap(pmaps [t], pmaps [t] -> fileName, argc, argv);
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
      
      if (albedo > FTINY && ray -> rlvl > 0)
         /* Add to volume photon map */
         newPhoton(volumePmap, ray);
         
      /* Absorbed? */
      if (pmapRandom(rouletteState) > albedo) 
         return 0;
      
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
#ifdef PMAP_RUNAWAY_WARN   
      error(WARNING, "runaway photon!");
#endif      
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
   
   if (photonRepTime) {
      sprintf(errmsg, "Precomputing irradiance for %ld global photons...\n", 
              numPreComp);
      eputs(errmsg);
      fflush(stderr);
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
       * avoid duplicates and clutering */
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
   
   if (photonRepTime) {
      eputs("Rebuilding precomputed photon map...\n");
      fflush(stderr);
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
   char           errmsg2 [128], shmFname [255];
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
      }

   initPhotonEmissionFuncs();
   initPhotonScatterFuncs();
   
   /* Get photon ports if specified */
   if (ambincl == 1) 
      getPhotonPorts();

   /* Get photon sensor modifiers */
   getPhotonSensors(photonSensorList);
   
   /* Set up shared mem for photon counters (zeroed by ftruncate) */
#if 0   
   snprintf(shmFname, 255, PMAP_SHMFNAME, getpid());
   shmFile = shm_open(shmFname, O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
#else
   strcpy(shmFname, PMAP_SHMFNAME);
   shmFile = mkstemp(shmFname);
#endif      

   if (shmFile < 0)
      error(SYSTEM, "failed opening shared memory file in distribPhotons");

   if (ftruncate(shmFile, sizeof(*photonCnt)) < 0)
      error(SYSTEM, "failed setting shared memory size in distribPhotons");

   photonCnt = mmap(NULL, sizeof(*photonCnt), PROT_READ | PROT_WRITE, 
                    MAP_SHARED, shmFile, 0);
                     
   if (photonCnt == MAP_FAILED)
      error(SYSTEM, "failed mapping shared memory in distribPhotons");

   if (photonRepTime)
      eputs("\n");
   
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

   /* MAIN LOOP */   
   for (proc = 0; proc < numProc; proc++) {
      if (!(pid = fork())) {
         /* SUBPROCESS ENTERS HERE.
            All opened and memory mapped files are inherited */
         unsigned       passCnt = 0, prePassCnt = 0;
         unsigned long  lastNumPhotons [NUM_PMAP_TYPES];
         unsigned long  localNumEmitted = 0; /* Num photons emitted by this
                                                subprocess alone */
         
         /* Seed RNGs from PID for decorellated photon distribution */
         pmapSeed(randSeed + proc, partState);
         pmapSeed(randSeed + proc, emitState);
         pmapSeed(randSeed + proc, cntState);
         pmapSeed(randSeed + proc, mediumState);
         pmapSeed(randSeed + proc, scatterState);
         pmapSeed(randSeed + proc, rouletteState);
                  
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
                  sprintf(errmsg, 
                          "proc %d, source %s: too many prepasses",
                          proc, source [srcIdx].so -> oname);               

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

            /* Update shared completion counter for prog.report by parent */
            photonCnt -> numComplete += numEmit;                             

            /* PHOTON DISTRIBUTION LOOP */
            for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
               unsigned portCnt = 0;
               emap.src = source + srcIdx;

               do {  /* Need at least one iteration if no ports! */
                  emap.port = emap.src -> sflags & SDISTANT 
                              ? photonPorts + portCnt : NULL;
                  photonPartition [emap.src -> so -> otype] (&emap);

                  if (photonRepTime && !proc) {
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
                     
                     sprintf(errmsg2, "(%lu partitions)\n",
                             emap.numPartitions);
                     strcat(errmsg, errmsg2);
                     eputs(errmsg);
                     fflush(stderr);
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
                      * flux are uniform (the latter only varying in RGB).
                      * */ 
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
                        tracePhoton(&photonRay);
                     }                                          
                     
                     /* Update shared global photon count for each pmap */
                     for (t = 0; t < NUM_PMAP_TYPES; t++)
                        if (pmaps [t]) {
                           photonCnt -> numPhotons [t] += 
                              pmaps [t] -> numPhotons - lastNumPhotons [t];
                           lastNumPhotons [t] = pmaps [t] -> numPhotons;
                        }
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
            
            if (t >= NUM_PMAP_TYPES) {
               /* No empty photon maps found; now do pass 2 */
               passCnt++;
#if 0
               if (photonRepTime)
                  eputs("\n");
#endif
            }
         } while (passCnt < 2);
         
         /* Unmap shared photon counters */
#if 0         
         munmap(photonCnt, sizeof(*photonCnt));
         close(shmFile);
#endif
         
         /* Flush heap buffa for every pmap one final time; this is required
          * to prevent data corruption! */
         for (t = 0; t < NUM_PMAP_TYPES; t++)
            if (pmaps [t]) {
#if 0            
               eputs("Final flush\n");
#endif               
               flushPhotonHeap(pmaps [t]);
               fclose(pmaps [t] -> heap);
#ifdef DEBUG_PMAP               
               sprintf(errmsg, "Proc %d: total %ld photons\n", getpid(),
                       pmaps [t] -> numPhotons);
               eputs(errmsg);
#endif               
            }

         exit(0);
      }
      else if (pid < 0)
         error(SYSTEM, "failed to fork subprocess in distribPhotons");         
   }

   /* PARENT PROCESS CONTINUES HERE */
   /* Record start time and enable progress report signal handler */
   repStartTime = time(NULL);
#ifdef SIGCONT
   signal(SIGCONT, pmapDistribReport);
#endif
   
   if (photonRepTime)
      eputs("\n");
   
   /* Wait for subprocesses to complete while reporting progress */
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
      /* Update progress report from shared subprocess counters */
      repEmitted = repProgress = photonCnt -> numEmitted;
      repComplete = photonCnt -> numComplete;

      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if ((pm = pmaps [t])) {
#if 0         
            /* Get photon count from heapfile size for progress update */
            fseek(pm -> heap, 0, SEEK_END);
            pm -> numPhotons = ftell(pm -> heap) / sizeof(Photon); */
#else            
            /* Get global photon count from shmem updated by subprocs */
            pm -> numPhotons = photonCnt -> numPhotons [t];
#endif            
         }

      if (photonRepTime > 0 && time(NULL) >= repLastTime + photonRepTime)
         pmapDistribReport();
#ifdef SIGCONT
      else signal(SIGCONT, pmapDistribReport);
#endif
   }

   /* ===================================================================
    * POST-DISTRIBUTION - Set photon flux and build data struct for photon
    * storage, etc.
    * =================================================================== */
#ifdef SIGCONT    
   signal(SIGCONT, SIG_DFL);
#endif
   free(emap.samples);
   
   /* Set photon flux (repProgress is total num emitted) */
   totalFlux /= photonCnt -> numEmitted;
   
   /* Photon counters no longer needed, unmap shared memory */
   munmap(photonCnt, sizeof(*photonCnt));
   close(shmFile);
#if 0   
   shm_unlink(shmFname);
#else
   unlink(shmFname);
#endif      
   
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (pmaps [t]) {
         if (photonRepTime) {
            sprintf(errmsg, "\nBuilding %s photon map...\n", pmapName [t]);
            eputs(errmsg);
            fflush(stderr);
         }
         
         /* Build underlying data structure; heap is destroyed */
         buildPhotonMap(pmaps [t], &totalFlux, NULL, numProc);
      }

   /* Precompute photon irradiance if necessary */
   if (preCompPmap)
      preComputeGlobal(preCompPmap);
}



void photonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon density estimate. Returns irradiance at ray -> rop. */
{
   unsigned                      i;
   float                         r;
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
      
      /* Average radius between furthest two photons to improve accuracy */      
      r = max(sqn -> dist2, (sqn + 1) -> dist2);
      r = 0.25 * (pmap -> maxDist2 + r + 2 * sqrt(pmap -> maxDist2 * r));   
      
      /* Skip the extra photon */
      for (i = 1 ; i < pmap -> squeue.tail; i++, sqn++) {
         photon = getNearestPhoton(&pmap -> squeue, sqn -> idx);
         getPhotonFlux(photon, flux);         
#ifdef PMAP_EPANECHNIKOV
         /* Apply Epanechnikov kernel to photon flux based on photon dist */
         scalecolor(flux, 2 * (1 - sqn -> dist2 / r));
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
   Photon p;
   
   setcolor(irrad, 0, 0, 0);

   /* Ignore sources */
   if (r -> ro && islight(objptr(r -> ro -> omod) -> otype)) 
      return;
      
   find1Photon(preCompPmap, r, &p);
   getPhotonFlux(&p, irrad);
}



void volumePhotonDensity (PhotonMap *pmap, RAY *ray, COLOR irrad)
/* Photon volume density estimate. Returns irradiance at ray -> rop. */
{
   unsigned                      i;
   float                         r, gecc2, ph;
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
      
      /* Average radius between furthest two photons to improve accuracy */      
      r = max(sqn -> dist2, (sqn + 1) -> dist2);
      r = 0.25 * (pmap -> maxDist2 + r + 2 * sqrt(pmap -> maxDist2 * r));   
      
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
      scalecolor(irrad, 3 / (16 * PI * PI * r * sqrt(r)));
      return;
   }
#if 0
   else 
      /* Apply bias compensation to density estimate */
      volumeBiasComp(pmap, ray, irrad);
#endif      
}
