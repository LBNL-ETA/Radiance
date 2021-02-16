#ifndef lint
static const char RCSid[] = "$Id: pmapcontrib.c,v 2.20 2021/02/16 20:06:06 greg Exp $";
#endif

/* 
   ======================================================================
   Photon map for light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapcontrib.c,v 2.20 2021/02/16 20:06:06 greg Exp $
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
#if NIX
   #include <sys/mman.h>
   #include <sys/wait.h>   
#endif


static PhotonPrimaryIdx newPhotonPrimary (PhotonMap *pmap, 
                                          const RAY *primRay, 
                                          FILE *primHeap)
/* Add primary ray for emitted photon and save light source index, origin on
 * source, and emitted direction; used by contrib photons. The current
 * primary is stored in pmap -> lastPrimary.  If the previous primary
 * contributed photons (has srcIdx >= 0), it's appended to primHeap.  If
 * primRay == NULL, the current primary is still flushed, but no new primary
 * is set.  Returns updated primary counter pmap -> numPrimary.  */
{
   if (!pmap || !primHeap)
      return 0;
      
   /* Check if last primary ray has spawned photons (srcIdx >= 0, see
    * newPhoton()), in which case we save it to the primary heap file
    * before clobbering it */
   if (pmap -> lastPrimary.srcIdx >= 0) {
      if (!fwrite(&pmap -> lastPrimary, sizeof(PhotonPrimary), 1, primHeap))
         error(SYSTEM, "failed writing photon primary in newPhotonPrimary");
         
      pmap -> numPrimary++;
      if (pmap -> numPrimary > PMAP_MAXPRIMARY)
         error(INTERNAL, "photon primary overflow in newPhotonPrimary");
   }

   /* Mark unused with negative source index until path spawns a photon (see
    * newPhoton()) */
   pmap -> lastPrimary.srcIdx = -1;
     
   if (primRay) { 
      FVECT dvec;

#ifdef PMAP_PRIMARYDIR            
      /* Reverse incident direction to point to light source */
      dvec [0] = -primRay -> rdir [0];
      dvec [1] = -primRay -> rdir [1];
      dvec [2] = -primRay -> rdir [2];
      pmap -> lastPrimary.dir = encodedir(dvec);
#endif      
#ifdef PMAP_PRIMARYPOS      
      VCOPY(pmap -> lastPrimary.pos, primRay -> rop);
#endif      
   }
   
   return pmap -> numPrimary;
}



#ifdef DEBUG_PMAP
static int checkPrimaryHeap (FILE *file)
/* Check heap for ordered primaries */
{
   Photon   p, lastp;
   int      i, dup;
   
   rewind(file);
   memset(&lastp, 0, sizeof(lastp));
   
   while (fread(&p, sizeof(p), 1, file)) {
      dup = 1;
      
      for (i = 0; i <= 2; i++) {
         if (p.pos [i] < thescene.cuorg [i] || 
             p.pos [i] > thescene.cuorg [i] + thescene.cusize) { 
             
            sprintf(errmsg, "corrupt photon in heap at [%f, %f, %f]\n", 
                    p.pos [0], p.pos [1], p.pos [2]);
            error(WARNING, errmsg);
         }
         
         dup &= p.pos [i] == lastp.pos [i];
      }
      
      if (dup) {
         sprintf(errmsg,
                 "consecutive duplicate photon in heap at [%f, %f, %f]\n", 
                 p.pos [0], p.pos [1], p.pos [2]);
         error(WARNING, errmsg);
      }
   }

   return 0;
}
#endif



static PhotonPrimaryIdx buildPrimaries (PhotonMap *pmap, FILE **primaryHeap,
                                        char **primaryHeapFname,
                                        PhotonPrimaryIdx *primaryOfs, 
                                        unsigned numHeaps)
/* Consolidate per-subprocess photon primary heaps into the primary array
 * pmap -> primaries. Returns offset for primary index linearisation in
 * numPrimary. The heap files in primaryHeap are closed on return. */
{
   PhotonPrimaryIdx  heapLen;
   unsigned          heap;
   
   if (!pmap || !primaryHeap || !primaryOfs || !numHeaps)
      return 0;
      
   pmap -> numPrimary = 0;
   
   for (heap = 0; heap < numHeaps; heap++) {
      primaryOfs [heap] = pmap -> numPrimary;
      
      if (fseek(primaryHeap [heap], 0, SEEK_END) < 0)
         error(SYSTEM, "failed photon primary seek in buildPrimaries");
      pmap -> numPrimary += heapLen = ftell(primaryHeap [heap]) / 
                                      sizeof(PhotonPrimary);      

      pmap -> primaries = realloc(pmap -> primaries, 
                                  pmap -> numPrimary * 
                                  sizeof(PhotonPrimary));  
      if (!pmap -> primaries)
         error(SYSTEM, "failed photon primary alloc in buildPrimaries");

      rewind(primaryHeap [heap]);
      if (fread(pmap -> primaries + primaryOfs [heap], sizeof(PhotonPrimary),
                heapLen, primaryHeap [heap]) != heapLen)
         error(SYSTEM, "failed reading photon primaries in buildPrimaries");
      
      fclose(primaryHeap [heap]);
      unlink(primaryHeapFname [heap]);
   }
   
   return pmap -> numPrimary;
}      



/* Defs for photon emission counter array passed by sub-processes to parent
 * via shared memory */
typedef  unsigned long  PhotonContribCnt;

/* Indices for photon emission counter array: num photons stored and num
 * emitted per source */
#define  PHOTONCNT_NUMPHOT    0
#define  PHOTONCNT_NUMEMIT(n) (1 + n)






void distribPhotonContrib (PhotonMap* pm, unsigned numProc)
{
   EmissionMap       emap;
   char              errmsg2 [128], shmFname [PMAP_TMPFNLEN];
   unsigned          srcIdx, proc;
   int               shmFile, stat, pid;
   double            *srcFlux,         /* Emitted flux per light source */
                     srcDistribTarget; /* Target photon count per source */
   PhotonContribCnt  *photonCnt;       /* Photon emission counter array */
   unsigned          photonCntSize = sizeof(PhotonContribCnt) * 
                                     PHOTONCNT_NUMEMIT(nsources);
   FILE              **primaryHeap = NULL;
   char              **primaryHeapFname = NULL;
   PhotonPrimaryIdx  *primaryOfs = NULL;
                                    
   if (!pm)
      error(USER, "no photon map defined in distribPhotonContrib");
      
   if (!nsources)
      error(USER, "no light sources in distribPhotonContrib");

   /* Allocate photon flux per light source; this differs for every 
    * source as all sources contribute the same number of distributed
    * photons (srcDistribTarget), hence the number of photons emitted per
    * source does not correlate with its emitted flux. The resulting flux
    * per photon is therefore adjusted individually for each source. */
   if (!(srcFlux = calloc(nsources, sizeof(double))))
      error(SYSTEM, "can't allocate source flux in distribPhotonContrib");

   /* ===================================================================
    * INITIALISATION - Set up emission and scattering funcs
    * =================================================================== */
   emap.samples = NULL;
   emap.src = NULL;
   emap.maxPartitions = MAXSPART;
   emap.partitions = (unsigned char*)malloc(emap.maxPartitions >> 1);
   if (!emap.partitions)
      error(USER, "can't allocate source partitions in distribPhotonContrib");

   /* Initialise contrib photon map */   
   initPhotonMap(pm, PMAP_TYPE_CONTRIB);
   initPhotonHeap(pm);
   initPhotonEmissionFuncs();
   initPhotonScatterFuncs();
   
   /* Per-subprocess / per-source target counts */
   pm -> distribTarget /= numProc;
   srcDistribTarget = nsources ? (double)pm -> distribTarget / nsources : 0;   
   
   if (!pm -> distribTarget)
      error(INTERNAL, "no photons to distribute in distribPhotonContrib");
   
   /* Get photon ports from modifier list */
   getPhotonPorts(photonPortList);
      
   /* Get photon sensor modifiers */
   getPhotonSensors(photonSensorList);      

#if NIX   
   /* Set up shared mem for photon counters (zeroed by ftruncate) */
   strcpy(shmFname, PMAP_TMPFNAME);
   shmFile = mkstemp(shmFname);
   
   if (shmFile < 0 || ftruncate(shmFile, photonCntSize) < 0)
      error(SYSTEM, "failed shared mem init in distribPhotonContrib");

   photonCnt = mmap(NULL, photonCntSize, PROT_READ | PROT_WRITE, 
                    MAP_SHARED, shmFile, 0);
                     
   if (photonCnt == MAP_FAILED)
      error(SYSTEM, "failed shared mem mapping in distribPhotonContrib");
#else
   /* Allocate photon counters statically on Windoze */
   if (!(photonCnt = malloc(photonCntSize)))
      error(SYSTEM, "failed trivial malloc in distribPhotonContrib");
   
   for (srcIdx = 0; srcIdx < PHOTONCNT_NUMEMIT(nsources); srcIdx++)
      photonCnt [srcIdx] = 0;
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

   /* =============================================================
    * FLUX INTEGRATION - Get total flux emitted from sources/ports
    * ============================================================= */   
   for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
      unsigned portCnt = 0;      
      srcFlux [srcIdx] = 0;
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
            srcFlux [srcIdx] += colorAvg(emap.partFlux);
         }
         
         portCnt++;
      } while (portCnt < numPhotonPorts);         
      
      if (srcFlux [srcIdx] < FTINY) {
         sprintf(errmsg, "source %s has zero emission", 
                 source [srcIdx].so -> oname);
         error(WARNING, errmsg);
      }
   }   
   
   /* Allocate & init per-subprocess primary heap files */
   primaryHeap = calloc(numProc, sizeof(FILE*));
   primaryHeapFname = calloc(numProc, sizeof(char*));
   primaryOfs = calloc(numProc, sizeof(PhotonPrimaryIdx));
   if (!primaryHeap || !primaryHeapFname || !primaryOfs)
      error(SYSTEM, "failed primary heap allocation in "
            "distribPhotonContrib");
      
   for (proc = 0; proc < numProc; proc++) {
      primaryHeapFname [proc] = malloc(PMAP_TMPFNLEN);
      if (!primaryHeapFname [proc])
         error(SYSTEM, "failed primary heap file allocation in "
               "distribPhotonContrib");
               
      mktemp(strcpy(primaryHeapFname [proc], PMAP_TMPFNAME));
      if (!(primaryHeap [proc] = fopen(primaryHeapFname [proc], "w+b")))
         error(SYSTEM, "failed opening primary heap file in "
               "distribPhotonContrib");
   }               

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
         /* SUBPROCESS ENTERS HERE; opened and mmapped files inherited */
#else
      if (1) {
         /* No subprocess under Windoze */
#endif   
         /* Local photon counters for this subprocess */
         unsigned long  lastNumPhotons = 0, localNumEmitted = 0;
         double         photonFluxSum = 0;   /* Accum. photon flux */

         /* Seed RNGs from PID for decorellated photon distribution */
         pmapSeed(randSeed + proc, partState);
         pmapSeed(randSeed + (proc + 1) % numProc, emitState);
         pmapSeed(randSeed + (proc + 2) % numProc, cntState);
         pmapSeed(randSeed + (proc + 3) % numProc, mediumState);
         pmapSeed(randSeed + (proc + 4) % numProc, scatterState);
         pmapSeed(randSeed + (proc + 5) % numProc, rouletteState);

#ifdef PMAP_SIGUSR                       
   double partNumEmit;
   unsigned long partEmitCnt;
   double srcPhotonFlux, avgPhotonFlux;
   unsigned       portCnt, passCnt, prePassCnt;
   float          srcPreDistrib;
   double         srcNumEmit;     /* # to emit from source */
   unsigned long  srcNumDistrib;  /* # stored */

   void sigUsrDiags()
   /* Loop diags via SIGUSR1 */
   {
      sprintf(errmsg, 
              "********************* Proc %d Diags *********************\n"
              "srcIdx = %d (%s)\nportCnt = %d (%s)\npassCnt = %d\n"
              "srcFlux = %f\nsrcPhotonFlux = %f\navgPhotonFlux = %f\n"
              "partNumEmit = %f\npartEmitCnt = %lu\n\n",
              proc, srcIdx, findmaterial(source [srcIdx].so) -> oname, 
              portCnt, photonPorts [portCnt].so -> oname,
              passCnt, srcFlux [srcIdx], srcPhotonFlux, avgPhotonFlux,
              partNumEmit, partEmitCnt);
      eputs(errmsg);
      fflush(stderr);
   }
#endif
         
#ifdef PMAP_SIGUSR
         signal(SIGUSR1, sigUsrDiags);
#endif         

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

         /* =============================================================
          * 2-PASS PHOTON DISTRIBUTION
          * Pass 1 (pre):  emit fraction of target photon count
          * Pass 2 (main): based on outcome of pass 1, estimate remaining 
          *                number of photons to emit to approximate target 
          *                count
          * ============================================================= */
         for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
#ifndef PMAP_SIGUSR         
            unsigned       portCnt, passCnt = 0, prePassCnt = 0;
            float          srcPreDistrib = preDistrib;
            double         srcNumEmit = 0;       /* # to emit from source */
            unsigned long  srcNumDistrib = pm -> numPhotons;  /* # stored */
#else
            passCnt = prePassCnt = 0;
            srcPreDistrib = preDistrib;
            srcNumEmit = 0;       /* # to emit from source */
            srcNumDistrib = pm -> numPhotons;  /* # stored */
#endif            

            if (srcFlux [srcIdx] < FTINY)
               continue;
                        
            while (passCnt < 2) {
               if (!passCnt) {   
                  /* INIT PASS 1 */
                  if (++prePassCnt > maxPreDistrib) {
                     /* Warn if no photons contributed after sufficient
                      * iterations; only output from subprocess 0 to reduce
                      * console clutter */
                     if (!proc) {
                        sprintf(errmsg, 
                                "source %s: too many prepasses, skipped",
                                source [srcIdx].so -> oname);
                        error(WARNING, errmsg);
                     }

                     break;
                  }
                  
                  /* Num to emit is fraction of target count */
                  srcNumEmit = srcPreDistrib * srcDistribTarget;
               }
               else {
                  /* INIT PASS 2 */
#ifndef PMAP_SIGUSR
                  double srcPhotonFlux, avgPhotonFlux;
#endif
                  
                  /* Based on the outcome of the predistribution we can now
                   * figure out how many more photons we have to emit from
                   * the current source to meet the target count,
                   * srcDistribTarget. This value is clamped to 0 in case
                   * the target has already been exceeded in pass 1.
                   * srcNumEmit and srcNumDistrib is the number of photons
                   * emitted and distributed (stored) from the current
                   * source in pass 1, respectively. */
                  srcNumDistrib = pm -> numPhotons - srcNumDistrib;
                  srcNumEmit *= srcNumDistrib 
                                ? max(srcDistribTarget/srcNumDistrib, 1) - 1
                                : 0;

                  if (!srcNumEmit)
                     /* No photons left to distribute in main pass */
                     break;
                     
                  srcPhotonFlux = srcFlux [srcIdx] / srcNumEmit;
                  avgPhotonFlux = photonFluxSum / (srcIdx + 1);
                  
                  if (avgPhotonFlux > FTINY && 
                      srcPhotonFlux / avgPhotonFlux < FTINY) {
                     /* Skip source if its photon flux is grossly below the
                      * running average, indicating negligible contributions
                      * at the expense of excessive distribution time; only
                      * output from subproc 0 to reduce console clutter */
                     if (!proc) {
                        sprintf(errmsg, 
                                "source %s: itsy bitsy photon flux, skipped",
                                source [srcIdx].so -> oname);                     
                        error(WARNING, errmsg);
                     }

                     srcNumEmit = 0;   /* Or just break??? */
                  }
                        
                  /* Update sum of photon flux per light source */
                  photonFluxSum += srcPhotonFlux;
               }
                              
               portCnt = 0;
               do {    /* Need at least one iteration if no ports! */
                  emap.src = source + srcIdx;
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
#ifndef PMAP_SIGUSR                       
                     double partNumEmit;
                     unsigned long partEmitCnt;
#endif
                     
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
                        
                     /* Update local and shared global emission counter */
                     photonCnt [PHOTONCNT_NUMEMIT(srcIdx)] += partEmitCnt;
                     localNumEmitted += partEmitCnt;                                    
                     
                     /* Integer counter avoids FP rounding errors during
                      * iteration */
                     while (partEmitCnt--) {
                        RAY photonRay;
                     
                        /* Emit photon according to PDF (if any), allocate
                         * associated primary ray, and trace through scene
                         * until absorbed/leaked; emitPhoton() sets the
                         * emitting light source index in photonRay */
                        emitPhoton(&emap, &photonRay);
#if 1
                        if (emap.port)
                           /* !!!  PHOTON PORT REJECTION SAMPLING HACK: set
                            * !!!  photon port as fake hit object for
                            * !!!  primary ray to check for intersection in
                            * !!!  tracePhoton() */                        
                           photonRay.ro = emap.port -> so;
#endif
                        newPhotonPrimary(pm, &photonRay, primaryHeap[proc]);
                        /* Set subprocess index in photonRay for post-
                         * distrib primary index linearisation; this is
                         * propagated with the primary index in photonRay
                         * and set for photon hits by newPhoton() */
                        PMAP_SETRAYPROC(&photonRay, proc);
                        tracePhoton(&photonRay);
                     }
                     
                     /* Update shared global photon count */                     
                     photonCnt [PHOTONCNT_NUMPHOT] += pm -> numPhotons - 
                                                      lastNumPhotons;
                     lastNumPhotons = pm -> numPhotons;
#if !NIX
                     /* Synchronous progress report on Windoze */
                     if (!proc && photonRepTime > 0 && 
                           time(NULL) >= repLastTime + photonRepTime) {
                        unsigned s;                        
                        repComplete = pm -> distribTarget * numProc;
                        repProgress = photonCnt [PHOTONCNT_NUMPHOT];
                        
                        for (repEmitted = 0, s = 0; s < nsources; s++)
                           repEmitted += photonCnt [PHOTONCNT_NUMEMIT(s)];

                        pmapDistribReport();
                     }
#endif
                  }

                  portCnt++;
               } while (portCnt < numPhotonPorts);                  

               if (pm -> numPhotons == srcNumDistrib) {
                  /* Double predistrib factor in case no photons were stored
                   * for this source and redo pass 1 */
                  srcPreDistrib *= 2;
               }
               else {
                  /* Now do pass 2 */
                  passCnt++;
               }
            }
         }
                        
         /* Flush heap buffa one final time to prevent data corruption */
         flushPhotonHeap(pm);         
         /* Flush final photon primary to primary heap file */
         newPhotonPrimary(pm, NULL, primaryHeap [proc]);
         /* Heap files closed automatically on exit
            fclose(pm -> heap);
            fclose(primaryHeap [proc]); */
                  
#ifdef DEBUG_PMAP
         sprintf(errmsg, "Proc %d total %ld photons\n", proc, 
                 pm -> numPhotons);
         eputs(errmsg);
         fflush(stderr);
#endif

#ifdef PMAP_SIGUSR
         signal(SIGUSR1, SIG_DFL);
#endif

#if NIX
         /* Terminate subprocess */
         exit(0);
#endif
      }
      else if (pid < 0)
         error(SYSTEM, "failed to fork subprocess in distribPhotonContrib");
   }

#if NIX
   /* PARENT PROCESS CONTINUES HERE */
#ifdef SIGCONT
   /* Enable progress report signal handler */
   signal(SIGCONT, pmapDistribReport);
#endif
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

      /* Asynchronous progress report from shared subprocess counters */      
      repComplete = pm -> distribTarget * numProc;
      repProgress = photonCnt [PHOTONCNT_NUMPHOT];      
      
      for (repEmitted = 0, srcIdx = 0; srcIdx < nsources; srcIdx++)
         repEmitted += photonCnt [PHOTONCNT_NUMEMIT(srcIdx)];

      /* Get global photon count from shmem updated by subprocs */
      pm -> numPhotons = photonCnt [PHOTONCNT_NUMPHOT];

      if (photonRepTime > 0 && time(NULL) >= repLastTime + photonRepTime)
         pmapDistribReport();
#ifdef SIGCONT
      else signal(SIGCONT, pmapDistribReport);
#endif
   }
#endif /* NIX */

   /* ================================================================
    * POST-DISTRIBUTION - Set photon flux and build kd-tree, etc.
    * ================================================================ */
#ifdef SIGCONT    
   /* Reset signal handler */
   signal(SIGCONT, SIG_DFL);
#endif   
   free(emap.samples);

   if (!pm -> numPhotons)
      error(USER, "empty contribution photon map");

   /* Load per-subprocess primary rays into pm -> primary array */
   /* Dumb compilers apparently need the char** cast */
   pm -> numPrimary = buildPrimaries(pm, primaryHeap, 
                                     (char**)primaryHeapFname,
                                     primaryOfs, numProc);
   if (!pm -> numPrimary)
      error(INTERNAL, "no primary rays in contribution photon map");
   
   /* Set photon flux per source */
   for (srcIdx = 0; srcIdx < nsources; srcIdx++)
      srcFlux [srcIdx] /= photonCnt [PHOTONCNT_NUMEMIT(srcIdx)];
#if NIX
   /* Photon counters no longer needed, unmap shared memory */
   munmap(photonCnt, sizeof(*photonCnt));
   close(shmFile);
   unlink(shmFname);
#else
   free(photonCnt);   
#endif      
   
   if (verbose) {
      eputs("\nBuilding contribution photon map...\n");
#if NIX      
      fflush(stderr);
#endif      
   }
   
   /* Build underlying data structure; heap is destroyed */
   buildPhotonMap(pm, srcFlux, primaryOfs, numProc);
   
   /* Free per-subprocess primary heap files */
   for (proc = 0; proc < numProc; proc++)
      free(primaryHeapFname [proc]);
      
   free(primaryHeapFname);
   free(primaryHeap);
   free(primaryOfs);
   
   if (verbose)
      eputs("\n");
}
