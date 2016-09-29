#ifndef lint
static const char RCSid[] = "$Id: pmapcontrib.c,v 2.13 2016/09/29 21:51:58 greg Exp $";
#endif

/* 
   ======================================================================
   Photon map support for building light source contributions

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapcontrib.c,v 2.13 2016/09/29 21:51:58 greg Exp $
*/


#include "pmapcontrib.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmaprand.h"
#include "pmapio.h"
#include "pmapdiag.h"
#include "rcontrib.h"
#include "otypes.h"
#include <sys/mman.h>
#include <sys/wait.h>



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
    * newPhoton()), in which case we write it to the primary heap file
    * before overwriting it */
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
      
      /* Reverse incident direction to point to light source */
      dvec [0] = -primRay -> rdir [0];
      dvec [1] = -primRay -> rdir [1];
      dvec [2] = -primRay -> rdir [2];
      pmap -> lastPrimary.dir = encodedir(dvec);
#ifdef PMAP_PRIMARYPOS      
      VCOPY(pmap -> lastPrimary.pos, primRay -> rop);
#endif      
   }
   
   return pmap -> numPrimary;
}



#ifdef DEBUG_PMAP_CONTRIB
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
      
      if (fseek(primaryHeap [heap], 0, SEEK_END))
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
   char              errmsg2 [128], shmFname [255];
   unsigned          srcIdx, proc;
   int               shmFile, stat, pid;
   double            *srcFlux,         /* Emitted flux per light source */
                     srcDistribTarget; /* Target photon count per source */
   PhotonContribCnt  *photonCnt;       /* Photon emission counter array */
   const unsigned    photonCntSize = sizeof(PhotonContribCnt) * 
                                     PHOTONCNT_NUMEMIT(nsources);
   FILE              *primaryHeap [numProc];
   PhotonPrimaryIdx  primaryOfs [numProc];
                                    
   if (!pm)
      error(USER, "no photon map defined in distribPhotonContrib");
      
   if (!nsources)
      error(USER, "no light sources in distribPhotonContrib");

   if (nsources > MAXMODLIST)
      error(USER, "too many light sources in distribPhotonContrib");
      
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
   
   if (shmFile < 0 || ftruncate(shmFile, photonCntSize) < 0)
      error(SYSTEM, "failed shared mem init in distribPhotonContrib");

   photonCnt = mmap(NULL, photonCntSize, PROT_READ | PROT_WRITE, 
                    MAP_SHARED, shmFile, 0);
                     
   if (photonCnt == MAP_FAILED)
      error(SYSTEM, "failed shared mem mapping in distribPhotonContrib");

   /* =============================================================
    * FLUX INTEGRATION - Get total flux emitted from light source
    * ============================================================= */   
   for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
      unsigned portCnt = 0;
      
      srcFlux [srcIdx] = 0;
      emap.src = source + srcIdx;
      
      if (photonRepTime) 
         eputs("\n");

      do {  /* Need at least one iteration if no ports! */      
         emap.port = emap.src -> sflags & SDISTANT ? photonPorts + portCnt 
                                                   : NULL;
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
            
            sprintf(errmsg2, "(%lu partitions)...\n", emap.numPartitions);
            strcat(errmsg, errmsg2);
            eputs(errmsg);
            fflush(stderr);
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

   if (photonRepTime) 
      eputs("\n");

   /* Init per-subprocess primary heap files */
   for (proc = 0; proc < numProc; proc++)
      if (!(primaryHeap [proc] = tmpfile()))
         error(SYSTEM, "failed opening primary heap file in "
               "distribPhotonContrib");
               
   /* MAIN LOOP */
   for (proc = 0; proc < numProc; proc++) {
      if (!(pid = fork())) {
         /* SUBPROCESS ENTERS HERE; 
          * all opened and memory mapped files are inherited */
          
         /* Local photon counters for this subprocess */
         unsigned long  lastNumPhotons = 0, localNumEmitted = 0;
         double         photonFluxSum = 0;   /* Running photon flux sum */

         /* Seed RNGs from PID for decorellated photon distribution */
         pmapSeed(randSeed + proc, partState);
         pmapSeed(randSeed + proc, emitState);
         pmapSeed(randSeed + proc, cntState);
         pmapSeed(randSeed + proc, mediumState);
         pmapSeed(randSeed + proc, scatterState);
         pmapSeed(randSeed + proc, rouletteState);
         
         /* =============================================================
          * 2-PASS PHOTON DISTRIBUTION
          * Pass 1 (pre):  emit fraction of target photon count
          * Pass 2 (main): based on outcome of pass 1, estimate remaining 
          *                number of photons to emit to approximate target 
          *                count
          * ============================================================= */         
         for (srcIdx = 0; srcIdx < nsources; srcIdx++) {
            unsigned       portCnt, passCnt = 0, prePassCnt = 0;
            float          srcPreDistrib = preDistrib;
            double         srcNumEmit = 0;       /* # to emit from source */
            unsigned long  srcNumDistrib = pm -> numPhotons;  /* # stored */

            if (srcFlux [srcIdx] < FTINY)
               continue;
                        
            while (passCnt < 2) {
               if (!passCnt) {   
                  /* INIT PASS 1 */
                  if (++prePassCnt > maxPreDistrib) {
                     /* Warn if no photons contributed after sufficient
                      * iterations */
                     sprintf(errmsg, "proc %d, source %s: "
                             "too many prepasses, skipped",
                             proc, source [srcIdx].so -> oname);
                     error(WARNING, errmsg);
                     break;
                  }
                  
                  /* Num to emit is fraction of target count */
                  srcNumEmit = srcPreDistrib * srcDistribTarget;
               }
               else {
                  /* INIT PASS 2 */
                  double srcPhotonFlux, avgPhotonFlux;
                  
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
                  
                  if (avgPhotonFlux > 0 && 
                      srcPhotonFlux / avgPhotonFlux < FTINY) {
                     /* Skip source if its photon flux is grossly below the
                      * running average, indicating negligible contribs at
                      * the expense of excessive distribution time */
                     sprintf(errmsg, "proc %d, source %s: "
                             "itsy bitsy photon flux, skipped",
                             proc, source [srcIdx].so -> oname);
                     error(WARNING, errmsg);
                     srcNumEmit = 0;
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
                  
                  if (photonRepTime && !proc) {
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
                     localNumEmitted += partEmitCnt;                                    
                     photonCnt [PHOTONCNT_NUMEMIT(srcIdx)] += partEmitCnt;
                     
                     /* Integer counter avoids FP rounding errors */
                     while (partEmitCnt--) {
                        RAY photonRay;
                     
                        /* Emit photon according to PDF (if any), allocate
                         * associated primary ray, and trace through scene
                         * until absorbed/leaked; emitPhoton() sets the
                         * emitting light source index in photonRay */
                        emitPhoton(&emap, &photonRay);
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
                  }

                  portCnt++;
               } while (portCnt < numPhotonPorts);                  

               if (pm -> numPhotons == srcNumDistrib) 
                  /* Double predistrib factor in case no photons were stored
                   * for this source and redo pass 1 */
                  srcPreDistrib *= 2;
               else {
                  /* Now do pass 2 */
                  passCnt++;
/*                if (photonRepTime)
                     eputs("\n"); */
               }
            }
         }
                        
         /* Flush heap buffa one final time to prevent data corruption */
         flushPhotonHeap(pm);
         fclose(pm -> heap);
         
         /* Flush final photon primary to primary heap file */
         newPhotonPrimary(pm, NULL, primaryHeap [proc]);
         fclose(primaryHeap [proc]);
                  
#ifdef DEBUG_PMAP
         sprintf(errmsg, "Proc %d exited with total %ld photons\n", proc, 
                 pm -> numPhotons);
         eputs(errmsg);
#endif

         exit(0);
      }
      else if (pid < 0)
         error(SYSTEM, "failed to fork subprocess in distribPhotonContrib");
   }

   /* PARENT PROCESS CONTINUES HERE */
   /* Record start time and enable progress report signal handler */
   repStartTime = time(NULL);
#ifdef SIGCONT
   signal(SIGCONT, pmapDistribReport);
#endif
/*
   if (photonRepTime)
      eputs("\n"); */
   
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

   /* ================================================================
    * POST-DISTRIBUTION - Set photon flux and build kd-tree, etc.
    * ================================================================ */
#ifdef SIGCONT    
   signal(SIGCONT, SIG_DFL);
#endif   
   free(emap.samples);

   if (!pm -> numPhotons)
      error(USER, "empty photon map");

   /* Load per-subprocess primary rays into pm -> primary array */
   pm -> numPrimary = buildPrimaries(pm, primaryHeap, primaryOfs, numProc);
   if (!pm -> numPrimary)
      error(INTERNAL, "no primary rays in contribution photon map");
   
   /* Set photon flux per source */
   for (srcIdx = 0; srcIdx < nsources; srcIdx++)
      srcFlux [srcIdx] /= photonCnt [PHOTONCNT_NUMEMIT(srcIdx)];

   /* Photon counters no longer needed, unmap shared memory */
   munmap(photonCnt, sizeof(*photonCnt));
   close(shmFile);
#if 0   
   shm_unlink(shmFname);
#else
   unlink(shmFname);
#endif      
   
   if (photonRepTime) {
      eputs("\nBuilding contrib photon map...\n");
      fflush(stderr);
   }
   
   /* Build underlying data structure; heap is destroyed */
   buildPhotonMap(pm, srcFlux, primaryOfs, numProc);   
}
