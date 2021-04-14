#ifndef lint
static const char RCSid[] = "$Id: pmapdata.c,v 2.23 2021/04/14 11:26:25 rschregle Exp $";
#endif

/* 
   =========================================================================
   Photon map types and interface to nearest neighbour lookups in underlying
   point cloud data structure.
   
   The default data structure is an in-core kd-tree (see pmapkdt.{h,c}).
   This can be overriden with the PMAP_OOC compiletime switch, which enables
   an out-of-core octree (see oococt.{h,c}).

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme" (FARESYS)) 
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   ==========================================================================
   
   $Id: pmapdata.c,v 2.23 2021/04/14 11:26:25 rschregle Exp $
*/



#include "pmapdata.h"
#include "pmaprand.h"
#include "pmapmat.h"
#include "otypes.h"
#include "source.h"
#include "rcontrib.h"
#include "random.h"



PhotonMap *photonMaps [NUM_PMAP_TYPES] = {
   NULL, NULL, NULL, NULL, NULL, NULL
};



/* Include routines to handle underlying point cloud data structure */
#ifdef PMAP_OOC
   #include "pmapooc.c"
#else
   #include "pmapkdt.c"
#endif

/* Ambient include/exclude set (from ambient.c) */
#ifndef  MAXASET
   #define MAXASET  4095
#endif
extern OBJECT ambset [MAXASET+1];

/* Callback to print photon attributes acc. to user defined format */
int (*printPhoton)(RAY *r, Photon *p, PhotonMap *pm);



void initPhotonMap (PhotonMap *pmap, PhotonMapType t)
/* Init photon map 'n' stuff... */
{
   if (!pmap) 
      return;
      
   pmap -> numPhotons = 0;
   pmap -> biasCompHist = NULL;
   pmap -> maxPos [0] = pmap -> maxPos [1] = pmap -> maxPos [2] = -FHUGE;
   pmap -> minPos [0] = pmap -> minPos [1] = pmap -> minPos [2] = FHUGE;
   pmap -> minGathered = pmap -> maxGathered = pmap -> totalGathered = 0;
   pmap -> gatherTolerance = gatherTolerance;
   pmap -> minError = pmap -> maxError = pmap -> rmsError = 0;
   pmap -> numDensity = 0;
   pmap -> distribRatio = 1;
   pmap -> type = t;
   pmap -> squeue.node = NULL;
   pmap -> squeue.len = 0;   

   /* Init local RNG state */
   pmap -> randState [0] = 10243;
   pmap -> randState [1] = 39829;
   pmap -> randState [2] = 9433;
   pmapSeed(randSeed, pmap -> randState);
   
   /* Set up type-specific photon lookup callback */
   pmap -> lookup = pmapLookup [t];

   /* Mark primary photon ray as unused */
   pmap -> lastPrimary.srcIdx = -1;
   pmap -> numPrimary = 0;  
   pmap -> primaries = NULL;
   
   /* Init storage */
   pmap -> heap = NULL;
   pmap -> heapBuf = NULL;
   pmap -> heapBufLen = 0;
#ifdef PMAP_OOC
   OOC_Null(&pmap -> store);
#else
   kdT_Null(&pmap -> store);
#endif
}



void initPhotonHeap (PhotonMap *pmap)
{
   int fdFlags;
   
   if (!pmap)
      error(INTERNAL, "undefined photon map in initPhotonHeap");
      
   if (!pmap -> heap) {
      /* Open heap file */
      mktemp(strcpy(pmap -> heapFname, PMAP_TMPFNAME));
      if (!(pmap -> heap = fopen(pmap -> heapFname, "w+b")))
         error(SYSTEM, "failed opening heap file in initPhotonHeap");

#ifdef F_SETFL	/* XXX is there an alternate needed for Windows? */
      fdFlags = fcntl(fileno(pmap -> heap), F_GETFL);
      fcntl(fileno(pmap -> heap), F_SETFL, fdFlags | O_APPEND);
#endif/*      ftruncate(fileno(pmap -> heap), 0); */
   }
}



void flushPhotonHeap (PhotonMap *pmap)
{
   int                  fd;
   const unsigned long  len = pmap -> heapBufLen * sizeof(Photon);
   
   if (!pmap)
      error(INTERNAL, "undefined photon map in flushPhotonHeap");

   if (!pmap -> heap || !pmap -> heapBuf) {
      /* Silently ignore undefined heap 
      error(INTERNAL, "undefined heap in flushPhotonHeap"); */
      return;
   }

   /* Atomically seek and write block to heap */
   /* !!! Unbuffered I/O via pwrite() avoids potential race conditions
    * !!! and buffer corruption which can occur with lseek()/fseek()
    * !!! followed by write()/fwrite(). */   
   fd = fileno(pmap -> heap);
   
#ifdef DEBUG_PMAP
   sprintf(errmsg, "Proc %d: flushing %ld photons from pos %ld\n", getpid(), 
           pmap -> heapBufLen, lseek(fd, 0, SEEK_END) / sizeof(Photon)); 
   eputs(errmsg);
#endif

   /*if (pwrite(fd, pmap -> heapBuf, len, lseek(fd, 0, SEEK_END)) != len) */
   if (write(fd, pmap -> heapBuf, len) != len)
      error(SYSTEM, "failed append to heap file in flushPhotonHeap");

#if NIX
   if (fsync(fd))
      error(SYSTEM, "failed fsync in flushPhotonHeap");
#endif

   pmap -> heapBufLen = 0;
}



#ifdef DEBUG_PMAP
static int checkPhotonHeap (FILE *file)
/* Check heap for nonsensical or duplicate photons */
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



int newPhoton (PhotonMap* pmap, const RAY* ray)
{
   unsigned i;
   Photon photon;
   COLOR photonFlux;
   
   /* Account for distribution ratio */
   if (!pmap || pmapRandom(pmap -> randState) > pmap -> distribRatio) 
      return -1;
      
   /* Don't store on sources */
   if (ray -> robj > -1 && islight(objptr(ray -> ro -> omod) -> otype)) 
      return -1;

   /* Ignore photon if modifier in/outside exclude/include set */
   if (ambincl != -1 && ray -> ro && 
       ambincl != inset(ambset, ray -> ro -> omod))
      return -1;

   if (pmapNumROI && pmapROI) {
      unsigned inROI = 0;
      FVECT    photonDist;
      
      /* Store photon if within a region of interest (for ze Ecksperts!) 
         Note size of spherical ROI is squared. */
      for (i = 0; !inROI && i < pmapNumROI; i++) {
         VSUB(photonDist, ray -> rop, pmapROI [i].pos);
         
         inROI = (
            PMAP_ROI_ISSPHERE(pmapROI + i) 
            ?  DOT(photonDist, photonDist) <= pmapROI [i].siz [0] 
            :  fabs(photonDist [0]) <= pmapROI [i].siz [0] &&
               fabs(photonDist [1]) <= pmapROI [i].siz [1] && 
               fabs(photonDist [2]) <= pmapROI [i].siz [2]
         );
      }
      if (!inROI)
         return -1;
   }
   
   /* Adjust flux according to distribution ratio and ray weight */
   copycolor(photonFlux, ray -> rcol);         
#if 0
   /* Factored out ray -> rweight as deprecated (?) for pmap, and infact
      erroneously attenuates volume photon flux based on extinction,
      which is already factored in by photonParticipate() */
   scalecolor(photonFlux, 
              ray -> rweight / (pmap -> distribRatio ? pmap -> distribRatio
                                                     : 1));
#else
   scalecolor(photonFlux, 
              1.0 / (pmap -> distribRatio ? pmap -> distribRatio : 1));
#endif
   setPhotonFlux(&photon, photonFlux);

   /* Set photon position and flags */
   VCOPY(photon.pos, ray -> rop);
   photon.flags = 0;
   photon.caustic = PMAP_CAUSTICRAY(ray);

   /* Set contrib photon's primary ray and subprocess index (the latter
    * to linearise the primary ray indices after photon distribution is
    * complete). Also set primary ray's source index, thereby marking it
    * as used. */
   if (isContribPmap(pmap)) {
      photon.primary = pmap -> numPrimary;
      photon.proc = PMAP_GETRAYPROC(ray);
      pmap -> lastPrimary.srcIdx = ray -> rsrc;
   }
   else photon.primary = 0;
   
   /* Set normal */
   for (i = 0; i <= 2; i++)
      photon.norm [i] = 127.0 * (isVolumePmap(pmap) ? ray -> rdir [i] 
                                                    : ray -> ron [i]);

   if (!pmap -> heapBuf) {
      /* Lazily allocate heap buffa */
#if NIX
      /* Randomise buffa size to temporally decorellate flushes in
       * multiprocessing mode */
      srandom(randSeed + getpid());
      pmap -> heapBufSize = PMAP_HEAPBUFSIZE * (0.5 + frandom());
#else
      /* Randomisation disabled for single processes on WIN; also useful
       * for reproducability during debugging */         
      pmap -> heapBufSize = PMAP_HEAPBUFSIZE;
#endif         
      if (!(pmap -> heapBuf = calloc(pmap -> heapBufSize, sizeof(Photon))))
         error(SYSTEM, "failed heap buffer allocation in newPhoton");
      pmap -> heapBufLen = 0;      
   }

   /* Photon initialised; now append to heap buffa */
   memcpy(pmap -> heapBuf + pmap -> heapBufLen, &photon, sizeof(Photon));
               
   if (++pmap -> heapBufLen >= pmap -> heapBufSize)
      /* Heap buffa full, flush to heap file */
      flushPhotonHeap(pmap);

   pmap -> numPhotons++;
         
   /* Print photon attributes */
   if (printPhoton)
      /* Non-const kludge */
      printPhoton((RAY*)ray, &photon, pmap);
            
   return 0;
}



void buildPhotonMap (PhotonMap *pmap, double *photonFlux, 
                     PhotonPrimaryIdx *primaryOfs, unsigned nproc)
{
   unsigned long  n, nCheck = 0;
   unsigned       i;
   Photon         *p;
   COLOR          flux;
   char           nuHeapFname [sizeof(PMAP_TMPFNAME)];
   FILE           *nuHeap;
   /* Need double here to reduce summation errors */
   double         avgFlux [3] = {0, 0, 0}, CoG [3] = {0, 0, 0}, CoGdist = 0;
   FVECT          d;
   
   if (!pmap)
      error(INTERNAL, "undefined photon map in buildPhotonMap");
      
   /* Get number of photons from heapfile size */
   if (fseek(pmap -> heap, 0, SEEK_END) < 0)
      error(SYSTEM, "failed seek to end of photon heap in buildPhotonMap");
   pmap -> numPhotons = ftell(pmap -> heap) / sizeof(Photon);
   
   if (!pmap -> numPhotons)
      error(INTERNAL, "empty photon map in buildPhotonMap");   

   if (!pmap -> heap)
      error(INTERNAL, "no heap in buildPhotonMap");

#ifdef DEBUG_PMAP
   eputs("Checking photon heap consistency...\n");
   checkPhotonHeap(pmap -> heap);
   
   sprintf(errmsg, "Heap contains %ld photons\n", pmap -> numPhotons);
   eputs(errmsg);
#endif

   /* Allocate heap buffa */
   if (!pmap -> heapBuf) {
      pmap -> heapBufSize = PMAP_HEAPBUFSIZE;
      pmap -> heapBuf = calloc(pmap -> heapBufSize, sizeof(Photon));
      if (!pmap -> heapBuf)
         error(SYSTEM, "failed to allocate postprocessed photon heap in" 
               "buildPhotonMap");
   }

   /* We REALLY don't need yet another @%&*! heap just to hold the scaled
    * photons, but can't think of a quicker fix... */
   mktemp(strcpy(nuHeapFname, PMAP_TMPFNAME));
   if (!(nuHeap = fopen(nuHeapFname, "w+b")))
      error(SYSTEM, "failed to open postprocessed photon heap in "
            "buildPhotonMap");
            
   rewind(pmap -> heap);

#ifdef DEBUG_PMAP 
   eputs("Postprocessing photons...\n");
#endif
   
   while (!feof(pmap -> heap)) {   
#ifdef DEBUG_PMAP 
      printf("Reading %lu at %lu... ", pmap -> heapBufSize, ftell(pmap->heap));
#endif      
      pmap -> heapBufLen = fread(pmap -> heapBuf, sizeof(Photon), 
                                 pmap -> heapBufSize, pmap -> heap);
#ifdef DEBUG_PMAP                                 
      printf("Got %lu\n", pmap -> heapBufLen);
#endif      

      if (ferror(pmap -> heap))
         error(SYSTEM, "failed to read photon heap in buildPhotonMap");

      for (n = pmap -> heapBufLen, p = pmap -> heapBuf; n; n--, p++) {
         /* Update min and max pos and set photon flux */
         for (i = 0; i <= 2; i++) {
            if (p -> pos [i] < pmap -> minPos [i]) 
               pmap -> minPos [i] = p -> pos [i];
            else if (p -> pos [i] > pmap -> maxPos [i]) 
               pmap -> maxPos [i] = p -> pos [i];   

            /* Update centre of gravity with photon position */                 
            CoG [i] += p -> pos [i];                  
         }  
         
         if (primaryOfs)
            /* Linearise photon primary index from subprocess index using the
             * per-subprocess offsets in primaryOfs */
            p -> primary += primaryOfs [p -> proc];
         
         /* Scale photon's flux (hitherto normalised to 1 over RGB); in
          * case of a contrib photon map, this is done per light source,
          * and photonFlux is assumed to be an array */
         getPhotonFlux(p, flux);

         if (photonFlux) {
            scalecolor(flux, photonFlux [isContribPmap(pmap) ? 
                                            photonSrcIdx(pmap, p) : 0]);
            setPhotonFlux(p, flux);
         }

         /* Update average photon flux; need a double here */
         addcolor(avgFlux, flux);
      }
         
      /* Write modified photons to new heap */
      fwrite(pmap -> heapBuf, sizeof(Photon), pmap -> heapBufLen, nuHeap);
                
      if (ferror(nuHeap))
         error(SYSTEM, "failed postprocessing photon flux in "
               "buildPhotonMap");
      
      nCheck += pmap -> heapBufLen;
   }

#ifdef DEBUG_PMAP
   if (nCheck < pmap -> numPhotons)
      error(INTERNAL, "truncated photon heap in buildPhotonMap");
#endif
   
   /* Finalise average photon flux */
   scalecolor(avgFlux, 1.0 / pmap -> numPhotons);
   copycolor(pmap -> photonFlux, avgFlux);

   /* Average photon positions to get centre of gravity */
   for (i = 0; i < 3; i++)
      pmap -> CoG [i] = CoG [i] /= pmap -> numPhotons;
      
   rewind(pmap -> heap);
   
   /* Compute average photon distance to centre of gravity */
   while (!feof(pmap -> heap)) {
      pmap -> heapBufLen = fread(pmap -> heapBuf, sizeof(Photon), 
                                 pmap -> heapBufSize, pmap -> heap);
      
      for (n = pmap -> heapBufLen, p = pmap -> heapBuf; n; n--, p++) {
         VSUB(d, p -> pos, CoG);
         CoGdist += DOT(d, d);
      }
   }   

   pmap -> CoGdist = CoGdist /= pmap -> numPhotons;

   /* Swap heaps, discarding unscaled photons */
   fclose(pmap -> heap);
   unlink(pmap -> heapFname);
   pmap -> heap = nuHeap;
   strcpy(pmap -> heapFname, nuHeapFname);
   
#ifdef PMAP_OOC
   OOC_BuildPhotonMap(pmap, nproc);
#else
   kdT_BuildPhotonMap(pmap);
#endif

   /* Trash heap and its buffa */
   free(pmap -> heapBuf);
   fclose(pmap -> heap);
   unlink(pmap -> heapFname);
   pmap -> heap = NULL;
   pmap -> heapBuf = NULL;
}



/* Dynamic max photon search radius increase and reduction factors */
#define PMAP_MAXDIST_INC 4
#define PMAP_MAXDIST_DEC 0.9

/* Num successful lookups before reducing in max search radius */
#define PMAP_MAXDIST_CNT 1000

/* Threshold below which we assume increasing max radius won't help */
#define PMAP_SHORT_LOOKUP_THRESH 1

/* Coefficient for adaptive maximum search radius */
#define PMAP_MAXDIST_COEFF 100

void findPhotons (PhotonMap* pmap, const RAY* ray)
{
   int redo = 0;
   
   if (!pmap -> squeue.len) {
      /* Lazy init priority queue */
#ifdef PMAP_OOC
      OOC_InitFindPhotons(pmap);
#else
      kdT_InitFindPhotons(pmap);      
#endif      
      pmap -> minGathered = pmap -> maxGather;
      pmap -> maxGathered = pmap -> minGather;
      pmap -> totalGathered = 0;
      pmap -> numLookups = pmap -> numShortLookups = 0;
      pmap -> shortLookupPct = 0;
      pmap -> minError = FHUGE;
      pmap -> maxError = -FHUGE;
      pmap -> rmsError = 0;
      /* SQUARED max search radius limit is based on avg photon distance to
       * centre of gravity, unless fixed by user (maxDistFix > 0) */
      pmap -> maxDist0 = pmap -> maxDist2Limit = 
         maxDistFix > 0 ? maxDistFix * maxDistFix
                        : PMAP_MAXDIST_COEFF * pmap -> squeue.len * 
                          pmap -> CoGdist / pmap -> numPhotons;
   }

   do {
      pmap -> squeue.tail = 0;
      pmap -> maxDist2 = pmap -> maxDist0;
         
      /* Search position is ray -> rorg for volume photons, since we have no 
         intersection point. Normals are ignored -- these are incident 
         directions). */
      /* NOTE: status returned by XXX_FindPhotons() is currently ignored; 
         if no photons are found, an empty queue is returned under the
         assumption all photons are too distant to contribute significant
         flux. */
      if (isVolumePmap(pmap)) {
#ifdef PMAP_OOC
         OOC_FindPhotons(pmap, ray -> rorg, NULL);
#else
         kdT_FindPhotons(pmap, ray -> rorg, NULL);
#endif                  
      }
      else {
#ifdef PMAP_OOC
         OOC_FindPhotons(pmap, ray -> rop, ray -> ron);         
#else
         kdT_FindPhotons(pmap, ray -> rop, ray -> ron);
#endif
      }

#ifdef PMAP_LOOKUP_INFO
      fprintf(stderr, "%d/%d %s photons found within radius %.3f " 
              "at (%.2f,%.2f,%.2f) on %s\n", pmap -> squeue.tail, 
              pmap -> squeue.len, pmapName [pmap -> type], sqrt(pmap -> maxDist2),
              ray -> rop [0], ray -> rop [1], ray -> rop [2], 
              ray -> ro ? ray -> ro -> oname : "<null>");
#endif      
            
      if (pmap -> squeue.tail < pmap -> squeue.len * pmap -> gatherTolerance) {
         /* Short lookup; too few photons found */
         if (pmap -> squeue.tail > PMAP_SHORT_LOOKUP_THRESH) {
            /* Ignore short lookups which return fewer than
             * PMAP_SHORT_LOOKUP_THRESH photons under the assumption there
             * really are no photons in the vicinity, and increasing the max
             * search radius therefore won't help */
#ifdef PMAP_LOOKUP_WARN
            sprintf(errmsg, 
                    "%d/%d %s photons found at (%.2f,%.2f,%.2f) on %s",
                    pmap -> squeue.tail, pmap -> squeue.len, 
                    pmapName [pmap -> type], 
                    ray -> rop [0], ray -> rop [1], ray -> rop [2], 
                    ray -> ro ? ray -> ro -> oname : "<null>");
            error(WARNING, errmsg);         
#endif             

            /* Bail out after warning if maxDist is fixed */
            if (maxDistFix > 0)
               return;
            
            if (pmap -> maxDist0 < pmap -> maxDist2Limit) {
               /* Increase max search radius if below limit & redo search */
               pmap -> maxDist0 *= PMAP_MAXDIST_INC;
#ifdef PMAP_LOOKUP_REDO
               redo = 1;
#endif               
#ifdef PMAP_LOOKUP_WARN
               sprintf(errmsg, 
                       redo ? "restarting photon lookup with max radius %.1e"
                            : "max photon lookup radius adjusted to %.1e",
                       pmap -> maxDist0);
               error(WARNING, errmsg);
#endif
            }
#ifdef PMAP_LOOKUP_REDO
            else {
               sprintf(errmsg, "max photon lookup radius clamped to %.1e",
                       pmap -> maxDist0);
               error(WARNING, errmsg);
            }
#endif
         }
         
         /* Reset successful lookup counter */
         pmap -> numLookups = 0;
      }
      else {
         /* Bail out after warning if maxDist is fixed */
         if (maxDistFix > 0)
            return;
            
         /* Increment successful lookup counter and reduce max search radius if
          * wraparound */
         pmap -> numLookups = (pmap -> numLookups + 1) % PMAP_MAXDIST_CNT;
         if (!pmap -> numLookups)
            pmap -> maxDist0 *= PMAP_MAXDIST_DEC;
            
         redo = 0;
      }
      
   } while (redo);
}



Photon *find1Photon (PhotonMap *pmap, const RAY* ray, Photon *photon)
{
   /* Init (squared) search radius to avg photon dist to centre of gravity */
   float maxDist2_0 = pmap -> CoGdist;   
   int found = 0;   
#ifdef PMAP_LOOKUP_REDO
   #define REDO 1
#else
   #define REDO 0
#endif
   
   do {
      pmap -> maxDist2 = maxDist2_0;
#ifdef PMAP_OOC
      found = OOC_Find1Photon(pmap, ray -> rop, ray -> ron, photon);
#else
      found = kdT_Find1Photon(pmap, ray -> rop, ray -> ron, photon);
#endif
      if (found < 0) {
         /* Expand search radius to retry */
         maxDist2_0 *= 2;      
#ifdef PMAP_LOOKUP_WARN
         sprintf(errmsg, "failed 1-NN photon lookup"
#ifdef PMAP_LOOKUP_REDO
            ", retrying with search radius %.2f", maxDist2_0
#endif
         );
         error(WARNING, errmsg);
#endif
      }
   } while (REDO && found < 0);

   /* Return photon buffer containing valid photon, else NULL */
   return found < 0 ? NULL : photon;
}



void getPhoton (PhotonMap *pmap, PhotonIdx idx, Photon *photon)
{
#ifdef PMAP_OOC
   if (OOC_GetPhoton(pmap, idx, photon))      
#else
   if (kdT_GetPhoton(pmap, idx, photon))
#endif
      error(INTERNAL, "failed photon lookup");
}



Photon *getNearestPhoton (const PhotonSearchQueue *squeue, PhotonIdx idx)
{
#ifdef PMAP_OOC
   return OOC_GetNearestPhoton(squeue, idx);
#else
   return kdT_GetNearestPhoton(squeue, idx);
#endif
}



PhotonIdx firstPhoton (const PhotonMap *pmap)
{
#ifdef PMAP_OOC
   return OOC_FirstPhoton(pmap);
#else
   return kdT_FirstPhoton(pmap);
#endif
}



void deletePhotons (PhotonMap* pmap)
{
#ifdef PMAP_OOC
   OOC_Delete(&pmap -> store);
#else
   kdT_Delete(&pmap -> store);
#endif

   free(pmap -> squeue.node);
   free(pmap -> biasCompHist);
   
   pmap -> numPhotons = pmap -> minGather = pmap -> maxGather =
      pmap -> squeue.len = pmap -> squeue.tail = 0;
}
