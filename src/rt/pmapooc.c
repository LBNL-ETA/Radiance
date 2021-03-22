/* 
   ======================================================================
   Photon map interface to out-of-core octree

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmapooc.c,v 1.7 2021/03/22 23:00:00 rschregle Exp $
*/



#include "pmapdata.h"   /* Includes pmapooc.h */
#include "source.h"
#include "otspecial.h"
#include "oocsort.h"
#include "oocbuild.h"



/* Returns photon position as sorting key for OOC_Octree & friends (notably
 * for Morton code generation).
 * !!! Uses type conversion from float via TEMPORARY storage; 
 * !!! THIS IS NOT THREAD SAFE!
 * !!! RETURNED KEY PERSISTS ONLY IF COPIED BEFORE NEXT CALL! */
RREAL *OOC_PhotonKey (const void *p)
{
   static FVECT photonPos; /* Temp storage for type conversion */
   
   VCOPY(photonPos, ((Photon*)p) -> pos);
   return photonPos;
}
           


#ifdef DEBUG_OOC
static int OOC_CheckKeys (FILE *file, const OOC_Octree *oct)
{
   Photon         p, lastp;
   RREAL          *key;
   OOC_MortonIdx  zIdx, lastzIdx = 0;
   
   rewind(file);
   memset(&lastp, 0, sizeof(lastp));
   
   while (fread(&p, sizeof(p), 1, file) > 0) {
      key = OOC_PhotonKey(&p);
      zIdx = OOC_KEY2MORTON(key, oct);
      
      if (zIdx < lastzIdx) {
         error(INTERNAL, "photons not sorted");
         return -1;
      }
      
      if (zIdx == lastzIdx) {
         sprintf(errmsg, "identical key %021ld "
                 "for [%.9f, %.9f, %.9f]\tand [%.9f, %.9f, %.9f]",
                 zIdx, lastp.pos [0], lastp.pos [1], lastp.pos [2], 
                 p.pos [0], p.pos [1], p.pos [2]);
         error(WARNING, errmsg);
      } 
            
      lastzIdx = zIdx;
      memcpy(&lastp, &p, sizeof(p));
   }

   return 0;
}
#endif



void OOC_BuildPhotonMap (struct PhotonMap *pmap, unsigned numProc)
{
   FILE        *leafFile;
   char        leafFname [1024];
   FVECT       d, octOrg;
   int         i;
   RREAL       octSize = 0;
   
   /* Determine octree size and origin from pmap extent and init octree */
   VCOPY(octOrg, pmap -> minPos);
   VSUB(d, pmap -> maxPos, pmap -> minPos);
   for (i = 0; i < 3; i++)
      if (octSize < d [i])
         octSize = d [i];
         
   if (octSize < FTINY)
      error(INTERNAL, "zero octree size in OOC_BuildPhotonMap");
            
   /* Derive leaf filename from photon map and open file */
   strncpy(leafFname, pmap -> fileName, sizeof(leafFname));
   strncat(leafFname, PMAP_OOC_LEAFSUFFIX, 
           sizeof(leafFname) - strlen(leafFname) - 1);   
   if (!(leafFile = fopen(leafFname, "w+b")))
      error(SYSTEM, "failed to open leaf file in OOC_BuildPhotonMap");

#ifdef DEBUG_OOC
   eputs("Sorting photons by Morton code...\n");
#endif
   
   /* Sort photons in heapfile by Morton code and write to leaf file */
   if (OOC_Sort(pmap -> heap, leafFile, PMAP_OOC_NUMBLK, PMAP_OOC_BLKSIZE,
                numProc, sizeof(Photon), octOrg, octSize, OOC_PhotonKey))
      error(INTERNAL, "failed out-of-core photon sort in OOC_BuildPhotonMap");

   /* Init and build octree */   
   OOC_Init(&pmap -> store, sizeof(Photon), octOrg, octSize, OOC_PhotonKey, 
            leafFile);
            
#ifdef DEBUG_OOC
   eputs("Checking leaf file consistency...\n");
   OOC_CheckKeys(leafFile, &pmap -> store);
   
   eputs("Building out-of-core octree...\n");
#endif   
            
   if (!OOC_Build(&pmap -> store, PMAP_OOC_LEAFMAX, PMAP_OOC_MAXDEPTH))
      error(INTERNAL, "failed out-of-core octree build in OOC_BuildPhotonMap");
      
#ifdef DEBUG_OOC
   eputs("Checking out-of-core octree consistency...\n");
   if (OOC_Check(&pmap -> store, OOC_ROOT(&pmap -> store), 
                 octOrg, octSize, 0))
      error(INTERNAL, "inconsistent out-of-core octree; Time4Harakiri");
#endif
}



int OOC_SavePhotons (const struct PhotonMap *pmap, FILE *out)
{
   return OOC_SaveOctree(&pmap -> store, out);
}



int OOC_LoadPhotons (struct PhotonMap *pmap, FILE *nodeFile)
{
   FILE  *leafFile;
   char  leafFname [1024];
   
   /* Derive leaf filename from photon map and open file */
   strncpy(leafFname, pmap -> fileName, sizeof(leafFname));
   strncat(leafFname, PMAP_OOC_LEAFSUFFIX, 
           sizeof(leafFname) - strlen(leafFname) - 1);   
   
   if (!(leafFile = fopen(leafFname, "r")))
      error(SYSTEM, "failed to open leaf file in OOC_LoadPhotons");

   if (OOC_LoadOctree(&pmap -> store, nodeFile, OOC_PhotonKey, leafFile))
      return -1;

#ifdef DEBUG_OOC
   /* Check octree for consistency */
   if (OOC_Check(
      &pmap -> store, OOC_ROOT(&pmap -> store),
      pmap -> store.org, pmap -> store.size, 0
   ))
      return -1;
#endif

   return 0;
}



void OOC_InitFindPhotons (struct PhotonMap *pmap)
{
   if (OOC_InitNearest(&pmap -> squeue, pmap -> maxGather + 1,
                       sizeof(Photon)))
      error(SYSTEM, "can't allocate photon search queue");
}



static void OOC_InitPhotonCache (struct PhotonMap *pmap)
/* Initialise OOC photon cache */
{
   static char warn = 1;
   
   if (!pmap -> store.cache && !pmap -> numDensity) {
      if (pmapCacheSize > 0) {
         const unsigned pageSize = pmapCachePageSize * pmap -> maxGather,
                        numPages = pmapCacheSize / pageSize;
         /* Allocate & init I/O cache in octree */
         pmap -> store.cache = malloc(sizeof(OOC_Cache));
         
         if (!pmap -> store.cache || 
             OOC_CacheInit(pmap -> store.cache, numPages, pageSize, 
                           sizeof(Photon))) {
            error(SYSTEM, "failed OOC photon map cache init");
         }
      }
      else if (warn) {
         error(WARNING, "OOC photon map cache DISABLED");
         warn = 0;
      }
   }   
}



typedef struct {
   const PhotonMap   *pmap;
   const float       *norm;
} OOC_FilterData;



int OOC_FilterPhoton (void *p, void *fd)
/* Filter callback for photon kNN search, used by OOC_FindNearest() */
{
   const Photon         *photon = p;
   const OOC_FilterData *filtData = fd; 
   const PhotonMap      *pmap = filtData -> pmap;

   /* Reject photon if normal faces away (ignored for volume photons) with
    * tolerance to account for perturbation; note photon normal is coded
    * in range [-127,127], hence we factor this in */
   if (filtData -> norm && 
       DOT(filtData->norm, photon->norm) <= PMAP_NORM_TOL * 127 * frandom())
      return 0;
      
   if (isContribPmap(pmap)) {
      /* Lookup in contribution photon map; filter according to emitting
       * light source if contrib list set, else accept all */
       
      if (pmap -> srcContrib) {
         OBJREC *srcMod; 
         const int srcIdx = photonSrcIdx(pmap, photon);
      
         if (srcIdx < 0 || srcIdx >= nsources)
            error(INTERNAL, "invalid light source index in photon map");
      
         srcMod = findmaterial(source [srcIdx].so);

         /* Reject photon if contributions from light source which emitted
          * it are not sought */
         if (!lu_find(pmap -> srcContrib, srcMod -> oname) -> data)
            return 0;
      }

      /* Reject non-caustic photon if lookup for caustic contribs */
      if (pmap -> lookupCaustic && !photon -> caustic)
         return 0;
   }
   
   /* Accept photon */
   return 1;   
}



int OOC_FindPhotons (struct PhotonMap *pmap, const FVECT pos, const FVECT norm)
{
   OOC_SearchFilter  filt;
   OOC_FilterData    filtData;
   float             n [3];
   
   /* Lazily init OOC cache */
   if (!pmap -> store.cache)
      OOC_InitPhotonCache(pmap);
         
   /* Set up filter callback */
   filtData.pmap = pmap;
   if (norm)
      VCOPY(n, norm);
   filtData.norm = norm ? n : NULL;
   filt.data = &filtData;
   filt.func = OOC_FilterPhoton;

   pmap -> maxDist2 = OOC_FindNearest(&pmap -> store, 
                                      OOC_ROOT(&pmap -> store), 0, 
                                      pmap -> store.org, pmap -> store.size,
                                      pos, &filt, &pmap -> squeue, 
                                      pmap -> maxDist2);

   if (pmap -> maxDist2 < 0)
      error(INTERNAL, "failed k-NN photon lookup in OOC_FindPhotons");
   
   /* Return success or failure (empty queue => none found) */
   return pmap -> squeue.tail ? 0 : -1; 
}



int OOC_Find1Photon (struct PhotonMap* pmap, const FVECT pos, 
                     const FVECT norm, Photon *photon)
{
   OOC_SearchFilter     filt;
   OOC_FilterData       filtData;
   float                n [3], maxDist2;
   
   /* Lazily init OOC cache */
   if (!pmap -> store.cache)
      OOC_InitPhotonCache(pmap);
   
   /* Set up filter callback */
   filtData.pmap = pmap;
   if (norm)
      VCOPY(n, norm);
   filtData.norm = norm ? n : NULL;   
   filt.data = &filtData;
   filt.func = OOC_FilterPhoton;
   
   maxDist2 = OOC_Find1Nearest(&pmap -> store, 
                               OOC_ROOT(&pmap -> store), 0, 
                               pmap -> store.org, pmap -> store.size,
                               pos, &filt, photon, pmap -> maxDist2);
                               
   if (maxDist2 < 0)
      error(INTERNAL, "failed 1-NN photon lookup in OOC_Find1Photon");      
      
   if (maxDist2 >= pmap -> maxDist2)
      /* No photon found => failed */
      return -1;
   else {
      /* Set photon distance => success */
      pmap -> maxDist2 = maxDist2;
      return 0;
   }
}



int OOC_GetPhoton (struct PhotonMap *pmap, PhotonIdx idx, 
                   Photon *photon)
{
   return OOC_GetData(&pmap -> store, idx, photon);
}



Photon *OOC_GetNearestPhoton (const PhotonSearchQueue *squeue, PhotonIdx idx)
{
   return OOC_GetNearest(squeue, idx);
}



PhotonIdx OOC_FirstPhoton (const struct PhotonMap* pmap)
{
   return 0;
}
