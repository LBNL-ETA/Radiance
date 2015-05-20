/* 
   ==================================================================
   Photon map data structures and kd-tree handling   

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapdata.h,v 2.3 2015/05/20 12:58:31 greg Exp $
*/


#ifndef PMAPDATA_H
   #define PMAPDATA_H

   #include "ray.h"   
   #include "pmaptype.h"
   #include "lookup.h"


   /* Define PMAP_FLOAT_FLUX to store photon flux as floats instead of
    * compact RGBE, which was found to improve accuracy in analytical
    * validation. */
   #ifdef PMAP_FLOAT_FLUX
      #define setPhotonFlux(p,f) copycolor((p) -> flux, f)
      #define getPhotonFlux(p,f) copycolor(f, (p) -> flux)
   #else
      #define setPhotonFlux(p,f) setcolr((p)->flux, (f)[0], (f)[1], (f)[2])
      #define getPhotonFlux(p,f) colr_color(f, (p) -> flux)
   #endif


   /* Primary photon ray for light source contributions */
   typedef struct {
      short srcIdx;               /* Index of emitting light source */
      float dir [3], pos [3];     /* Incident dir & hit point */
   } PhotonPrimary;
      
   #define photonSrcIdx(pm, p) ((pm) -> primary [(p) -> primary].srcIdx)

   
   typedef struct {
      float pos [3];                 /* Photon position */
      signed char norm [3];          /* Surface normal at pos (incident
                                        direction for volume photons) */
      char flags;                    /* Bit 0-1: kd-tree discriminator axis,
                                        Bit 2:   caustic photon */
      #ifdef PMAP_FLOAT_FLUX
         COLOR flux;
      #else
         COLR flux;                  /* Photon flux */
      #endif

      unsigned primary;              /* Index to primary ray */
   } Photon;
   
   /* Photon flag bitmasks and manipulation macros */
   #define PMAP_DISCR_BIT 3
   #define PMAP_CAUST_BIT 4
   #define photonDiscr(p)       ((p).flags & PMAP_DISCR_BIT)
   #define setPhotonDiscr(p, d) ((p).flags = ((p).flags & ~PMAP_DISCR_BIT) | \
                                              ((d) & PMAP_DISCR_BIT))

   /* Photon map type tests */
   #define isGlobalPmap(p)    ((p) -> type == PMAP_TYPE_GLOBAL)
   #define isCausticPmap(p)   ((p) -> type == PMAP_TYPE_CAUSTIC)
   #define isContribPmap(p)   ((p) -> type == PMAP_TYPE_CONTRIB)
   #define isVolumePmap(p)    ((p) -> type == PMAP_TYPE_VOLUME)

   
   /* Queue node for photon lookups */
   typedef struct {                  
      Photon *photon;
      float dist;
   } PhotonSQNode;

   /* Bias compensation history node */
   typedef struct {                  
      COLOR irrad;
      float weight;
   } PhotonBCNode;


   struct PhotonMap;                 /* Forward declarations */

   typedef struct PhotonMap {
      PhotonMapType type;            /* See pmaptype.h */
      char *fileName;                /* Photon map file */
      Photon *heap;                  /* Photon k-d tree as linear array */
      PhotonSQNode *squeue;          /* Search queue */
      PhotonBCNode *biasCompHist;    /* Bias compensation history */      
      char lookupFlags;              /* Flags passed to findPhotons() */
      
      unsigned long distribTarget,   /* Num stored specified by user */
                    heapSize,
                    heapSizeInc, 
                    heapEnd,         /* Num actually stored in heap */
                    numDensity,      /* Num density estimates */
                    totalGathered,   /* Total photons gathered */
                    numLookups,      /* Counters for short photon lookups */
                    numShortLookups;
                    
      unsigned minGather,            /* Specified min/max photons per */
               maxGather,            /* density estimate */
               squeueSize,
               squeueEnd,
               minGathered,          /* Min/max photons actually gathered */
               maxGathered,          /* per density estimate */
               shortLookupPct;       /* Percentage of short lookups (for
                                        statistics output */
               
                                     /* NOTE: All distances are SQUARED */
      float maxDist,                 /* Max search radius during NN search */
            maxDist0,                /* Initial value for above */
            maxDistLimit,            /* Hard limit for above */
            CoGdist,                 /* Avg distance to centre of gravity */
            maxPos [3], minPos [3],  /* Max & min photon positions */
            distribRatio,            /* Probability of photon storage */
            gatherTolerance,         /* Fractional deviation from minGather/
                                        maxGather for short lookup */
            minError, maxError,      /* Min/max/rms density estimate error */
            rmsError;
            
      FVECT CoG;                     /* Centre of gravity (avg photon pos) */
      COLOR photonFlux;              /* Average photon flux */
      unsigned short randState [3];  /* Local RNG state */
      
      void (*lookup)(struct PhotonMap*, RAY*, COLOR);   
                                     /* Callback for type-specific photon
                                      * lookup (usually density estimate) */

      PhotonPrimary *primary;        /* Primary photon rays & associated
                                        counters */
      unsigned primarySize, primaryEnd;
               
      LUTAB *srcContrib;             /* lookup table for source contribs */
   } PhotonMap;



   /* Photon maps by type (see PhotonMapType) */
   extern PhotonMap *photonMaps [];

   /* Macros for specific photon map types */
   #define globalPmap            (photonMaps [PMAP_TYPE_GLOBAL])
   #define preCompPmap           (photonMaps [PMAP_TYPE_PRECOMP])
   #define causticPmap           (photonMaps [PMAP_TYPE_CAUSTIC])
   #define volumePmap            (photonMaps [PMAP_TYPE_VOLUME])
   #define directPmap            (photonMaps [PMAP_TYPE_DIRECT])
   #define contribPmap           (photonMaps [PMAP_TYPE_CONTRIB])



   void initPhotonMap (PhotonMap *pmap, PhotonMapType t);
   /* Initialise empty photon map of specified type */

   const PhotonPrimary* addPhotonPrimary (PhotonMap *pmap, const RAY *ray);
   /* Add primary ray for emitted photon and save light source index, origin
    * on source, and emitted direction; used by contrib photons */

   const Photon* addPhoton (PhotonMap *pmap, const RAY *ray);
   /* Create new photon with ray's direction, intersection point, and flux,
      and add to photon map subject to acceptance probability pmap ->
      distribRatio for global density control; if the photon is rejected,
      NULL is returned. The flux is scaled by ray -> rweight and
      1 / pmap -> distribRatio. */

   void balancePhotons (PhotonMap *pmap, double *photonFlux);
   /* Build a balanced kd-tree as heap to guarantee logarithmic search
    * times. This must be called prior to performing photon search with
    * findPhotons(). 
    * PhotonFlux is the flux scaling factor per photon averaged over RGB. In
    * the case of a contribution photon map, this is an array with a
    * separate factor specific to each light source due to non-uniform
    * photon emission; Otherwise it is referenced as a scalar value. Flux is
    * not scaled if photonFlux == NULL. */

   void buildHeap (Photon *heap, unsigned long *heapIdx,
                   unsigned long *heapXdi, const float min [3], 
                   const float max [3], unsigned long left, 
                   unsigned long right, unsigned long root);
   /* Recursive part of balancePhotons(..); builds heap from subarray
    * defined by indices left and right. */
      
   void findPhotons (PhotonMap* pmap, const RAY *ray);
   /* Find pmap -> squeueSize closest photons to ray -> rop with similar 
      normal. For volume photons ray -> rorg is used and the normal is 
      ignored (being the incident direction in this case). Found photons 
      are placed in pmap -> squeue, with pmap -> squeueEnd being the number
      actually found. */

   Photon *find1Photon (PhotonMap *pmap, const RAY *ray);
   /* Finds single closest photon to ray -> rop with similar normal. 
      Returns NULL if none found. */

   void deletePhotons (PhotonMap*);
   /* Free dem mammaries... */

#endif
