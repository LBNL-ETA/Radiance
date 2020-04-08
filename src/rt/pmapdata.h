/* RCSid $Id: pmapdata.h,v 2.14 2020/04/08 15:14:21 rschregle Exp $ */

/* 
   =========================================================================
   Photon map types and interface to nearest neighbour lookups in underlying
   point cloud data structure.

   The default data structure is an in-core kd-tree (see pmapkdt.{h,c}).
   This can be overriden with the PMAP_OOC compiletime switch, which enables
   an out-of-core octree (see oococt.{h,c}).
   
   Defining PMAP_FLOAT_FLUX stores photon flux as floats rather than packed
   RGBE for greater precision; this may be necessary when the flux differs
   significantly in individual colour channels, e.g. with highly saturated
   colours.

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   =========================================================================
   
   $Id: pmapdata.h,v 2.14 2020/04/08 15:14:21 rschregle Exp $
*/



#ifndef PMAPDATA_H
   #define PMAPDATA_H

   #ifndef NIX
      #if defined(_WIN32) || defined(_WIN64)
         #define NIX 0
      #else
         #define NIX 1
      #endif
   #endif

   #if (defined(PMAP_OOC) && !NIX)
      #error "OOC currently only supported on NIX -- tuff luck."
   #endif

   #ifdef PMAP_CBDM
      /* Enable photon primary hitpoints and incident directions (see struct
       * PhotonPrimary below).  Note this will increase the size of photon
       * primaries 9-fold (10-fold after alignment)!!! */
      #define  PMAP_PRIMARYPOS
      #define  PMAP_PRIMARYDIR
   #endif

   #include "ray.h"
   #include "pmaptype.h"
   #include "paths.h"
   #include "lookup.h"
   #include <stdint.h>


   
   /* Primary photon ray for light source contributions */
   typedef struct {
      int16    srcIdx;              /* Index of emitting light source */
                                    /* !!! REDUCED FROM 32 BITS !!! */
#ifdef PMAP_PRIMARYDIR
      int32    dir;                 /* Encoded ray direction */
#endif
#ifdef PMAP_PRIMARYPOS
      float    pos [3];             /* Hit point */
#endif
   } PhotonPrimary;
      
   #define photonSrcIdx(pm, p)      ((pm)->primaries[(p)->primary].srcIdx)



   /* Photon primary ray index type and limit */
   typedef  uint32            PhotonPrimaryIdx;
   #define  PMAP_MAXPRIMARY   UINT32_MAX

   /* Macros for photon's generating subprocess field */
#ifdef PMAP_OOC
   #define  PMAP_PROCBITS  7
#else
   #define  PMAP_PROCBITS  5
#endif
   #define  PMAP_MAXPROC         (1 << PMAP_PROCBITS)
   #define  PMAP_GETRAYPROC(r)   ((r) -> crtype >> 8)
   #define  PMAP_SETRAYPROC(r,p) ((r) -> crtype |= p << 8)

   /* Tolerance for photon normal check during lookups */
   #define  PMAP_NORM_TOL  0.02
   

   
   typedef struct {
      float                pos [3];       /* Photon position */
      signed char          norm [3];      /* Surface normal at pos (incident
                                             direction for volume photons) */
      union {
         struct {
#ifndef PMAP_OOC
            unsigned char  discr    : 2;  /* kd-tree discriminator axis */
#endif            
            unsigned char  caustic  : 1;  /* Specularly scattered (=caustic) */
            
            /* Photon's generating subprocess index, used for primary ray
             * index linearisation when building contrib pmaps; note this is
             * reduced for kd-tree to accommodate the discriminator field */
            unsigned char  proc  : PMAP_PROCBITS; 
         };
         
         unsigned char     flags;
      };
      
#ifdef PMAP_FLOAT_FLUX
      COLOR                flux;
#else
      COLR                 flux;
#endif
      PhotonPrimaryIdx     primary;       /* Index to primary ray */
   } Photon;


   
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

   
   /* Bias compensation history node */
   typedef struct {
      COLOR irrad;
      float weight;
   } PhotonBiasCompNode;


   /* Forward declaration */
   struct PhotonMap;
   
   
   /* Define search queue and underlying data struct types */
#ifdef PMAP_OOC
   #include "pmapooc.h"
#else
   #include "pmapkdt.h"
#endif


   /* Mean size of heapfile write buffer, in number of photons */
   #define PMAP_HEAPBUFSIZE   1e6
   
   /* Mean idle time between heap locking attempts, in usec */
   #define PMAP_HEAPBUFSLEEP  2e6
   
   /* Temporary filename for photon heaps */
   #define PMAP_TMPFNAME      TEMPLATE
   #define PMAP_TMPFNLEN      (TEMPLEN + 1)


   typedef struct PhotonMap {
      PhotonMapType  type;             /* See pmaptype.h */
      char           *fileName;        /* Photon map file */
      
      
      /* ================================================================
       * PRE/POST-BUILD STORAGE
       * ================================================================ */
      FILE           *heap;            /* Unsorted photon heap prior to
                                          construction of store */
      char           heapFname [sizeof(PMAP_TMPFNAME)];       
      Photon         *heapBuf;         /* Write buffer for above */
      unsigned long  heapBufLen,       /* Current & max size of heapBuf */
                     heapBufSize;
      PhotonStorage  store;            /* Photon storage in space
                                          subdividing data struct */
       
      /* ================================================================
       * PHOTON DISTRIBUTION STUFF
       * ================================================================ */
      unsigned long  distribTarget,    /* Num stored specified by user */
                     numPhotons;       /* Num actually stored */
      float          distribRatio;     /* Probability of photon storage */
      COLOR          photonFlux;       /* Average photon flux */
      unsigned short randState [3];    /* Local RNG state */


      /* ================================================================
       * PHOTON LOOKUP STUFF 
       * ================================================================ */
      union {                          /* Flags passed to findPhotons() */
         char        lookupCaustic : 1;
         char        lookupFlags;
      };
      
      PhotonSearchQueue squeue;        /* Search queue for photon lookups */
      unsigned       minGather,        /* Specified min/max photons per */
                     maxGather;        /* density estimate */
                                       
                                       /* NOTE: All distances are SQUARED */
      float          maxDist2,         /* Max search radius */
                     maxDist0,         /* Initial value for above */
                     maxDist2Limit,    /* Hard limit for above */
                     gatherTolerance;  /* Fractional deviation from minGather/
                                          maxGather for short lookup */
      void (*lookup)(struct PhotonMap*, 
                     RAY*, COLOR);     /* Callback for type-specific photon
                                        * lookup (usually density estimate) */                                          

                     
      /* ================================================================
       * BIAS COMPENSATION STUFF
       * ================================================================ */
      PhotonBiasCompNode   *biasCompHist;    /* Bias compensation history */
      

      /* ================================================================
       * STATISTIX
       * ================================================================ */
      unsigned long  totalGathered,    /* Total photons gathered */
                     numDensity,       /* Num density estimates */
                     numLookups,       /* Counters for short photon lookups */
                     numShortLookups;
                     
      unsigned       minGathered,      /* Min/max photons actually gathered */
                     maxGathered,      /* per density estimate */
                     shortLookupPct;   /* % of short lookups for stats */                                          
                     
      float          minError,         /* Min/max/rms density estimate error */
                     maxError,
                     rmsError,
                     CoGdist,          /* Avg distance to centre of gravity */
                     maxPos [3],       /* Max & min photon positions */
                     minPos [3];
                     
      FVECT          CoG;              /* Centre of gravity (avg photon pos) */
                     
       
      /* ================================================================
       * PHOTON CONTRIB/COEFF STUFF
       * ================================================================ */
      PhotonPrimary  *primaries,    /* Photon primary array for rendering */
                     lastPrimary;   /* Current primary for photon distrib */
      PhotonPrimaryIdx  numPrimary;    /* Number of primary rays */
      LUTAB             *srcContrib;   /* Lookup table for source contribs */
   } PhotonMap;



   /* Photon maps by type (see PhotonMapType) */
   extern PhotonMap  *photonMaps [];

   /* Macros for specific photon map types */
   #define globalPmap         (photonMaps [PMAP_TYPE_GLOBAL])
   #define preCompPmap        (photonMaps [PMAP_TYPE_PRECOMP])
   #define causticPmap        (photonMaps [PMAP_TYPE_CAUSTIC])
   #define volumePmap         (photonMaps [PMAP_TYPE_VOLUME])
   #define directPmap         (photonMaps [PMAP_TYPE_DIRECT])
   #define contribPmap        (photonMaps [PMAP_TYPE_CONTRIB])

   /* Photon map type tests */
   #define isGlobalPmap(p)    ((p) -> type == PMAP_TYPE_GLOBAL)
   #define isCausticPmap(p)   ((p) -> type == PMAP_TYPE_CAUSTIC)
   #define isContribPmap(p)   ((p) -> type == PMAP_TYPE_CONTRIB)
   #define isVolumePmap(p)    ((p) -> type == PMAP_TYPE_VOLUME)



   void initPhotonMap (PhotonMap *pmap, PhotonMapType t);
   /* Initialise empty photon map of specified type */

   int newPhoton (PhotonMap *pmap, const RAY *ray);
   /* Create new photon with ray's direction, intersection point, and flux,
      and append to unsorted photon heap pmap -> heap. The photon is
      accepted with probability pmap -> distribRatio for global density
      control; if the photon is rejected, -1 is returned, else 0.  The flux
      is scaled by ray -> rweight and 1 / pmap -> distribRatio. */

   void initPhotonHeap (PhotonMap *pmap);
   /* Open photon heap file */

   void flushPhotonHeap (PhotonMap *pmap);
   /* Flush photon heap buffa pmap -> heapBuf to heap file pmap -> heap;
    * used by newPhoton() and to finalise heap in distribPhotons(). */

   void buildPhotonMap (PhotonMap *pmap, double *photonFlux,
                        PhotonPrimaryIdx *primaryOfs, unsigned nproc);
   /* Postpress unsorted photon heap pmap -> heap and build underlying data
    * structure pmap -> store.  This is prerequisite to photon lookups with
    * findPhotons(). */
    
   /* PhotonFlux is the flux per photon averaged over RGB; this is
    * multiplied with each photon's flux during the postprocess.  In the
    * case of a contribution photon map, this is an array with a separate
    * flux specific to each light source due to non-uniform photon emission;
    * Otherwise it is referenced as a scalar value.  Flux is not scaled if
    * photonFlux == NULL.  */

   /* Photon map construction may be parallelised if nproc > 1, if
    * supported. The heap is destroyed on return.  */
    
   /* PrimaryOfs is an array of index offsets for the primaries in pmap ->
    * primaries generated by each of the nproc subprocesses during contrib
    * photon distribution (see distribPhotonContrib()).  These offsets are
    * used to linearise the photon primary indices in the postprocess.  This
    * linearisation is skipped if primaryOfs == NULL.  */

   void findPhotons (PhotonMap* pmap, const RAY *ray);
   /* Find pmap -> squeue.len closest photons to ray -> rop with similar 
      normal. For volume photons ray -> rorg is used and the normal is 
      ignored (being the incident direction in this case). Found photons 
      are placed search queue starting with the furthest photon at pmap ->
      squeue.node, and pmap -> squeue.tail being the number actually found. */

   Photon *find1Photon (PhotonMap *pmap, const RAY *ray, Photon *photon);
   /* Find single closest photon to ray -> rop with similar normal. 
      Return NULL if none found, else the supplied Photon* buffer, 
      indicating that it contains a valid photon. */

   void getPhoton (PhotonMap *pmap, PhotonIdx idx, Photon *photon);
   /* Retrieve photon referenced by idx from pmap -> store */

   Photon *getNearestPhoton (const PhotonSearchQueue *squeue, PhotonIdx idx);
   /* Retrieve photon from NN search queue after calling findPhotons() */
   
   PhotonIdx firstPhoton (const PhotonMap *pmap);
   /* Index to first photon, to be passed to getPhoton(). Indices to
    * subsequent photons can be optained via increment operator (++) */
    
   void deletePhotons (PhotonMap*);
   /* Free dem mammaries... */

#endif
