#ifndef lint
static const char RCSid[] = "$Id: pmapsrc.c,v 2.20 2021/04/09 17:42:34 rschregle Exp $";
#endif
/* 
   ======================================================================
   Photon map support routines for emission from light sources

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation (DFG) 
       under the FARESYS project.
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF #147053).
   (c) Tokyo University of Science,
       supported by the JSPS KAKENHI Grant Number JP19KK0115.
   ======================================================================

   $Id: pmapsrc.c,v 2.20 2021/04/09 17:42:34 rschregle Exp $"
*/



#include "pmapsrc.h"
#include "pmap.h"
#include "pmaprand.h"
#include "otypes.h"
#include "otspecial.h"



/* List of photon port modifier names */
char *photonPortList [MAXSET + 1] = {NULL};
/* Photon port objects (with modifiers in photonPortMods) */
SRCREC *photonPorts = NULL;
unsigned numPhotonPorts = 0;

void (*photonPartition [NUMOTYPE]) (EmissionMap*);
void (*photonOrigin [NUMOTYPE]) (EmissionMap*);



/* PHOTON PORT SUPPORT ROUTINES ------------------------------------------ */



/* Get/set photon port orientation flags from/in source flags.
 * HACK: the port orientation flags are embedded in the source flags and
 * shifted so they won't clobber the latter, since these are interpreted 
 * by the *PhotonPartition() and *PhotonOrigin() routines! */
#define PMAP_SETPORTFLAGS(portdir) ((int)(portdir) << 14)
#define PMAP_GETPORTFLAGS(sflags)  ((sflags) >> 14)

/* Set number of source partitions.
 * HACK: this is doubled if the source acts as a bidirectionally
 * emitting photon port, resulting in alternating front/backside partitions,
 * although essentially each partition is just used twice with opposing
 * normals. */
#define PMAP_SETNUMPARTITIONS(emap)   ( \
   (emap) -> numPartitions <<= ( \
      (emap) -> port && \
      PMAP_GETPORTFLAGS((emap) -> port -> sflags) == PMAP_PORTBI \
   ) \
)

/* Get current source partition and numer of partitions
 * HACK: halve the number partitions if the source acts as a bidrectionally
 * emitting photon port, since each partition is used twice with opposing 
 * normals. */
#define PMAP_GETNUMPARTITIONS(emap) (\
   (emap) -> numPartitions >> ( \
      (emap) -> port && \
      PMAP_GETPORTFLAGS((emap) -> port -> sflags) == PMAP_PORTBI \
   ) \
)
#define PMAP_GETPARTITION(emap)  ( \
   (emap) -> partitionCnt >> ( \
      (emap) -> port && \
      PMAP_GETPORTFLAGS((emap) -> port -> sflags) == PMAP_PORTBI \
   ) \
)



void getPhotonPorts (char **portList)
/* Find geometry declared as photon ports from modifiers in portList */
{
   OBJECT   i;
   OBJREC   *obj, *mat;
   int      mLen;
   char     **lp;   
   
   /* Init photon port objects */
   photonPorts = NULL;
   
   if (!portList [0])
      return;
   
   for (i = numPhotonPorts = 0; i < nobjects; i++) {
      obj = objptr(i);
      mat = findmaterial(obj);
      
      /* Check if object is a surface and NOT a light source (duh) and 
       * resolve its material (if any) via any aliases, then check for 
       * inclusion in modifier list; note use of strncmp() to ignore port 
       * flags */
      if (issurface(obj -> otype) && mat && !islight(mat -> otype)) {
         mLen = strlen(mat -> oname);
         for (lp = portList; *lp && strncmp(mat -> oname, *lp, mLen); lp++);
         
         if (*lp) {
            /* Add photon port */
            photonPorts = (SRCREC*)realloc(
               photonPorts, (numPhotonPorts + 1) * sizeof(SRCREC)
            );
            if (!photonPorts) 
               error(USER, "can't allocate photon ports");
            
            photonPorts [numPhotonPorts].so = obj;
            /* Extract port orientation flags and embed in source flags.
             * Note setsource() combines (i.e. ORs) these with the actual 
             * source flags below. */
            photonPorts [numPhotonPorts].sflags = 
               PMAP_SETPORTFLAGS((*lp) [mLen]);
         
            if (!sfun [obj -> otype].of || !sfun[obj -> otype].of -> setsrc)
               objerror(obj, USER, "illegal photon port");
         
            setsource(photonPorts + numPhotonPorts, obj);
            numPhotonPorts++;
         }
      }
   }
   if (!numPhotonPorts)
      error(USER, "no valid photon ports found");
}



static void setPhotonPortNormal (EmissionMap* emap)
/* Set normal for current photon port partition (if defined) based on its
 * orientation */
{

   int i, portFlags;
   
   if (emap -> port) {
      /* Extract photon port orientation flags, set surface normal as follows:
         -- Port oriented forwards --> flip surface normal to point outwards, 
            since normal points inwards per mkillum convention) 
         -- Port oriented backwards --> surface normal is NOT flipped, since 
            it already points inwards.
         -- Port is bidirectionally/bilaterally oriented --> flip normal based
            on the parity of the current partition emap -> partitionCnt. In 
            this case, photon emission alternates between port front/back 
            faces for consecutive partitions.
      */   
      portFlags = PMAP_GETPORTFLAGS(emap -> port -> sflags);

      if (
         portFlags == PMAP_PORTFWD || 
         portFlags == PMAP_PORTBI && !(emap -> partitionCnt & 1)
      )
         for (i = 0; i < 3; i++)     
            emap -> ws [i] = -emap -> ws [i];
   }
}



/* SOURCE / PHOTON PORT PARTITIONING ROUTINES----------------------------- */
   


static int flatPhotonPartition2 (
   EmissionMap* emap, unsigned long mp, FVECT cent, FVECT u, FVECT v, 
   double du2, double dv2
)
/* Recursive part of flatPhotonPartition(..) */
{
   FVECT newct, newax;
   unsigned long npl, npu;

   if (mp > emap -> maxPartitions) {
      /* Enlarge partition array */
      emap -> maxPartitions <<= 1;
      emap -> partitions = (unsigned char*)realloc(
         emap -> partitions, emap -> maxPartitions >> 1
      );

      if (!emap -> partitions) 
         error(USER, "can't allocate source partitions");

      memset(
         emap -> partitions + (emap -> maxPartitions >> 2), 0,
         emap -> maxPartitions >> 2
      );
   }
   
   if (du2 * dv2 <= 1) {                /* hit limit? */
      setpart(emap -> partitions, emap -> partitionCnt, S0);
      emap -> partitionCnt++;
      return 1;      
   }
   
   if (du2 > dv2) {                     /* subdivide in U */
      setpart(emap -> partitions, emap -> partitionCnt, SU);
      emap -> partitionCnt++;
      newax [0] = 0.5 * u [0];
      newax [1] = 0.5 * u [1];
      newax [2] = 0.5 * u [2];
      u = newax;
      du2 *= 0.25;
   } 
   
   else {                               /* subdivide in V */
      setpart(emap -> partitions, emap -> partitionCnt, SV);
      emap -> partitionCnt++;
      newax [0] = 0.5 * v [0];
      newax [1] = 0.5 * v [1];
      newax [2] = 0.5 * v [2];
      v = newax;
      dv2 *= 0.25;
   }
   
   /* lower half */
   newct [0] = cent [0] - newax [0];
   newct [1] = cent [1] - newax [1];
   newct [2] = cent [2] - newax [2];
   npl = flatPhotonPartition2(emap, mp << 1, newct, u, v, du2, dv2);
   
   /* upper half */
   newct [0] = cent [0] + newax [0];
   newct [1] = cent [1] + newax [1];
   newct [2] = cent [2] + newax [2];
   npu = flatPhotonPartition2(emap, mp << 1, newct, u, v, du2, dv2);
   
   /* return total */
   return npl + npu;
}



static void flatPhotonPartition (EmissionMap* emap)
/* Partition flat source for photon emission */
{
   RREAL    *vp;
   double   du2, dv2;

   memset(emap -> partitions, 0, emap -> maxPartitions >> 1);
   emap -> partArea = srcsizerat * thescene.cusize;
   emap -> partArea *= emap -> partArea;
   vp = emap -> src -> ss [SU];
   du2 = DOT(vp, vp) / emap -> partArea;
   vp = emap -> src -> ss [SV];
   dv2 = DOT(vp, vp) / emap -> partArea;
   emap -> partitionCnt = 0;
   emap -> numPartitions = flatPhotonPartition2(
      emap, 1, emap -> src -> sloc, 
      emap -> src -> ss [SU], emap -> src -> ss [SV], du2, dv2
   );
   emap -> partitionCnt = 0;
   emap -> partArea = emap -> src -> ss2 / emap -> numPartitions;
}



static void sourcePhotonPartition (EmissionMap* emap)
/* Partition scene cube faces or photon port for photon emission from 
   distant source */
{   
   if (emap -> port) {
      /* Relay partitioning to photon port */
      SRCREC *src = emap -> src;
      emap -> src = emap -> port;
      photonPartition [emap -> src -> so -> otype] (emap);
      PMAP_SETNUMPARTITIONS(emap);
      emap -> src = src;
   }
   
   else {
      /* No photon ports defined; partition scene cube faces (SUBOPTIMAL) */
      memset(emap -> partitions, 0, emap -> maxPartitions >> 1);
      setpart(emap -> partitions, 0, S0);
      emap -> partitionCnt = 0;
      emap -> numPartitions = 1 / srcsizerat;
      emap -> numPartitions *= emap -> numPartitions;
      emap -> partArea = sqr(thescene.cusize) / emap -> numPartitions;
      emap -> numPartitions *= 6;  
   }
}



static void spherePhotonPartition (EmissionMap* emap)
/* Partition spherical source into equal solid angles using uniform
   mapping for photon emission */
{
   unsigned numTheta, numPhi;
   
   memset(emap -> partitions, 0, emap -> maxPartitions >> 1);
   setpart(emap -> partitions, 0, S0);
   emap -> partArea = 4 * PI * sqr(emap -> src -> srad);
   emap -> numPartitions = 
      emap -> partArea / sqr(srcsizerat * thescene.cusize);

   numTheta = max(sqrt(2 * emap -> numPartitions / PI) + 0.5, 1);
   numPhi = 0.5 * PI * numTheta + 0.5;
   
   emap -> numPartitions = (unsigned long)numTheta * numPhi;
   emap -> partitionCnt = 0;
   emap -> partArea /= emap -> numPartitions;
}



static int cylPhotonPartition2 (
   EmissionMap* emap, unsigned long mp, FVECT cent, FVECT axis, double d2
)
/* Recursive part of cyPhotonPartition(..) */
{
   FVECT newct, newax;
   unsigned long npl, npu;

   if (mp > emap -> maxPartitions) {
      /* Enlarge partition array */
      emap -> maxPartitions <<= 1;
      emap -> partitions = (unsigned char*)realloc(
         emap -> partitions, emap -> maxPartitions >> 1
      );
      if (!emap -> partitions) 
         error(USER, "can't allocate source partitions");
         
      memset(
         emap -> partitions + (emap -> maxPartitions >> 2), 0,
         emap -> maxPartitions >> 2
      );
   }
   
   if (d2 <= 1) {
      /* hit limit? */
      setpart(emap -> partitions, emap -> partitionCnt, S0);
      emap -> partitionCnt++;
      return 1;
   }
   
   /* subdivide */
   setpart(emap -> partitions, emap -> partitionCnt, SU);
   emap -> partitionCnt++;
   newax [0] = 0.5 * axis [0];
   newax [1] = 0.5 * axis [1];
   newax [2] = 0.5 * axis [2];
   d2 *= 0.25;
   
   /* lower half */
   newct [0] = cent [0] - newax [0];
   newct [1] = cent [1] - newax [1];
   newct [2] = cent [2] - newax [2];
   npl = cylPhotonPartition2(emap, mp << 1, newct, newax, d2);
   
   /* upper half */
   newct [0] = cent [0] + newax [0];
   newct [1] = cent [1] + newax [1];
   newct [2] = cent [2] + newax [2];
   npu = cylPhotonPartition2(emap, mp << 1, newct, newax, d2);
   
   /* return total */
   return npl + npu;
}



static void cylPhotonPartition (EmissionMap* emap)
/* Partition cylindrical source for photon emission */
{
   double d2;

   memset(emap -> partitions, 0, emap -> maxPartitions >> 1);
   d2 = srcsizerat * thescene.cusize;
   d2 = PI * emap -> src -> ss2 / (2 * emap -> src -> srad * sqr(d2));
   d2 *= d2 * DOT(emap -> src -> ss [SU], emap -> src -> ss [SU]);

   emap -> partitionCnt = 0;
   emap -> numPartitions = cylPhotonPartition2(
      emap, 1, emap -> src -> sloc, emap -> src -> ss [SU], d2
   );
   emap -> partitionCnt = 0;
   emap -> partArea = PI * emap -> src -> ss2 / emap -> numPartitions;
}



/* PHOTON ORIGIN ROUTINES ------------------------------------------------ */



static void flatPhotonOrigin (EmissionMap* emap)
/* Init emission map with photon origin and associated surface axes on 
   flat light source surface. Also sets source aperture and sampling
   hemisphere axes at origin */
{
   int i, cent[3], size[3], parr[2];
   FVECT vpos;

   cent [0] = cent [1] = cent [2] = 0;
   size [0] = size [1] = size [2] = emap -> maxPartitions;
   parr [0] = 0; 
   parr [1] = PMAP_GETPARTITION(emap);
   
   if (!skipparts(cent, size, parr, emap -> partitions))
      error(CONSISTENCY, "bad source partition in flatPhotonOrigin");
      
   vpos [0] = (1 - 2 * pmapRandom(partState)) * size [0] / 
              emap -> maxPartitions;
   vpos [1] = (1 - 2 * pmapRandom(partState)) * size [1] / 
              emap -> maxPartitions;
   vpos [2] = 0;
   
   for (i = 0; i < 3; i++) 
      vpos [i] += (double)cent [i] / emap -> maxPartitions;
      
   /* Get origin */
   for (i = 0; i < 3; i++) 
      emap -> photonOrg [i] = 
         emap -> src -> sloc [i] + 
         vpos [SU] * emap -> src -> ss [SU][i] +
         vpos [SV] * emap -> src -> ss [SV][i] + 
         vpos [SW] * emap -> src -> ss [SW][i];
                              
   /* Get surface axes */
   VCOPY(emap -> us, emap -> src -> ss [SU]);
   normalize(emap -> us);
   VCOPY(emap -> ws, emap -> src -> ss [SW]);
   /* Flip normal emap -> ws if port and required by its orientation */
   setPhotonPortNormal(emap);
   fcross(emap -> vs, emap -> ws, emap -> us);
   
   /* Get hemisphere axes & aperture */
   if (emap -> src -> sflags & SSPOT) {
      VCOPY(emap -> wh, emap -> src -> sl.s -> aim);
      i = 0;
      
      do {
         emap -> vh [0] = emap -> vh [1] = emap -> vh [2] = 0;
         emap -> vh [i++] = 1;
         fcross(emap -> uh, emap -> vh, emap -> wh);
      } while (normalize(emap -> uh) < FTINY);
      
      fcross(emap -> vh, emap -> wh, emap -> uh);
      emap -> cosThetaMax = 1 - emap -> src -> sl.s -> siz / (2 * PI);
   }
   
   else {
      VCOPY(emap -> uh, emap -> us);
      VCOPY(emap -> vh, emap -> vs);
      VCOPY(emap -> wh, emap -> ws);
      emap -> cosThetaMax = 0;
   }
}



static void spherePhotonOrigin (EmissionMap* emap)
/* Init emission map with photon origin and associated surface axes on 
   spherical light source. Also sets source aperture and sampling
   hemisphere axes at origin */
{
   int i = 0;
   unsigned numTheta, numPhi, t, p;
   RREAL cosTheta, sinTheta, phi;

   /* Get current partition */
   numTheta = max(sqrt(2 * PMAP_GETNUMPARTITIONS(emap) / PI) + 0.5, 1);
   numPhi = 0.5 * PI * numTheta + 0.5;

   t = PMAP_GETPARTITION(emap) / numPhi;
   p = PMAP_GETPARTITION(emap) - t * numPhi;

   emap -> ws [2] = cosTheta = 1 - 2 * (t + pmapRandom(partState)) / numTheta;
   sinTheta = sqrt(1 - sqr(cosTheta));
   phi = 2 * PI * (p + pmapRandom(partState)) / numPhi;
   emap -> ws [0] = cos(phi) * sinTheta;
   emap -> ws [1] = sin(phi) * sinTheta;
   /* Flip normal emap -> ws if port and required by its orientation */
   setPhotonPortNormal(emap);

   /* Get surface axes us & vs perpendicular to ws */
   do {
      emap -> vs [0] = emap -> vs [1] = emap -> vs [2] = 0;
      emap -> vs [i++] = 1;
      fcross(emap -> us, emap -> vs, emap -> ws);
   } while (normalize(emap -> us) < FTINY);

   fcross(emap -> vs, emap -> ws, emap -> us);
   
   /* Get origin */
   for (i = 0; i < 3; i++)
      emap -> photonOrg [i] = emap -> src -> sloc [i] + 
                              emap -> src -> srad * emap -> ws [i];
                              
   /* Get hemisphere axes & aperture */
   if (emap -> src -> sflags & SSPOT) {
      VCOPY(emap -> wh, emap -> src -> sl.s -> aim);
      i = 0;

      do {
         emap -> vh [0] = emap -> vh [1] = emap -> vh [2] = 0;
         emap -> vh [i++] = 1;
         fcross(emap -> uh, emap -> vh, emap -> wh);
      } while (normalize(emap -> uh) < FTINY);

      fcross(emap -> vh, emap -> wh, emap -> uh);
      emap -> cosThetaMax = 1 - emap -> src -> sl.s -> siz / (2 * PI);
   }

   else {
      VCOPY(emap -> uh, emap -> us);
      VCOPY(emap -> vh, emap -> vs);
      VCOPY(emap -> wh, emap -> ws);
      emap -> cosThetaMax = 0;
   }
}



static void sourcePhotonOrigin (EmissionMap* emap)
/* Init emission map with photon origin and associated surface axes 
   on scene cube face for distant light source. Also sets source 
   aperture (solid angle) and sampling hemisphere axes at origin */
{  
   unsigned long i, partsPerDim, partsPerFace;
   unsigned face;
   RREAL du, dv;                            
      
   if (emap -> port) {
      /* Relay to photon port; get origin on its surface */
      SRCREC *src = emap -> src;
      emap -> src = emap -> port;
      photonOrigin [emap -> src -> so -> otype] (emap);
      emap -> src = src;
   }
   
   else {
      /* No ports defined, so get origin on scene cube face (SUBOPTIMAL) */
      /* Get current face from partition number */
      partsPerDim = 1 / srcsizerat;
      partsPerFace = sqr(partsPerDim);
      face = emap -> partitionCnt / partsPerFace;
      if (!(emap -> partitionCnt % partsPerFace)) {
         /* Skipped to a new face; get its normal */
         emap -> ws [0] = emap -> ws [1] = emap -> ws [2] = 0;
         emap -> ws [face >> 1] = face & 1 ? 1 : -1;
         
         /* Get surface axes us & vs perpendicular to ws */
         face >>= 1;
         emap -> vs [0] = emap -> vs [1] = emap -> vs [2] = 0;
         emap -> vs [(face + (emap -> ws [face] > 0 ? 2 : 1)) % 3] = 1;
         fcross(emap -> us, emap -> vs, emap -> ws);
      }
      
      /* Get jittered offsets within face from partition number 
         (in range [-0.5, 0.5]) */
      i = emap -> partitionCnt % partsPerFace;
      du = (i / partsPerDim + pmapRandom(partState)) / partsPerDim - 0.5;
      dv = (i % partsPerDim + pmapRandom(partState)) / partsPerDim - 0.5;

      /* Jittered destination point within partition */
      for (i = 0; i < 3; i++)
         emap -> photonOrg [i] = thescene.cuorg [i] + thescene.cusize * (
            0.5 + du * emap -> us [i] + dv * emap -> vs [i] +
            0.5 * emap -> ws [i]
         );
   }
   
   /* Get hemisphere axes & aperture */
   VCOPY(emap -> wh, emap -> src -> sloc);
   i = 0;
   
   do {
      emap -> vh [0] = emap -> vh [1] = emap -> vh [2] = 0;
      emap -> vh [i++] = 1;
      fcross(emap -> uh, emap -> vh, emap -> wh);
   } while (normalize(emap -> uh) < FTINY);
   
   fcross(emap -> vh, emap -> wh, emap -> uh);
   
   /* Get aperture */
   emap -> cosThetaMax = 1 - emap -> src -> ss2 / (2 * PI);
   emap -> cosThetaMax = min(1, max(-1, emap -> cosThetaMax));
}



static void cylPhotonOrigin (EmissionMap* emap)
/* Init emission map with photon origin and associated surface axes 
   on cylindrical light source surface. Also sets source aperture 
   and sampling hemisphere axes at origin */
{
   int i, cent[3], size[3], parr[2];
   FVECT v;

   cent [0] = cent [1] = cent [2] = 0;
   size [0] = size [1] = size [2] = emap -> maxPartitions;
   parr [0] = 0; 
   parr [1] = PMAP_GETPARTITION(emap);

   if (!skipparts(cent, size, parr, emap -> partitions))
      error(CONSISTENCY, "bad source partition in cylPhotonOrigin");

   v [SU] = 0;
   v [SV] = (1 - 2 * pmapRandom(partState)) * (double)size [1] / 
            emap -> maxPartitions;
   v [SW] = (1 - 2 * pmapRandom(partState)) * (double)size [2] / 
            emap -> maxPartitions;
   normalize(v);
   v [SU] = (1 - 2 * pmapRandom(partState)) * (double)size [1] / 
            emap -> maxPartitions; 

   for (i = 0; i < 3; i++) 
      v [i] += (double)cent [i] / emap -> maxPartitions;

   /* Get surface axes */
   for (i = 0; i < 3; i++) 
      emap -> photonOrg [i] = emap -> ws [i] = (
         v [SV] * emap -> src -> ss [SV][i] + 
         v [SW] * emap -> src -> ss [SW][i]
      ) / 0.8559;

   /* Flip normal emap -> ws if port and required by its orientation */
   setPhotonPortNormal(emap);

   normalize(emap -> ws);
   VCOPY(emap -> us, emap -> src -> ss [SU]);
   normalize(emap -> us);
   fcross(emap -> vs, emap -> ws, emap -> us);

   /* Get origin */
   for (i = 0; i < 3; i++) 
      emap -> photonOrg [i] += 
         v [SU] * emap -> src -> ss [SU][i] + emap -> src -> sloc [i];

   /* Get hemisphere axes & aperture */
   if (emap -> src -> sflags & SSPOT) {
      VCOPY(emap -> wh, emap -> src -> sl.s -> aim);
      i = 0;

      do {
         emap -> vh [0] = emap -> vh [1] = emap -> vh [2] = 0;
         emap -> vh [i++] = 1;
         fcross(emap -> uh, emap -> vh, emap -> wh);
      } while (normalize(emap -> uh) < FTINY);

      fcross(emap -> vh, emap -> wh, emap -> uh);
      emap -> cosThetaMax = 1 - emap -> src -> sl.s -> siz / (2 * PI);
   }

   else {
      VCOPY(emap -> uh, emap -> us);
      VCOPY(emap -> vh, emap -> vs);
      VCOPY(emap -> wh, emap -> ws);
      emap -> cosThetaMax = 0;
   }
}



/* PHOTON EMISSION ROUTINES ---------------------------------------------- */



static void defaultEmissionFunc (EmissionMap* emap)
/* Default behaviour when no emission funcs defined for this source type */
{
   objerror(emap -> src -> so, INTERNAL, 
            "undefined photon emission function");
}



void initPhotonEmissionFuncs ()
/* Init photonPartition[] and photonOrigin[] dispatch tables */
{
   int i;
   
   for (i = 0; i < NUMOTYPE; i++) 
      photonPartition [i] = photonOrigin [i] = defaultEmissionFunc;
   
   photonPartition [OBJ_FACE] = photonPartition [OBJ_RING] = flatPhotonPartition;
   photonPartition [OBJ_SOURCE] = sourcePhotonPartition;
   photonPartition [OBJ_SPHERE] = spherePhotonPartition;
   photonPartition [OBJ_CYLINDER] = cylPhotonPartition;
   photonOrigin [OBJ_FACE] = photonOrigin [OBJ_RING] = flatPhotonOrigin;
   photonOrigin [OBJ_SOURCE] = sourcePhotonOrigin;
   photonOrigin [OBJ_SPHERE] = spherePhotonOrigin;
   photonOrigin [OBJ_CYLINDER] = cylPhotonOrigin;
}      



void initPhotonEmission (EmissionMap *emap, float numPdfSamples)
/* Initialise photon emission from partitioned light source emap -> src;
 * this involves integrating the flux emitted from the current photon
 * origin emap -> photonOrg and setting up a PDF to sample the emission
 * distribution with numPdfSamples samples */
{
   unsigned i, t, p;
   double phi, cosTheta, sinTheta, du, dv, dOmega, thetaScale;
   EmissionSample* sample;
   const OBJREC* mod =  findmaterial(emap -> src -> so);
   static RAY r;   

   photonOrigin [emap -> src -> so -> otype] (emap);
   setcolor(emap -> partFlux, 0, 0, 0);   
   emap -> cdf = 0;
   emap -> numSamples = 0;
   cosTheta = DOT(emap -> ws, emap -> wh);         
   
   if (cosTheta <= 0 && sqrt(1-sqr(cosTheta)) <= emap -> cosThetaMax + FTINY)
      /* Aperture completely below surface; no emission from current origin */
      return;
   
   if (
      mod -> omod == OVOID && !emap -> port && (
         cosTheta >= 1 - FTINY || (
            emap -> src -> sflags & SDISTANT && 
            acos(cosTheta) + acos(emap -> cosThetaMax) <= 0.5 * PI
         )
      )
   ) {
      /* Source is unmodified and has no port (which requires testing for 
         occlusion), and is either local with normal aligned aperture or 
         distant with aperture above surface
         --> get flux & PDF via analytical solution */
      setcolor(
         emap -> partFlux, mod -> oargs.farg [0], mod -> oargs.farg [1], 
         mod -> oargs.farg [2]
      );

      /* Multiply radiance by projected Omega * dA to get flux */
      scalecolor(
         emap -> partFlux, 
         PI * cosTheta * (1 - sqr(max(emap -> cosThetaMax, 0))) * 
            emap -> partArea
      );
   }
   
   else {
      /* Source is either modified, has a port, is local with off-normal 
         aperture, or distant with aperture partly below surface
         --> get flux via numerical integration */
      thetaScale = (1 - emap -> cosThetaMax);
      
      /* Figga out numba of aperture samples for integration;
         numTheta / numPhi ratio is 1 / PI */
      du = sqrt(pdfSamples * 2 * thetaScale);
      emap -> numTheta = max(du + 0.5, 1);
      emap -> numPhi = max(PI * du + 0.5, 1);

      dOmega = 2 * PI * thetaScale / (emap -> numTheta * emap -> numPhi);
      thetaScale /= emap -> numTheta;
      
      /* Allocate PDF, baby */
      sample = emap -> samples = (EmissionSample*)realloc(
         emap -> samples, 
         sizeof(EmissionSample) * emap -> numTheta * emap -> numPhi
      );
      if (!emap -> samples) 
         error(USER, "can't allocate emission PDF");
         
      VCOPY(r.rorg, emap -> photonOrg);
      VCOPY(r.rop, emap -> photonOrg);
      r.rmax = 0;
      
      for (t = 0; t < emap -> numTheta; t++) {
         for (p = 0; p < emap -> numPhi; p++) {
            /* This uniform mapping handles 0 <= thetaMax <= PI */
            cosTheta = 1 - (t + pmapRandom(emitState)) * thetaScale;
            sinTheta = sqrt(1 - sqr(cosTheta));
            phi = 2 * PI * (p + pmapRandom(emitState)) / emap -> numPhi;
            du = cos(phi) * sinTheta;
            dv = sin(phi) * sinTheta;
            rayorigin(&r, PRIMARY, NULL, NULL);
            
            for (i = 0; i < 3; i++)
               r.rdir [i] = (
                  du * emap -> uh [i] + dv * emap -> vh [i] +
                  cosTheta * emap -> wh [i]
               );            
                            
            /* Sample behind surface? */
            VCOPY(r.ron, emap -> ws);
            if ((r.rod = DOT(r.rdir, r.ron)) <= 0) 
               continue;
               
            /* Get radiance emitted in this direction; to get flux we 
               multiply by cos(theta_surface), dOmega, and dA. Ray 
               is directed towards light source for raytexture(). */
            if (!(emap -> src -> sflags & SDISTANT)) 
               for (i = 0; i < 3; i++) 
                  r.rdir [i] = -r.rdir [i];
                  
            /* Port occluded in this direction? */
            if (emap -> port && localhit(&r, &thescene)) 
               continue; 
               
            raytexture(&r, mod -> omod);
            setcolor(
               r.rcol, mod -> oargs.farg [0], mod -> oargs.farg [1], 
               mod -> oargs.farg [2]
            );
            multcolor(r.rcol, r.pcol);
            
            /* Multiply by cos(theta_surface) */
            scalecolor(r.rcol, r.rod);
            
            /* Add PDF sample if nonzero; importance info for photon emission
             * could go here... */
            if (colorAvg(r.rcol)) {
               copycolor(sample -> pdf, r.rcol);
               sample -> cdf = emap -> cdf += colorAvg(sample -> pdf);
               sample -> theta = t;
               sample++ -> phi = p;
               emap -> numSamples++;
               addcolor(emap -> partFlux, r.rcol);
            }
         }
      }
      
      /* Multiply by dOmega * dA */
      scalecolor(emap -> partFlux, dOmega * emap -> partArea);
   }
}



#define vomitPhoton     emitPhoton
#define bluarrrghPhoton vomitPhoton

void emitPhoton (const EmissionMap* emap, RAY* ray)
/* Emit photon from current partition emap -> partitionCnt based on
   emission distribution. Returns new photon ray. */
{
   unsigned long i, lo, hi;
   const EmissionSample* sample = emap -> samples;
   RREAL du, dv, cosTheta, cosThetaSqr, sinTheta, phi;  
   const OBJREC* mod = findmaterial(emap -> src -> so);
   
   /* Choose a new origin within current partition for every 
      emitted photon to break up clustering artifacts */
   photonOrigin [emap -> src -> so -> otype] ((EmissionMap*)emap);
   /* If we have a local glow source with a maximum radius, then
      restrict our photon to the specified distance, otherwise we set
      the limit imposed by photonMaxDist (or no limit if 0) */
   if (
      mod -> otype == MAT_GLOW && 
      !(emap -> src -> sflags & SDISTANT) && mod -> oargs.farg[3] > FTINY
   )
      ray -> rmax = mod -> oargs.farg[3];
   else
      ray -> rmax = photonMaxDist;
   rayorigin(ray, PRIMARY, NULL, NULL);
   
   if (!emap -> numSamples) {
      /* Source is unmodified and has no port, and either local with 
         normal aligned aperture, or distant with aperture above surface
         --> use cosine weighted distribution */
      cosThetaSqr = (1 - 
         pmapRandom(emitState) * (1 - sqr(max(emap -> cosThetaMax, 0)))
      );
      cosTheta = sqrt(cosThetaSqr);
      sinTheta = sqrt(1 - cosThetaSqr);
      phi = 2 * PI * pmapRandom(emitState);
      setcolor(
         ray -> rcol, mod -> oargs.farg [0], mod -> oargs.farg [1], 
         mod -> oargs.farg [2]
      );
   }
   
   else {
      /* Source is either modified, has a port, is local with off-normal 
         aperture, or distant with aperture partly below surface
         --> choose  direction from constructed cumulative distribution 
         function with Monte Carlo inversion using binary search. */
      du = pmapRandom(emitState) * emap -> cdf;
      lo = 1;
      hi = emap -> numSamples;
      
      while (hi > lo) {
         i = (lo + hi) >> 1;
         sample = emap -> samples + i - 1;
         
         if (sample -> cdf >= du) 
            hi = i;
         if (sample -> cdf < du) 
            lo = i + 1;
      }

      /* Finalise found sample */
      i = (lo + hi) >> 1;
      sample = emap -> samples + i - 1;

      /* This is a uniform mapping, mon */
      cosTheta = (1 - 
         (sample -> theta + pmapRandom(emitState)) *
         (1 - emap -> cosThetaMax) / emap -> numTheta
      );
      sinTheta = sqrt(1 - sqr(cosTheta));
      phi = 2 * PI * (sample -> phi + pmapRandom(emitState)) / emap -> numPhi;
      copycolor(ray -> rcol, sample -> pdf);
   }

   /* Normalize photon flux so that average over RGB is 1 */
   colorNorm(ray -> rcol);   
   
   VCOPY(ray -> rorg, emap -> photonOrg);   
   du = cos(phi) * sinTheta;
   dv = sin(phi) * sinTheta;
   
   for (i = 0; i < 3; i++)
      ray -> rdir [i] = du * emap -> uh [i] + dv * emap -> vh [i] + 
                        cosTheta * emap -> wh [i];
                        
   if (emap -> src -> sflags & SDISTANT)
      /* Distant source; reverse ray direction to point into the scene. */     
      for (i = 0; i < 3; i++) 
         ray -> rdir [i] = -ray -> rdir [i];
      
   if (emap -> port)
      /* Photon emitted from port; move origin just behind port so it 
         will be scattered */
      for (i = 0; i < 3; i++) 
         ray -> rorg [i] -= 2 * FTINY * ray -> rdir [i];    
   
   /* Assign emitting light source index */
   ray -> rsrc = emap -> src - source;
}



/* SOURCE CONTRIBS FROM DIRECT / VOLUME PHOTONS -------------------------- */



void multDirectPmap (RAY *r)
/* Factor irradiance from direct photons into r -> rcol; interface to
 * direct() */
{
   COLOR photonIrrad;
   
   /* Lookup direct photon irradiance */
   (directPmap -> lookup)(directPmap, r, photonIrrad);
   
   /* Multiply with coefficient in ray */
   multcolor(r -> rcol, photonIrrad);
   
   return;
}



void inscatterVolumePmap (RAY *r, COLOR inscatter)
/* Add inscattering from volume photon map; interface to srcscatter() */
{
   /* Add ambient in-scattering via lookup callback */
   (volumePmap -> lookup)(volumePmap, r, inscatter);
}

