/* 
   ==================================================================
   Photon map support routines for emission from light sources

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmapsrc.h,v 2.2 2015/05/08 13:20:23 rschregle Exp $
*/



#ifndef PMAPSRC_H
   #define PMAPSRC_H

   #include "ray.h"
   #include "source.h"



   /* Data structures for photon emission */
   typedef struct {
      unsigned theta, phi;               /* Direction indices */
      COLOR pdf;                         /* Probability of emission in this
                                            direction per RGB */
      float cdf;                         /* Cumulative probability up to
                                            this sample */
   } EmissionSample;

   typedef struct {
      unsigned numTheta, numPhi;         /* Num divisions in theta & phi */
      RREAL cosThetaMax;                 /* cos(source aperture) */
      FVECT uh, vh, wh;                  /* Emission aperture axes at origin, 
                                            w is normal*/
      FVECT us, vs, ws;                  /* Source surface axes at origin, 
                                            w is normal */
      FVECT photonOrg;                   /* Current photon origin */
      EmissionSample *samples;           /* PDF samples for photon emission 
                                            directions */
      unsigned long numPartitions;       /* Number of surface partitions */
      RREAL partArea;                    /* Area covered by each partition */
      SRCREC *src, *port;                /* Current source and port */
      unsigned long partitionCnt,        /* Current surface partition */
                    maxPartitions,       /* Max allocated partitions */
                    numSamples;          /* Number of PDF samples */
      unsigned char* partitions;         /* Source surface partitions */
      COLOR partFlux;                    /* Total flux emitted by partition */
      float cdf;                         /* Cumulative probability */
   } EmissionMap;



   /* Photon port list for emission from geometry en lieu of light sources */
   extern SRCREC *photonPorts;
   extern unsigned numPhotonPorts;

   /* Dispatch tables for partitioning light source surfaces and generating
      an origin for photon emission */
   extern void (*photonPartition []) (EmissionMap*);
   extern void (*photonOrigin []) (EmissionMap*);



   void getPhotonPorts ();
   /* Find geometry declared as photon ports */

   void initPhotonEmissionFuncs ();
   /* Init photonPartition[] and photonOrigin[] dispatch tables with source
      specific emission routines */

   void initPhotonEmission (EmissionMap *emap, float numPdfSamples);
   /* Initialize photon emission from partitioned light source emap -> src;
    * this involves integrating the flux emitted from the current photon
    * origin emap -> photonOrg and setting up a PDF to sample the emission
    * distribution with numPdfSamples samples */

   void emitPhoton (const EmissionMap*, RAY* ray);
   /* Emit photon from current source partition based on emission 
      distribution and return new photon ray */

   void multDirectPmap (RAY *r);
   /* Factor irradiance from direct photons into r -> rcol; interface to
    * direct() */

   void inscatterVolumePmap (RAY *r, COLOR inscatter);
   /* Add inscattering from volume photon map; interface to srcscatter() */

#endif
