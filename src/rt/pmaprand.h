/* RCSid $Id: pmaprand.h,v 2.7 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ==================================================================
   Random number generators for photon distribution

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
   $Id: pmaprand.h,v 2.7 2016/05/17 17:39:47 rschregle Exp $
*/



#ifndef PMAPRAND_H
   #define PMAPRAND_H

   /* According to the analytical validation, skipping numbers in the sequence
      introduces bias in scenes with high reflectance. We therefore use
      erand48() with separate states for photon emission, scattering, and
      russian roulette.  The pmapSeed() and pmapRandom() macros can be
      adapted to other (better?) RNGs.  */

#if defined(_WIN32) || defined(BSD)
   /* Assume no erand48(), so use standard RNG without explicit multistate 
      control; the resulting sequences will be suboptimal */      
   #include "random.h"
      
   #define pmapSeed(seed, state) (srandom(seed))
   #define pmapRandom(state)     (frandom())
#else
   #define pmapSeed(seed, state) (state [0] += seed, state [1] += seed, \
                                  state [2] += seed)
   #define pmapRandom(state) erand48(state)
#endif

      
   extern unsigned short partState [3], emitState [3], cntState [3],
                         mediumState [3], scatterState [3], rouletteState [3],
                         randSeed;
#endif

