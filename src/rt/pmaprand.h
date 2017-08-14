/* RCSid $Id: pmaprand.h,v 2.8 2017/08/14 21:12:10 rschregle Exp $ */

/* 
   ======================================================================
   Random number generators for photon distribution

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation (SNSF, #147053)
   ======================================================================
   
   $Id: pmaprand.h,v 2.8 2017/08/14 21:12:10 rschregle Exp $
*/



#ifndef PMAPRAND_H
   #define PMAPRAND_H

   /* According to the analytical validation, skipping sequential samples
    * when sharing a single RNG among multiple sampling domains introduces
    * overlapping photon rays (reported as 'duplicate keys' when building
    * the underlying data structure) and therefore bias.  This is aggravated
    * when running parallel photon distribution subprocesses, where photon
    * rays from different subprocesses may correlate.    
    * We therefore maintain a separate RNG state for each sampling domain
    * (e.g. photon emission, scattering, and russian roulette).  With
    * multiprocessing, each subprocess has its own instance of the RNG
    * state, which is independently seeded for decorellation -- see
    * distribPhotons() and distribPhotonContrib().    
    * The pmapSeed() and pmapRandom() macros below can be adapted to
    * platform-specific RNGs if necessary.  */
#if defined(_WIN32) || defined(_WIN64)
   /* Use standard RNG without state management; the generated sequences
    * will be suboptimal */
   #include "random.h"      
   #define pmapSeed(seed, state) srandom(seed)
   #define pmapRandom(state)     frandom()
#else
   /* Assume NIX and manage RNG state via erand48() */
   #define pmapSeed(seed, state) (state [0] += seed, state [1] += seed, \
                                  state [2] += seed)
   #define pmapRandom(state)     erand48(state)
#endif

      
   extern unsigned short partState [3], emitState [3], cntState [3],
                         mediumState [3], scatterState [3], rouletteState [3],
                         randSeed;
#endif

