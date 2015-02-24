#ifndef PMAPRAND_H
#define PMAPRAND_H

/* 
   ==================================================================
   Random number generators for photon distribution

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: pmaprand.h,v 2.1 2015/02/24 19:39:27 greg Exp $
*/



/* According to the analytical validation, skipping numbers in the sequence
   introduces bias in scenes with high reflectance. We therefore use
   erand48() with separate states for photon emission, scattering, and
   russian roulette. The pmapSeed() and pmapRandom() macros can be adapted
   to other (better?) RNGs. */   
   
#define pmapSeed(seed, state) (state [0] += seed, state [1] += seed, \
                               state [2] += seed)
#define pmapRandom(state) erand48(state)


   
extern unsigned short partState [3], emitState [3], cntState [3],
                      mediumState [3], scatterState [3], rouletteState [3],
                      randSeed;

#endif
