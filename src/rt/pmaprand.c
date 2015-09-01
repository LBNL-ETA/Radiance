#ifndef lint
static const char RCSid[] = "$Id: pmaprand.c,v 2.4 2015/09/01 16:27:52 greg Exp $";
#endif
/* 
   ==================================================================
   Random number generators for photon distribution

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================
   
*/



/* 
   Separate RNG states are used for the following variates during photon 
   distribution:
   
   - source partition
   - emission direction
   - emission counter fraction
   - mean free distance in medium
   - scattering direction
   - russian roulette.
   
   Each photon map also has a local state randState used for distribRatio 
   during distribution and for bias compensation during gathering. 
   
   The random seed randSeed can added to each initial state with pmapSeed()
   so the RNGs can be externally seeded if necessary.
*/
   
unsigned short partState [3] = {47717, 5519, 21521},
               emitState [3] = {33997, 59693, 11003},
               cntState [3] = {17077, 4111, 48907},
               mediumState [3] = {25247, 7507, 33797},
               scatterState [3] = {21863, 45191, 5099},
               rouletteState [3] = {10243, 39829, 9433},
               randSeed = 0;
