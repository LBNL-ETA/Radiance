/* RCSid $Id: pmapopt.h,v 2.6 2024/09/17 16:36:05 greg Exp $ */

/* 
   ==================================================================
   Photon map interface to RADIANCE render options

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapopt.h,v 2.6 2024/09/17 16:36:05 greg Exp $
*/



#ifndef PMAPOPT_H
   #define PMAPOPT_H

#ifdef __cplusplus
extern "C" {
#endif

   int getPmapRenderOpt (int ac, char *av []);
   /* Parse next render option for photon map; interface to getrenderopt();
    * return -1 if parsing failed, else number of parameters consumed */

   void printPmapDefaults ();
   /* Print defaults for photon map render options */

#ifdef __cplusplus
}
#endif

#endif
