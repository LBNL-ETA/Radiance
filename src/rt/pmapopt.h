/* RCSid $Id: pmapopt.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $ */

/* 
   ==================================================================
   Photon map interface to RADIANCE render options

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapopt.h,v 2.5 2016/05/17 17:39:47 rschregle Exp $
*/



#ifndef PMAPOPT_H
   #define PMAPOPT_H

   int getPmapRenderOpt (int ac, char *av []);
   /* Parse next render option for photon map; interface to getrenderopt();
    * return -1 if parsing failed, else number of parameters consumed */

   void printPmapDefaults ();
   /* Print defaults for photon map render options */
#endif
