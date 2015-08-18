/* RCSid $Id: pmapopt.h,v 2.3 2015/08/18 18:45:55 greg Exp $ */
/* 
   ==================================================================
   Photon map interface to RADIANCE render options

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapopt.h,v 2.3 2015/08/18 18:45:55 greg Exp $
*/



int getPmapRenderOpt (int ac, char *av []);
/* Parse next render option for photon map; interface to getrenderopt();
 * return -1 if parsing failed, else number of parameters consumed */

void printPmapDefaults ();
/* Print defaults for photon map render options */
