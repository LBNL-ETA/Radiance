/* 
   ==================================================================
   Photon map interface to RADIANCE render options

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================   
   
   $Id: pmapopt.h,v 2.1 2015/02/24 19:39:27 greg Exp $
*/



int getPmapRenderOpt (int ac, char *av []);
/* Parse next render option for photon map; interface to getrenderopt();
 * return -1 if parsing failed, else number of parameters consumed */

void printPmapDefaults ();
/* Print defaults for photon map render options */
