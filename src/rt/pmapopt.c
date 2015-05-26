/* 
   ==================================================================
   Photon map interface to RADIANCE render options

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
   (c) Lucerne University of Applied Sciences and Arts,
   supported by the Swiss National Science Foundation (SNSF, #147053)
   ==================================================================   
   
   $Id: pmapopt.c,v 2.5 2015/05/26 13:31:19 rschregle Exp $
*/



#include "ray.h"
#include "pmapparm.h"



int getPmapRenderOpt (int ac, char *av [])
/* Parse next render option for photon map; interface to getrenderopt();
 * return -1 if parsing failed, else number of parameters consumed */
{
   #define check(ol,al) (av[0][ol] || badarg(ac-1,av+1,al))
   
   static int t = -1;	/* pmap parameter index */
                        
   if (ac < 1 || !av [0] || av [0][0] != '-')
      return -1;
       
   switch (av [0][1]) {
      case 'a':
         switch (av [0][2]) {
            case 'p': /* photon map */
	       /* Asking for photon map, ergo ambounce != 0 */
	       ambounce += (ambounce == 0);
               if (!check(3, "s")) {
                  /* File -> assume bwidth = 1 or precomputed pmap */
                  if (++t >= NUM_PMAP_TYPES)
                     error(USER, "too many photon maps");
                     
                  pmapParams [t].fileName = savqstr(av [1]);
                  pmapParams [t].minGather = pmapParams [t].maxGather =
                     defaultGather;
               }
               else return -1;
                  
               if (!check(3, "si")) {
                  /* File + bandwidth */
                  pmapParams [t].minGather = pmapParams [t].maxGather = 
                     atoi(av [2]);

                  if (!pmapParams [t].minGather)
                     return -1;
               }
               else {
                  sprintf(errmsg, "no photon lookup bandwidth specified, "
                          "using default %d", defaultGather);
                  error(WARNING, errmsg);
                  return 1;
               }
               
               if (!check(3, "sii")) {
                  /* File + min bwidth + max bwidth -> bias compensation */
                  pmapParams [t].maxGather = atoi(av [3]);

                  if (pmapParams [t].minGather >= pmapParams [t].maxGather)
                     return -1;

                  return 3;
               }
               else return 2;
               
            case 'm': /* Fixed max photon search radius */
               if (check(3, "f") || (maxDistFix = atof(av [1])) <= 0)
                  error(USER, "invalid max photon search radius coefficient");

               return 1;
         }      
   }
#undef check
   
   /* Unknown option */
   return -1;
}



void printPmapDefaults ()
/* Print defaults for photon map render options */
{
   printf("-am %.1f\t\t\t\t# max photon search radius\n", maxDistFix);
}
