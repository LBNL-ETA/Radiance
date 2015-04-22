/* 
   ==================================================================
   Photon map diagnostic output and progress reports

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts
   ==================================================================
   
   $Id: pmapdiag.c,v 2.3 2015/04/22 20:28:16 rschregle Exp $
*/



#include "pmapdiag.h"
#include "pmapdata.h"
#include "standard.h"



time_t repStartTime, repLastTime = 0;   /* Time at start & last report */
unsigned long repProgress,              /* Report progress counter */
              repComplete;              /* Report completion count */



static char* biasCompStats (const PhotonMap *pmap, PhotonMapType type, 
                            char *stats)
/* Dump bias compensation statistics */
{
   unsigned avgBwidth;
   float avgErr;
   
   /* Check if photon map is valid and bias compensated */
   if (pmap && pmap -> maxGather > pmap -> minGather) {      
      avgBwidth = pmap -> numDensity 
                  ? pmap -> totalGathered / pmap -> numDensity : 0,
      avgErr = pmap -> numDensity 
               ? sqrt(pmap -> rmsError / pmap -> numDensity) : 0;
                                              
      sprintf(stats, "%d/%d/%d %s (%.1f/%.1f/%.1f%% error), ",
              pmap -> minGathered, pmap -> maxGathered, avgBwidth, 
              pmapName [type],
              100 * pmap -> minError, 100 * pmap -> maxError, 100 * avgErr);
                 
      return stats;
   }
   
   return NULL;
}



void pmapBiasCompReport (char *stats)
/* Append full bias compensation statistics to stats; interface to rpict's
 * report() */
{
   char tmp [512];
   unsigned t;
   
   stats [0] = '\0';
   
   for (t = 0; t < NUM_PMAP_TYPES; t++)
      if (biasCompStats(photonMaps [t], t, tmp))
         strcat(stats, tmp);
}



#ifndef NON_POSIX
   void pmapPreCompReport()
   /* Report global photon precomputation progress */
   {
      extern char *myhostname();
      float u, s;
      char tmp [512];
      
      #ifdef BSD
         struct rusage rubuf;
      #else
         struct tms tbuf;
         float period;
      #endif

      repLastTime = time(NULL);
      
      #ifdef BSD
         getrusage(RUSAGE_SELF, &rubuf);
         u = rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s = rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
         getrusage(RUSAGE_CHILDREN, &rubuf);
         u += rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s += rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
      #else
         times(&tbuf);

         #ifdef _SC_CLK_TCK
            period = 1.0 / sysconf(_SC_CLK_TCK);
         #else
            period = 1.0 / 60;
         #endif
      
         u = (tbuf.tms_utime + tbuf.tms_cutime) * period;
         s = (tbuf.tms_stime + tbuf.tms_cstime) * period;
      #endif

      sprintf(errmsg, "%lu precomputed, ", repProgress);

      /* Append bias compensation stats */      
      biasCompStats(preCompPmap, PMAP_TYPE_PRECOMP, tmp);
      strcat(errmsg, tmp);
      
      sprintf(tmp, "%4.2f%% after %.3fu %.3fs %.3fr hours on %s\n",
              100.0 * repProgress / repComplete, u / 3600, s / 3600, 
              (repLastTime - repStartTime) / 3600.0, myhostname());
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
      
      #ifdef SIGCONT
         signal(SIGCONT, pmapPreCompReport);
      #endif
   }



   void pmapDistribReport()
   /* Report photon distribution progress */
   {
      extern char *myhostname();
      float u, s;
      unsigned t;
      char tmp [512];
      
      #ifdef BSD
         struct rusage rubuf;
      #else
         struct tms tbuf;
         float period;
      #endif

      repLastTime = time(NULL);
      
      #ifdef BSD
         getrusage(RUSAGE_SELF, &rubuf);
         u = rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s = rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
         getrusage(RUSAGE_CHILDREN, &rubuf);
         u += rubuf.ru_utime.tv_sec + rubuf.ru_utime.tv_usec / 1e6;
         s += rubuf.ru_stime.tv_sec + rubuf.ru_stime.tv_usec / 1e6;
      #else
         times(&tbuf);
      
         #ifdef _SC_CLK_TCK
            period = 1.0 / sysconf(_SC_CLK_TCK);
         #else
            period = 1.0 / 60;
         #endif
      
         u = (tbuf.tms_utime + tbuf.tms_cutime) * period;
         s = (tbuf.tms_stime + tbuf.tms_cstime) * period;
      #endif

      sprintf(errmsg, "%lu emitted, ", repProgress);
      
      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (photonMaps [t]) {
            sprintf(tmp, "%lu %s, ", 
                    photonMaps [t] -> heapEnd, pmapName [t]);
            strcat(errmsg, tmp);
         }
      
      sprintf(tmp, "%4.2f%% after %.3fu %.3fs %.3fr hours on %s\n", 
              100.0 * repProgress / repComplete, u / 3600, s / 3600, 
              (repLastTime - repStartTime) / 3600.0, myhostname());
              
      strcat(errmsg, tmp);
      eputs(errmsg);
      fflush(stderr);
      
      #ifndef SIGCONT
         signal(SIGCONT, pmapDistribReport);
      #endif
   }
   
#else /* POSIX */

   void pmapPreCompReport()
   /* Report global photon precomputation progress */
   {
      char tmp [512];
      
      repLastTime = time(NULL);
      sprintf(errmsg, "%lu precomputed, ", repProgress);

      /* Append bias compensation stats */
      biasCompStats(preCompPmap, PMAP_TYPE_PRECOMP, tmp);
      strcat(errmsg, tmp);      

      sprintf(tmp, "%4.2f%% after %5.4f hours\n", 
              100.0 * repProgress / repComplete, 
              (repLastTime - repStartTime) / 3600.0);
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
   }



   void pmapDistribReport()
   /* Report photon distribution progress */
   {
      char tmp [512];
      unsigned t;
      
      repLastTime = time(NULL);
      sprintf(errmsg, "%lu emitted, ", repProgress);

      for (t = 0; t < NUM_PMAP_TYPES; t++)
         if (photonMaps [t]) {
            sprintf(tmp, "%lu %s, ", 
                    photonMaps [t] -> heapEnd, pmapName [t]);
            strcat(errmsg, tmp);
         }      
      
      sprintf(tmp, "%4.2f%% after %5.4f hours\n", 
              100.0 * repProgress / repComplete, 
              (repLastTime - repStartTime) / 3600.0);
      strcat(errmsg, tmp);
      
      eputs(errmsg);
      fflush(stderr);
   }
#endif



