#ifndef lint
static const char RCSid[] = "$Id: pmapio.c,v 2.13 2021/02/18 17:08:50 rschregle Exp $";
#endif

/* 
   =======================================================================
   Photon map portable file I/O

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme FARESYS")
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components",
        SNSF #179067, "Light Fields for Spatio-Temporal Glare Assessment")
   =======================================================================
   
   $Id: pmapio.c,v 2.13 2021/02/18 17:08:50 rschregle Exp $    
*/



#include "pmapio.h"
#include "pmapdiag.h"
#include "resolu.h"



void savePhotonMap (const PhotonMap *pmap, const char *fname,
                    int argc, char **argv)
{
   unsigned long  i, j;
   FILE           *file;

   if (!pmap || !pmap -> numPhotons || !validPmapType(pmap -> type)) {
      error(INTERNAL, "attempt to save empty or invalid photon map");
      return;
   }
      
   if (verbose) {
      if (pmap -> numPrimary)
         sprintf(errmsg, "Saving %s (%ld photons, %d primaries)\n", 
                 fname, pmap -> numPhotons, pmap -> numPrimary);
      else sprintf(errmsg, "Saving %s (%ld photons)\n", fname, 
                   pmap -> numPhotons);
                   
      eputs(errmsg);
      fflush(stderr);
   }
   
   if (!(file = fopen(fname, "wb"))) {
      sprintf(errmsg, "can't open photon map file %s", fname);
      error(SYSTEM, errmsg);
   }
      
   /* Write header */
   newheader("RADIANCE", file);
   
   /* Write command line */
   printargs(argc, argv, file);
   
   /* Include statistics in info text */
   fprintf(
      file, 
      "NumPhotons\t= %ld\n"
      "AvgFlux\t\t= [%.2e, %.2e, %.2e]\n"
      "Bbox\t\t= [%.3f, %.3f, %.3f] [%.3f, %.3f, %.3f]\n"
      "CoG\t\t= [%.3f, %.3f, %.3f]\n"
      "MaxDist^2\t= %.3f\n",
      pmap -> numPhotons, 
      pmap -> photonFlux [0], pmap -> photonFlux [1], pmap -> photonFlux [2],
      pmap -> minPos [0], pmap -> minPos [1], pmap -> minPos [2],
      pmap -> maxPos [0], pmap -> maxPos [1], pmap -> maxPos [2],
      pmap -> CoG [0], pmap -> CoG [1], pmap -> CoG [2], 
      pmap -> CoGdist
   );
           
   if (pmap -> primaries)
      fprintf(file, "%d primary rays\n", pmap -> numPrimary);
   
   /* Write format, including human-readable file version */
   fputformat((char*)pmapFormat [pmap -> type], file);
   fprintf(file, "VERSION=%s\n", PMAP_FILEVER);
   
   /* Empty line = end of header */
   putc('\n', file);
   
   /* Write machine-readable file format version */
   putstr(PMAP_FILEVER, file); 
   
   /* Write number of photons as fixed size, which possibly results in
    * padding of MSB with 0 on some platforms.  Unlike sizeof() however,
    * this ensures portability since this value may span 32 or 64 bits
    * depending on platform.  */
   putint(pmap -> numPhotons, PMAP_LONGSIZE, file);
   
   /* Write average photon flux */
   for (j = 0; j < 3; j++) 
      putflt(pmap -> photonFlux [j], file);

   /* Write max and min photon positions */
   for (j = 0; j < 3; j++) {
      putflt(pmap -> minPos [j], file);
      putflt(pmap -> maxPos [j], file);
   }

   /* Write centre of gravity */
   for (j = 0; j < 3; j++)
      putflt(pmap -> CoG [j], file);
      
   /* Write avg distance to centre of gravity */
   putflt(pmap -> CoGdist, file);
   
   /* Write out primary photon rays (or just zero count if none) */
   if (pmap -> primaries) {
      putint(pmap -> numPrimary, sizeof(pmap -> numPrimary), file);
      
      for (i = 0; i < pmap -> numPrimary; i++) {
         PhotonPrimary *prim = pmap -> primaries + i;      

         putint(prim -> srcIdx, sizeof(prim -> srcIdx), file);
#ifdef PMAP_PRIMARYDIR         
         putint(prim -> dir, sizeof(prim -> dir), file);
#endif         
#ifdef PMAP_PRIMARYPOS         
         for (j = 0; j < 3; j++)
            putflt(prim -> pos [j], file);
#endif
         if (ferror(file))
            error(SYSTEM, "error writing primary photon rays");
      }
   }
   else putint(0, sizeof(pmap -> numPrimary), file);
   
   /* Save photon storage */
#ifdef PMAP_OOC
   if (OOC_SavePhotons(pmap, file)) {
#else
   if (kdT_SavePhotons(pmap, file)) {
#endif
      sprintf(errmsg, "error writing photon map file %s", fname);
      error(SYSTEM, errmsg);
   }

   fclose(file);
}



PhotonMapType loadPhotonMap (PhotonMap *pmap, const char *fname)
{
   PhotonMapType  ptype = PMAP_TYPE_NONE;
   char           format [MAXFMTLEN];
   unsigned long  i, j;
   FILE           *file;

   if (!pmap)
      return PMAP_TYPE_NONE;
   
   if ((file = fopen(fname, "rb")) == NULL) {
      sprintf(errmsg, "can't open photon map file %s", fname);
      error(SYSTEM, errmsg);
   }

   /* Get format string */
   strcpy(format, PMAP_FORMAT_GLOB);
   if (checkheader(file, format, NULL) != 1) {
      sprintf(errmsg, "photon map file %s has unknown format %s", 
              fname, format);
      error(USER, errmsg);
   }
   
   /* Identify photon map type from format string */
   for (ptype = 0; 
        ptype < NUM_PMAP_TYPES && strcmp(pmapFormat [ptype], format);
        ptype++);
      
   if (!validPmapType(ptype)) {
      sprintf(errmsg, "file %s contains an unknown photon map type", fname);
      error(USER, errmsg);
   }

   initPhotonMap(pmap, ptype);
   
   /* Get file format version and check for compatibility */
   if (strcmp(getstr(format, file), PMAP_FILEVER))
      error(USER, "incompatible photon map file format");
      
   /* Get number of photons as fixed size, which possibly results in
    * padding of MSB with 0 on some platforms.  Unlike sizeof() however,
    * this ensures portability since this value may span 32 or 64 bits
    * depending on platform.  */
   pmap -> numPhotons = getint(PMAP_LONGSIZE, file);
      
   /* Get average photon flux */
   for (j = 0; j < 3; j++) 
      pmap -> photonFlux [j] = getflt(file);
   
   /* Get max and min photon positions */
   for (j = 0; j < 3; j++) {
      pmap -> minPos [j] = getflt(file);
      pmap -> maxPos [j] = getflt(file);
   }
   
   /* Get centre of gravity */
   for (j = 0; j < 3; j++)
      pmap -> CoG [j] = getflt(file);
      
   /* Get avg distance to centre of gravity */
   pmap -> CoGdist = getflt(file);

   /* Get primary photon rays */
   pmap -> numPrimary = getint(sizeof(pmap -> numPrimary), file);
   
   if (pmap -> numPrimary) {
      pmap -> primaries = calloc(pmap -> numPrimary, sizeof(PhotonPrimary));      
      if (!pmap -> primaries)
         error(INTERNAL, "can't allocate primary photon rays");
         
      for (i = 0; i < pmap -> numPrimary; i++) {
         PhotonPrimary *prim = pmap -> primaries + i;
         
         prim -> srcIdx = getint(sizeof(prim -> srcIdx), file);
#ifdef PMAP_PRIMARYDIR         
         prim -> dir = getint(sizeof(prim -> dir), file);
#endif         
#ifdef PMAP_PRIMARYPOS            
         for (j = 0; j < 3; j++)
            prim -> pos [j] = getflt(file);
#endif            
         if (feof(file))
            error(SYSTEM, "error reading primary photon rays");
      }
   }      
   
   /* Load photon storage */
#ifdef PMAP_OOC
   if (OOC_LoadPhotons(pmap, file)) {
#else
   if (kdT_LoadPhotons(pmap, file)) {
#endif
      sprintf(errmsg, "error reading photon map file %s", fname);
      error(SYSTEM, errmsg);
   }

   fclose(file);
   return ptype;
}
