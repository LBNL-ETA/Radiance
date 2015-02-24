/* 
   ===============================================================   
   Dump photon maps as RADIANCE scene description to stdout

   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       Lucerne University of Applied Sciences & Arts   
   ==================================================================
   
   $Id: pmapdump.c,v 2.1 2015/02/24 19:39:27 greg Exp $
*/



#include "pmapio.h"
#include "pmapparm.h"
#include "pmaptype.h"
#include "rtio.h"
#include "resolu.h"
#include "random.h"


/* Defaults */
/*    Sphere radius as fraction of avg. intersphere dist */
/*    Relative scale for sphere radius (fudge factor) */
/*    Number of spheres */
#define RADCOEFF 0.05
#define RADSCALE 1.0
#define NSPHERES 10000


/* RADIANCE material and object defs for each photon type */
typedef struct {
   char *mat, *obj;
} RadianceDef;

   
static char header [] = "$Revision: 2.1 $";


/* Colour code is as follows:    global         = blue
                                 precomp global = cyan
                                 caustic        = red
                                 volume         = green
                                 direct         = magenta 
                                 contrib        = yellow */   
const RadianceDef radDefs [] = {
   {  "void plastic mat.global\n0\n0\n5 0 0 1 0 0\n",
      "mat.global sphere obj.global\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void plastic mat.pglobal\n0\n0\n5 0 1 1 0 0\n",
      "mat.pglobal sphere obj.global\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void plastic mat.caustic\n0\n0\n5 1 0 0 0 0\n",
      "mat.caustic sphere obj.caustic\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void plastic mat.volume\n0\n0\n5 0 1 0 0 0\n",
      "mat.volume sphere obj.volume\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void plastic mat.direct\n0\n0\n5 1 0 1 0 0\n",
      "mat.direct sphere obj.direct\n0\n0\n4 %g %g %g %g\n"
   },
   {  "void plastic mat.contrib\n0\n0\n5 1 1 0 0 0\n",
      "mat.contrib sphere obj.contrib\n0\n0\n4 %g %g %g %g\n"
   }
};



int main (int argc, char** argv)
{
   char format [128];
   RREAL rad, radScale = RADSCALE, vol, dumpRatio;
   FVECT minPos, maxPos;
   unsigned arg, j, ptype;
   long numPhotons, numSpheres = NSPHERES;
   FILE *pmapFile;
   Photon p;
   
   if (argc < 2) {
      puts("Dump photon maps as RADIANCE scene description\n");
      printf("Usage: %s [-r radscale1] [-n nspheres1] pmap1 "
             "[-r radscale2] [-n nspheres2] pmap2 ...\n", argv [0]);
      return 1;
   }
   
   for (arg = 1; arg < argc; arg++) {
      /* Parse options */
      if (argv [arg][0] == '-') {
         switch (argv [arg][1]) {
            case 'r':
               if ((radScale = atof(argv [++arg])) <= 0)
                  error(USER, "invalid radius scale");
               break;
               
            case 'n':
               if ((numSpheres = parseMultiplier(argv [++arg])) <= 0)
                  error(USER, "invalid number of spheres");
               break;
               
            default:
               sprintf(errmsg, "unknown option %s", argv [arg]);
               error(USER, errmsg);
               return -1;
         }
         
         continue;
      }
      
      /* Dump photon map */
      if (!(pmapFile = fopen(argv [arg], "rb"))) {
         sprintf(errmsg, "can't open %s", argv [arg]);
         error(SYSTEM, errmsg);
      }
         
      /* Get format string */
      strcpy(format, PMAP_FORMAT_GLOB);
      if (checkheader(pmapFile, format, NULL) != 1) {
         sprintf(errmsg, "photon map file %s has unknown format %s", 
                 argv [arg], format);
         error(USER, errmsg);
      }
      
      /* Identify photon map type from format string */
      for (ptype = 0; 
           strcmp(pmapFormat [ptype], format) && ptype < NUM_PMAP_TYPES; 
           ptype++);
      
      if (!validPmapType(ptype)) {
         sprintf(errmsg, "file %s contains an unknown photon map type", 
                argv [arg]);
         error(USER, errmsg);
      }

      /* Get file format version and check for compatibility */
      if (getint(sizeof(PMAP_FILEVER), pmapFile) != PMAP_FILEVER)
         error(USER, "incompatible photon map file format");
         
      /* Dump command line as comment */
      fputs("# ", stdout);
      printargs(argc, argv, stdout);
      fputc('\n', stdout);
      
      /* Dump material def */   
      fputs(radDefs [ptype].mat, stdout);
      fputc('\n', stdout);
      
      /* Get number of photons (is this sizeof() hack portable?) */
      numPhotons = getint(sizeof(((PhotonMap*)NULL) -> heapSize), pmapFile);
      
      /* Skip avg photon flux */ 
      for (j = 0; j < 3; j++) 
         getflt(pmapFile);
      
      /* Get distribution extent (min & max photon positions) */
      for (j = 0; j < 3; j++) {
         minPos [j] = getflt(pmapFile);
         maxPos [j] = getflt(pmapFile);
      }
      
      /* Skip centre of gravity, and avg photon dist to it */
      for (j = 0; j < 4; j++)
         getflt(pmapFile);
      
      /* Sphere radius based on avg intersphere dist 
         (= sphere distrib density ^-1/3) */
      vol = (maxPos [0] - minPos [0]) * (maxPos [1] - minPos [1]) * 
            (maxPos [2] - minPos [2]);
      rad = radScale * RADCOEFF * cbrt(vol / numSpheres);
      
      /* Photon dump probability to satisfy target sphere count */
      dumpRatio = numSpheres < numPhotons ? (float)numSpheres / numPhotons
                                          : 1;
      
      while (numPhotons-- > 0) {
         /* Get photon position */            
         for (j = 0; j < 3; j++) 
            p.pos [j] = getflt(pmapFile);

         /* Dump photon probabilistically acc. to target sphere count */
         if (frandom() <= dumpRatio) {
            printf(radDefs [ptype].obj, p.pos [0], p.pos [1], p.pos [2], rad);
            fputc('\n', stdout);
         }
         
         /* Skip photon normal and flux */
         for (j = 0; j < 3; j++) 
            getint(sizeof(p.norm [j]), pmapFile);
            
         #ifdef PMAP_FLOAT_FLUX
            for (j = 0; j < 3; j++) 
               getflt(pmapFile);
         #else      
            for (j = 0; j < 4; j++) 
               getint(1, pmapFile);
         #endif

         /* Skip primary ray index */
         getint(sizeof(p.primary), pmapFile);

         /* Skip flags */
         getint(sizeof(p.flags), pmapFile);
         
         if (feof(pmapFile)) {
            sprintf(errmsg, "error reading %s", argv [arg]);
            error(USER, errmsg);
         }
      }
      
      fclose(pmapFile);
      
      /* Reset defaults for next dump */
      radScale = RADSCALE;
      numSpheres = NSPHERES;
   }
   
   return 0;
}
