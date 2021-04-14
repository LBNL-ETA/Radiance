#ifndef lint
static const char RCSid[] = "$Id: mkpmap.c,v 2.11 2021/04/14 11:26:25 rschregle Exp $";
#endif


/* 
   ======================================================================
   Photon map generator
   
   Roland Schregle (roland.schregle@{hslu.ch, gmail.com})
   (c) Fraunhofer Institute for Solar Energy Systems,
       supported by the German Research Foundation 
       (DFG LU-204/10-2, "Fassadenintegrierte Regelsysteme" (FARESYS)) 
   (c) Lucerne University of Applied Sciences and Arts,
       supported by the Swiss National Science Foundation 
       (SNSF #147053, "Daylight Redirecting Components")
   (c) Tokyo University of Science,
       supported by the JSPS Grants-in-Aid for Scientific Research 
       (KAKENHI JP19KK0115, "Three-Dimensional Light Flow")
   ======================================================================
   
   $Id: mkpmap.c,v 2.11 2021/04/14 11:26:25 rschregle Exp $    
*/


#include "pmap.h"
#include "pmapmat.h"
#include "pmapsrc.h"
#include "pmapcontrib.h"
#include "pmaprand.h"
#include "paths.h"
#include "ambient.h"
#include "resolu.h"
#include "source.h"
#include <ctype.h>
#include <string.h>
#include <sys/stat.h>


/* Enable options for Ze Ekspertz only! */
#define PMAP_EKSPERTZ


extern char VersionID [];


char*    progname;                  /* argv[0] */
int      dimlist [MAXDIM];          /* sampling dimensions */
int      ndims = 0;                 /* number of sampling dimenshunns */
char*    octname = NULL;            /* octree name */
CUBE     thescene;                  /* scene top-level octree */
OBJECT   nsceneobjs;                /* number of objects in scene */
double   srcsizerat = 0.01;         /* source partition size ratio */
int      backvis = 1;               /* back face visibility */
int      clobber = 0;               /* overwrite output */
COLOR    cextinction = BLKCOLOR;    /* global extinction coefficient */
COLOR    salbedo = BLKCOLOR;        /* global scattering albedo */
double   seccg = 0;                 /* global scattering eccentricity */
char     *amblist [AMBLLEN + 1];    /* ambient include/exclude list */
int      ambincl = -1;              /* include == 1, exclude == 0 */
char     *diagFile = NULL;          /* diagnostics output file */
int      rand_samp = 1;             /* uncorrelated random sampling */
unsigned nproc = 1;                 /* number of parallel processes */
#ifdef EVALDRC_HACK
   char  *angsrcfile = NULL;        /* angular source file for EvalDRC */
#endif


/* Dummies for linkage */

COLOR ambval = BLKCOLOR;
double   shadthresh = .05, ambacc = 0.2, shadcert = .5, minweight = 5e-3, 
         ssampdist = 0, dstrsrc = 0.0, specthresh = 0.15, specjitter = 1.0,
         avgrefl = 0.5;
int      ambvwt = 0, ambssamp = 0, ambres = 32, ambounce = 0, 
         directrelay = 1, directvis = 1, samplendx, do_irrad = 0, 
         ambdiv = 128, vspretest = 512, maxdepth = 6, contrib = 0;
char     *shm_boundary = NULL, *ambfile = NULL, *RCCONTEXT = NULL;
void     (*trace)() = NULL, (*addobjnotify [])() = {ambnotify, NULL};

  
void printdefaults()
/* print default values to stdout */
{
#ifdef EVALDRC_HACK
   /* EvalDRC support */
   puts("-A\t\t\t\t# angular source file");
#endif    
   puts("-ae  mod\t\t\t\t# exclude modifier");
   puts("-aE  file\t\t\t\t# exclude modifiers from file");
   puts("-ai  mod\t\t\t\t# include modifier");
   puts("-aI  file\t\t\t\t# include modifiers from file");
#ifdef PMAP_EKSPERTZ
   puts("-api xmin ymin zmin xmax ymax zmax\t# rectangular region of interest");
   puts("-apI xpos ypos zpos radius\t\t# spherical region of interest");
#endif
   puts("-apg file nPhotons\t\t\t# global photon map");
   puts("-apc file nPhotons\t\t\t# caustic photon map");
   puts("-apd file nPhotons\t\t\t# direct photon map");
   puts("-app file nPhotons bwidth\t\t# precomputed global photon map");
   puts("-apv file nPhotons\t\t\t# volume photon map");
   puts("-apC file nPhotons\t\t\t# contribution photon map");   
   printf("-apD %f\t\t\t\t# predistribution factor\n", preDistrib);
   printf("-apM %d\t\t\t\t\t# max predistrib passes\n", maxPreDistrib);
#if 1
   /* Kept for backwards compat, will be gradually phased out by -ld, -lr */
   printf("-apm %ld\t\t\t\t# limit photon bounces\n", photonMaxBounce);
#endif
   puts("-apo+ mod\t\t\t\t# photon port modifier");
   puts("-apO+ file\t\t\t\t# photon ports from file");
   printf("-apP %f\t\t\t\t# precomputation factor\n", finalGather);
   printf("-apr %d\t\t\t\t\t# random seed\n", randSeed);
   puts("-aps mod\t\t\t\t# antimatter sensor modifier");
   puts("-apS file\t\t\t\t# antimatter sensors from file");

   printf(backvis 
      ? "-bv+\t\t\t\t\t# back face visibility on\n"
      : "-bv-\t\t\t\t\t# back face visibility off\n"
   );
   printf("-dp  %.1f\t\t\t\t# PDF samples / sr\n", pdfSamples);
   printf("-ds  %f\t\t\t\t# source partition size ratio\n", srcsizerat);
   printf("-e   %s\t\t\t\t# diagnostics output file\n", diagFile);
   printf(clobber 
      ?  "-fo+\t\t\t\t\t# force overwrite\n" 
      :  "-fo-\t\t\t\t\t# do not overwrite\n"
   );
#ifdef PMAP_EKSPERTZ
   /* (Formerly) NU STUFF for Ze Exspertz! */      
   printf("-ld %.1f\t\t\t\t\t# limit photon distance\n", photonMaxDist);
   printf("-lr %ld\t\t\t\t# limit photon bounces\n", photonMaxBounce);   
#endif   
   printf(
      "-ma  %.2f %.2f %.2f\t\t\t# scattering albedo\n", 
      colval(salbedo,RED), colval(salbedo,GRN), colval(salbedo,BLU)
   );
   printf(
      "-me  %.2e %.2e %.2e\t\t# extinction coefficient\n", 
      colval(cextinction,RED), colval(cextinction,GRN), 
      colval(cextinction,BLU)
   );
   printf("-mg  %.2f\t\t\t\t# scattering eccentricity\n", seccg);
#if NIX   
   /* Multiprocessing on NIX only; so tuff luck, Windoze Weenies! */
   printf("-n   %d\t\t\t\t\t# number of parallel processes\n", nproc);
#endif   
   printf("-t   %-9d\t\t\t\t# time between reports\n", photonRepTime);
   printf(verbose 
      ?  "-v+\t\t\t\t\t# verbose console output\n"
      :  "-v-\t\t\t\t\t# terse console output\n"
   );
}


int main (int argc, char* argv [])
{
   #define check(ol, al) if ( \
      argv [i][ol] || badarg(argc - i - 1,argv + i + 1, al) \
   ) goto badopt
   
   /* Evaluate boolean option, setting var accordingly */                         
   #define check_bool(olen, var) switch (argv [i][olen]) { \
      case '\0': \
         var = !var; break; \
      case 'y': case 'Y': case 't': case 'T': case '+': case '1': \
         var = 1; break; \
      case 'n': case 'N': case 'f': case 'F': case '-': case '0': \
         var = 0; break; \
      default: \
         goto badopt; \
   }
   
   /* Evaluate trinary option, setting bits v1 and v2 in var accordingly */
   #define check_tri(olen, v1, v2, var) switch (argv [i][olen]) { \
      case '\0': case '+': \
         var = v1; break; \
      case '-': \
         var = v2; break;\
      case '0': \
         var = v1 | v2; break; \
      default: \
         goto badopt; \
   }

   int      loadflags = IO_CHECK | IO_SCENE | IO_TREE | IO_BOUNDS, 
            rval, i, j, n;
   char     **portLp = photonPortList, **sensLp = photonSensorList, 
            **amblp = NULL, sbuf [MAXSTR], portFlags [2] = "\0\0";
   struct   stat pmstat;

   /* Global program name */
   progname = fixargv0(argv [0]);
   /* Initialize object types */
   initotypes();
   
   /* Parse options */
   for (i = 1; i < argc; i++) {
      /* Eggs-pand arguments */
      while ((rval = expandarg(&argc, &argv, i)))
         if (rval < 0) {
            sprintf(errmsg, "cannot eggs-pand '%s'", argv [i]);
            error(SYSTEM, errmsg);
         }
         
      if (argv[i] == NULL) 
         break;
         
      if (!strcmp(argv [i], "-version")) {
         puts(VersionID);
         quit(0);
      }
      
      if (!strcmp(argv [i], "-defaults") || !strcmp(argv [i], "-help")) {
         printdefaults();
         quit(0);
      }
      
      /* Get octree */
      if (i == argc - 1) {
         octname = argv [i];
         break;
      }
            
      switch (argv [i][1]) {
         case 'a': /* Ambient */
            switch (argv [i][2]) {
               case 'i': /* Ambient include */
               case 'I':
                  check(3, "s");
                  if (ambincl != 1) {
                     ambincl = 1;
                     amblp = amblist;
                  }
                  if (isupper(argv [i][2])) {
                     /* Add modifiers from file */
                     rval = wordfile(
                        amblp, AMBLLEN - (amblp - amblist),
                        getpath(argv [++i], getrlibpath(), R_OK)
                     );
                     if (rval < 0) {
                        sprintf(
                           errmsg, "cannot open ambient include file \"%s\"",
                           argv [i]
                        );
                        error(SYSTEM, errmsg);
                     }
                     amblp += rval;
                  } 
                  else {
                     /* Add modifier from next arg */
                     *amblp++ = savqstr(argv [++i]);
                     *amblp = NULL;
                  }
                  break;

               case 'e': /* Ambient exclude */
               case 'E':
                  check(3, "s");
                  if (ambincl != 0) {
                     ambincl = 0;
                     amblp = amblist;
                  }
                  if (isupper(argv [i][2])) { 
                     /* Add modifiers from file */
                     rval = wordfile(
                        amblp, AMBLLEN - (amblp - amblist),
                        getpath(argv [++i], getrlibpath(), R_OK)
                     );
                     if (rval < 0) {
                        sprintf(
                           errmsg, "cannot open ambient exclude file \"%s\"", 
                           argv [i]
                        );
                        error(SYSTEM, errmsg);
                     }
                     amblp += rval;
                  } 
                  else {
                     /* Add modifier from next arg */ 
                     *amblp++ = savqstr(argv [++i]);
                     *amblp = NULL;
                  }
                  break;
            
               case 'p': /* Pmap-specific */
                  switch (argv [i][3]) {
                     case 'g': /* Global photon map */
                        check(4, "ss");
                        globalPmapParams.fileName = argv [++i];
                        globalPmapParams.distribTarget = 
                           parseMultiplier(argv [++i]);
                        if (!globalPmapParams.distribTarget) 
                           goto badopt;
                        globalPmapParams.minGather = 
                           globalPmapParams.maxGather = 0;
                        break;
            
                     case 'p': /* Precomputed global photon map */
                        check(4, "ssi");
                        preCompPmapParams.fileName = argv [++i];
                        preCompPmapParams.distribTarget = 
                           parseMultiplier(argv [++i]);
                        if (!preCompPmapParams.distribTarget) 
                           goto badopt;
                        preCompPmapParams.minGather = 
                           preCompPmapParams.maxGather = atoi(argv [++i]);
                        if (!preCompPmapParams.maxGather) 
                           goto badopt;
                        break;
            
                     case 'c': /* Caustic photon map */
                        check(4, "ss");
                        causticPmapParams.fileName = argv [++i];
                        causticPmapParams.distribTarget = 
                           parseMultiplier(argv [++i]);
                        if (!causticPmapParams.distribTarget) 
                           goto badopt;
                        break;
                  
                     case 'v': /* Volume photon map */
                        check(4, "ss");
                        volumePmapParams.fileName = argv [++i];
                        volumePmapParams.distribTarget = 
                           parseMultiplier(argv [++i]);
                        if (!volumePmapParams.distribTarget) 
                           goto badopt;
                        break;
                     
                     case 'd': /* Direct photon map */
                        check(4, "ss");
                        directPmapParams.fileName = argv [++i];
                        directPmapParams.distribTarget = 
                           parseMultiplier(argv [++i]);
                        if (!directPmapParams.distribTarget) 
                           goto badopt;
                        break;
                     
                     case 'C': /* Contribution photon map */
                        check(4, "ss");
                        contribPmapParams.fileName = argv [++i];
                        contribPmapParams.distribTarget =
                           parseMultiplier(argv [++i]);
                        if (!contribPmapParams.distribTarget)
                           goto badopt;
                        break;

                     case 'D': /* Predistribution factor */
                        check(4, "f");
                        preDistrib = atof(argv [++i]);
                        if (preDistrib <= 0)
                           error(USER, "predistrib factor must be > 0");
                        break;

                     case 'M': /* Max predistribution passes */
                        check(4, "i");
                        maxPreDistrib = atoi(argv [++i]);
                        if (maxPreDistrib <= 0)
                           error(USER, "max predistrib passes must be > 0");
                        break;

#if 1
                     /* Kept for backwards compat, to be phased out by -lr */
                     case 'm': /* Max photon bounces */
                        check(4, "i");
                        photonMaxBounce = atol(argv [++i]);
                        if (photonMaxBounce <= 0) 
                           error(USER, "max photon bounces must be > 0");
                        break;
#endif

#ifdef PMAP_EKSPERTZ
                     case 'i': /* Add rectangular region of interest */
                     case 'I': /* Add spherical region of interest */
                        check(4, isupper(argv [j=i][3]) ? "ffff" : "ffffff");
                        n = pmapNumROI;
                        
                        pmapROI = realloc(
                           pmapROI, ++pmapNumROI * sizeof(PhotonMapROI)
                        );
                        if (!pmapROI)
                           error(SYSTEM, "failed to allocate ROI");
                        
                        pmapROI [n].pos [0] = atof(argv [++i]);
                        pmapROI [n].pos [1] = atof(argv [++i]);
                        pmapROI [n].pos [2] = atof(argv [++i]);
                        pmapROI [n].siz [0] = atof(argv [++i]);

                        if (isupper(argv [j][3])) {
                           /* Spherical ROI; radius^2 */
                           pmapROI [n].siz [0] *= pmapROI [n].siz [0];
                           PMAP_ROI_SETSPHERE(pmapROI + n);
                           if (pmapROI [n].siz [0] <= FTINY)
                              error(
                                 USER, 
                                 "region of interest has invalid radius"
                              );
                        }
                        else {
                           /* Rectangular ROI */
                           pmapROI [n].siz [1] = atof(argv [++i]);
                           pmapROI [n].siz [2] = atof(argv [++i]);

                           for (j = 0; j < 3; j++) {
                              /* Pos at rectangle centre, siz symmetric */
                              pmapROI [n].pos [j] = 0.5 * (
                                 pmapROI [n].pos [j] + pmapROI [n].siz [j]
                              );
                              pmapROI [n].siz [j] = fabs(
                                 pmapROI [n].siz [j] - pmapROI [n].pos [j]
                              );
                              if (pmapROI [n].siz [j] <= FTINY)
                                 error(
                                    USER, 
                                    "region of interest has invalid size"
                                 );
                           }
                        }
                        break;
#endif

                     case 'P': /* Global photon precomp ratio */
                        check(4, "f");
                        finalGather = atof(argv [++i]);
                        if (finalGather <= 0 || finalGather > 1)
                           error(
                              USER, "global photon precomputation ratio "
                              "must be in range ]0, 1]"
                           );
                        break;
                     
                     case 'o': /* Photon port */ 
                     case 'O':
                        /* Check for bad arg and length, taking into account
                         * default forward orientation if none specified, in
                         * order to maintain previous behaviour */
                        check(argv [i][4] ? 5 : 4, "s");
                        /* Get port orientation flags */
                        check_tri(
                           4, PMAP_PORTFWD, PMAP_PORTBWD, portFlags [0]
                        );
                        
                        if (isupper(argv [i][3])) {
                           /* Add port modifiers from file */
                           rval = wordfile(
                              portLp, MAXSET - (portLp - photonPortList),
                              getpath(argv [++i], getrlibpath(), R_OK)
                           );
                           if (rval < 0) {
                               sprintf(
                                  errmsg, "cannot open photon port file %s", 
                                  argv [i]
                               );
                               error(SYSTEM, errmsg);
                           }
                           /* HACK: append port orientation flags to every 
                            * modifier; note this requires reallocation */
                           for (; rval--; portLp++) {
                              j = strlen(*portLp);
                              if (!(*portLp = realloc(*portLp, j + 2))) {
                                  sprintf(
                                     errmsg, 
                                     "cannot allocate photon port modifiers"
                                     " from file %s", argv [i]
                                  );
                                  error(SYSTEM, errmsg);
                              }
                              strcat(*portLp, portFlags);
                           }
                        } 
                        else {
                           /* Append port flags to port modifier, add to 
                            * port list and mark of end list with NULL */
                           strcpy(sbuf, argv [++i]);
                           strcat(sbuf, portFlags);
                           *portLp++ = savqstr(sbuf);
                           *portLp = NULL;
                        }
                        break;
                     
                     case 'r': /* Random seed */
                        check(4, "i");
                        randSeed = atoi(argv [++i]);
                        break;                   

                     case 's': /* Antimatter sensor */ 
                     case 'S':
                        check(4, "s");
                        if (isupper(argv[i][3])) {
                           /* Add sensor modifiers from file */
                           rval = wordfile(
                              sensLp, MAXSET - (sensLp - photonSensorList),
                              getpath(argv [++i], getrlibpath(), R_OK)
                           );
                           if (rval < 0) {
                               sprintf(
                                 errmsg, 
                                 "cannot open antimatter sensor file %s",
                                 argv [i]
                              );
                              error(SYSTEM, errmsg);
                           }
                           sensLp += rval;
                        } 
                        else {
                           /* Append modifier to sensor list, mark end with
                            * NULL */
                           *sensLp++ = savqstr(argv [++i]);
                           *sensLp = NULL;
                        }
                        break;

                     default: goto badopt;
                  }
                  break;
                  
               default: goto badopt;
            }
            break;
                
         case 'b': /* Back face visibility */
            if (argv [i][2] == 'v') {
               check_bool(3, backvis);
            }
            else goto badopt;
            break;
                   
         case 'd': /* Direct */
            switch (argv [i][2]) {
               case 'p': /* PDF samples */
                  check(3, "f");
                  pdfSamples = atof(argv [++i]);
                  break;
                  
               case 's': /* Source partition size ratio */
                  check(3, "f");
                  srcsizerat = atof(argv [++i]);
                  break;
                  
               default: goto badopt;
            }                   
            break;
                   
         case 'e': /* Diagnostics file */
            check(2, "s");
            diagFile = argv [++i];
            break;
                  
         case 'f': /* Force overwrite */
            if (argv [i][2] == 'o') {
               check_bool(3, clobber);
            }
            else goto badopt;
            break; 

#ifdef PMAP_EKSPERTZ
         case 'l': /* Limits */
            switch (argv [i][2]) {
               case 'd': /* Limit photon path distance */
                  check(3, "f");
                  photonMaxDist = atof(argv [++i]);
                  if (photonMaxDist <= 0)
                     error(USER, "max photon distance must be > 0");
                  break;
                 
               case 'r': /* Limit photon bounces */               
                  check(3, "i");               
                  photonMaxBounce = atol(argv [++i]);
                  if (photonMaxBounce <= 0) 
                     error(USER, "max photon bounces must be > 0");
                  break;
               
               default: goto badopt;
            }
            break;
#endif

         case 'm': /* Medium */
            switch (argv[i][2]) {
               case 'e': /* Eggs-tinction coefficient */
                  check(3, "fff");
                  setcolor(
                     cextinction, atof(argv [i + 1]),
                     atof(argv [i + 2]), atof(argv [i + 3])
                  );
                  i += 3;
                  break;
                                
               case 'a':	/* Albedo */
                  check(3, "fff");
                  setcolor(
                     salbedo, atof(argv [i + 1]), 
                     atof(argv [i + 2]), atof(argv [i + 3])
                  );
                  i += 3;
                  break;
                                
               case 'g':	/* Scattering eccentricity */
                  check(3, "f");
                  seccg = atof(argv [++i]);
                  break;
                                
               default: goto badopt;
            }                   
            break;
            
#if NIX
         case 'n': /* Num parallel processes (NIX only) */
            check(2, "i");
            nproc = atoi(argv [++i]);
            
            if (nproc > PMAP_MAXPROC) {
               nproc = PMAP_MAXPROC;
               sprintf(
                  errmsg, "too many parallel processes, clamping to %d\n", 
                  nproc
               );
               error(WARNING, errmsg);
            }
            break;
#endif

         case 't': /* Timer */
            check(2, "i");
            photonRepTime = atoi(argv [++i]);
            break;
            
         case 'v':   /* Verbosity */
            check_bool(2, verbose);
            break;
            
#ifdef EVALDRC_HACK
         case 'A':   /* Angular source file */
            check(2,"s");
            angsrcfile = argv[++i];
            break;
#endif

        default: goto badopt;
      }
   }
   
   /* Open diagnostics file */
   if (diagFile) {
      if (!freopen(diagFile, "a", stderr)) 
         quit(2);
      fprintf(stderr, "**************\n*** PID %5d: ", getpid());
      printargs(argc, argv, stderr);
      putc('\n', stderr);
      fflush(stderr);
   }
   
#ifdef NICE
   /* Lower priority */
   nice(NICE);
#endif

   if (octname == NULL) 
      error(USER, "missing octree argument");
      
   /* Allocate photon maps and set parameters */
   for (i = 0; i < NUM_PMAP_TYPES; i++) {
      setPmapParam(photonMaps + i, pmapParams + i);
      
      /* Don't overwrite existing photon map unless clobbering enabled */
      if (
         photonMaps [i] && !stat(photonMaps [i] -> fileName, &pmstat) &&
         !clobber
      ) {
         sprintf(
            errmsg, "photon map file %s exists, not overwritten",
           photonMaps [i] -> fileName
         );
         error(USER, errmsg);
      }
   }
       
   for (i = 0; i < NUM_PMAP_TYPES && !photonMaps [i]; i++);
   if (i >= NUM_PMAP_TYPES)
      error(USER, "no photon maps specified");
   
   readoct(octname, loadflags, &thescene, NULL);
#ifdef EVALDRC_HACK   
   if (angsrcfile)
      readobj(angsrcfile);    /* load angular sources */
#endif         
   nsceneobjs = nobjects;
   
   /* Get sources */
   marksources();

   /* Do forward pass and build photon maps */
   if (contribPmap)
      /* Just build contrib pmap, ignore others */
      distribPhotonContrib(contribPmap, nproc);
   else
      distribPhotons(photonMaps, nproc);
   
   /* Save photon maps; no idea why GCC needs an explicit cast here... */
   savePmaps((const PhotonMap**)photonMaps, argc, argv);
   cleanUpPmaps(photonMaps);
   
   quit(0);

badopt:
   sprintf(errmsg, "command line error at '%s'", argv[i]);
   error(USER, errmsg);

   #undef check
   #undef check_bool
   return 0;
}

