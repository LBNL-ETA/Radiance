set(rtrad_SOURCES
  addobjnotify.c badarg.c biggerlib.c bmalloc.c bmpfile.c bsdf.c bsdf_m.c
  bsdf_t.c byteswap.c caldefn.c calexpr.c calfunc.c calprnt.c ccolor.c ccyrgb.c
  chanvalue.c clip.c color.c colrops.c cone.c dircode.c disk2square.c ealloc.c
  eputs.c erf.c error.c expandarg.c ezxml.c face.c falsecolor.c fdate.c
  fgetline.c fgetval.c fgetword.c fixargv0.c font.c fputword.c free_os.c
  fropen.c fvect.c gethomedir.c getlibpath.c getpath.c header.c hilbert.c
  image.c instance.c invmat4.c lamps.c linregr.c loadbsdf.c loadvars.c lookup.c
  mat4.c mesh.c modobject.c multisamp.c myhostname.c objset.c octree.c otypes.c
  paths.c plocate.c portio.c process.c quit.c readfargs.c readmesh.c readobj.c
  readoct.c resolu.c rexpr.c savqstr.c savestr.c sceneio.c spec_rgb.c tcos.c
  tmap16bit.c tmapcolrs.c tmapluv.c timegm.c tmaptiff.c tmesh.c tonemap.c
  urand.c urind.c words.c wordfile.c wputs.c xf.c zeroes.c)

if(UNIX)
  list(APPEND rtrad_SOURCES unix_process.c strcmp.c)
else()
  list(APPEND rtrad_SOURCES win_process.c win_popen.c)
endif()

add_library(rtrad ${rtrad_SOURCES})
target_link_libraries(rtrad ${LIB_M})
if(WIN32)
  target_link_libraries(rtrad ws2_32)
endif()

add_library(mgf mgf_parser.c mgf_object.c mgf_xf.c mgf_context.c)

install(TARGETS rtrad mgf
  LIBRARY DESTINATION bin
  ARCHIVE DESTINATION bin)
