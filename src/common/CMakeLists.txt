set(rtrad_SOURCES
  addobjnotify.c
  badarg.c
  biggerlib.c
  bmalloc.c
  bmpfile.c
  bsdf.c
  bsdf_m.c
  bsdf_t.c
  byteswap.c
  caldefn.c
  calexpr.c
  calfunc.c
  calprnt.c
  ccolor.c
  ccyrgb.c
  chanvalue.c
  clip.c
  color.c
  colrops.c
  cone.c
  cvtcmd.c
  depthcodec.c
  dircode.c
  disk2square.c
  ealloc.c
  eputs.c
  erf.c
  error.c
  expandarg.c
  ezxml.c
  face.c
  falsecolor.c
  fdate.c
  fgetline.c
  fgetval.c
  fgetword.c
  fixargv0.c
  fltdepth.c
  font.c
  fputword.c
  free_os.c
  fropen.c
  fvect.c
  gethomedir.c
  getlibpath.c
  getpath.c
  header.c
  hilbert.c
  idmap.c
  image.c
  instance.c
  interp2d.c
  invmat4.c
  lamps.c
  linregr.c
  loadbsdf.c
  loadvars.c
  lookup.c
  mat4.c
  mesh.c
  modobject.c
  multisamp.c
  myhostname.c
  normcodec.c
  objset.c
  octree.c
  otypes.c
  paths.c
  plocate.c
  portio.c
  process.c
  quit.c
  readfargs.c
  readmesh.c
  readobj.c
  readoct.c
  resolu.c
  rexpr.c
  savestr.c
  savqstr.c
  sceneio.c
  spec_rgb.c
  tcos.c
  timegm.c
  tmap16bit.c
  tmapcolrs.c
  tmapluv.c
  tmaptiff.c
  tmesh.c
  tonemap.c
  triangulate.c
  urand.c
  urind.c
  wordfile.c
  words.c
  wputs.c
  xf.c
  zeroes.c
)

if(UNIX)
  list(APPEND rtrad_SOURCES unix_process.c)
  if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
    list(APPEND rtrad_SOURCES strlcpy.c)
  endif()
else()
  list(APPEND rtrad_SOURCES win_process.c win_popen.c win_usleep.c strlcpy.c)
endif()

add_library(rtrad ${rtrad_SOURCES})
target_link_libraries(rtrad ${LIB_M})
if(WIN32)
  target_link_libraries(rtrad ws2_32)
endif()

add_library(rgl
  rgldomat.c
  rglfile.c
  rglinst.c
  rglmat.c
  rglsrc.c
  rglsurf.c
)

add_library(mgf
  mgf_context.c
  mgf_object.c
  mgf_parser.c
  mgf_xf.c
)

add_library(wfobj
  objutil.c
  readwfobj.c
  convertobj.c
  writewfobj.c
  objtriangulate.c
)

add_library(libreetz
  g3affine.c
  g3affine.h
  g3flist.c
  g3flist.h
  g3list.h
  g3nlist.h
  g3sphere.c
  g3sphere.h
  g3vector.c
  g3vector.h
  gbasic.c
  gbasic.h
  muc_randvar.c
  muc_randvar.h
)

add_executable(testBSDF
  bsdf.c
  rtio.h
  testBSDF.c
)
target_link_libraries(testBSDF radiance rtrad)

install(TARGETS rtrad mgf
  LIBRARY DESTINATION "bin"
  ARCHIVE DESTINATION "bin"
)

file(COPY tmesh.cal
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)
