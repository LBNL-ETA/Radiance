set(VERSION_FILE "${radiance_BINARY_DIR}/src/util/Version.c")
create_version_file("${VERSION_FILE}")

add_executable(dctimestep dctimestep.c cmbsdf.c cmatrix.c)
target_link_libraries(dctimestep rtrad)

add_executable(findglare findglare.c glareval.c glaresrc.c setscan.c)
target_link_libraries(findglare rtrad)

add_executable(getinfo getinfo.c)
target_link_libraries(getinfo rtrad)

add_executable(glarendx glarendx.c)
target_link_libraries(glarendx rtrad)

add_executable(rad rad.c)
target_link_libraries(rad rtrad)

if(WIN32)
  set(netproc_SOURCES win_netproc.c)
else()
  set(netproc_SOURCES netproc.c)
endif()

add_executable(ranimate ranimate.c ${netproc_SOURCES})
target_link_libraries(ranimate rtrad)

add_executable(rpiece rpiece.c "${VERSION_FILE}")
target_link_libraries(rpiece rtrad)

add_executable(vwrays vwrays.c)
target_link_libraries(vwrays rtrad)

add_executable(vwright vwright.c)
target_link_libraries(vwright rtrad)

add_executable(dcglare dcglare.c cmbsdf.c cmatrix.c cmglare.c)
target_link_libraries(dcglare rtrad)

#the next few libraries all need to include header files from the rt directory
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../rt")

add_executable(ranimove ranimove.c ranimove1.c ranimove2.c)
target_link_libraries(ranimove raycalls radiance rtrad)

add_executable(rsensor rsensor.c)
target_link_libraries(rsensor raycalls radiance rtrad)

add_executable(rttree_reduce rttree_reduce.c)
target_link_libraries(rttree_reduce rtrad)

add_executable(rcollate rcollate.c)
target_link_libraries(rcollate rtrad)

add_executable(rmtxop rmtxop.c rmatrix.c cmbsdf.c cmatrix.c)
target_link_libraries(rmtxop rtrad)

add_executable(wrapBSDF wrapBSDF.c)
target_link_libraries(wrapBSDF rtrad)

add_executable(eplus_adduvf eplus_adduvf.c eplus_idf.c)
target_link_libraries(eplus_adduvf rtrad)

add_executable(rfluxmtx rfluxmtx.c)
target_link_libraries(rfluxmtx rtrad)

add_executable(evalglare evalglare.c pictool.c)
target_link_libraries(evalglare rtrad libreetz)

add_executable(radcompare radcompare.c)
target_link_libraries(radcompare rtrad)

add_executable(rcode_depth rcode_depth.c)
target_link_libraries(rcode_depth rtrad)

add_executable(rcode_norm rcode_norm.c)
target_link_libraries(rcode_norm rtrad)

add_executable(rcode_ident rcode_ident.c)
target_link_libraries(rcode_ident rtrad)

if(X11_FOUND)
  add_executable(xglaresrc
    xglaresrc.c
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/x11findwind.c"
  )
  target_link_libraries(xglaresrc rtrad ${X11_LIBRARIES})
  install(TARGETS xglaresrc
    RUNTIME DESTINATION "bin"
  )
endif()

#glrad needs to link against OpenGL libraries
if(APPLE)
  find_library(OPENGL_gl_LIBRARY
    NAMES GL MesaGL
    PATHS ${X11_LIB_SEARCH_PATH}
  )
  find_library(OPENGL_glu_LIBRARY
    NAMES GLU MesaGLU
    PATHS ${X11_LIB_SEARCH_PATH}
  )
  set(OPENGL_LIBRARIES ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
else()
  find_package(OpenGL REQUIRED)
  include_directories(${OPENGL_INCLUDE_DIR})
endif()

#check if we have SGIStereo.h
find_file(HAVE_SGISTEREO_H "SGIStereo.h" ${OPENGL_INCLUDE_DIR})
if(NOT HAVE_SGISTEREO_H)
  add_definitions(-DNOSTEREO)
endif()
if(UNIX AND X11_FOUND)
  add_executable(glrad glrad.c)
  target_link_libraries(glrad rgl rtrad ${OPENGL_LIBRARIES} ${X11_LIBRARIES})
  install(TARGETS glrad
    RUNTIME DESTINATION "bin"
  )
  if(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
    target_link_libraries(glrad rt)
  endif()
endif()

if(UNIX)
  install(TARGETS rpiece
    RUNTIME DESTINATION "bin"
  )
endif()

install(TARGETS
  dcglare
  dctimestep
  eplus_adduvf
  evalglare
  findglare
  getinfo
  glarendx
  rad
  radcompare
  ranimate
  ranimove
  rcode_depth
  rcode_ident
  rcode_norm
  rcollate
  rfluxmtx
  rmtxop
  rsensor
  rttree_reduce
  vwrays
  vwright
  wrapBSDF
  RUNTIME DESTINATION "bin"
)

file(COPY
  ambpos.cal
  disk2square.cal
  klems_full.cal
  klems_half.cal
  klems_quarter.cal
  minimalBSDFt.xml
  rambpos.cal
  reinhartb.cal
  tregsamp.dat
  WINDOW6BSDFt.xml
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)

if(WIN32)
   add_custom_target(bsdfview ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/bsdfview.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/bsdfview.exe -d
  )
  add_custom_target(genBSDF ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/genBSDF.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/genBSDF.exe -d
  )
  add_custom_target(genklemsamp ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/genklemsamp.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/genklemsamp.exe -d
  )
  add_custom_target(genskyvec ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/genskyvec.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/genskyvec.exe -d
  )
  add_custom_target(ltpict ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/ltpict.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/ltpict.exe -d
  )
  add_custom_target(ltview ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/ltview.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/ltview.exe -d
  )
  add_custom_target(objpict ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/objpict.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/objpict.exe -d
  )
  add_custom_target(objview ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/objview.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/objview.exe -d
  )
  add_custom_target(genambpos ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/genambpos.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/genambpos.exe -d
  )
  add_custom_target(rcode2bmp ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/rcode2bmp.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/rcode2bmp.exe -d
  )
  add_custom_target(rtpict ALL
    COMMAND pp ${CMAKE_CURRENT_SOURCE_DIR}/rtpict.pl -o ${CMAKE_CURRENT_SOURCE_DIR}/rtpict.exe -d
  )
  install(PROGRAMS
    bsdfview.exe
    genBSDF.exe
    genklemsamp.exe
    genskyvec.exe
    ltpict.exe
    ltview.exe
    objpict.exe
    objview.exe
    genambpos.exe
    rcode2bmp.exe
    rtpict.exe
    DESTINATION "bin"
  )
else()
  install_without_extension(
    bsdfview.pl
    compamb.csh
    dayfact.csh
    fieldcomb.csh
    genBSDF.pl
    genambpos.pl
    genklemsamp.pl
    genskyvec.pl
    glare.csh
    ltpict.pl
    ltview.pl
    objline.csh
    objpict.pl
    objview.pl
    raddepend.csh
    rcode2bmp.pl
    rlux.csh
    rtpict.pl
    vinfo.csh
  )
  add_custom_target(trad ALL
    COMMAND csh -f ${CMAKE_CURRENT_SOURCE_DIR}/tradinstall.csh ../../build/bin ../../build/lib
  )
  install(FILES
    ${CMAKE_BINARY_DIR}/bin/trad
    DESTINATION "bin"
    PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ
  )
endif()
