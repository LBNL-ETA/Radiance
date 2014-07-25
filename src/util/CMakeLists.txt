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

#the next few libraries all need to include header files from the rt directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../rt)

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

add_executable(eplus_adduvf eplus_adduvf.c eplus_idf.c)
target_link_libraries(eplus_adduvf rtrad)

add_executable(rfluxmtx rfluxmtx.c)
target_link_libraries(rfluxmtx rtrad)

if(X11_FOUND)
  add_executable(xglaresrc
    xglaresrc.c ${CMAKE_CURRENT_SOURCE_DIR}/../common/x11findwind.c)
  target_link_libraries(xglaresrc rtrad ${X11_LIBRARIES})
  install(TARGETS xglaresrc RUNTIME DESTINATION bin)
endif()

#glrad needs to link against OpenGL libraries
if(APPLE)
  find_library(OPENGL_gl_LIBRARY
    NAMES GL MesaGL
    PATHS /usr/X11R6/lib
  )
  find_library(OPENGL_glu_LIBRARY
    NAMES GLU MesaGLU
    PATHS /usr/X11R6/lib
  )
  set(OPENGL_LIBRARIES ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY})
else()
  find_package(OpenGL REQUIRED)
  include_directories(${OPENGL_INCLUDE_DIR})
endif()

#check if	we have	SGIStereo.h
find_file(HAVE_SGISTEREO_H "SGIStereo.h" ${OPENGL_INCLUDE_DIR})
if(NOT HAVE_SGISTEREO_H)
	add_definitions(-DNOSTEREO)
endif()
if(UNIX AND X11_FOUND)
	add_executable(glrad glrad.c)
	target_link_libraries(glrad rgl rtrad ${OPENGL_LIBRARIES} ${X11_LIBRARIES})
	install(TARGETS rpiece RUNTIME DESTINATION bin)
  if(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
    target_link_libraries(glrad rt)
  endif()
endif()
install(TARGETS dctimestep findglare getinfo glarendx rad ranimate
    vwrays vwright ranimove rsensor rttree_reduce rcollate
		eplus_adduvf rfluxmtx rmtxop 
    RUNTIME DESTINATION bin)

install(FILES rambpos.cal ambpos.cal tregsamp.dat reinhartb.cal 
    klems_full.cal klems_half.cal klems_quarter.cal disk2square.cal
    DESTINATION lib)

if(WIN32)
  install(PROGRAMS genBSDF.pl DESTINATION bin)
  install(PROGRAMS genklemsamp.pl DESTINATION bin)
  install(PROGRAMS genskyvec.pl DESTINATION bin)
  install(PROGRAMS ltview.pl DESTINATION bin)
  install(PROGRAMS ltpict.pl DESTINATION bin)
  install(PROGRAMS objview.pl DESTINATION bin)
  install(PROGRAMS objpict.pl DESTINATION bin)
else()
  install(PROGRAMS genBSDF.pl DESTINATION bin RENAME genBSDF)
  install(PROGRAMS genklemsamp.pl DESTINATION bin RENAME genklemsamp)
  install(PROGRAMS genskyvec.pl DESTINATION bin RENAME genskyvec)
  install(PROGRAMS genambpos.csh DESTINATION bin RENAME genambpos)
  install(PROGRAMS ltview.pl DESTINATION bin)
  install(PROGRAMS ltpict.pl DESTINATION bin)
  install(PROGRAMS objview.pl DESTINATION bin RENAME objview)
  install(PROGRAMS objpict.pl DESTINATION bin RENAME objpict)
  install(PROGRAMS glare.csh DESTINATION bin RENAME glare)
  install(PROGRAMS dayfact.csh DESTINATION bin RENAME dayfact)
  install(PROGRAMS rlux.csh DESTINATION bin RENAME rlux)
  install(PROGRAMS raddepend.csh DESTINATION bin RENAME raddepend)
  install(PROGRAMS trad.wsh DESTINATION bin RENAME trad)
  install(PROGRAMS objline.csh DESTINATION bin RENAME objline)
  install(PROGRAMS compamb.csh DESTINATION bin RENAME compamb)
  install(PROGRAMS vinfo.csh DESTINATION bin RENAME vinfo)
  install(PROGRAMS fieldcomb.csh DESTINATION bin RENAME fieldcomb)
endif()
