set(VERSION_FILE "${radiance_BINARY_DIR}/src/util/Version.c")
add_custom_command(
  OUTPUT "${VERSION_FILE}"
  DEPENDS "${radiance_SOURCE_DIR}/src/rt/VERSION"
  COMMAND "${CMAKE_COMMAND}"
  -DVERSION_OUT_FILE="${VERSION_FILE}"
  -DVERSION_IN_FILE="${radiance_SOURCE_DIR}/src/rt/VERSION"
  -P "${radiance_SOURCE_DIR}/src/rt/create_version.cmake"
)

add_executable(dctimestep dctimestep.c)
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

add_executable(rtcontrib
  rtcontrib.c "${VERSION_FILE}")
target_link_libraries(rtcontrib rtrad)

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

if(X11_FOUND)
  add_executable(xglaresrc
    xglaresrc.c ${CMAKE_CURRENT_SOURCE_DIR}/../common/x11findwind.c)
  target_link_libraries(xglaresrc rtrad ${X11_LIBRARIES})
  install(TARGETS xglaresrc RUNTIME DESTINATION bin)
endif()

find_library(LIB_RGL "rgl" DOC "Path to librgl")
if(LIB_RGL)
  #glrad needs to be able to include OpenGL headers
  include_directories(${OPENGL_INCLUDE_DIR})

  #check if we have SGIStereo.h
  find_file(HAVE_SGISTEREO_H "SGIStereo.h" ${OPENGL_INCLUDE_DIR})
  if(NOT HAVE_SGISTEREO_H)
    add_definitions(-DNOSTEREO)
  endif()

  add_executable(glrad glrad.c)
  target_link_libraries(glrad
    rtrad ${LIB_RGL} ${OPENGL_LIBRARIES} ${X11_LIBRARIES})
  install(TARGETS glrad RUNTIME DESTINATION bin)
endif()

install(TARGETS dctimestep findglare getinfo glarendx rad ranimate rpiece
                rtcontrib vwrays vwright ranimove rsensor
        RUNTIME DESTINATION bin)
install(FILES objview.csh objpict.csh glare.csh dayfact.csh rlux.csh
              raddepend.csh trad.wsh objline.csh compamb.csh vinfo.csh
              genambpos.csh fieldcomb.csh genklemsamp.pl genskyvec.pl
              genBSDF.pl tregsamp.dat
        DESTINATION bin)
