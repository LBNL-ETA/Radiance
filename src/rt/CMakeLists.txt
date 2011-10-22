set(VERSION_FILE "${radiance_BINARY_DIR}/src/rt/Version.c")
add_custom_command(
  OUTPUT "${VERSION_FILE}"
  COMMAND "${CMAKE_COMMAND}"
  -DVERSION_OUT_FILE="${VERSION_FILE}"
  -DVERSION_IN_FILE="${radiance_SOURCE_DIR}/src/util/VERSION"
  -P "${radiance_SOURCE_DIR}/src/common/create_version.cmake"
)

add_library(radiance
  ambcomp.c ambient.c ambio.c aniso.c data.c dielectric.c fprism.c freeobjmem.c 
  func.c glass.c initotypes.c m_alias.c m_brdf.c m_bsdf.c m_clip.c m_direct.c
  m_mirror.c m_mist.c mx_func.c mx_data.c noise3.c normal.c o_face.c o_cone.c
  o_instance.c o_mesh.c p_data.c p_func.c preload.c raytrace.c renderopts.c
  source.c sphere.c srcobstr.c srcsupp.c srcsamp.c t_data.c t_func.c text.c
  "${VERSION_FILE}" virtuals.c)

if(WIN32)
  set(rayp_SOURCES raypwin.c)
else()
  set(rayp_SOURCES raypcalls.c)
endif()

add_library(raycalls raycalls.c ${rayp_SOURCES} rayfifo.c)

add_executable(rtrace  rtmain.c rtrace.c duphead.c persist.c)
target_link_libraries(rtrace raycalls radiance rtrad)

add_executable(rpict rpmain.c rpict.c srcdraw.c duphead.c persist.c)
target_link_libraries(rpict radiance rtrad)
if(X11_FOUND)
  set(X11_SOURCES x11.c x11twind.c)
  add_definitions("-DHAS_X11")
endif()
if(QT4_FOUND)
  add_definitions("-DHAS_QT")
endif()
add_executable(rvu
  rvmain.c rview.c rv2.c rv3.c devtable.c devcomm.c editline.c ${X11_SOURCES}
  colortab.c)
target_link_libraries(rvu raycalls radiance rtrad ${X11_LIBRARIES})
# only link to rvu-qt-lib if qt is found
if(QT4_FOUND)
  target_link_libraries(rvu qtrvu)
endif()

add_executable(lookamb lookamb.c ambio.c)
target_link_libraries(lookamb rtrad)

add_test(rtrace_Version ${EXECUTABLE_OUTPUT_PATH}/rtrace -version)
set_tests_properties(rtrace_Version PROPERTIES
  PASS_REGULAR_EXPRESSION "RADIANCE")

install(TARGETS radiance raycalls rtrace rpict rvu
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION bin
        ARCHIVE DESTINATION bin)
install(FILES rayinit.cal DESTINATION bin)
