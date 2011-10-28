set(VERSION_FILE "${radiance_BINARY_DIR}/src/hd/Version.c")
add_custom_command(
  OUTPUT "${VERSION_FILE}"
  COMMAND "${CMAKE_COMMAND}"
  -DVERSION_OUT_FILE="${VERSION_FILE}"
  -DVERSION_IN_FILE="${radiance_SOURCE_DIR}/src/rt/VERSION"
  -P "${radiance_SOURCE_DIR}/src/rt/create_version.cmake"
)

add_executable(genrhgrid genrhgrid.c holo.c)
target_link_libraries(genrhgrid rtrad)

add_executable(rhcopy rhcopy.c clumpbeams.c holo.c holofile.c)
target_link_libraries(rhcopy rtrad)

add_executable(rhinfo rhinfo.c holo.c holofile.c)
target_link_libraries(rhinfo rtrad)
if(NOT WIN32)
  set(VERSION_FILE "${radiance_BINARY_DIR}/src/hd/Version.c")
  create_version_file("${VERSION_FILE}")

  add_executable(rholo rholo.c rholo2.c rholo2l.c rholo3.c rholo4.c holo.c
    holofile.c viewbeams.c ${VERSION_FILE})
  target_link_libraries(rholo rtrad)
  add_executable(rhoptimize rhoptimize.c clumpbeams.c holo.c holofile.c)
  target_link_libraries(rhoptimize rtrad)
  add_executable(rhpict
    rhpict.c rhpict2.c holo.c holofile.c viewbeams.c ${VERSION_FILE})
  target_link_libraries(rhpict rtrad)
  install(TARGETS rholo rhoptimize rhpict RUNTIME DESTINATION bin) 
endif()

install(TARGETS genrhgrid rhcopy RUNTIME DESTINATION bin)
