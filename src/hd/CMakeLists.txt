add_executable(genrhgrid genrhgrid.c holo.c)
target_link_libraries(genrhgrid rtrad)

add_executable(rhcopy rhcopy.c clumpbeams.c holo.c holofile.c)
target_link_libraries(rhcopy rtrad)

add_executable(rhinfo rhinfo.c holo.c holofile.c)
target_link_libraries(rhinfo rtrad)
if(UNIX)
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


	add_executable(x11.hdi rhdisp.c rhdisp2.c viewbeams.c rhd_x11.c rhd_qtree.c rhd_qtree2r.c rhd_ctab.c holo.c)
	target_link_libraries (x11.hdi rtrad ${X11_LIBRARIES})
  if(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
    target_link_libraries(x11.hdi rt)
  endif()
	install(TARGETS x11.hdi RUNTIME DESTINATION bin)
endif()


install(TARGETS genrhgrid rhcopy rhinfo RUNTIME DESTINATION bin)
