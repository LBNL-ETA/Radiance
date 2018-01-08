add_executable(genrhgrid genrhgrid.c holo.c)
target_link_libraries(genrhgrid rtrad)

add_executable(rhcopy rhcopy.c clumpbeams.c holo.c holofile.c)
target_link_libraries(rhcopy rtrad)

add_executable(rhinfo rhinfo.c holo.c holofile.c)
target_link_libraries(rhinfo rtrad)

if(UNIX)
  set(VERSION_FILE "${radiance_BINARY_DIR}/src/hd/Version.c")
  create_version_file("${VERSION_FILE}")

  add_executable(rholo
    holo.c
    holofile.c
    rholo.c
    rholo2.c
    rholo2l.c
    rholo3.c
    rholo4.c
    viewbeams.c
    ${VERSION_FILE}
  )
  target_link_libraries(rholo rtrad)

  add_executable(rhoptimize rhoptimize.c clumpbeams.c holo.c holofile.c)
  target_link_libraries(rhoptimize rtrad)

  add_executable(rhpict
    holo.c
    holofile.c
    rhpict.c
    rhpict2.c
    viewbeams.c
    ${VERSION_FILE}
  )
  target_link_libraries(rhpict rtrad)

  install(TARGETS rholo rhoptimize rhpict
    RUNTIME DESTINATION "bin"
  )

  if(X11_FOUND)
    add_executable(x11.hdi
      holo.c
      rhd_ctab.c
      rhd_qtree.c
      rhd_qtree2r.c
      rhd_x11.c
      rhdisp.c
      rhdisp2.c
      viewbeams.c
    )
    target_link_libraries(x11.hdi rtrad ${X11_LIBRARIES})

    if(CMAKE_SYSTEM_NAME STREQUAL "SunOS")
      target_link_libraries(x11.hdi rt)
    endif()
    install(TARGETS x11.hdi
      RUNTIME DESTINATION "bin"
    )
  endif()
endif()

install(TARGETS genrhgrid rhcopy rhinfo
  RUNTIME DESTINATION "bin"
)
