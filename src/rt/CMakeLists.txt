set(VERSION_FILE "${radiance_BINARY_DIR}/src/rt/Version.c")
create_version_file("${VERSION_FILE}")

add_library(radiance
  "${VERSION_FILE}"
  ../common/paths.c
  ../common/platform.h
  ../common/random.h
  ambcomp.c
  ambient.c
  ambio.c
  aniso.c
  ashikhmin.c
  data.c
  dielectric.c
  fprism.c
  freeobjmem.c
  func.c
  glass.c
  initotypes.c
  m_alias.c
  m_brdf.c
  m_bsdf.c
  m_clip.c
  m_direct.c
  m_mirror.c
  m_mist.c
  mx_data.c
  mx_func.c
  noise3.c
  normal.c
  o_cone.c
  o_face.c
  o_instance.c
  o_mesh.c
  p_data.c
  p_func.c
  pmap.c
  pmapamb.c
  pmapbias.c
  pmapcontrib.c
  pmapdata.c
  pmapdiag.c
  pmapio.c
  pmapmat.c
  pmapopt.c
  pmapparm.c
  pmaprand.c
  pmapray.c
  pmapsrc.c
  pmaptype.c
  pmcontrib2.c
  pmutil.c
  preload.c
  raytrace.c
  renderopts.c
  source.c
  sphere.c
  srcobstr.c
  srcsamp.c
  srcsupp.c
  t_data.c
  t_func.c
  text.c
  virtuals.c
)

if(WIN32)
  set(rayp_SOURCES raypwin.c)
else()
  set(rayp_SOURCES raypcalls.c)
endif()

add_library(raycalls raycalls.c ${rayp_SOURCES} rayfifo.c)

add_executable(rtrace rtmain.c rtrace.c duphead.c persist.c)
target_link_libraries(rtrace raycalls radiance rtrad)

add_executable(rpict rpmain.c rpict.c srcdraw.c duphead.c persist.c)
target_link_libraries(rpict radiance rtrad)

add_executable(lookamb lookamb.c ambio.c)
target_link_libraries(lookamb rtrad)

add_executable(rcontrib rcmain.c rcontrib.c rc2.c rc3.c)
target_link_libraries(rcontrib radiance rtrad)

add_executable(mkpmap mkpmap.c)
target_link_libraries(mkpmap radiance rtrad)
add_executable(pmapdump pmapdump.c pmaptype.c pmapparm.c)
target_link_libraries(pmapdump raycalls radiance rtrad)

set(targets_to_install
  lookamb
  mkpmap
  pmapdump
  radiance
  raycalls
  rcontrib
  rpict
  rtrace
)

if(X11_FOUND)
  set(X11_SOURCES x11.c x11twind.c)
  add_definitions(-DHAS_X11)
endif()
if(Qt5Widgets_DIR)
  include_directories(${Qt5Widgets_INCLUDE_DIRS})
  add_definitions(-DHAS_QT)
endif()
if(NOT BUILD_HEADLESS)
  if(X11_FOUND OR Qt5Widgets_DIR)
    add_executable(rvu
      colortab.c
      devcomm.c
      devtable.c
      editline.c
      rv2.c
      rv3.c
      rview.c
      rvmain.c
      ${X11_SOURCES}
    )
    target_link_libraries(rvu raycalls radiance rtrad ${X11_LIBRARIES})
    list(APPEND targets_to_install rvu)

    # only link to rvu-qt-lib if qt is found
    if(Qt5Widgets_DIR)
      target_link_libraries(rvu qtrvu)
    endif()
  endif()
endif()

#add_test(rtrace_version "${EXECUTABLE_OUTPUT_PATH}/rtrace" -version)
#set_tests_properties(rtrace_version PROPERTIES
#  PASS_REGULAR_EXPRESSION "RADIANCE"
#)

install(TARGETS ${targets_to_install}
  RUNTIME DESTINATION "bin"
  LIBRARY DESTINATION "bin"
  ARCHIVE DESTINATION "bin"
)
file(COPY rayinit.cal
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)

#install qt.conf & Qt plugins loaded by rvu at runtime. Ignore if we are building as headless
if(NOT BUILD_HEADLESS)
  get_target_property(GIF_PLUGIN Qt5::QGifPlugin LOCATION)
  get_target_property(ICO_PLUGIN Qt5::QICOPlugin LOCATION)
  get_target_property(JPEG_PLUGIN Qt5::QJpegPlugin LOCATION)
  get_target_property(TIFF_PLUGIN Qt5::QTiffPlugin LOCATION)
  install(PROGRAMS "${GIF_PLUGIN}" "${ICO_PLUGIN}" "${JPEG_PLUGIN}" "${TIFF_PLUGIN}"
    DESTINATION "bin/plugins/imageformats"
  )
  install(FILES "${CMAKE_SOURCE_DIR}/InstallRules/qt.conf"
    DESTINATION "bin"
  )

  if(WIN32)
    get_target_property(QPA_PLUGIN Qt5::QWindowsIntegrationPlugin LOCATION)
  elseif(APPLE)
    get_target_property(QPA_PLUGIN Qt5::QCocoaIntegrationPlugin LOCATION)
  elseif(UNIX)
    get_target_property(QPA_PLUGIN Qt5::QXcbIntegrationPlugin LOCATION)
  endif()

  install(PROGRAMS "${QPA_PLUGIN}"
    DESTINATION "bin/plugins/platforms"
  )
endif()
