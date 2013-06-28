if(X11_FOUND)
  get_filename_component(X11_LIBRARY_DIR ${X11_X11_LIB} PATH)
endif()
if(TIFF_LIBRARY)
  get_filename_component(TIFF_LIBRARY_DIR ${TIFF_LIBRARY} PATH)
  if(WIN32)
    #also need the location of libtiff3.dll on Windows.
    #By default this is ../bin relative to TIFF_LIBRARY_DIR.
    set(TIFF_BINARY_DIR "${TIFF_LIBRARY_DIR}/../bin")
    get_filename_component(TIFF_BINARY_DIR "${TIFF_BINARY_DIR}" ABSOLUTE)
    set(TIFF_LIBRARY_DIRS "${TIFF_LIBRARY_DIR};${TIFF_BINARY_DIR}")
  else()
    set(TIFF_LIBRARY_DIRS "${TIFF_LIBRARY_DIR}")
  endif()
endif()

#plugins loaded by rvu at runtime
get_filename_component(qgif_plugin "${QT_QGIF_${qt_plugin_term}_RELEASE}" NAME)
get_filename_component(qico_plugin "${QT_QICO_${qt_plugin_term}_RELEASE}" NAME)
get_filename_component(qjpeg_plugin "${QT_QJPEG_${qt_plugin_term}_RELEASE}" NAME)
get_filename_component(qmng_plugin "${QT_QMNG_${qt_plugin_term}_RELEASE}" NAME)
get_filename_component(qsvg_plugin "${QT_QSVG_${qt_plugin_term}_RELEASE}" NAME)
get_filename_component(qtiff_plugin "${QT_QTIFF_${qt_plugin_term}_RELEASE}" NAME)

configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/dependencies.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake"
   @ONLY
   )
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake")

