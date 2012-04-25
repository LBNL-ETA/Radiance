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
configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/dependencies.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake"
   @ONLY
   )
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake")

