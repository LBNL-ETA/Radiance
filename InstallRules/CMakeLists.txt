configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/dependencies.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake"
   @ONLY
   )
install(SCRIPT "${CMAKE_CURRENT_BINARY_DIR}/dependencies.cmake")

