if(UNIX)
  find_library(LIB_M m DOC "Path to libm")
  if(NOT LIB_M)
    message(FATAL_ERROR "Cannot build radiance without libm.  Please set LIB_M")
  endif()
else()
  set(LIB_M)
endif()

include_directories(common)
# create a macro to generate
macro(create_version_file version_file)
  add_custom_command(
    OUTPUT "${version_file}"
    COMMAND "${CMAKE_COMMAND}"
      -DRADIANCE_VERSION="${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH}.${PROJECT_VERSION_TWEAK}"
      -DVERSION_OUT_FILE="${version_file}"
      -DVERSION_IN_FILE="${radiance_SOURCE_DIR}/src/rt/VERSION"
      -DVERSION_GOLD="${radiance_SOURCE_DIR}/src/rt/Version.c"
      -P "${radiance_SOURCE_DIR}/src/rt/create_version.cmake"
  )
endmacro()
macro(install_without_extension)
  foreach(file IN ITEMS ${ARGN})
    string(REGEX REPLACE "\\.[^.]*$" "" file_we ${file})
    install(PROGRAMS "${file}"
      DESTINATION "bin"
      RENAME "${file_we}"
    )
  endforeach()
endmacro()
macro(install_and_copy_without_extension)
  foreach(file IN ITEMS ${ARGN})
    string(REGEX REPLACE "\\.[^.]*$" "" file_we ${file})
    install(PROGRAMS "${file}"
      DESTINATION "bin"
      RENAME ${file_we}
    )
    file(COPY "${file}" DESTINATION "${CMAKE_BINARY_DIR}/bin")
    file(RENAME "${CMAKE_BINARY_DIR}/bin/${file}" "${CMAKE_BINARY_DIR}/bin/${file_we}")
  endforeach()
endmacro()
add_subdirectory(common)
add_subdirectory(cal)
add_subdirectory(cv)
add_subdirectory(gen)
add_subdirectory(hd)
add_subdirectory(meta)
add_subdirectory(ot)
add_subdirectory(px)
add_subdirectory(rt)
add_subdirectory(util)
if(NOT BUILD_HEADLESS AND Qt5Widgets_DIR)
  # only use qtrvu if qt5 is found
  add_subdirectory(qtrvu)
endif()
