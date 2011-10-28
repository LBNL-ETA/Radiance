if(UNIX)
  find_library(LIB_M m DOC "Path to libm")
  if(NOT LIB_M)
    message(FATAL_ERROR "Cannot build radiance without libm.  Please set LIB_M")
  endif()
else()
  set(LIB_M )
endif()

find_package(OpenGL)
find_package(X11)

include_directories(common)
# create a macro to generate
macro(create_version_file version_file)
  add_custom_command(
    OUTPUT "${version_file}"
    COMMAND "${CMAKE_COMMAND}"
    -DVERSION_OUT_FILE="${version_file}"
    -DVERSION_IN_FILE="${radiance_SOURCE_DIR}/src/rt/VERSION"
    -DVERSION_GOLD="${radiance_SOURCE_DIR}/src/rt/Version.c"
    -P "${radiance_SOURCE_DIR}/src/rt/create_version.cmake"
    )
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
# only use rvu if qt4 is found
if(QT4_FOUND)
  add_subdirectory(qtrvu)
endif()

