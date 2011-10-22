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

