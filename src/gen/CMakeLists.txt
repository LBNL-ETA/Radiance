# Locate the threads package
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../rt)

add_executable(genbeads genbeads.c hermite3.c)
target_link_libraries(genbeads ${LIB_M})

add_executable(genbox genbox.c)
target_link_libraries(genbox wfobj rtrad ${LIB_M})

add_executable(genmarble genmarble.c ../common/random.h)
target_link_libraries(genmarble ${LIB_M})

add_executable(gensky gensky.c sun.c)
target_link_libraries(gensky rtrad ${LIB_M})


add_executable(genssky genssky.c sun.c atmos.c ../rt/data.c)
target_link_libraries(genssky rtrad ${LIB_M} Threads::Threads)

add_executable(gendaylit gendaylit.c sun.c)
target_link_libraries(gendaylit rtrad ${LIB_M})

add_executable(gendaymtx gendaymtx.c sun.c)
target_link_libraries(gendaymtx rtrad ${LIB_M})

add_executable(genblinds genblinds.c)
target_link_libraries(genblinds rtrad ${LIB_M})

add_executable(genprism genprism.c)
target_link_libraries(genprism rtrad ${LIB_M})

add_executable(xform xform.c)
target_link_libraries(xform rtrad)

add_executable(replmarks replmarks.c)
target_link_libraries(replmarks rtrad)

add_executable(mkillum mkillum.c mkillum2.c mkillum3.c)
target_link_libraries(mkillum raycalls radiance rtrad)

add_executable(mksource mksource.c)
target_link_libraries(mksource raycalls radiance rtrad)

add_executable(genrev genrev.c)
target_link_libraries(genrev rtrad)

add_executable(gencatenary gencat.c)
target_link_libraries(gencatenary rtrad)

add_executable(genworm genworm.c)
target_link_libraries(genworm rtrad)

add_executable(gensurf gensurf.c)
target_link_libraries(gensurf rtrad)

add_executable(genclock genclock.c)
target_link_libraries(genclock ${LIB_M})

add_executable(genbranch genbranch.c)
target_link_libraries(genbranch ${LIB_M})

install(TARGETS
  genblinds
  genbox
  genclock
  gendaylit
  gendaymtx
  genprism
  genrev
  gensky
  genssky
  gensurf
  genworm
  mkillum
  mksource
  replmarks
  xform
  DESTINATION "bin"
)

file(COPY
  clockface.hex
  glaze1.cal
  glaze2.cal
  illum.cal
  perezlum.cal
  perezlum_c.cal
  rev.cal
  skybright.cal
  surf.cal
  mie_ca.dat
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)

if(NOT WIN32)
  install(PROGRAMS glaze.csh
    DESTINATION "bin"
    RENAME glaze
  )
endif()
