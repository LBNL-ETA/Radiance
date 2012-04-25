include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../rt)

add_executable(genbeads genbeads.c hermite3.c)
target_link_libraries(genbeads ${LIB_M})

add_executable(genblinds genblinds.c)
target_link_libraries(genblinds ${LIB_M})

add_executable(genbox genbox.c)

add_executable(genbranch genbranch.c)
target_link_libraries(genbranch ${LIB_M})

add_executable(gencatenary gencat.c)
target_link_libraries(gencatenary rtrad)

add_executable(genclock genclock.c)
target_link_libraries(genclock ${LIB_M})

add_executable(gendaylit gendaylit.c sun.c)
target_link_libraries(gendaylit rtrad)

add_executable(genmarble genmarble.c ../common/random.h)
target_link_libraries(genmarble ${LIB_M})

add_executable(genprism genprism.c)
target_link_libraries(genprism ${LIB_M})

add_executable(genrev genrev.c )
target_link_libraries(genrev rtrad)

add_executable(gensky gensky.c sun.c)
target_link_libraries(gensky ${LIB_M})

add_executable(gensurf gensurf.c)
target_link_libraries(gensurf rtrad)

add_executable(genworm genworm.c)
target_link_libraries(genworm rtrad)

add_executable(mkillum mkillum.c mkillum2.c mkillum3.c)
target_link_libraries(mkillum raycalls radiance rtrad)

add_executable(mksource mksource.c)
target_link_libraries(mksource raycalls radiance rtrad)

add_executable(replmarks replmarks.c)
target_link_libraries(replmarks rtrad)

add_executable(xform xform.c)
target_link_libraries(xform rtrad)

install(TARGETS genbox gensky xform genblinds genrev genworm gendaylit gensurf
                genprism replmarks mkillum mksource genclock
        DESTINATION bin)

install(FILES illum.cal rev.cal skybright.cal surf.cal clockface.hex glaze1.cal
              glaze2.cal perezlum.cal coeff_perez.dat defangle.dat
        DESTINATION lib)

if (NOT WIN32)
  install(PROGRAMS glaze.csh DESTINATION bin RENAME glaze)
endif()