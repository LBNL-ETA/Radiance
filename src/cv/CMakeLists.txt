add_executable(ies2rad ies2rad.c)
target_link_libraries(ies2rad rtrad)

add_executable(arch2rad arch2rad.c trans.c)
target_link_libraries(arch2rad rtrad)

add_executable(nff2rad nff2rad.c)

add_executable(lampcolor lampcolor.c)
target_link_libraries(lampcolor rtrad)

add_executable(tmesh2rad tmesh2rad.c)
target_link_libraries(tmesh2rad rtrad)

add_executable(obj2rad obj2rad.c trans.c)
target_link_libraries(obj2rad rtrad)

add_executable(mgf2rad mgf2rad.c)
target_link_libraries(mgf2rad mgf rtrad)

add_executable(rad2mgf rad2mgf.c)
target_link_libraries(rad2mgf rtrad)

add_executable(mgfilt mgfilt.c)
target_link_libraries(mgfilt mgf rtrad)

add_executable(mgf2inv mgf2inv.c)
target_link_libraries(mgf2inv mgf rtrad)

add_executable(3ds2mgf 3ds2mgf.c rayopt3ds.c vect3ds.c)
target_link_libraries(3ds2mgf rtrad)

if(BUILD_PABOPTO_UTILS)
  add_executable(pabopto2bsdf pabopto2bsdf.c bsdfrep.c bsdfrbf.c bsdfmesh.c bsdftrans.cpp)
  target_link_libraries(pabopto2bsdf rtrad)

  #add_executable(pabopto2xyz pabopto2xyz.c)
  #target_link_libraries(pabopto2xyz rtrad)

  add_executable(pabopto2rad pabopto2bsdf.c bsdfrep.c bsdfrbf.c bsdfmesh.c bsdftrans.cpp)
  target_link_libraries(pabopto2rad rtrad)
  set_target_properties(pabopto2rad PROPERTIES COMPILE_DEFINITIONS TEST_MAIN)
endif()

add_executable(bsdf2rad bsdf2rad.c bsdfrep.c bsdfinterp.c)
target_link_libraries(bsdf2rad rtrad)

add_executable(bsdfquery bsdfquery.c bsdfrep.c bsdfinterp.c)
target_link_libraries(bsdfquery rtrad)

add_executable(bsdf2ttree bsdf2ttree.c bsdfrep.c bsdfinterp.c)
target_link_libraries(bsdf2ttree rtrad)

add_executable(bsdf2klems bsdf2klems.c bsdfrep.c bsdfinterp.c)
target_link_libraries(bsdf2klems rtrad)

add_executable(pkgBSDF pkgBSDF.c trans.c)
target_link_libraries(pkgBSDF rtrad)

add_executable(epw2wea epw2wea.c)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../meta)
add_executable(mgf2meta mgf2meta.c)
target_link_libraries(mgf2meta mgf meta rtrad)

install(TARGETS
  3ds2mgf
  bsdf2klems
  bsdf2ttree
  epw2wea
  ies2rad
  lampcolor
  mgf2inv
  mgf2meta
  mgf2rad
  mgfilt
  nff2rad
  obj2rad
  pkgBSDF
  rad2mgf
  tmesh2rad
  RUNTIME DESTINATION "bin"
)

file(COPY
  bsdf2rad.cal
  lamp.tab
  source.cal
  tilt.cal
  window.cal
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)

if(NOT WIN32)
  install(PROGRAMS optics2rad.csh
    DESTINATION "bin"
    RENAME optics2rad
  )
endif()

if(BUILD_PABOPTO_UTILS)
  install(TARGETS pabopto2bsdf pabopto2rad
    DESTINATION "bin"
  )
endif()
