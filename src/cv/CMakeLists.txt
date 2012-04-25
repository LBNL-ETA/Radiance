add_executable(arch2rad arch2rad.c trans.c)
target_link_libraries(arch2rad rtrad)

add_executable(ies2rad ies2rad.c)
target_link_libraries(ies2rad rtrad)

add_executable(lampcolor lampcolor.c)
target_link_libraries(lampcolor rtrad)

add_executable(mgf2rad mgf2rad.c)
target_link_libraries(mgf2rad mgf rtrad)

add_executable(nff2rad nff2rad.c)

add_executable(obj2rad obj2rad.c trans.c)
target_link_libraries(obj2rad rtrad)

add_executable(pkgBSDF pkgBSDF.c trans.c)
target_link_libraries(pkgBSDF rtrad)

add_executable(rad2mgf rad2mgf.c)
target_link_libraries(rad2mgf rtrad)

add_executable(tmesh2rad tmesh2rad.c)
target_link_libraries(tmesh2rad rtrad)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../meta)
add_executable(mgf2meta mgf2meta.c)
target_link_libraries(mgf2meta mgf meta rtrad)

add_executable(mgfilt mgfilt.c)
target_link_libraries(mgfilt mgf rtrad)

add_executable(3ds2mgf 3ds2mgf.c rayopt3ds.c vect3ds.c)
target_link_libraries(3ds2mgf rtrad) 

add_executable(mgf2inv mgf2inv.c)
target_link_libraries(mgf2inv mgf rtrad) 

install(TARGETS 
    ies2rad nff2rad lampcolor tmesh2rad obj2rad mgf2rad
    rad2mgf mgf2meta mgfilt mgf2inv pkgBSDF 3ds2mgf
  RUNTIME DESTINATION bin)

INSTALL(FILES source.cal tilt.cal lamp.tab window.cal
  DESTINATION lib)

if (NOT WIN32)
  INSTALL(PROGRAMS optics2rad.csh DESTINATION bin RENAME optics2rad)
endif()