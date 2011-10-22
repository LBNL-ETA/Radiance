add_executable(getbbox getbbox.c readobj2.c bbox.c init2otypes.c)
target_link_libraries(getbbox rtrad)

add_executable(obj2mesh obj2mesh.c cvmesh.c wfconv.c o_face.c writemesh.c)
target_link_libraries(obj2mesh rtrad)

add_executable(oconv oconv.c sphere.c writeoct.c o_face.c o_cone.c o_instance.c bbox.c initotypes.c)
target_link_libraries(oconv rtrad)

install(TARGETS oconv getbbox obj2mesh RUNTIME DESTINATION bin)
