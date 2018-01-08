add_executable(getbbox
  bbox.c
  getbbox.c
  init2otypes.c
  readobj2.c
)
target_link_libraries(getbbox rtrad)

add_executable(obj2mesh
  cvmesh.c
  o_face.c
  obj2mesh.c
  wfconv.c
  writemesh.c
)
target_link_libraries(obj2mesh rtrad)

add_executable(oconv
  bbox.c
  initotypes.c
  o_cone.c
  o_face.c
  o_instance.c
  oconv.c
  sphere.c
  writeoct.c
)
target_link_libraries(oconv rtrad)

install(TARGETS getbbox obj2mesh oconv
  RUNTIME DESTINATION "bin"
)
