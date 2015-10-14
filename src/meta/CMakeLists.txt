set(common_src mfio.c syscalls.c misc.c)

add_executable(bgraph bgraph.c mgvars.c mgraph.c)
target_link_libraries(bgraph meta rtrad)

add_executable(cv cv.c mfio.c cvhfio.c syscalls.c misc.c)
target_link_libraries(cv rtrad)

add_executable(dgraph dgraph.c cgraph.c mgvars.c)
target_link_libraries(dgraph rtrad)

add_executable(gcomp gcomp.c gcalc.c mgvars.c)
target_link_libraries(gcomp rtrad)

add_executable(igraph igraph.c mgvars.c mgraph.c gcalc.c cgraph.c)
target_link_libraries(igraph meta rtrad)

add_library(meta metacalls.c primout.c progname.c ${common_src})

add_executable(meta2bmp meta2bmp.c rplot.c plot.c palloc.c ${common_src})
target_link_libraries(meta2bmp rtrad)

add_executable(meta2tga meta2tga.c rplot.c plot.c palloc.c ${common_src})
target_link_libraries(meta2tga rtrad)

add_executable(mtext mtext.c primout.c ${common_src})
target_link_libraries(mtext rtrad)

add_executable(pexpand pexpand.c expand.c segment.c palloc.c ${common_src})
target_link_libraries(pexpand rtrad)

add_executable(plot4 plot4.c primout.c ${common_src})
target_link_libraries(plot4 rtrad)

add_executable(plotin plotin.c primout.c ${common_src})
target_link_libraries(plotin rtrad)

add_executable(psmeta psmeta.c psplot.c ${common_src})
target_link_libraries(psmeta rtrad)

add_executable(psort psort.c sort.c palloc.c ${common_src})
target_link_libraries(psort rtrad)

#add_executable(tbar tbar.c tgraph.c primout.c)
#target_link_libraries(tbar rtrad)

#add_executable(tcurve tcurve.c tgraph.c primout.c)
#target_link_libraries(tcurve ${LIB_M})

#add_executable(tscat tscat.c tgraph.c primout.c)
#target_link_libraries(tscat rtrad)

if(X11_FOUND)
  add_executable(x11meta xmeta.c x11plot.c plot.c palloc.c ${common_src})
  target_link_libraries(x11meta rtrad ${X11_LIBRARIES})
  install(TARGETS x11meta RUNTIME DESTINATION bin)
endif()

add_custom_target(mta ALL
  COMMAND cv "${CMAKE_CURRENT_SOURCE_DIR}/symbols.met" > "${CMAKE_CURRENT_BINARY_DIR}/symbols.mta"
  COMMAND cv "${CMAKE_CURRENT_SOURCE_DIR}/vchars.met" > "${CMAKE_CURRENT_BINARY_DIR}/vchars.mta"
  SOURCES symbols.met vchars.met
)

install(FILES 
  "${CMAKE_CURRENT_BINARY_DIR}/symbols.mta"
  "${CMAKE_CURRENT_BINARY_DIR}/vchars.mta"
  DESTINATION lib
)

install(FILES symbols.met vchars.met boxw.plt curve.plt line.plt scatter.plt 
				cartesian.plt function.plt polar.plt standard.plt
				DESTINATION lib)

install(TARGETS meta2tga pexpand psort cv psmeta plotin bgraph igraph dgraph
        gcomp plot4 meta2bmp
        RUNTIME DESTINATION bin)

