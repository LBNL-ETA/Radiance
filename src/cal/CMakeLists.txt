add_executable(cnt cnt.c)

add_executable(ev ev.c)
target_link_libraries(ev rtrad)

add_executable(histo histo.c)
target_link_libraries(histo ${LIB_M})

add_executable(icalc calc.c)
target_link_libraries(icalc rtrad)

add_executable(neaten neat.c)

add_executable(rcalc rcalc.c)
target_link_libraries(rcalc rtrad)

add_executable(rlam lam.c)
target_link_libraries(rlam rtrad)

add_executable(tabfunc tabfunc.c)
target_link_libraries(tabfunc rtrad)

add_executable(total total.c)
target_link_libraries(total ${LIB_M})

install(TARGETS cnt ev histo icalc neaten rcalc rlam tabfunc total
  RUNTIME DESTINATION bin)
