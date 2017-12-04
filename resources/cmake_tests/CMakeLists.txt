include(FindPerl)

set(office_dir
  "${radiance_BINARY_DIR}/Downloads/Source/radiance_support/obj/office")
set(test_dir "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(setup_paths.cmake.in setup_paths.cmake @ONLY)
# configure these files so that they will be in the binary
# directory.  This will allow them to include setup_paths.cmake
# from the binary directory and just use regular variables in the
# configured files.
configure_file(test_setup.cmake test_setup.cmake COPYONLY)
configure_file(test_rtrace.cmake test_rtrace.cmake COPYONLY)
configure_file(test_rpict.cmake test_rpict.cmake COPYONLY)
configure_file(test_getinfo.cmake test_getinfo.cmake COPYONLY)
configure_file(test_pcond.cmake test_pcond.cmake COPYONLY)
configure_file(test_falsecolor.cmake test_falsecolor.cmake COPYONLY)
configure_file(test_DC.cmake test_DC.cmake COPYONLY)

add_test(test_setup ${CMAKE_COMMAND}
  -P ${CMAKE_CURRENT_BINARY_DIR}/test_setup.cmake)
set_tests_properties(test_setup PROPERTIES FAIL_REGULAR_EXPRESSION
  "No such file|not found|cannot find file|unexpected character")

add_test(test_rtrace
  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_rtrace.cmake)
set_tests_properties(test_rtrace PROPERTIES PASS_REGULAR_EXPRESSION "passed"
                                            FAIL_REGULAR_EXPRESSION "failed"
                                            DEPENDS "setup")
add_test(test_rpict
  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_rpict.cmake)
set_tests_properties(test_rpict PROPERTIES PASS_REGULAR_EXPRESSION "passed"
                                            FAIL_REGULAR_EXPRESSION "failed"
                                            DEPENDS "setup")
add_test(test_getinfo
  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_getinfo.cmake)
add_test(test_pcond
  ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_pcond.cmake)
set_tests_properties(test_pcond PROPERTIES PASS_REGULAR_EXPRESSION "passed"
                                            FAIL_REGULAR_EXPRESSION "failed"
                                            DEPENDS "setup")
if(PERL_FOUND)
  add_test(test_falsecolor
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_falsecolor.cmake)
  set_tests_properties(test_falsecolor PROPERTIES PASS_REGULAR_EXPRESSION "passed"
                                              FAIL_REGULAR_EXPRESSION "failed"
                                              DEPENDS "setup")
  add_test(test_DC
    ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/test_DC.cmake)
  set_tests_properties(test_DC PROPERTIES PASS_REGULAR_EXPRESSION "passed"
                                              FAIL_REGULAR_EXPRESSION "failed"
                                              DEPENDS "setup")
endif()
