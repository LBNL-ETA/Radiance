include(setup_paths.cmake)

execute_process(
  WORKING_DIRECTORY ${office_dir}
  COMMAND ${perl} ${px_dir}/falsecolor.pl -e -i  test/raytest_rpict.hdr
  OUTPUT_FILE test/raytest_falsecolor.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from falsecolor, res = ${res}")
endif()

file(READ ${office_dir}/test/raytest_falsecolor.hdr test_output)
if(test_output MATCHES "minv")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
