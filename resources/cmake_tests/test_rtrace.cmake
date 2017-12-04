include(setup_paths.cmake)
file(WRITE ${office_dir}/test/raytest_rtrace.in
"-23 24.6 5.8 0 1 0
")
execute_process(
  COMMAND rtrace${CMAKE_EXECUTABLE_SUFFIX} -ab 1 test/raytest_model.oct
  WORKING_DIRECTORY ${office_dir}
  INPUT_FILE test/raytest_rtrace.in
  OUTPUT_FILE test/raytest_rtrace.out
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from rtrace, res = ${res}")
endif()

file(READ ${office_dir}/test/raytest_rtrace.out test_output)
if(test_output MATCHES "SOFTWARE= RADIANCE")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
