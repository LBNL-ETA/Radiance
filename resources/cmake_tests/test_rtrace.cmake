include(setup_paths.cmake)
file(WRITE ${test_output_dir}/rtrace.in
"0 0 0 0 0 -1
")
execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND rtrace${CMAKE_EXECUTABLE_SUFFIX} -ab 1 cornell_box.oct
  INPUT_FILE rtrace.in
  OUTPUT_FILE rtrace.out
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from rtrace, res = ${res}")
endif()

file(READ ${test_output_dir}/rtrace.out test_output)
if(test_output MATCHES "SOFTWARE= RADIANCE")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
