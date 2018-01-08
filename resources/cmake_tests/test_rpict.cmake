include(setup_paths.cmake)

# run rpict on the Cornell box
execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND rpict${CMAKE_EXECUTABLE_SUFFIX} -ab 1 -vf ${resources_dir}/cornell_box/cornell.vf ${test_output_dir}/cornell_box.oct
  OUTPUT_FILE ${test_output_dir}/rpict_out.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from rpict, res = ${res}")
endif()

# scan output file for valid header
file(READ ${test_output_dir}/rpict_out.hdr test_output)
if(test_output MATCHES "SOFTWARE= RADIANCE")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
