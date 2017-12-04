include(setup_paths.cmake)

# run rpict on office model
execute_process(
  COMMAND rpict${CMAKE_EXECUTABLE_SUFFIX} -ab 0 -vf model.vp test/raytest_model.oct
  WORKING_DIRECTORY ${office_dir}
  OUTPUT_FILE test/raytest_rpict.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from rpict, res = ${res}")
endif()

# scan output file for valid header
file(READ ${office_dir}/test/raytest_rpict.hdr test_output)
if(test_output MATCHES "SOFTWARE= RADIANCE")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
