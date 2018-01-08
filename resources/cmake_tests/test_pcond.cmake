include(setup_paths.cmake)

execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND pcond${CMAKE_EXECUTABLE_SUFFIX} -h rpict_out.hdr
  OUTPUT_FILE pcond_out.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from pcond, res = ${res}")
endif()

file(READ ${test_output_dir}/pcond_out.hdr test_output)
if(test_output MATCHES "pcond")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
