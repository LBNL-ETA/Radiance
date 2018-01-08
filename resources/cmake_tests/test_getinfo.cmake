include(setup_paths.cmake)

execute_process(
  COMMAND getinfo${CMAKE_EXECUTABLE_SUFFIX} ${resources_dir}/evalglare/testimage.hdr
  OUTPUT_FILE test_output
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from getinfo, res = ${res}")
endif()

file(READ ${resources_dir}/evalglare/testimage.hdr test_output)
if(test_output MATCHES "RADIANCE")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
