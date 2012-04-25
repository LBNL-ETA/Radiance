include(setup_paths.cmake)

execute_process(
  WORKING_DIRECTORY ${office_dir}
  COMMAND rpict${CMAKE_EXECUTABLE_SUFFIX} -ab 0 -vf model.vp test/raytest_model.oct
  OUTPUT_FILE test/raytest_rpict.hdr
  COMMAND pcond${CMAKE_EXECUTABLE_SUFFIX} -h test/raytest_rpict.hdr
  OUTPUT_FILE test/raytest_rpict_pcd.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from pcond, res = ${res}")
endif()

file(READ ${office_dir}/test/raytest_rpict_pcd.hdr test_output)
if(test_output MATCHES "pcond")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()