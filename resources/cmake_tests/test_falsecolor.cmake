include(setup_paths.cmake)

execute_process(
  COMMAND ${perl} ${rpath}/falsecolor${CMAKE_PERL_SUFFIX} -e -i ${resources_dir}/evalglare/testimage.hdr
  OUTPUT_FILE ${test_output_dir}/falsecolor_out.hdr
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from falsecolor, res = ${res}")
endif()

file(READ ${test_output_dir}/falsecolor_out.hdr test_output)
if(test_output MATCHES "pcomb")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
