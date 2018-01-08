include(setup_paths.cmake)

file(REMOVE_RECURSE ${test_output_dir})
file(MAKE_DIRECTORY ${test_output_dir})

execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND oconv${CMAKE_EXECUTABLE_SUFFIX} -f ${resources_dir}/cornell_box/cornell.rad
  OUTPUT_FILE cornell_box.oct
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from oconv, res = ${res}")
endif()
