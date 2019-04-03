include(setup_paths.cmake)

execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND ./resources/cmake_tests/test_lib.rb test_getinfo
  RESULT_VARIABLE res
)
