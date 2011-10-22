include(setup_paths.cmake) # configured file found in the same dir as this file

file(REMOVE_RECURSE ${test_output_dir})
file(MAKE_DIRECTORY ${test_output_dir})

message("Run: [oconv${CMAKE_EXECUTABLE_SUFFIX} model.b90 desk misc]")
execute_process(
  COMMAND oconv${CMAKE_EXECUTABLE_SUFFIX} model.b90 desk misc
  WORKING_DIRECTORY ${office_dir}
  OUTPUT_FILE test/raytest_modelb.oct
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from oconv, res = ${res}")
endif()

message("Run:
  [oconv${CMAKE_EXECUTABLE_SUFFIX} -f -i test/raytest_modelb.oct
  window blinds lights lamp picture]")
execute_process(
  COMMAND oconv${CMAKE_EXECUTABLE_SUFFIX}
  -f -i test/raytest_modelb.oct window blinds lights lamp picture
  WORKING_DIRECTORY ${office_dir}
  OUTPUT_FILE test/raytest_model.oct
  RESULT_VARIABLE res
)

if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from oconv, res = ${res}")
endif()
