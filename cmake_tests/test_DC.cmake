include(setup_paths.cmake)

file(MAKE_DIRECTORY "${test_output_dir}/dc")

set(view_def -x 50 -y 50 -vf model.vp)
message(STATUS "view_def=${view_def}")

file(WRITE ${office_dir}/test/raytest_dc.sky 
"void glow skyglow\n0\n0\n4\n1 1 1 0\n\nskyglow source sky\n0\n0\n4\n0 0 1 360
")

message("Run: oconv${CMAKE_EXECUTABLE_SUFFIX} -f -i test/raytest_model.oct test/raytest_dc.sky")
execute_process(
  WORKING_DIRECTORY ${office_dir}
  COMMAND oconv${CMAKE_EXECUTABLE_SUFFIX} -f -i test/raytest_model.oct test/raytest_dc.sky
  OUTPUT_FILE test/raytest_model_dc.oct
)

execute_process(
  WORKING_DIRECTORY ${office_dir}
  COMMAND vwrays -d ${view_def}
  OUTPUT_VARIABLE vwrays_out
)

message(STATUS "vwrays_out: ${vwrays_out}")
string(STRIP "${vwrays_out}" vwrays_out)
separate_arguments(vwrays_out)

message("Run: vwrays${CMAKE_EXECUTABLE_SUFFIX} -ff ${view_def} | rtcontrib${CMAKE_EXECUTABLE_SUFFIX} -ab 1 ${rtcontrib_threads}-V+ -fo -ffc $\(vwrays -d ${view_def}\) -f tregenza.cal -bn 146 -b tbin -o ${test_output_dir}/dc/treg%03d.hdr -m skyglow test/raytest_model_dc.oct")
execute_process(
  WORKING_DIRECTORY ${office_dir}
  COMMAND vwrays${CMAKE_EXECUTABLE_SUFFIX} -ff ${view_def} 
  COMMAND rtcontrib${CMAKE_EXECUTABLE_SUFFIX} -ab 1 ${rtcontrib_threads}-V+ -fo -ffc ${vwrays_out} -f tregenza.cal -bn 146 -b tbin -o ${test_output_dir}/dc/treg%03d.hdr -m skyglow test/raytest_model_dc.oct
  RESULT_VARIABLE res
)
if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value from vwrays, res = ${res}")
endif() 

execute_process( 
  WORKING_DIRECTORY ${office_dir}
  COMMAND gensky${CMAKE_EXECUTABLE_SUFFIX} 03 21 12 -a 40 -o 105 -m 105 +s
  COMMAND ${perl} ${util_dir}/genskyvec.pl -m 1
  COMMAND dctimestep${CMAKE_EXECUTABLE_SUFFIX} test/dc/treg%03d.hdr
  COMMAND pfilt${CMAKE_EXECUTABLE_SUFFIX} -1 -x /1 -y /1
  OUTPUT_FILE test/raytest_timestep.hdr
  RESULT_VARIABLE res
)

if(NOT ${res} EQUAL 0)
  message(FATAL_ERROR "Bad return value, res = ${res}")
endif()

file(READ ${office_dir}/test/raytest_timestep.hdr test_output)
if(test_output MATCHES "dctimestep${CMAKE_EXECUTABLE_SUFFIX} test/dc/treg%03d.hdr")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
