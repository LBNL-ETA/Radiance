include(setup_paths.cmake)

execute_process(
  WORKING_DIRECTORY ${test_output_dir}
  COMMAND evalglare${CMAKE_EXECUTABLE_SUFFIX} -L 400 500 0.6 1.2 -B 0.4 -b 2500 -d -c outputimage.hdr ${resources_dir}/evalglare/testimage.hdr
  OUTPUT_FILE evalglare_out.txt
  RESULT_VARIABLE res
)

file(READ ${test_output_dir}/evalglare_out.txt test_output)
if(test_output MATCHES "band:band_omega,band_av_lum,band_median_lum,band_std_lum,band_perc_75,band_perc_95,band_lum_min,band_lum_max: 2.423399 2353.496003 1339.703125 2218.441248 3993.937463 5750.375098 68.348633 37142.499489")
  message(STATUS "passed")
else()
  message(STATUS "failed")
endif()
