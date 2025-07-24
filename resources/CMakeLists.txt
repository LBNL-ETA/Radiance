# RCSid $Id: CMakeLists.txt,v 1.3 2025/07/24 01:24:01 greg Exp $

add_subdirectory(cmake_tests)

file(COPY fonts/helvet.fnt
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)
file(COPY test_lib.rb
  DESTINATION "${CMAKE_BINARY_DIR}/resources"
)
