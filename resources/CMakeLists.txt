add_subdirectory(cmake_tests)

file(COPY fonts/helvet.fnt
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)
file(COPY test_lib.rb
  DESTINATION "${CMAKE_BINARY_DIR}/resources"
)
