add_subdirectory(cmake_tests)

file(COPY fonts/helvet.fnt
  DESTINATION "${CMAKE_BINARY_DIR}/lib"
)
