#${CMAKE_COMMAND} -DRADIANCE_VERSION=v -DVERSION_OUT_FILE=v -DVERSION_IN_FILE=src/rt/VERSION -DVERSION_GOLD=src/rt/Version.c -P src/common/create_version.cmake

# if the gold version exists then use that instead
if(EXISTS "${VERSION_GOLD}")
  configure_file("${VERSION_GOLD}" "${VERSION_OUT_FILE}" COPYONLY)
  return()
endif()

find_program(DATE date)
if(DATE)
  execute_process(COMMAND ${DATE} "+%F"
    OUTPUT_VARIABLE DATE_STR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()
file(READ "${VERSION_IN_FILE}" VERSION)
string(STRIP "${VERSION}" VERSION)
message("${VERSION}")
file(WRITE "${VERSION_OUT_FILE}"
  "char VersionID[]=\"${VERSION} ${DATE_STR} LBNL (${RADIANCE_VERSION})\";\n"
)
