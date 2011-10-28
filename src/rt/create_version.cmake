#${CMAKE_COMMAND} -DVERSION_OUT_FILE=v  -DVERSION_IN_FILE=src/rt/VERSION -DVERSION_GOLD=src/rt/Version.c -P src/common/create_version.cmake

# if the gold version exists then use that instead
if(EXISTS "${VERSION_GOLD}")
  configure_file("${VERSION_GOLD}" "${VERSION_OUT_FILE}" COPYONLY)
  return()
endif()

find_program(DATE date)
if(DATE)
  execute_process(COMMAND ${DATE} OUTPUT_VARIABLE DATE_STR)
  string(STRIP "${DATE_STR}" DATE_STR)
endif()
find_program(WHO whoami)
if(WHO)
  execute_process(COMMAND ${WHO} OUTPUT_VARIABLE WHO_STR)
  string(STRIP "${WHO_STR}" WHO_STR)
endif()
find_program(HOSTNAME hostname)
if(HOSTNAME)
  execute_process(COMMAND ${HOSTNAME} OUTPUT_VARIABLE HOST_STR)
  message("DATE= ${DATE_STR}")
  string(STRIP "${HOST_STR}" HOST_STR)
endif()
file(READ "${VERSION_IN_FILE}" VERSION)
string(STRIP "${VERSION}" VERSION)
message("${VERSION}")
file(WRITE "${VERSION_OUT_FILE}"
"char VersionID[]=\"${VERSION} lastmod ${DATE_STR}"
" by ${WHO_STR} on ${HOST_STR}\";\n"
)

# todo actually get the date, user name and hostname in this script
