add_executable(pfilt pf2.c pf3.c pfilt.c)
target_link_libraries(pfilt rtrad)

add_executable(pcond pcond.c pcond2.c pcond3.c pcond4.c warp3d.c)
target_link_libraries(pcond rtrad)

add_executable(pcwarp pcwarp.c warp3d.c)
target_link_libraries(pcwarp rtrad)

add_executable(ttyimage ttyimage.c)
target_link_libraries(ttyimage rtrad)

add_executable(pvalue pvalue.c)
target_link_libraries(pvalue rtrad)

add_executable(pcompos pcompos.c)
target_link_libraries(pcompos rtrad)

add_executable(psign psign.c)
target_link_libraries(psign rtrad)

add_executable(ra_gif ra_gif.c clrtab.c neuclrtab.c)
target_link_libraries(ra_gif rtrad)

add_executable(ra_ps ra_ps.c)
target_link_libraries(ra_ps rtrad)

add_executable(ra_ppm ra_ppm.c)
target_link_libraries(ra_ppm rtrad)

add_executable(ra_bmp ra_bmp.c)
target_link_libraries(ra_bmp rtrad)

add_executable(ra_t8 ra_t8.c clrtab.c neuclrtab.c)
target_link_libraries(ra_t8 rtrad)

add_executable(ra_t16 ra_t16.c)
target_link_libraries(ra_t16 rtrad)

add_executable(ra_rgbe ra_rgbe.c)
target_link_libraries(ra_rgbe rtrad)

add_executable(ra_pict ra_pict.c)
target_link_libraries(ra_pict rtrad)

add_executable(ra_hexbit ra_hexbit.c)
target_link_libraries(ra_hexbit rtrad)

if(X11_FOUND)
  add_executable(ximage x11image.c x11raster.c clrtab.c)
  target_link_libraries(ximage rtrad ${X11_LIBRARIES})
  add_executable(xshowtrace
    xshowtrace.c
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/x11findwind.c"
  )
  target_link_libraries(xshowtrace rtrad ${X11_LIBRARIES})
  install(TARGETS ximage xshowtrace
    RUNTIME DESTINATION "bin"
  )
endif()

add_executable(pflip pflip.c)
target_link_libraries(pflip rtrad)

add_executable(protate protate.c)
target_link_libraries(protate rtrad)

add_executable(pextrem pextrem.c)
target_link_libraries(pextrem rtrad)

add_executable(pcomb pcomb.c)
target_link_libraries(pcomb rtrad)

add_executable(pinterp pinterp.c)
target_link_libraries(pinterp rtrad)

add_executable(ra_xyze ra_xyze.c)
target_link_libraries(ra_xyze rtrad)

add_executable(macbethcal macbethcal.c pmapgen.c mx3.c warp3d.c)
target_link_libraries(macbethcal rtrad)

add_executable(pmblur2 pmblur2.c)
target_link_libraries(pmblur2 rtrad)

add_executable(psketch psketch.c)
target_link_libraries(psketch rtrad)

if(NOT TIFF_LIBRARY)
  find_package(TIFF)
endif()
mark_as_advanced(CLEAR TIFF_INCLUDE_DIR TIFF_LIBRARY)

#special case for manually-built libtiff where find_package fails
if(NOT TIFF_FOUND)
  find_library(TIFF_LIBRARY
    NAMES libtiff.a tiff tiff3 libtiff
    PATHS
      "${CMAKE_CURRENT_SOURCE_DIR}/tiff/libtiff/.libs"
      "${radiance_BINARY_DIR}/Downloads/Source/radiance_support/src/px/tiff/libtiff/.libs"
      "C:/Program Files/GnuWin32/lib"
      "C:/Program Files (x86)/GnuWin32/lib"
  )
  #if we found the library, search for tiff.h in <tiff-lib-dir>/../include
  if(TIFF_LIBRARY)
    get_filename_component(TIFF_LIBRARY_DIR ${TIFF_LIBRARY} PATH)
    set(TIFF_INCLUDE_HINT "${TIFF_LIBRARY_DIR}/../include")
    get_filename_component(TIFF_INCLUDE_HINT "${TIFF_INCLUDE_HINT}" ABSOLUTE)
    find_file(TIFF_HEADER tiff.h PATHS ${TIFF_INCLUDE_HINT})
    mark_as_advanced(TIFF_HEADER)
    if(TIFF_HEADER)
      get_filename_component(TIFF_INCLUDE_DIR "${TIFF_HEADER}" PATH)
      set(TIFF_INCLUDE_DIR "${TIFF_INCLUDE_DIR}" CACHE PATH "Directory containing tiff.h" FORCE)
    endif()
  endif()

  find_library(Z_LIBRARY z)
  if(TIFF_LIBRARY AND Z_LIBRARY)
    set(TIFF_FOUND TRUE)
    get_filename_component(Z_LIBRARY_RESOLVED ${Z_LIBRARY} REALPATH)
    set(TIFF_LIBRARIES ${TIFF_LIBRARY} ${Z_LIBRARY_RESOLVED})
  else()
    #try find_package again
    #it should succeed now if TIFF_LIBRARY and TIFF_INCLUDE_DIR are set
    find_package(TIFF)
    mark_as_advanced(CLEAR TIFF_INCLUDE_DIR TIFF_LIBRARY)
  endif()
endif()

if(TIFF_FOUND)
  add_executable(ra_tiff ra_tiff.c)
  target_link_libraries(ra_tiff rtrad ${TIFF_LIBRARIES})
  add_executable(normtiff normtiff.c)
  target_link_libraries(normtiff rtrad ${TIFF_LIBRARIES})
  if(BUILD_LIBTIFF)
    add_dependencies(ra_tiff libtiff)
    add_dependencies(normtiff libtiff)
  endif()
  install(TARGETS ra_tiff normtiff
    RUNTIME DESTINATION "bin"
  )
  # install(FILES ${TIFF_LIBRARIES} DESTINATION "bin")
  install(PROGRAMS ran2tiff.csh
    DESTINATION "bin"
    RENAME ran2tiff
  )
#TODO: else case
endif()

install(TARGETS
  macbethcal
  pcomb
  pcompos
  pcond
  pcwarp
  pextrem
  pfilt
  pflip
  pinterp
  pmblur2
  protate
  psign
  psketch
  pvalue
  ra_bmp
  ra_gif
  ra_hexbit
  ra_pict
  ra_ppm
  ra_ps
  ra_rgbe
  ra_t16
  ra_t8
  ra_xyze
  ttyimage
  RUNTIME DESTINATION "bin"
)

if(WIN32)
  install(PROGRAMS falsecolor.pl
    DESTINATION "bin"
  )
  file(COPY falsecolor.pl
    DESTINATION "${CMAKE_BINARY_DIR}/bin/Release"
  )
else()
  install_and_copy_without_extension(
    falsecolor.pl
    normpat.csh
    pbilat.csh
    pdfblur.csh
    pgblur.csh
    phisto.csh
    pmblur.csh
    pmdblur.csh
    ra_pfm.csh
    xyzimage.csh
  )
endif()
