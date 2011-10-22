add_executable(pfilt pfilt.c pf2.c pf3.c)
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
  add_executable(xshowtrace xshowtrace.c
    ${CMAKE_CURRENT_SOURCE_DIR}/../common/x11findwind.c)
  target_link_libraries(xshowtrace rtrad ${X11_LIBRARIES})
  install(TARGETS ximage xshowtrace RUNTIME DESTINATION bin)
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

find_package(TIFF)
mark_as_advanced(CLEAR TIFF_INCLUDE_DIR TIFF_LIBRARY)
if(TIFF_FOUND)
  add_executable(ra_tiff ra_tiff.c)
  target_link_libraries(ra_tiff rtrad ${TIFF_LIBRARIES})
  install(TARGETS ra_tiff RUNTIME DESTINATION bin)
  install(FILES ${TIFF_LIBRARIES} DESTINATION bin)
#TODO: else case
endif()

install(TARGETS pfilt pflip ttyimage psign pvalue pcompos protate
  ra_hexbit ra_bmp ra_t8 ra_t16 pcomb pinterp ra_ppm ra_rgbe ra_pict ra_ps
  pextrem ra_gif ra_xyze macbethcal pcond pcwarp
  RUNTIME DESTINATION bin)

install(FILES normpat.csh falsecolor.pl pdfblur.csh pmblur.csh pmdblur.csh
  xyzimage.csh pgblur.csh ra_pfm.csh pbilat.csh
  DESTINATION bin)
