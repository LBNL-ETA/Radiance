
include(${QT_USE_FILE})

qt4_wrap_cpp(rvu_moc rvuwidget.h mainwindow.h)
qt4_wrap_ui(rvu_ui mainwindow.ui)
qt4_wrap_ui(exposure_dialog_ui exposuredialog.ui)
qt4_wrap_ui(parameter_dialog_ui parameterdialog.ui)
qt4_wrap_ui(view_dialog_ui viewdialog.ui)
qt4_wrap_ui(increments_dialog_ui incrementsdialog.ui)
qt4_wrap_ui(commands_dialog_ui commandsdialog.ui)

include_directories(${radiance_SOURCE_DIR}/src/rt
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_library(qtrvu  qt.c rvuwidget.cxx qt_rvu_main.cxx mainwindow.cxx
  ${rvu_moc} ${rvu_ui} ${exposure_dialog_ui} ${parameter_dialog_ui}
  ${view_dialog_ui} ${increments_dialog_ui} ${commands_dialog_ui})
target_link_libraries(qtrvu ${QT_LIBRARIES})

install(TARGETS qtrvu
        LIBRARY DESTINATION bin
        ARCHIVE DESTINATION bin)
