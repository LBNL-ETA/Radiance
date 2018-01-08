qt5_wrap_cpp(rvu_moc rvuwidget.h mainwindow.h)
qt5_wrap_ui(rvu_ui mainwindow.ui)
qt5_wrap_ui(exposure_dialog_ui exposuredialog.ui)
qt5_wrap_ui(parameter_dialog_ui parameterdialog.ui)
qt5_wrap_ui(view_dialog_ui viewdialog.ui)
qt5_wrap_ui(increments_dialog_ui incrementsdialog.ui)
qt5_wrap_ui(commands_dialog_ui commandsdialog.ui)

include_directories(
  ${CMAKE_SOURCE_DIR}/src/rt
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

add_library(qtrvu
  mainwindow.cxx
  qt.c
  qt_rvu_main.cxx
  rvuwidget.cxx
  ${commands_dialog_ui}
  ${exposure_dialog_ui}
  ${increments_dialog_ui}
  ${parameter_dialog_ui}
  ${rvu_moc}
  ${rvu_ui}
  ${view_dialog_ui}
)
target_link_libraries(qtrvu Qt5::Widgets)
