
#ifndef __MainWindow_H
#define __MainWindow_H

#include <QtGui/QMainWindow>
#include "fvect.h"
#include "rpaint.h"

class QAction;
class QDialog;
class RvuWidget;

namespace Ui
{
class MainWindow;
class exposureDialog;
class parameterDialog;
class viewDialog;
class incrementsDialog;
class commandsDialog;
}

class MainWindow : public QMainWindow
{
  Q_OBJECT

public:
  MainWindow(int width, int height, const char* name, const char* id);
  ~MainWindow();

  RvuWidget* getRvuWidget() const;
  void resizeImage(int width, int height);
  void setStatus(const char*);
  void setProgress(int);
  void pick( int* x, int* y);
protected:
  void closeEvent(QCloseEvent *event);
  void doSubmit();
  void updatePositionLabels();
  void runCommand(const char *command);
  void move(int direction, float amount);
  void refreshView(VIEW *nv);

protected slots:
  void buttonPressed();
  void redraw();
  void traceRay();
  void saveImage();
  void saveCurrentImage();
  void loadView();
  void toggleToolBar();
  void toggleTranslateTool();
  void toggleStatusBar();
  void moveXPlusOne();
  void moveXPlusTwo();
  void moveYPlusOne();
  void moveYPlusTwo();
  void moveZPlusOne();
  void moveZPlusTwo();
  void moveXMinusOne();
  void moveXMinusTwo();
  void moveYMinusOne();
  void moveYMinusTwo();
  void moveZMinusOne();
  void moveZMinusTwo();
  void showExposureDialog();
  void hideExposureDialog();
  void showParameterDialog();
  void hideParameterDialog();
  void adjustExposure();
  void updatePointRadio();
  void adjustParameters();
  void showViewDialog();
  void hideViewDialog();
  void adjustView();
  void toggleBackfaceVisibility();
  void toggleIrradiance();
  void toggleGrayscale();
  void appendToRif();
  void appendToView();
  void showIncrementsDialog();
  void hideIncrementsDialog();
  void adjustIncrements();
  void showCommandsDialog();
  void showAbout();
  /** Enable/disable elements of the UI while rendering. */
  void enableInterface(bool enable);

private:
  void createActions();
  void connectSlots();

  /** Our MainWindow GUI. */
  Ui::MainWindow *m_ui;
  QAction *m_minimizeAction;
  QAction *m_maximizeAction;
  QAction *m_restoreAction;
  QAction *m_quitAction;

  /** Exposure Dialog */
  QDialog *m_exposureDialog;
  Ui::exposureDialog *m_exposureDialogUi;

  /** Parameter Dialog */
  QDialog *m_parameterDialog;
  Ui::parameterDialog *m_parameterDialogUi;

  /** View Dialog */
  QDialog *m_viewDialog;
  Ui::viewDialog *m_viewDialogUi;
  
  /** Increments Dialog */
  QDialog *m_incrementsDialog;
  Ui::incrementsDialog *m_incrementsDialogUi;
  float smallIncrement;
  float largeIncrement;
  
  /** Commands Dialog */
  QDialog *m_commandsDialog;
  Ui::commandsDialog *m_commandsDialogUi;
  
  /** save image */
  QString currentImageName;
};

#endif
