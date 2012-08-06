#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "ui_exposuredialog.h"
#include "ui_parameterdialog.h"
#include "ui_viewdialog.h"
#include "ui_incrementsdialog.h"
#include "ui_commandsdialog.h"

#include <QtGui/QMessageBox>
#include <QtGui/QLineEdit>
#include <QtGui/QCloseEvent>
#include <QtGui/QFileDialog>
#include <QtGui/QInputDialog>
#include <QtCore/QDebug>
#include <QtCore/QTime>
#include <QtCore/QTextStream>

#include <iostream>

// Process a command
extern "C" void qt_process_command(const char*);
// set the abort flag to stop a render in progress
extern "C" void qt_set_abort(int );
extern "C"
{
//externs for exposure dialog
extern VIEW ourview;
extern double exposure;
//externs for parameters dialog
extern COLOR  ambval;   /* ambient value */
extern int  ambvwt;   /* initial weight for ambient value */
extern double ambacc;   /* ambient accuracy */
extern int  ambres;   /* ambient resolution */
extern int  ambdiv;   /* ambient divisions */
extern int  ambssamp; /* ambient super-samples */
extern int  ambounce; /* ambient bounces */
extern double shadcert; /* shadow testing certainty */
extern double shadthresh; /* shadow threshold */
extern int  directvis;  /* light sources visible to eye? */
extern COLOR  cextinction;  /* global extinction coefficient */
extern COLOR  salbedo;  /* global scattering albedo */
extern double seccg;    /* global scattering eccentricity */
extern double ssampdist;  /* scatter sampling distance */
extern double specthresh; /* specular sampling threshold */
extern double specjitter; /* specular sampling jitter */
extern int  maxdepth; /* maximum recursion depth */
extern double minweight;  /* minimum ray weight */
extern double dstrsrc;  /* square source distribution */
extern double srcsizerat; /* maximum source size/dist. ratio */
extern int  psample;    /* pixel sample size */
extern double  maxdiff;   /* max. sample difference */
void quit(int code);
}
MainWindow::MainWindow(int width, int height, const char* name, const char* id)
{
  m_ui = new Ui::MainWindow;
  m_ui->setupUi(this);
  m_exposureDialog = new QDialog();
  m_exposureDialogUi = new Ui::exposureDialog;
  m_exposureDialogUi->setupUi(m_exposureDialog);
  m_parameterDialog = new QDialog();
  m_parameterDialogUi = new Ui::parameterDialog;
  m_parameterDialogUi->setupUi(m_parameterDialog);
  m_viewDialog = new QDialog();
  m_viewDialogUi = new Ui::viewDialog;
  m_viewDialogUi->setupUi(m_viewDialog);
  this->smallIncrement = 0.1;
  this->largeIncrement = 0.5;
  m_incrementsDialog = new QDialog();
  m_incrementsDialogUi = new Ui::incrementsDialog;
  m_incrementsDialogUi->setupUi(m_incrementsDialog);
  m_commandsDialog = new QDialog();
  m_commandsDialogUi = new Ui::commandsDialog;
  m_commandsDialogUi->setupUi(m_commandsDialog);
  this->currentImageName  = "";
  createActions();
  connectSlots();
  m_ui->progressBar->setMinimum(0);
  m_ui->progressBar->setMaximum(100);
  setWindowTitle(tr(id));
  resize(width, height);
  updatePositionLabels();
}

MainWindow::~MainWindow()
{
  delete m_ui;
  delete m_exposureDialog;
  delete m_parameterDialog;
  delete m_viewDialog;
  delete m_incrementsDialog;
  delete m_commandsDialog;
}

void MainWindow::setProgress(int p)
{
  m_ui->progressBar->setValue(p);
}

RvuWidget* MainWindow::getRvuWidget() const
{
  return m_ui->rvuWidget;
}

void MainWindow::closeEvent(QCloseEvent *event)
{
  // Here is a good spot to shut down the rest of the code before exit.
  this->m_exposureDialog->close();
  this->m_parameterDialog->close();
  this->m_viewDialog->close();
  this->m_incrementsDialog->close();
  this->m_commandsDialog->close();
}

void MainWindow::createActions()
{
  m_minimizeAction = new QAction(tr("Mi&nimize"), this);
  connect(m_minimizeAction, SIGNAL(triggered()), this, SLOT(hide()));

  m_maximizeAction = new QAction(tr("Ma&ximize"), this);
  connect(m_maximizeAction, SIGNAL(triggered()), this, SLOT(showMaximized()));

  m_restoreAction = new QAction(tr("&Restore"), this);
  connect(m_restoreAction, SIGNAL(triggered()), this, SLOT(showNormal()));
}

void MainWindow::connectSlots()
{
  connect(m_ui->exit, SIGNAL(triggered()), qApp, SLOT(quit()));
  connect(m_ui->pushButton, SIGNAL(clicked()), this, SLOT(buttonPressed()));
  connect(m_ui->lineEdit, SIGNAL(returnPressed()),
          this, SLOT(buttonPressed()));
  connect(m_ui->redraw, SIGNAL(triggered()), this, SLOT(redraw()));
  connect(m_ui->trace, SIGNAL(triggered()), this, SLOT(traceRay()));
  connect(m_ui->toggleToolBar, SIGNAL(triggered()),
          this, SLOT(toggleToolBar()));
  connect(m_ui->toggleStatusBar, SIGNAL(triggered()),
          this, SLOT(toggleStatusBar()));
  connect(m_ui->toggleTranslateTool, SIGNAL(triggered()),
          this, SLOT(toggleTranslateTool()));
  connect(m_ui->exposure, SIGNAL(triggered()),
          this, SLOT(showExposureDialog()));
  connect(m_ui->parameters, SIGNAL(triggered()),
          this, SLOT(showParameterDialog()));
  connect(m_ui->view, SIGNAL(triggered()),
          this, SLOT(showViewDialog()));
  connect(m_ui->saveImage, SIGNAL(triggered()),
          this, SLOT(saveCurrentImage()));
  connect(m_ui->saveImageAs, SIGNAL(triggered()),
          this, SLOT(saveImage()));
  connect(m_ui->loadView, SIGNAL(triggered()),
          this, SLOT(loadView()));
  connect(m_ui->loadRif, SIGNAL(triggered()),
          this, SLOT(loadRif()));
  connect(m_ui->backfaceVisibility, SIGNAL(triggered()),
          this, SLOT(toggleBackfaceVisibility()));
  connect(m_ui->grayscale, SIGNAL(triggered()),
          this, SLOT(toggleGrayscale()));
  connect(m_ui->irradiance, SIGNAL(triggered()),
          this, SLOT(toggleIrradiance()));
  connect(m_ui->appendToRif, SIGNAL(triggered()),
          this, SLOT(appendToRif()));
  connect(m_ui->appendToView, SIGNAL(triggered()),
          this, SLOT(appendToView()));
  connect(m_ui->GUI_Increments, SIGNAL(triggered()),
          this, SLOT(showIncrementsDialog()));
  connect(m_ui->commands, SIGNAL(triggered()),
          this, SLOT(showCommandsDialog()));
  connect(m_ui->aboutRvu, SIGNAL(triggered()),
          this, SLOT(showAbout()));

  //movement buttons
  connect(m_ui->xPlus1Button, SIGNAL(clicked()),
          this, SLOT(moveXPlusOne()));
  connect(m_ui->xPlus2Button, SIGNAL(clicked()),
          this, SLOT(moveXPlusTwo()));
  connect(m_ui->xMinus1Button, SIGNAL(clicked()),
          this, SLOT(moveXMinusOne()));
  connect(m_ui->xMinus2Button, SIGNAL(clicked()),
          this, SLOT(moveXMinusTwo()));
  connect(m_ui->yPlus1Button, SIGNAL(clicked()),
          this, SLOT(moveYPlusOne()));
  connect(m_ui->yPlus2Button, SIGNAL(clicked()),
          this, SLOT(moveYPlusTwo()));
  connect(m_ui->yMinus1Button, SIGNAL(clicked()),
          this, SLOT(moveYMinusOne()));
  connect(m_ui->yMinus2Button, SIGNAL(clicked()),
          this, SLOT(moveYMinusTwo()));
  connect(m_ui->zPlus1Button, SIGNAL(clicked()),
          this, SLOT(moveZPlusOne()));
  connect(m_ui->zPlus2Button, SIGNAL(clicked()),
          this, SLOT(moveZPlusTwo()));
  connect(m_ui->zMinus1Button, SIGNAL(clicked()),
          this, SLOT(moveZMinusOne()));
  connect(m_ui->zMinus2Button, SIGNAL(clicked()),
          this, SLOT(moveZMinusTwo()));

  //exposure dialog
  connect(m_exposureDialogUi->exposureButtonBox, SIGNAL(accepted()),
          this, SLOT(adjustExposure()));
  connect(m_exposureDialogUi->exposureButtonBox, SIGNAL(rejected()),
          this, SLOT(hideExposureDialog()));
  connect(m_exposureDialogUi->natural, SIGNAL(toggled(bool)),
          this, SLOT(updatePointRadio()));
  connect(m_exposureDialogUi->relative, SIGNAL(toggled(bool)),
          this, SLOT(updatePointRadio()));

  //parameter dialog
  connect(m_parameterDialogUi->buttonBox, SIGNAL(accepted()),
          this, SLOT(adjustParameters()));
  connect(m_parameterDialogUi->buttonBox, SIGNAL(rejected()),
          this, SLOT(hideParameterDialog()));

  //view dialog
  connect(m_viewDialogUi->buttonBox, SIGNAL(accepted()),
          this, SLOT(adjustView()));
  connect(m_viewDialogUi->buttonBox, SIGNAL(rejected()),
          this, SLOT(hideViewDialog()));

  //increments dialog
  connect(m_incrementsDialogUi->buttonBox, SIGNAL(accepted()),
          this, SLOT(adjustIncrements()));
  connect(m_incrementsDialogUi->buttonBox, SIGNAL(rejected()),
          this, SLOT(hideIncrementsDialog()));
}

void MainWindow::doSubmit()
{
  QString command = m_ui->lineEdit->text();
  bool centerCrossHairs = false;
  if(command == "aim")
    {
    centerCrossHairs = true;
    }
  m_ui->lineEdit->setText("");
  m_ui->progressBar->setEnabled(true);
  this->setProgress(0);
  QTime t;
  t.start();
  if (command.isEmpty())
    {
    this->m_ui->pushButton->setText("Submit");
    enableInterface(true);
    return;
    }
  qt_process_command(command.toAscii());
  QString msg;
  QTextStream(&msg) << "Render Time: " << t.elapsed() << " ms";
  m_ui->messageBox->appendPlainText(msg);
  this->m_ui->pushButton->setText("Submit");
  enableInterface(true);
  if(centerCrossHairs)
    {
    QPoint p = this->m_ui->rvuWidget->geometry().center();
    this->m_ui->rvuWidget->setPosition(p.x(), p.y());
    }
  updatePositionLabels();
}


void MainWindow::buttonPressed()
{
  if( m_ui->pushButton->text() == QString("Submit") )
    {
    qt_set_abort(0);
    this->m_ui->pushButton->setText("Abort");
    enableInterface(false);
    this->doSubmit();
    }
  else if( m_ui->pushButton->text() == QString("Abort") )
    {
    qt_set_abort(1);
    this->m_ui->pushButton->setText("Submit");
    enableInterface(true);
    }
}

void MainWindow::resizeImage(int width, int height)
{
  int pad = m_ui->lineEdit->size().height();
  pad += m_ui->menubar->size().height();
  pad += m_ui->messageBox->size().height();
  this->resize(width+4, height+pad+4);
  this->m_ui->rvuWidget->resizeImage(width, height);
}

void MainWindow::setStatus(const char* m)
{
  m_ui->messageBox->appendPlainText(QString(m).trimmed());
}

void MainWindow::redraw()
{
  runCommand("redraw");
}

void MainWindow::traceRay()
{
  runCommand("t");
}

void MainWindow::toggleToolBar()
{
  m_ui->toolBar->setVisible(m_ui->toggleToolBar->isChecked());
}

void MainWindow::toggleStatusBar()
{
  m_ui->messageBox->setVisible(m_ui->toggleStatusBar->isChecked());
}

void MainWindow::toggleTranslateTool()
{
  m_ui->xWidget->setVisible(m_ui->toggleTranslateTool->isChecked());
  m_ui->yWidget->setVisible(m_ui->toggleTranslateTool->isChecked());
  m_ui->zWidget->setVisible(m_ui->toggleTranslateTool->isChecked());
}

void MainWindow::refreshView(VIEW *nv)
{
	if (waitrays() < 0)			/* clear ray queue */
		quit(1);
  newview(nv);
  this->redraw();
}

void MainWindow::move(int direction, float amount)
{
  if( m_ui->pushButton->text() == QString("Abort") )
    {
    qt_set_abort(1);
    return;
    } 
  VIEW nv = ourview;
  nv.vp[direction] = nv.vp[direction] + amount;
  this->refreshView(&nv);
  updatePositionLabels();
}

void MainWindow::moveXPlusOne()
{
  this->move(0, this->smallIncrement);
}

void MainWindow::moveXPlusTwo()
{
  this->move(0, this->largeIncrement);
}

void MainWindow::moveXMinusOne()
{
  this->move(0, -1 * this->smallIncrement);
}

void MainWindow::moveXMinusTwo()
{
  this->move(0, -1 * this->largeIncrement);
}

void MainWindow::moveYPlusOne()
{
  this->move(1, this->smallIncrement);
}

void MainWindow::moveYPlusTwo()
{
  this->move(1, this->largeIncrement);
}

void MainWindow::moveYMinusOne()
{
  this->move(1, -1 * this->smallIncrement);
}

void MainWindow::moveYMinusTwo()
{
  this->move(1, -1 * this->largeIncrement);
}

void MainWindow::moveZPlusOne()
{
  this->move(2, this->smallIncrement);
}

void MainWindow::moveZPlusTwo()
{
  this->move(2, this->largeIncrement);
}

void MainWindow::moveZMinusOne()
{
  this->move(2, -1 * this->smallIncrement);
}

void MainWindow::moveZMinusTwo()
{
  this->move(2, -1 * this->largeIncrement);
}

void MainWindow::updatePositionLabels()
{
  m_ui->xValue->setText(QString::number(ourview.vp[0]));
  m_ui->yValue->setText(QString::number(ourview.vp[1]));
  m_ui->zValue->setText(QString::number(ourview.vp[2]));
}

void MainWindow::showExposureDialog()
{
  m_exposureDialogUi->exposureValue->setText(QString::number(exposure));
  m_exposureDialog->show();
}

void MainWindow::hideExposureDialog()
{
  m_exposureDialog->hide();
}

void MainWindow::showParameterDialog()
{
  //ambient tab
  m_parameterDialogUi->valueX->setText(QString::number(ambval[0]));
  m_parameterDialogUi->valueY->setText(QString::number(ambval[1]));
  m_parameterDialogUi->valueZ->setText(QString::number(ambval[2]));
  m_parameterDialogUi->weight->setText(QString::number(ambvwt));
  m_parameterDialogUi->accuracy->setText(QString::number(ambacc));
  m_parameterDialogUi->divisions->setText(QString::number(ambdiv));
  m_parameterDialogUi->supersamples->setText(QString::number(ambssamp));
  m_parameterDialogUi->bounces->setText(QString::number(ambounce));
  m_parameterDialogUi->resolution->setText(QString::number(ambres));

  //direct tab
  m_parameterDialogUi->certainty->setText(QString::number(shadcert));
  m_parameterDialogUi->threshold->setText(QString::number(shadthresh));
  m_parameterDialogUi->visibility->setText(QString::number(directvis));
  m_parameterDialogUi->jitter->setText(QString::number(dstrsrc));
  m_parameterDialogUi->sampling->setText(QString::number(srcsizerat));
  //values on direct tab that I can't find: reflection

  //limit tab
  m_parameterDialogUi->limit_weight->setText(QString::number(minweight));
  m_parameterDialogUi->reflection->setText(QString::number(maxdepth));

  //medium tab
  m_parameterDialogUi->albedoX->setText(QString::number(salbedo[0]));
  m_parameterDialogUi->albedoY->setText(QString::number(salbedo[1]));
  m_parameterDialogUi->albedoZ->setText(QString::number(salbedo[2]));
  m_parameterDialogUi->extinctionX->setText(QString::number(cextinction[0]));
  m_parameterDialogUi->extinctionY->setText(QString::number(cextinction[1]));
  m_parameterDialogUi->extinctionZ->setText(QString::number(cextinction[2]));
  m_parameterDialogUi->eccentricity->setText(QString::number(seccg));
  m_parameterDialogUi->mistDistance->setText(QString::number(ssampdist));

  //pixel tab
  m_parameterDialogUi->sample->setText(QString::number(psample));
  m_parameterDialogUi->px_threshold->setText(QString::number(maxdiff));

  //specular tab
  m_parameterDialogUi->sp_jitter->setText(QString::number(specjitter));
  m_parameterDialogUi->sp_threshold->setText(QString::number(specthresh));

  m_parameterDialog->show();
}

void MainWindow::hideParameterDialog()
{
  m_parameterDialog->hide();
}

void MainWindow::adjustExposure()
{
  bool checkForPoint = false;

  QString command = "exposure ";

  if(m_exposureDialogUi->absolute->isChecked())
    {
    command += "=";
    }
  else if(m_exposureDialogUi->fstop->isChecked())
    {
    //TODO: check if the value is positive or not here
    command += "+";
    }
  else if(m_exposureDialogUi->natural->isChecked())
    {
    command += "@";
    checkForPoint = true;
    }
  else if(m_exposureDialogUi->relative->isChecked())
    {
    checkForPoint = true;
    }

  if(!(checkForPoint && m_exposureDialogUi->point->isChecked()))
    {
    command += m_exposureDialogUi->exposureSetting->text();
    }

  runCommand(command.toAscii());
}

void MainWindow::updatePointRadio()
{
  if(m_exposureDialogUi->relative->isChecked() ||
     m_exposureDialogUi->natural->isChecked())
    {
    m_exposureDialogUi->point->setEnabled(true);
    m_exposureDialogUi->average->setEnabled(true);
    }
  else
    {
    m_exposureDialogUi->point->setEnabled(false);
    m_exposureDialogUi->average->setEnabled(false);
    }
}

void MainWindow::adjustParameters()
{
  if(m_parameterDialogUi->valueX->isModified() ||
     m_parameterDialogUi->valueY->isModified() ||
     m_parameterDialogUi->valueZ->isModified())
    {
    QString command = "set av ";
    command += m_parameterDialogUi->valueX->text();
    command += " ";
    command += m_parameterDialogUi->valueY->text();
    command += " ";
    command += m_parameterDialogUi->valueZ->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->weight->isModified())
    {
    QString command = "set aw ";
    command += m_parameterDialogUi->weight->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->accuracy->isModified())
    {
    QString command = "set aa ";
    command += m_parameterDialogUi->accuracy->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->divisions->isModified())
    {
    QString command = "set ad ";
    command += m_parameterDialogUi->divisions->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->supersamples->isModified())
    {
    QString command = "set as ";
    command += m_parameterDialogUi->supersamples->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->bounces->isModified())
    {
    QString command = "set ab ";
    command += m_parameterDialogUi->bounces->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->resolution->isModified())
    {
    QString command = "set ar ";
    command += m_parameterDialogUi->resolution->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->certainty->isModified())
    {
    QString command = "set dc ";
    command += m_parameterDialogUi->certainty->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->threshold->isModified())
    {
    QString command = "set dt ";
    command += m_parameterDialogUi->threshold->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->visibility->isModified())
    {
    QString command = "set dv ";
    command += m_parameterDialogUi->visibility->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->jitter->isModified())
    {
    QString command = "set dj ";
    command += m_parameterDialogUi->jitter->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->sampling->isModified())
    {
    QString command = "set ds ";
    command += m_parameterDialogUi->sampling->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->limit_weight->isModified())
    {
    QString command = "set lw ";
    command += m_parameterDialogUi->limit_weight->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->reflection->isModified())
    {
    QString command = "set lr ";
    command += m_parameterDialogUi->reflection->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->albedoX->isModified() ||
     m_parameterDialogUi->albedoY->isModified() ||
     m_parameterDialogUi->albedoZ->isModified())
    {
    QString command = "set ma ";
    command += m_parameterDialogUi->albedoX->text();
    command += " ";
    command += m_parameterDialogUi->albedoY->text();
    command += " ";
    command += m_parameterDialogUi->albedoZ->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->extinctionX->isModified() ||
     m_parameterDialogUi->extinctionY->isModified() ||
     m_parameterDialogUi->extinctionZ->isModified())
    {
    QString command = "set me ";
    command += m_parameterDialogUi->extinctionX->text();
    command += " ";
    command += m_parameterDialogUi->extinctionY->text();
    command += " ";
    command += m_parameterDialogUi->extinctionZ->text();
    }
  if(m_parameterDialogUi->eccentricity->isModified())
    {
    QString command = "set mg ";
    command += m_parameterDialogUi->eccentricity->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->mistDistance->isModified())
    {
    QString command = "set ms ";
    command += m_parameterDialogUi->mistDistance->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->sample->isModified())
    {
    QString command = "set ps ";
    command += m_parameterDialogUi->sample->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->px_threshold->isModified())
    {
    QString command = "set pt ";
    command += m_parameterDialogUi->px_threshold->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->sp_jitter->isModified())
    {
    QString command = "set ss ";
    command += m_parameterDialogUi->sp_jitter->text();
    runCommand(command.toAscii());
    }
  if(m_parameterDialogUi->sp_threshold->isModified())
    {
    QString command = "set st ";
    command += m_parameterDialogUi->sp_threshold->text();
    runCommand(command.toAscii());
    }
  m_parameterDialog->hide();
}

void MainWindow::showViewDialog()
{
  if(ourview.type==VT_PER)
    {
    m_viewDialogUi->perspective->setChecked(true);
    }
  else if(ourview.type==VT_PAR)
    {
    m_viewDialogUi->parallel->setChecked(true);
    }
  else if(ourview.type==VT_HEM)
    {
    m_viewDialogUi->hemispherical->setChecked(true);
    }
  else if(ourview.type==VT_ANG)
    {
    m_viewDialogUi->angular->setChecked(true);
    }
  else if(ourview.type==VT_PLS)
    {
    m_viewDialogUi->planispheric->setChecked(true);
    }
  else if(ourview.type==VT_CYL)
    {
    m_viewDialogUi->cylindrical->setChecked(true);
    }

  m_viewDialogUi->viewX->setText(QString::number(ourview.vp[0]));
  m_viewDialogUi->viewY->setText(QString::number(ourview.vp[1]));
  m_viewDialogUi->viewZ->setText(QString::number(ourview.vp[2]));

  m_viewDialogUi->dirX->setText(QString::number(ourview.vdir[0], 'f', 4));
  m_viewDialogUi->dirY->setText(QString::number(ourview.vdir[1], 'f', 4));
  m_viewDialogUi->dirZ->setText(QString::number(ourview.vdir[2], 'f', 4));

  m_viewDialogUi->upX->setText(QString::number(ourview.vup[0]));
  m_viewDialogUi->upY->setText(QString::number(ourview.vup[1]));
  m_viewDialogUi->upZ->setText(QString::number(ourview.vup[2]));

  m_viewDialogUi->sizeHoriz->setText(QString::number(ourview.horiz));
  m_viewDialogUi->sizeVert->setText(QString::number(ourview.vert));
  m_viewDialogUi->fore->setText(QString::number(ourview.vfore));
  m_viewDialogUi->aft->setText(QString::number(ourview.vaft));
  m_viewDialogUi->offsetHoriz->setText(QString::number(ourview.hoff));
  m_viewDialogUi->offsetVert->setText(QString::number(ourview.voff));

  m_viewDialog->show();
}

void MainWindow::hideViewDialog()
{
  m_viewDialog->hide();
}

void MainWindow::adjustView()
{
  bool changed = false;
  VIEW nv = ourview;

  //type
  if(m_viewDialogUi->perspective->isChecked() && nv.type != VT_PER)
    {
    nv.type = VT_PER;
    changed = true;
    }
  if(m_viewDialogUi->parallel->isChecked() && nv.type != VT_PAR)
    {
    nv.type = VT_PAR;
    changed = true;
    }
  if(m_viewDialogUi->hemispherical->isChecked() && nv.type != VT_HEM)
    {
    nv.type = VT_HEM;
    changed = true;
    }
  if(m_viewDialogUi->angular->isChecked() && nv.type != VT_ANG)
    {
    nv.type = VT_ANG;
    changed = true;
    }
  if(m_viewDialogUi->planispheric->isChecked() && nv.type != VT_PLS)
    {
    nv.type = VT_PLS;
   changed = true;
    }
  if(m_viewDialogUi->cylindrical->isChecked() && nv.type != VT_CYL)
    {
    nv.type = VT_CYL;
    changed = true;
    }

  //viewpoint
  if(m_viewDialogUi->viewX->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->viewX->text().toFloat(&okay);
    if(okay)
      {
      nv.vp[0] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->viewY->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->viewY->text().toFloat(&okay);
    if(okay)
      {
      nv.vp[1] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->viewZ->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->viewZ->text().toFloat(&okay);
    if(okay)
      {
      nv.vp[2] = f;
      changed = true;
      }
    }

  //direction
  if(m_viewDialogUi->dirX->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->dirX->text().toFloat(&okay);
    if(okay)
      {
      nv.vdir[0] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->dirY->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->dirY->text().toFloat(&okay);
    if(okay)
      {
      nv.vdir[1] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->dirZ->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->dirZ->text().toFloat(&okay);
    if(okay)
      {
      nv.vdir[2] = f;
      changed = true;
      }
    }

  //up
  if(m_viewDialogUi->upX->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->upX->text().toFloat(&okay);
    if(okay)
      {
      nv.vup[0] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->upY->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->upY->text().toFloat(&okay);
    if(okay)
      {
      nv.vup[1] = f;
      changed = true;
      }
    }
  if(m_viewDialogUi->upZ->isModified())
    {
    bool okay = true;
    float f = m_viewDialogUi->upZ->text().toFloat(&okay);
    if(okay)
      {
      nv.vup[2] = f;
      changed = true;
      }
    }

  //size
  if(m_viewDialogUi->sizeHoriz->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->sizeHoriz->text().toDouble(&okay);
    if(okay)
      {
      nv.horiz = d;
      changed = true;
      }
    }
  if(m_viewDialogUi->sizeVert->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->sizeVert->text().toDouble(&okay);
    if(okay)
      {
      nv.vert = d;
      changed = true;
      }
    }

  //clipping plane
  if(m_viewDialogUi->fore->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->fore->text().toDouble(&okay);
    if(okay)
      {
      nv.vfore = d;
      changed = true;
      }
    }
  if(m_viewDialogUi->aft->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->aft->text().toDouble(&okay);
    if(okay)
      {
      nv.vaft = d;
      changed = true;
      }
    }

  //clipping plane
  if(m_viewDialogUi->offsetHoriz->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->offsetHoriz->text().toDouble(&okay);
    if(okay)
      {
      nv.hoff = d;
      changed = true;
      }
    }
  if(m_viewDialogUi->offsetVert->isModified())
    {
    bool okay = true;
    double d = m_viewDialogUi->offsetVert->text().toDouble(&okay);
    if(okay)
      {
      nv.voff = d;
      changed = true;
      }
    }

  if(changed)
    {
    this->refreshView(&nv);
    }
}

void MainWindow::enableInterface(bool enable)
{
  m_ui->toolBar->setEnabled(enable);
  m_ui->progressBar->setEnabled(!enable);
}

void MainWindow::runCommand(const char *command)
{
  qt_set_abort(0);
  this->m_ui->pushButton->setText("Abort");
  enableInterface(false);
  qt_process_command(command);
  qt_set_abort(1);
  this->m_ui->pushButton->setText("Submit");
  enableInterface(true);
}

void MainWindow::pick(int* x, int* y)
{
  QCursor oldCursor = this->cursor();
  this->setCursor(Qt::CrossCursor);
  this->getRvuWidget()->getPosition(x, y);
  // restore state
  this->setCursor(oldCursor);
}
  
void MainWindow::saveImage()
{
  this->currentImageName = QFileDialog::getSaveFileName(this, tr("Save File"),
    "", tr(".hdr images (*.hdr)"));
  QString command = "write " + this->currentImageName;
  this->runCommand(command.toAscii());
}

void MainWindow::saveCurrentImage()
{
  if(this->currentImageName == "")
    {
	this->saveImage();
	return;
	}
  QString command = "write " + this->currentImageName;
  this->runCommand(command.toAscii());
}

void MainWindow::loadView()
{
 QString viewFileName = QFileDialog::getOpenFileName(this, tr("Open View File"),
   "", tr("View Files (*.vp *.vf)"));
 QString command = "last " + viewFileName;
 this->runCommand(command.toAscii());
}
          
void MainWindow::toggleGrayscale()
{
  QString command = "set b ";
  if(m_ui->grayscale->isChecked())
    {
    command += " 1";
    }
  else
    {
    command += " 0";
    }
  this->runCommand(command.toAscii());
  this->runCommand("new");
}

void MainWindow::toggleBackfaceVisibility()
{
  QString command = "set bv";
  if(m_ui->backfaceVisibility->isChecked())
    {
    command += " 1";
    }
  else
    {
    command += " 0";
    }
  this->runCommand(command.toAscii());
  this->runCommand("new");
}

void MainWindow::toggleIrradiance()
{
  QString command = "set i";
  if(m_ui->irradiance->isChecked())
    {
    command += " 1";
    }
  else
    {
    command += " 0";
    }
  this->runCommand(command.toAscii());
  this->runCommand("new");
}
  
void MainWindow::loadRif()
{
  bool ok;
  QString viewName = QInputDialog::getText(this, tr("Input view name"),
                                      tr("Name of view:"), QLineEdit::Normal,
                                      "", &ok);
  if (ok && !viewName.isEmpty())
    {
    QString radFileName = QFileDialog::getOpenFileName(this, tr("Open rad File"),
      "", tr("rad files (*.rif)"));
    QString command = "L " + viewName + " " + radFileName;
    this->runCommand(command.toAscii());
    }
}
  
void MainWindow::appendToRif()
{
  bool ok;
  QString viewName = QInputDialog::getText(this, tr("Input view name"),
                                      tr("Name of view:"), QLineEdit::Normal,
                                      "", &ok);
  if (ok && !viewName.isEmpty())
    {
    QString radFileName = QFileDialog::getSaveFileName(this, tr("Save File"),
      "", tr("rad files (*.rif)"), 0, QFileDialog::DontConfirmOverwrite);
    QString command = "V " + viewName + " " + radFileName;
    this->runCommand(command.toAscii());
    }
}

void MainWindow::appendToView()
{
 QString viewFileName = QFileDialog::getSaveFileName(this, tr("Save File"),
   "", tr("view files (*.vp *.vf)"), 0, QFileDialog::DontConfirmOverwrite);
 if(viewFileName != "")
   {
   QString command = "view " + viewFileName;
   this->runCommand(command.toAscii());
   }
}

void MainWindow::showIncrementsDialog()
{
  m_incrementsDialogUi->small->setText(QString::number(this->smallIncrement));
  m_incrementsDialogUi->large->setText(QString::number(this->largeIncrement));
  m_incrementsDialog->show();
}

void MainWindow::hideIncrementsDialog()
{
  m_incrementsDialog->hide();
}

void MainWindow::adjustIncrements()
{
  bool okay = true;
  float f = m_incrementsDialogUi->small->text().toFloat(&okay);
  if(okay)
    {
    this->smallIncrement = f;
    }
  f = m_incrementsDialogUi->large->text().toFloat(&okay);
  if(okay)
    {
    this->largeIncrement = f;
    }
}

void MainWindow::showCommandsDialog()
{
  m_commandsDialog->show();
}

void MainWindow::showAbout()
{
  QString aboutText =
    "rvu was written mainly by Greg Ward.\nThe Qt-based GUI was developed by Kitware, Inc. in 2011";
  QMessageBox::about(this, "About rvu", aboutText);
}
