#include <QtGui/QApplication>
#include <QtGui/QInputDialog>
#include <QtGui/QColor>
#include <QtCore/QDir>
#include <string>
#include <iostream>
#include "rvuwidget.h"
#include "mainwindow.h"

// file level globals
static RvuWidget* WidgetInstance = 0;
static MainWindow* MainWindowInstance = 0;
static QApplication* ApplicationInstance = 0;

// default size of window
static int XSize = 1024;
static int YSize = 768;

extern "C"
int qt_rvu_paint(int r,int g,int b,
                 int xmin, int ymin, int xmax, int ymax
  )
{
  WidgetInstance->drawRect(xmin, YSize - ymax -1,
                           abs(xmax-xmin),
                           abs(ymax-ymin),
                           QColor(r, g, b));
  return 0;
}

extern "C" void qt_set_progress(int p)
{
  MainWindowInstance->setProgress(p);
}
extern "C" void qt_flush_display()
{
  WidgetInstance->update();
  ApplicationInstance->processEvents();
}

extern "C" int qt_rvu_run()
{
  MainWindowInstance->show();
  int ret = ApplicationInstance->exec();
  return ret;
}

extern "C" void qt_rvu_get_position(int *x, int *y)
{
  MainWindowInstance->pick(x, y);
  *y = YSize - *y -1;
}

extern "C"
int qt_rvu_init(char* name, char* id,
                int* xsize, int* ysize)
{
  static int argc = 1;
  static char * argv[] = {"rvu"};
  ApplicationInstance = new QApplication(argc, argv);
  *xsize = XSize;
  *ysize = YSize;
  MainWindowInstance = new MainWindow(*xsize, *ysize, name, id);
  WidgetInstance = MainWindowInstance->getRvuWidget();
  return 0;
}

extern "C"
void qt_window_comout(const char* m)
{
  MainWindowInstance->setStatus(m);
}

extern "C"
void qt_resize(int X, int Y)
{
  XSize = X;
  YSize = Y;
  MainWindowInstance->resizeImage(X, Y);
  MainWindowInstance->show();
}

// Return 1 if inp was changed
// Return 0 if no change was made to inp
extern "C"
int qt_open_text_dialog(char* inp, const char* prompt)
{
  bool ok;
  QString text = QInputDialog::getText(MainWindowInstance,
                                       "QInputDialog::getText()",
                                       prompt, QLineEdit::Normal,
                                       "", &ok);
  if (ok && !text.isEmpty())
    {
    if(text != inp)
      {
      strcpy(inp, text.toAscii());
      return 1;
      }
    }
  return 0;
}
