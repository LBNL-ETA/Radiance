#include <QtGui/QApplication>
#include <QtGui/QColor>

extern "C" int qt_rvu_run();
extern "C" int qt_rvu_init(char* name, char* id,
                           int* xsize, int* ysize);

int main(int argc, char* argv[])
{
  int x, y;
  qt_rvu_init("test", "id", &x, &y);
  return qt_rvu_run();
}
