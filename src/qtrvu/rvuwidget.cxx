#include "rvuwidget.h"

#include <QtGui/QApplication>
#include <QtGui/QCursor>
#include <QtGui/QColor>
#include <QtGui/QPainter>
#include <QtGui/QPixmap>
#include <QtGui/QMouseEvent>

RvuWidget::RvuWidget(QWidget* parent) : QWidget(parent), m_x(0), m_y(0)
{
  m_image = new QPixmap(200, 200);
  m_image->fill(QColor(0, 0, 0));
  m_painter = new QPainter(m_image);
  m_do_pick = false;
}

RvuWidget::~RvuWidget()
{
  m_painter->end();
  delete m_painter;
  delete m_image;
}

void RvuWidget::resizeImage(int X, int Y)
{
  m_painter->end();
  delete m_image;
  m_image =  new QPixmap(X, Y);
  m_image->fill(QColor(0, 0, 0));
  m_painter->begin(m_image);
  m_painter->setPen(Qt::NoPen);
}

void RvuWidget::drawRect(int x, int y, int width, int height,
                         const QColor &color)
{
  m_painter->fillRect(x, y, width, height, color);
}

void RvuWidget::getPosition(int *x, int *y)
{
  this->m_do_pick = true;
  while(this->m_do_pick)
    {
    QApplication::processEvents();
    }
  *x = m_x;
  *y = m_y;
}

void RvuWidget::paintEvent(QPaintEvent *)
{
  QPainter painter(this);
  // Draw QImage
  painter.drawPixmap(0, 0, *m_image);

  // Draw the crosshairs
  painter.setPen(QColor(0, 255, 0));
  painter.drawLine(m_x - 10, m_y, m_x + 10, m_y);
  painter.drawLine(m_x, m_y - 10, m_x, m_y + 10);
}

void RvuWidget::mousePressEvent(QMouseEvent *event)
{
  if (event->button() == Qt::LeftButton)
    {
    this->m_do_pick = false;
    // Set the cursor position
    m_x = event->x();
    m_y = event->y();
    }
  this->repaint();
}

void RvuWidget::setPosition(int x, int y)
{
  m_x = x;
  m_y = y;
  this->repaint();
}
