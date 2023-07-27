from PySide6.QtWidgets import  QWidget, QSizePolicy, QGraphicsView, QGraphicsScene
from PySide6.QtCore import QSize, QRect , Qt, QPointF, QPoint
from PySide6.QtGui import QPixmap, QIcon, QImage, QPaintEvent, QColor, \
    QFont, QPainter, QBrush, QFont, QPen, QDragMoveEvent, QDragLeaveEvent, \
    QMouseEvent, QTransform

class DrawingCanvas(QGraphicsScene):

    def __init__(self,
                 drawing_colour : QColor=None,
                 font_type: QFont=None,
                 font_size: int=None) -> None:
        super(DrawingCanvas, self).__init__()

        self.drawing_color = drawing_colour
        self.drawing_font = font_type
        self.font_size = font_size

        self.black_brush = QBrush(Qt.black)
        self.drawing_pen = QPen(Qt.black)

        self.draw_ellipses()

    def draw_text(self):
        """
            draws text to the active viewport
        """
        pass

    def draw_ellipses(self):
        """
            draws particle to screen
        """
        num_particles_to_draw_x = 100
        num_particles_to_draw_y = 100

        spacing = 1

        for i in range(0, num_particles_to_draw_x, spacing):
            for j in range(0, num_particles_to_draw_y, spacing):
                self.addEllipse(i, j, i + 1, j + 1, self.drawing_pen,
                                self.black_brush)