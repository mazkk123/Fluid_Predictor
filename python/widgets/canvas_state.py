from PySide6.QtWidgets import  QWidget, QSizePolicy, QGraphicsView, QGraphicsScene
from PySide6.QtCore import QSize, QRect , Qt, QPointF, QPoint
from PySide6.QtGui import QPixmap, QIcon, QImage, QPaintEvent, QColor, \
    QFont, QPainter, QBrush, QFont, QPen, QDragMoveEvent, QDragLeaveEvent, \
    QMouseEvent, QTransform

class DrawingCanvas(QGraphicsView, QGraphicsScene):

    def __init__(self,
                 canvas_size : QRect=None,
                 drawing_colour : QColor=None,
                 font_type: QFont=None,
                 font_size: int=None) -> None:
        super(DrawingCanvas, self).__init__()

        self.main_scene = QGraphicsScene()
        self.main_frame_w = QGraphicsView(self.main_scene)

        self.setGeometry(QRect(0, 0, 1000, 400))
        
        self.drawing_color = drawing_colour
        self.drawing_font = font_type
        self.font_size = font_size

    def paintEvent(self, event: QPaintEvent) -> None:
        """
            controls main painting to the widget window
        """
        self.main_painter = QPainter()

        self.render_points()

        self.main_painter.end()

    def render_points(self, point_radius: int=None, point_colour: QColor=None,
                      pos_x: int=None, pos_y=None, point_spacing: int=1,
                      num_points: int=100) -> None :
        """
            main render points drawing calls
        """
        self.drawing_pen = QPen(QColor(0,0,0,255))
        self.drawing_pen.setWidth(3)
        self.main_painter.setPen(self.drawing_pen)

        # draws a uniform grid scale of spacing specified by parameter
        num_x_points = num_points
        num_y_points = num_points

        for i in range(1,num_x_points,5):
            for j in range(1,num_y_points,5):
                self.main_painter.drawPoint(QPointF(i,j))
    
    def dragMoveEvent(self, event: QDragMoveEvent) -> None:
        pass
    
    def dragLeaveEvent(self, event: QDragLeaveEvent) -> None:
        pass
    
    def mousePressEvent(self, event: QMouseEvent) -> None:
        pass
    
    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        pass