from PySide6.QtWidgets import  QWidget, QSizePolicy
from PySide6.QtCore import QSize, QRect , Qt, QPointF, QPoint
from PySide6.QtGui import QPixmap, QIcon, QImage, QPaintEvent, QColor, \
    QFont, QPainter, QBrush, QFont, QPen, QDragMoveEvent, QDragLeaveEvent, \
    QMouseEvent, QTransform

class DrawingCanvas(QWidget):

    def __init__(self,
                 canvas_size : QRect,
                 drawing_colour : QColor,
                 font_type: QFont,
                 font_size: int) -> None:
        super().__init__()

        """ self.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                       QSizePolicy.Policy.Expanding))
        self.setFixedWidth(canvas_size.width()) """
        self.setGeometry(canvas_size)
        self.drawing_color = drawing_colour
        self.drawing_font = font_type
        self.font_size = font_size

    def paintEvent(self, event: QPaintEvent) -> None:
        """
            controls main painting to the widget window
        """
        self.draw_feature_text()

        self.render_points()

        self.main_painter.end()

    def draw_feature_text(self) -> None:
        """
            draws feature text on top right hand side of screen        
        """
        self.fps_text_to_draw = "FPS: None"
        self.pt_num_text_to_draw = "No. Points: None"
        self.projection_text_to_draw = "Projection: persp"

        self.font_to_draw = QFont()
        self.font_to_draw.setFamily("Arial")
        self.font_to_draw.setPixelSize(self.font_size)

        self.main_painter = QPainter(self)
        self.main_painter.setFont(self.font_to_draw)

        self.main_painter.drawText(800,20, self.fps_text_to_draw)
        self.main_painter.drawText(800,35, self.pt_num_text_to_draw)
        self.main_painter.drawText(800,50, self.projection_text_to_draw)

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