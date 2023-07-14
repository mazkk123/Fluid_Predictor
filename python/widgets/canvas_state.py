from PySide6.QtWidgets import  QWidget
from PySide6.QtCore import QSize, QRect , Qt, QPointF
from PySide6.QtGui import QPixmap, QIcon, QImage, QPaintEvent, QColor, \
    QFont, QPainter, QBrush, QFont, QPen, QDragMoveEvent, QDragLeaveEvent, \
    QMouseEvent

class DrawingCanvas(QWidget):

    def __init__(self,
                 canvas_size : QRect,
                 drawing_colour : QColor,
                 font_type: QFont,
                 font_size: int) -> None:
        super().__init__()

        self.setGeometry(canvas_size)
        self.drawing_color = drawing_colour
        self.drawing_font = font_type
        self.font_size = font_size

    def paintEvent(self, event: QPaintEvent) -> None:
        """
            controls main painting to the widget window
        """
        self.fps_text_to_draw = "FPS: None"
        self.pt_num_text_to_draw = "No. Points: None"
        self.projection_text_to_draw = "Projection: persp"

        self.font_to_draw = QFont()
        self.font_to_draw.setFamily("Arial")
        self.font_to_draw.setPixelSize(self.font_size)

        self.main_painter = QPainter(self)
        self.main_painter.setFont(self.font_to_draw)

        self.main_painter.drawText(825,20, self.fps_text_to_draw)
        self.main_painter.drawText(825,35, self.pt_num_text_to_draw)
        self.main_painter.drawText(825,50, self.projection_text_to_draw)

        self.main_painter.end()
    
    def dragMoveEvent(self, event: QDragMoveEvent) -> None:
        pass
    
    def dragLeaveEvent(self, event: QDragLeaveEvent) -> None:
        pass
    
    def mousePressEvent(self, event: QMouseEvent) -> None:
        pass
    
    def mouseReleaseEvent(self, event: QMouseEvent) -> None:
        pass