from PySide6.QtWidgets import QVBoxLayout, QWidget,  QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QFrame, QMenu, QMenuBar, QDockWidget, QScrollArea, \
    QStyleOption
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage, QMouseEvent
from typing import Union, Optional

class AppearancePushButton(QPushButton):

    STYLE = """
    """

    def __init__(self) -> object: 
        """
            constructor for Appearance slider widget
        """
        super(AppearancePushButton, self).__init__()


class AppearanceSlider(QSlider):

    def __init__(self) -> object: 
        """
            constructor for Appearance slider widget
        """
        super(AppearanceSlider, self).__init__()


class AppearanceDoubleSpinBox(QDoubleSpinBox):
    
    def __init__(self) -> object:
        super(AppearanceDoubleSpinBox, self).__init__()


class AppearanceSpinBox(QSpinBox):

    def __init__(self) -> object:
        
        super(AppearanceSpinBox, self).__init__()


class AppearanceLineEdit(QLineEdit):

    def __init__(self) -> object:
        super(AppearanceLineEdit, self).__init__()


class AppearanceLabel(QLabel):

    def __init__(self) -> object:
        super(AppearanceLabel, self).__init__()

class AppearanceGroupBox(QGroupBox):

    def __init__(self) -> object:
        super(AppearanceGroupBox, self).__init__()

class AppearanceComboBox(QComboBox):
    
    def __init__(self) -> object:
        super(AppearanceComboBox, self).__init__()

# ----------------------------------------- MENU WIDGETS ----------------------------------

class AppearanceTabWidget(QTabWidget):

    def __init__(self) -> object:
        
        super(AppearanceTabWidget, self).__init__(self)

class AppearanceScrollArea(QScrollArea):
    
    def __init__(self):
        super(AppearanceScrollArea, self).__init__()

class AppearanceToolBar(QToolBar, QToolButton):

    
    def __init__(self):
        super(AppearanceToolBar, self).__init__()

class AppearanceMenu(QMenu, QMenuBar):

    def __init__(self):
        super(AppearanceMenu, self).__init__()


class AppearanceDockableWidget(QDockWidget):

    def __init__(self):
        super(AppearanceDockableWidget, self).__init__()
