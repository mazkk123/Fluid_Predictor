from PySide6.QtWidgets import QVBoxLayout, QWidget,  QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QScrollArea, QGraphicsScene, QGraphicsView
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage

class Utilities:

    def __init__(self):
        pass

    def add_icons_to_actions(self, icon: QToolBar=None, 
                             sub_menu: str ="Toolbar", 
                             image_name : str=None) -> None:
        """
            adds icons based on their name to local list
        """
        if icon is not None:
            if image_name is not None:
                path = "images/" + sub_menu + "/" + str(image_name)
                icon.setIcon(QIcon(path))

    def add_icons_to_widgets(self, widget : QWidget=None, 
                             sub_menu: str="Bottom_Bar", 
                             image_name : str=None) -> None:
        """
            adds icons based on their name to local list
        """
        if widget is not None:
            if image_name is not None:
                path = "images/" + sub_menu + "/" + str(image_name)
                widget.setIcon(QIcon(path))

    def configure_sliders(self, 
                          slider : QSlider,
                          start_value : float,
                          increment : float) -> None:
        """
            configures the values of sliders    
        """
        slider.setTickPosition(start_value)
        slider.setSizeIncrement(increment)
        slider.setValue(start_value)

    def connect_slider_with_sBox(self, 
                               slider : QSlider,
                               spin_box : QSpinBox,
                               ) -> None:
        """
            connects slider values to spin box values
        """
        if slider.sliderMoved():
            spin_box.setValue(slider.value)
        if spin_box.valueChanged():
            slider.setValue(spin_box.value)

    def set_default_state(self, group_box : QGroupBox=None) -> None:
        """
            set default group box children states
        """
        if group_box is not None:
            for child in group_box.children():
                if child.isWidgetType():
                    child.setHidden(True)
            group_box.setChecked(False)

    def set_fixed_size_policy(self, group_box: QGroupBox=None) -> None:
        """
            set all group size policies to expanding fixed.
        """
        if group_box is not None:
            group_box.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Expanding, 
                                                QSizePolicy.Policy.Fixed))

    def hide_group_box_widgets(self, group_box : QGroupBox=None) -> None:
        """
            generic hide function for group box widget children
        """
        if group_box is not None:
            for child in group_box.children():
                if child.isWidgetType():
                    if group_box.isChecked():
                        child.setHidden(False)
                    else:
                        child.setHidden(True)