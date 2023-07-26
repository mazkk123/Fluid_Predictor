from PySide6.QtWidgets import QVBoxLayout, QWidget,  QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QFrame, QMenu, QMenuBar, QDockWidget, QScrollArea, \
    QStyleOption, QDialog
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage
from typing import Union, Optional

# ----------------------------------------- SEPARATOR OPTIONS ------------------------------------

class VerticalSeparator(QFrame):
    
    def __init__(self) -> None:
        super(VerticalSeparator, self).__init__()
        self.setFrameShape(QFrame.Shape.VLine)
        self.setFrameShadow(QFrame.Shadow.Sunken)

class HorizontalSeparator(QFrame):
    
    def __init__(self) -> None:
        super(VerticalSeparator, self).__init__()
        self.setFrameShape(QFrame.Shape.VLine)
        self.setFrameShadow(QFrame.Shadow.Sunken)

class UtilFuncs(VerticalSeparator, HorizontalSeparator):

    ABOUT_UTILITY = """
                    provides abstract utility functions to main   
                    widget class calls from the main widgets module.
                    """

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
                             tool_button: QToolButton=None, 
                             sub_menu: str="Bottom_Bar", 
                             image_name : str=None) -> None:
        """
            adds icons based on their name to local list
        """
        if widget is not None:
            if image_name is not None:
                path = "images/" + sub_menu + "/" + str(image_name)
                widget.setIcon(QIcon(path))
        if tool_button is not None:
            if image_name is not None:
                path = "images/" + sub_menu + "/" + str(image_name)
                tool_button.setIcon(QIcon(path))

    def add_multiple_layouts(self, 
                             parent : QGroupBox=None,
                             parent_layout : QVBoxLayout or QHBoxLayout=None,
                             child_layouts : list=None) -> None:
        """
            adds multiple child layouts into a given vbox
            or hbox layout clas
        """
        parent.setLayout(parent_layout)
        for child in child_layouts:
            if child.isWidgetType() is False:
                parent_layout.addLayout(child)

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

    def create_frame_bars(self, layout: QHBoxLayout=None,
                          start: int=None, end: int=None,
                          line_edit_top: QLineEdit=None,
                          line_edit_side: QLineEdit=None,
                          scale_factor: float=6) -> None:
        """
            creates vertical separator corresponding to the number of frames
            between the spacers
        """
        separation_length = end - start

        try:
            line_edit_top_val = int(line_edit_top.text())
            line_edit_side_val = int(line_edit_side.text())
            
            num_frames = line_edit_side_val - line_edit_top_val
            spacing_between_frames = separation_length // num_frames

            for i in range(num_frames-1):
                separator_obj = VerticalSeparator()
                separator_obj.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                        QSizePolicy.Policy.Expanding))
                separator_obj.setFixedWidth(spacing_between_frames-((separator_obj.lineWidth() +
                                                                    separator_obj.midLineWidth())*
                                                                    scale_factor))
                separator_obj.setLineWidth(3)
                separator_obj.setMidLineWidth(3)
                layout.addWidget(separator_obj)

        except (ValueError, TypeError):
            pass


    def adjust_label_spacing(self, box_widget: QGroupBox=None,
                             scale_factor : int=8) -> None:
        """
            adjusts label spacing by group box proportional
            amount
        """
        widget_pos_list = []

        for child in box_widget.children():
                if isinstance(child, QLabel):
                    widget_pos_list.append(child.pos().x() + child.width())

        maximum_dist = max(widget_pos_list)/scale_factor

        for child in box_widget.children():
            if child.isWidgetType():
                if isinstance(child, QLabel):
                    child.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                    QSizePolicy.Policy.Fixed))
                    child.setFixedWidth(maximum_dist)
                    child.setAlignment(Qt.AlignmentFlag.AlignHCenter)

    def adjust_label_spacing_layout(self, box_widget: QWidget=None,
                                    scale_factor : int=8) -> None:
        """
            adjusts label spacing by group box proportional
            amount
        """
        widget_pos_list = []

        for child in box_widget.children():
                if isinstance(child, QLabel):
                    widget_pos_list.append(child.pos().x() + child.width())

        maximum_dist = max(widget_pos_list)/scale_factor

        for child in box_widget.children():
            if child.isWidgetType():
                if isinstance(child, QLabel):
                    child.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                    QSizePolicy.Policy.Fixed))
                    child.setFixedWidth(maximum_dist)
                    child.setAlignment(Qt.AlignmentFlag.AlignHCenter)

    def __repr__(self) -> str:
        return self.ABOUT_UTILITY

