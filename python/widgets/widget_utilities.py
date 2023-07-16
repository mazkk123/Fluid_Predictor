from PySide6.QtWidgets import QVBoxLayout, QWidget,  QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QFrame, QMenu, QMenuBar, QDockWidget, QScrollArea, \
    QStyleOption
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage
from typing import Union, Optional
from utility_functions import UtilFuncs

# ------------------------------------------ WIDGET ITEMS ------------------------------------------

class CustomPushButton(QPushButton):
    
     STYLE_SHEET = """
                """
    
     ABOUT_CUSTOM_SLIDER = """
        This is a custom slider with its own style widgets                        
        and utility functions overloaded from the QSlider 
        widget class
    """

     def __init__(self, title: str=None, style_sheet : str=None,
                  v_size_policy : QSizePolicy.Policy=None, 
                  h_size_policy : QSizePolicy.Policy=None,
                  size_policy : QSizePolicy=None,
                  minimum_size : QSize=None,maximum_size : QSize=None,
                  width : QSize=None, height : QSize=None, 
                  linked_btn : QPushButton=None, linked_spin_box: QSpinBox=None,
                  linked_double_spin_box : QDoubleSpinBox = None) -> object: 
        """
            constructor for custom slider widget
        """
        super(CustomPushButton).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

        self.setText(title)

        self.setMaximumSize(maximum_size)
        self.setMinimumSize(minimum_size)
        self.setFixedWidth(width)
        self.setFixedHeight(height)

        if style_sheet is not None:
            self.setStyleSheet(self.STYLE_SHEET)

        self.setStyleSheet(style_sheet)

class CustomSlider(QSlider):

    STYLE_SHEET = """
                """
    
    ABOUT_CUSTOM_SLIDER = """
        This is a custom slider with its own style widgets                        
        and utility functions overloaded from the QSlider 
        widget class
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum : int = None, maximum : int = None,
                 increment : float=None, 
                 minimum_size : QSize=None,maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_btn : QPushButton=None, linked_spin_box: QSpinBox=None,
                 linked_double_spin_box : QDoubleSpinBox = None) -> object: 
        """
            constructor for custom slider widget
        """
        super(CustomSlider).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

        self.styleSheet = style_sheet
        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setSizeIncrement(increment)

        self.min = minimum
        self.max = maximum
        self.val = self.value

        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setSingleStep(1)
        self.setMaximumSize(maximum_size)
        self.setMinimumSize(minimum_size)
        self.setFixedWidth(width)
        self.setFixedHeight(height)

        if style_sheet is not None:
            self.setStyleSheet(self.STYLE_SHEET)

        self.setStyleSheet(style_sheet)

        self.linked_btn = linked_btn
        self.linked_spin_box = linked_spin_box
        self.linked_double_spin_box = linked_double_spin_box

        self.activate_slider_links()

    # ---------------------------------- BUTTON SLIDERS ---------------------------------
    def activate_slider_links(self):
        """
            activates slots for slider
        """
        self.valueChanged.connect(self.value_changed)
        self.sliderMoved.connect(self.slider_moved)
        self.sliderPressed.connect(self.slider_pressed)
        self.sliderReleased.connect(self.slider_released)

    # -------------------------------------- BUTTON CALLBACKS ----------------------------
    def value_changed(self):
        """
            callback function on linked button
        """
        self.linked_spin_box.setValue(int(self.value))
        self.linked_double_spin_box.setValue(self.value)

    def slider_moved(self):
        """
            callback function on linked button
        """
        self.linked_spin_box.setValue(int(self.value))
        self.linked_double_spin_box.setValue(self.value)

    def slider_pressed(self):
        """
            callback function on linked button
        """
        pass

    def slider_released(self):
        """
            callback function on linked button
        """
        pass

    # -------------------------------- CLASS DUNDER METHODS ----------------------------
    def __repr__(self):
        print(self.ABOUT_CUSTOM_SLIDER)
    
    def __str__(self):
        return self.ABOUT_CUSTOM_SLIDER

class CustomDoubleSpinBox(QDoubleSpinBox):
    
    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_spin_box: QSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomDoubleSpinBox).__init__()
        
        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setSizeIncrement(increment)

        self.setMinimumSize(minimum_size)
        self.setMaximumSize(maximum_size)

        self.setFixedWidth(width)
        self.setFixedHeight(height)
        self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(self.STYLE_SHEET)

        self.setStyleSheet(style_sheet)

        self.linked_spin_box = linked_spin_box
        self.linked_slider = linked_slider

        self.activate_spin_box_links()

    def activate_spin_box_links(self):
        """
            activates all callback slot links to different
            widgets    
        """
        self.valueChanged.connect(self.value_changed)

    def value_changed(self):
        """
            value changed callback slot 
        """
        self.linked_spin_box.setValue(self.value)
        self.linked_slider.setValue(int(self.value))

class CustomSpinBox(QSpinBox):

    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomSpinBox).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

        self.setMinimum(minimum)
        self.setMaximum(maximum)
        self.setSizeIncrement(increment)

        self.setMinimumSize(minimum_size)
        self.setMaximumSize(maximum_size)

        self.setFixedWidth(width)
        self.setFixedHeight(height)
        self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(self.STYLE_SHEET)

        self.setStyleSheet(style_sheet)

        self.linked_double_spin_box = linked_double_spin_box
        self.linked_slider = linked_slider

        self.activate_all_slots()

    def activate_all_slots(self):
        """
            activates all slot callback methods
        """
        self.valueChanged.connect(self.value_changed)

    def value_changed(self):
        """
            value changed slot
        """
        self.linked_double_spin_box.setValue(float(self.value))
        self.linked_slider.setValue(self.value)

class CustomLineEdit(QLineEdit):

    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomLineEdit).__init__()
    
        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

class CustomLabel(QLabel):
    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None, title : str=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomLabel).__init__()

        self.setText(title)

class CustomGroupBox(QGroupBox, UtilFuncs):
    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None, title : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None) -> object:
        super(CustomGroupBox).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        self.setSizePolicy(v_size_policy, h_size_policy)

        self.setTitle(title)

        self.setMinimumSize(minimum_size)   
        self.setMaximumSize(maximum_size)

        self.setFixedWidth(width)
        self.setFixedHeight(height)
        self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(self.STYLE_SHEET)

        self.setStyleSheet(style_sheet)

        self.activate_all_slots()

    def activate_all_slots(self):
        """
            activate all slots associated with group box widget
        """
        self.toggled.connect(self.toggled_slot)

    def toggled_slot(self):
        """
            toggled group box slot
        """
        self.hide_group_box_widgets(self)

class CustomComboBox(QComboBox):
    STYLE_SHEET = """
    """

    ABOUT = """
    """

    def __init__(self, style_sheet : str=None, add_item : bool=None,
                 add_items : bool=None, item : str=None, items_to_add : list=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomComboBox).__init__()

        self.all_items_to_add = items_to_add
        if add_items is True:
            self.add_items_from_list(self.all_items_to_add)
        elif add_item is True:
            self.addItem(item)
        
    def add_items_from_list(self):
        """
            add items from list helper method
        """
        for item in self.all_items_to_add:
            self.addItem(item)

# ----------------------------------------- MENU WIDGETS ----------------------------------

class CustomTabWidget(QTabWidget):
    pass

class CustomScrollArea(QScrollArea):
    pass

class CustomToolBar(QToolBar, QToolButton):
    pass

class CustomMenu(QMenu, QMenuBar):
    pass

class CustomDockableWidget(QDockWidget):
    pass

