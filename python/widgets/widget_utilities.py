from PySide6.QtWidgets import QVBoxLayout, QWidget,  QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QFrame, QMenu, QMenuBar, QDockWidget, QScrollArea, \
    QStyleOption
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage, QMouseEvent
from typing import Union, Optional
from utility_functions import UtilFuncs, VerticalSeparator, HorizontalSeparator

# ------------------------------------------ WIDGET ITEMS ------------------------------------------

class CustomPushButton(QPushButton):
    
     STYLE_SHEET = """
                """
    
     ABOUT_CUSTOM_SLIDER = """
        This is a custom push button with its own style widgets                        
        and utility functions overloaded from the QPushButton 
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
        super(CustomPushButton, self).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)

        if title is not None:
            self.setText(title)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(self.STYLE_SHEET)

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
                 orientation : Qt.Orientation=None,
                 minimum : int = None, maximum : int = None,
                 increment : float=None, multiplier: float=None,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_btn : QPushButton=None, linked_spin_box: QSpinBox=None,
                 linked_double_spin_box : QDoubleSpinBox = None) -> object: 
        """
            constructor for custom slider widget
        """
        super(CustomSlider, self).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)
        
        if orientation is not None:
            self.setOrientation(orientation)

        self.min = minimum
        self.max = maximum
        self.val = self.value()

        if multiplier is not None:
            self.multiplier = multiplier
        if minimum is not None:
            self.setMinimum(minimum)
        if maximum is not None:
            self.setMaximum(maximum)
        if increment is not None:
            self.setSizeIncrement(increment)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)
        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(self.STYLE_SHEET)

        self.btn = linked_btn
        self.spin_box = linked_spin_box
        self.double_spin_box = linked_double_spin_box

        self.activate_slider_links()

    # ----------------------------------- SLIDER GETTERS AND SETTERS --------------------------
    @property
    def spin_box(self):
        return self._spin_box
    
    @spin_box.setter
    def spin_box(self, value : QSpinBox=None):
        self._spin_box = value

    @property
    def btn(self):
        return self._btn
    
    @btn.setter
    def btn(self, value : QPushButton=None):
        self._btn = value

    @property
    def double_spin_box(self):
        return self._double_spin_box
    
    @double_spin_box.setter
    def double_spin_box(self, value : QDoubleSpinBox=None):
        self._double_spin_box = value

    # ---------------------------------- BUTTON SLIDERS ---------------------------------
    def activate_slider_links(self):
        """
            activates slots for slider
        """
        self.sliderMoved.connect(self.slider_moved)
        self.sliderPressed.connect(self.slider_pressed)
        self.sliderReleased.connect(self.slider_released)

    # -------------------------------------- BUTTON CALLBACKS ----------------------------

    def slider_moved(self):
        """
            callback function on linked button
        """
        if self.spin_box is not None:
            self.spin_box.setValue(int(self.value()*self.multiplier))
        if self.double_spin_box is not None:
            self.double_spin_box.setValue(self.value()*self.multiplier)

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

    ABOUT_DOUBLE_SPIN_BOX = """
        This is a custom double spin box with its own style widgets                        
        and utility functions overloaded from the QDoubleSpinBox 
        widget class
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, value : float= None, multiplier: float=None,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_spin_box: QSpinBox=None, linked_slider : QSlider = None) -> object:
        super(CustomDoubleSpinBox, self).__init__()
        
        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)

        if multiplier is not None:
            self.multiplier = multiplier
        if minimum is not None:
            self.setMinimum(minimum)
        if maximum is not None:
            self.setMaximum(maximum)
        if increment is not None:
            self.setSingleStep(increment)
        if value is not None:
            self.setValue(value)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(self.STYLE_SHEET)

        self.spin_box = linked_spin_box
        self.slider = linked_slider

        self.activate_spin_box_links()

    # -------------------------- DOUBLE SPIN BOX GETTERS AND SETTERS ---------------------
    @property
    def spin_box(self):
        return self._spin_box
    
    @spin_box.setter
    def spin_box(self, value: QSpinBox=None):
        self._spin_box = value

    @property
    def slider(self):
        return self._slider
    
    @slider.setter
    def slider(self, value: QSlider=None):
        self._slider = value

    # --------------------------------- SLOTS AND CALLBACKS -------------------------------

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
        if self.spin_box is not None:
            self.spin_box.setValue(int(self.value()*self.multiplier))
        if self.slider is not None:
            self.slider.setValue(int(self.value()*self.multiplier))

    def __str__(self):
        return self.ABOUT_DOUBLE_SPIN_BOX
    
    def __repr__(self):
        print(self.ABOUT_DOUBLE_SPIN_BOX)

class CustomSpinBox(QSpinBox):

    STYLE_SHEET = """
    """

    ABOUT_SPIN_BOX = """
        This is a custom spin box with its own style widgets                        
        and utility functions overloaded from the QSpinBox 
        widget class
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None, multiplier : int =None,
                 minimum : float = None, maximum : float = None,
                 increment : float=None, value : float = None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        
        super(CustomSpinBox, self).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)
        
        if multiplier is not None:
            self.multiplier = multiplier
        if minimum is not None:
            self.setMinimum(minimum)
        if maximum is not None:
            self.setMaximum(maximum)
        if increment is not None:
            self.setSingleStep(increment)
        if value is not None:
            self.setValue(value)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(self.STYLE_SHEET)

        self.double_spin_box = linked_double_spin_box
        self.slider = linked_slider

        self.activate_all_slots()

    # ------------------------------ SPIN BOX GETTERS AND SETTERS ------------------------------
    @property
    def double_spin_box(self):
        return self._double_spin_box
    
    @double_spin_box.setter
    def double_spin_box(self, value: QDoubleSpinBox=None):
        self._double_spin_box = value

    @property
    def slider(self):
        return self._slider
    
    @slider.setter
    def slider(self, value: QSlider=None):
        self._slider = value

    # ------------------------------ SLOTS AND CALLBACKS ------------------------------

    def activate_all_slots(self):
        """
            activates all slot callback methods
        """
        self.valueChanged.connect(self.value_changed)

    def value_changed(self):
        """
            value changed slot
        """
        if self.double_spin_box is not None:
            self.double_spin_box.setValue(float(self.value()*self.multiplier))
        if self.slider is not None:
            self.slider.setValue(self.value()*self.multiplier)

class CustomLineEdit(QLineEdit):

    STYLE_SHEET = """
    """

    ABOUT_LINE_EDIT = """
        This is a custom line edit with its own style widgets                        
        and utility functions overloaded from the QLineEdit 
        widget class
    """

    def __init__(self, style_sheet : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None,
                 size_policy : QSizePolicy=None, text : str=None,
                 minimum : float = None, maximum : float = None,
                 increment : QSize=None, 
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomLineEdit, self).__init__()
    
        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)

        if text is not None:
            self.setText(text)

    def __str__(self):
        return self.ABOUT_LINE_EDIT
    
    def __repr__(self):
        print(self.ABOUT_LINE_EDIT)

class CustomLabel(QLabel):
    STYLE_SHEET = """
    """

    ABOUT_LABEL = """
        This is a custom label with its own style widgets                        
        and utility functions overloaded from the QLabel 
        widget class
    """

    def __init__(self, style_sheet : str=None, title : str=None,
                 minimum : float = None, maximum : float = None,
                 increment : int=None, label_control : bool = False,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 size_policy : QSizePolicy.Policy = None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_spin : QSpinBox=None, linked_double_spin : QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomLabel, self).__init__()

        if title is not None:
            self.setText(title)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if size_policy is not None:
            self.setSizePolicy(size_policy)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)

        self.mouse_pressed = False
        self.slider = linked_slider
        self.spin = linked_spin
        self.double_spin = linked_double_spin

        self.label_control = label_control
        self.increment_size = increment
        self.curr_tick = 0
        self.start_increment = 0
    
    @property
    def slider(self):
        return self._slider

    @slider.setter
    def slider(self, value : QSlider=None) -> None:
        """
            sets internal slider object instance
        """
        self._slider = value

    @property
    def spin(self):
        return self._spin
    
    @spin.setter
    def spin(self, value : QSpinBox=None) -> None:
        """
            sets double spin box object instance
        """
        self._spin = value

    @property
    def double_spin(self):
        return self._double_spin
    
    @double_spin.setter
    def double_spin(self, value : QDoubleSpinBox=None) -> None:
        """
            sets double spin box object instance
        """
        self._double_spin = value

    def mousePressEvent(self, ev: QMouseEvent) -> None:
        """
            callback event when mouse is pressed down
            by the user
        """
        if self.underMouse() is True:
            if self.label_control:
                self.mouse_pressed = True
    
    def mouseMoveEvent(self, ev: QMouseEvent) -> None:
        """
            callback event when mouse is being 
            moved by the user
        """
        if self.mouse_pressed:
            if self.double_spin is not None:
                self.double_spin.setValue(float(self.start_increment))
                self.start_increment += self.increment_size
            if self.slider is not None:
                self.slider.setValue(self.curr_tick)
            if self.spin is not None:
                self.spin.setValue(self.start_increment)
                self.start_increment += self.increment_size

    def mouseReleaseEvent(self, ev: QMouseEvent) -> None:
        """
            callback event when mouse is released by
            the user
        """
        self.mouse_pressed = False

class CustomGroupBox(QGroupBox, UtilFuncs):
    STYLE_SHEET = """
    """

    ABOUT_GROUP_BOX = """
        This is a custom group box with its own style widgets                        
        and utility functions overloaded from the QGroupBox 
        widget class
    """

    def __init__(self, style_sheet : str=None, title : str=None,
                 v_size_policy : QSizePolicy.Policy=None, 
                 h_size_policy : QSizePolicy.Policy=None, default_state : bool = False,
                 fixed_size_policy : bool=False,
                 size_policy : QSizePolicy=None, checkable : bool = False,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None) -> object:
        super(CustomGroupBox, self).__init__()

        if size_policy is not None:
            self.setSizePolicy(size_policy)
            
        if v_size_policy is not None and h_size_policy is not None:
            self.setSizePolicy(v_size_policy, h_size_policy)
        self.setCheckable(checkable)

        if default_state is True:
            self.set_default_state(self)
        if fixed_size_policy is True:
            self.set_fixed_size_policy(self)

        if title is not None:
            self.setTitle(title)

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(self.STYLE_SHEET)

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

    ABOUT_COMBO_BOX = """
        This is a custom combo box with its own style widgets                        
        and utility functions overloaded from the QComboBox 
        widget class
    """

    def __init__(self, style_sheet : str=None, add_item : bool=None,
                 add_items : bool=None, item : str=None, items_to_add : list =None,
                 minimum_size : QSize=None, maximum_size : QSize=None,
                 width : QSize=None, height : QSize=None, size : QSize=None,
                 linked_double_spin_box: QDoubleSpinBox=None,
                 linked_slider : QSlider = None) -> object:
        super(CustomComboBox, self).__init__()

        self.all_items_to_add = items_to_add
        if add_items is True:
            self.add_items_from_list()
        if add_item is True:
            self.addItem(item)
        if style_sheet is not None:
            self.setStyleSheet(style_sheet)
        
        self.currentIndexChanged.connect(self.index_changed)
        
        self.setStyleSheet(self.STYLE_SHEET)
            
    def add_items_from_list(self):
        """
            add items from list helper method
        """
        for item in self.all_items_to_add:
            self.addItem(item)

    def index_changed(self):
        """
            slot callback when index item changed
        """
        pass

# ----------------------------------------- MENU WIDGETS ----------------------------------

class CustomTabWidget(QTabWidget):
    
    STYLE_SHEET = """

    """

    ABOUT_TAB_WIDGET = """
        This is a custom tab widget with its own style widgets                        
        and utility functions overloaded from the QTabWidgets 
        widget class
    """

    def __init__(self, text : str=None,
                 style_sheet : str=None, closable: bool=False,
                 dockable : bool = False, size_policy : QSizePolicy.Policy=None,
                 minimum_size : int=None, maximum_size : int=None,
                 width : int=None, height : int=None,
                 size : QSize=None) -> object:
        
        super(CustomTabWidget, self).__init__(self)

        self.setTabsClosable(closable)
        
        if size_policy is not None:
            self.setSizePolicy(size_policy)
        if size is not None:
            self.setFixedSize(size)
        if text is not None:
            self.tabText(text)
        if dockable is True:
            self.isDockable = True
        if style_sheet is not None:
            self.setStyleSheet(style_sheet)
        
        self.isDockable = False

        if minimum_size is not None:
            self.setMinimumSize(minimum_size)   
        if maximum_size is not None:
            self.setMaximumSize(maximum_size)

        if width is not None:
            self.setFixedWidth(width)
        if height is not None:
            self.setFixedHeight(height)
        if size is not None:
            self.setFixedSize(size)

        self.setStyleSheet(self.STYLE_SHEET)

class CustomScrollArea(QScrollArea):
    
    STYLE_SHEET = """

    """

    ABOUT_SCROLL_AREA = """
        This is a custom scroll area with its own style widgets                        
        and utility functions overloaded from the QScrollArea 
        widget class
    """
    
    def __init__(self):
        super(CustomScrollArea, self).__init__()

    def __str__(self):
        return self.ABOUT_SCROLL_AREA
    
    def __repr__(self):
        print(self.ABOUT_SCROLL_AREA)

class CustomToolBar(QToolBar, QToolButton):
    
    STYLE_SHEET = """

    """

    ABOUT_TOOL_BAR = """
        This is a custom tool bar with its own style widgets                        
        and utility functions overloaded from the QToolBar and  
        QToolButton widget class
    """
    
    def __init__(self):
        super(CustomToolBar, self).__init__()

    def __str__(self):
        return self.ABOUT_TOOL_BAR
    
    def __repr__(self):
        print(self.ABOUT_TOOL_BAR)

class CustomMenu(QMenu, QMenuBar):
    
    STYLE_SHEET = """

    """

    ABOUT_MENU = """
        This is a custom combo box with its own style widgets                        
        and utility functions overloaded from the QComboBox 
        widget class
    """

    def __init__(self):
        super(CustomMenu, self).__init__()

    def __str__(self):
        return self.ABOUT_MENU
    
    def __repr__(self):
        print(self.ABOUT_MENU)

class CustomDockableWidget(QDockWidget):
    
    STYLE_SHEET = """

    """

    ABOUT_DOCK_WIDGET = """
        This is a custom combo box with its own style widgets                        
        and utility functions overloaded from the QComboBox 
        widget class
    """

    def __init__(self):
        super(CustomDockableWidget, self).__init__()

    def __str__(self):
        return self.ABOUT_DOCK_WIDGET
    
    def __repr__(self):
        print(self.ABOUT_DOCK_WIDGET)
