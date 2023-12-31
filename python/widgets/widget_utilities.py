from PySide6.QtWidgets import QPushButton, \
    QToolButton, QToolBar, QTabWidget, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox,  \
    QSizePolicy,  QMenu, QMenuBar, QDockWidget, QScrollArea
from PySide6.QtCore import QSize, Qt
from PySide6.QtGui import QMouseEvent, QAction
from utility_functions import UtilFuncs
from Fluid_Utilities.system import FluidSystem
from appearance import *

# ------------------------------------------ WIDGET ITEMS ------------------------------------------

class CustomPushButton(QPushButton):
     
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
                  linked_double_spin_box : QDoubleSpinBox = None,
                  alias:str = None, system_obj:FluidSystem = None) -> object: 
        """
            constructor for custom slider widget
        """
        super(CustomPushButton, self).__init__()

        styling = AppearancePushButton()
        
        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

        self.setStyleSheet(styling.STYLE)

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

        if alias is not None:
            self.alias = alias
        if system_obj is not None:
            self.system_obj = system_obj

        self.clicked.connect(self.clicked_fn)

     def clicked_fn(self):
        
        if self.alias is not None:
            if isinstance(self.alias, list):
                if self.alias == "gizmos":
                    pass
            else:
                if self.alias == "play":
                    if self.system_obj.finished_caching != True:
                        self.system_obj.start_play = True
                        self.system_obj.stop = False
                        self.system_obj.update_stored_positions()
                    else:
                        self.system_obj.stop = False
                        self.system_obj.start_playforward = True
                        self.system_obj.start_playback = False
                elif self.alias == "stop":
                    self.system_obj.stop = True
                    self.system_obj.start_play = False
                elif self.alias == "backward":
                    self.system_obj.start_playback = True
                    self.system_obj.start_playforward = False
                    self.system_obj.stop = False

class CustomSlider(QSlider):
    
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
                 linked_double_spin_box : QDoubleSpinBox = None,
                 alias:str = None, system_obj:FluidSystem = None) -> object: 
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

        self.btn = linked_btn
        self.spin_box = linked_spin_box
        self.double_spin_box = linked_double_spin_box

        if alias is not None:
            self.alias = alias
        if system_obj is not None:
            self.system_obj = system_obj

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

        self.system_obj.finished_caching = False
        self.alias_callback()

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
    
    def alias_callback(self):
        """
            change attributes in system class object
        """
        if self.alias is not None:
            if self.alias == "particle number":
                self.system_obj.num_particles = self.value()
            elif self.alias == "particle separation":
                self.system_obj.USER_PARAMETERS["grid_separation"] = self.value()
            elif self.alias == "cell size":
                self.system_obj.USER_PARAMETERS["cell_size"] = self.value()
            elif self.alias == "mass":
                self.system_obj.USER_PARAMETERS["mass"] = self.value()
            elif self.alias == "gravity":
                self.system_obj.USER_PARAMETERS["gravity"] = self.value()
            elif self.alias == "buoyancy":
                self.system_obj.USER_PARAMETERS["buoyancy"] = self.value()
            elif self.alias == "viscosity":
                self.system_obj.USER_PARAMETERS["viscosity"] = self.value()
            elif self.alias == "pressure":
                self.system_obj.USER_PARAMETERS["pressure"] = self.value()
            elif self.alias == "mass density":
                self.system_obj.USER_PARAMETERS["mass_density"] = self.value()
            elif self.alias == "speed loss":
                self.system_obj.USER_PARAMETERS["speed_loss"] = self.value()
            elif self.alias == "delta time":
                self.system_obj.USER_PARAMETERS["delta_time"] = self.value() 

    # -------------------------------- CLASS DUNDER METHODS ----------------------------
    def __repr__(self):
        print(self.ABOUT_CUSTOM_SLIDER)
    
    def __str__(self):
        return self.ABOUT_CUSTOM_SLIDER

class CustomDoubleSpinBox(QDoubleSpinBox):

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
                 linked_spin_box: QSpinBox=None, linked_slider : QSlider = None,
                 alias:str = None, system_obj:FluidSystem = None,
                 direction:int = None) -> object:
        
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

        if alias is not None:
            self.alias = alias
        if system_obj is not None:
            self.system_obj = system_obj
        if direction is not None:
            self.direction = direction

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

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

        self.system_obj.finished_caching = False
        self.alias_callback()

    def alias_callback(self):
        if self.alias is not None:
            if self.alias == "particle separation":
                self.system_obj.USER_PARAMETERS["grid_separation"] = self.value()
            elif self.alias == "cell size":
                self.system_obj.USER_PARAMETERS["cell_size"] = self.value()
            elif self.alias == "velocity":
                if self.direction is not None:
                    self.system_obj.USER_PARAMETERS["initial_velocity"][self.direction] = self.value()
            elif self.alias == "acceleration":
                if self.direction is not None:
                    self.system_obj.USER_PARAMETERS["initial_acceleration"][self.direction] = self.value()
            elif self.alias == "tank radius":
                if self.direction is not None:
                    self.system_obj.TANK_ATTRS["dimensions"]["size"][self.direction] = self.value()
            elif self.alias == "tank position":
                if self.direction is not None:
                    self.system_obj.TANK_ATTRS["dimensions"]["location"][self.direction] = self.value()
            elif self.alias == "mass":
                self.system_obj.USER_PARAMETERS["mass"] = self.value()
            elif self.alias == "gravity":
                self.system_obj.USER_PARAMETERS["gravity"] = self.value()
            elif self.alias == "buoyancy":
                self.system_obj.USER_PARAMETERS["buoyancy"] = self.value()
            elif self.alias == "viscosity":
                self.system_obj.USER_PARAMETERS["viscosity"] = self.value()
            elif self.alias == "pressure":
                self.system_obj.USER_PARAMETERS["pressure"] = self.value()
            elif self.alias == "mass density":
                self.system_obj.USER_PARAMETERS["mass_density"] = self.value()
            elif self.alias == "speed loss":
                self.system_obj.USER_PARAMETERS["speed_loss"] = self.value()
            elif self.alias == "delta time":
                self.system_obj.USER_PARAMETERS["delta_time"] = self.value() 

    def __str__(self):
        return self.ABOUT_DOUBLE_SPIN_BOX
    
    def __repr__(self):
        print(self.ABOUT_DOUBLE_SPIN_BOX)

class CustomSpinBox(QSpinBox):

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
                 linked_slider : QSlider = None,
                 alias:str=None, system_obj:FluidSystem=None) -> object:
        
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

        if alias is not None:
            self.alias = alias
        if system_obj is not None:
            self.system_obj = system_obj

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)

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

        self.system_obj.finished_caching = False
        self.alias_callback()

    def alias_callback(self):
        if self.alias is not None:
            if self.alias == "particle number":
                self.system_obj.num_particles = self.value()

class CustomLineEdit(QLineEdit):

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

class CustomAction(QAction):

    def __init__(self,
                 text= None) -> None:

        super(CustomAction, self).__init__()

        if text is not None:
            self.setText(text)

        self.trigger_signals()

    def trigger_signals(self):

        self.changed.connect(self.changed_fn)
        self.checkableChanged.connect(self.checkable_changed)

    def changed_fn(self):
        pass

    def checkable_changed(self):
        pass

class CustomComboBox(QComboBox):

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
                 linked_slider : QSlider = None, alias:str = None,
                 system_obj:FluidSystem = None) -> object:
        super(CustomComboBox, self).__init__()

        self.all_items_to_add = items_to_add
        if add_items is True:
            self.add_items_from_list()
        if add_item is True:
            self.addItem(item)

        if alias is not None:
            self.alias = alias
        if system_obj is not None:
            self.system_obj = system_obj  

        if style_sheet is not None:
            self.setStyleSheet(style_sheet)
        
        self.currentIndexChanged.connect(self.index_changed)
            
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
        self.alias_callback()
        self.system_obj.start_playback = False

    def alias_callback(self):
        if self.alias is not None:
            if self.alias == "solver type":
                self.system_obj.simulation_type = self.currentText()
            elif self.alias == "distribution type":
                self.system_obj.orientation_type = self.currentText()
            elif self.alias == "search method":
                self.system_obj.search_method = self.currentText()
            elif self.alias == "time integrator":
                self.system_obj.time_stepping = self.currentText()

class CustomGroupBox(QGroupBox, UtilFuncs):

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
    
# ----------------------------------------- MENU WIDGETS ----------------------------------

class CustomTabWidget(QTabWidget):

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

class CustomScrollArea(QScrollArea):

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
