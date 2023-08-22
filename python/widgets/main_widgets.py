from PySide6.QtWidgets import QVBoxLayout, QWidget, QMainWindow, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, \
    QSpacerItem, QSizePolicy,  QScrollArea, QSplitter, QStatusBar, \
    QDockWidget
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QColor, QFont, QKeyEvent
from render_scene import RenderScene
from utility_functions import *
from widget_utilities import *
import sys

sys.path.append("C:\\Users\\Student\\OneDrive - Bournemouth University\\Desktop\\Personal\\Python\\Fluid_Predictor\\python\\")

from Fluid_Utilities.system import FluidSystem

class MainWindow(QMainWindow, UtilFuncs):

    def __init__(self, app):
        super().__init__()
        self.app = app
        
        self.setGeometry(QRect(150, 30, 1280, 700))
        self.setWindowTitle("Fluid Application")
        self.window_icon = QPixmap("images/Toolbar/fluid.png")
        self.setWindowIcon(self.window_icon)
        self.system_obj = FluidSystem()

        self.initUI()

        self.setCentralWidget(self.main_splitter_widget)

    def initUI(self):
        """
            initializes calls to most UI elements on screen
        """
        self.add_menu_bars()
        self.add_tool_bars()
        self.add_left_tool_bar()
        self.add_status_bar()

        self.main_options_widget = QWidget()
        self.main_options_widget.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Expanding,
                                                           QSizePolicy.Policy.Fixed))
        self.main_layout = QVBoxLayout()
        self.tab_v_layout = QVBoxLayout()

        self.create_scroll_area_widget()
        self.initialize_main_group_boxes()

        self.main_options_widget.setLayout(self.main_layout)
        self.create_main_canvas()
        self.create_frame_controls()

        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.fluid_options_dock_w)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.fluid_pred_dock_w)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.animation_controls_dock_w)

        self.main_splitter_widget = QSplitter(Qt.Orientation.Horizontal)

        self.main_splitter_widget.addWidget(self.frame_and_graphics_splitter)
        self.main_splitter_widget.addWidget(self.fluid_options_splitter_w)

        self.main_splitter_widget.setChildrenCollapsible(False)

        self.fluid_options_dock_w.setMinimumWidth(300)
        self.fluid_options_dock_w.setMinimumHeight(400)

    def add_status_bar(self):
        """
            add state bar text to show current simulation steps
        """
        self.main_state_bar = QStatusBar()

        self.simulation_text = CustomLabel(title="Simulation is running....")
        self.main_state_bar.addWidget(self.simulation_text)

        self.setStatusBar(self.main_state_bar)

    def add_left_tool_bar(self):
        """
            add movement and cursor viewport handles on the left
        """

        self.left_tool_dock_w = QDockWidget("Tools")
        self.left_tool_dock_w.setSizePolicy(QSizePolicy.Policy.Minimum,
                                            QSizePolicy.Policy.Expanding)
        
        self.left_tool_bar_w = QWidget()
        self.left_bar_v_layout = QVBoxLayout()
        self.left_tool_bar_w.setLayout(self.left_bar_v_layout)
        self.left_tool_dock_w.setWidget(self.left_tool_bar_w)

        self.cursor_gizmo = QToolButton()
        self.cursor_gizmo.setMinimumSize(QSize(50,50))
        self.add_icons_to_widgets(tool_button=self.cursor_gizmo,
                                  sub_menu="Gizmos", image_name="select.png")
        self.translate_gizmo = QToolButton()
        self.translate_gizmo.setMinimumSize(QSize(50,50))
        self.add_icons_to_widgets(tool_button=self.translate_gizmo,
                                  sub_menu="Gizmos", image_name="move.png")
        self.rotate_gizmo = QToolButton()
        self.rotate_gizmo.setMinimumSize(QSize(50,50))
        self.add_icons_to_widgets(tool_button=self.rotate_gizmo,
                                  sub_menu="Gizmos", image_name="rotate.png")
        self.scale_gizmo = QToolButton()
        self.scale_gizmo.setMinimumSize(QSize(50,50))
        self.add_icons_to_widgets(tool_button=self.scale_gizmo,
                                  sub_menu="Gizmos", image_name="scale.png")
        self.current_gizmo = QToolButton()
        self.current_gizmo.setMinimumSize(QSize(50,50))
        self.add_icons_to_widgets(tool_button=self.current_gizmo,
                                  sub_menu="Gizmos", image_name="select.png")

        self.left_bar_v_layout.addWidget(self.cursor_gizmo)
        self.left_bar_v_layout.addWidget(self.translate_gizmo)
        self.left_bar_v_layout.addWidget(self.rotate_gizmo)
        self.left_bar_v_layout.addWidget(self.scale_gizmo)
        self.left_bar_v_layout.addStretch(300)
        self.left_bar_v_layout.addWidget(self.current_gizmo)

        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.left_tool_dock_w)

    def create_scroll_area_widget(self):
        """
            creates and configures the initial scroll area widget binded to each
            tab
        """

        self.fluid_options_dock_w = QDockWidget("Fluid Options")
        self.fluid_options_dock_w.setFloating(True)
        self.fluid_options_dock_w.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)

        self.fluid_sim_scroll_area = QScrollArea()
        self.fluid_sim_scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.fluid_sim_scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.fluid_sim_scroll_area.setWidgetResizable(True)

        self.fluid_sim_scroll_area.setWidget(self.main_options_widget)

        self.fluid_options_dock_w.setWidget(self.fluid_sim_scroll_area)

        self.fluid_pred_dock_w = QDockWidget("Fluid Prediction")
        
        self.fluid_pred_options_widget = QWidget()

        self.fluid_pred_w()

        self.fluid_pred_dock_w.setWidget(self.fluid_pred_options_widget)

        self.fluid_options_splitter_w = QSplitter(Qt.Orientation.Vertical)

        self.fluid_options_splitter_w.addWidget(self.fluid_options_dock_w)
        self.fluid_options_splitter_w.addWidget(self.fluid_pred_dock_w)

    def add_tool_bars(self):
        """
            creates tool bars in the main window on the top underneath menu bars.
        """
        
        self.main_tool_bar = QToolBar()
        self.addToolBar(self.main_tool_bar)
        self.main_tool_bar.setMovable(True)

        self.create_shelves()
        self.main_tool_bar.addWidget(self.main_shelf_w)

    def initialize_main_group_boxes(self) -> None:
        """
            creates main group box widgets within tabs
        """
        self.create_reset_button()
        self.create_fluid_solver_widget()
        self.create_appearance_widgets()
        self.create_motion_widget()
        self.create_distribution_widgets()
        self.create_physical_widget()
        self.create_fluid_tank_widgets()
        self.create_system_widgets()
        self.main_layout.addSpacerItem(QSpacerItem(0,250))

# ------------------------------------------------------------- CREATE SHELVES ----------------------------------------------------------------------------
    def create_shelves(self):
        """
            widget to create shelf bars on top
        """
        self.main_shelf_w = QTabWidget()
        self.main_shelf_w.setMovable(True)
        self.main_shelf_w.setTabShape(QTabWidget.TabShape.Rounded)

        self.shelf_v_layout = QVBoxLayout()

        self.create_tank_shelf()
        self.create_colliders_shelf()

        self.main_shelf_w.setLayout(self.shelf_v_layout)

        self.main_shelf_w.addTab(self.tank_shelf, "Tank")
        self.main_shelf_w.addTab(self.colliders_shelf, "Colliders")

    def create_tank_shelf(self):
        """
            create collider shapes shelf tool
        """   
        
        self.tank_shelf = QWidget()
        self.tank_h_layout = QHBoxLayout()

        self.tank_shelf.setLayout(self.tank_h_layout)

        self.rectangular_collider = QToolButton()
        self.rectangular_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.rectangular_collider,
                                  sub_menu="Shelf", image_name="cube.png")
        self.sphere_collider = QToolButton()
        self.sphere_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.sphere_collider,
                                  sub_menu="Shelf", image_name="sphere.png")
        self.capsule_collider = QToolButton()
        self.capsule_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.capsule_collider,
                                  sub_menu="Shelf", image_name="capsule.png")
        self.cylinder_collider = QToolButton()
        self.cylinder_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.cylinder_collider,
                                  sub_menu="Shelf", image_name="cylinder.png")

        self.tank_h_layout.addWidget(self.rectangular_collider)
        self.tank_h_layout.addWidget(self.sphere_collider)
        self.tank_h_layout.addWidget(self.capsule_collider)
        self.tank_h_layout.addWidget(self.cylinder_collider)
        self.tank_h_layout.addStretch(200)

    def create_colliders_shelf(self):
        """
            create tank types shelf tool
        """
        self.colliders_shelf = QWidget()
        self.colliders_h_layout = QHBoxLayout()

        self.colliders_shelf.setLayout(self.colliders_h_layout)

        self.ellipse_collider = QToolButton()
        self.ellipse_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.ellipse_collider,
                                  sub_menu="Shelf", image_name="ellipse.png")
        self.rectangle_collider = QToolButton()
        self.rectangle_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.rectangle_collider,
                                  sub_menu="Shelf", image_name="rectangle.png")
        self.polygon_collider = QToolButton()
        self.polygon_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.polygon_collider,
                                  sub_menu="Shelf", image_name="polygon.png")
        self.triangle_collider = QToolButton()
        self.triangle_collider.setMinimumSize(QSize(40,40))
        self.add_icons_to_widgets(tool_button=self.triangle_collider,
                                  sub_menu="Shelf", image_name="triangle.png")

        self.colliders_h_layout.addWidget(self.ellipse_collider)
        self.colliders_h_layout.addWidget(self.rectangle_collider)
        self.colliders_h_layout.addWidget(self.polygon_collider)
        self.colliders_h_layout.addWidget(self.triangle_collider)
        self.colliders_h_layout.addStretch(200)

# -------------------------------------------------------------- CREATE DRAWING AREAS ---------------------------------------------------------------------
    def create_main_canvas(self):
        """
            creates and attaches the bottom toolbar to the main window
            widget
        """
        self.graphics_view_v_layout = QVBoxLayout()
        self.main_canvas = RenderScene(transparent=False,
                                       system_obj=self.system_obj)

    def create_frame_controls(self):
        """
            create frame toolbar at bottom of window
        """
        self.play_bar_h_layout = QHBoxLayout()
        self.play_simulation_btn = CustomPushButton()
        self.add_icons_to_widgets(widget=self.play_simulation_btn, image_name="play_backward.png")
        self.stop_simulation_btn = CustomPushButton()
        self.add_icons_to_widgets(widget=self.stop_simulation_btn, image_name="stop.png")
        self.play_back_simulation_btn = CustomPushButton()
        self.add_icons_to_widgets(widget=self.play_back_simulation_btn, image_name="play_forward.png")

        self.frame_spacing_item = QSpacerItem(700, 0)
        self.prev_frame_btn = CustomPushButton()
        self.add_icons_to_widgets(widget=self.prev_frame_btn, image_name="frame_backward.png")
        self.curr_frame_lbl = CustomLineEdit()
        self.curr_frame_lbl.setSizePolicy(QSizePolicy(QSizePolicy.Policy.MinimumExpanding,
                                                        QSizePolicy.Policy.Fixed))
        self.curr_frame_lbl.setFixedSize(QSize(40,20))
        self.next_frame_btn = CustomPushButton()
        self.add_icons_to_widgets(widget=self.next_frame_btn, image_name="frame_forward.png")

        self.frame_range_h_layout = QHBoxLayout()
        self.start_frame_field = CustomLineEdit(text="0")
        self.start_frame_field.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                         QSizePolicy.Policy.Fixed))
        self.start_frame_field.setFixedSize(QSize(40,20))
        self.end_frame_field = CustomLineEdit(text="100")
        self.end_frame_field.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                        QSizePolicy.Policy.Fixed))
        self.end_frame_field.setFixedSize(QSize(40,20))

        self.frame_range_h_layout.addStretch(800)
        self.frame_range_h_layout.addWidget(self.start_frame_field)
        self.frame_range_h_layout.addWidget(self.end_frame_field)

        # adding widgets to play bar after frame range hierarchies to have 
        # access to the frame range value class instances

        self.play_bar_h_layout.addWidget(self.play_simulation_btn)
        self.play_bar_h_layout.addWidget(self.stop_simulation_btn)
        self.play_bar_h_layout.addWidget(self.play_back_simulation_btn)

        self.create_frame_bars(layout=self.play_bar_h_layout,
                               start = 0, end = 1000,
                               line_edit_top=self.start_frame_field,
                               line_edit_side=self.end_frame_field)

        self.play_bar_h_layout.addWidget(self.prev_frame_btn)
        self.play_bar_h_layout.addWidget(self.curr_frame_lbl)
        self.play_bar_h_layout.addWidget(self.next_frame_btn)

        self.frame_control_widget = QWidget()
        self.frame_control_v_layout = QVBoxLayout()
        self.frame_control_widget.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                            QSizePolicy.Policy.Fixed))
        self.frame_control_v_layout.addLayout(self.frame_range_h_layout)
        self.frame_control_v_layout.addLayout(self.play_bar_h_layout)
        self.frame_control_widget.setLayout(self.frame_control_v_layout)

        self.animation_controls_dock_w = QDockWidget("Frame Controls")
        self.animation_controls_dock_w.setWidget(self.frame_control_widget)

        self.animation_controls_dock_w.setSizePolicy(QSizePolicy.Policy.Minimum,
                                                     QSizePolicy.Policy.Minimum)

        self.create_bottom_splitter()

    def create_bottom_splitter(self):
        """
            frame splitter and dock widget creation
        """
        self.frame_and_graphics_splitter = QSplitter(Qt.Orientation.Vertical)
        self.frame_and_graphics_splitter.setChildrenCollapsible(False)

        self.frame_and_graphics_splitter.addWidget(self.main_canvas)
        self.frame_and_graphics_splitter.addWidget(self.animation_controls_dock_w)
   
    def create_reset_button(self):
        """
            creates a reset to defaults button in main menu widget
        """
        self.reset_btn = CustomPushButton(title="Reset to Defaults", 
                                          size_policy=QSizePolicy(QSizePolicy.Policy.Expanding,
                                                                QSizePolicy.Policy.Fixed))
        self.main_layout.addWidget(self.reset_btn)

    def create_fluid_solver_widget(self):
        """
            creates the widget for fluid solvers
        """
        self.solver_vertical_layout = QVBoxLayout()
        self.solver_horizontal_layout = QHBoxLayout()

        self.fluid_solver_group_box = CustomGroupBox(title="Solver Type")
        self.set_fixed_size_policy(self.fluid_solver_group_box)

        self.fluid_solver_label = CustomLabel(title="Choose Solver")
        fluid_solver_lbls = ["SPH", "Multi SPH", "WCSPH", "PCSPH", "IISPH", "Eulerian", "FLIP", "PIC"]
        
        self.fluid_solver_combo_box = CustomComboBox(add_items=True, items_to_add=fluid_solver_lbls)

        self.solver_horizontal_layout.addWidget(self.fluid_solver_label)
        self.solver_horizontal_layout.addWidget(self.fluid_solver_combo_box)

        self.fluid_solver_group_box.setLayout(self.solver_vertical_layout)
        self.solver_vertical_layout.addLayout(self.solver_horizontal_layout)
        self.main_layout.addWidget(self.fluid_solver_group_box)

# -------------------- PARTICLE WIDGETS -------------------------------------

    def create_appearance_widgets(self):
        """
            creates widgets responsible for fluid appearance
        """
        self.appearance_gBox = CustomGroupBox(title="Appearance",
                                              fixed_size_policy=True,
                                              checkable=True)
        self.appearance_v_layout = QVBoxLayout()

        self.nbr_particles_h_layout = QHBoxLayout()
        self.nbr_of_particles_lbl = CustomLabel(title="No. of particles",
                                                label_control=True,
                                                increment=1)
        self.nbr_of_particle_sBox = CustomSpinBox(maximum=150000, minimum=100, 
                                                  increment=1, value=100000, 
                                                  alias="particle number", 
                                                  system_obj=self.system_obj)
        self.nbr_of_particles_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                                    linked_spin_box=self.nbr_of_particle_sBox,
                                                    minimum=100, maximum=150000, alias="particle number",
                                                    system_obj=self.system_obj)
        self.nbr_of_particle_sBox.slider = self.nbr_of_particles_slider

        self.nbr_of_particles_lbl.spin = self.nbr_of_particle_sBox
        self.nbr_of_particles_lbl.slider = self.nbr_of_particles_slider

        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_lbl)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particle_sBox)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_slider)

        self.particle_size_h_layout = QHBoxLayout()
        self.particle_size_lbl = CustomLabel(title="Particle Size",
                                             label_control=True,
                                             increment=1)
        self.particle_size_sBox = CustomDoubleSpinBox(maximum=10,
                                                    minimum=0,
                                                    increment=0.1, value=0.1)
        self.particle_size_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                                   linked_double_spin_box=self.particle_size_sBox,
                                                   minimum=0, maximum=10)
        self.particle_size_sBox.slider = self.particle_size_slider_w

        self.particle_size_lbl.spin = self.particle_size_sBox
        self.particle_size_lbl.slider = self.particle_size_slider_w
        
        self.particle_size_h_layout.addWidget(self.particle_size_lbl)
        self.particle_size_h_layout.addWidget(self.particle_size_sBox)
        self.particle_size_h_layout.addWidget(self.particle_size_slider_w)

        self.appearance_gBox.setLayout(self.appearance_v_layout)
        self.appearance_v_layout.addLayout(self.nbr_particles_h_layout)
        self.appearance_v_layout.addLayout(self.particle_size_h_layout)
        self.set_default_state(self.appearance_gBox)
        self.adjust_label_spacing(box_widget=self.appearance_gBox,
                                  scale_factor=8)

        self.main_layout.addWidget(self.appearance_gBox)

    def create_motion_widget(self):
        """
            creates the widget for the motion of the particles
        """
        self.motion_v_layout = QVBoxLayout()
        self.motion_gBox = CustomGroupBox(title="Motion",
                                          fixed_size_policy=True,
                                          checkable=True)


        self.velocity_h_layout = QHBoxLayout()
        self.velocity_lbl = CustomLabel(title="Velocity",
                                        label_control=True,
                                        increment = 0.1)
        self.velocity_x_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                                   value=0.5, alias="velocity", direction=0,
                                                   system_obj=self.system_obj)
        self.velocity_y_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                                   value=0.2, alias="velocity", direction=1,
                                                   system_obj=self.system_obj)
        self.velocity_z_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                                   value=0.7, alias="velocity", direction=2,
                                                   system_obj=self.system_obj)

        self.velocity_lbl.double_spin = self.velocity_x_sBox
        self.velocity_lbl.double_spin = self.velocity_y_sBox
        self.velocity_lbl.double_spin = self.velocity_z_sBox

        self.velocity_h_layout.addWidget(self.velocity_lbl)
        self.velocity_h_layout.addWidget(self.velocity_x_sBox)
        self.velocity_h_layout.addWidget(self.velocity_y_sBox)
        self.velocity_h_layout.addWidget(self.velocity_z_sBox)

        self.acc_h_layout = QHBoxLayout()
        self.acc_lbl = CustomLabel(title="Acceleration",
                                   label_control=True,
                                   increment = 0.1)
        self.acc_x_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                              value=0.2, alias="accaeleration", direction=0,
                                              system_obj=self.system_obj)
        self.acc_y_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                              value=0.4, alias="acceleration", direction=1,
                                              system_obj=self.system_obj)
        self.acc_z_sBox = CustomDoubleSpinBox(maximum=20, minimum=-20, increment=0.1,
                                              value=0.7, alias="acceleration", direction=2,
                                              system_obj=self.system_obj)

        self.acc_lbl.double_spin = self.acc_x_sBox
        self.acc_lbl.double_spin = self.acc_y_sBox
        self.acc_lbl.double_spin = self.acc_z_sBox

        self.acc_h_layout.addWidget(self.acc_lbl)
        self.acc_h_layout.addWidget(self.acc_x_sBox)
        self.acc_h_layout.addWidget(self.acc_y_sBox)
        self.acc_h_layout.addWidget(self.acc_z_sBox)

        self.motion_gBox.setLayout(self.motion_v_layout)
        self.motion_v_layout.addLayout(self.velocity_h_layout)
        self.motion_v_layout.addLayout(self.acc_h_layout)
        self.set_default_state(self.motion_gBox)
        self.adjust_label_spacing(box_widget=self.motion_gBox,
                                  scale_factor=8)

        self.main_layout.addWidget(self.motion_gBox)

    def create_physical_widget(self):
        """
            creates physical forces widget appearance
        """
        self.physical_v_layout = QVBoxLayout()
        self.physical_gBox = CustomGroupBox(title="Physical Forces",
                                            fixed_size_policy=True,
                                            checkable=True)

        self.gravity_h_layout = QHBoxLayout()
        self.physical_g_lbl = CustomLabel(title="Gravity",
                                          label_control=True,
                                          increment=0.1)
        self.physical_g_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                   value=-9.81, multiplier=10, alias="gravity",
                                                   system_obj=self.system_obj)
        self.physical_g_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                              linked_double_spin_box=self.physical_g_sBox,
                                              minimum=-300, maximum=300, multiplier=0.1,
                                              alias = "gravity", system_obj=self.system_obj)

        self.gravity_h_layout.addWidget(self.physical_g_lbl)
        self.gravity_h_layout.addWidget(self.physical_g_sBox)
        self.gravity_h_layout.addWidget(self.physical_g_slider)

        self.physical_g_lbl.double_spin = self.physical_g_sBox
        self.physical_g_lbl.slider = self.physical_g_slider
        self.physical_g_sBox.slider = self.physical_g_slider

        self.buoyancy_h_layout = QHBoxLayout()
        self.physical_b_lbl = CustomLabel(title="Buoyancy",
                                          label_control=True,
                                          increment=0.1)
        self.physical_b_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                   value=0, multiplier=10, alias="buoyancy",
                                                   system_obj=self.system_obj)
        self.physical_b_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                              linked_double_spin_box=self.physical_b_sBox,
                                              minimum=-300, maximum=300, multiplier=0.1, 
                                              alias = "buoyancy", system_obj=self.system_obj)
        
        self.buoyancy_h_layout.addWidget(self.physical_b_lbl)
        self.buoyancy_h_layout.addWidget(self.physical_b_sBox)
        self.buoyancy_h_layout.addWidget(self.physical_b_slider)

        self.physical_b_lbl.double_spin = self.physical_b_sBox
        self.physical_b_lbl.slider = self.physical_b_slider
        self.physical_b_sBox.slider = self.physical_b_slider
        
        self.viscosity_h_layout = QHBoxLayout()
        self.physical_v_lbl = CustomLabel(title="Viscosity",
                                          label_control=True,
                                          increment=0.1)
        self.physical_v_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                   value=3.5, multiplier=10, alias="viscosity",
                                                   system_obj=self.system_obj)
        self.physical_v_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                              linked_double_spin_box=self.physical_v_sBox,
                                              minimum=-300, maximum=300, multiplier=0.1,
                                              alias = "viscosity", system_obj=self.system_obj)

        self.viscosity_h_layout.addWidget(self.physical_v_lbl)
        self.viscosity_h_layout.addWidget(self.physical_v_sBox)
        self.viscosity_h_layout.addWidget(self.physical_v_slider)

        self.physical_v_lbl.double_spin = self.physical_v_sBox
        self.physical_v_lbl.slider = self.physical_v_slider
        self.physical_v_sBox.slider = self.physical_v_slider

        self.pressure_h_layout = QHBoxLayout()
        self.physical_p_lbl = CustomLabel(title="Pressure",
                                          label_control=True,
                                          increment=0.1)
        self.physical_p_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                   value=5, multiplier=10, alias="pressure",
                                                   system_obj=self.system_obj)
        self.physical_p_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                              linked_spin_box=self.physical_p_sBox,
                                              minimum=-300, maximum=300, multiplier=0.1,
                                              alias="pressure", system_obj=self.system_obj)

        self.pressure_h_layout.addWidget(self.physical_p_lbl)
        self.pressure_h_layout.addWidget(self.physical_p_sBox)
        self.pressure_h_layout.addWidget(self.physical_p_slider)
        
        self.physical_p_lbl.double_spin = self.physical_p_sBox
        self.physical_p_lbl.slider = self.physical_p_slider
        self.physical_p_sBox.slider = self.physical_p_slider
        
        self.mass_h_layout = QHBoxLayout()
        self.physical_m_lbl = CustomLabel(title="Mass",
                                          label_control=True,
                                          increment=0.1)
        self.physical_m_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                   value=0.1, multiplier=10, alias="mass",
                                                   system_obj=self.system_obj)
        self.physical_m_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                              linked_double_spin_box=self.physical_m_sBox,
                                              minimum=-300, maximum=300, multiplier=0.1,
                                              alias="mass", system_obj=self.system_obj)

        self.mass_h_layout.addWidget(self.physical_m_lbl)
        self.mass_h_layout.addWidget(self.physical_m_sBox)
        self.mass_h_layout.addWidget(self.physical_m_slider)

        self.physical_m_lbl.double_spin = self.physical_m_sBox
        self.physical_m_lbl.slider = self.physical_m_slider
        self.physical_m_sBox.slider = self.physical_m_slider
        
        self.massD_h_layout = QHBoxLayout()
        self.physical_md_lbl = CustomLabel(title="Mass Density",
                                          label_control=True,
                                          increment=0.1)
        self.physical_md_sBox = CustomDoubleSpinBox(minimum=0, maximum=2000, increment=1,
                                                    value=998.2, multiplier=1, alias="mass density",
                                                    system_obj=self.system_obj)
        self.physical_md_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                               linked_double_spin_box=self.physical_md_sBox,
                                               minimum=0, maximum=2000, multiplier=1,
                                               alias="mass density", system_obj=self.system_obj)

        self.massD_h_layout.addWidget(self.physical_md_lbl)
        self.massD_h_layout.addWidget(self.physical_md_sBox)
        self.massD_h_layout.addWidget(self.physical_md_slider)
        
        self.physical_md_lbl.double_spin = self.physical_md_sBox
        self.physical_md_lbl.slider = self.physical_md_slider
        self.physical_md_sBox.slider = self.physical_md_slider

        self.speedLoss_h_layout = QHBoxLayout()
        self.physical_sL_lbl = CustomLabel(title="Speed Loss",
                                          label_control=True,
                                          increment=0.1)
        self.physical_sL_sBox = CustomDoubleSpinBox(minimum=-30, maximum=30, increment=0.1,
                                                    value=0.1, multiplier=10, alias="speed loss",
                                                    system_obj=self.system_obj)
        self.physical_sL_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                               linked_double_spin_box=self.physical_sL_sBox,
                                               minimum=-300, maximum=300, multiplier=0.1,
                                               alias="speed loss", system_obj=self.system_obj)

        self.speedLoss_h_layout.addWidget(self.physical_sL_lbl)
        self.speedLoss_h_layout.addWidget(self.physical_sL_sBox)
        self.speedLoss_h_layout.addWidget(self.physical_sL_slider)

        self.physical_sL_lbl.double_spin = self.physical_sL_sBox
        self.physical_sL_lbl.slider = self.physical_sL_slider
        self.physical_sL_sBox.slider = self.physical_sL_slider

        self.add_multiple_layouts(self.physical_gBox, self.physical_v_layout,
                                [self.mass_h_layout, self.gravity_h_layout,
                                self.buoyancy_h_layout, self.viscosity_h_layout, 
                                self.pressure_h_layout, self.massD_h_layout, 
                                self.speedLoss_h_layout])

        self.set_default_state(self.physical_gBox)
        self.adjust_label_spacing(box_widget=self.physical_gBox,
                                  scale_factor=8)

        self.main_layout.addWidget(self.physical_gBox)

# ----------------------- TANK WIDGET ---------------------------------------

    def create_fluid_tank_widgets(self):
        """
            creates the widgets for the fluid tank session of the main window UI
        """
        self.tank_v_layout = QVBoxLayout()
        self.tank_gBox = CustomGroupBox(title="Tank Control",
                                        fixed_size_policy=True,
                                        checkable=True)

        self.tank_type_h_layout = QHBoxLayout()
        self.tank_type_lbl = CustomLabel(title="Tank Type")
        tank_items = ["Circular", "Cubic", "Cylindrical", "Capsule"]
        self.tank_combo_box = CustomComboBox(add_items=True,
                                             items_to_add=tank_items)

        self.tank_type_h_layout.addWidget(self.tank_type_lbl)
        self.tank_type_h_layout.addWidget(self.tank_combo_box)

        self.tank_radius_h_layout = QHBoxLayout()
        self.tank_radius_lbl = CustomLabel(title="Tank Radius",
                                          label_control=True,
                                          increment=0.1)
        self.tank_x_radius_sBox = CustomDoubleSpinBox(minimum=1, maximum=10, increment=0.1,
                                                      alias = "tank radius", direction=0,
                                                      system_obj=self.system_obj)
        self.tank_y_radius_sBox = CustomDoubleSpinBox(minimum=1, maximum=10, increment=0.1,
                                                      alias = "tank radius", direction=1,
                                                      system_obj=self.system_obj)
        self.tank_z_radius_sBox = CustomDoubleSpinBox(minimum=1, maximum=10, increment=0.1,
                                                      alias = "tank radius", direction=2,
                                                      system_obj=self.system_obj)

        self.tank_radius_lbl.double_spin = self.tank_x_radius_sBox
        self.tank_radius_lbl.double_spin = self.tank_y_radius_sBox
        self.tank_radius_lbl.double_spin = self.tank_z_radius_sBox

        self.tank_radius_h_layout.addWidget(self.tank_radius_lbl)
        self.tank_radius_h_layout.addWidget(self.tank_x_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_y_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_z_radius_sBox)

        self.tank_pos_h_layout = QHBoxLayout()
        self.tank_pos_lbl = CustomLabel(title="Tank Position",
                                          label_control=True,
                                          increment=0.1)
        self.tank_x_pos_sBox = CustomDoubleSpinBox(minimum=-10, maximum=10, increment=0.1,
                                                      alias = "tank position", direction=0,
                                                      system_obj=self.system_obj)
        self.tank_y_pos_sBox = CustomDoubleSpinBox(minimum=-10, maximum=10, increment=0.1,
                                                      alias = "tank position", direction=1,
                                                      system_obj=self.system_obj)
        self.tank_z_pos_sBox = CustomDoubleSpinBox(minimum=-10, maximum=10, increment=0.1,
                                                      alias = "tank position", direction=2,
                                                      system_obj=self.system_obj)

        self.tank_pos_lbl.double_spin = self.tank_x_pos_sBox
        self.tank_pos_lbl.double_spin = self.tank_y_pos_sBox
        self.tank_pos_lbl.double_spin = self.tank_z_pos_sBox

        self.tank_pos_h_layout.addWidget(self.tank_pos_lbl)
        self.tank_pos_h_layout.addWidget(self.tank_x_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_y_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_z_pos_sBox)

        self.tank_gBox.setLayout(self.tank_v_layout)
        self.tank_v_layout.addLayout(self.tank_radius_h_layout)
        self.tank_v_layout.addLayout(self.tank_pos_h_layout)
        self.set_default_state(self.tank_gBox)
        self.adjust_label_spacing(box_widget=self.tank_gBox,
                                  scale_factor=8)

        self.main_layout.addWidget(self.tank_gBox)

# -------------------- MAIN SYSTEM WIDGETS ----------------------------------

    def create_system_widgets(self):
        """
            responsible for the system of particle controls
        """
        self.timeSteps_v_layout = QVBoxLayout()
        self.timeSteps_gBox = CustomGroupBox(title="Simulation Steps",
                                             fixed_size_policy=True,
                                             checkable=True)

        self.substep_type_h_layout = QHBoxLayout()
        self.substep_type_lbl = CustomLabel(title="Time Integrator")
        substep_type = ["Default", "Verlet", "Leap_Frog"]
        self.substep_comboBox = CustomComboBox(add_items=True,
                                               items_to_add=substep_type)

        self.substep_type_h_layout.addWidget(self.substep_type_lbl)
        self.substep_type_h_layout.addWidget(self.substep_comboBox)

        self.delta_t_h_layout = QHBoxLayout()
        self.delta_t_lbl = CustomLabel(title="Delta Time",
                                          label_control=True,
                                          increment=0.01)
        self.delta_t_slider = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                           minimum=0, maximum=100, multiplier=0.01,
                                           alias = "delta time", system_obj=self.system_obj)
        self.delta_t_sBox = CustomDoubleSpinBox(minimum=0, maximum=0.2, increment=0.01,
                                                alias = "delta time", system_obj=self.system_obj,
                                                value=0.01, multiplier=100)

        self.delta_t_h_layout.addWidget(self.delta_t_lbl)
        self.delta_t_h_layout.addWidget(self.delta_t_slider)
        self.delta_t_h_layout.addWidget(self.delta_t_sBox)

        self.delta_t_lbl.slider = self.delta_t_slider
        self.delta_t_lbl.double_spin = self.delta_t_sBox

        self.timeSteps_gBox.setLayout(self.timeSteps_v_layout)
        self.timeSteps_v_layout.addLayout(self.substep_type_h_layout)
        self.timeSteps_v_layout.addLayout(self.delta_t_h_layout)
        self.set_default_state(self.timeSteps_gBox)

        self.main_layout.addWidget(self.timeSteps_gBox)
    
    def create_distribution_widgets(self):
        """
            controls widgets for neighbouring between particle solver types
            as well as parameters to control particle separation and grid
            scale
        """
        self.distr_v_layout = QVBoxLayout()
        self.distr_gBox = CustomGroupBox(title="Particle Distribution",
                                         fixed_size_policy=True,
                                         checkable=True,
                                         default_state=True)

        self.particle_distr_h_layout = QHBoxLayout()
        self.particle_distr_lbl = CustomLabel(title="Particle Distribution")
        particle_distr_items = ["Uniform", "Random"]
        self.particle_distr_comboBox = CustomComboBox(add_items=True,
                                                      items_to_add=particle_distr_items)

        self.particle_distr_h_layout.addWidget(self.particle_distr_lbl)
        self.particle_distr_h_layout.addWidget(self.particle_distr_comboBox)

        self.neighbr_solver_h_layout = QHBoxLayout()
        neighbr_solver_items = ["Distance", "K d Trees", "Octree",
                                "Quadtree", "Spatial Hashing"]
        self.neighbr_solver_lbl = CustomLabel(title="Neighbour Solver")
        self.neighbr_solver_comboBox = CustomComboBox(add_items=True,
                                                      items_to_add=neighbr_solver_items)

        self.neighbr_solver_h_layout.addWidget(self.neighbr_solver_lbl)
        self.neighbr_solver_h_layout.addWidget(self.neighbr_solver_comboBox)

        self.particle_sep_h_layout = QHBoxLayout()
        self.particle_sep_lbl = CustomLabel(title="Particle Separation",
                                          label_control=True,
                                          increment=0.01)
        self.particle_sep_spinBox = CustomDoubleSpinBox(minimum=0.0001, maximum=0.1, increment=0.01,
                                                        multiplier=1000, alias = "particle separation",
                                                        system_obj=self.system_obj)
        self.particle_sep_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                                  linked_double_spin_box=self.particle_sep_spinBox,
                                                  multiplier=0.0001, alias = "particle separation", 
                                                  system_obj=self.system_obj)
        self.particle_sep_spinBox.slider = self.particle_sep_slider_w

        self.particle_sep_h_layout.addWidget(self.particle_sep_lbl)
        self.particle_sep_h_layout.addWidget(self.particle_sep_spinBox)
        self.particle_sep_h_layout.addWidget(self.particle_sep_slider_w)

        self.particle_sep_lbl.double_spin = self.particle_sep_spinBox
        self.particle_sep_lbl.slider = self.particle_sep_slider_w

        self.cell_size_h_layout = QHBoxLayout()
        self.cell_size_lbl = CustomLabel(title="Cell Size",
                                          label_control=True,
                                          increment=0.01)
        self.cell_size_spinBox = CustomDoubleSpinBox(minimum=0.001, maximum=0.5, increment=0.01,
                                                     multiplier=1000, alias="cell size",
                                                     system_obj=self.system_obj)
        self.cell_size_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal,
                                               linked_double_spin_box=self.cell_size_spinBox,
                                               multiplier=0.01, alias="cell size",
                                               system_obj=self.system_obj)
        self.cell_size_spinBox.slider = self.cell_size_slider_w

        self.cell_size_lbl.double_spin = self.cell_size_spinBox
        self.cell_size_lbl.slider = self.cell_size_slider_w

        self.cell_size_h_layout.addWidget(self.cell_size_lbl)
        self.cell_size_h_layout.addWidget(self.cell_size_spinBox)
        self.cell_size_h_layout.addWidget(self.cell_size_slider_w)

        self.distr_gBox.setLayout(self.distr_v_layout)
        self.distr_v_layout.addLayout(self.particle_distr_h_layout)
        self.distr_v_layout.addLayout(self.neighbr_solver_h_layout)
        self.distr_v_layout.addLayout(self.particle_sep_h_layout)
        self.distr_v_layout.addLayout(self.cell_size_h_layout)
        self.set_default_state(self.distr_gBox)
        self.adjust_label_spacing(box_widget=self.distr_gBox,
                                  scale_factor=6)

        self.main_layout.addWidget(self.distr_gBox)
#-----------------------------------------------------------------------------

    def add_menu_bars(self):
        """
            responsible for menu bar widgets and dockers in central main window
        """
        menu_bar = self.menuBar()
        file_bar = menu_bar.addMenu("File")

        new_action = file_bar.addAction("New")
        open_action = file_bar.addAction("Open")
        open_recent = file_bar.addMenu("Open Recent...")

        file_bar.addSeparator()

        set_project = file_bar.addAction("Set Project")

        file_bar.addSeparator()

        save_action = file_bar.addAction("Save")
        save_as_action = file_bar.addAction("Save As")
        import_action = file_bar.addAction("Import")
        export_action = file_bar.addAction("Export")

        file_bar.addSeparator()
        quit_action = file_bar.addAction("Quit")

        # creating the edit toolbar
        edit_bar = menu_bar.addMenu("Edit")

        copy_action = edit_bar.addAction("Copy")
        paste_action = edit_bar.addAction("Paste")

        edit_bar.addSeparator()

        undo_action = edit_bar.addAction("Undo")
        redo_action = edit_bar.addAction("Redo")

        edit_bar.addSeparator()

        delete_action = edit_bar.addAction("Delete")

        window_bar = menu_bar.addMenu("Window")

        full_screen = window_bar.addAction("FullScreen")

        window_bar.addSeparator()

        resizable = window_bar.addAction("Resizable")
        resizable.setCheckable(True)

        expand_window = window_bar.addAction("Expand")
        shrink_window = window_bar.addAction("Shrink")

        display_bar = menu_bar.addMenu("Display")

        display_all = display_bar.addAction("Display All")
        display_all.setCheckable(True)

        hide_menu = display_bar.addMenu("Hide")
        
        hide_tank_display = hide_menu.addAction("Tank")
        hide_tank_display.setCheckable(True)
        hide_colliders_display = hide_menu.addAction("Collider")
        hide_colliders_display.setCheckable(True)
        hide_fps_control = hide_menu.addAction("FPS")
        hide_fps_control.setCheckable(True)
        hide_display_particle_count = hide_menu.addAction("No. of particles")
        hide_display_particle_count.setCheckable(True)

        display_bar.addSeparator()

        tank_display = display_bar.addAction("Tank")
        tank_display.setCheckable(True)
        colliders_display = display_bar.addAction("Collider")
        colliders_display.setCheckable(True)
        fps_control = display_bar.addAction("FPS")
        fps_control.setCheckable(True)
        display_particle_count = display_bar.addAction("No. of particles")
        display_particle_count.setCheckable(True)
        reset_simulation = display_bar.addAction("Reset Simulation")

        display_bar.addSeparator()

        dockers_options = display_bar.addMenu("Dockers")
        fluid_options = dockers_options.addAction("Fluid Options")
        fluid_options.setCheckable(True)
        fluid_pred_options = dockers_options.addAction("Fluid Prediction")
        fluid_pred_options.setCheckable(True)
        frame_controls = dockers_options.addAction("Frame Controls")
        frame_controls.setCheckable(True)
        viewport_tools = dockers_options.addAction("Tools")
        viewport_tools.setCheckable(True)

        all_viewport_tools = display_bar.addMenu("Tools")

        select_action = all_viewport_tools.addAction("Select")
        translate_action = all_viewport_tools.addAction("Translate")
        rotate_action = all_viewport_tools.addAction("Rotate")
        scale_action = all_viewport_tools.addAction("Scale")

        shelves = display_bar.addMenu("Shelves")

        tank_shelf = shelves.addAction("Tank")
        collider_shelf = shelves.addAction("Collider")

        predictor_bar = menu_bar.addMenu("Predictor")

        #creating the settings bar
        settings_bar = menu_bar.addMenu("Settings")
        restore_values = settings_bar.addAction("Restore Session")
        reset_values = settings_bar.addAction("Reset Values")

        settings_bar.addSeparator()

        more_settings = settings_bar.addAction("More Settings...")

        information_bar = menu_bar.addMenu("About")

        about_app = information_bar.addAction("About Predictor")
        learn_more = information_bar.addAction("Learn More")
        help = information_bar.addAction("Help")

    def fluid_pred_w(self):
        """
            creates fluid prediction controls
        """
        
        self.fluid_pred_v_layout = QVBoxLayout()

        self.fluid_pred_options_widget.setLayout(self.fluid_pred_v_layout)

        self.attrib_pred_h_layout = QHBoxLayout()
        
        self.attribs_pred_lbl = CustomLabel(title="Attrs to Predict")
        attribs_pred_list = ["Default", "Position", "Velocity", "Colour", "Pressure", "Mass Density"]
        self.attribs_pred_combo_box = CustomComboBox(add_items=True, items_to_add=attribs_pred_list)

        self.attrib_pred_h_layout.addWidget(self.attribs_pred_lbl)
        self.attrib_pred_h_layout.addWidget(self.attribs_pred_combo_box)

        self.pred_model_h_layout = QHBoxLayout()

        self.pred_model_lbl = CustomLabel(title="Choose Model")
        pred_model_list = ["Linear", "Recmoid", "Laplacian"]
        self.pred_model_combo_box = CustomComboBox(add_items=True, items_to_add=pred_model_list)

        self.pred_model_h_layout.addWidget(self.pred_model_lbl)
        self.pred_model_h_layout.addWidget(self.pred_model_combo_box)

        self.epoch_h_layout = QHBoxLayout()

        self.num_epochs_lbl = CustomLabel(title="No. Epochs")
        self.num_epochs_spin_box = CustomSpinBox()
        self.num_epochs_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.epoch_h_layout.addWidget(self.num_epochs_lbl)
        self.epoch_h_layout.addWidget(self.num_epochs_spin_box)
        self.epoch_h_layout.addWidget(self.num_epochs_slider)

        self.export_loc_h_layout = QHBoxLayout()

        self.export_file_path = CustomLabel(title="File Path")
        self.export_line_edit = CustomLineEdit(text="Export location specify here...")

        self.export_loc_h_layout.addWidget(self.export_file_path)
        self.export_loc_h_layout.addWidget(self.export_line_edit)

        self.reset_default_btn = QPushButton("Reset")

        self.fluid_pred_v_layout.addWidget(self.reset_default_btn)
        self.fluid_pred_v_layout.addLayout(self.attrib_pred_h_layout)
        self.fluid_pred_v_layout.addLayout(self.pred_model_h_layout)
        self.fluid_pred_v_layout.addLayout(self.epoch_h_layout)
        self.fluid_pred_v_layout.addLayout(self.export_loc_h_layout)
        self.fluid_pred_v_layout.addStretch(100)

        self.adjust_label_spacing_layout(box_widget=self.fluid_pred_options_widget,
                                         scale_factor=8)
