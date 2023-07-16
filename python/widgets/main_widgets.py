from PySide6.QtWidgets import QVBoxLayout, QWidget, QMainWindow, QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox,  \
    QGroupBox, QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap,  \
    QSpacerItem, QSizePolicy,  QScrollArea, QSplitter, QDockWidget
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon, QImage, QColor, QFont
from canvas_state import DrawingCanvas
from utility_functions import *
from widget_utilities import *

class MainWindow(QMainWindow, UtilFuncs):

    def __init__(self, app):
        super().__init__()
        self.app = app
        
        self.setGeometry(QRect(100, 50, 1280, 720))
        self.setWindowTitle("Fluid Application")
        self.window_icon = QPixmap("images/Toolbar/fluid.png")
        self.setWindowIcon(self.window_icon)

        self.initUI()

        self.setCentralWidget(self.main_splitter_widget)

    def initUI(self):
        """
            initializes calls to most UI elements on screen
        """
        self.add_menu_bars()
        self.add_tool_bars()

        self.main_tabs_widget = QTabWidget()

        self.main_options_widget = QWidget()
        self.main_options_widget.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Expanding,
                                                           QSizePolicy.Policy.Fixed))
        self.main_layout = QVBoxLayout()

        self.fluid_pred_options_widget = QWidget()
        self.tab_v_layout = QVBoxLayout()

        self.create_scroll_area_widget()

        self.initialize_main_group_boxes()

        self.horizontal_spacer_layout = QHBoxLayout()

        self.main_options_widget.setLayout(self.main_layout)
        self.create_bottom_toolbar()

        self.main_tabs_widget.setMinimumSize(QSize(225,400))
        self.main_tabs_widget.setMaximumSize(QSize(325,720))

        self.main_splitter_widget = QSplitter(Qt.Horizontal)
        self.main_splitter_widget.addWidget(self.main_graphics_widget)
        self.main_splitter_widget.addWidget(self.main_tabs_widget)
        self.main_splitter_widget.setChildrenCollapsible(False)

    def create_scroll_area_widget(self):
        """
            creates and configures the initial scroll area widget binded to each
            tab
        """
        self.fluid_sim_scroll_area = QScrollArea()
        self.fluid_sim_scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.fluid_sim_scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.fluid_sim_scroll_area.setWidgetResizable(True)

        self.fluid_sim_scroll_area.setWidget(self.main_options_widget)
        self.main_tabs_widget.addTab(self.fluid_sim_scroll_area, "Fluid Simulation")
        self.main_tabs_widget.addTab(self.fluid_pred_options_widget, "Fluid Prediction")
        self.main_tabs_widget.setLayout(self.tab_v_layout)

    def add_tool_bars(self):
        """
            creates tool bars in the main window on the top underneath menu bars.
        """
        
        self.main_tool_bar = QToolBar()
        self.addToolBar(self.main_tool_bar)
        self.main_tool_bar.setMovable(True)

        self.undo_tool_action = self.main_tool_bar.addAction("Undo")
        self.redo_tool_action = self.main_tool_bar.addAction("Redo")

        self.add_icons_to_actions(icon=self.undo_tool_action, image_name="undo.png")
        self.add_icons_to_actions(icon=self.redo_tool_action, image_name="redo.png")

        self.main_tool_bar.addSeparator()

        self.save_tool_action = self.main_tool_bar.addAction("Save")
        self.copy_tool_action = self.main_tool_bar.addAction("Copy")
        self.paste_tool_action = self.main_tool_bar.addAction("Paste")

        self.add_icons_to_actions(icon=self.save_tool_action, image_name="save.png")
        self.add_icons_to_actions(icon=self.copy_tool_action, image_name="copy.png")
        self.add_icons_to_actions(icon=self.paste_tool_action, image_name="paste.png")

        self.main_tool_bar.addSeparator()

        self.settings_tool_action = self.main_tool_bar.addAction("Settings")
        self.add_icons_to_actions(icon=self.settings_tool_action, image_name="settings.png")

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
    
    def create_bottom_toolbar(self):
        """
            creates and attaches the bottom toolbar to the main window
            widget
        """
        self.graphics_view_v_layout = QVBoxLayout()
        self.main_canvas = DrawingCanvas(QRect(0, 0, 700, 400), 
                                         QColor(255,0,0), 
                                         QFont("Arial"), 14)

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

        self.play_bar_h_layout.addWidget(self.play_simulation_btn)
        self.play_bar_h_layout.addWidget(self.stop_simulation_btn)
        self.play_bar_h_layout.addWidget(self.play_back_simulation_btn)

        self.play_bar_h_layout.addSpacerItem(self.frame_spacing_item)

        self.play_bar_h_layout.addWidget(self.prev_frame_btn)
        self.play_bar_h_layout.addWidget(self.curr_frame_lbl)
        self.play_bar_h_layout.addWidget(self.next_frame_btn)

        self.frame_range_h_layout = QHBoxLayout()
        self.start_frame_field = CustomLineEdit()
        self.start_frame_field.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                         QSizePolicy.Policy.Fixed))
        self.start_frame_field.setFixedSize(QSize(40,20))
        self.end_frame_field = CustomLineEdit()
        self.end_frame_field.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                        QSizePolicy.Policy.Fixed))
        self.end_frame_field.setFixedSize(QSize(40,20))

        self.frame_range_h_layout.addSpacerItem(QSpacerItem(850,0))
        self.frame_range_h_layout.addWidget(self.start_frame_field)
        self.frame_range_h_layout.addWidget(self.end_frame_field)

        self.frame_control_widget = QWidget()
        self.frame_control_v_layout = QVBoxLayout()
        self.frame_control_widget.setSizePolicy(QSizePolicy(QSizePolicy.Policy.Fixed,
                                                            QSizePolicy.Policy.Fixed))
        self.frame_control_v_layout.addLayout(self.frame_range_h_layout)
        self.frame_control_v_layout.addLayout(self.play_bar_h_layout)
        self.frame_control_widget.setLayout(self.frame_control_v_layout)

        self.frame_and_graphics_splitter = QSplitter(Qt.Vertical)
        self.frame_and_graphics_splitter.addWidget(self.main_canvas)
        self.frame_and_graphics_splitter.addWidget(self.frame_control_widget)

        self.graphics_view_v_layout.addWidget(self.frame_and_graphics_splitter)

        self.main_graphics_widget = QWidget()
        self.main_graphics_widget.setLayout(self.graphics_view_v_layout)
   
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
        fluid_solver_lbls = ["SPH", "IISPH", "WCSPH", "Multi SPH"]
        
        self.fluid_solver_combo_box = CustomComboBox(add_items=True, items_to_add=fluid_solver_lbls)

        self.solver_horizontal_layout.addWidget(self.fluid_solver_label)
        self.solver_horizontal_layout.addWidget(self.fluid_solver_combo_box)

        self.fluid_solver_group_box.setLayout(self.solver_vertical_layout)
        self.solver_vertical_layout.addLayout(self.solver_horizontal_layout)
        self.main_layout.addWidget(self.fluid_solver_group_box)

# -------------------- PARTICLE WIDGETS ------------------------------

    def create_appearance_widgets(self):
        """
            creates widgets responsible for fluid appearance
        """
        self.appearance_gBox = CustomGroupBox(title="Appearance",
                                              fixed_size_policy=True,
                                              checkable=True)
        self.appearance_v_layout = QVBoxLayout()

        self.nbr_particles_h_layout = QHBoxLayout()
        self.nbr_of_particles_lbl = CustomLabel(title="No. of particles")
        self.nbr_of_particle_sBox = CustomSpinBox()
        self.nbr_of_particles_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_lbl)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particle_sBox)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_slider)

        self.particle_size_h_layout = QHBoxLayout()
        self.particle_size_lbl = CustomLabel(title="Particle Size")
        self.particle_size_sBox = CustomSpinBox()
        self.particle_size_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal)
        
        self.particle_size_h_layout.addWidget(self.particle_size_lbl)
        self.particle_size_h_layout.addWidget(self.particle_size_sBox)
        self.particle_size_h_layout.addWidget(self.particle_size_slider_w)

        self.appearance_gBox.setLayout(self.appearance_v_layout)
        self.appearance_v_layout.addLayout(self.nbr_particles_h_layout)
        self.appearance_v_layout.addLayout(self.particle_size_h_layout)
        self.set_default_state(self.appearance_gBox)

        self.main_layout.addWidget(self.appearance_gBox)

    def create_motion_widget(self):
        """
            creates the widget for the motion of the particles
        """
        self.motion_v_layout = QVBoxLayout()
        self.motion_gBox = CustomGroupBox(title="Motion",
                                          fixed_size_policy=True,
                                          checkable=True)

        self.position_h_layout = QHBoxLayout()
        self.position_lbl = CustomLabel(title="Position")
        self.position_x_sBox = CustomDoubleSpinBox()
        self.position_y_sBox = CustomDoubleSpinBox()
        self.position_z_sBox = CustomDoubleSpinBox()

        self.position_h_layout.addWidget(self.position_lbl)
        self.position_h_layout.addWidget(self.position_x_sBox)
        self.position_h_layout.addWidget(self.position_y_sBox)
        self.position_h_layout.addWidget(self.position_z_sBox)

        self.velocity_h_layout = QHBoxLayout()
        self.velocity_lbl = CustomLabel(title="Velocity")
        self.velocity_x_sBox = CustomDoubleSpinBox()
        self.velocity_y_sBox = CustomDoubleSpinBox()
        self.velocity_z_sBox = CustomDoubleSpinBox()

        self.velocity_h_layout.addWidget(self.velocity_lbl)
        self.velocity_h_layout.addWidget(self.velocity_x_sBox)
        self.velocity_h_layout.addWidget(self.velocity_y_sBox)
        self.velocity_h_layout.addWidget(self.velocity_z_sBox)

        self.acc_h_layout = QHBoxLayout()
        self.acc_lbl = CustomLabel(title="Acceleration")
        self.acc_x_sBox = CustomDoubleSpinBox()
        self.acc_y_sBox = CustomDoubleSpinBox()
        self.acc_z_sBox = CustomDoubleSpinBox()

        self.acc_h_layout.addWidget(self.acc_lbl)
        self.acc_h_layout.addWidget(self.acc_x_sBox)
        self.acc_h_layout.addWidget(self.acc_y_sBox)
        self.acc_h_layout.addWidget(self.acc_z_sBox)

        self.motion_gBox.setLayout(self.motion_v_layout)
        self.motion_v_layout.addLayout(self.position_h_layout)
        self.motion_v_layout.addLayout(self.velocity_h_layout)
        self.motion_v_layout.addLayout(self.acc_h_layout)
        self.set_default_state(self.motion_gBox)

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
        self.physical_g_lbl = CustomLabel(title="Gravity")
        self.physical_g_sBox = CustomDoubleSpinBox()
        self.physical_g_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.gravity_h_layout.addWidget(self.physical_g_lbl)
        self.gravity_h_layout.addWidget(self.physical_g_sBox)
        self.gravity_h_layout.addWidget(self.physical_g_slider)

        self.buoyancy_h_layout = QHBoxLayout()
        self.physical_b_lbl = CustomLabel(title="Buoyancy")
        self.physical_b_sBox = CustomDoubleSpinBox()
        self.physical_b_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.buoyancy_h_layout.addWidget(self.physical_b_lbl)
        self.buoyancy_h_layout.addWidget(self.physical_b_sBox)
        self.buoyancy_h_layout.addWidget(self.physical_b_slider)

        self.viscosity_h_layout = QHBoxLayout()
        self.physical_v_lbl = CustomLabel(title="Viscosity")
        self.physical_v_sBox = CustomDoubleSpinBox()
        self.physical_v_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.viscosity_h_layout.addWidget(self.physical_v_lbl)
        self.viscosity_h_layout.addWidget(self.physical_v_sBox)
        self.viscosity_h_layout.addWidget(self.physical_v_slider)

        self.pressure_h_layout = QHBoxLayout()
        self.physical_p_lbl = CustomLabel(title="Pressure")
        self.physical_p_sBox = CustomDoubleSpinBox()
        self.physical_p_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.pressure_h_layout.addWidget(self.physical_p_lbl)
        self.pressure_h_layout.addWidget(self.physical_p_sBox)
        self.pressure_h_layout.addWidget(self.physical_p_slider)

        self.mass_h_layout = QHBoxLayout()
        self.physical_m_lbl = CustomLabel(title="Mass")
        self.physical_m_sBox = CustomDoubleSpinBox()
        self.physical_m_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.mass_h_layout.addWidget(self.physical_m_lbl)
        self.mass_h_layout.addWidget(self.physical_m_sBox)
        self.mass_h_layout.addWidget(self.physical_m_slider)

        self.massD_h_layout = QHBoxLayout()
        self.physical_md_lbl = CustomLabel(title="Mass Density")
        self.physical_md_sBox = CustomDoubleSpinBox()
        self.physical_md_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.massD_h_layout.addWidget(self.physical_md_lbl)
        self.massD_h_layout.addWidget(self.physical_md_sBox)
        self.massD_h_layout.addWidget(self.physical_md_slider)

        self.speedLoss_h_layout = QHBoxLayout()
        self.physical_sL_lbl = CustomLabel(title="Speed Loss")
        self.physical_sL_sBox = CustomDoubleSpinBox()
        self.physical_sL_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.speedLoss_h_layout.addWidget(self.physical_sL_lbl)
        self.speedLoss_h_layout.addWidget(self.physical_sL_sBox)
        self.speedLoss_h_layout.addWidget(self.physical_sL_slider)

        self.physical_gBox.setLayout(self.physical_v_layout)
        self.physical_v_layout.addLayout(self.mass_h_layout)
        self.physical_v_layout.addLayout(self.gravity_h_layout)
        self.physical_v_layout.addLayout(self.buoyancy_h_layout)
        self.physical_v_layout.addLayout(self.viscosity_h_layout)
        self.physical_v_layout.addLayout(self.pressure_h_layout)
        self.physical_v_layout.addLayout(self.massD_h_layout)
        self.physical_v_layout.addLayout(self.speedLoss_h_layout)
        self.set_default_state(self.physical_gBox)

        self.main_layout.addWidget(self.physical_gBox)

# ----------------------- TANK WIDGET -------------------------------

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
        self.tank_radius_lbl = CustomLabel(title="Tank Radius")
        self.tank_x_radius_sBox = CustomDoubleSpinBox()
        self.tank_y_radius_sBox = CustomDoubleSpinBox()
        self.tank_z_radius_sBox = CustomDoubleSpinBox()

        self.tank_radius_h_layout.addWidget(self.tank_radius_lbl)
        self.tank_radius_h_layout.addWidget(self.tank_x_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_y_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_z_radius_sBox)

        self.tank_pos_h_layout = QHBoxLayout()
        self.tank_pos_lbl = CustomLabel(title="Tank Position")
        self.tank_x_pos_sBox = CustomDoubleSpinBox()
        self.tank_y_pos_sBox = CustomDoubleSpinBox()
        self.tank_z_pos_sBox = CustomDoubleSpinBox()

        self.tank_pos_h_layout.addWidget(self.tank_pos_lbl)
        self.tank_pos_h_layout.addWidget(self.tank_x_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_y_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_z_pos_sBox)

        self.tank_gBox.setLayout(self.tank_v_layout)
        self.tank_v_layout.addLayout(self.tank_radius_h_layout)
        self.tank_v_layout.addLayout(self.tank_pos_h_layout)
        self.set_default_state(self.tank_gBox)

        self.main_layout.addWidget(self.tank_gBox)


# -------------------- MAIN SYSTEM WIDGETS ---------------------------------

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
        self.delta_t_lbl = CustomLabel(title="Delta Time")
        self.delta_t_slider = CustomSlider(orientation=Qt.Orientation.Horizontal)
        self.delta_t_sBox = CustomDoubleSpinBox()

        self.delta_t_h_layout.addWidget(self.delta_t_lbl)
        self.delta_t_h_layout.addWidget(self.delta_t_slider)
        self.delta_t_h_layout.addWidget(self.delta_t_sBox)

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
        self.particle_sep_lbl = CustomLabel(title="Particle Separation")
        self.particle_sep_spinBox = CustomDoubleSpinBox()
        self.particle_sep_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.particle_sep_h_layout.addWidget(self.particle_sep_lbl)
        self.particle_sep_h_layout.addWidget(self.particle_sep_spinBox)
        self.particle_sep_h_layout.addWidget(self.particle_sep_slider_w)

        self.cell_size_h_layout = QHBoxLayout()
        self.cell_size_lbl = CustomLabel(title="Cell Size")
        self.cell_size_spinBox = CustomDoubleSpinBox()
        self.cell_size_slider_w = CustomSlider(orientation=Qt.Orientation.Horizontal)

        self.cell_size_h_layout.addWidget(self.cell_size_lbl)
        self.cell_size_h_layout.addWidget(self.cell_size_spinBox)
        self.cell_size_h_layout.addWidget(self.cell_size_slider_w)

        self.distr_gBox.setLayout(self.distr_v_layout)
        self.distr_v_layout.addLayout(self.particle_distr_h_layout)
        self.distr_v_layout.addLayout(self.neighbr_solver_h_layout)
        self.distr_v_layout.addLayout(self.particle_sep_h_layout)
        self.distr_v_layout.addLayout(self.cell_size_h_layout)
        self.set_default_state(self.distr_gBox)

        self.main_layout.addWidget(self.distr_gBox)
#-----------------------------------------------------------------------------

    def add_menu_bars(self):
        """
            responsible for menu bar widgets and dockers in central main window
        """
        menu_bar = self.menuBar()
        file_bar = menu_bar.addMenu("File")

        quit_action = file_bar.addAction("Quit")
        undo_action = file_bar.addAction("Undo")
        redo_action = file_bar.addAction("Redo")

        self.add_icons_to_actions(quit_action, image_name="quit.png")
        self.add_icons_to_actions(undo_action, image_name="undo.png")
        self.add_icons_to_actions(redo_action, image_name="redo.png")

        # creating the edit toolbar
        edit_bar = menu_bar.addMenu("Edit")

        save_action = edit_bar.addAction("Save")
        copy_action = edit_bar.addAction("Copy")
        paste_action = edit_bar.addAction("Paste")
        
        self.add_icons_to_actions(icon=save_action, image_name="save.png")
        self.add_icons_to_actions(icon=copy_action, image_name="copy.png")
        self.add_icons_to_actions(icon=paste_action, image_name="paste.png")

        #creating the import/export menu bars
        import_bar = menu_bar.addMenu("Import/Export")

        json_import = import_bar.addAction("Import JSON")
        txt_import = import_bar.addAction("Import txt")
        import_bar.addSeparator()
        json_export = import_bar.addAction("Export JSON")
        txt_export = import_bar.addAction("Export txt")

        self.add_icons_to_actions(icon=json_import, image_name="import.png")
        self.add_icons_to_actions(icon=txt_import, image_name="import.png")
        self.add_icons_to_actions(icon=json_export, image_name="export.png")
        self.add_icons_to_actions(icon=txt_export, image_name="export.png")

        #creating the settings bar
        settings_bar = menu_bar.addMenu("Settings")

        fps_control = settings_bar.addAction("Display FPS")
        display_particle_count = settings_bar.addAction("Display No. of particles")
        reset_simulation = settings_bar.addAction("Reset Simulation")

        self.add_icons_to_actions(icon=fps_control, image_name="fps.png")
        self.add_icons_to_actions(icon=display_particle_count, image_name="count.png")
        self.add_icons_to_actions(icon=reset_simulation, image_name="reset.png")