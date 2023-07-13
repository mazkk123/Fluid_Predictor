from PySide6.QtWidgets import QVBoxLayout, QWidget, QMainWindow, QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox, QGroupBox, \
    QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap, QRadioButton, QSpacerItem
from PySide6.QtCore import QSize, QRect , Qt
from PySide6.QtGui import QPixmap, QIcon

class MainWindow(QMainWindow):

    def __init__(self, app):
        super().__init__()
        self.app = app

        self.setGeometry(QRect(0, 0, 1280, 720))
        self.setWindowTitle("Fluid Application")
        self.window_icon = QPixmap("images/Toolbar/fluid.png")
        self.setWindowIcon(self.window_icon)

        self.add_menu_bars()
        self.add_tool_bars()

        self.main_tabs_widget = QTabWidget()

        self.main_options_widget = QWidget()
        self.main_layout = QVBoxLayout()

        self.fluid_pred_options_widget = QWidget()
        self.tab_v_layout = QVBoxLayout()

        self.main_tabs_widget.addTab(self.main_options_widget, "Fluid Simulation")
        self.main_tabs_widget.addTab(self.fluid_pred_options_widget, "Fluid Prediction")
        self.main_tabs_widget.setLayout(self.tab_v_layout)

        # creates main group box widgets
        self.create_reset_button()
        self.create_fluid_solver_widget()
        self.create_appearance_widgets()
        """ self.create_motion_widget()
        self.create_distribution_widgets()
        self.create_physical_widget() """
        self.create_fluid_tank_widgets()
        self.create_system_widgets()

        self.horizontal_spacer_layout = QHBoxLayout()

        self.main_options_widget.setLayout(self.main_layout)
        self.horizontal_spacer_layout.addSpacerItem(QSpacerItem(1000, 400))
        self.horizontal_spacer_layout.addWidget(self.main_tabs_widget)

        self.main_window = QWidget()
        self.main_window.setLayout(self.horizontal_spacer_layout)

        self.setCentralWidget(self.main_window)

    def add_tool_bars(self):
        """
            creates tool bars in the main window on the top underneath menu bars.
        """
        
        self.main_tool_bar = QToolBar()
        self.addToolBar(self.main_tool_bar)
        self.main_tool_bar.setMovable(False)

        self.quit_tool_action = self.main_tool_bar.addAction("Quit")
        self.undo_tool_action = self.main_tool_bar.addAction("Undo")
        self.redo_tool_action = self.main_tool_bar.addAction("Redo")

        self.add_icons_to_actions(icon=self.quit_tool_action, image_name="quit.png")
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

    def create_reset_button(self):
        """
            creates a reset to defaults button in main menu widget
        """
        self.reset_btn = QPushButton("Reset to Defaults")
        self.main_layout.addWidget(self.reset_btn)

    def create_fluid_solver_widget(self):
        """
            creates the widget for fluid solvers
        """
        self.solver_vertical_layout = QVBoxLayout()
        self.solver_horizontal_layout = QHBoxLayout()

        self.fluid_solver_group_box = QGroupBox("Solver Type")
        self.fluid_solver_group_box.setCheckable(True)

        self.fluid_solver_label = QLabel("Choose Solver")
        self.fluid_solver_combo_box = QComboBox()
        self.fluid_solver_combo_box.addItem("SPH")
        self.fluid_solver_combo_box.addItem("IISPH")
        self.fluid_solver_combo_box.addItem("WCSPH")
        self.fluid_solver_combo_box.addItem("Multi SPH")

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
        self.appearance_group_box = QGroupBox("Appearance")
        self.appearance_group_box.setCheckable(True)
        self.appearance_v_layout = QVBoxLayout()

        self.nbr_particles_h_layout = QHBoxLayout()
        self.nbr_of_particles_lbl = QLabel("No. of particles")
        self.nbr_of_particle_sBox = QSpinBox()
        self.nbr_of_particles_slider = QSlider(Qt.Horizontal)

        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_lbl)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particle_sBox)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_slider)

        self.particle_size_h_layout = QHBoxLayout()
        self.particle_size_lbl = QLabel("Particle Size")
        self.particle_size_sBox = QDoubleSpinBox()
        self.particle_size_slider = QSlider(Qt.Horizontal)
        
        self.particle_size_h_layout.addWidget(self.particle_size_lbl)
        self.particle_size_h_layout.addWidget(self.particle_size_sBox)
        self.particle_size_h_layout.addWidget(self.particle_size_slider)

        self.appearance_group_box.setLayout(self.appearance_v_layout)
        self.appearance_v_layout.addLayout(self.nbr_particles_h_layout)
        self.appearance_v_layout.addLayout(self.particle_size_h_layout)

        self.main_layout.addWidget(self.appearance_group_box)

    def create_motion_widget(self):
        """
            creates the widget for the motion of the particles
        """
        self.motion_v_layout = QVBoxLayout()
        self.motion_gBox = QGroupBox("Motion")
        self.motion_gBox.setCheckable(True)

        self.position_h_layout = QHBoxLayout()
        self.position_lbl = QLabel("Position")
        self.position_x_sBox = QDoubleSpinBox()
        self.position_y_sBox = QDoubleSpinBox()
        self.position_z_sBox = QDoubleSpinBox()

        self.position_h_layout.addWidget(self.position_lbl)
        self.position_h_layout.addWidget(self.position_x_sBox)
        self.position_h_layout.addWidget(self.position_y_sBox)
        self.position_h_layout.addWidget(self.position_z_sBox)

        self.velocity_h_layout = QHBoxLayout()
        self.velocity_lbl = QLabel("Velocity")
        self.velocity_x_sBox = QDoubleSpinBox()
        self.velocity_y_sBox = QDoubleSpinBox()
        self.velocity_z_sBox = QDoubleSpinBox()

        self.velocity_h_layout.addWidget(self.velocity_lbl)
        self.velocity_h_layout.addWidget(self.velocity_x_sBox)
        self.velocity_h_layout.addWidget(self.velocity_y_sBox)
        self.velocity_h_layout.addWidget(self.velocity_z_sBox)

        self.acc_h_layout = QHBoxLayout()
        self.acc_lbl = QLabel("Acceleration")
        self.acc_x_sBox = QDoubleSpinBox()
        self.acc_y_sBox = QDoubleSpinBox()
        self.acc_z_sBox = QDoubleSpinBox()

        self.acc_h_layout.addWidget(self.acc_lbl)
        self.acc_h_layout.addWidget(self.acc_x_sBox)
        self.acc_h_layout.addWidget(self.acc_y_sBox)
        self.acc_h_layout.addWidget(self.acc_z_sBox)

        self.motion_gBox.setLayout(self.motion_v_layout)
        self.motion_v_layout.addLayout(self.position_h_layout)
        self.motion_v_layout.addLayout(self.velocity_h_layout)
        self.motion_v_layout.addLayout(self.acc_h_layout)

        self.main_layout.addWidget(self.motion_gBox)

    def create_physical_widget(self):
        """
            creates physical forces widget appearance
        """
        self.physical_v_layout = QVBoxLayout()
        self.physical_gBox = QGroupBox("Physical Forces")
        self.physical_gBox.setCheckable(True)

        self.gravity_h_layout = QHBoxLayout()
        self.physical_g_lbl = QLabel("Gravity")
        self.physical_g_sBox = QDoubleSpinBox()
        self.physical_g_slider = QSlider(Qt.Horizontal)

        self.gravity_h_layout.addWidget(self.physical_g_lbl)
        self.gravity_h_layout.addWidget(self.physical_g_sBox)
        self.gravity_h_layout.addWidget(self.physical_g_slider)

        self.buoyancy_h_layout = QHBoxLayout()
        self.physical_b_lbl = QLabel("Buoyancy")
        self.physical_b_sBox = QDoubleSpinBox()
        self.physical_b_slider = QSlider(Qt.Horizontal)

        self.buoyancy_h_layout.addWidget(self.physical_b_lbl)
        self.buoyancy_h_layout.addWidget(self.physical_b_sBox)
        self.buoyancy_h_layout.addWidget(self.physical_b_slider)

        self.viscosity_h_layout = QHBoxLayout()
        self.physical_v_lbl = QLabel("Viscosity")
        self.physical_v_sBox = QDoubleSpinBox()
        self.physical_v_slider = QSlider(Qt.Horizontal)

        self.viscosity_h_layout.addWidget(self.physical_v_lbl)
        self.viscosity_h_layout.addWidget(self.physical_v_sBox)
        self.viscosity_h_layout.addWidget(self.physical_v_slider)

        self.pressure_h_layout = QHBoxLayout()
        self.physical_p_lbl = QLabel("Pressure")
        self.physical_p_sBox = QDoubleSpinBox()
        self.physical_p_slider = QSlider(Qt.Horizontal)

        self.pressure_h_layout.addWidget(self.physical_p_lbl)
        self.pressure_h_layout.addWidget(self.physical_p_sBox)
        self.pressure_h_layout.addWidget(self.physical_p_slider)

        self.mass_h_layout = QHBoxLayout()
        self.physical_m_lbl = QLabel("Mass")
        self.physical_m_sBox = QDoubleSpinBox()
        self.physical_m_slider = QSlider(Qt.Horizontal)

        self.mass_h_layout.addWidget(self.physical_m_lbl)
        self.mass_h_layout.addWidget(self.physical_m_sBox)
        self.mass_h_layout.addWidget(self.physical_m_slider)

        self.massD_h_layout = QHBoxLayout()
        self.physical_md_lbl = QLabel("Mass Density")
        self.physical_md_sBox = QDoubleSpinBox()
        self.physical_md_slider = QSlider(Qt.Horizontal)

        self.massD_h_layout.addWidget(self.physical_md_lbl)
        self.massD_h_layout.addWidget(self.physical_md_sBox)
        self.massD_h_layout.addWidget(self.physical_md_slider)

        self.speedLoss_h_layout = QHBoxLayout()
        self.physical_sL_lbl = QLabel("Speed Loss")
        self.physical_sL_sBox = QDoubleSpinBox()
        self.physical_sL_slider = QSlider(Qt.Horizontal)

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

        self.main_layout.addWidget(self.physical_gBox)

# ----------------------- TANK WIDGET -------------------------------

    def create_fluid_tank_widgets(self):
        """
            creates the widgets for the fluid tank session of the main window UI
        """
        self.tank_v_layout = QVBoxLayout()
        self.tank_gBox = QGroupBox("Tank Control")
        self.tank_gBox.setCheckable(True)

        self.tank_type_h_layout = QHBoxLayout()
        self.tank_type_lbl = QLabel("Tank Type")
        self.tank_combo_box = QComboBox()

        self.tank_combo_box.addItem("Circular")
        self.tank_combo_box.addItem("Cubic")
        self.tank_combo_box.addItem("Cylindrical")
        self.tank_combo_box.addItem("Capsule")

        self.tank_type_h_layout.addWidget(self.tank_type_lbl)
        self.tank_type_h_layout.addWidget(self.tank_combo_box)

        self.tank_radius_h_layout = QHBoxLayout()
        self.tank_radius_lbl = QLabel("Tank Radius")
        self.tank_x_radius_sBox = QDoubleSpinBox()
        self.tank_y_radius_sBox = QDoubleSpinBox()
        self.tank_z_radius_sBox = QDoubleSpinBox()

        self.tank_radius_h_layout.addWidget(self.tank_radius_lbl)
        self.tank_radius_h_layout.addWidget(self.tank_x_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_y_radius_sBox)
        self.tank_radius_h_layout.addWidget(self.tank_z_radius_sBox)

        self.tank_pos_h_layout = QHBoxLayout()
        self.tank_pos_lbl = QLabel("Tank Position")
        self.tank_x_pos_sBox = QDoubleSpinBox()
        self.tank_y_pos_sBox = QDoubleSpinBox()
        self.tank_z_pos_sBox = QDoubleSpinBox()

        self.tank_pos_h_layout.addWidget(self.tank_pos_lbl)
        self.tank_pos_h_layout.addWidget(self.tank_x_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_y_pos_sBox)
        self.tank_pos_h_layout.addWidget(self.tank_z_pos_sBox)

        self.tank_gBox.setLayout(self.tank_v_layout)
        self.tank_v_layout.addLayout(self.tank_radius_h_layout)
        self.tank_v_layout.addLayout(self.tank_pos_h_layout)

        self.main_layout.addWidget(self.tank_gBox)


# -------------------- MAIN SYSTEM WIDGETS ---------------------------------

    def create_system_widgets(self):
        """
            responsible for the system of particle controls
        """
        self.timeSteps_v_layout = QVBoxLayout()
        self.timeSteps_gBox = QGroupBox("Simulation Steps")
        self.timeSteps_gBox.setCheckable(True)

        self.substep_type_h_layout = QHBoxLayout()
        self.substep_type_lbl = QLabel("Time Integrator")
        self.substep_comboBox = QComboBox()

        self.substep_comboBox.addItem("Default")
        self.substep_comboBox.addItem("Verlet")
        self.substep_comboBox.addItem("Leap_Frog")

        self.substep_type_h_layout.addWidget(self.substep_type_lbl)
        self.substep_type_h_layout.addWidget(self.substep_comboBox)

        self.delta_t_h_layout = QHBoxLayout()
        self.delta_t_lbl = QLabel("Delta Time")
        self.delta_t_slider = QSlider(Qt.Horizontal)
        self.delta_t_sBox = QDoubleSpinBox()

        self.delta_t_h_layout.addWidget(self.delta_t_lbl)
        self.delta_t_h_layout.addWidget(self.delta_t_slider)
        self.delta_t_h_layout.addWidget(self.delta_t_sBox)

        self.timeSteps_gBox.setLayout(self.timeSteps_v_layout)
        self.timeSteps_v_layout.addLayout(self.substep_type_h_layout)
        self.timeSteps_v_layout.addLayout(self.delta_t_h_layout)

        self.main_layout.addWidget(self.timeSteps_gBox)
    
    def create_distribution_widgets(self):
        """
            controls widgets for neighbouring between particle solver types
            as well as parameters to control particle separation and grid
            scale
        """
        self.distr_v_layout = QVBoxLayout()
        self.distr_gBox = QGroupBox("Particle Distribution")
        self.distr_gBox.setCheckable(True)

        self.particle_distr_h_layout = QHBoxLayout()
        self.particle_distr_lbl = QLabel("Particle Distribution")
        self.particle_distr_comboBox = QComboBox()

        self.particle_distr_comboBox.addItem("Uniform")
        self.particle_distr_comboBox.addItem("Random")

        self.particle_distr_h_layout.addWidget(self.particle_distr_lbl)
        self.particle_distr_h_layout.addWidget(self.particle_distr_comboBox)

        self.neighbr_solver_h_layout = QHBoxLayout()
        self.neighbr_solver_lbl = QLabel("Neighbour Solver")
        self.neighbr_solver_comboBox = QComboBox()

        self.neighbr_solver_comboBox.addItem("Distance")
        self.neighbr_solver_comboBox.addItem("K d trees")
        self.neighbr_solver_comboBox.addItem("Octree")
        self.neighbr_solver_comboBox.addItem("Quadtree")
        self.neighbr_solver_comboBox.addItem("Spatial Hashing")

        self.neighbr_solver_h_layout.addWidget(self.neighbr_solver_lbl)
        self.neighbr_solver_h_layout.addWidget(self.neighbr_solver_comboBox)

        self.particle_sep_h_layout = QHBoxLayout()
        self.particle_sep_lbl = QLabel("Particle Separation")
        self.particle_sep_spinBox = QDoubleSpinBox()
        self.particle_sep_slider_w = QSlider(orientation=Qt.Horizontal)

        self.particle_sep_h_layout.addWidget(self.particle_sep_lbl)
        self.particle_sep_h_layout.addWidget(self.particle_sep_spinBox)
        self.particle_sep_h_layout.addWidget(self.particle_sep_slider_w)

        self.cell_size_h_layout = QHBoxLayout()
        self.cell_size_lbl = QLabel("Cell Size")
        self.cell_size_spinBox = QDoubleSpinBox()
        self.cell_size_slider_w = QSlider(orientation=Qt.Horizontal)

        self.cell_size_h_layout.addWidget(self.cell_size_lbl)
        self.cell_size_h_layout.addWidget(self.cell_size_spinBox)
        self.cell_size_h_layout.addWidget(self.cell_size_slider_w)

        self.distr_gBox.setLayout(self.distr_v_layout)
        self.distr_v_layout.addLayout(self.particle_distr_h_layout)
        self.distr_v_layout.addLayout(self.neighbr_solver_h_layout)
        self.distr_v_layout.addLayout(self.particle_sep_h_layout)
        self.distr_v_layout.addLayout(self.cell_size_h_layout)

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
    
    def add_icons_to_actions(self, icon=None, sub_menu="Toolbar", image_name=None):
        """
            adds icons based on their name to local list
        """
        if icon is not None:
            if image_name is not None:
                path = "images/" + sub_menu + "/" + str(image_name)
                icon.setIcon(QIcon(path))

# ------------------------------ MAIN SLOT ACTIVATION ------------------------------

    def activate_all_slots(self):
        """
            activates all slot widgets
        """
        pass


    def activate_reset_slot(self):
        """
            activates the reset slot calls
        """
        self.reset_btn.clicked.connect(self.reset_btn_func)

    def activate_appearance_slots(self):
        """
            activates appearance slot callbacks
        """
        self.nbr_of_particles_slider.valueChanged.connect(self.particle_num_slider)
        self.nbr_of_particle_sBox.valueChanged.connect(self.particle_num_spin_box)

        self.particle_size_sBox.valueChanged.connect(self.particle_size_spin_box)

    def activate_motion_slots(self):
        """
            activates motion slot callbacks
        """
        self.position_x_sBox.valueChanged.connect(self.mPosition_spin_box)
        self.position_y_sBox.valueChanged.connect(self.mPosition_spin_box)
        self.position_z_sBox.valueChanged.connect(self.mPosition_spin_box)

        self.velocity_x_sBox.valueChanged.connect(self.mVel_spin_box)
        self.velocity_y_sBox.valueChanged.connect(self.mVel_spin_box)
        self.velocity_z_sBox.valueChanged.connect(self.mVel_spin_box)

        self.acc_x_sBox.valueChanged.connect(self.mAccel_spin_box)
        self.acc_y_sBox.valueChanged.connect(self.mAccel_spin_box)
        self.acc_z_sBox.valueChanged.connect(self.mAccel_spin_box)

    def activate_physical_slots(self):
        """
            activates physical slot callbacks
        """
        self.physical_m_slider.valueChanged.connect(self.mass_slider)
        self.physical_m_sBox.valueChanged.connect(self.mass_spin_box)

        self.physical_g_sBox.valueChanged.connect(self.gravity_spin_box)
        self.physical_g_slider.valueChanged.connect(self.gravity_slider)

        self.physical_b_sBox.valueChanged.connect(self.buoyancy_spin_box)
        self.physical_b_slider.valueChanged.connect(self.buoyancy_slider)

        self.physical_v_sBox.valueChanged.connect(self.viscosity_spin_box)
        self.physical_v_slider.valueChanged.connect(self.viscosity_slider)

        self.physical_p_sBox.valueChanged.connect(self.pressure_spin_box)
        self.physical_p_slider.valueChanged.connect(self.pressure_slider)

        self.physical_md_sBox.valueChanged.connect(self.massD_spin_box)
        self.physical_md_slider.valueChanged.connect(self.massD_slider)

        self.physical_sL_sBox.valueChanged.connect(self.speed_loss_spin_box)
        self.physical_sL_slider.valueChanged.connect(self.speed_loss_slider)

    def activate_distr_slots(self):
        """
            activates distribution slot callbacks
        """
        self.particle_distr_comboBox.currentIndexChanged.connect(self.distr_combo_box)
        self.neighbr_solver_comboBox.currentIndexChanged.connect(self.solver_combo_box)

        self.particle_sep_spinBox.valueChanged.connect(self.particle_sep_spin_box)
        self.particle_sep_slider_w.valueChanged.connect(self.particle_sep_slider)

        self.cell_size_spinBox.valueChanged.connect(self.cell_size_spin_box)
        self.cell_size_slider_w.valueChanged.connect(self.cell_size_slider)

    def activate_tank_slots(self):
        """
            activates tank slot callbacks
        """
        self.tank_x_radius_sBox.valueChanged.connect(self.tank_radius_spin_box)
        self.tank_y_radius_sBox.valueChanged.connect(self.tank_radius_spin_box)
        self.tank_z_radius_sBox.valueChanged.connect(self.tank_radius_spin_box)

        self.tank_combo_box.currentIndexChanged.connect(self.tank_type_combo_box)
        
        self.tank_x_pos_sBox.valueChanged.connect(self.tank_position_spin_box)
        self.tank_y_pos_sBox.valueChanged.connect(self.tank_position_spin_box)
        self.tank_z_pos_sBox.valueChanged.connect(self.tank_position_spin_box)

    def activate_time_slots(self):
        """
            activates time 
        """
        self.delta_t_sBox.valueChanged.connect(self.delta_time_spin_box)
        self.delta_t_slider.valueChanged.connect(self.delta_time_slider)
        
        self.substep_comboBox.currentIndexChanged.connect(self.time_integrator_combo_box)

        
# ------------------------------- GROUPBOX SLOTS ------------------------------------

    def appearance_gBox(self):
        """
            sets appearance slot functions
        """
        pass

    def motion_gBox(self):
        """
            sets motion slot functions
        """
        pass

    def distr_gBox(self):
        """
            sets distribution slot functions
        """
        pass

    def physical_gBox(self):
        """
            sets physical slot functions
        """
        pass

    def tank_gBox(self):
        """
            sets tank slot functions
        """
        pass

    def time_gBox(self):
        """
            sets time variation slot functions
        """
        pass

# ---------------------------------- SLIDER SLOTS -----------------------------------------

    def particle_num_slider(self):
        """
            slider template code
        """
        pass

    def particle_size_slider(self):
        """
            slider template code
        """
        pass

    def particle_sep_slider(self):
        """
            slider template code
        """
        pass
    
    def cell_size_slider(self):
        """
            slider template code
        """
        pass

    def mass_slider(self):
        """
            slider template code
        """
        pass

    def gravity_slider(self):
        """
            slider template code
        """
        pass

    def buoyancy_slider(self):
        """
            slider template code
        """
        pass

    def pressure_slider(self):
        """
            slider template code
        """
        pass

    def viscosity_slider(self):
        """
            slider template code
        """
        pass

    def massD_slider(self):
        """
            slider template code
        """
        pass

    def speed_loss_slider(self):
        """
            slider template code
        """
        pass

    def delta_time_slider(self):
        """
            slider template code
        """
        pass
        
# ---------------------------------- BUTTON SLOTS -----------------------------------------

    def reset_btn_func(self):
        """
            controls reset button functionality
        """
        pass

# ---------------------------------- SPINBOX SLOTS ----------------------------------------
    def particle_num_spin_box(self):
        pass

    def particle_size_spin_box(self):
        pass
    
    def mPosition_spin_box(self):
        """
            spin box template code...
        """
        pass

    def mVel_spin_box(self):
        """
            spin box template code...
        """
        pass

    def mAccel_spin_box(self):
        """
            spin box template code...
        """
        pass

    def particle_sep_spin_box(self):
        """
            spin box template code...
        """
        pass

    def cell_size_spin_box(self):
        """
            spin box template code...
        """
        pass

    def createSpinBox(self):
        """
            spin box template code...
        """
        pass

    def mass_spin_box(self):
        """
            spin box template code...
        """
        pass

    def gravity_spin_box(self):
        """
            spin box template code...
        """
        pass


    def buoyancy_spin_box(self):
        """
            spin box template code...
        """
        pass

    def viscosity_spin_box(self):
        """
            spin box template code...
        """
        pass

    def massD_spin_box(self):
        """
            spin box template code...
        """
        pass

    def pressure_spin_box(self):
        """
            spin box template code...
        """
        pass

    def speed_loss_spin_box(self):
        """
            spin box template code...
        """
        pass

    def tank_radius_spin_box(self):
        """
            spin box template code...
        """
        pass

    def tank_position_spin_box(self):
        """
            spin box template code...
        """
        pass

    def delta_time_spin_box(self):
        """
            spin box template code...
        """
        pass

    
# ---------------------------------- COMBOBOX SLOTS ---------------------------------------

    def tank_type_combo_box(self):
        """
            combo box template code
        """
        pass

    def solver_combo_box(self):
        """
            combo box template code
        """
        pass

    def distr_combo_box(self):
        """
            combo box template code
        """
        pass

    def time_integrator_combo_box(self):
        """
            combo box template code
        """
        pass
# ------------------------------------ TAB SLOTS ------------------------------------------