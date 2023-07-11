from PySide6.QtWidgets import QVBoxLayout, QWidget, QMainWindow, QPushButton, \
    QToolButton, QToolBar, QTabWidget, QHBoxLayout, QVBoxLayout, QComboBox, QGroupBox, \
    QLabel, QLineEdit, QSlider, QDoubleSpinBox, QSpinBox, QColormap, QRadioButton, QGridLayout 
from PySide6.QtCore import QSize, QRect
from PySide6.QtGui import QPixmap, QIcon

class MainWindow(QMainWindow):

    def __init__(self, app):
        super().__init__()
        self.app = app

        self.setGeometry(QRect(0, 0, 1280, 720))
        self.setWindowTitle("Fluid Predictor Application")
        self.window_icon = QPixmap("images/Toolbar/fluid.png")
        self.setWindowIcon(self.window_icon)

        self.add_menu_bars()
        self.add_tool_bars()

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

    def create_fluid_solver_widget(self):
        """
            creates the widget for fluid solvers
        """
        self.solver_vertical_layout = QVBoxLayout()
        self.solver_horizontal_layout = QHBoxLayout()

        self.fluid_solver_group_box = QGroupBox("Solver Type")

        self.fluid_solver_label = QLabel("Choose Solver")
        self.fluid_solver_combo_box = QComboBox()
        self.fluid_solver_combo_box.addItem("SPH")
        self.fluid_solver_combo_box.addItem("Eulerian")

        self.solver_horizontal_layout.addWidget(self.fluid_solver_label)
        self.solver_horizontal_layout.addWidget(self.fluid_solver_combo_box)

        self.fluid_solver_group_box.setLayout(self.solver_vertical_layout)
        self.main_layout.addWidget(self.fluid_solver_group_box)

# -------------------- PARTICLE WIDGETS ------------------------------

    def create_appearance_widgets(self):
        """
            creates widgets responsible for fluid appearance
        """
        self.appearance_group_box = QGroupBox("Appearance")
        self.appearance_v_layout = QVBoxLayout()

        self.nbr_particles_h_layout = QHBoxLayout()
        self.nbr_of_particles_lbl = QLabel("No. of particles")
        self.nbr_of_particles_slider = QSlider(orientation="horizontal")

        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_lbl)
        self.nbr_particles_h_layout.addWidget(self.nbr_of_particles_slider)

        self.particle_size_h_layout = QHBoxLayout()
        self.particle_size_lbl = QLabel("Particle Size")
        self.particle_size_slider = QSlider(orientation="horizontal")
        
        self.particle_size_h_layout.addWidget(self.particle_size_lbl)
        self.particle_size_h_layout.addWidget(self.particle_size_slider)

        self.particle_colour_h_layout = QHBoxLayout()
        self.particle_colour_label = QLabel("Particle colour")
        self.particle_colour = QColormap()

        self.particle_colour_h_layout.addWidget(self.particle_colour_label)
        self.particle_colour_h_layout.addWidget(self.particle_colour)

        self.appearance_v_layout.addWidget(self.appearance_group_box)

    def create_motion_widget(self):
        """
            creates the widget for the motion of the particles
        """
        pass

    def create_physical_widget(self):
        """
            creates physical forces widget appearance
        """
        pass


# ----------------------- TANK WIDGET -------------------------------

    def create_fluid_tank_widgets(self):
        """
            creates the widgets for the fluid tank session of the main window UI
        """
        pass

# -------------------- MAIN SYSTEM WIDGETS ---------------------------------

    def create_system_widgets(self):
        """
            resposible for the system of particle controls
        """
        pass


    def add_menu_bars(self):
        """
            responsible for menu bar widgets and dockers in central main window
        """

        # creating the file menu toolbar
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
