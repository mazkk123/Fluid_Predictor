from PySide6.QtWidgets import QVBoxLayout, QWidget, QMainWindow, QPushButton
import sys

class MainWindow(QMainWindow):

    def __init__(self, app):
        super().__init__()
        self.app = app

        menuBar = self.menuBar()
        file = menuBar.addMenu("File")
        quit_action = file.addAction("Quit")
        quit_action.triggered.connect(self.quit_app)

    def quit_app(self):
        self.app.quit()
