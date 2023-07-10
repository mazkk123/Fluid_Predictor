from PySide6.QtWidgets import QApplication, QMainWindow, QPushButton
from widgets import MainWindow
import sys

app = QApplication(sys.argv)
window = MainWindow(app)
window.show()
app.exec()
