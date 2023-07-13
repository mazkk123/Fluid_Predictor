from PySide6.QtWidgets import QApplication, QMainWindow
from main_window import MainWindow
import sys

app = QApplication(sys.argv)
window = MainWindow(app)
window.show()
app.exec()
