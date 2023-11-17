
# ------------------------------------------------- ----- 
# --------------------- mplwidget.py -------------------- 
# ---------------------------------------------------------------- ---- 
from PyQt5.QtWidgets import QWidget, QVBoxLayout

from matplotlib.backends.backend_qt5agg import FigureCanvas

from matplotlib.figure import Figure

class MplWidget(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        self.setLayout(vertical_layout)
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        #self.canvas.axes.set_xlabel("Energy (eV)")
        #self.canvas.axes.set_ylabel("Intensity (a.u.)")
           