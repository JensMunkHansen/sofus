import sys
import numpy as np

#--- QT imports
from PyQt4 import QtGui
from PyQt4 import QtCore

from PyQt4.QtGui import *
from PyQt4.QtCore import *

class MyForm(QtGui.QMainWindow):

    k_aspect = {'limits':  None,
                'auto':    'auto'}
    
    # ------------------------------------------------------------------------
    def __init__(self, parent=None):
        super(MyForm, self).__init__(parent)

        central = QWidget(self)
        self.setCentralWidget(central)
        
        self.label_a = QLabel("a");
        self.label_b = QLabel("b");
        self.slider_a = QSlider(Qt.Horizontal);
        self.slider_b = QSlider(Qt.Horizontal);
        self.slider_a.blockSignals(True)
        self.slider_b.blockSignals(True)
        self.slider_a.setMinimum(0)
        self.slider_a.setMaximum(99)
        self.slider_b.setMinimum(100)
        self.slider_b.setMaximum(199)
        
        grid = QGridLayout(central)
        grid.addWidget(self.label_a, 0, 0)
        grid.addWidget(self.slider_a, 1, 0)
        grid.addWidget(self.label_b, 2, 0)
        grid.addWidget(self.slider_b, 3, 0)

        self.connect(self.slider_a, SIGNAL("valueChanged(int)"),
                     self.slider_changed)
#        self.connect(self.slider_b, SIGNAL("valueChanged(int)"),
#                     self.slider_changed)
        self.slider_b.valueChanged.connect(self.slider_changed)
        self.slider_a.blockSignals(False)
        self.slider_b.blockSignals(False)

    def slider_changed(self,value):
        self.slider_a.blockSignals(True)
        self.slider_b.blockSignals(True)
        if (value in range(self.slider_a.minimum(),self.slider_a.maximum()+1)):
            if (self.slider_a.value() + 100 > self.slider_b.value()):
                self.slider_b.setValue(self.slider_a.value() + 100)
        if (value in range(self.slider_b.minimum(),self.slider_b.maximum()+1)):
            if (self.slider_b.value() - 100 < self.slider_a.value()):
                self.slider_a.setValue(self.slider_b.value() - 100)
        self.slider_a.blockSignals(False)
        self.slider_b.blockSignals(False)

        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)

    win = MyForm(None)
    win.show()
    app.exec_()
        
