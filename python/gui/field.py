from PyQt4 import (QtCore, QtGui)

from PyQt4.QtCore import (Qt, QVariant, pyqtSignal, pyqtSlot)

import addpaths

from dicts import dotdict

class PropagatorModel(QtCore.QAbstractListModel):

  #select = pyqtSignal(str)
  
  def __init__(self,data=[],parent=None):
    QtCore.QAbstractListModel.__init__(self, parent)
    self.listdata = data

  def rowCount(self,parent):
    return len(self.listdata)

  def data(self, index, role):
    if index.isValid() and role == Qt.DisplayRole:
      return QVariant(self.listdata[index.row()])
    else:
      return QVariant()

  def test(self, selection):
    if len(selection.indexes()) > 0:
      index = selection.indexes()[0]
      shape = index.data().toString()
      #self.select.emit(shape)


class GridTableModel(QtGui.QStandardItemModel):
  def __init__(self, data, horzHeader, vertHeader, parent=None):
    """
    Args:
        datain: a list of lists\n
        horHeader: a list of strings
        verHeader: a list of strings
    """
    
    QtGui.QStandardItemModel.__init__(self, len(data), len(data[0]), parent)
    self.arraydata = data
    self.horzHeader = horzHeader
    self.vertHeader = vertHeader

    types = [int, float, float, float, float]
    
    for row in range(self.rowCount()):
      for col in range(self.columnCount()):
        index = self.index(row,col)
        self.setData(index,types[row](data[row][col]))
    self.setHorizontalHeaderLabels(horzHeader)
    self.setVerticalHeaderLabels(vertHeader)
    self.itemChanged.connect(self.on_itemChanged)

  def flags(self, index):
    """
    Only the first 3 rows are enabled
    """
    f = Qt.ItemIsSelectable
    if (index.row() < 3):
      f = f | Qt.ItemIsEditable | Qt.ItemIsEnabled
    return f

  def on_itemChanged(self, item):
    """
    Update min/max values
    """
    col = item.index().column()

    N,retval = self.data(self.index(0,col),Qt.DisplayRole).toInt()
    offset,_ = self.data(self.index(1,col),Qt.DisplayRole).toFloat()
    delta,_  = self.data(self.index(2,col),Qt.DisplayRole).toFloat()
    colmin = (-(N-1.0)/2 + offset)*delta
    colmin = float("{:.3f}".format(colmin))
    colmax = ((N-1.0) - (N-1.0)/2 + offset)*delta
    colmax = float("{:.3f}".format(colmax))

    self.blockSignals(True)
    self.setData(self.index(3,col),colmin)
    self.setData(self.index(4,col),colmax)
    self.blockSignals(False)

  def getData(self):
    opt = dotdict({})
    # Consider returning full data
    opt['nx'], _    = self.data(self.index(0,0),Qt.DisplayRole).toInt()
    opt['dx']       = float(self.data(self.index(2,0),Qt.DisplayRole).toString())
    opt['offset_x'] = float(self.data(self.index(1,0),Qt.DisplayRole).toString())
    opt['nz'], _    = self.data(self.index(0,2),Qt.DisplayRole).toInt()
    opt['dz']       = float(self.data(self.index(2,2),Qt.DisplayRole).toString())
    opt['offset_z'] = float(self.data(self.index(1,2),Qt.DisplayRole).toString())
    opt['ny'], _    = self.data(self.index(0,1),Qt.DisplayRole).toInt()
    opt['dy']       = float(self.data(self.index(2,1),Qt.DisplayRole).toString())
    opt['offset_y'] = float(self.data(self.index(1,1),Qt.DisplayRole).toString())
    return opt
    
class GridElementDelegate(QtGui.QItemDelegate):
  """
  Uses type alone for distinction
  """
  def __init__(self, parent=None):
    QtGui.QItemDelegate.__init__(self, parent)
    
  def createEditor(self, parent, option, index):
    if index.model().data(index,Qt.DisplayRole).type() == QVariant.Int:
      editor = QtGui.QSpinBox(parent)
      editor.setMinimum(1)
      editor.setMaximum(1024)
    elif index.row() == 1:
      editor = QtGui.QDoubleSpinBox(parent)
      editor.setSingleStep(0.5)
      editor.setMaximum(1024.0)
      editor.setMinimum(-1024.0)
      editor.setDecimals(1)
    else:
      editor = QtGui.QLineEdit(parent)
    return editor
        
  def setEditorData(self, editor, index):
    if index.model().data(index,Qt.DisplayRole).type() == QVariant.Int:
      value, retval = index.model().data(index,Qt.EditRole).toInt()
      spinbox = editor
      spinbox.setValue(value)
    elif index.row() == 1:
      value, retval = index.model().data(index,Qt.EditRole).toDouble()
      spinbox = editor
      spinbox.setValue(value)
    else:
      value = index.model().data(index,Qt.EditRole).toString()
      lineedit = editor
      lineedit.setText(value)

  def setModelData(self, editor, model, index):
    if index.model().data(index, Qt.DisplayRole).type() == QVariant.Int:
      spinbox = editor
      spinbox.interpretText()
      value = int(spinbox.value())
    elif index.row() == 1:
      spinbox = editor
      spinbox.interpretText()
      value = float(spinbox.value())
    else:
      lineedit = editor
      value = lineedit.text()
    model.setData(index, value, Qt.EditRole)

  def updateEditorGeometry(self, editor, option, index):
    editor.setGeometry(option.rect)

class ReadOnlyDelegate(QtGui.QItemDelegate):
  def __init__(self, parent=None):
    QtGui.QItemDelegate.__init__(self, parent)
    
  def createEditor(self, parent, option, index):
    return None
    
