from PyQt4 import (QtCore, QtGui)
from PyQt4.QtCore import (Qt,QVariant, pyqtSignal, pyqtSlot)

import addpaths
from dicts import dotdict

xdc_data = dotdict({'Linear'                : dotdict({'N'      : int,
                                                       'width'  : float,
                                                       'kerf'   : float,
                                                       'pitch'  : float,
                                                       'height' : float,}),
                    'Linear (focused)'      : dotdict({'N'      : int,
                                                       'width'  : float,
                                                       'kerf'   : float,
                                                       'pitch'  : float,
                                                       'M'      : int,
                                                       'height' : float,
                                                       'efocus'  : float}),
                    'Curvelinear'           : dotdict({'N'      : int,
                                                       'width'  : float,
                                                       'kerf'   : float,
                                                       'pitch'  : float,
                                                       'radius' : float,
                                                       'height' : float,}),
                    'Curvelinear (focused)' : dotdict({'N'      : int,
                                                       'width'  : float,
                                                       'kerf'   : float,
                                                       'pitch'  : float,
                                                       'radius' : float,
                                                       'M'      : int,
                                                       'height' : float,
                                                       'efocus'  : float,})})

class ElevationPlacementModel(QtCore.QAbstractListModel):

  updated = pyqtSignal(str)
  
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

  def update(self, selection):
    if len(selection.indexes()) > 0:
      index = selection.indexes()[0]
      placement = index.data().toString()
      self.dataChanged.emit(index, index, ())
      #self.updated.emit(placement)

class TransducerShapeModel(QtCore.QAbstractListModel):

  select = pyqtSignal(str)
  
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
      self.select.emit(shape)

class TransducerElementModelProxy(QtGui.QSortFilterProxyModel):
  def __init__(self,parent=None):
    super(TransducerElementModelProxy,self).__init__(parent)
    self.enabledColumns = range(8)
    for i in [4,5,7]:
      self.enabledColumns.remove(i)
    self.setFilterRole(Qt.DisplayRole)

  def filterAcceptsColumn(self, column, parent):
    if ( column in self.enabledColumns): return True
    return False

  @QtCore.pyqtSlot(str)
  def on_transducer_shape_select(self, shape):
    self.layoutAboutToBeChanged.emit()
    enabledColumns = range(8)
    if not(shape.count('focused')):
      enabledColumns.remove(7)
      enabledColumns.remove(5)
    if not(shape.count('Curvelinear')):
      enabledColumns.remove(4)
    self.enabledColumns = enabledColumns
    self.setFilterRegExp("")
    self.layoutChanged.emit()

  def getData(self):
    opt = dotdict({})
    for i in range(self.columnCount()):
      exec("opt['%s']" % (str(self.headerData(i, Qt.Horizontal, Qt.DisplayRole).toString())) + "=%s" %(self.data(self.index(0,i),Qt.DisplayRole).toString()))
    return opt

class TransducerElementModel(QtGui.QStandardItemModel):
  def __init__(self,defaults,parent=None):
    QtGui.QStandardItemModel.__init__(self, 1, 8, parent)

    labels = ["N",
              "width", 
              "kerf",      
              "pitch", 
              "radius",
              "M",     
              "height",
              "efocus"]
    types = [int,
             float,
             float,
             float,
             float,
             int,
             float,
             float]

    for row in range(self.rowCount()):
      for col in range(self.columnCount()):
        index = self.index(row,col)
        self.setData(index,types[col](defaults[col]))

    self.setHorizontalHeaderLabels(labels)
    self.setVerticalHeaderLabels(["value"])    


class TransducerElementDelegate(QtGui.QItemDelegate):
  """
  Uses type alone for distinction
  """
  def __init__(self, parent=None):
    QtGui.QItemDelegate.__init__(self, parent)
    
  def createEditor(self, parent, option, index):
    if index.model().data(index,QtCore.Qt.DisplayRole).type() == QVariant.Int:
      editor = QtGui.QSpinBox(parent)
      editor.setMinimum(0)
      editor.setMaximum(192)
    else:
      editor = QtGui.QLineEdit(parent)
    return editor
        
  def setEditorData(self, editor, index):
    if index.model().data(index,QtCore.Qt.DisplayRole).type() == QVariant.Int:
      value, retval = index.model().data(index,QtCore.Qt.EditRole).toInt()
      spinbox = editor
      spinbox.setValue(value)
    else:
      value = index.model().data(index,QtCore.Qt.EditRole).toString()
      lineedit = editor
      lineedit.setText(value)

  def setModelData(self, editor, model, index):
    if index.model().data(index,QtCore.Qt.DisplayRole).type() == QVariant.Int:
      spinbox = editor
      spinbox.interpretText()
      value = int(spinbox.value())
    else:
      lineedit = editor
      value = lineedit.text()
    model.setData(index, value, QtCore.Qt.EditRole)

  def updateEditorGeometry(self, editor, option, index):
    editor.setGeometry(option.rect)

# Edit QTableView from code
# QAbstractItemView::edit ( const QModelIndex & index ) [slot]
# void QAbstractItemView::closeEditor ( QWidget * editor, QAbstractItemDelegate::EndEditHint hint ) [virtual protected slot]
# void QAbstractItemView::editorDestroyed ( QObject * editor )   [virtual protected slot]


# QModelIndex index = ui->tableView->model()->index(0, 0, QModelIndex());
# ui->tableView->edit(index);

# Select items
# tableView->selectionModel()->select(tabelView->model()->index(row,colum), QItemSelectionModel::Select);

# Select and edit next cell in a row
# QModelIndex index = tableViewPowerDegree->currentIndex();
#	int row = index.row() + 1;
#	int column = 1;
#	QModelIndex newIndex  = tableViewPowerDegree->model()->index(row,column);
#	tableViewPowerDegree->selectionModel()->select(newIndex, QItemSelectionModel::Select);
#	tableViewPowerDegree->setCurrentIndex(newIndex);
#	tableViewPowerDegree->setFocus();
#	tableViewPowerDegree->edit(newIndex);

# Disable editing of QTableView
# inherit and override edit(const QModelIndex & index, EditTrigger trigger, QEvent * event)
#
# or
#
# call TableView->setEditTriggers(QAbstractItemView::NoEditTriggers)
#
# or
#
# use TableView->setItemDelegateForColumn() and make read-only delegate

class GridTableModel(QtCore.QAbstractTableModel):
    def __init__(self, data, horzHeader, vertHeader, parent=None):
        """
        Args:
            datain: a list of lists\n
            horHeader: a list of strings
            verHeader: a list of strings
        """
        QtCore.QAbstractTableModel.__init__(self, parent)
        self.arraydata = data
        self.horzHeader = horzHeader
        self.vertHeader = vertHeader

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        if len(self.arraydata) > 0: 
            return len(self.arraydata[0]) 
        return 0

    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        elif role != Qt.DisplayRole:
            return QVariant()
        return QVariant(self.arraydata[index.row()][index.column()])

    def setData(self, index, value, role = Qt.EditRole):
      if role == Qt.EditRole:
        setattr(self.arraydata[index.row()], self.columns[index.column()], value)
        self.dataChanged.emit(index, index, ())
        return True
      else:
        return False 
      
    def headerData(self, index, orientation, role):
        if role == Qt.DisplayRole:
          if orientation == Qt.Horizontal:
            return QVariant(self.horzHeader[index])
          else:
            return QVariant(self.vertHeader[index])
        return QVariant()

# Use QStandardItem and setFlags(!Qt.ItemIsEditable);

# for d in data:
#     row = []
#     for name in d:
#         item = QStandardItem(name)
#         item.setEditable(False)
#         row.append(item)
#     model.appendRow(row)
# setHorizontalHeaderItem
# setItem

    
      
# class TransducerElementModel(QtGui.QAbstractProxyModel):
#   def __init__(self, parent=None):
#     QtGui.QAbstractProxyModel.__init__(self, parent)

    

# class TransducerElementModels(QtGui.QStandardItemModel):
#   def __init__(self, parent=None, data=None):
#     QtGui.QStandardItemModel.__init__(self, parent)
#     self.modelIDs = data.keys()
#     for key in self.modelIDs:
#       parent = self.invisibleRootItem()
#       item   = QStandardItem(QString(key))
#       parent.appendRow(parent)
#       parent = item
#       parent.appendRows([QStandardItem(QString(param)) for param in data[key].keys()])
      
#      if index.isValid():
#        self.items[index.row()] = value.toString()
#      return True


# Write a proxy model derived from QIdentityProxyModel, and implement the QIdentityProxyModel::headerData() function and return what ever data is required in in header.

# In Qt::Vertical header, return nothing for Qt::ItemDataRole and return icon for Qt::DecorationRole
# In Qt::Horizontal header return the desired label in Qt::ItemDataRole






    
      
# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
