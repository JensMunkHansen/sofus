import sys
import os
from PyQt4 import QtGui, QtCore
from PyQt4.uic import loadUiType
from PyQt4.QtCore import (QString,QRect,QObject,QRunnable, QThreadPool, Qt)
from PyQt4.Qt import (SIGNAL, SLOT, QObject)
import weakref
import multiprocessing

# TODO: Make changing the elevationPlacement trigger update
#       On propagator changed: Set available order

from PyQt4.QtCore import (pyqtSignal)

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import (art3d, Axes3D)

import numpy as np

from transducer import (TransducerShapeModel,
                        TransducerElementModel,
                        TransducerElementModelProxy,
                        TransducerElementDelegate,
                        ElevationPlacementModel)

from field import (GridTableModel,
                   GridElementDelegate,
                   ReadOnlyDelegate,
                   PropagatorModel)

import addpaths
from dicts import dotdict

import swig_fnm as fnm
from fnm_arrays import (linear_array, linear_array3, convex_array, convex_array3)

multithreaded = True

_emitterCache = weakref.WeakKeyDictionary()

def log_compress(pressure,dBrange=60):
  logp = np.abs(pressure)
  logp = logp / logp.max()
  logp = np.maximum(logp,np.finfo(np.float32).eps)
  logp = 20*np.log10(logp)
  logp[logp < -dBrange] = -dBrange
  #logp[0,0] = 0
  #logp[-1,-1] = -dBrange
  return logp

ui, QMainWindow = loadUiType('design.ui')

class WorkerSignals(QObject):
    result = pyqtSignal(int)
        
class Worker(QRunnable):
  def __init__(self,task,data=None):
    super(Worker, self).__init__()
    self.data = data
    self.task = task
    self.signals = WorkerSignals()
      
  def run(self):
    # Update graphics
    try:
      if not(self.data.has_key('M')):
        self.data.M = 1
      
      if not(self.data.has_key('efocus')):
        self.data.efocus = None
          
      if self.data.has_key('radius'):
        a = convex_array3(nElements = self.data.N,
                          nSubH     = self.data.M,
                          pitch     = self.data.pitch,
                          kerf      = self.data.kerf,
                          height    = self.data.height,
                          efocus    = self.data.efocus,
                          radius    = self.data.radius)
      elif self.data.M > 1:
        a = linear_array3(nElements = self.data.N,
                          nSubH     = self.data.M,
                          pitch     = self.data.pitch,
                          kerf      = self.data.kerf,
                          height    = self.data.height,
                          efocus    = self.data.efocus,
                          elePlacement = self.data.elePlacement)
      else:
        a = linear_array(nElements = self.data.N,
                         nSubH     = self.data.M,
                         pitch     = self.data.pitch,
                         kerf      = self.data.kerf,
                         height    = self.data.height,
                         efocus    = self.data.efocus)
      
      self.data.ax.clear()
      a.show(ax=self.data.ax)
    except Exception as inst:
      exc_type, exc_obj, exc_tb = sys.exc_info()
      print(str(type(inst)))
      print(str(inst.args))
      fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
      print(str(exc_type) + fname + ' ' + str(exc_tb.tb_lineno))
    finally:
      self.signals.result.emit(self.task)

class Tasks(QObject):
  def __init__(self):
    super(Tasks, self).__init__()
    self.pool = QThreadPool.globalInstance()
    self.pool.setMaxThreadCount(2)

  def process_result(self, task):
    print 'Receiving', task

  def start(self):
    for task in range(1):
      worker = Worker(task)
      worker.signals.result.connect(self.process_result)
      self.pool.start(worker)
      self.pool.waitForDone()

class ComputeSignals(QObject):
    result = pyqtSignal(int)


def emitter(ob):
  if ob not in _emitterCache:
    _emitterCache[ob] = QObject()
  return _emitterCache[ob]
    
    
class Callback(fnm.ProgressBarInterface):
  def show(self,f):
    QObject.emit(emitter(self), SIGNAL('test'), f)
        
# Move elsewhere
class ComputeField(QRunnable):
  def __init__(self,task,data=None):
    super(ComputeField, self).__init__()
    self.task = task
    self.data = data
    self.signals = ComputeSignals()
      
  def run(self):
    if not(self.data.has_key('M')):
        self.data.M = 1

    nx = self.data.nx
    ny = self.data.ny
    nz = self.data.nz
    
    N  = self.data.N
    nCycles = self.data.nCycles
    f0 = self.data.f0
    fs = self.data.fs
    if not(self.data.has_key('efocus')):
      self.data.efocus = 0.0

    if (self.data.has_key('radius')):
      print('convex')
      a = convex_array3(nElements = self.data.N,
                        nSubH     = self.data.M,
                        pitch     = self.data.pitch,
                        kerf      = self.data.kerf,
                        height    = self.data.height,
                        efocus    = self.data.efocus,
                        radius    = self.data.radius)
    else:
      a = linear_array3(nElements = self.data.N,
                        nSubH     = self.data.M,
                        pitch     = self.data.pitch,
                        kerf      = self.data.kerf,
                        height    = self.data.height,
                        efocus    = self.data.efocus,
                        elePlacement = self.data.elePlacement)
    a.f0 = f0
    a.c  = 1500
    a.nDivW = self.data.order
    a.nDivH = self.data.order
    a.nthreads = multiprocessing.cpu_count()
    focus   = [0,0, self.data.focus]
    a.focus = focus 

    a.att_enabled = self.data.use_att
    a.alpha       = self.data.alpha
    a.beta        = self.data.beta

    a.focus_type = fnm.FocusingType.Pythagorean
    if (self.data.focus < 0.0):
      delays = np.sqrt(np.sum((a.pos - np.r_[self.data.N*[focus]])**2,axis=1)) / a.c
      dmin   = np.min(delays)
      delays = -(dmin - delays)
      a.delays = delays.astype(np.float32)
      a.focus_type = fnm.FocusingType.Delays
    
    apodization = a.apodization
    nActive = int((self.data.focus / self.data.f) / self.data.pitch)
    apodization = np.zeros(apodization.shape,np.float32)

    if (N % 2) == 1:
      apodization[N/2] = 1.0
    apodization[N/2:int(N/2 + nActive/2)]      = 1.0
    apodization[max(N/2-int(nActive/2),0):N/2] = 1.0
    a.apodization = apodization

      
    # TODO: Support any plane
    xs1 = (np.r_[0:nx] - (nx-1.0)/2.0 + self.data.offset_x) * self.data.dx
    zs1 = (np.r_[0:nz] - (nz-1.0)/2.0 + self.data.offset_z) * self.data.dz

    xs,zs = np.meshgrid(xs1,zs1,indexing='ij')

    pos = np.c_[xs.flatten(), np.zeros(nx*nz), zs.flatten()].astype(np.float32)

    a.focus = [0,0,self.data.focus]

    if self.data.simtype == 0:
      a.focus_type = fnm.FocusingType.Rayleigh
      out = a.CalcCwFast(pos)[1]
      result3 = np.squeeze(out.reshape((nx,ny,nz)))

      rho = 10000
      omega = 2*np.pi*1e6#f0
      pressure = -1j * omega * rho * np.exp(1j * omega * 0) * result3
    else:
      a.fs = fs
      a.c  = 1540.0

      # Set order
      sysparm   = fnm.SysParmFloat()
      sysparm.c = 1540.0
      sysparm.pulseWaveIntOrder  = int(self.data.order)
      sysparm.timeDomainCalcType = int(self.data.proptype)  
      a.SysParmSet(sysparm)
      
      impulse_response = np.sin(np.linspace(0, nCycles * 2 * np.pi, round(nCycles * fs / f0)))
      impulse_response = impulse_response * np.hamming(len(impulse_response))
      excitation       = np.sin(np.linspace(0, nCycles * 2 * np.pi, round(nCycles * fs / f0)))
      if not(len(excitation) < 2):
        a.excitation = excitation.astype(np.float32)
        a.impulse    = impulse_response.astype(np.float32)

      # Consider setting callback outside
      #QObject.connect(emitter(cb), SIGNAL('test'), main, SLOT('on_progress(float)'), Qt.QueuedConnection)
      # AutoConnection, DirectConnection, QueuedConnection, BlockingQueuedConnection, UniqueConnection, AutoCompatConnection

      if multithreaded:
        tstart, hp1 = a.CalcPwFieldThreaded(pos)
      else:
        cb = Callback()
        a.ProgressBarSet(cb)
        QObject.connect(emitter(cb), SIGNAL('test'), main.on_progress)
        tstart, hp1 = a.CalcPwField(pos)        
      hp = np.sum(hp1**2,axis=1)
      hp = np.sqrt(hp)
      pressure = np.squeeze(hp.reshape((nx,ny,nz)))

    # Remove NaN
    pressure[~np.isfinite(pressure)] = np.finfo(np.float32).eps
    pressure = np.abs(pressure)
    logp = log_compress(pressure)
    self.data.ax.clear()

    nDim = 3 - [nx,ny,nz].count(1)

    # Fix - find dimension which are not singleton
    indices = [i for i, x in enumerate([nx,ny,nz]) if x == 1]
    
    if nDim == 2:
      self.data.ax.imshow(logp,extent=[zs.min(),zs.max(),xs.min(),xs.max()],aspect='auto')
      self.data.ax.set_xlabel('Depth [m]')
      self.data.ax.set_ylabel('Width [m]')
    elif nDim == 1:
      if nx == 1:
        self.data.ax.plot(zs1, logp.flatten())
        self.data.ax.set_xlabel('Depth [m]')
        self.data.ax.set_ylabel('Gain [dB]')

    # Plot shit
    self.signals.result.emit(self.task)
      
class Main(QMainWindow, ui):

  orderCW = 20
  orderPW = 3
  
  def __init__(self, ):
    # Consider using multiple QStringListModels
    # QVXYModelMapper can map columns or rows to an xy-model
    
    super(Main, self).__init__()
    self.setupUi(self)
    self.action_Exit.triggered.connect(self.close)

    self.modelShape = TransducerShapeModel(["Linear",
                                            "Linear (focused)",
                                            "Curvelinear",
                                            "Curvelinear (focused)"])

    # TODO: connect to model, which emits updated
    #self.cboxElePlacement.view().selectionModel().selectionChanged.connect(self.modelElevationPlacement.update)

    # TODO: update parameter, which is used for selecting inner or output
    #self.modelElevationPlacement.updated.connect(self.proxy.on_elevation_placement_select)
    
    self.cboxShape.setModel(self.modelShape)
    self.cboxShape.view().selectionModel().selectionChanged.connect(self.modelShape.test) # model emits select

    # TODO: How to use a model for a comboBox???
    #self.modelElevationPlacement = ElevationPlacementModel(["Inner", "Outer"])
    #self.cboxElePlacement.setModel(self.modelElevationPlacement)
    #self.cboxElePlacement.view().selectionModel().selectionChanged.connect(self.update_image)

    # Temporary solution without model
    self.cboxElePlacement.setCurrentIndex(0)
    self.cboxElePlacement.currentIndexChanged[QString].connect(self.on_cbox_ele_placement_changed)
    
    self.delegate = TransducerElementDelegate(self)

    # Baffle is disabled for CW (default)
    self.cboxBaffle.setDisabled(1)
    
    # CW defaults
    defaults = [128,
                3e-3,
                5.0e-4,
                3.5e-3,
                61e-3,
                5,
                50e-3,
                8.5e-2]
    # PW defaults
    defaults = [64,
                1.8e-4,
                2e-5,
                2e-4,
                61e-3,
                5,
                20e-3,
                1.5e-2]
    self.model = TransducerElementModel(defaults, self)

    self.proxy = TransducerElementModelProxy(self)
    self.proxy.setSourceModel(self.model)

    self.modelShape.select.connect(self.proxy.on_transducer_shape_select) # select triggers update 

    self.tblElements.setModel(self.proxy)
    self.tblElements.setItemDelegate(self.delegate)
    
    self.tblElements.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch) #QtGui.QHeaderView.ResizeToContents
    self.tblElements.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)
    
    self.tblElements.setVisible(False)
    self.tblElements.resizeColumnsToContents()
    self.tblElements.resizeRowsToContents()
    self.tblElements.setVisible(True)

    self.pool = QThreadPool.globalInstance()
    self.pool.setMaxThreadCount(4)
    self.active = False

    self.proxy.dataChanged.connect(self.update_image)
    #self.modelElevationPlacement.dataChanged.connect(self.update_image)
    #self.modelElevationPlacement.updated.connect(self.update_image)
    
    self.cboxShape.view().selectionModel().selectionChanged.connect(self.update_image)

    # TODO: Make model work
    #self.propModel = PropagatorModel(["Point", "Plane", "Sphere"], self)
    #self.propModel = PropagatorModel(["Sphere"], self)
    #self.cboxPropagator.setModel(self.propModel)
    
    self.sboxOrder.setValue(2)
    self.sboxOrder.setMinimum(2)
    self.sboxOrder.setMaximum(500)

    self.cboxDomain.setCurrentIndex(0)
    self.cboxDomain.currentIndexChanged[QString].connect(self.on_cbox_domain_changed)
    self.leSampleFrq.setDisabled(True)
    self.leExcitationCycles.setDisabled(True)

    self.sboxOrder.valueChanged[int].connect(self.on_sbox_order_changed)
    
    # CW grid
    data = [[170, 1, 250],
            [0, 0, 125],
            [ 0.004,  0.1,  0.0036],
            [-0.338,  0.00, 0.0018],
            [ 0.3380, 0.00, 0.8982]]
    # PW grid
    data = [[100, 1, 100],
            [0, 0, 60],
            [ 6e-5,   0.10, 4e-4],
            [-0.003,  0.00, 0.00],
            [ 0.003,  0.00, 0.04]]
    
    self.gridXYZ = GridTableModel(data, ['X','Y','Z'], ['N','offset', 'delta','min','max'], self)

    self.tblGrid.setModel(self.gridXYZ)
    self.tblGrid.setItemDelegate(self.delegate)
    
    self.tblGrid.horizontalHeader().setResizeMode(QtGui.QHeaderView.Stretch)
    self.tblGrid.verticalHeader().setResizeMode(QtGui.QHeaderView.Fixed)

    self.tblGrid.setVisible(False)
    self.tblGrid.resizeColumnsToContents()
    self.tblGrid.resizeRowsToContents()

    self.gridDelegate = GridElementDelegate(self)
    self.readOnlyDelegate = ReadOnlyDelegate(self)
    self.tblGrid.setItemDelegateForRow(0,self.gridDelegate)
    self.tblGrid.setItemDelegateForRow(1,self.gridDelegate)
    self.tblGrid.setItemDelegateForRow(2,self.gridDelegate)
    self.tblGrid.setItemDelegateForRow(3,self.readOnlyDelegate)
    self.tblGrid.setItemDelegateForRow(4,self.readOnlyDelegate)
    self.tblGrid.setVisible(True)

    self.progressBar.setRange(0,100)

    self.leFNumber.setAlignment(Qt.AlignRight)
    self.btnField.clicked.connect(self.on_calc_clicked)
    #self.proxy.dataChanged.connect(self.tblGrid.adjust)
    
    # TODO: Add signal and slot to event loop responsible for updating graph
    # TODO: Use TableView's sizechanged signal to model, emit dataChanged

    # Implement TableView::sizeHint() {return minimumSizeHint()}
    # MyTableView::mimimumSizeHint() { // compute width}

  @QtCore.pyqtSlot(float)
  def on_progress(self, f):
    self.progressBar.setValue(int(f))
    self.progressBar.repaint()
    #QtGui.QApplication.processEvents()

  def on_cbox_ele_placement_changed(self, index):
    inner = index.count('Inner') > 0
    self.update_image()

  def on_sbox_order_changed(self, index):
    index = self.cboxDomain.currentIndex()
    if (index == 0):
      Main.orderCW = index
    else:
      Main.orderPW = index
    
  def on_cbox_domain_changed(self, index):
    cw = index.count('Continous') > 0
    if cw:
      self.leSampleFrq.setDisabled(1)
      self.leExcitationCycles.setDisabled(1)
      self.sboxOrder.blockSignals(True)
      self.sboxOrder.setMinimum(2)
      self.sboxOrder.setMaximum(300)
      self.sboxOrder.setValue(Main.orderCW)
      self.sboxOrder.blockSignals(False)

      self.cboxPropagator.blockSignals(True)
      # Populate model with Sphere only
      index = self.cboxPropagator.findText("Sphere", QtCore.Qt.MatchFixedString)
      if ( index != -1 ):
        self.cboxPropagator.setCurrentIndex(index);
      self.cboxBaffle.setDisabled(1)
      self.cboxPropagator.blockSignals(False)
    else:
      self.leSampleFrq.setEnabled(1)
      self.leExcitationCycles.setEnabled(1)
      self.sboxOrder.blockSignals(True)
      value = self.sboxOrder.value()
      _min = 1
      _max = 5
      value = max(value,_min)
      value = min(value,_max)
      self.sboxOrder.setValue(Main.orderPW)
      self.sboxOrder.setMinimum(_min)
      self.sboxOrder.setMaximum(_max)
      self.cboxBaffle.setEnabled(True)
      self.sboxOrder.blockSignals(False)
  
  def addmpl(self, fig):
    self.fig = fig
    self.canvas = FigureCanvas(fig)
    self.mplvl.addWidget(self.canvas)
    self.canvas.draw()
    self.toolbar = NavigationToolbar(self.canvas, 
                                     self.mplwindow, coordinates=True)
    self.mplvl.addWidget(self.toolbar)
    ax = fig.get_axes()[0]
    Axes3D.mouse_init(ax)
    
  def addfield(self, fig):
    self.fig1 = fig
    self.canvas1 = FigureCanvas(fig)
    self.mplvl1.addWidget(self.canvas1)
    self.canvas1.draw()
    self.toolbar1 = NavigationToolbar(self.canvas1, 
                                      self.mplwindow1, coordinates=True)
    self.mplvl1.addWidget(self.toolbar1)
    ax = fig.get_axes()[0]

  def done(self,task):
    self.canvas.draw()
    self.active = False

  def update_image(self):
    options = {}
    options = self.proxy.getData()
    options.ax = self.fig.gca()

    # TODO: Use model instead
    options.elePlacement = str(self.cboxElePlacement.currentText())

    # TODO: 
    worker = Worker(0,options)
    worker.signals.result.connect(self.done)

    if not(self.active):
      self.active = True
      self.pool.start(worker)
    #self.pool.waitForDone()

  def on_calc_clicked(self):
    self.btnField.setEnabled(False)
    self.progressBar.setValue(0)
    # Get transducer model
    options = self.proxy.getData()

    options.ax = self.fig1.gca()
    options.update(self.gridXYZ.getData())

    options.focus  = float(self.leFocusDist.text())
    options.order  = int(self.sboxOrder.value())
    options.f0     = float(self.leExcitationFrq.text())
    options.f      = float(self.leFNumber.text())

    # Hack: Find propagator using text instead
    options.proptype = self.cboxPropagator.currentIndex()

    options.use_att = self.checkAttenuation.isChecked()
    dBcmMHz = float(self.leAttenuation.text())

    options.alpha = dBcmMHz * 100/1e6 * options.f0 * fnm.ApertureFloat.Neper_dB
    options.beta  = options.alpha / options.f0
    
    # TODO: Use model instead
    options.elePlacement = str(self.cboxElePlacement.currentText())
    
    # Hack: Get simulation type
    index = self.cboxDomain.currentIndex()
    options.simtype = index

    options.fs = float(self.leSampleFrq.text())
    # Fractional number of cycles
    options.nCycles = float(self.leExcitationCycles.text())

    # TEST Update
    options.aperture = fnm.ApertureFloat()

    
    
    # Consider thr = QThread(self)
    # worker.moveToThread(thr)
    # connect(worker, SIGNAL(updateprogress), self, updatemyprogress)
    # connect(thr, SIGNAL(destroyed), worker, SLOT(deleteLater()))
    # worker->run()
    
    # Transducer data
    worker = ComputeField(0,options)
    worker.signals.result.connect(self.on_calc_done)
    self.pool.start(worker)

    
  def on_calc_done(self,task):
    self.canvas1.draw()
    self.progressBar.setValue(100)
    self.btnField.setEnabled(True)
                                  

    
if __name__ == '__main__':

  fig = Figure()
  ax = fig.add_subplot(111, projection='3d')
  #fig.tight_layout()
  Axes3D.mouse_init(ax)

  fig1 = Figure()
  ax1 = fig1.add_subplot(111)
  #fig1.tight_layout()

  app = QtGui.QApplication(sys.argv)
  app.setStyle(QtGui.QStyleFactory.create('Fusion'))
  main = Main()
  main.addmpl(fig)
  main.addfield(fig1)
  main.show()
  sys.exit(app.exec_())


# In Qt 5.2, we can do
# view->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContentsOnFirstShow);
# or
# view->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContents);

#self.myTable.verticalHeader().sectionResized.connect(
#  self.row_resized)
#def row_resized(self, index, old_size, new_size):
#    for i in range(self.myTable.verticalHeader().count()):
#        self.myTable.setRowHeight(i, new_size)
    

  
# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
