%include "attribute.i"

#pragma message("hello")

%insert("python") %{
from pprint import pprint
from inspect import getmembers
from types import FunctionType

def attributes(obj):
  disallowed_names = {
    name for name, value in getmembers(type(obj))
      if isinstance(value, FunctionType)}
  return {
    name: getattr(obj, name) for name in dir(obj)
      if name[0] != '_' and name not in disallowed_names and hasattr(obj, name)}

from mpl_toolkits.mplot3d import art3d
from sys import version_info
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import itertools
import numpy as np

from dicts import dotdict
from euler import (euler2rot, rot2euler)
from pyutils import nargout
plt.ion()
%}

/* %fnm_typemaps() macro
 *
 * This macro applies a long list of typemaps for templated functions. It is meant
 * to be executed for float and double
 *
 */
%define %fnm_typemaps(DATA_TYPE)

/* Fixed-length 1-dim arrays */
%apply (DATA_TYPE ARGOUT_ARRAY1[ANY]) {(DATA_TYPE oFocus[3])}
%apply (DATA_TYPE IN_ARRAY1[ANY])     {(const DATA_TYPE iFocus[3])}
%apply (DATA_TYPE IN_ARRAY1[ANY])     {(const DATA_TYPE iEuler[3])}


/* One-dimensional arrays */
%apply (DATA_TYPE* IN_ARRAY1, int DIM1) \
{(const DATA_TYPE* data, const size_t nData)}
%apply (std::complex<DATA_TYPE>* IN_ARRAY1, int DIM1)    \
{(const std::complex<DATA_TYPE>* pFieldValues, const size_t nComplexValues)}
%apply (DATA_TYPE* INPLACE_ARRAY1, size_t DIM1) \
{(DATA_TYPE* ioData, size_t nIOdata)}

/* Two-dimensional arrays */
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* pos, const size_t nPositions, const size_t nDim)};
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* iData, const size_t nChannels, const size_t nSamples)};
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* p0, const size_t nx, const size_t ny)};

%apply (DATA_TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3)             \
{(const DATA_TYPE* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim)};

%apply(DATA_TYPE** ARGOUTVIEW_ARRAY1, size_t* DIM1) {(DATA_TYPE** data, size_t* nData)}
%apply(DATA_TYPE** ARGOUTVIEW_ARRAY1, size_t* DIM1) {(DATA_TYPE** pRefData, size_t* nRefData)}


%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(DATA_TYPE** odata, size_t* nData)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(DATA_TYPE** ppData, size_t* pnData)}
%apply (int** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(int** offsets, size_t* nOffsets)}
%apply (int** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(int** lengths, size_t* nLengths)}
%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(std::complex<DATA_TYPE>** odata, size_t* nOutPositions)};

%apply (int** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(int** sizeAndOffsets, size_t* nFilters, size_t* nTwo)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** out, size_t* nElements, size_t* nParams)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** coordinates, size_t* nDim, size_t* nLimits)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** odata, size_t* nSignals, size_t* nSamples)}
%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(std::complex<DATA_TYPE>** odata, size_t* nOutPositions, size_t* nSamples)};

%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(DATA_TYPE** out, size_t* nElements, size_t* nSubElements, size_t* nParams)}
%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY3, size_t* DIM1,
 size_t* DIM2, size_t* DIM3) {(std::complex<DATA_TYPE>** p1, size_t*
 onx, size_t* ony, size_t* onz)};


/********************************
 *  Used for static constructors
 ********************************/
%typemap(in, numinputs=0) fnm::Aperture<DATA_TYPE> **obj (fnm::Aperture<DATA_TYPE> *temp) {
  $1 = &temp;
}

%typemap(argout) fnm::Aperture<DATA_TYPE> ** {
  PyObject* temp = NULL;
  if (!PyList_Check($result)) {
    temp = $result;
    $result = PyList_New(1);
    PyList_SetItem($result, 0, temp);
  }
  
  // Create shadow object (do not use SWIG_POINTER_NEW)
  temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1),
			    $descriptor(fnm::Aperture<DATA_TYPE>*),
			    SWIG_POINTER_OWN | 0);

  PyList_Append($result, temp);
  Py_DECREF(temp);
}

%typemap(in, numinputs=0) sofus::FocusLine<DATA_TYPE> **pFocusLine (sofus::FocusLine<DATA_TYPE> *temp) {
  $1 = &temp;
}

%typemap(argout) sofus::FocusLine<DATA_TYPE> ** {
  PyObject* temp = NULL;
  if (!PyList_Check($result)) {
    temp = $result;
    $result = PyList_New(1);
    PyList_SetItem($result, 0, temp);
  }

  // Create shadow object (do not use SWIG_POINTER_NEW)
  temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1),
			    $descriptor(sofus::FocusLine<DATA_TYPE>*),
			    SWIG_POINTER_OWN | 0);

  PyList_Append($result, temp);
  Py_DECREF(temp);
}

%typemap(in, numinputs=0) sofus::FocusLineList<DATA_TYPE> **obj (sofus::FocusLineList<DATA_TYPE> *temp) {
  $1 = &temp;
}

%typemap(argout) sofus::FocusLineList<DATA_TYPE> ** {
  PyObject* temp = NULL;
  if (!PyList_Check($result)) {
    temp = $result;
    $result = PyList_New(1);
    PyList_SetItem($result, 0, temp);
  }
  // Create shadow object (do not use SWIG_POINTER_NEW)
  temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1),
			    $descriptor(sofus::FocusLineList<DATA_TYPE>*),
			    SWIG_POINTER_OWN | 0);

  PyList_Append($result, temp);
  Py_DECREF(temp);
}

// Ensures no double deletion
%delobject sofus::FocusLineListDestroy<DATA_TYPE>;

%extend sofus::FocusLineList<DATA_TYPE> {
  const sofus::FocusLine<DATA_TYPE>& __getitem__(int index) const {
    return (*($self))[index];
  }
  int __len__() const {
    return (int) ($self)->count;
  }
// Causes double definition
//  ~FocusLineList<DATA_TYPE>() {
//    sofus::FocusLineListDestroy<DATA_TYPE>($self);
//  }
}

%enddef    /* %fnm_typemaps() */

/* %fnm_extensions() macro
 *
 * Macro for extending the Python classes Aperture<DATA_TYPE>
 */
%define %fnm_extensions(DATA_TYPE)

%extend sofus::FocusLine<DATA_TYPE> {
%pythoncode %{
# Implement __dir__() to return it plus all of the other members
def __dir__(self):
  if version_info >= (3, 3, 0):
    return FocusLineFloat.__swig_dir__
  else:
    return self.__dict__.keys() + FocusLineFloat.__swig_dir__ + FocusLineFloat.__dict__.keys()
%}
};

%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, fs, FsGet, FsSet);
%attribute(fnm::Aperture<DATA_TYPE>, int, focus_type, FocusingTypeGet, FocusingTypeSet);
%attribute(fnm::Aperture<DATA_TYPE>, size_t, nthreads, NThreadsGet, NThreadsSet);

%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, alpha, AlphaGet, AlphaSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, beta, BetaGet, BetaSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, w, WGet, WSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, f0, F0Get, F0Set);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, rho, DensityGet, DensitySet);

%attribute(fnm::Aperture<DATA_TYPE>, size_t, nDivH, NDivHGet, NDivHSet);
%attribute(fnm::Aperture<DATA_TYPE>, size_t, nDivW, NDivWGet, NDivWSet);
%attribute(fnm::Aperture<DATA_TYPE>, bool, att_enabled, AttenuationEnabledGet, AttenuationEnabledSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, bandwidth, BandWidthGet, BandWidthSet);
%attribute(fnm::Aperture<DATA_TYPE>, size_t, nelements, NElementsGet);

%attribute(fnm::Aperture<DATA_TYPE>, int, excitation_type, ExcitationTypeGet, ExcitationTypeSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, fc, FCGet, FCSet);
%attribute(fnm::Aperture<DATA_TYPE>, DATA_TYPE, c, CGet, CSet);

%attribute(fnm::Aperture<DATA_TYPE>, int, apodization_type, ApodizationTypeGet, ApodizationTypeSet);

#if FNM_PULSED_WAVE
  %attribute(fnm::Aperture<DATA_TYPE>, bool, normalize, NormalizeGet, NormalizeSet);
  %attribute(fnm::Aperture<DATA_TYPE>, int, impulse_type, ImpulseTypeGet, ImpulseTypeSet);
#endif

%extend fnm::Aperture<DATA_TYPE> {
%pythoncode %{

focus = property(_swig_fnm.ApertureFloat_FocusGet, _swig_fnm.ApertureFloat_FocusSet)
center_focus = property(_swig_fnm.ApertureFloat_CenterFocusGet, _swig_fnm.ApertureFloat_CenterFocusSet)
elements = property(_swig_fnm.ApertureFloat_ElementsGet, _swig_fnm.ApertureFloat_ElementsSet)
subelements = property(_swig_fnm.ApertureFloat_SubElementsGet, _swig_fnm.ApertureFloat_SubElementsSet)
pos = property(_swig_fnm.ApertureFloat_PositionsGet, _swig_fnm.ApertureFloat_PositionsSet)
apodization = property(_swig_fnm.ApertureFloat_ApodizationGet, _swig_fnm.ApertureFloat_ApodizationSet)
delays = property(_swig_fnm.ApertureFloat_DelaysGet, _swig_fnm.ApertureFloat_DelaysSet)

if 'FNM_PULSED_WAVE' in globals():
  focus2 = property(_swig_fnm.ApertureFloat_Focus2Get, _swig_fnm.ApertureFloat_Focus2Set)
  excitation = property(_swig_fnm.ApertureFloat_ExcitationGet, _swig_fnm.ApertureFloat_ExcitationSet)
  impulse = property(_swig_fnm.ApertureFloat_ImpulseGet, _swig_fnm.ApertureFloat_ImpulseSet)
  sysparm = property(_swig_fnm.ApertureFloat_SysParmGet, _swig_fnm.ApertureFloat_SysParmSet)

# Read-only properties
phases = property(_swig_fnm.ApertureFloat_PhasesGet)
subphases = property(_swig_fnm.ApertureFloat_SubPhasesGet)
rectangles =property(_swig_fnm.ApertureFloat_RectanglesGet)
extent =property(_swig_fnm.ApertureFloat_ExtentGet)
area =property(_swig_fnm.ApertureFloat_AreaGet)

#Implement __dir__() to return it plus all of the other members
def __dir__(self):
  if version_info >= (3, 3, 0):
     return list(ApertureFloat.__dict__.keys())
  else:
    return self.__dict__.keys() + ApertureFloat.__swig_dir__ + ApertureFloat.__dict__.keys()

if 'FNM_CLOSURE_FUNCTIONS' in globals():
  # An error was introduced in SWIG 3.0.10
  def RwFloatParamSetMulti(self, *args): return _swig_fnm.ApertureFloat_RwFloatParamSet(self, *args)

def show_rect(self, ax, corners, color, trans):
  if (trans):
    rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)],alpha=0.1)
  else:
    rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)],alpha=0.9)
  if color == None:
    rect.set_color(colors.rgb2hex([0.60,0.05,0.65]))
  else:
    rect.set_color(color)
  rect.set_edgecolor('k')
  ax.add_collection3d(rect)
  return rect

def show_box(self,**kwargs):
  # TODO: Introduce function in fnm_data returning rectangles of box

  opt = dotdict({'ax' : None,
                 'color' : None})
  opt.update(**kwargs)
  ax = opt.ax
  if (ax==None):
    if len(plt.get_fignums()) > 0:
      ax = plt.gca()
      if (type(ax).name != '3d'):
        print('Current axis is not 3d')
        return None
    else:
      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      self.show(ax=ax)

  e = self.extent
  _min = e[:,0]
  _max = e[:,1]

  allcorners = np.array(list(itertools.product(*zip(_min,_max))))

  # Positive orientation
  orient = [1, 0, 2, 1]

  for iXYZ in range(e.shape[0]):
    for iMinMax in range(e.shape[1]):
      indices = np.where(allcorners[:,iXYZ] == e[iXYZ,iMinMax])[0]
      corners = sorted(allcorners[indices,:],
                       key=lambda x: np.arctan2(x[orient[2-iXYZ]],
                                                x[orient[3-iXYZ]]))
      self.show_rect(ax, corners, None, True)
  plt.draw()
def rotate(self, euler, point):
  subelements = self.subelements.copy()
  rotm = euler2rot(euler[0], euler[1], euler[2], conv='yxy')
  for i in range(subelements.shape[0]):
    for j in range(subelements.shape[1]):
      rp = subelements[i,j,2:5] - point
      angles = subelements[i,j,5:8]
      rotm0 = euler2rot(angles[0],angles[1],angles[2],conv='yxy')
      crot = np.dot(rotm,rotm0)
      subelements[i,j,5:8] = rot2euler(crot,conv='yxy')
      subelements[i,j,2:5] = point + np.dot(rotm, rp)

  self.SubElementsSet(subelements)
def show(self,**kwargs):
  cmap = cm.get_cmap('rainbow')
  opt = dotdict({'ax' : None,
                 'color' : None})
  opt.update(**kwargs)
  ax = opt.ax
  if (ax==None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
  nElements = self.NElementsGet()
  nSub = self.NSubElementsGet()

  _rectangles = np.reshape(self.rectangles,(nElements,nSub,4,3))
  _min = [1,1,1]
  _max = [0,0,0]
  e = []

  for i in range(nElements):
    for j in range(nSub):
      corners = _rectangles[i,j,:,:]
      _min = np.minimum(np.min(corners,axis=0),_min)
      _max = np.maximum(np.max(corners,axis=0),_max)
      e.append(self.show_rect(ax, corners, cmap(self.apodization[i]), None))

  _min = _min - np.finfo(np.float32).eps
  _max = _max + np.finfo(np.float32).eps
  ax.set_xlim([_min[0],_max[0]])
  ax.set_ylim([_min[1],_max[1]])
  ax.set_zlim([_min[2],_max[2]])
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  ax.set_aspect('auto')
  return nargout(e)
def __getstate__(self):
  args = (self.FsGet(),
          self.CGet(),
          self.SubElementsGet(),
          self.ApodizationGet(),
          self.ExcitationGet(),
          self.ImpulseGet(),
          self.FocusGet(),
          self.CenterFocusGet(),
          self.DelaysGet())
  return args
def __setstate__(self, state):
  self.__init__()
  (fs, c, subelements, apodization, excitation, impulse, focus, center_focus, delays) = state
  self.FsSet(fs)
  self.CSet(c)
  self.SubElementsSet(subelements)
  self.ApodizationSet(apodization)
  if 'FNM_PULSED_WAVE' in globals():
    self.ExcitationSet(excitation)
  self.ImpulseSet(impulse)
  self.FocusSet(focus)
  self.CenterFocusSet(center_focus)
  self.DelaysSet(delays)
%}
};
%enddef    /* %fnm_extensions() */

/* %fnm_circular() macro
 *
 * Macro for extending the Python classes CircularAperture<DATA_TYPE>
 */
%define %fnm_circular(DATA_TYPE)

%extend fnm::CircularAperture<DATA_TYPE> {
%pythoncode %{
def __dir__(self):
  if version_info >= (3, 3, 0):
    return CircularApertureFloat.__swig_dir__ + list(CircularApertureFloat.__dict__.keys())
  else:
    return self.__dict__.keys() + CircularApertureFloat.__swig_dir__ + CircularApertureFloat.__dict__.keys()
%}

};
%enddef    /* %fnm_circular() */

/* Local variables: */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* End: */
