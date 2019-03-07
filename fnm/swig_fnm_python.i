%insert("python") %{
from mpl_toolkits.mplot3d import art3d
from sys import version_info
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import itertools
import numpy as np

from dicts import dotdict

plt.ion()
%}


/* %fnm_typemaps() macro
 *
 * This macro applies a long list of typemaps for templated functions. It is meant
 * to be executed for float and double
 *
 */
%define %fnm_typemaps(DATA_TYPE)

%apply (DATA_TYPE ARGOUT_ARRAY1[ANY]) {(DATA_TYPE oFocus[3])}
%apply (DATA_TYPE IN_ARRAY1[ANY])     {(const DATA_TYPE iFocus[3])}

%apply (DATA_TYPE* INPLACE_ARRAY1, size_t DIM1) \
{(DATA_TYPE* ioData, size_t nIOdata)}

%apply (DATA_TYPE* IN_ARRAY1, int DIM1) \
{(const DATA_TYPE* data, const size_t nData)}

%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* pos, const size_t nPositions, const size_t nDim)};
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* iData, const size_t nChannels, const size_t nSamples)};
%apply (DATA_TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3) \
{(const DATA_TYPE* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim)};

%apply(DATA_TYPE** ARGOUTVIEW_ARRAY1, size_t* DIM1) {(DATA_TYPE** data, size_t* nData)}

%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(DATA_TYPE** odata, size_t* nData)}

%apply (int** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(int** offsets, size_t* nOffsets)}

%apply (int** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(int** lengths, size_t* nLengths)}

%apply (int** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(int** sizeAndOffsets, size_t* nFilters, size_t* nTwo)}

%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** out, size_t* nElements, size_t* nParams)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** coordinates, size_t* nDim, size_t* nLimits)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** odata, size_t* nSignals, size_t* nSamples)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(DATA_TYPE** out, size_t* nElements, size_t* nSubElements, size_t* nParams)}

%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(std::complex<DATA_TYPE>** odata, size_t* nOutPositions)};

%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(std::complex<DATA_TYPE>** odata, size_t* nOutPositions, size_t* nSamples)};


// TEST (for angular spectrum approach)
%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY3, size_t* DIM1,
 size_t* DIM2, size_t* DIM3) {(std::complex<DATA_TYPE>** p1, size_t*
 onx, size_t* ony, size_t* onz)};
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* p0, const size_t nx, const size_t ny)};

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

%enddef    /* %fnm_typemaps() */

/* %fnm_extensions() macro
 *
 * Macro for extending the Python classes Aperture<DATA_TYPE>
 */
%define %fnm_extensions(DATA_TYPE)

%extend fnm::Aperture<DATA_TYPE> {
%pythoncode %{
# Read-write properties
__swig_getmethods__["f0"]           = F0Get
__swig_setmethods__["f0"]           = F0Set
__swig_getmethods__["c"]            = CGet
__swig_setmethods__["c"]            = CSet
__swig_getmethods__["nDivW"]        = NDivWGet
__swig_setmethods__["nDivW"]        = NDivWSet
__swig_getmethods__["nDivH"]        = NDivHGet
__swig_setmethods__["nDivH"]        = NDivHSet
__swig_getmethods__["focus"]        = FocusGet
__swig_setmethods__["focus"]        = FocusSet
__swig_getmethods__["center_focus"] = CenterFocusGet
__swig_setmethods__["center_focus"] = CenterFocusSet
__swig_getmethods__["focus_type"]   = FocusingTypeGet
__swig_setmethods__["focus_type"]   = FocusingTypeSet
__swig_getmethods__["att_enabled"]  = AttenuationEnabledGet
__swig_setmethods__["att_enabled"]  = AttenuationEnabledSet
__swig_getmethods__["alpha"]        = AlphaGet
__swig_setmethods__["alpha"]        = AlphaSet
__swig_getmethods__["beta"]         = BetaGet
__swig_setmethods__["beta"]         = BetaSet
__swig_getmethods__["w"]            = WGet
__swig_setmethods__["w"]            = WSet

__swig_getmethods__["nthreads"]     = NThreadsGet
__swig_setmethods__["nthreads"]     = NThreadsSet

__swig_getmethods__["subelements"]  = SubElementsGet
__swig_setmethods__["subelements"]  = SubElementsSet
__swig_getmethods__["elements"]     = ElementsGet
__swig_setmethods__["elements"]     = ElementsSet
__swig_getmethods__["pos"]          = PositionsGet
__swig_setmethods__["pos"]          = PositionsSet
__swig_getmethods__["apodization"]  = ApodizationGet
__swig_setmethods__["apodization"]  = ApodizationSet
__swig_getmethods__["delays"]       = DelaysGet
__swig_setmethods__["delays"]       = DelaysSet

# Branch in python on presence of pulsed-waves
if 'FNM_PULSED_WAVE' in globals():
  __swig_getmethods__["fc"]           = FCGet
  __swig_setmethods__["fc"]           = FCSet
  __swig_getmethods__["impulse_type"] = ImpulseTypeGet
  __swig_setmethods__["impulse_type"] = ImpulseTypeSet
  __swig_getmethods__["excitation_type"] = ExcitationTypeGet
  __swig_setmethods__["excitation_type"] = ExcitationTypeSet
  __swig_getmethods__["fs"]           = FsGet
  __swig_setmethods__["fs"]           = FsSet
  __swig_getmethods__["normalize"]    = NormalizeGet
  __swig_setmethods__["normalize"]    = NormalizeSet
  __swig_getmethods__["excitation"]   = ExcitationGet
  __swig_setmethods__["excitation"]   = ExcitationSet
  __swig_getmethods__["impulse"]      = ImpulseGet
  __swig_setmethods__["impulse"]      = ImpulseSet
  __swig_getmethods__["sysparm"]      = SysParmGet
  __swig_setmethods__["sysparm"]      = SysParmSet
  __swig_getmethods__["bandwidth"]    = BandWidthGet
  __swig_setmethods__["bandwidth"]    = BandWidthSet
  

# Read-only properties
__swig_getmethods__["phases"]       = PhasesGet
__swig_getmethods__["subphases"]    = SubPhasesGet
__swig_getmethods__["nelements"]    = NElementsGet
__swig_getmethods__["rectangles"]   = RectanglesGet
__swig_getmethods__["extent"]       = ExtentGet
__swig_getmethods__["area"]         = AreaGet

if _newclass:
    # Read-write properties
  f0       = property(F0Get, F0Set)

  if 'FNM_PULSED_WAVE' in globals():
    fs       = property(FsGet, FsSet)

  c        = property(CGet, CSet)
  w        = property(WGet, WSet)
  nDivW    = property(NDivWGet, NDivWSet)
  nDivH    = property(NDivHGet, NDivHSet)
  nthreads = property(NThreadsGet, NThreadsSet)
    # Not necessary to keep as class property
  focus    = property(FocusGet, FocusSet)
  center_focus = property(CenterFocusGet, CenterFocusSet)
  focus_type = property(FocusingTypeGet, FocusingTypeSet)
  rectangles = property(RectanglesGet)
  extent     = property(ExtentGet)
  area       = property(AreaGet)
      
  # Merge the two method dictionaries, and get the keys
  # __swig_dir__ = dict(__swig_getmethods__.items() + __swig_setmethods__.items()).keys() # Works in 2.7
  # __swig_dir__ = list({**__swig_getmethods__, **__swig_setmethods__}.keys())            # Works in 3.5+ (syntax error in 2.7)
  __swig_dir__ = __swig_getmethods__.copy()
  __swig_dir__.update(__swig_setmethods__)
  __swig_dir__ = list(__swig_dir__.keys())

#Implement __dir__() to return it plus all of the other members
def __dir__(self):
  if version_info >= (3, 3, 0):
    return ApertureFloat.__swig_dir__
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

def show(self,**kwargs):
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

  for i in range(nElements):
    for j in range(nSub):
      corners = _rectangles[i,j,:,:]
      _min = np.minimum(np.min(corners,axis=0),_min)
      _max = np.maximum(np.max(corners,axis=0),_max)
      self.show_rect(ax,corners,opt.color,None)

  _min = _min - np.finfo(np.float32).eps
  _max = _max + np.finfo(np.float32).eps
  ax.set_xlim([_min[0],_max[0]])
  ax.set_ylim([_min[1],_max[1]])
  ax.set_zlim([_min[2],_max[2]])
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  ax.set_aspect('auto')
  return ax
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
# Read-write properties
__swig_getmethods__["f0"]           = F0Get
__swig_setmethods__["f0"]           = F0Set
__swig_getmethods__["fs"]           = FsGet
__swig_setmethods__["radius"]       = RadiusSet
__swig_getmethods__["radius"]       = RadiusGet
__swig_setmethods__["fs"]           = FsSet
__swig_getmethods__["w"]            = WGet
__swig_setmethods__["w"]            = WSet
__swig_getmethods__["nDivA"]        = NMaxDivAGet
__swig_setmethods__["nDivA"]        = NMaxDivASet
__swig_getmethods__["gridSectorScale"] = GridSectorScaleGet
__swig_setmethods__["gridSectorScale"] = GridSectorScaleSet
if _newclass:
    # Read-write properties
  f0       = property(F0Get, F0Set)
  nDivA    = property(NMaxDivAGet, NMaxDivASet)
  radius   = property(RadiusGet, RadiusSet)
  gridSectorScale = property(GridSectorScaleGet, GridSectorScaleSet)
%}
};
%enddef    /* %fnm_circular() */


/* Local variables: */
/* mode: text */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* End: */
