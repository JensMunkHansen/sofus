/* %fnm_typemaps() macro
 *
 * This macro applies a long list of typemaps for templated functions. It is meant
 * to be executed for float and double
 *
 */
%define %fnm_typemaps(DATA_TYPE)

%apply (DATA_TYPE ARGOUT_ARRAY1[ANY]) {(DATA_TYPE oFocus[3])}
%apply (DATA_TYPE IN_ARRAY1[ANY])     {(const DATA_TYPE iFocus[3])}

%apply (DATA_TYPE* IN_ARRAY1, int DIM1) \
{(const DATA_TYPE* data, const size_t nData)}
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* pos, const size_t nPositions, const size_t nDim)};
%apply (DATA_TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3) \
{(const DATA_TYPE* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim)};

%apply(DATA_TYPE** ARGOUTVIEW_ARRAY1, size_t* DIM1) {(DATA_TYPE** data, size_t* nData)}

%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(DATA_TYPE** odata, size_t* nData)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** out, size_t* nElements, size_t* nParams)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) {(DATA_TYPE** coordinates, size_t* nDim, size_t* nLimits)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) {(DATA_TYPE** odata, size_t* nSignals, size_t* nSamples)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(DATA_TYPE** out, size_t* nElements, size_t* nSubElements, size_t* nParams)}

%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(std::complex<DATA_TYPE>** odata, size_t* nOutPositions)};

// TEST
%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY3, size_t* DIM1,
 size_t* DIM2, size_t* DIM3) {(std::complex<DATA_TYPE>** p1, size_t*
 onx, size_t* ony, size_t* onz)};
%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* p0, const size_t nx, const size_t ny)};


%enddef    /* %fnm_typemaps() */

/* %fnm_extensions() macro
 *
 * Macro for extending the Python class fnm
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
__swig_getmethods__["focus_type"]   = FocusingTypeGet
__swig_setmethods__["focus_type"]   = FocusingTypeSet

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
__swig_getmethods__["fs"]           = FsGet
__swig_setmethods__["fs"]           = FsSet
__swig_getmethods__["excitation"]   = ExcitationGet
__swig_setmethods__["excitation"]   = ExcitationSet
__swig_getmethods__["impulse"]      = ImpulseGet
__swig_setmethods__["impulse"]      = ImpulseSet

# Read-only properties
__swig_getmethods__["phases"]       = PhasesGet
__swig_getmethods__["delays"]       = DelaysGet
__swig_getmethods__["nelements"]    = NElementsGet
__swig_getmethods__["rectangles"]   = RectanglesGet
__swig_getmethods__["extent"]       = ExtentGet
__swig_getmethods__["area"]         = AreaGet

if _newclass:
    # Read-write properties
    f0       = property(F0Get, F0Set)
    fs       = property(FsGet, FsSet)
    c        = property(CGet, CSet)
    nDivW    = property(NDivWGet, NDivWSet)
    nDivH    = property(NDivHGet, NDivHSet)
    nthreads = property(NThreadsGet, NThreadsSet)
    # Not necessary to keep as class property
    focus    = property(FocusGet, FocusSet)
    focus_type = property(FocusingTypeGet, FocusingTypeSet)
    rectangles = property(RectanglesGet)
    extent     = property(ExtentGet)
    area       = property(AreaGet)
      
# Merge the two method dictionaries, and get the keys
__swig_dir__ = dict(__swig_getmethods__.items() + __swig_setmethods__.items()).keys()
# Implement __dir__() to return it plus all of the other members
def __dir__(self):
  return self.__dict__.keys() + ApertureFloat.__swig_dir__

def show_rect(self,ax,corners,color):
  from mpl_toolkits.mplot3d import art3d
  import matplotlib.colors as colors
  import numpy as np
  rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)])
  if color == None:
    rect.set_color(colors.rgb2hex([0.60,0.05,0.65]))
  else:
    rect.set_color(color)
  rect.set_edgecolor('k')
  ax.add_collection3d(rect)

def show(self,**kwargs):
  import numpy as np
  import matplotlib.pyplot as plt
  from mpl_toolkits.mplot3d import art3d
  from dicts import dotdict
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
      corners = corners[[0,1,3,2],:]
      _min = np.minimum(np.min(corners,axis=0),_min)
      _max = np.maximum(np.max(corners,axis=0),_max)
      self.show_rect(ax,corners,opt.color)

  _min = _min - np.finfo(np.float32).eps
  _max = _max + np.finfo(np.float32).eps
  ax.set_xlim([_min[0],_max[0]])
  ax.set_ylim([_min[1],_max[1]])
  ax.set_zlim([_min[2],_max[2]])
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_zlabel('z')
  ax.set_aspect('auto')
%}
};
%enddef    /* %fnm_extensions() */

/* Local variables: */
/* mode: text */
/* indent-tabs-mode: nil */
/* tab-width: 2 */
/* End: */
