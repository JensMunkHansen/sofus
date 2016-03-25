#define _SINGLE_LIBRARY
%module(docstring="This is a Python wrapper for Sofus") swig_fnm
#pragma SWIG nowarn=320
%{

  #define SWIG_FILE_WITH_INIT
  // Try to ignore cast between pointer-to-function and pointer-to-object
  #include "swig_system.h"
  #include "FnmMath.hpp"
  #include "fnm.hpp"
%}

#ifdef _SWIG_WIN32
%include "windows.i"
#endif

// %feature ("flatnested");

%include "carrays.i"
%array_class(float, floatArrayClass);
%array_class(double, doubleArrayClass);
%array_class(size_t, sizetArrayClass);
%array_functions(float, floatArray);
%array_functions(double, doubleArray);

%array_functions(unsigned char, uint8Array);
%array_class(unsigned char, uint8ArrayClass);

%array_functions(int, int32Array);
%array_class(int, int32ArrayClass);

%include "std_complex.i"

%include "numpy.i"
%include "numpy_ex.i"

%init {
  import_array();
}

//%include "numpy_std_complex.i"

#define STATIC_INLINE_BEGIN
#define STATIC_INLINE_END

// Not defined unless sps/export.h is included
#define SPS_EXPORT
#define FNM_EXPORT

%apply (float* IN_ARRAY2, int DIM1, int DIM2) {(const float* pos, const size_t nPositions, const size_t nDim)};
%apply (std::complex<float>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(std::complex<float>** odata, size_t* nOutPositions)};

%apply (float ARGOUT_ARRAY1[ANY]) {(float oFocus[3])}
%apply (float IN_ARRAY1[ANY])     {(const float iFocus[3])}

%apply (float** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(float** phases, size_t* nPhases)}
%apply (float** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) {(float** pos, size_t* nElements, size_t* nParams)}

%apply (float** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) {(float** pos, size_t* nElements, size_t* nSubElements, size_t* nParams)}
%apply (float** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(float** outRectangles, size_t* nElements, size_t* nSubElements, size_t* nCornerCoordinates)}


%apply (double* IN_ARRAY2, int DIM1, int DIM2) {(const double* pos, const size_t nPositions, const size_t nDim)};
%apply (std::complex<double>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(std::complex<double>** odata, size_t* nOutPositions)};

%apply (double ARGOUT_ARRAY1[ANY]) {(double oFocus[3])}
%apply (double IN_ARRAY1[ANY])     {(const double iFocus[3])}

%apply (double** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(double** phases, size_t* nPhases)}
%apply (double** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) {(double** pos, size_t* nElements, size_t* nParams)}

%apply (double** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) {(double** pos, size_t* nElements, size_t* nSubElements, size_t* nParams)}
%apply (double** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(double** outRectangles, size_t* nElements, size_t* nSubElements, size_t* nCornerCoordinates)}

// Individual modules
%include "fnm.hpp"

namespace fnm {
  %template (ApertureFloat) Aperture<float>;
  %template (ApertureDouble) Aperture<double>;
}



%clear (float oFocus[3]);
%clear (const float iFocus[3]);

%extend fnm::Aperture<float> {
  %pythoncode %{
    # Read-write properties
    __swig_getmethods__["f0"]           = F0Get
    __swig_setmethods__["f0"]           = F0Set
    __swig_getmethods__["c"]            = CGet
    __swig_setmethods__["c"]            = CSet
    __swig_getmethods__["nDivW"]         = NDivWGet
    __swig_setmethods__["nDivW"]         = NDivWSet
    __swig_getmethods__["nDivH"]         = NDivHGet
    __swig_setmethods__["nDivH"]         = NDivHSet
    __swig_getmethods__["focus"]        = FocusGet
    __swig_setmethods__["focus"]        = FocusSet
    __swig_getmethods__["nthreads"]     = NThreadsGet
    __swig_setmethods__["nthreads"]     = NThreadsSet

    # Read-only properties
    __swig_getmethods__["pos"]           = PositionsGet
    __swig_getmethods__["phases"]        = PhasesGet
    __swig_getmethods__["rectangles"]   = RectanglesGet
    if _newclass:
        # Read-write properties
        f0    = property(F0Get, F0Set)
        c     = property(CGet, CSet)
        nDivW    = property(NDivWGet, NDivWSet)
        nDivH    = property(NDivHGet, NDivHSet)
        focus    = property(FocusGet, FocusSet)
        nthreads = property(NThreadsGet, NThreadsSet)

        rectangles = property(RectanglesGet)
          
    # Merge the two method dictionaries, and get the keys
    __swig_dir__ = dict(__swig_getmethods__.items() + __swig_setmethods__.items()).keys()
    # Implement __dir__() to return it plus all of the other members
    def __dir__(self):
      return self.__dict__.keys() + ApertureFloat.__swig_dir__

    def show_rect(self,ax,corners):
      from mpl_toolkits.mplot3d import art3d
      import matplotlib.colors as colors
      import numpy as np
      rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)])
      rect.set_color(colors.rgb2hex([0.60,0.05,0.65]))
      rect.set_edgecolor('k')
      ax.add_collection3d(rect)

    def show(self):
      import numpy as np
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import art3d

      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      nElements = self.NElementsGet()
      nSub = self.NSubElementsGet()

      bum = np.reshape(self.rectangles,(nElements,nSub,4,3))
      _min = [1,1,1]
      _max = [0,0,0]

      for i in range(nElements):
        for j in range(nSub):
          corners = bum[i,j,:,:]
          corners = corners[[0,1,3,2],:]
          _min = np.minimum(np.min(corners,axis=0),_min)
          _max = np.maximum(np.max(corners,axis=0),_max)
          self.show_rect(ax,corners)

      ax.set_xlim([_min[0],_max[0]])
      ax.set_ylim([_min[1],_max[1]])
      ax.set_zlim([_min[2],_max[2]])
      ax.set_aspect('auto')
  %}
   
};

%extend fnm::Aperture<double> {
  %pythoncode %{
    # Read-write properties
    __swig_getmethods__["f0"]           = F0Get
    __swig_setmethods__["f0"]           = F0Set
    __swig_getmethods__["c"]            = CGet
    __swig_setmethods__["c"]            = CSet
    __swig_getmethods__["nDivW"]         = NDivWGet
    __swig_setmethods__["nDivW"]         = NDivWSet
    __swig_getmethods__["nDivH"]         = NDivHGet
    __swig_setmethods__["nDivH"]         = NDivHSet
    __swig_getmethods__["focus"]        = FocusGet
    __swig_setmethods__["focus"]        = FocusSet
    __swig_getmethods__["nthreads"]     = NThreadsGet
    __swig_setmethods__["nthreads"]     = NThreadsSet

    # Read-only properties
    __swig_getmethods__["positions"]        = PositionsGet
    __swig_getmethods__["phases"]        = PhasesGet
    if _newclass:
        # Read-write properties
        f0    = property(F0Get, F0Set)
        c     = property(CGet, CSet)
        nDivW    = property(NDivWGet, NDivWSet)
        nDivH    = property(NDivHGet, NDivHSet)
        focus    = property(FocusGet, FocusSet)
        nthreads = property(NThreadsGet, NThreadsSet)
    # Merge the two method dictionaries, and get the keys
    __swig_dir__ = dict(__swig_getmethods__.items() + __swig_setmethods__.items()).keys()
    # Implement __dir__() to return it plus all of the other members
    def __dir__(self):
      return self.__dict__.keys() + ApertureDouble.__swig_dir__

    def show_rect(self,ax,corners):
      from mpl_toolkits.mplot3d import art3d
      import matplotlib.colors as colors
      import numpy as np
      rect = art3d.Poly3DCollection([np.roll(corners,-2,axis=0)])
      rect.set_color(colors.rgb2hex([0.60,0.05,0.65]))
      rect.set_edgecolor('k')
      ax.add_collection3d(rect)

    def show(self):
      import numpy as np
      import matplotlib.pyplot as plt
      from mpl_toolkits.mplot3d import art3d

      fig = plt.figure()
      ax = fig.add_subplot(111, projection='3d')
      nElements = self.NElementsGet()
      nSub = self.NSubElementsGet()

      bum = np.reshape(self.rectangles,(nElements,nSub,4,3))
      _min = [1,1,1]
      _max = [0,0,0]

      for i in range(nElements):
        for j in range(nSub):
          corners = bum[i,j,:,:]
          corners = corners[[0,1,3,2],:]
          _min = np.minimum(np.min(corners,axis=0),_min)
          _max = np.maximum(np.max(corners,axis=0),_max)
          self.show_rect(ax,corners)

      ax.set_xlim([_min[0],_max[0]])
      ax.set_ylim([_min[1],_max[1]])
      ax.set_zlim([_min[2],_max[2]])
      ax.set_aspect('auto')
  %}
};


// Instantiate templates

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
