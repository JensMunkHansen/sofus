// TODO: Crazy constant error when including export files

#define _SINGLE_LIBRARY
%module(docstring="This is a wrapper for FNM") swig_fnm
#pragma SWIG nowarn=320
%{

  #define SWIG_FILE_WITH_INIT
  #include "swig_system.h"
  #include <fnm/FnmMath.hpp>
  #include <fnm/fnm.hpp>

%}

%include "windows.i"

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

#ifdef SWIGPYTHON
%include "numpy.i"

%init {
  import_array();
}
#endif

// _WIN32 is not defined, when running SWIG on Windows
#define STATIC_INLINE_BEGIN
#define STATIC_INLINE_END

#define ALIGN16_BEGIN
#define ALIGN16_END

// Define FNM_EXPORT if you skip DEFINE_NO_DEPRECATED for GenerateExportHeader
   
#ifdef SWIG_INCLUDE_DOCUMENTATION
%import "documentation.i"
#endif

%include "typemaps.i"



#ifdef SWIGPYTHON

/* %fnm_typemaps() macro
 *
 * This macro applies a long list of typemaps for templated functions. It is meant
 * to be executed for float and double
 *
 */
%define %fnm_typemaps(DATA_TYPE)

%apply (DATA_TYPE ARGOUT_ARRAY1[ANY]) {(DATA_TYPE oFocus[3])}
%apply (DATA_TYPE IN_ARRAY1[ANY])     {(const DATA_TYPE iFocus[3])}

%apply (DATA_TYPE* IN_ARRAY2, int DIM1, int DIM2) \
{(const DATA_TYPE* pos, const size_t nPositions, const size_t nDim)};
%apply (DATA_TYPE* IN_ARRAY3, int DIM1, int DIM2, int DIM3) \
{(const DATA_TYPE* pos, const size_t nElements, const size_t nSubElementsPerElement, const size_t nDim)};

%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY1, size_t* DIM1) \
{(DATA_TYPE** phases, size_t* nPhases)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY2, size_t* DIM1, size_t* DIM2) \
{(DATA_TYPE** out, size_t* nElements, size_t* nParams)}
%apply (DATA_TYPE** ARGOUTVIEWM_ARRAY3, size_t* DIM1, size_t* DIM2, size_t* DIM3) \
{(DATA_TYPE** out, size_t* nElements, size_t* nSubElements, size_t* nParams)}

%apply (std::complex<DATA_TYPE>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(std::complex<DATA_TYPE>** odata, size_t* nOutPositions)};

%enddef    /* %fnm_typemaps() */

/* %fnm_extensions() macro
 *
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

__swig_getmethods__["nthreads"]     = NThreadsGet
__swig_setmethods__["nthreads"]     = NThreadsSet

__swig_getmethods__["subelements"]  = SubElementsGet
__swig_setmethods__["subelements"]  = SubElementsSet
__swig_getmethods__["elements"]     = ElementsGet
__swig_setmethods__["elements"]     = ElementsSet

# Read-only properties
__swig_getmethods__["pos"]          = PositionsGet
__swig_getmethods__["phases"]       = PhasesGet
__swig_getmethods__["nelements"]    = NElementsGet
__swig_getmethods__["rectangles"]   = RectanglesGet

if _newclass:
    # Read-write properties
    f0       = property(F0Get, F0Set)
    c        = property(CGet, CSet)
    nDivW    = property(NDivWGet, NDivWSet)
    nDivH    = property(NDivHGet, NDivHSet)
    nthreads = property(NThreadsGet, NThreadsSet)
    # Not necessary to keep as class property
    focus    = property(FocusGet, FocusSet)
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

  _min = _min - np.finfo(np.float32).eps
  _max = _max + np.finfo(np.float32).eps
  ax.set_xlim([_min[0],_max[0]])
  ax.set_ylim([_min[1],_max[1]])
  ax.set_zlim([_min[2],_max[2]])
  ax.set_aspect('auto')
%}
};
%enddef    /* %fnm_extensions() */

%fnm_typemaps(float)
%fnm_typemaps(double)

#endif


%typemap(in,numinputs=0) (int* COUNT, Zone **OUTPUT) (int dim, Zone* temp) %{
  $1 = &dim;
  $2 = &temp;
%}

// Build a list of tuples containing the two entries from the struct.
// Append the new Python list object to the existing "int" result.
%typemap(argout) (int* COUNT, Zone **OUTPUT) {
  PyObject* out = NULL;
  if (!PyList_Check($result)) {
    out = $result;
    $result = PyList_New(1);
    PyList_SetItem($result, 0, out);
  }

  printf("count: %d\n", (*$1));
  int i = 0;
  for (i = 0 ; i < (*$1) ; i++)
  {
    // No shadow dtor never called SWIG_POINTER_NOSHADOW
    PyObject* s = SWIG_NewPointerObj(*$2+i, $descriptor(Zone*), SWIG_POINTER_NEW); // SWIG_POINTER_OWN
    PyList_Append($result,s);
    Py_DECREF(s);
  }
  free($2);
}

// Apply the OUTPUT typemap set to the "results" parameter.
%apply (int* COUNT, Zone **OUTPUT) {(int* count, Zone** results)};

/*
%delobject _Zone::~_Zone;

struct _Zone {};
%extend _Zone {
  ~_Zone() {
    printf("extend _Zone::~_Zone()\n");
    //    Zone* z = (Zone*) ($self);
    //    z.~Zone();
  }
}
%ignore _Zone;
*/

// Individual modules
%include <sps/cenv.h>
%include <fnm/fnm_export.h>
%include <fnm/FnmMath.hpp>
%include <fnm/fnm.hpp>

namespace fnm {
  %template (ApertureFloat) Aperture<float>;
  %template (ApertureDouble) Aperture<double>;
  %template (SysParmFloat) sysparm_t<float>;
}

#ifdef SWIGPYTHON
  // Include matplotlib for display of apertures 
  %pythoncode %{
    import matplotlib.pyplot as plt
    plt.ion()
  %}
  
  %fnm_extensions(float)
  %fnm_extensions(double)
#endif


// Instantiate templates

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
