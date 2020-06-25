#define _SINGLE_LIBRARY

%module(docstring="This is a wrapper for FNM", directors="1") swig_fnm
#pragma SWIG nowarn=320

%{
  #define SWIG_FILE_WITH_INIT
  #include <sps/config.h>
  #include <sps/progress_if.hpp>
  #include <sps/progress.hpp>
  #include <fnm/config.h>
  #include <sps/cenv.h>
#if FNM_PULSED_WAVE
  #include <sofus/sofus_types.h>
  #include <sofus/sofus_types.hpp>
#endif
  #include <fnm/fnm_export.h>
  #include <fnm/fnm_types.h>
  #include <fnm/fnm_types.hpp>
  #include <fnm/fnm.hpp>
  #include <fnm/fnm_profiling.hpp>
  #include <fnm/circular.hpp>
  #include <gl/gl_export.h>
  #include <gl/gl.hpp>
%}

%include "windows.i"

%include "carrays.i"
%array_class(size_t, sizetArrayClass);

%include "cstring.i"
%cstring_output_allocate(char **ostring, free(*$1));

%array_class(float, floatArrayClass);
%array_functions(float, floatArray);

%array_class(double, doubleArrayClass);
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

// Consider using
// %insert("header") %{
// #define STATIC_INLINE_BEGIN
// %}

#define STATIC_INLINE_BEGIN
#define STATIC_INLINE_END

#define ALIGN16_BEGIN
#define ALIGN16_END

// Define exports or use DEFINE_NO_DEPRECATED when calling GenerateExportHeader
#define GL_EXPORT
#define FNM_EXPORT
#define SOFUS_EXPORT
#define SPS_EXPORT

// This must be included before #include's
#ifdef SWIG_INCLUDE_DOCUMENTATION
  %include "documentation.i"
#endif

%include "typemaps.i"

#ifdef SWIGPYTHON
 %include "swig_fnm_python.i"
 %fnm_typemaps(float)
 %fnm_typemaps(double)
#endif

// TODO: Make work using swig 3.0.10 (works using 3.0.8) 
%typecheck(SWIG_TYPECHECK_DOUBLE_ARRAY) (const float* iMultiData)
{
  $1 = is_array($input) || PySequence_Check($input);
}

/* We could use numinputs=1, skip nDim and rely on dimensions from input array */
%typemap(in) (const float* iMultiData, size_t nDim, ...) (PyArrayObject* array=NULL, int is_new_object=0, size_t vargs[3]) {
  size_t i;
  size_t argc;
  bool singleton = true;
  size_t nDim, dim;
  nDim = 0;
  for (i = 0; i < 3; i++) vargs[i] = 0;

  /* Note: nDim is parsed as a variable argument - for a 3-dimensional array, argc = 4 */
  argc = PyTuple_Size(varargs);
  if (argc > 0) {
    nDim = PyInt_AsSsize_t(PyTuple_GetItem(varargs,0));
  }
  if (argc > 4) {
    PyErr_SetString(PyExc_ValueError, "Too many arguments");
    SWIG_fail;
  }
  for (i = 1; i < argc; i++) {
    PyObject *o = PyTuple_GetItem(varargs, i);
    if (!PyInt_Check(o)) {
      PyErr_SetString(PyExc_ValueError, "Expected an integer");
      SWIG_fail;
    }
    vargs[i-1] = PyInt_AsSsize_t(o);
    if (vargs[i-1] > 1)
      singleton = false;
  }

  array = obj_to_array_contiguous_allow_conversion($input, NPY_FLOAT, &is_new_object);
  $1 = (float*) array_data(array);

  if (singleton) {
    // Scalars and 1-dimensional arrays with length 0 or 1
    $2 = nDim;
  }
  else {
    $2 = (size_t) array_numdims(array);
    // Verify nDim matches dimension of input
    if ($2 != nDim) {
      PyErr_Format(PyExc_ValueError, "Dimensions mismatch: nDim: %zu, array_numdims: %zu", nDim, $2);
      SWIG_fail;
    }
  }
  
  // Verify each dimension from the NumPy array
  for (i=0 ; i < $2 ; i++) {
    dim = (size_t) array_dimensions(array)[i];
    if (dim != vargs[i]) {
      PyErr_Format(PyExc_ValueError, "Dimensions mismatch: varg[%zu]: %zu, dim%zu: %zu", i, vargs[i], i, dim);
      SWIG_fail;
    }
  }
  $3 = (void *)vargs;
}

%typemap(freearg)
(const float* iMultiData, size_t nDim)
{
  if (is_new_object$argnum && array$argnum)
    { Py_DECREF(array$argnum); }
}

%feature("action") fnm::Aperture<float>::RwFloatParamSet {
  size_t **vargs = (size_t **) arg5;
  result = arg1->RwFloatParamSet(arg2, arg3, arg4, vargs[0], vargs[1], vargs[2], NULL); 
}

// Works using swig 3.0.10
%typemap(in)
(float** oMultiData, size_t nDim, ...)
(float* data_temp=NULL, size_t vargs[3])
{
  int i;
  for (i = 0; i < 3; i++) vargs[i] = 0;
  
  /* Could use $symname to branch between input and */
  $1 = &data_temp;            // arg3
  $2 = PyInt_AsSsize_t(obj2); // arg4
  $3 = (void *) vargs;        // arg5

  ((void)(varargs));
}

%feature("action") fnm::Aperture<float>::RwFloatParamGet {

  size_t* vargs = (size_t*) arg5;
  result = arg1->RwFloatParamGet(arg2, arg3, arg4, &(vargs[0]), &(vargs[1]), &(vargs[2]), NULL);
  
}
 
%typemap(argout)
(float** oMultiData, size_t nDim, ...)
{
  size_t *vargs = (size_t *) arg5;
  npy_intp dims[3] = { (npy_intp) vargs[0], (npy_intp) vargs[1], (npy_intp) vargs[2] };

  PyObject* obj;

  int int_result = 0;
  if ($result != SWIG_Py_Void()) {
    int_result = PyInt_AsLong($result);
    if (int_result) {
      PyErr_Format(PyExc_TypeError,
                   "Error calling function: %s", "$symname");
      SWIG_fail;
    }
  }

  if (arg4 == 0) {
    // For scalars, an array is returned
    dims[0] = 1; 
    obj = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, (void*)(*$1));
    // return PyFloat_FromDouble(**$1);
  }
  else {
    obj = PyArray_SimpleNewFromData((int)arg4, dims, NPY_FLOAT, (void*)(*$1));
  }

  PyArrayObject* array = (PyArrayObject*) obj;

  if (!array) SWIG_fail;

%#ifdef SWIGPY_USE_CAPSULE
    PyObject* cap = PyCapsule_New((void*)(*$1), SWIGPY_CAPSULE_NAME, free_cap);
%#else
    PyObject* cap = PyCObject_FromVoidPtr((void*)(*$1), free);
%#endif
  
%#if NPY_API_VERSION < 0x00000007
  PyArray_BASE(array) = cap;
%#else
  PyArray_SetBaseObject(array,cap);
%#endif

  $result = SWIG_Python_AppendOutput($result,obj);
}

#ifdef USE_PROGRESS_BAR
namespace sps {
  %feature("director") ProgressBarInterface;
}
#endif

%include <sps/config.h>
%include <sps/progress_if.hpp>
%include <sps/progress.hpp>

// Individual modules
%include <fnm/config.h>
%include <sps/cenv.h>
%include <fnm/fnm_export.h>

#if FNM_PULSED_WAVE
%include <sofus/sofus_export.h>
#endif

%include <gl/gl_export.h>

%rename(FocusingType) FNM_FocusingTypeNS;

%rename(ImpulseType) SOFUS_ImpulseTypeNS;

%rename(ExcitationType) SOFUS_ExcitationTypeNS;

%rename(TimeDomainCalcType) SOFUS_TimeDomainCalcTypeNS;
%rename(TimeDomainIntOrder) SOFUS_TimeDomainIntOrderNS;

%rename(RwParamType) RwParamTypeNS;

%rename(SuperType) FNM_TypeNS;

#if FNM_PULSED_WAVE
  %include <sofus/sofus_types.h>
  %include <sofus/sofus_types.hpp>
#endif
%include <fnm/fnm_types.h>
%include <fnm/fnm_types.hpp>
%include <fnm/fnm.hpp>
%include <fnm/circular.hpp>

%include <gl/gl.hpp>

// Instantiate templates
namespace fnm {
  %template (ApertureFloat) Aperture<float>;
  %template (SysParmFloat) sysparm_t<float>;
  %template (CircularApertureFloat) CircularAperture<float>;
#if FNM_DOUBLE_SUPPORT
  %template (ApertureDouble) Aperture<double>;
  %template (SysParmDouble) sysparm_t<double>;
  %template (CircularApertureDouble) CircularAperture<double>;
#endif
}

#ifdef SWIGPYTHON
  // TODO: Separate typemaps from class extensions
  %fnm_extensions(float)
  %fnm_circular(float)
# if FNM_DOUBLE_SUPPORT
  %fnm_extensions(double)
  %fnm_circular(double)
# endif
#endif

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
