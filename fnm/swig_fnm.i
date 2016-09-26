// TODO: Crazy constant error when including export files

#define _SINGLE_LIBRARY
%module(docstring="This is a wrapper for FNM") swig_fnm
#pragma SWIG nowarn=320
%{

  #define SWIG_FILE_WITH_INIT
  #include <fnm/fnm.hpp>
  #include <gl/gl.hpp>
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

// TODO: Use (cenv.h)
#define STATIC_INLINE_BEGIN
#define STATIC_INLINE_END

#define ALIGN16_BEGIN
#define ALIGN16_END

// Define exports or use DEFINE_NO_DEPRECATED when calling GenerateExportHeader
#define GL_EXPORT
#define FNM_EXPORT
   
#ifdef SWIG_INCLUDE_DOCUMENTATION
%import "documentation.i"
#endif

%include "typemaps.i"

#ifdef SWIGPYTHON
%include "swig_fnm_python.i"

%fnm_typemaps(float)
%fnm_typemaps(double)
#endif

// Individual modules
%include <sps/cenv.h>
%include <fnm/fnm_export.h>
%include <gl/gl_export.h>

%include <fnm/fnm.hpp>
%include <gl/gl.hpp>

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
