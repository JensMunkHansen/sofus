// TODO: Crazy constant error when including export files

#define _SINGLE_LIBRARY
%module(docstring="This is a wrapper for FNM", directors="1") swig_fnm
#pragma SWIG nowarn=320

#if 0
  // Experimental stuff for setting structures 
  %pythoncode %{
  def StructArgs(type_name):
    def wrap(f):
      def _wrapper(*args, **kwargs):
        ty=globals()[type_name]
        arg=(ty(),) if kwargs else tuple()
        for it in kwargs.iteritems():
          setattr(arg[0], *it)
        return f(*(args+arg))
      return _wrapper
    return wrap
  %}
  
  %define %StructArgs(func, ret, type)
  %pythoncode %{ @StructArgs(#type) %} // *very* position sensitive
  %pythonprepend func %{ %} // Hack to workaround problem with #3
  ret func(const type*);
  %ignore func;
  %enddef
#endif

%{

  #define SWIG_FILE_WITH_INIT
  #include <sps/progress.hpp>
  #include <fnm/config.h>
  #include <sps/cenv.h>
  #include <fnm/fnm_types.hpp>
  #include <fnm/fnm.hpp>
  #include <gl/gl.hpp>
%}

%include "windows.i"

%include "carrays.i"
%array_class(size_t, sizetArrayClass);

%array_class(float, floatArrayClass);
%array_functions(float, floatArray);

%array_functions(double, doubleArray);
%array_class(double, doubleArrayClass);

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

// TODO: Consider including cenv.h
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
 // %fnm_typemaps(double)
#endif

%apply (std::complex<float>** ARGOUTVIEWM_ARRAY1, size_t* DIM1) {(std::complex<float>** outTest, size_t* nOutTest)};

namespace sps {
  %feature("director") ProgressBarInterface;
}
%include <sps/progress.hpp>

// Individual modules
%include <fnm/config.h>
%include <sps/cenv.h>
%include <fnm/fnm_export.h>
%include <gl/gl_export.h>

%rename(FocusingType) FocusingTypeNS;

%include <fnm/fnm_types.hpp>
%include <fnm/fnm.hpp>
%include <gl/gl.hpp>

namespace fnm {
  %template (ApertureFloat) Aperture<float>;
  %template (SysParmFloat) sysparm_t<float>;
#ifdef FNM_DOUBLE_SUPPORT
   //  %template (ApertureDouble) Aperture<double>;
   //  %template (SysParmDouble) sysparm_t<double>;
#endif

  // %StructArgs(ApertureFloat::SysParmSet, void, SysParmFloat) // Not working must include definition right after
}

#ifdef SWIGPYTHON
  // Include matplotlib for display of apertures 
  %pythoncode %{
    import matplotlib.pyplot as plt
    plt.ion()
  %}
  
  %fnm_extensions(float)
# ifdef FNM_DOUBLE_SUPPORT
   //  %fnm_extensions(double)
# endif
#endif



  
// Instantiate templates

/* Local variables: */
/* indent-tab-mode: nil */
/* tab-width: 2 */
/* c-basic-offset: 2 */
/* indent-tabs-mode: nil */
/* End: */
