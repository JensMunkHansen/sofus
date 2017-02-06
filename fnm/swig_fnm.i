// TODO: Crazy constant error when including export files

#define _SINGLE_LIBRARY
%module(docstring="This is a wrapper for FNM", directors="1") swig_fnm // Gives trouble debugging
//%module(docstring="This is a wrapper for FNM") swig_fnm
#pragma SWIG nowarn=320

%{
  #define SWIG_FILE_WITH_INIT
  #include <sps/config.h>
  #include <sps/progress_if.hpp>
  #include <sps/progress.hpp>
  #include <fnm/config.h>
  #include <sps/cenv.h>
#if FNM_PULSED_WAVE
  #include <sofus/sofus_types.hpp>
#endif
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
#define SOFUS_EXPORT
#define SPS_EXPORT

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

%rename(FocusingType) FocusingTypeNS;
%rename(TimeDomainCalcType) TimeDomainCalcTypeNS;
%rename(PulsedWaveIntOrder) PulsedWaveIntOrderNS;

#if FNM_PULSED_WAVE
%include <sofus/sofus_types.hpp>
#endif
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
