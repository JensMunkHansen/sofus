#ifdef SWIGPYTHON

%include "numpy.i"

%include <std_complex.i>

%numpy_typemaps(std::complex<float>,  NPY_CFLOAT , int)
%numpy_typemaps(std::complex<double>, NPY_CDOUBLE, int)

%numpy_typemaps(std::complex<float>,  NPY_CFLOAT , size_t)
%numpy_typemaps(std::complex<double>, NPY_CDOUBLE, size_t)

#endif /* SWIGPYTHON */
