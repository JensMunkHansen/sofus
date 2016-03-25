/**
 * @file   numpy_ex.i
 * @author Jens Munk Hansen <jmh@jmhlaptop.parknet.dk>
 * @date   Thu Apr  3 23:49:11 2014
 *
 * @brief  Extensions for SWIG for Numerical Python
 *
 *
 */
#ifdef SWIGPYTHON
%{
#ifndef SWIG_FILE_WITH_INIT
#  define NO_IMPORT_ARRAY
#endif
#include "stdio.h"
#include <numpy/arrayobject.h>
%}

%include "numpy.i"

/* %numpy_ex_typemaps() macro
 *
 * This macro defines a family of 5 typemaps that allow C arguments
 * of the form
 *
 * (DATA_TYPE**  DYNARGOUT_ARRAY1, DIM_TYPE* DIM1)                     // Under development
 * (DATA_TYPE*   DYNARGOUT_ARRAY2_FIXED[ANY], DIM_TYPE* DIM1)          // Under development
 * (DATA_TYPE**  DYNARGOUTVIEW_ARRAY1, DIM_TYPE* DIM1)                 // Working
 * (DATA_TYPE*   DYNARGOUTVIEW_ARRAY2_FIXED[ANY], DIM_TYPE* DIM1)      // Working
 * (DATA_TYPE*** DYNARGOUTVIEW_ARRAY2, DIM_TYPE* DIM1, DIM_TYPE* DIM2) // Working
 *
 * where "DATA_TYPE" is any type supported by the NumPy module, and
 * "DIM_TYPE" is any int-like type suitable for specifying dimensions.
 * In python, the dimensions will not need to be specified. For all
 * typemaps both an input and output typemap must be supplied. With
 * this construction, the output dimension is supplied using the input
 * typemap
 *
 * The typemaps can be applied to existing functions using the
 * %apply directive.  For example:
 *
 *     %apply (double** DYNARGOUT_ARRAY1, size_t* DIM1) {(double** coefs, size_t* length)};
 *     %apply void DYNARGOUT_ARRAY1_NPY_DOUBLE {void getFilter}
 *     void getFilter(double** coefs, size_t* length);
 *
 * The C-function must provide the length and a reference to the data.
 * The data is copied by using DYNARGOUT_ARRAY1 and referenced, when
 * using DYNARGOUTVIEW_ARRAY1.
 *
 * or directly with
 *
 *     void DYNARGOUT_ARRAY1_NPY_DOUBLE(double** DYNARGOUT_ARRAY1, size_t* DIM1);
 *
 */
%define %numpy_ex_typemaps(DATA_TYPE, DATA_TYPECODE, DIM_TYPE)

/*********************************/
/* Dynamic Output Array Typemaps */
/*********************************/


// Apparently this segfaults

/* Typemap suite for (DATA_TYPE** DYNARGOUT_ARRAY1, DIM_TYPE* DIM1)
 */
%typemap(in,numinputs=0,noblock=1)
(DATA_TYPE** DYNARGOUT_ARRAY1, DIM_TYPE* DIM1)
{
  DIM_TYPE templen;
  DATA_TYPE* tempData;
  $1 = &tempData;
  $2 = &templen;
}
%typemap(out)
void DYNARGOUT_ARRAY1_##DATA_TYPECODE
{
  npy_intp dims[1] = { (npy_intp) templen };
  PyObject *array = PyArray_SimpleNew(1, dims, DATA_TYPECODE);
  if (!array) SWIG_fail;
  DATA_TYPE* ff = (DATA_TYPE*) (((PyArrayObject *)array)->data);
  memcpy(ff,tempData,sizeof(DATA_TYPE)*dims[0]);
  $result = SWIG_Python_AppendOutput($result,array);
}

/* Typemap suite for (DATA_TYPE** DYNARGOUTVIEW_ARRAY1, DIM_TYPE* DIM1)
 */
%typemap(in,numinputs=0,noblock=1)
(DATA_TYPE** DYNARGOUTVIEW_ARRAY1, DIM_TYPE* DIM1)
{
  DIM_TYPE templen;
  DATA_TYPE* tempData;
  $1 = &tempData;
  $2 = &templen;
}
%typemap(out)
void DYNARGOUTVIEW_ARRAY1_##DATA_TYPECODE
{
  npy_intp dims[1] = { (npy_intp) templen };
  PyObject *array = PyArray_SimpleNewFromData(1, dims, DATA_TYPECODE, (void*)(tempData));
  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,array);
}


/* Typemap suite for (DATA_TYPE** DYNARGOUTVIEW_ARRAY2_FIXED[ANY], DIM_TYPE* DIM1)
 */
%typemap(in,numinputs=0,noblock=1)
(DATA_TYPE* DYNARGOUTVIEW_ARRAY2_FIXED[ANY], DIM_TYPE* DIM1)
{
  size_t secondDim = $1_dim0;
  size_t firstDim;
  DATA_TYPE* tempData;
  $1 = &tempData;
  $2 = &firstDim;
}
%typemap(out)
void DYNARGOUTVIEW_ARRAY2_FIXED_##DATA_TYPECODE
{
  npy_intp dims[2] = {(npy_intp) firstDim, (npy_intp) secondDim};
%#ifdef __GNUC__
__extension__
%#endif
  PyObject * array = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE, (void*)(tempData));
  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,array);
}


/* Typemap suite for (DATA_TYPE*** DYNARGOUTVIEW_ARRAY2, DIM_TYPE* DIM1)
 */
%typemap(in,numinputs=0,noblock=1)
(DATA_TYPE** DYNARGOUTVIEW_ARRAY2, DIM_TYPE* DIM1, DIM_TYPE* DIM2)
{
  DIM_TYPE dim1,dim2;
  DATA_TYPE* tempData;
  $1 = &tempData;
  $2 = &dim1;
  $3 = &dim2;
}
%typemap(out)
void DYNARGOUTVIEW_ARRAY2_##DATA_TYPECODE
(PyObject * array = NULL)
{
  npy_intp dims[2] = { (npy_intp) dim1, (npy_intp) dim2 };
%#ifdef __GNUC__
__extension__
%#endif
  void* pData = (void*)(tempData);
%#ifdef __GNUC__
__extension__
%#endif
  array = PyArray_SimpleNewFromData(2, dims, DATA_TYPECODE, pData);
  if (!array) SWIG_fail;
  $result = SWIG_Python_AppendOutput($result,array);
}

%enddef    /* %numpy_ex_typemaps() macro */

%numpy_ex_typemaps(double            , NPY_DOUBLE   , size_t)
%numpy_ex_typemaps(float             , NPY_FLOAT    , size_t)
%numpy_ex_typemaps(signed char       , NPY_BYTE     , size_t)
%numpy_ex_typemaps(unsigned char     , NPY_UBYTE    , size_t)
%numpy_ex_typemaps(short             , NPY_SHORT    , size_t)
%numpy_ex_typemaps(unsigned short    , NPY_USHORT   , size_t)
%numpy_ex_typemaps(int               , NPY_INT      , size_t)
%numpy_ex_typemaps(unsigned int      , NPY_UINT     , size_t)
%numpy_ex_typemaps(long              , NPY_LONG     , size_t)
%numpy_ex_typemaps(unsigned long     , NPY_ULONG    , size_t)
%numpy_ex_typemaps(long long         , NPY_LONGLONG , size_t)
%numpy_ex_typemaps(unsigned long long, NPY_ULONGLONG, size_t)
%numpy_ex_typemaps(float             , NPY_FLOAT    , size_t)

#endif
