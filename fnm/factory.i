/* %factory() macro
 *
 * Macro for wrapping factory constructors
 */
%define %factory_typemaps(CLASS)

%typemap(in, numinputs=0)
    (CLASS** OUTPUT_INST)
    (CLASS* temp)
{
  $1 = &temp;
}

%typemap(argout)
    (CLASS ** OUTPUT_INST)
{
  PyObject* temp = NULL;
  if (!PyList_Check($result)) {
    temp = $result;
    $result = PyList_New(1);
    PyList_SetItem($result, 0, temp);
  }

  temp = SWIG_NewPointerObj(SWIG_as_voidptr(*$1),
                            $descriptor(CLASS*),
                            SWIG_POINTER_OWN | 0);

  PyList_Append($result, temp);
  Py_DECREF(temp);
}

%enddef    /* %factory_typemaps() */

/*
 *     %apply (bftx::Aperture<DATA_TYPE>** OUTPUT_INST) {(bftx::Aperture<DATA_TYPE>** obj)};
*/
