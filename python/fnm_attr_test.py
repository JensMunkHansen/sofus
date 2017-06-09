#!/usr/bin/env python

import sys
import re
import numpy as np

import unittest
import subprocess

import types

import addpaths

import swig_fnm as fnm

class FnmTest(unittest.TestCase):
  def setUp(self):
    """
    For some reason, the unittest framework returns an error
    
    SyntaxError: unqualified exec is not allowed in function
    'test_properties_set' because it contains a nested function
    with free variables
    
    Workaround is to introduce the member nSetMethods
    """
    setmethods = fnm.ApertureFloat.__swig_setmethods__.keys()
    self.nSetMethods = len(filter(lambda i: not re.compile('^_').search(i), setmethods))
  def test_properties_get(self):

    # Filter out builtin_function_type
    getmethods = {k: v for k,v in fnm.ApertureFloat.__swig_getmethods__.iteritems() if type(v) != types.BuiltinFunctionType}
    # Filter out lambda functions (constructors)
    getmethods = {k: v for k,v in getmethods.iteritems() if v.func_name != '<lambda>'}.keys()

    # Update without pulsed waves also
    if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
      assert(len(getmethods)==28)
    else:
      assert(len(getmethods)==23)
      
    argss = [[0,0.1,0.1,0.1],[5,0.1,0.1,0.1],[100,0.1,0.1,0.1]]
    nGetSuccess = 0
    testMe = set()
    try:
      for args in argss:
        a = fnm.ApertureFloat(*args)
        for method in getmethods:
          try:
            value = eval('a.'+method)
            nGetSuccess = nGetSuccess + 1
          except:
            raise(Exception('error getting '+method))
    except Exception as e:
      print(e.message)
  
    self.assertEqual(nGetSuccess/3,len(getmethods))
  def test_properties_set(self):
    setmethods = fnm.ApertureFloat.__swig_setmethods__.keys()

    if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
      assert(len(setmethods)==22)
    else:
      assert(len(setmethods)==17)

    argss = [[0,0.1,0.1,0.1],[1,0.1,0.1,0.1],[100,0.1,0.1,0.1]]

    nSetSuccess = 0
    testMe = set()
    method = setmethods[0]
    try:
      for args in argss:
        a = fnm.ApertureFloat(*args)
        for method in setmethods:
          if not(method.find('_') == 0):
            try:
              value = eval('a.'+method)
              exec('a.'+method+'=value')
              nSetSuccess = nSetSuccess + 1
            except:
              print('error setting '+method)
              raise(Exception('error setting '+method))
    except:
      print('error setting '+method)
      raise(Exception('error setting '+method))
    self.assertEqual(nSetSuccess/3,self.nSetMethods)

  def test_closures_scalar(self):
    # TODO: Expose list of properties (C)
    props = [[fnm.RwParamType.Alpha, 'alpha'],
             [fnm.RwParamType.Beta,  'beta'],
             [fnm.RwParamType.C,     'c'],
             [fnm.RwParamType.F0,    'f0'],
             [fnm.RwParamType.Fs,    'fs'],
             [fnm.RwParamType.W,     'w']]

    a = fnm.ApertureFloat(1,1.0,0.0,1.0)
    for prop in props:
      a.RwFloatParamSet0D(prop[0],7.0)
      value = eval('a.'+prop[1])
      self.assertEqual(value,7.0)
    
  def test_setting_elements(self):
    a = fnm.ApertureFloat()
    elements0 = np.random.rand(32,8).astype(np.float32)
    a.elements = elements0
    elements1 = a.elements
    self.assertTrue(np.all(elements0==elements1))

  def test_setting_subelements(self):
    a = fnm.ApertureFloat()
    subelements0 = np.random.rand(32,4,8).astype(np.float32)
    a.subelements = subelements0
    subelements1 = a.subelements
    self.assertTrue(np.all(subelements0==subelements1))
        
if __name__ == '__main__':
  unittest.main()

# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
