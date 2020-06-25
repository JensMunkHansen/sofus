#!/usr/bin/env python
# -*- coding: utf-8; tab-width: 2; python-indent: 2; indent-tabs-mode: nil -*-

import sys
import os
import re
import numpy as np

from sys import version_info

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
    self.nSetMethods = len(tuple(filter(lambda i: not re.compile('^_').search(i), setmethods)))
  def test_properties_get(self):

    if version_info < (3, 0, 0):
      # Filter out builtin_function_type
      getmethods = {k: v for k,v in fnm.ApertureFloat.__swig_getmethods__.iteritems() if type(v) != types.BuiltinFunctionType}
      # Filter out lambda functions (constructors)
      getmethods = {k: v for k,v in getmethods.iteritems() if v.func_name != '<lambda>'}.keys()

      # Update without pulsed waves also
      if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
        self.assertEqual(len(getmethods), 34)
      else:
        self.assertEqual(len(getmethods), 23)
    else:
      getmethods = {k: v for k,v in fnm.ApertureFloat.__swig_getmethods__.items()}

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

    nSetMethods = 0
    if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
      nSetMethods=28
    else:
      nSetMethods=17

    #if 'FNM_CLOSURE_FUNCTIONS' in fnm.__dict__.keys():
    #  nSetMethods = nSetMethods + 3

    self.assertEqual(len(setmethods),nSetMethods)

    argss = [[0,0.1,0.1,0.1],[1,0.1,0.1,0.1],[100,0.1,0.1,0.1]]

    nSetSuccess = 0
    testMe = set()
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

  def test_closures_dimensions(self):
    if tuple(fnm.ApertureFloat.__dict__.keys()).count('RwFloatParamSetMulti') > 0:

      a = fnm.ApertureFloat(5,1.0,0.0,1.0)

      if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
        params = [fnm.RwParamType.Excitation,
                  fnm.RwParamType.Impulse]

        # Variable dimensions
        for param in params:
          refValue = np.random.rand(10).astype(np.float32)
          a.RwFloatParamSetMulti(param, refValue, 1, 10)
          testValue = a.RwFloatParamGet(param,1)[1]
          self.assertTrue(np.all(refValue == testValue))

      params = [fnm.RwParamType.ElementDelays,
                fnm.RwParamType.Apodization,
                fnm.RwParamType.Focus,
                fnm.RwParamType.CenterFocus]

      # Fixed dimensions
      for param in params:
        refValue = a.RwFloatParamGet(param,1)[1]
        refValue = refValue + np.random.rand(len(refValue)).astype(np.float32)
        a.RwFloatParamSetMulti(param, refValue, 1, len(refValue))
        testValue = a.RwFloatParamGet(param,1)[1]
        self.assertTrue(np.all(refValue == testValue))

  def test_closures_scalar_float(self):
    # Not compiling using MSVC 2013
    if tuple(fnm.ApertureFloat.__dict__.keys()).count('RwFloatParamSetMulti') > 0:
      # TODO: Expose list of properties from C/C++
      props = [[fnm.RwParamType.Alpha, 'alpha'],
               [fnm.RwParamType.Beta,  'beta'],
               [fnm.RwParamType.C,     'c'],
               [fnm.RwParamType.F0,    'f0'],
               [fnm.RwParamType.Fs,    'fs'],
               [fnm.RwParamType.W,     'w']]

      a = fnm.ApertureFloat(1,1.0,0.0,1.0)
      refValue = 0.0
      for prop in props:
        refValue = refValue + 1.0
        a.RwFloatParamSetMulti(prop[0], refValue)
        value = eval('a.'+prop[1])
        self.assertEqual(value,refValue)
        value1 = a.RwFloatParamGet(prop[0],0)[1]
        self.assertEqual(value1,refValue)

  def test_setting_elements(self):
    a = fnm.ApertureFloat()
    elements0 = np.random.rand(32,8).astype(np.float32)
    a.elements = elements0
    elements1 = a.elements
    self.assertTrue(np.all(elements0==elements1))

    if tuple(fnm.ApertureFloat.__dict__.keys()).count('RwFloatParamSetMulti') > 0:
      elements2 = np.random.rand(32,8).astype(np.float32)
      a.RwFloatParamSetMulti(fnm.RwParamType.Elements, elements2, 2, 32, 8)
      elements3 = a.elements
      self.assertTrue(np.all(elements2==elements3))

  def test_setting_subelements(self):
    a = fnm.ApertureFloat()
    subelements0 = np.random.rand(32,4,8).astype(np.float32)
    a.subelements = subelements0
    subelements1 = a.subelements
    self.assertTrue(np.all(subelements0==subelements1))

if __name__ == '__main__':
  unittest.main()

# Local variables: #
# tab-width: 2 #
# python-indent: 2 #
# indent-tabs-mode: nil #
# End: #
