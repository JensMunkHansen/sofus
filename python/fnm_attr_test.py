#!/usr/bin/env python

import sys
import re
import numpy as np

import unittest
import subprocess

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
    getmethods = fnm.ApertureFloat.__swig_getmethods__.keys()

    if 'FNM_PULSED_WAVE' in fnm.__dict__.keys():
      assert(len(getmethods)==25)
    else:
      assert(len(getmethods)==22)
      
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
      assert(len(setmethods)==19)
    else:
      assert(len(setmethods)==16)

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
        
if __name__ == '__main__':
  unittest.main()

# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
    
