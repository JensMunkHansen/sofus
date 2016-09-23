#!/usr/bin/env python
import sys
import re
import numpy as np

import unittest
import subprocess

import addpaths

import swig_fnm as fnm

print('Hello')
class FnmTest(unittest.TestCase):
    def test_properties_get(self):
        getmethods = fnm.ApertureFloat.__swig_getmethods__.keys()

        assert(len(getmethods)==13)

        argss = [[0,0.1,0.1,0.1],[1,0.1,0.1,0.1],[100,0.1,0.1,0.1]]
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
            for args in argss:
                a = fnm.ApertureDouble(*args)
                for method in getmethods:
                    try:
                        value = eval('a.'+method)
                        nGetSuccess = nGetSuccess + 1
                    except:
                        raise(Exception('error getting '+method))
        except Exception as e:
            print(e.message)
      
        self.assertEqual(nGetSuccess/6,len(getmethods))
        
if __name__ == '__main__':
    unittest.main()

