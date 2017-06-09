#!/usr/bin/env python

import sys
import re
import numpy as np

import unittest
import subprocess

import matplotlib.pyplot as plt

import addpaths

import swig_fnm as fnm

from euler import euler2rot

plt.ion()

class FnmFourTest(unittest.TestCase):
  def test_edge_sym(self):
    """
    Scatter outside four edges
    """
    width  = 1.0
    kerf   = 0.0
    height = 1.0
    a = fnm.ApertureFloat(1, width, kerf, height)
    a.c          = 1.0
    a.f0         = 1.0
    a.nSubH      = 20
    a.nSubW      = 20

    pos = np.r_[1.5,0,5]
    pos = pos.reshape([1,3]).astype(np.float32)

    hp = a.CalcCwFieldFourRef(pos)[1]

    hpr = np.ones(4)*hp
    hps = []

    for iEdge in range(4):
      rotm = euler2rot(iEdge*np.pi/2.0,0,0,conv='zxz',intrinsic=True)
      pos1 = np.dot(rotm,pos.T).astype(np.float32).T
      hp = a.CalcCwFieldFourRef(pos1)[1]
      hps.append(hp)
    self.assertTrue(np.allclose(hps,hpr))
  def test_superposition(self):
    """
    Test superposition with all signs
    """
    hw = 2.0
    hh = 1.0
    kerf = 0.0
    z = 20.0
    dx = 0.5
    dy = 0.5
    
    a = fnm.ApertureFloat(1, 2*hw, kerf, 2*hh)
    a.f0 = 5
    a.c  = 1
    a.nDivH = 4
    a.nDivW = 4
    
    method = a.CalcCwFieldFourRef
    #method = a.CalcCwFieldRef
    #method = a.CalcCwFast
    
    pos0 = np.r_[hw+dx,hh+dy,z]
    pos0 = pos0.reshape([1,3]).astype(np.float32)

    hpr = method(pos0)[1]

    pos0 = pos0.flatten()

    hps = []

    # Large rectangle
    w = abs(pos0[0])+hw
    h = abs(pos0[1])+hh

    a.elements = np.r_[np.c_[abs(w)/2.0,abs(h)/2.0,0,0,0,0,0,0]].astype(np.float32)
    pos = np.r_[abs(w)/2.0,abs(h)/2.0,z]
    pos = pos.reshape([1,3]).astype(np.float32)
    hp = method(pos)[1]
    hps.append(np.sign(w)*np.sign(h)*hp)
    

    # Consider adding sign
    w = hw - abs(pos0[0])
    h = abs(pos0[1])+hh
    a.elements = np.r_[np.c_[abs(w)/2.0,h/2.0,0,0,0,0,0,0]].astype(np.float32)
    pos = np.r_[abs(w)/2.0,abs(h)/2.0,z]
    pos = pos.reshape([1,3]).astype(np.float32)
    hp = method(pos)[1]
    hps.append(np.sign(w)*hp)
    

    # Consider adding sign
    w = abs(pos0[0])+hw
    h = hh - abs(pos0[1])
    a.elements = np.r_[np.c_[w/2.0,abs(h)/2.0,0,0,0,0,0,0]].astype(np.float32)
    pos = np.r_[abs(w)/2.0,abs(h)/2.0,z]
    pos = pos.reshape([1,3]).astype(np.float32)
    hp = method(pos)[1]
    hps.append(np.sign(h)*hp)
    
    
    # Consider adding sign
    w = hw - abs(pos0[0])
    h = hh - abs(pos0[1])
    a.elements = np.r_[np.c_[abs(w)/2.0,abs(h)/2.0,0,0,0,0,0,0]].astype(np.float32)
    pos = np.r_[abs(w)/2.0,abs(h)/2.0,z]
    pos = pos.reshape([1,3]).astype(np.float32)
    hp = method(pos)[1]
    hps.append(np.sign(w)*np.sign(h)*hp)
    

    hp = np.sum(hps)
    diff = np.abs(hp - hpr)
    print(diff)
    print(abs(hpr))
if __name__ == '__main__':
  unittest.main()
    
# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
    
