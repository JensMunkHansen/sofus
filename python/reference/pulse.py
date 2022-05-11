# -*- coding: utf-8 -*-
"""Python reference implementation for computation of arc length of sphere
intersecting a rectange

Example
-------
TODO: How to run this::

    $ python example_numpy.py

Notes
-----
    No notes for this file

"""
import os, sys
import numpy as np

filedir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(filedir, '../../fnm'))

from dicts import dotdict

class pulse(dotdict):
  def __init__(self,*args, **kwargs):
    opt = dotdict({'fs' : 10e6,
                   'f0' : 1e6,
                   'c'  : 1500,
                   'nCycles' : 3})
    opt.update(*args,**kwargs)
    super(pulse,self).update(**opt)
    self.w = self.nCycles / self.f0
  def eval(self,t):
    """
    """
    return (np.fabs(t/self.w - 0.5) < 0.5) * np.sin(2.0*np.pi*self.f0*t)
  def sbf(self, i, tau):
    if i==0:
      return np.cos(2*np.pi*self.f0*tau)
    else:
      return np.sin(2*np.pi*self.f0*tau)
  def tbf(self, i, t):
    if i==0:
      return np.sin(2*np.pi*self.f0*t)
    else:
      return -np.cos(2*np.pi*self.f0*t)

  def adjust(self, sigma, weight, adjacent, z, E):
    eps = np.finfo(np.float32).eps
    # Only works for rectangular response
    denom = adjacent**2 + sigma**2
    if (denom > 0.0):
      factor = weight / denom
    else:
      factor = 0.0
    sqrtArg = z**2 + denom
    tau = np.sqrt(max(0,sqrtArg)) / self.c
    for i in range(2):
      E[i] = E[i] +  factor * self.sbf(i,tau)
