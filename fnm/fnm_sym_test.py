#!/usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt

import addpaths
from dicts import dotdict
import swig_fnm as fnm
plt.ion()

from fnm import rect

def compareWithPython(**kwargs):
  opt = dotdict({'show'     : False,
                 'method'   : 'CalcCwFieldRef',
                 'location' : 'inside'})
  opt.update(**kwargs)
  
  areas = [2.0,3.0,4.0,5.0]
  widths  = np.array([areas[0],1.0,areas[2],1.0],dtype=np.float32)
  heights = np.array([1.0,areas[1],1.0,areas[3]],dtype=np.float32)

  xsign = np.array([1.0,-1.0,-1.0,1.0],dtype=np.float32)
  ysign = np.array([1.0,1.0,-1.0,-1.0],dtype=np.float32)

  show = opt.show
  # Ensure that we either slightly inside or outside
  scale = 40.0 * np.finfo(np.float32).eps

  epss = dict({'inside'  : -scale,
               'outside' :  scale})
  
  ndiv = 2
  iCorners = [0,1,2,3]

  eps    = epss[opt.location]
  method = opt.method

  scatters = np.c_[(widths+eps)/2.0 * xsign,
                   (heights+eps)/2.0 * ysign,
                   np.ones(4,dtype=np.float32)]

  print('Testing: %s vs %s: %s' % (opt.method,'H_accurate (Python)',opt.location))
  for iCorner in iCorners:
    a = fnm.ApertureFloat(1,float(widths[iCorner]),
                          float(0.0),
                          float(heights[iCorner]))
    a.nthreads = 1
    a.nDivH = ndiv
    a.nDivW = ndiv
    a.f0 = 1
    a.c  = 1
    if show:
      a.show()
      ax = plt.gca()
      ax.scatter([scatters[iCorner][0]],
                 [scatters[iCorner][1]],
                 [scatters[iCorner][2]],color='black',marker='o',s=20)
      
      ax.set_xlim([-widths.max(),widths.max()])
      ax.set_ylim([-heights.max(),heights.max()])
      ax.set_zlim([-5,5])
      ax.set_aspect('equal')
      plt.show()
    
    exec('result = a.'+method+'([scatters[iCorner]])')

    r = rect(hw=widths[iCorner]/2.0,hh=heights[iCorner]/2.0,center=[0,0,0],nAbcissa=ndiv)
          
    k = 2*np.pi / (a.c / a.f0)
  
    xs = np.ones((1,1))*scatters[iCorner][0]
    ys = np.ones((1,1))*scatters[iCorner][1]
    zs = np.ones((1,1))*scatters[iCorner][2]
  
    reference = r.H_accurate(xs,ys,zs,k)
    assert(np.abs(reference) - np.abs(result[1]) < np.finfo(np.float32).eps)
    diff = np.abs(reference)-np.abs(result[1])
    print('Relative difference: %f' % (diff / np.mean(np.abs(reference)+np.abs(result[1]))))


def compareWithReference(**kwargs):
  opt = dotdict({'method0'   : 'CalcCwFieldRef',
                 'method1'   : 'CalcCwFieldRef',
                 'location' : 'inside',
                 'ndiv'     : 2})

  opt.update(**kwargs)
  
  areas = [2.0,3.0,4.0,5.0]
  widths  = np.array([areas[0],1.0,areas[2],1.0],dtype=np.float32)
  heights = np.array([1.0,areas[1],1.0,areas[3]],dtype=np.float32)

  xsign = np.array([1.0,-1.0,-1.0,1.0],dtype=np.float32)
  ysign = np.array([1.0,1.0,-1.0,-1.0],dtype=np.float32)

  show = opt.show
  # Ensure that we either slightly inside or outside
  scale = 40.0 * np.finfo(np.float32).eps

  epss = dict({'inside'  : -scale,
               'outside' :  scale})
  
  ndiv = opt.ndiv
  iCorners = [0,1,2,3]

  eps    = epss[opt.location]
  method0 = opt.method0
  method1 = opt.method1

  scatters = np.c_[(widths+eps)/2.0 * xsign,
                   (heights+eps)/2.0 * ysign,
                   np.ones(4,dtype=np.float32)]
  
  print('Testing: %s vs %s: %s' % (opt.method0,opt.method1, opt.location))
  for iCorner in iCorners:
    a = fnm.ApertureFloat(1,
                          float(widths[iCorner]),
                          float(0.0),
                          float(heights[iCorner]))
    a.nthreads = 1
    a.nDivH = ndiv
    a.nDivW = ndiv
    a.f0 = 1
    a.c  = 1
    
    exec('result0 = a.'+method0+'([scatters[iCorner]])[1]')
    exec('result1 = a.'+method1+'([scatters[iCorner]])[1]')
    diff = np.abs(result0)-np.abs(result1)
    assert(diff < 100*np.finfo(np.float32).eps)
    print('Relative difference: %f' % (diff / np.mean(np.abs(result0)+np.abs(result1))))

if __name__ == "__main__":
  compareWithPython(location='inside', method='CalcCwFieldRef')
  compareWithPython(location='outside',method='CalcCwFieldRef')

  # SSE with C refernece
  compareWithReference(location='outside',method1='CalcCwFast')
  compareWithReference(location='inside', method1='CalcCwFast')
  # When inside the errors should be identical (they are)
  compareWithReference(location='inside', method0='CalcCwField',method1='CalcCwField2')
