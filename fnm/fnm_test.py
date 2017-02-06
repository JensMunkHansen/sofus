#!/usr/bin/env python

# Work in progress
import addpaths

import swig_fnm as fnm
import sys
import re

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

getmethods = fnm.ApertureFloat.__swig_getmethods__.keys()
setmethods = fnm.ApertureFloat.__swig_setmethods__.keys()

assert(len(getmethods)==22)
assert(len(setmethods)==15)

argss = [[0,0.1,0.1,0.1],[1,0.1,0.1,0.1],[100,0.1,0.1,0.1]]
nGetSuccess = 0
nSetSuccess = 0

doubleSupported = True

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
        for method in setmethods:
            if not(method.find('_') == 0):
                try:
                    value = eval('a.'+method)
                    exec('a.'+method+'=value')
                    nSetSuccess = nSetSuccess + 1
                except:
                    print('error setting '+method)
                    raise(Exception('error setting '+method))
    
    for args in argss:
        a = fnm.ApertureDouble(*args)
        for method in getmethods:
            try:
                value = eval('a.'+method)
                nGetSuccess = nGetSuccess + 1
            except:
                raise(Exception('error getting '+method))
        for method in setmethods:
            if not(method.find('_') == 0):
                try:
                    value = eval('a.'+method)
                    exec('a.'+method+'=value')
                    nSetSuccess = nSetSuccess + 1
                except:
                    print('error setting '+method)
                    raise(Exception('error setting '+method))
except Exception as e:
    doubleSupported = False
    print(e.message)

factor = 3
if doubleSupported:
    factor = 6
  
print('Getting %d out of %d properties' % (nGetSuccess/factor, len(getmethods)))
print('Setting %d out of %d properties' % (nSetSuccess/factor, len(setmethods)))

assert(nGetSuccess/factor == len(getmethods))
# Not possible to set properties starting with _
assert(nSetSuccess/factor == len(filter(lambda i: not re.compile('^_').search(i), setmethods)))

a = fnm.ApertureFloat(1,1.0,0.0,1.0)

pos = np.meshgrid(np.r_[-1:2],np.r_[-1:2])
pos = np.c_[pos[0].flatten(),pos[1].flatten(),np.ones(9)].astype(np.float32)

a = fnm.ApertureFloat(1,1.0,0.0,1.0)
for method in [a.CalcCwFast, a.CalcCwField, a.CalcCwField2, a.CalcCwFieldRef]:
    result = method(pos)[1].reshape((3,3))
    result = result.flatten()
    r0 = result[1:8:2]
    r1 = result[[0,2,6,8]]
    # Check symmetry
    #assert(np.allclose(r0,r0[0]))
    #assert(np.allclose(r1,r1[0]))
    print(r0[0])
    print(r1[0])
    print(result[4])

pos = np.meshgrid(np.r_[-0.25:0.5:0.25],np.r_[-0.25:0.5:0.25])
pos = np.c_[pos[0].flatten(),pos[1].flatten(),np.ones(9)].astype(np.float32)

# TODO: Use python reference and use thin elements to track error

print('Inside')
reference = a.CalcCwFast(pos)[1].reshape((3,3)).flatten()
for method in [a.CalcCwField, a.CalcCwField2, a.CalcCwFieldRef]:
    result = method(pos)[1].reshape((3,3)).flatten()
    r0 = result[1:8:2]
    r1 = result[[0,2,6,8]]
    # Check symmetry
    #assert(np.allclose(r0,r0[0]))
    #assert(np.allclose(r1,r1[0]))
    if not(np.all(reference == result)):
        print('method ' + str(method.im_func) + ' not working properly inside')
    print(r0[0])
    print(r1[0])
    print(result[4])

# TODO: Gives NaN over corner.
eps = np.finfo(np.float32).eps
a = fnm.ApertureFloat(1,1.0-eps,0.0,1.0-eps)
cornerPositions = np.array([-0.5,0.5,4.0,0.5,0.5,4.0,0.5,-0.5,4.0,-0.5,-0.5,4.0],dtype=np.float32).reshape((4,3))
for iCorner in range(4):
    result = a.CalcCwFast([cornerPositions[iCorner]])
    print(np.abs(result[1]))

print('Just outside')
# Gives NaN over corner, if dz=0
eps = 4*eps
a = fnm.ApertureFloat(1,1.0-eps,0.0,1.0-eps)
# Just outside each corner (almost equal, when eps is larger, they are equal)
cornerPositions = np.array([-0.5,0.5,0.0,0.5,0.5,0.0,0.5,-0.5,0.0,-0.5,-0.5,0.0],dtype=np.float32).reshape((4,3))
for iCorner in range(4):
    a = fnm.ApertureFloat()
    center = cornerPositions[iCorner % 4].copy()
    a.elements = [np.array([(1.0-eps)/2.0,(1.0-eps)/2.0,center[0],center[1],center[2],0,0,0],dtype=np.float32)]
    scatter = cornerPositions[(iCorner+2) % 4].copy()
    scatter[2] = 4.0
    exec('result = a.CalcCwFast([scatter])')
#    exec('result = a.CalcCwFieldRef([scatter])') # Should be equal (it is not)
#    exec('result = a.CalcCwField([scatter])')    # Should be different
#    exec('result = a.CalcCwField2([scatter])')   # Should equal CalcCwField (it is)
    print(np.abs(result[1]))
    
# Local variables: #
# indent-tab-mode: nil #
# tab-width: 2 #
# python-indent: 2 #
# py-indent-offset: 2 #
# indent-tabs-mode: nil #
# End: #
