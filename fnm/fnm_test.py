#!/usr/bin/env python
import addpaths

import swig_fnm as fnm
import sys
import re

import numpy as np
import matplotlib.pyplot as plt
plt.ion()

getmethods = fnm.ApertureFloat.__swig_getmethods__.keys()
setmethods = fnm.ApertureFloat.__swig_setmethods__.keys()

assert(len(getmethods)==18)
assert(len(setmethods)==12)

argss = [[0,0.1,0.1,0.1],[1,0.1,0.1,0.1],[100,0.1,0.1,0.1]]
nGetSuccess = 0
nSetSuccess = 0

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
    print(e.message)

print('Getting %d out of %d properties' % (nGetSuccess/6, len(getmethods)))
print('Setting %d out of %d properties' % (nSetSuccess/6, len(setmethods)))

assert(nGetSuccess/6 == len(getmethods))
# Not possible to set properties starting with _
assert(nSetSuccess/6 == len(filter(lambda i: not re.compile('^_').search(i), setmethods)))

            
sys.exit(0)
