from __future__ import division
from random import random
from math import floor
import numpy as np

def fractionize(R, n, d):
    error = abs(n/d - R)
    return (n, d, error)  # (numerator, denominator, absolute difference to R)

def better(a, b):
    return a if a[2] < b[2] else b

def approximate(R, n, m):
    best = (0, 1, R)
    for d in xrange(1, m+1):
        n1 = min(n, int(floor(R * d)))
        n2 = min(n, n1 + 1) # ceil(R*d)
        best = better(best, fractionize(R, n1, d))
        best = better(best, fractionize(R, n2, d))
    return best

def calc_incr_int(startDepth, stopDepth, startValue=0, stopValue=-6):
    fs = 12e6
    ds = 1540.0 / (2*fs)
    a = 10**(startValue / 20.0)
    b = 10**(stopValue / 20.0)
    u = int(startDepth / ds)
    v = int(stopDepth / ds)
    maxStep = 255
    increment = 1.0e-6
    maxIncr   = 2**16 - 1
    R = (np.round((a-b) / increment) -1) / (v-u)
    return (approximate(R, maxIncr, maxStep), (int(R),1,R-int(R)))



if __name__ == '__main__': 
    def main():
        R = random()
        n = 30
        m = 100
        print R, approximate(R, n, m)
    main()
