import numpy as np

def f(x):
    return x**5 + x**2

# Accurate upto 2*n - 1

order = 3
xs = np.polynomial.legendre.leggauss(order)[0]

result = np.sum(f(xs) * np.polynomial.legendre.leggauss(order)[1])

# Test integrating x**2 from [-2,6]

xs = np.polynomial.legendre.leggauss(order)[0]
xs = xs * (6 - -2)/2.0 + (6 + -2)/2.0
result = np.sum(f(xs) * (6 - -2)/2.0 * np.polynomial.legendre.leggauss(order)[1])

a = -2
b = 2
order = 3

xs = np.polynomial.legendre.leggauss(order)[0]
xs = xs * (b - a)/2.0 + (b + a)/2.0
result1 = np.sum(f(xs) * (b - a)/2.0 * np.polynomial.legendre.leggauss(order)[1])

a = 2
b = 6

xs = np.polynomial.legendre.leggauss(order)[0]
xs = xs * (b - a)/2.0 + (b + a)/2.0
result2 = np.sum(f(xs) * (b - a)/2.0 * np.polynomial.legendre.leggauss(order)[1])

