from sympy import *
[t,tau,f0,W] = symbols(['t','tau','f0','W'])

g = []
g.append(0.5*sin(2*pi*f0*t))
g.append(-0.5*cos(2*pi*f0*t))
g.append(-0.5*cos(2*pi*t/W)*sin(2*pi*f0*t))
g.append(0.5*cos(2*pi*t/W)*cos(2*pi*f0*t))
g.append(-0.5*sin(2*pi*t/W)*sin(2*pi*f0*t))
g.append(0.5*sin(2*pi*t/W)*cos(2*pi*f0*t))

f = []
f.append(cos(2*pi*f0*tau))
f.append(sin(2*pi*f0*tau))
f.append(cos(2*pi*tau/W)*cos(2*pi*f0*tau))
f.append(cos(2*pi*tau/W)*sin(2*pi*f0*tau))
f.append(sin(2*pi*tau/W)*cos(2*pi*f0*tau))
f.append(sin(2*pi*tau/W)*sin(2*pi*f0*tau))
         
