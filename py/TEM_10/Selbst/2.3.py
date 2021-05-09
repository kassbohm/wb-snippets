from sympy.physics.units import *
from sympy import *

sqrt2 = sqrt(2)

# Given symbols:
q0, l = var("q0, l")

x = var("x")

q = q0 * x / (sqrt2 * l)

# R:
R = integrate(q, (x, 0, sqrt2*l))
pprint("\nR:")
pprint(R)

I = integrate(q * x, (x, 0, sqrt2*l))
xR = I / R

pprint(" \nxʀ:")
pprint(xR)
