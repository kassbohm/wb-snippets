from sympy import *

a3, x, l = var("a3, x, l")

w = a3*x*x*(x - l)

tmp = integrate(w, (x,0,l))
pprint(tmp)

#      4 
# -a₃⋅l
# ───────
#    12
