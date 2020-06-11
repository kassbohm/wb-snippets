# https://colab.research.google.com/github/kassbohm/wb-snippets/blob/master/ipynb/TEM_10/ESA/a1.ipynb

from sympy import *

F,l = var("F,l")
R = 3*F/2
lu = l/sqrt(3)

# Kai
Ah,Av,Bh,Bv,Ch,Cv = var("Ah,Av,Bh,Bv,Ch,Cv")

e1 = Eq(Ah + Bh + F)
e2 = Eq(Av + Bv - R)
e3 = Eq(Bv*l - Bh*l - F*l/2 - R*7/18*l)

e4 = Eq(Ch - Bh)
e5 = Eq(Cv - F - Bv)
e6 = Eq(F*lu/2 + Bv*lu + Bh*l)

eqs = [e1,e2,e3,e4,e5,e6]
unknowns = [Ah,Av,Bh,Bv,Ch,Cv]

# # Armin:
# Ax,Ay,Bx,By,Gx,Gy = var("Ax,Ay,Bx,By,Gx,Gy")
#
# e1 = Eq(Ay + Gy - R)
# e2 = Eq(Ax + F - Gx)
# e3 = Eq(F/2 + 7*R/18 - Gy - Gx)
# e4 = Eq(-Gy -F + By)
# e5 = Eq(Gx - Bx)
# e6 = Eq(Gx - sqrt(3)*F/6 - Gy/sqrt(3))
#
# eqs = [e1,e2,e3,e4,e5,e6]
# unknowns = [Ax,Ay,Bx,By,Gx,Gy]


sol = solve(eqs,unknowns)

pprint("Reactions / F:")
for v in  sorted(sol,key=default_sort_key):
    pprint("\n\n")
    s = sol[v]
    sF = (s/F).simplify()
    pprint([v, sF, N(sF,2)])

# Kai:
# ⎡      43   19⋅√3       ⎤
# ⎢Ah, - ── + ─────, -0.42⎥
# ⎣      24     24        ⎦
#
#
# ⎡      3   19⋅√3     ⎤
# ⎢Av, - ─ + ─────, 1.0⎥
# ⎣      8     24      ⎦
#
#
# ⎡      19⋅√3   19       ⎤
# ⎢Bh, - ───── + ──, -0.58⎥
# ⎣        24    24       ⎦
#
#
# ⎡      19⋅√3   15     ⎤
# ⎢Bv, - ───── + ──, 0.5⎥
# ⎣        24    8      ⎦
#
#
# ⎡      19⋅√3   19       ⎤
# ⎢Ch, - ───── + ──, -0.58⎥
# ⎣        24    24       ⎦
#
#
# ⎡      19⋅√3   23     ⎤
# ⎢Cv, - ───── + ──, 1.5⎥
# ⎣        24    8      ⎦


# Armin:
#
# ⎡      43   19⋅√3       ⎤
# ⎢Ax, - ── + ─────, -0.42⎥
# ⎣      24     24        ⎦
#
#
# ⎡      3   19⋅√3     ⎤
# ⎢Ay, - ─ + ─────, 1.0⎥
# ⎣      8     24      ⎦
#
#
# ⎡      19   19⋅√3      ⎤
# ⎢Bx, - ── + ─────, 0.58⎥
# ⎣      24     24       ⎦
#
#
# ⎡      19⋅√3   23     ⎤
# ⎢By, - ───── + ──, 1.5⎥
# ⎣        24    8      ⎦
#
#
# ⎡      19   19⋅√3      ⎤
# ⎢Gx, - ── + ─────, 0.58⎥
# ⎣      24     24       ⎦
#
#
# ⎡      19⋅√3   15     ⎤
# ⎢Gy, - ───── + ──, 0.5⎥
# ⎣        24    8      ⎦
