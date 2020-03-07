from sympy import *

F,l = var("F,l")

Ah,Av,Bh,Bv,Ch,Cv = var("Ah,Av,Bh,Bv,Ch,Cv")

# Correct:
R = 3*F/2

lu = l/sqrt(3)

# left:
e1 = Eq(Ah + Bh + F)
e2 = Eq(Av + Bv - R)
e3 = Eq(Bv*l - Bh*l - F*l/2 - R*7/18*l)

e4 = Eq(Ch - Bh)
e5 = Eq(Cv - F - Bv)
e6 = Eq(F*lu/2 + Bv*lu + Bh*l)

eqs = [e1,e2,e3,e4,e5,e6]
unknowns = [Ah,Av,Bh,Bv,Ch,Cv]
sol = solve(eqs,unknowns)

pprint("Reactions / F:")
for v in  sorted(sol,key=default_sort_key):
    pprint("\n\n")
    s = sol[v]
    sF = (s/F).simplify()
    pprint([v, sF, N(sF,2)])

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
