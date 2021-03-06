# https://colab.research.google.com/github/kassbohm/wb-snippets/blob/master/ipynb/TEM_10/ESA/a1_cc.ipynb

F,l = var("F,l")
R = 3*F/2
lu = l/sqrt(3)

Ah,Av,Bh,Bv,Ch,Cv = var("Ah,Av,Bh,Bv,Ch,Cv")

e1 = Eq(Ah + Bh + F)
e2 = Eq(Av + Bv - R)
e3 = Eq(Bv*l - Bh*l - F*l/2 - R*7/18*l)

e4 = Eq(Ch - Bh)
e5 = Eq(Cv - F - Bv)
e6 = Eq(F*lu/2 + Bv*lu + Bh*l)

eqs = [e1,e2,e3,e4,e5,e6]
unknowns = [Ah,Av,Bh,Bv,Ch,Cv]


pprint("\nEquations:")
for e in eqs:
    pprint(e)
    pprint("\n")


# Alternative Solution (also correct):
# Ah,Av,Bh,Bv,Gh,Gv = var("Ah,Av,Bh,Bv,Gh,Gv")
#
# e1 = Eq(Av + Gv - R)
# e2 = Eq(Ah + F - Gh)
# e3 = Eq(F/2 + 7*R/18 - Gv - Gh)
# e4 = Eq(-Gv -F + Bv)
# e5 = Eq(Gh - Bh)
# e6 = Eq(Gh - sqrt(3)*F/6 - Gv/sqrt(3))
#
# eqs = [e1,e2,e3,e4,e5,e6]
# unknowns = [Ah,Av,Bh,Bv,Gh,Gv]

sol = solve(eqs,unknowns)

pprint("\nReactions:")
pprint(sol)

pprint("\nReactions / F (rounded to 0.01):")
for v in  sorted(sol,key=default_sort_key):
    pprint("\n\n")
    s = sol[v]
    tmp = (s/F)
    tmp = tmp.simplify()
    # pprint(tmp)
    pprint([v, tmp, iso_round(tmp,0.01)])

# Reactions / F:
#
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
