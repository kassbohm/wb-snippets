F, E, A, l = var("F, E, A, l")

# Quantities:
sub_list = [
    ( F,        5 *kN    ),
    ( E, 200*1000 *MPa   ),
    ( A,       25 *mm**2 ),
    ( l,     1707 *mm    ),
    ]

c = sqrt(2)/2
EA = E*A
S1, S2, dl1, dl2, u = var("S1, S2, dl1, dl2, u")

eq1 = Eq( F/2 + S2 + S1*c )
eq2 = Eq( S1, EA*dl1 / (sqrt(2)*l) )
eq3 = Eq( S2, EA/2*dl2/l )
eq4 = Eq(dl1, -u*c)
eq5 = Eq(dl2, -u)

eqns = [eq1, eq2, eq3, eq4, eq5]
unks = [S1, S2, dl1, dl2, u]
sol = solve(eqns, unks)
pprint(sol)

u = sol[u]
u = u.subs(sub_list)

pprint("\nu / mm:")
tmp = u
tmp /= mm
tmp = iso_round(tmp,0.001)
pprint(tmp)

# ⎧      2⋅F   √2⋅F        4⋅F   √2⋅F         -√2⋅F⋅l            -2⋅F⋅l            2⋅F⋅l    ⎫
# ⎨S₁: - ─── + ────, S₂: - ─── + ────, dl₁: ────────────, dl₂: ────────────, u: ────────────⎬
# ⎩       7     14          7     7         A⋅E⋅(√2 + 4)       A⋅E⋅(√2 + 4)     A⋅E⋅(√2 + 4)⎭
#
# u / mm:
# 0.631
