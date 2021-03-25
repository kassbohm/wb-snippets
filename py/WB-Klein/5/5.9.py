# https://colab.research.google.com/github/kassbohm/wb-snippets/blob/master/ipynb/WB-Klein/5/5.9_cc.ipynb

pprint("\nSolution 1:")
A, F = var("A, F")
EI, l = var("EI, l")
sub_list = [
    (F,  1*Newton),
    (EI, 1*Newton*m**2),
    (l, 1*m),
    ]

l2 = l*l
l3 = l*l*l

pl, w, pr, A = var("psi_l, w, psi_r, A")

eq1 = Eq( EI/l3   * (  4*l2 *pl +  6*l * w) , 0 )
eq2 = Eq( EI/l3   * (  6*l  *pl + 12   * w) , + A + F/2 )
eq3 = Eq( EI/2/l3 * (  4*l2 *pr -  6*l * w) , 0 )
eq4 = Eq( EI/2/l3 * ( -6*l  *pr + 12   * w) , - A + F/2 )
eqns = [eq1, eq2, eq3, eq4]
unks = [pl, w, pr, A]
sol = solve(eqns, unks)
pprint(sol)
pl = sol[pl]
w = sol[w]
pr = sol[pr]
A = sol[A]

pprint("\nψₗ, ψᵣ / rad:")
for s in [pl, pr]:
    s = s.subs(sub_list)
    pprint(s)

pprint("\nw / m:")
w = w.subs(sub_list)
w /= m
pprint(w)

pprint("\nSolution 2:")

EI = 1 *Newton*m**2
l = 1 *m
l2 = l*l
l3 = l*l*l

K = EI/l**3
K *= Matrix([
[  18 ,  6*l, -3*l  ],
[  6*l, 4*l2,  0    ],
[ -3*l,   0,  2*l2  ],
])

# pprint("\nStiffness matrix / N:")
# pprint(K/Newton)

w2, p2m, p2p = var("w₂, ψ₂⁻, ψ₂⁺")
u = Matrix([w2, p2m, p2p])
f = Matrix([1*Newton,0,0])

eq = Eq(K*u , f)
sol = solve(eq, [w2, p2m, p2p])

pprint(sol)

# Solution 1:
# ⎧              2          2          3⎫
# ⎪   F      -F⋅l        F⋅l      2⋅F⋅l ⎪
# ⎨A: ─, φₗ: ──────, φᵣ: ────, w: ──────⎬
# ⎪   6       3⋅EI       3⋅EI      9⋅EI ⎪
# ⎩                                     ⎭
#
# φL, φR / rad:
# -1/3
# 1/3
#
# w / m:
# 2/9
#
# Solution 2:
#
# Stiffness matrix / N:
# ⎡  18                   ⎤
# ⎢─────     6       -3   ⎥
# ⎢meter                  ⎥
# ⎢                       ⎥
# ⎢  6    4⋅meter     0   ⎥
# ⎢                       ⎥
# ⎣ -3       0     2⋅meter⎦
# ⎧                         2⋅meter⎫
# ⎨p2m: -1/3, p2p: 1/3, w₂: ───────⎬
# ⎩                            9   ⎭
