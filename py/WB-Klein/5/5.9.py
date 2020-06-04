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

pprint("\nStiffness matrix / N:")
pprint(K/Newton)

w2, p2m, p2p = var("w2, p2m, p2p")
u = Matrix([w2, p2m, p2p])
f = Matrix([1*Newton,0,0])

eq = Eq(K*u , f)
sol = solve(eq, [w2, p2m, p2p])

pprint("\nSolution 1:")
pprint(sol)

# Stiffness matrix / N:
# ⎡18          ⎤
# ⎢──   6   -3 ⎥
# ⎢m           ⎥
# ⎢            ⎥
# ⎢6   4⋅m   0 ⎥
# ⎢            ⎥
# ⎣-3   0   2⋅m⎦
#
# Solution:
# ⎧                         2⋅m⎫
# ⎨p2m: -1/3, p2p: 1/3, w₂: ───⎬
# ⎩                          9 ⎭


pprint("\nSolution 2:")
A, F = var("A, F")
EI, l = var("EI, l")
sub_list = [
    (F,  1*Newton),
    (EI, 1*Newton*m**2),
    (l, 1*m),
    ]

l2 = l*l
l3 = l*l*l

pl, w, pr, A = var("phi_l, w, phi_r, A")

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

pprint("\nφL, φR / rad:")
for s in [pl, pr]:
    s = s.subs(sub_list)
    pprint(s)

pprint("\nw / m:")
w = w.subs(sub_list)
w /= m
pprint(w)
