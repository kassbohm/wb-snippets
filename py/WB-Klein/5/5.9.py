EI = 1 *Newton*m**2
l = 1 *m
l2 = l*l

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

pprint("\nSolution:")
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
