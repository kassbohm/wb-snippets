llEA, l, F1, F2 = var("EA, l, F1, F2")

sub_list = [
    ( EA,  2 *Pa*m**2 ),
    ( l,   1 *m       ),
    ( F1,  1 *Newton  ),
    ( F2,  2 *Newton  ),
    ]

def k(phi):
    """ element stiffness matrix """
    # phi is angle between:
    # 1. vector along global x axis
    # 2. vector along 1-2-axis of truss
    # phi is counted positively about z.
    # pprint("phi / deg:")
    # pprint(N(deg(phi),3))
    (c, s) = ( cos(phi), sin(phi) )
    (cc, ss, sc) = ( c*c, s*s, s*c)
    return Matrix(
        [
        [ cc,  sc, -cc, -sc],
        [ sc,  ss, -sc, -ss],
        [-cc, -sc,  cc,  sc],
        [-sc, -ss,  sc,  ss],
        ])


(p1, p2, p3) = (135 *pi/180, 0 *pi/180, 45 *pi/180)
(k1, k2, k3) = (EA/l*k(p1), EA/l*k(p2), EA/l*k(p3))

pprint("\nk1 / (EA / l): ")
pprint(k1 / (EA/l) )
pprint("\nk2 / (EA / l): ")
pprint(k2 / (EA/l) )
pprint("\nk3 / (EA / l): ")
pprint(k3 / (EA/l) )

A = EA/l*Matrix([
    [ S(3)/2 ,  -1 ],
    [ -1,   S(3)/2 ]
    ])

u2x, u3x = var("u2x, u3x")
u = Matrix([u2x  , u3x  ])
f = Matrix([F1/2 , F2/2 ])

u2x, u3x = var("u2x, u3x")

eq = Eq(A*u , f)
sol = solve(eq, [u2x, u3x])
pprint("\nSolution:")
pprint(sol)

u2x, u3x = sol[u2x], sol[u3x]

pprint("\nu2x:")
tmp = u2x.subs(sub_list)
pprint(N(tmp,4))

pprint("\nu3x:")
tmp = u3x.subs(sub_list)
pprint(N(tmp,4))

pprint("\nF1x:")
tmp = - EA/l * u2x/2
tmp = tmp.subs(sub_list)
tmp /= Newton
pprint(tmp)


# k1 / (EA / l):
# ⎡1/2   -1/2  -1/2  1/2 ⎤
# ⎢                      ⎥
# ⎢-1/2  1/2   1/2   -1/2⎥
# ⎢                      ⎥
# ⎢-1/2  1/2   1/2   -1/2⎥
# ⎢                      ⎥
# ⎣1/2   -1/2  -1/2  1/2 ⎦
#
# k2 / (EA / l):
# ⎡1   0  -1  0⎤
# ⎢            ⎥
# ⎢0   0  0   0⎥
# ⎢            ⎥
# ⎢-1  0  1   0⎥
# ⎢            ⎥
# ⎣0   0  0   0⎦
#
# k3 / (EA / l):
# ⎡1/2   1/2   -1/2  -1/2⎤
# ⎢                      ⎥
# ⎢1/2   1/2   -1/2  -1/2⎥
# ⎢                      ⎥
# ⎢-1/2  -1/2  1/2   1/2 ⎥
# ⎢                      ⎥
# ⎣-1/2  -1/2  1/2   1/2 ⎦
#
# Solution:
# ⎧     l⋅(3⋅F₁ + 2⋅F₂)       l⋅(2⋅F₁ + 3⋅F₂)⎫
# ⎨u2x: ───────────────, u3x: ───────────────⎬
# ⎩           5⋅EA                  5⋅EA     ⎭
#
# u2x:
# 0.7⋅m
#
# u3x:
# 0.8⋅m
#
# F1x:
# -7/10
