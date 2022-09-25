EA, l, F1, F2 = var("EA, l, F1, F2")
Newton = kg*m/s**2
Pa     = Newton/m**2

sub_list = [
    ( EA,  2 *Pa*m**2 ),
    ( l,   1 *m       ),
    ( F1,  1 *Newton /2  ), # due to symmetry
    ( F2,  2 *Newton /2  ), # due to symmetry
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

(p1, p2, p3) = (315*pi/180, 0 *pi/180, 45 *pi/180)
# k2 uses only 1/2 A due to symmetry:
(k1, k2, k3) = (EA/l*k(p1), EA/2/l*k(p2), EA/l*k(p3))

pprint("\nk1 / (EA / l): ")
pprint(k1 / (EA/l) )
pprint("\nk2 / (EA / l): ")
pprint(k2 / (EA/l) )
pprint("\nk3 / (EA / l): ")
pprint(k3 / (EA/l) )

K = EA/l*Matrix([
    [  1     ,  -S(1)/2 ],
    [ -S(1)/2,       1  ]
    ])


u2x, u3x = var("u2x, u3x")
u = Matrix([u2x  , u3x  ])
f = Matrix([F1 , F2 ])

u2x, u3x = var("u2x, u3x")

eq = Eq(K*u , f)
sol = solve(eq, [u2x, u3x])
pprint("\nSolution:")
pprint(sol)

u2x, u3x = sol[u2x], sol[u3x]

pprint("\nu2x / m:")
tmp = u2x.subs(sub_list)
tmp /= m
pprint(tmp)

pprint("\nu3x / m:")
tmp = u3x.subs(sub_list)
tmp /= m
pprint(tmp)

pprint("\nF1x / N:")
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
# ⎡1/2   0  -1/2  0⎤
# ⎢                ⎥
# ⎢ 0    0   0    0⎥
# ⎢                ⎥
# ⎢-1/2  0  1/2   0⎥
# ⎢                ⎥
# ⎣ 0    0   0    0⎦
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
# ⎧     2⋅l⋅(2⋅F₁ + F₂)       2⋅l⋅(F₁ + 2⋅F₂)⎫
# ⎨u2x: ───────────────, u3x: ───────────────⎬
# ⎩           3⋅EA                  3⋅EA     ⎭
#
# u2x / m:
# 2/3
#
# u3x / m:
# 5/6
#
# F1x / N:
# -2/3
