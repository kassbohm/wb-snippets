EA, l, F1, F2 = var("EA, l, F1, F2")

sub_list = [
    ( EA,  2 *Pa*m**2 ),
    ( l,   1 *m       ),
    ( F1,  1 *Newton /2  ),
    ( F2,  2 *Newton /2  ),
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


(p1, p2, p3) = (0 *pi/180, 30 *pi/180, 45 *pi/180)
(k1, k2, k3) = (EA/sqrt(3)/l*k(p1), EA/2/l*k(p2), EA/sqrt(2)/l*k(p3))

pprint("\nk1 / (EA / l): ")
pprint(k1 / (EA/l) )
pprint("\nk2 / (EA / l): ")
pprint(k2 / (EA/l) )
pprint("\nk3 / (EA / l): ")
pprint(k3 / (EA/l) )

# A = EA/l*Matrix([
#     [  1     ,  -S(1)/2 ],
#     [ -S(1)/2,       1  ]
#     ])
#
# u2x, u3x = var("u2x, u3x")
# u = Matrix([u2x  , u3x  ])
# f = Matrix([F1 , F2 ])
#
# u2x, u3x = var("u2x, u3x")
#
# eq = Eq(A*u , f)
# sol = solve(eq, [u2x, u3x])
# pprint("\nSolution:")
# pprint(sol)
#
# u2x, u3x = sol[u2x], sol[u3x]
#
# pprint("\nu2x / m:")
# tmp = u2x.subs(sub_list)
# tmp /= m
# pprint(tmp)
#
# pprint("\nu3x / m:")
# tmp = u3x.subs(sub_list)
# tmp /= m
# pprint(tmp)
#
# pprint("\nF1x / N:")
# tmp = - EA/l * u2x/2
# tmp = tmp.subs(sub_list)
# tmp /= Newton
# pprint(tmp)
