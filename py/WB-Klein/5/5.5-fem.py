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

F, EA, l = var("F, EA, l")

sub_list = [
    (F,               5 *kN ),
    (EA, 200 *GPa*25 *mm**2 ),
    (l,            1707 *mm ),
    ]

(l1, l2) = (l*sqrt(2), l)
(p1, p2) = (45 *pi/180, 90*pi/180)
(k1, k2) = (EA/l1*k(p1), half*EA/l2*k(p2))

pprint("\nk1 / (EA / l): ")
pprint(k1 / (EA/l) )

pprint("\nk2 / (EA / l): ")
pprint(k2 / (EA/l) )

psi = sqrt(2)/4
c = EA/l
u3y = -(F/2) / (half+psi) / c

pprint("\nu3y:")
pprint(u3y)
pprint("\nu3y / mm:")
tmp = u3y
tmp = tmp.subs(sub_list)
tmp /= mm
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\nF1x / N:")
F1x = - c * psi * u3y
tmp = F1x
tmp = tmp.subs(sub_list)
tmp /= Newton
tmp = iso_round(tmp,0.1)
pprint(tmp)

# k1 / (EA / l):
# ⎡ √2    √2   -√2   -√2 ⎤
# ⎢ ──    ──   ────  ────⎥
# ⎢ 4     4     4     4  ⎥
# ⎢                      ⎥
# ⎢ √2    √2   -√2   -√2 ⎥
# ⎢ ──    ──   ────  ────⎥
# ⎢ 4     4     4     4  ⎥
# ⎢                      ⎥
# ⎢-√2   -√2    √2    √2 ⎥
# ⎢────  ────   ──    ── ⎥
# ⎢ 4     4     4     4  ⎥
# ⎢                      ⎥
# ⎢-√2   -√2    √2    √2 ⎥
# ⎢────  ────   ──    ── ⎥
# ⎣ 4     4     4     4  ⎦
#
# k2 / (EA / l):
# ⎡0   0    0   0  ⎤
# ⎢                ⎥
# ⎢0  1/2   0  -1/2⎥
# ⎢                ⎥
# ⎢0   0    0   0  ⎥
# ⎢                ⎥
# ⎣0  -1/2  0  1/2 ⎦
#
# u3y:
#     -F⋅l
# ─────────────
#      ⎛√2   1⎞
# 2⋅EA⋅⎜── + ─⎟
#      ⎝4    2⎠
#
# u3y / mm:
# -1.000
#
# F1x / N:
# 1035.5
