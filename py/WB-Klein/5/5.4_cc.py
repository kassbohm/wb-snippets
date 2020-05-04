from sympy.physics.units import *
from sympy import *

# Rounding:
import decimal
from decimal import Decimal as DX
def iso_round(obj, pv, rounding=decimal.ROUND_HALF_EVEN):
    import sympy
    """
    Rounding acc. to DIN EN ISO 80000-1:2013-08
    place value = Rundestellenwert
    """
    assert pv in set([
        # place value   #  round to:
        1,              #  1
        0.1,            #  1st digit after decimal
        0.01,           #  2nd
        0.001,          #  3rd
        0.0001,         #  4th
        0.00001,        #  5th
        0.000001,       #  6th
        0.0000001,      #  7th
        0.00000001,     #  8th
        0.000000001,    #  9th
        0.0000000001,   # 10th
        ])
    try:
        tmp = DX(str(float(obj)))
        obj = tmp.quantize(DX(str(pv)), rounding=rounding)
    except:
        for i in range(len(obj)):
            tmp = DX(str(float(obj[i])))
            obj[i] = tmp.quantize(DX(str(pv)), rounding=rounding)
    return obj

# LateX:
kwargs = {}
kwargs["mat_str"] = "bmatrix"
kwargs["mat_delim"] = ""
# kwargs["symbol_names"] = {FB: "F^{\mathsf B}", }

# Units:
(k, M, G ) = ( 10**3, 10**6, 10**9 )
(mm, cm, deg) = ( m/1000, m/100, pi/180)
Newton = kg*m/s**2
Pa     = Newton/m**2
MPa    = M*Pa
GPa    = G*Pa
kN     = k*Newton

half = S(1)/2

# ---

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
