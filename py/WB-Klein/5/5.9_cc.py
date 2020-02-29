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
