# Header starts here.
from sympy.physics.units import *
from sympy import *

# Rounding:
import decimal
from decimal import Decimal as DX
from copy import deepcopy
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
    objc = deepcopy(obj)
    try:
        tmp = DX(str(float(objc)))
        objc = tmp.quantize(DX(str(pv)), rounding=rounding)
    except:
        for i in range(len(objc)):
            tmp = DX(str(float(objc[i])))
            objc[i] = tmp.quantize(DX(str(pv)), rounding=rounding)
    return objc

# LateX:
kwargs = {}
kwargs["mat_str"] = "bmatrix"
kwargs["mat_delim"] = ""
# kwargs["symbol_names"] = {FB: "F^{\mathsf B}", }

# Units:
(k, M, G ) = ( 10**3, 10**6, 10**9 )
(mm, cm) = ( m/1000, m/100 )
Newton = kg*m/s**2
Pa     = Newton/m**2
MPa    = M*Pa
GPa    = G*Pa
kN     = k*Newton
deg    = pi/180

half = S(1)/2

# Header ends here.
#
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
