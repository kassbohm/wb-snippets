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
F, E, A, l = var("F, E, A, l")

# Quantities:
sub_list = [
    ( F,        5 *kN    ),
    ( E, 200*1000 *MPa   ),
    ( A,       25 *mm**2 ),
    ( l,     1707 *mm    ),
    ]

c = sqrt(2)/2
EA = E*A
S1, S2, dl1, dl2, u = var("S1, S2, dl1, dl2, u")

eq1 = Eq( F/2 + S2 + S1*c )
eq2 = Eq( S1, EA*dl1 / (sqrt(2)*l) )
eq3 = Eq( S2, EA/2*dl2/l )
eq4 = Eq(dl1, -u*c)
eq5 = Eq(dl2, -u)

eqns = [eq1, eq2, eq3, eq4, eq5]
unks = [S1, S2, dl1, dl2, u]
sol = solve(eqns, unks)
pprint(sol)

u = sol[u]
u = u.subs(sub_list)

pprint("\nu / mm:")
tmp = u
tmp /= mm
tmp = iso_round(tmp,0.001)
pprint(tmp)

# ⎧         √2⋅F           √2⋅F         -√2⋅F⋅l            -2⋅F⋅l            2⋅F⋅l    ⎫
# ⎨S₁: -F + ────, S₂: -F + ────, dl₁: ────────────, dl₂: ────────────, u: ────────────⎬
# ⎩          2              2         A⋅E⋅(√2 + 2)       A⋅E⋅(√2 + 2)     A⋅E⋅(√2 + 2)⎭
#
# u / mm:
# 1.000
