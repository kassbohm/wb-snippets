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
M,l,EI = var("M,l,EI")

sub_list=[
    ( M,   10 *Newton*m              ),
    ( l,    1 *m                     ),
    ( EI, 200*GPa * 2*mm*6*mm**3/ 12 ),
    ]
l2 = l*l
l3 = l*l*l

K = EI/l3
K *= Matrix(
[
[  4*l2 ,  -6*l ,  2*l2 ,   6*l ],
[ -6*l  ,  12   , -6*l  , -12   ],
[  2*l2 ,  -6*l ,  4*l2 ,   6*l ],
[  6*l  , -12   ,  6*l  ,  12   ],
]
)

p2 = var("psi2")
M1,F1,F2 = var("M1,F1,F2")

u = Matrix([0,0,p2,0])
f = Matrix([M1,F1,M,F2])

unks = [p2,M1,F1,F2]

eq = Eq(K*u , f)
sol = solve(eq, unks)

p2 = sol[p2]
M1 = sol[M1]
F1 = sol[F1]
F2 = sol[F2]

pprint("\nM1 / Nm:")
tmp = M1
pprint(tmp)
tmp = tmp.subs(sub_list)
tmp /= Newton*m
tmp = iso_round(tmp,1)
pprint(tmp)

pprint("\nF1 / N:")
tmp = F1
tmp = tmp.subs(sub_list)
tmp /= Newton
tmp = iso_round(tmp,1)
pprint(tmp)

pprint("\nF2 / N:")
tmp = F2
tmp = tmp.subs(sub_list)
tmp /= Newton
tmp = iso_round(tmp,1)
pprint(tmp)

pprint("\nψ₂:")
tmp = p2
pprint(tmp)

pprint("\nψ₂ / rad:")
tmp = p2
tmp = tmp.subs(sub_list)
tmp = iso_round(tmp,1)
pprint(tmp)

# M1 / Nm:
# M
# ─
# 2
# 5
#
# F1 / N:
# -15
#
# F2 / N:
# 15
#
# ψ₂:
# M⋅l
# ────
# 4⋅EI
#
# ψ₂ / rad:
# 12
