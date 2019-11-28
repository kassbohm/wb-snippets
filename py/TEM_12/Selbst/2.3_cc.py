from sympy.physics.units import *
from sympy import *

import decimal
from decimal import Decimal as DX
def iso_round(obj, pv, rounding=decimal.ROUND_HALF_EVEN):
    import sympy
    """
    Rounding according to: DIN EN ISO 80000-1:2013-08
    pv: place value = Rundestellenwert
    """
    # pv     :  Rounded digit
    # 1      :  last digit before decimal
    # 0.1    :  1st digit after decimal
    # 0.01   :  2nd digit
    # 0.001  :  3rd digit
    assert pv in set([ # round to:
        100,           #  3rd last digit before decimal
        10,            #  2nd last
        1,             #  last
        0.1,           #  1st digit after decimal
        0.01,          #  2nd
        0.001,         #  3rd
        0.0001,        #  4th
        0.00001,       #  5th
        0.000001,      #  6th
        0.0000001,     #  7th
        0.00000001,    #  8th
        0.000000001,   #  9th
        0.0000000001,  # 10th
        ])
    try:
        tmp = float(obj)
        tmp = DX(str(tmp))
        obj = tmp.quantize(DX(str(pv)), rounding=rounding)
    except:
        for i in range(len(obj)):
            tmp = float(obj[i])
            tmp = DX(str(tmp))
            obj[i] = tmp.quantize(DX(str(pv)), rounding=rounding)
    return obj

# LateX:
kwargs = {"mat_str": "bmatrix",  "mat_delim": ""}

kilo = 1000
mega = 1000*1000
giga = 1000*1000*1000
(mm, cm) = (m/1000, m/100)

deg = pi/180
Newton = kg*m/s**2
Pa = Newton/m**2
MPa = mega*Pa
GPa = giga*Pa
kN = kilo*Newton

###

w, l = var("w, l")

vB = Matrix([0,      w*l,   0])
vC = Matrix([-2*w*l,   0,   0])

vAx, vAy = var("vAx, vAy")
w3, w4 = var("w3, w4")

W3 = Matrix([0, 0, w3])
W4 = Matrix([0, 0, w4])

vA = Matrix([vAx, vAy, 0])

dAC = Matrix([-l,  0, 0])
dAB = Matrix([ 0, -l, 0])

eq1 = Eq(vC, vA + W3.cross(dAC))
eq2 = Eq(vB, vA + W4.cross(dAB))

sol = solve([eq1, eq2],[vAx, vAy, w3, w4])
pprint(sol)

# {vAx: -2⋅l⋅w, vAy: l⋅w, w₃: w, w₄: 2⋅w}
