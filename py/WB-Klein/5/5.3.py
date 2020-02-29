def k(phi):
    """ element stiffness matrix """
    # phi: Angle between two rays 1 and 2:
    # 1: ray along global x axis and
    # 2: ray along 1-2-axis of rod
    # phi is counted positively about z.
    (c, s) = ( cos(phi), sin(phi) )
    (cc, ss, sc) = ( c*c, s*s, s*c)
    return Matrix(
        [
        [ cc,  sc, -cc, -sc],
        [ sc,  ss, -sc, -ss],
        [-cc, -sc,  cc,  sc],
        [-sc, -ss,  sc,  ss],
        ])

F, c = var("F, c")

# Stiffness Matrix:
k1 = c*k(135 *pi/180)
pprint("\n\nk1 / c: ")
pprint(k1/c)

# Linear System:
u, F1x, F2x, F2y = var("u, F1x, F2x, F2y")
f_ = Matrix([F1x, -F/2, F2x, F2y])
u_ = Matrix([0, -u, 0, 0])

sol = solve(Eq(k1*u_,f_), [u, F1x, F2x, F2y], dict=True)
pprint("\n\nSolution:")
pprint(sol[0])

# k1 / c:
# ⎡1/2   -1/2  -1/2  1/2 ⎤
# ⎢                      ⎥
# ⎢-1/2  1/2   1/2   -1/2⎥
# ⎢                      ⎥
# ⎢-1/2  1/2   1/2   -1/2⎥
# ⎢                      ⎥
# ⎣1/2   -1/2  -1/2  1/2 ⎦
#
#
# Solution:
# ⎧     F       -F        F     F⎫
# ⎨F1x: ─, F2x: ───, F2y: ─, u: ─⎬
# ⎩     2        2        2     c⎭
