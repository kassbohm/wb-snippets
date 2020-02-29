c = var("c")

K = Matrix([
    [ 2*c , -c ],
    [  -c ,  c ]
    ])

M, p2, p3 = var("M, p2, p3")
u = Matrix([p2, p3 ])
f = Matrix([0, M ])

eq = Eq(K*u , f)
sol = solve(eq, [p2, p3])
pprint("\nSolution:")
pprint(sol)

# Solution:
# ⎧    M      2⋅M⎫
# ⎨p₂: ─, p₃: ───⎬
# ⎩    c       c ⎭
