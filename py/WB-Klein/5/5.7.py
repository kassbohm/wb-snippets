M,l,EI = var("M,l,EI")

sub_list=[
    ( M,   10 *Newton*m              ),
    ( l,    1 *m                     ),
    ( EI, 200*GPa * 2*mm*6*mm**2/ 12 ),
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

# M1 / Nm:
# 5
#
# F1 / N:
# -15
#
# F2 / N:
# 15
