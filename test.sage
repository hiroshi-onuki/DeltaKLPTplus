from quaternion import *

p = 103
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.maximal_order()
N = 10

print(RandomFixedNormIdeal(O0, N))


