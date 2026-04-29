from quaternion import *

p = 103
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.maximal_order()
N = 10007

I1, beta1 = RandomFixedNormIdeal(O0, N)
I2, beta2 = RandomFixedNormIdeal(O0, N)

C, D = IdealModConstraint(O0, qj, qk, beta1, beta2, N)

