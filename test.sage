from quaternion import *

p = 103
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.maximal_order()
N = 2003

I1, beta1 = RandomFixedNormIdeal(O0, N)
I2, beta2 = RandomFixedNormIdeal(O0, N)

C, D = IdealModConstraint(O0, qj, qk, beta1, beta2, N)
nu = StrongApproximation(O0, N, C, D, 2)
print("nu = %s, Nrd(nu) = %s" % (nu, factor(nu.reduced_norm())))
