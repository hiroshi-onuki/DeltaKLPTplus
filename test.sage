from quaternion import *

p = 103
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1,qi, (qi+qj)/2,(1+qk)/2])
N = 1000000007

I, beta1 = RandomFixedNormIdeal(O0, N)
J, beta2 = RandomFixedNormIdeal(O0, N)

C, D = IdealModConstraint(O0, qj, qk, beta1, beta2, N)
gamma = C * qj + D * qk

print(beta1 * gamma in J)

