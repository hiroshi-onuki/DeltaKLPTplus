from quaternion import *

p = 2^248 * 5 - 1
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
assert O0.is_maximal()
N = random_prime(ceil(p**(1/2)))
M = 2**128

beta1 = O0(0)
beta2 = O0(0)
while beta1 * (1 + qi)/2 in O0 or beta2 * (1 + qi)/2 in O0:
    I1, beta1 = RandomFixedNormIdeal(O0, N*M)
    I2, beta2 = RandomFixedNormIdeal(O0, N*M)
assert not beta1/2 in O0
assert not beta2/2 in O0

J1, J2, newN = EquivalentIdealsWithSameNormSmallN(I1, I2, N*M)
print(float(log(newN, 2)), float(1/2*log(p, 2) + 1/2*log(N*M, 2)))

"""
N = random_prime(ceil(p**(1/4)))
I, alpha = RandomFixedNormIdeal(O0, N)
Jd, _ = RandomFixedNormIdeal(O0, M)
J = (I.conjugate() * I.intersection(Jd)) * (1/N)
L, nu = newKLPT(I, J, 2, 1300)
print(factor(norm(L)))
print(norm(L) == 2**1300)
"""