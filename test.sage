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
print(f"Generated ideals I1, I2 with norm {N*M} and generators beta1, beta2")

L, nu = deltaKLPT(I1, I2, 2, 1500)
print(factor(norm(L)))

n1 = random_prime(ceil(p**(0.6)))
n2 = random_prime(ceil(p**(2.7)))
J = KLPT(L, n1, n2)
print(factor(norm(J)))


