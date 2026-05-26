from quaternion import *

proof.all(False)

p = 2^248 * 5 - 1
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
assert O0.is_maximal()
N1 = random_prime(ceil(p**4))
N2 = random_prime(ceil(p**4))

beta1 = O0(0)
beta2 = O0(0)
while beta1 * (1 + qi)/2 in O0 or beta2 * (1 + qi)/2 in O0:
    I1, beta1 = RandomFixedNormIdeal(O0, N1)
    I2, beta2 = RandomFixedNormIdeal(O0, N2)
assert not beta1/2 in O0
assert not beta2/2 in O0
print(f"Generated ideals I1, I2 with norm {N1}, {N2} and generators beta1, beta2")

n1 = random_prime(ceil(p**(0.7)))
n2 = random_prime(ceil(p**(2.8)))
J1 = KLPT(I1, n1, n2)
print(norm(J1) == n1*n2)
J2 = KLPT(I2, n1, n2)
print(norm(J2) == n1*n2)

N = n1*n2
cnt = 0
while norm(J1) > 2**255:
    J1, J2, N = EquivalentIdealsWithSameNorm(J1, J2, N)
    cnt += 1
print(norm(J1) == norm(J2) == N)
print(float(log(N, 2)), cnt)

L = deltaKLPT(J1, J2, 2, 260*4)


