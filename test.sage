from quaternion import *

p = 2^248 * 5 - 1
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
assert O0.is_maximal()
N = random_prime(p)
M = random_prime(p)

I1, beta1 = RandomFixedNormIdeal(O0, N*M)
I2, beta2 = RandomFixedNormIdeal(O0, N*M)

#beta1, beta2 = EquivalentIdealsWithSameNorm(I1, I2, N, M)
#print(float(log(beta1.reduced_norm()/(N*M), 2)), float(3/4*log(p, 2) + 1/4*log(N*M, 2)))

N = random_prime(ceil(p**(1/4)))
I, alpha = RandomFixedNormIdeal(O0, N)
Jd, _ = RandomFixedNormIdeal(O0, M)
J = (I.conjugate() * I.intersection(Jd)) * (1/N)
L, nu = newKLPT(I, J, 2, 1300)
print(factor(norm(L)))
print(norm(L) == 2**1300)