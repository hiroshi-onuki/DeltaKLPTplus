from quaternion import *

proof.all(False)

p = 2^248 * 5 - 1
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
assert O0.is_maximal()
Nsk = random_prime(ceil(p**4))
Ncom = random_prime(ceil(p**4))
Nchl = 2^128

beta1 = O0(0)
beta2 = O0(0)
Isk, beta1 = RandomFixedNormIdeal(O0, Nsk)
Icom, beta2 = RandomFixedNormIdeal(O0, Ncom)
Ichl, _ = RandomFixedNormIdeal(O0, Nchl)
print(f"Generated ideals Isk, Icom, Ichl.")

IskIchl = Isk.intersection(Ichl)
print(f"Computed Isk intersect Ichl.")
print(f"Norm of Isk intersect Ichl is {factor(IskIchl.norm())}")
assert IskIchl.left_order() == O0

L, nu = deltaKLPTforSign(Icom, IskIchl, 2, 265*4, 2^255)
beta = SmallGenerator(L)
I = L + O0 * 2**(265*4)
print(I.is_principal())


