from quaternion import *

p = 2^17-1
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.maximal_order()
N = random_prime(p**4)

I1, beta1 = RandomFixedNormIdeal(O0, N)
I2, beta2 = RandomFixedNormIdeal(O0, N)
IJ = I1.intersection(I2)
alpha = SmallGenerator(I1.conjugate() * I1)
J = IJ.right_order().right_ideal([alpha.inverse() * b for b in (I1.conjugate() * IJ).basis()])
print(I1.right_order() == J.left_order())
print(norm(I1), norm(I2), norm(J))

L, nu = newKLPT(I1, J, 2, 4*ceil(log(p, 2)))
