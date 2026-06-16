from id2iso import Qlapoti
from quaternion import RandomFixedNormIdeal

e = 248
f = 5
p = 2^e * f - 1  
Bpinf.<qi, qj, qk> = QuaternionAlgebra(-1, -p)
O0 = Bpinf.maximal_order(order_basis=(Bpinf(1), qi, (qi + qj)/2, (1 + qk)/2))
N = random_prime(10*p)
I, _ = RandomFixedNormIdeal(O0, N)
print(f"Generated ideal I with norm {N}")

beta1, beta2 = Qlapoti(I, e-2)
print(f"Found elements in I with norms {beta1.reduced_norm()}, {beta2.reduced_norm()}")
