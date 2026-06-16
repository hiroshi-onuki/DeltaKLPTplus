from id2iso import Qlapoti
from quaternion import RandomFixedNormIdeal

e = 248
f = 5
p = 2^e * f - 1  
D = 2^(e-2)
Bpinf.<ii, jj, kk> = QuaternionAlgebra(-1, -p)
O0 = Bpinf.maximal_order(order_basis=(Bpinf(1), ii, (ii + jj)/2, (1 + kk)/2))
N = random_prime(10*p)
I, _ = RandomFixedNormIdeal(O0, N)

I1, I2 = Qlapoti(I, D)
print(f"Found ideals I1, I2 with norm {norm(I1)}, {norm(I2)}")
