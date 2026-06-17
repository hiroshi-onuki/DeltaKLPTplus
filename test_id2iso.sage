from id2iso import IdealToIsogeny
from quaternion import RandomFixedNormIdeal
from special_curve import SpecialSuperSingularCurve

e = 248
f = 5
p = 2^e * f - 1  

E0withEnd = SpecialSuperSingularCurve(p, e, f)

for _ in range(10):
    N = random_prime(10*p)
    I, _ = RandomFixedNormIdeal(E0withEnd.order, N)
    print(f"Generated ideal I with norm {N}")

    IdealToIsogeny(E0withEnd, I)
    print("Isogeny computed successfully")