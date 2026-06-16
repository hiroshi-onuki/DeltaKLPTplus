from sage.all import (
    norm,
    vector,
    gcd,
    ZZ,
    log,
)

from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from quaternion import SmallestEquivalentIdeal, SmallGenerator, SumOf2Squares
from lattice import ClosestVectorDim2Euclidean

def Qlapoti(I, N, max_tries=1000):
    _, qi, _, _ = I.quaternion_algebra().basis()
    O = I.left_order()
    assert qi**2 == -1

    I = SmallestEquivalentIdeal(I)
    n = norm(I)

    while max_tries > 0:
        max_tries -= 1
        alpha = SmallGenerator(I)
        aa, ba, _, _ = alpha
        r = ZZ(alpha.reduced_norm()/n)
        if gcd(2*aa, n) > 1 and gcd(2*ba, n) > 1:
            continue
        if gcd(2*aa, n) == 1:
            ainv = ZZ(2*aa).inverse_mod(n)
            b1 = vector([-2*ba * ainv, 1])
            b2 = vector([n, 0])
            target = vector([(N - 2*r) * ainv % n, 0])
        else:
            binv = ZZ(2*ba).inverse_mod(n)
            b1 = vector([1, -2*aa * binv])
            b2 = vector([0, n])
            target = vector([0, (N - 2*r) * binv % n])
        v = ClosestVectorDim2Euclidean(b1, b2, target)
        v = target - v
        s, t = v
        print(float(log(s**2 + t**2, 2)), float(log(n, 2)))
        assert (2*aa * s + 2*ba * t - (N - 2*r)) % n == 0
        z = 2 * (N - 2*r - 2*aa*s - 2*ba*t) / n - s**2 - t**2
        if z < 0:
            continue
        if z % 4 == 0 and not (s % 2 == 0 and t % 2 == 0):
            continue
        if z % 4 == 1 and s % 2 == t % 2:
            continue
        if z % 4 == 2 and not(s % 2 == t % 2 == 1):
            continue
        if z % 4 == 3:
            continue
        try:
            print(f"Found suitable s, t, z = {s}, {t}, {z}")
            z0, z1 = SumOf2Squares(z)
        except:
            continue
        if not z0 % 2 == s % 2:
            z0, z1 = z1, z0
        a1 = (z0 + s) / 2
        b1 = (z1 + t) / 2
        a2 = s - a1
        b2 = t - b1
        beta1 = n*(a1 + b1 * ii) + alpha
        beta2 = n*(a2 + b2 * ii) + alpha
        I1 = O0.left_ideal([b * beta1.conjugate() / n for b in I.basis()])
        I2 = O0.left_ideal([b * beta2.conjugate() / n for b in I.basis()])
        assert norm(I1) + norm(I2) == N
        return I1, I2
    raise ValueError("No suitable ideals found")
