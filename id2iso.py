from lattice import EnumerateCloseVectorsDim2Euclidean, ShortBasisDim2Euclidean
from sage.all import (
    norm,
    vector,
    gcd,
    ZZ,
    ceil,
    inverse_mod,
)
from theta_structures.couple_point import CouplePoint
from theta_isogenies.product_isogeny import EllipticProductIsogeny
from quaternion import SmallestEquivalentIdeal, SmallGenerator, SumOf2Squares
from utilities.discrete_log import tate_pairing_pari

def Qlapoti(I, e, max_tries=10000):
    _, qi, _, _ = I.quaternion_algebra().basis()
    assert qi**2 == -1
    O = I.left_order()
    N = 2**e

    I, beta0 = SmallestEquivalentIdeal(I)
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
        rb1, rb2 = ShortBasisDim2Euclidean(b1, b2)
        vs = EnumerateCloseVectorsDim2Euclidean(rb1, rb2, target, 10000, ceil(2 * (N - 2*r) / n))
        for v in vs:
            v = target - v
            s, t = v
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
            z0, z1 = SumOf2Squares(z)
            if z0 is None or z1 is None:
                continue
            if not z0 % 2 == s % 2:
                z0, z1 = z1, z0
            a1 = (z0 + s) / 2
            b1 = (z1 + t) / 2
            a2 = s - a1
            b2 = t - b1
            beta1 = n*(a1 + b1 * qi) + alpha
            beta2 = n*(a2 + b2 * qi) + alpha
            assert beta1.reduced_norm() + beta2.reduced_norm() == N*n

            if beta1.reduced_norm() / n % 2 == 0:
                continue
            gamma = beta2*beta1.conjugate() / n
            v_gamma = vector(gamma) * O.basis_matrix().inverse()
            if (v_gamma - vector([1, 0, 0, 0])) % 2 == 0 or (v_gamma - vector([0, 1, 0, 0])) % 2 == 0:
                continue
            beta1 = beta1 * beta0 / n
            beta2 = beta2 * beta0 / n
            return beta1, beta2, gamma
    raise ValueError("No suitable ideals found")

def IdealToIsogeny(E0withEnd, I):
    E0 = E0withEnd.curve
    O = E0withEnd.order
    N = norm(I)
    assert I.left_order() == O
    assert N % 2 == 1
    e = E0withEnd.e

    """
    compute beta1, beta2, gamma s.t.
    beta1 = hat{phi_I} * phi_1, beta2 = hat{phi_I} * phi_2,
    gamma = hat{phi_2} * phi_1,
    where phi_I is the isogeny corresponding to I,
    deg(phi_1) + deg(phi_2) = 2^(e-2).
    """
    beta1, beta2, gamma = Qlapoti(I, e - 2)
    assert beta1 in I and beta2 in I
    assert beta2 * beta1.conjugate() / N == gamma
    d1 = beta1.reduced_norm() / N
    d2 = beta2.reduced_norm() / N
    assert d1 + d2 == 2**(e - 2)
    assert d1 % 2 == 1 and d2 % 2 == 1

    # compute (2^(e-2), 2^(e-2))-isogeny
    P1, Q1 = d1 * E0withEnd.P, d1 * E0withEnd.Q
    P2, Q2 = E0withEnd.quaternion_action(gamma)
    K1 = CouplePoint(P1, P2)
    K2 = CouplePoint(Q1, Q2)
    Phi = EllipticProductIsogeny((K1, K2), e-2)
    P0, Q0 = E0withEnd.quaternion_action(beta1.conjugate()) # P0 = hat{phi_1}*phi_I(P), Q0 = hat{phi_1}*phi_I(Q)
    image1 = Phi(CouplePoint(P0, E0(0)))
    image2 = Phi(CouplePoint(Q0, E0(0)))
    image_sum = Phi(CouplePoint(P0 + Q0, E0(0)))

    # determine the codomain EI by comparing the Tate pairings
    tP0Q0d1 = E0withEnd.tate_pairing_PQ**(d1**2*N)
    for idx in range(2):
        Pim, Qim, PQim = image1[idx], image2[idx], image_sum[idx]
        if not (Pim + Qim == PQim or Pim + Qim == -PQim):
            Qim = -Qim
        assert Pim + Qim == PQim or Pim + Qim == -PQim
        exp = (E0withEnd.p**2 - 1) / 2**e
        tPimQim = tate_pairing_pari(Pim, Qim, 2**e)**exp
        if tPimQim == tP0Q0d1:
            EI = Phi.codomain()[idx]
            PI = inverse_mod(d1, 2**e) * Pim
            QI = inverse_mod(d1, 2**e) * Qim
            return EI, PI, QI
    raise ValueError("Failed to determine the codomain of the isogeny")
