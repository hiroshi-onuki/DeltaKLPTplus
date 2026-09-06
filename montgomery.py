"""
x-only arithmetic on Montgomery curves and 2^n-isogeny chains over F_{p^2} = F_p(i), p = 3 mod 4.

Field elements are elements of Sage's GF(p^2) (modulus x^2 + 1). A point on the Kummer line of
E_A: y^2 = x^3 + A x^2 + x is a pair (X, Z) of field elements, Z == 0 meaning the point at
infinity; a curve is carried either as its affine coefficient A or projectively as
(A24p : C24) = (A + 2 : 4).

Everything here is deterministic (no randomness). The verifier relies on this: it has to reproduce
the signer's normalized curves (SQIsign._normalize_curve) and the deterministic torsion bases
exactly, and it only ever needs x-coordinates: the kernel of a 2^e-isogeny is <P + c*Q> or
<c*P + Q>, which the 3-point ladder computes from x(P), x(Q), x(P - Q); the cyclicity witness
2^(e-1) Q only has to be tested for being the point at infinity (Z == 0).

2-isogeny steps use the formulas of Renes (https://ia.cr/2017/1198, as in SIKE) for kernels
(x2, 0) with x2 != 0. A kernel <(0, 0)> can only occur at the first step of a chain (under Renes'
map the images of the other two 2-torsion points, i.e. the next kernels, are never (0, 0)); it is
handled by isogeny_2_special with the translated model, at the cost of one square root.
"""
from utilities.strategy import optimised_strategy_old
from util import fp2_order_key


def deterministic_sqrt(a):
    """
    Square root of a in F_{p^2} = F_p(i), p = 3 mod 4, by two exponentiations in F_p
    (Sage's GF(p^2).sqrt() is several hundred times slower); the root that is minimal under
    fp2_order_key, exactly like util.deterministic_sqrt. Raises ValueError if a is not a square.
    """
    F = a.parent()
    p = F.characteristic()
    x0, x1 = (a.list() + [0, 0])[:2]

    def sqrt_Fp(x):
        r = x ** ((p + 1) // 4)
        return r if r * r == x else None

    if x1 == 0:
        y0 = sqrt_Fp(x0)
        r = F(y0) if y0 is not None else F.gen() * sqrt_Fp(-x0)     # -1 is not a square in F_p
    else:
        sd = sqrt_Fp(x0 * x0 + x1 * x1)                               # norm of a, a square iff a is
        if sd is None:
            raise ValueError("element is not a square")
        y0 = sqrt_Fp((x0 + sd) / 2)
        if y0 is None:
            y0 = sqrt_Fp((x0 - sd) / 2)
        r = F(y0) + F(x1 / (2 * y0)) * F.gen()
    return min(r, -r, key=fp2_order_key)


# ---- Kummer line arithmetic ------------------------------------------------------------------

def A24_affine(A):
    """(A + 2) / 4"""
    return (A + 2) / 4


def A24_projective(A):
    """(A24p : C24) = (A + 2 : 4)"""
    return A + 2, 4


def xDBL(P, A24p, C24):
    """[2]P for (A24p : C24) = (A + 2 : 4). 4M + 2S."""
    X, Z = P
    t0 = (X - Z) * (X - Z)
    t1 = (X + Z) * (X + Z)
    Z2 = C24 * t0
    X2 = Z2 * t1
    t1 = t1 - t0                    # 4XZ
    return X2, (Z2 + A24p * t1) * t1


def xDBLe(P, A24p, C24, n):
    """[2^n]P"""
    for _ in range(n):
        P = xDBL(P, A24p, C24)
    return P


def xDBLADD(P, Q, PmQ, A24):
    """([2]P, P + Q) from P, Q and PmQ = P - Q (all projective) and the affine A24 = (A + 2) / 4."""
    XP, ZP = P
    XQ, ZQ = Q
    XD, ZD = PmQ
    t0 = XP + ZP
    t1 = XP - ZP
    X2P = t0 * t0
    Z2P = t1 * t1
    t0 = t0 * (XQ - ZQ)
    t1 = t1 * (XQ + ZQ)
    t2 = X2P - Z2P
    ZQP = (t0 - t1) * (t0 - t1)
    XQP = (t0 + t1) * (t0 + t1)
    return (X2P * Z2P, (A24 * t2 + Z2P) * t2), (ZD * XQP, XD * ZQP)


def ladder3pt(m, xP, xQ, xPQ, A24, nbits=None):
    """
    x(P + [m]Q) from the affine x-coordinates xP, xQ and xPQ = x(P - Q) (= x(Q - P)),
    with A24 = (A + 2) / 4 affine. nbits fixes the number of ladder steps (m < 2^nbits).
    """
    R0, R1, R2 = (xQ, 1), (xP, 1), (xPQ, 1)
    m = int(m)
    if nbits is None:
        nbits = m.bit_length()
    for i in range(nbits):
        if (m >> i) & 1:
            R0, R1 = xDBLADD(R0, R1, R2, A24)
        else:
            R0, R2 = xDBLADD(R0, R2, R1, A24)
    return R1


# ---- 2-isogenies -------------------------------------------------------------------------------

def isogeny_2_codomain(T):
    """(A24p : C24) of the codomain of the 2-isogeny with kernel T = (X2 : Z2), X2 != 0 (Renes)."""
    X2, Z2 = T
    C24 = Z2 * Z2
    return C24 - X2 * X2, C24


def isogeny_2_eval(T, P):
    """image of P = (X : Z) under the 2-isogeny with kernel T = (X2 : Z2), X2 != 0 (Renes). 4M."""
    X2, Z2 = T
    X, Z = P
    t0 = (X2 + Z2) * (X - Z)
    t1 = (X2 - Z2) * (X + Z)
    return X * (t0 + t1), Z * (t0 - t1)


def isogeny_2_special(A, R4):
    """
    The 2-isogeny with kernel <(0, 0)> on E_A, written with a Montgomery codomain.

    R4 = (X4 : Z4) is a point of order 4 above (0, 0), i.e. the next kernel point; it has
    x(R4) = +-1 and the translation is chosen so that its image is not (0, 0) again:
        x(R4) = -1 : x -> (x - 1)^2 / (2 s x),  s = sqrt(A + 2),  A' = (A + 6) / (2 s)
        x(R4) = +1 : x -> (x + 1)^2 / (2 s x),  s = sqrt(2 - A),  A' = (A - 6) / (2 s)
    Returns ((A24p : C24) of the codomain, evaluation function on (X : Z)).
    """
    X4, Z4 = R4
    if X4 == Z4:
        s = deterministic_sqrt(2 - A)
        A24p, sign = -(s - 2) * (s - 2), 1
    else:
        assert X4 == -Z4, "R4 is not a point of order 4 above (0, 0)"
        s = deterministic_sqrt(A + 2)
        A24p, sign = (s + 2) * (s + 2), -1

    def ev(P):
        X, Z = P
        t = X + sign * Z
        return t * t, 2 * s * X * Z

    return (A24p, 8 * s), ev


_STRATEGY_CACHE = {}

def strategy_2(n, mul_c=1.4):
    """
    Optimal strategy (SIDH spec Algorithm 60) for a 2^n chain with doubling cost mul_c relative
    to one isogeny evaluation (xDBL is 4M + 2S, isogeny_2_eval is 4M).
    """
    key = (n, mul_c)
    if key not in _STRATEGY_CACHE:
        _STRATEGY_CACHE[key] = optimised_strategy_old(n, mul_c=mul_c) if n >= 2 else []
    return _STRATEGY_CACHE[key]


def isogeny_chain_2e(A, xK, n, extra=(), strategy=None):
    """
    The 2^n-isogeny from E_A (A affine) with kernel <K>, where xK = (X : Z) is a kernel point of
    exact order 2^n, computed as a chain of 2-isogenies along an optimal strategy.

    Returns (A', images) with A' the affine Montgomery coefficient of the codomain (one of the six
    coefficients of its isomorphism class; the caller normalizes it) and images the list of the
    x-only images of the points in `extra`.
    """
    extra = list(extra)
    if n == 0:
        return A, extra
    A24p, C24 = A24_projective(A)
    strat = strategy_2(n) if strategy is None else strategy
    pts = []      # stack of kernel points still to be pushed through (as in SIDH)
    idx = []      # their "heights" (number of doublings applied)
    ind = 0
    P = xK
    sidx = 0
    R4 = None
    for row in range(n):
        target = n - 1 - row
        while ind < target:
            pts.append(P)
            idx.append(ind)
            m = strat[sidx]
            sidx += 1
            if row == 0 and ind + m == target:
                # keep the point of order 4 above the first kernel point for isogeny_2_special
                R4 = xDBLe(P, A24p, C24, m - 1)
                P = xDBL(R4, A24p, C24)
            else:
                P = xDBLe(P, A24p, C24, m)
            ind += m
        # P now has order 2
        if row == 0 and P[0] == 0:
            if R4 is None:
                R4 = (-1, 1)      # n == 1: there is no next kernel; any translation works
            (A24p, C24), ev = isogeny_2_special(A, R4)
        else:
            A24p, C24 = isogeny_2_codomain(P)
            ev = lambda Q, T=P: isogeny_2_eval(T, Q)
        pts = [ev(Q) for Q in pts]
        extra = [ev(Q) for Q in extra]
        if pts:
            P = pts.pop()
            ind = idx.pop()
    return 4 * A24p / C24 - 2, extra


# ---- Montgomery normalization on coefficients only ----------------------------------------------

def normalize_A(A):
    """
    Algorithm 1 (MontgomeryNormalize) of the SQIsign specification on the coefficient alone: the
    same computation as SQIsign._normalize_curve, returning (A', R, U2) where A' is the canonical
    coefficient of the isomorphism class of E_A and x -> U2 * (x + R) is the x-part of the
    isomorphism E_A -> E_A' used there (U2 = U^2, so no square root is needed for U).
    """
    A2 = A * A
    s = deterministic_sqrt(A2 - 4)
    u = (9 - A2) / 2
    t = (A2 * A - 3 * A) / (2 * s)
    Aprime = deterministic_sqrt(min((A2, u + t, u - t), key=fp2_order_key))
    if Aprime == A:
        return Aprime, 0, 1
    if Aprime == -A:
        return Aprime, 0, -1                    # U = i
    Ap2 = Aprime * Aprime
    R = (A2 + Ap2 - 6) * A / (A2 + 2 * Ap2 - 9)
    return Aprime, R, Aprime / (A - 3 * R)


def apply_isomorphism_x(P, R, U2):
    """(X : Z) -> (U2 * (X + R Z) : Z), the x-part of (x, y) -> (U^2 (x + R), U^3 y)"""
    X, Z = P
    return U2 * (X + R * Z), Z
