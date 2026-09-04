"""
x-only arithmetic on Montgomery curves and 2^n-isogeny chains over F_{p^2} = F_p(i), p = 3 mod 4.

Representation (chosen for speed: everything is plain gmpy2 integer arithmetic, no Sage objects):
  - an element of F_{p^2} is a pair (a0, a1) of gmpy2 integers meaning a0 + a1*i, with 0 <= a0, a1 < p;
  - a point on the Kummer line of E_A: y^2 = x^3 + A x^2 + x is a pair (X, Z) of field elements,
    Z == 0 meaning the point at infinity;
  - a curve is carried either as its affine coefficient A or projectively as (A24p : C24) = (A + 2 : 4).

Everything here is deterministic (no randomness, no PARI). The verifier relies on this: it has to
reproduce the signer's normalized curves (SQIsign._normalize_curve) and the deterministic torsion
bases exactly, and it only ever needs x-coordinates: the kernel of a 2^e-isogeny is <P + c*Q> or
<c*P + Q>, which the 3-point ladder computes from x(P), x(Q), x(P - Q); the cyclicity witness
2^(e-1) Q only has to be tested for being the point at infinity (Z == 0).

2-isogeny steps use the formulas of Renes (https://ia.cr/2017/1198, as in SIKE) for kernels
(x2, 0) with x2 != 0. A kernel <(0, 0)> can only occur at the first step of a chain (under Renes'
map the images of the other two 2-torsion points, i.e. the next kernels, are never (0, 0)); it is
handled by isogeny_2_special with the translated model, at the cost of one square root.
"""
import gmpy2
from gmpy2 import mpz
from utilities.strategy import optimised_strategy_old


class Fp2:
    """
    Arithmetic in F_{p^2} = F_p(i), i^2 = -1, on pairs of gmpy2 integers.

    The square root is the p = 3 mod 4 method of utilities/fast_sqrt.sqrt_Fp2 (two exponentiations
    in F_p instead of Tonelli-Shanks, whose cost is quadratic in v_2(p^2 - 1) = e + 1);
    deterministic_sqrt picks the root that is minimal in the lexicographic order of
    util.fp2_order_key, exactly like util.deterministic_sqrt.
    """
    __slots__ = ("p", "zero", "one", "two", "four", "_sqrt_exp", "_inv2", "_inv4")

    def __init__(self, p):
        p = mpz(p)
        assert p % 4 == 3
        self.p = p
        self.zero = (mpz(0), mpz(0))
        self.one = (mpz(1), mpz(0))
        self.two = (mpz(2), mpz(0))
        self.four = (mpz(4), mpz(0))
        self._sqrt_exp = (p + 1) // 4
        self._inv2 = (p + 1) // 2
        self._inv4 = gmpy2.invert(mpz(4), p)

    # ---- conversions -----------------------------------------------------------------------------
    def from_sage(self, z):
        """F_{p^2} element of Sage (a0 + a1*i, i = F.gen()) -> pair"""
        l = z.list()
        a0 = mpz(int(l[0])) if len(l) > 0 else mpz(0)
        a1 = mpz(int(l[1])) if len(l) > 1 else mpz(0)
        return (a0, a1)

    def to_sage(self, F, a):
        """pair -> element of the Sage field F (which must be F_p(i) with modulus x^2 + 1)"""
        return F([int(a[0]), int(a[1])])

    def from_int(self, n):
        return (mpz(n) % self.p, mpz(0))

    # ---- arithmetic ------------------------------------------------------------------------------
    def add(self, a, b):
        p = self.p
        return ((a[0] + b[0]) % p, (a[1] + b[1]) % p)

    def sub(self, a, b):
        p = self.p
        return ((a[0] - b[0]) % p, (a[1] - b[1]) % p)

    def neg(self, a):
        p = self.p
        return ((-a[0]) % p, (-a[1]) % p)

    def mul(self, a, b):
        p = self.p
        a0, a1 = a
        b0, b1 = b
        t0 = a0 * b0
        t1 = a1 * b1
        return ((t0 - t1) % p, ((a0 + a1) * (b0 + b1) - t0 - t1) % p)

    def sqr(self, a):
        p = self.p
        a0, a1 = a
        return (((a0 + a1) * (a0 - a1)) % p, (2 * a0 * a1) % p)

    def scale(self, a, n):
        """a * n for a small integer n"""
        p = self.p
        return ((a[0] * n) % p, (a[1] * n) % p)

    def inv(self, a):
        p = self.p
        a0, a1 = a
        n = (a0 * a0 + a1 * a1) % p
        ni = gmpy2.invert(n, p)
        return ((a0 * ni) % p, (-a1 * ni) % p)

    def is_zero(self, a):
        return a[0] == 0 and a[1] == 0

    def is_square(self, a):
        """Euler's criterion on the norm to F_p"""
        p = self.p
        n = (a[0] * a[0] + a[1] * a[1]) % p
        if n == 0:
            return True
        return gmpy2.powmod(n, (p - 1) // 2, p) == 1

    def _sqrt_Fp(self, x):
        """a square root of x in F_p (p = 3 mod 4), or None if x is not a square"""
        p = self.p
        r = gmpy2.powmod(x, self._sqrt_exp, p)
        if (r * r) % p != x % p:
            return None
        return r

    def sqrt(self, a):
        """some square root of a; raises ValueError if a is not a square"""
        p = self.p
        x0, x1 = a
        if x1 == 0:
            y0 = self._sqrt_Fp(x0)
            if y0 is not None:
                return (y0, mpz(0))
            y1 = self._sqrt_Fp((-x0) % p)
            if y1 is None:
                raise ValueError("element is not a square")
            return (mpz(0), y1)
        delta = (x0 * x0 + x1 * x1) % p       # norm of a, a square in F_p iff a is a square
        sd = self._sqrt_Fp(delta)
        if sd is None:
            raise ValueError("element is not a square")
        y02 = ((x0 + sd) * self._inv2) % p
        y0 = self._sqrt_Fp(y02)
        if y0 is None:
            y02 = (y02 - sd) % p
            y0 = self._sqrt_Fp(y02)
            if y0 is None:
                raise ValueError("element is not a square")
        y1 = (x1 * gmpy2.invert((2 * y0) % p, p)) % p
        return (y0, y1)

    def order_key(self, a):
        """the ordering of util.fp2_order_key: lexicographic on (a0, a1) in [0, p-1]"""
        return (int(a[0]), int(a[1]))

    def deterministic_sqrt(self, a):
        """min(r, -r) under order_key, identical to util.deterministic_sqrt"""
        r = self.sqrt(a)
        mr = self.neg(r)
        return min(r, mr, key=self.order_key)


# ---- Kummer line arithmetic ------------------------------------------------------------------

def A24_affine(K, A):
    """(A + 2) / 4"""
    a = K.add(A, K.two)
    return ((a[0] * K._inv4) % K.p, (a[1] * K._inv4) % K.p)


def A24_projective(K, A):
    """(A24p : C24) = (A + 2 : 4)"""
    return K.add(A, K.two), K.four


def x_affine(K, P):
    """X / Z (P must not be the point at infinity)"""
    X, Z = P
    return K.mul(X, K.inv(Z))


def same_x(K, P, Q):
    """whether (X:Z) and (X':Z') are the same point of the Kummer line"""
    return K.mul(P[0], Q[1]) == K.mul(Q[0], P[1])


def xDBL(K, P, A24p, C24):
    """[2]P for (A24p : C24) = (A + 2 : 4). 4M + 2S."""
    p = K.p
    X, Z = P
    x0, x1 = X
    z0, z1 = Z
    d0, d1 = x0 - z0, x1 - z1
    s0, s1 = x0 + z0, x1 + z1
    t0 = (((d0 + d1) * (d0 - d1)) % p, (2 * d0 * d1) % p)      # (X - Z)^2
    t1 = (((s0 + s1) * (s0 - s1)) % p, (2 * s0 * s1) % p)      # (X + Z)^2
    Z2 = K.mul(C24, t0)
    X2 = K.mul(Z2, t1)
    t1 = ((t1[0] - t0[0]) % p, (t1[1] - t0[1]) % p)            # 4XZ
    t0 = K.mul(A24p, t1)
    Z2 = ((Z2[0] + t0[0]) % p, (Z2[1] + t0[1]) % p)
    Z2 = K.mul(Z2, t1)
    return (X2, Z2)


def xDBLe(K, P, A24p, C24, n):
    """[2^n]P"""
    for _ in range(n):
        P = xDBL(K, P, A24p, C24)
    return P


def xDBLADD(K, P, Q, PmQ, A24):
    """
    ([2]P, P + Q) from P, Q and PmQ = P - Q (all projective) and the affine A24 = (A + 2) / 4.
    """
    mul, sqr = K.mul, K.sqr
    p = K.p
    XP, ZP = P
    XQ, ZQ = Q
    XD, ZD = PmQ
    t0 = ((XP[0] + ZP[0]) % p, (XP[1] + ZP[1]) % p)
    t1 = ((XP[0] - ZP[0]) % p, (XP[1] - ZP[1]) % p)
    X2P = sqr(t0)
    t2 = ((XQ[0] - ZQ[0]) % p, (XQ[1] - ZQ[1]) % p)
    XQP = ((XQ[0] + ZQ[0]) % p, (XQ[1] + ZQ[1]) % p)
    t0 = mul(t0, t2)
    Z2P = sqr(t1)
    t1 = mul(t1, XQP)
    t2 = ((X2P[0] - Z2P[0]) % p, (X2P[1] - Z2P[1]) % p)
    X2P = mul(X2P, Z2P)
    XQP = mul(A24, t2)
    ZQP = ((t0[0] - t1[0]) % p, (t0[1] - t1[1]) % p)
    Z2P = ((XQP[0] + Z2P[0]) % p, (XQP[1] + Z2P[1]) % p)
    XQP = ((t0[0] + t1[0]) % p, (t0[1] + t1[1]) % p)
    Z2P = mul(Z2P, t2)
    ZQP = sqr(ZQP)
    XQP = sqr(XQP)
    ZQP = mul(XD, ZQP)
    XQP = mul(ZD, XQP)
    return (X2P, Z2P), (XQP, ZQP)


def ladder3pt(K, m, xP, xQ, xPQ, A24, nbits=None):
    """
    x(P + [m]Q) from the affine x-coordinates xP, xQ and xPQ = x(P - Q) (= x(Q - P)),
    with A24 = (A + 2) / 4 affine. nbits fixes the number of ladder steps (m < 2^nbits).
    """
    one = K.one
    R0 = (xQ, one)
    R1 = (xP, one)
    R2 = (xPQ, one)
    m = int(m)
    if nbits is None:
        nbits = m.bit_length()
    for i in range(nbits):
        if (m >> i) & 1:
            R0, R1 = xDBLADD(K, R0, R1, R2, A24)
        else:
            R0, R2 = xDBLADD(K, R0, R2, R1, A24)
    return R1


# ---- 2-isogenies -------------------------------------------------------------------------------

def isogeny_2_codomain(K, T):
    """(A24p : C24) of the codomain of the 2-isogeny with kernel T = (X2 : Z2), X2 != 0 (Renes)."""
    X2, Z2 = T
    A24p = K.sqr(X2)
    C24 = K.sqr(Z2)
    A24p = K.sub(C24, A24p)
    return A24p, C24


def isogeny_2_eval(K, T, P):
    """image of P = (X : Z) under the 2-isogeny with kernel T = (X2 : Z2), X2 != 0 (Renes). 4M."""
    p = K.p
    X2, Z2 = T
    X, Z = P
    t0 = ((X2[0] + Z2[0]) % p, (X2[1] + Z2[1]) % p)
    t1 = ((X2[0] - Z2[0]) % p, (X2[1] - Z2[1]) % p)
    t2 = ((X[0] + Z[0]) % p, (X[1] + Z[1]) % p)
    t3 = ((X[0] - Z[0]) % p, (X[1] - Z[1]) % p)
    t0 = K.mul(t0, t3)
    t1 = K.mul(t1, t2)
    t2 = ((t0[0] + t1[0]) % p, (t0[1] + t1[1]) % p)
    t3 = ((t0[0] - t1[0]) % p, (t0[1] - t1[1]) % p)
    return (K.mul(X, t2), K.mul(Z, t3))


def isogeny_2_special(K, A, R4):
    """
    The 2-isogeny with kernel <(0, 0)> on E_A, written with a Montgomery codomain.

    R4 = (X4 : Z4) is a point of order 4 above (0, 0), i.e. the next kernel point; it has
    x(R4) = +-1 and the translation is chosen so that its image is not (0, 0) again:
        x(R4) = -1 : x -> (x - 1)^2 / (2 s x),  s = sqrt(A + 2),  A' = (A + 6) / (2 s)
        x(R4) = +1 : x -> (x + 1)^2 / (2 s x),  s = sqrt(2 - A),  A' = (A - 6) / (2 s)
    Returns ((A24p : C24) of the codomain, evaluation function on (X : Z)).
    """
    p = K.p
    X4, Z4 = R4
    if X4 == Z4:
        s = K.deterministic_sqrt(K.sub(K.two, A))
        A24p = K.neg(K.sqr(K.sub(s, K.two)))
        sign = 1
    else:
        assert X4 == K.neg(Z4), "R4 is not a point of order 4 above (0, 0)"
        s = K.deterministic_sqrt(K.add(A, K.two))
        A24p = K.sqr(K.add(s, K.two))
        sign = -1
    C24 = K.scale(s, 8)
    twos = K.scale(s, 2)

    def ev(P):
        X, Z = P
        t = K.add(X, Z) if sign == 1 else K.sub(X, Z)
        return (K.sqr(t), K.mul(twos, K.mul(X, Z)))

    return (A24p, C24), ev


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


def isogeny_chain_2e(K, A, xK, n, extra=(), strategy=None):
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
    A24p, C24 = A24_projective(K, A)
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
                P = xDBLe(K, P, A24p, C24, m - 1)
                R4 = P
                P = xDBL(K, P, A24p, C24)
            else:
                P = xDBLe(K, P, A24p, C24, m)
            ind += m
        # P now has order 2
        if row == 0 and K.is_zero(P[0]):
            if R4 is None:
                # n == 1: there is no next kernel; any translation works, pretend x(R4) = -1
                R4 = (K.neg(K.one), K.one)
            (A24p, C24), ev = isogeny_2_special(K, A, R4)
            pts = [ev(Q) for Q in pts]
            extra = [ev(Q) for Q in extra]
        else:
            A24p, C24 = isogeny_2_codomain(K, P)
            pts = [isogeny_2_eval(K, P, Q) for Q in pts]
            extra = [isogeny_2_eval(K, P, Q) for Q in extra]
        if pts:
            P = pts.pop()
            ind = idx.pop()
    Anew = K.sub(K.mul(K.four, K.mul(A24p, K.inv(C24))), K.two)      # A' = 4 A24p / C24 - 2
    return Anew, extra


# ---- Montgomery normalization on coefficients only ----------------------------------------------

def normalize_A(K, A):
    """
    Algorithm 1 (MontgomeryNormalize) of the SQIsign specification on the coefficient alone: the
    same computation as SQIsign._normalize_curve, returning (A', R, U2) where A' is the canonical
    coefficient of the isomorphism class of E_A and x -> U2 * (x + R) is the x-part of the
    isomorphism E_A -> E_A' used there (U2 = U^2, so no square root is needed for U).
    """
    p = K.p
    A2 = K.sqr(A)
    s = K.deterministic_sqrt(K.sub(A2, K.four))
    nine = K.from_int(9)
    u = K.sub(nine, A2)
    u = ((u[0] * K._inv2) % p, (u[1] * K._inv2) % p)                 # (9 - A^2) / 2
    t = K.mul(K.sub(K.mul(A2, A), K.scale(A, 3)), K.inv(K.scale(s, 2)))   # (A^3 - 3A) / (2s)
    Z = min((A2, K.add(u, t), K.sub(u, t)), key=K.order_key)
    Aprime = K.deterministic_sqrt(Z)
    if Aprime == A:
        R, U2 = K.zero, K.one
    elif Aprime == K.neg(A):
        R, U2 = K.zero, K.neg(K.one)                                 # U = i
    else:
        Ap2 = K.sqr(Aprime)
        R = K.mul(K.mul(K.sub(K.add(A2, Ap2), K.from_int(6)), A),
                  K.inv(K.sub(K.add(A2, K.add(Ap2, Ap2)), nine)))
        U2 = K.mul(Aprime, K.inv(K.sub(A, K.scale(R, 3))))
    return Aprime, R, U2


def apply_isomorphism_x(K, P, R, U2):
    """(X : Z) -> (U2 * (X + R Z) : Z), the x-part of (x, y) -> (U^2 (x + R), U^3 y)"""
    X, Z = P
    return (K.mul(U2, K.add(X, K.mul(R, Z))), Z)
