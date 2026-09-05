from sage.all import ( # type: ignore
    ZZ,
    pi,
    kronecker,
    ceil,
    floor,
    gcd,
    is_pseudoprime,
    norm,
    randint,
    sqrt,
    sum_of_k_squares,
    vector,
    matrix,
    IntegralLattice,
    log,
    RR,
    QQ,
)
from sage.rings.factorint import factor_trial_division # type: ignore
from util import LCG, deterministic_sqrt_mod
import lattice

# Trial-division bound used before the primality test of the remaining cofactor.
# A larger bound makes more integers recognizable as sums of two squares
# (100 -> 2^14 raises the success rate on ~185-bit inputs by about 1.8x) for a
# negligible cost (~0.06 ms per call).
_S2S_TRIAL_DIVISION_BOUND = 2**14

# return x, y s.t. n = x^2 + y^2, or None, None if no such x, y exist
def SumOf2Squares(n):
    if n < 0:
        return None, None
    if n == 0:
        return 0, 0
    if n == 1:
        return 1, 0

    factor = factor_trial_division(n, _S2S_TRIAL_DIVISION_BOUND)
    if factor and is_pseudoprime(factor[-1][0]):
        try:
            x, y = sum_of_k_squares(2, ZZ(n))
            return x, y
        except ValueError:
            return None, None
    return None, None

# LLL reduced basis of I as a left O-ideal
def LLLBasis(I):
    B = I.quaternion_algebra()
    O = I.left_order()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace()/2 for b1 in B.basis() for b2 in B.basis()])
    L = IntegralLattice(Gram, [list(b) for b in (O(2)*I).basis()])
    return [O(b/2) for b in L.LLL().basis()]

# return a s.t. I = O*a + O*nrd(I)
def SmallGenerator(I, bound=100, basis=None):
    # `basis` may be a precomputed LLLBasis(I); callers that sample many
    # generators of the same ideal (Qlapoti) pass it in to avoid recomputing LLL.
    if basis is None:
        basis = LLLBasis(I)
    a = 0
    n = 0
    N = norm(I)
    while gcd(n, N**2) != N:
        coeffs = [randint(-bound, bound) for _ in range(len(basis))]
        a = sum(c * b for c, b in zip(coeffs, basis))
        n = a.reduced_norm()
    return a

def SmallestGenerator(I):
    basis = LLLBasis(I)
    a = basis[0]
    assert gcd(a.reduced_norm(), norm(I)**2) == norm(I)
    return a

# return I*bar(beta)/norm(I)
def EquivalentIdeal(I, beta):
    assert beta in I
    return I * (beta.conjugate() / norm(I))

# return J ~ I with nrd(J) is prime
def EquivalentPrimeIdeal(I, omega, constraint=lambda N: True):
    O = I.left_order()
    N = norm(I)
    p = O.discriminant()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O.basis() for b2 in O.basis()])
    Q = O.basis_matrix()
    Qinv = Q.inverse()

    B = ceil(2*sqrt(2)/pi * sqrt(p * omega * log(2) * log(p)/2)) * N
    L = IntegralLattice(Gram, [vector(b) * Qinv for b in I.basis()])
    cond = lambda newN: is_pseudoprime(ZZ(newN/(2*N))) and constraint(ZZ(newN/(2*N)))
    vs = lattice.LatticeEnumeration(L, B, condition=cond, num_vectors=1)
    k = 0
    while not vs:
        # The search radius B is chosen for a "random" ideal class (Gaussian heuristic). If the class of I
        # contains an ideal of norm N1 < sqrt(p)/160 (probability ~ pi^2/(2*160^2) ~ 1/5000 for a random
        # class), the lattice is degenerate: b1 = i*b0, every vector inside the ball is (a+bi)*b0 with
        # N' = (a^2+b^2)*N1, and the next vectors have N' ~ p/(4*N1) > B/(2N). N1 composite then leaves no
        # prime candidate at all. Enlarge the radius and skip the sublattice Z[i]*b0 = <b0, b1>.
        k += 1
        vs = lattice.LatticeEnumeration(L, B * 2**k, condition=cond, num_vectors=1, skip_rank=2)
    v = vs[0]
    alpha = sum(c * b for c, b in zip(v, O.basis()))
    return EquivalentIdeal(I, alpha), alpha, ZZ(alpha.reduced_norm() // N)

# return J ~ I with small nrd(J)
def SmallestEquivalentIdeal(I):
    O0 = I.left_order()
    basis = LLLBasis(I)
    return EquivalentIdeal(I, basis[0]), basis[0]

# return gamma in O0 s.t. nrd(gamma) = n
def FullRepresentInteger(O0, n):
    p = O0.discriminant()
    assert n > p, "RepresentInteger requires n > p"

    B = floor(sqrt(4*n/p))
    while True:
        z = randint(-B, B)
        t = randint(-B, B)
        x, y = SumOf2Squares(4*n - p * (z**2 + t**2))
        if x is None or y is None:
            continue
        if (x - t) % 2 == (y - z) % 2 == 0:
            if gcd([(x-t)//2, (y-z)//2, z, t]) == 1:
                gamma = O0([x, y, z, t]) / 2
                assert gamma.reduced_norm() == n
                return gamma

# return gamma in O0 s.t. nrd(gamma) = n
def FullRepresentIntegerDeterministic(O0, n, seed):
    p = O0.discriminant()
    assert n > p, "RepresentInteger requires n > p"

    B = floor(sqrt(4*n/p))
    rng = LCG(seed)
    while True:
        z = rng.randint(-B, B)
        t = rng.randint(-B, B)
        x, y = SumOf2Squares(4*n - p * (z**2 + t**2))
        if x is None or y is None:
            continue
        if (x - t) % 2 == (y - z) % 2 == 0:
            if gcd([(x-t)//2, (y-z)//2, z, t]) == 1:
                gamma = O0([x, y, z, t]) / 2
                assert gamma.reduced_norm() == n
                return gamma

# return C, D s.t. gamma * (C*qj + D*qk) in O0*alpha + O0*N
def IdealModConstraint(O0, qj, qk, gamma, alpha, N):
    Q = O0.basis_matrix().inverse()
    v_gamma_qj = vector(gamma * qj) * Q % N
    v_gamma_qk = vector(gamma * qk) * Q % N
    M_alpha = (-alpha).matrix() * Q % N
    R = ZZ.quotient_ring(N)
    M = matrix(R, [v_gamma_qj, v_gamma_qk])
    M = M.stack(matrix(R, M_alpha))
    sol = M.left_kernel()[1]

    C, D = ZZ(sol[0]), ZZ(sol[1])
    if (C == 0 and D == 0) or (C**2 + D**2) % N == 0:
        raise ValueError("IdealModConstraint: no solution")
    assert (gamma * (C*qj + D*qk)) in O0.left_ideal([alpha, N])
    return C, D

# Return a quaternion nu such that nu = Cj + Dk mod N and Nrd(nu) = nrd.
def FullStrongApproximation(O0, N, C, D, nrd, max_cnt=1000, condition=lambda nu: True):
    p = O0.discriminant()
    Nrd_mu = ZZ(p * (C**2 + D**2))
    assert kronecker(Nrd_mu, N) == kronecker(nrd, N)
    lam = 2 * deterministic_sqrt_mod(nrd * Nrd_mu.inverse_mod(N), N)
    rhs = ZZ((4*nrd - lam**2 * Nrd_mu) / N)
    R = ZZ.quotient_ring(N)

    c = 1
    d = ZZ((R(rhs) / R(2*p*lam) - C*c) / R(D))

    x = ZZ(-R(C)/R(D))
    b0 = vector(ZZ, [N, N*x])
    b1 = vector(ZZ, [0, N**2])
    beta1, beta0 = lattice.ShortBasisDim2Euclidean(b0, b1)
    target = vector(ZZ, [-lam*C - N*c, -lam*D - N*d])
    bound = ZZ(floor(4*nrd / p))

    vs = lattice.EnumerateCloseVectorsDim2Euclidean(beta1, beta0, target, max_cnt, bound)
    for v in vs:
        Nc = N*c + v[0]
        Nd = N*d + v[1]

        tmp = ZZ((4*nrd - p*((lam*C + Nc)**2 + (lam*D + Nd)**2)) / N**2)
        a, b = SumOf2Squares(tmp)
        if a is not None and b is not None:
            nu = O0([N*a, N*b, lam*C + Nc, lam*D + Nd])
            if nu in O0*2:
                nu = nu / 2
                assert nu.reduced_norm() == nrd
                if condition(nu):
                    return nu, True
    return None, False

# return J ~ I with nrd(J) is l**2
def KLPT(I, l, e, omega):
    p = I.quaternion_algebra().discriminant()
    O = I.left_order()
    _, _, qj, qk = I.quaternion_algebra().basis()
    L, alpha, N = EquivalentPrimeIdeal(I, omega)
    beta = SmallestGenerator(L)

    B_FRI = ZZ(ceil(16 * omega/pi * log(2) * p * log(p)))
    e0 = ceil(log(B_FRI / N, l))
    e1 = e - e0
    num_vectors_SA = ceil(6*log(2) * omega * log(p))

    pCD = None
    seed = 0
    while pCD is None or kronecker(l**e1, N) != kronecker(pCD, N):
        gamma = FullRepresentIntegerDeterministic(L.left_order(), N * l**e0, seed)
        C, D = IdealModConstraint(L.left_order(), qj, qk, gamma, beta, N)
        pCD = p * (C**2 + D**2)
        seed += 1000
    nu, found = FullStrongApproximation(L.left_order(), N, C, D, l**e1, max_cnt=num_vectors_SA)
    assert found
    assert gamma * nu in L
    assert gamma * nu * alpha / N in I
    return EquivalentIdeal(L, gamma * nu), gamma * nu * alpha / N, True

# return a left O0-ideal of norm N
# Algorithm 3 in https://eprint.iacr.org/2024/760.pdf
def RandomFixedNormIdeal(O0, N):
    N = ZZ(N)
    if N <= 0:
        raise ValueError("N must be positive")

    if N == 1:
        return O0.left_ideal([O0(1)])

    M = ceil(O0.discriminant()**2 / N)
    while gcd(N, M) != 1:
        M += 1
    gamma = FullRepresentInteger(O0, N * M)

    basis = O0.basis()
    alpha = O0(0)

    while gcd(alpha.reduced_norm(), N) != 1:
        us = [randint(0, N - 1) for _ in range(4)]
        alpha = sum(u * b for u, b in zip(us, basis))
    I = O0.left_ideal([gamma*alpha, N])
    assert norm(I) == N
    return I, gamma*alpha

# ---- helpers for IdealNormReduce (pure basis-matrix arithmetic, no FreeModule/IntegralLattice objects) ----
_INR_CONST = RR(2*sqrt(2)/pi)
_order_cache = {}

# cached per-order data: (algebra, Gram matrix of the trace form on O0's basis, basis matrix, its inverse, basis, p)
def _order_data(O0):
    key = (O0.quaternion_algebra(), tuple(O0.basis()))
    d = _order_cache.get(key)
    if d is None:
        Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
        Q = O0.basis_matrix()
        d = (O0.quaternion_algebra(), Gram, Q, Q.inverse(), O0.basis(), ZZ(O0.discriminant()))
        _order_cache[key] = d
    return d

# intersection of two full-rank lattices given by 4x4 rational basis matrices (rows): A ∩ B = (A* + B*)*
def _lattice_intersection(A, B):
    Sz, d = A.inverse().transpose().stack(B.inverse().transpose())._clear_denom()
    H = Sz._hnf_pari(0, include_zero_rows=False)
    return (H / d).inverse().transpose()

# row HNF of a rational generating set (same basis as A.ideal(gens).basis_matrix())
def _hnf(V):
    Hz, d = V._clear_denom()
    H = Hz._hnf_pari(0, include_zero_rows=False)
    return H if d == 1 else H / d

# LLL-reduced basis (rows, O0-coordinates) of the lattice with integer basis M w.r.t. the reduced norm,
# together with its Gram matrix (entries tr(b_i * conj(b_j)) = 2 * <b_i, b_j>); the first row is identical to
# IntegralLattice(Gram, M.rows()).LLL().basis()[0]
def _lll_reduce(M, Gram):
    G = M * Gram * M.T
    U = G.LLL_gram().T
    Mred = U * M
    return Mred, Mred * Gram * Mred.T

# first vector of the LLL-reduced basis, as a quaternion
def _lll_shortest(M, Gram, basis):
    Mred, _ = _lll_reduce(M, Gram)
    return sum(c * b for c, b in zip(Mred[0], basis))

def IdealNormReduce(I1, I2, check=False, search_prime=True):
    """
    One norm-reduction step on a pair of left O0-ideals of the same norm N: returns (J1, J2, beta1, beta2, newN)
    with J1 = I1 * conj(beta1) / N, J2 = I2 * conj(beta2) / N of the same norm newN.

    beta1 is normally the LLL-shortest vector of the second lattice. With search_prime, when the Minkowski
    bound of that lattice already guarantees newN <= p (i.e. this is the last reduction step of
    DeltaKLPT_plus, which then needs newN to be prime), all lattice vectors with nrd/N <= p are enumerated
    and the shortest one with (pseudo)prime nrd/N is used instead; if there is none, the shortest vector is
    used as before.
    """
    N = norm(I1)
    assert norm(I2) == N
    O0 = I1.left_order()
    A, Gram, Q, Qinv, basis, p = _order_data(O0)
    B1 = I1.basis_matrix()
    B2 = I2.basis_matrix()

    # L = {x in O0 | I_2 x subseteq I_1}  (= (N*O0 ∩ conj(I2)*I1) / N)
    J = _hnf(matrix(QQ, [list(bc.conjugate() * b) for bc in I2.basis() for b in I1.basis()]))   # conj(I2)*I1
    L = _lattice_intersection(N * Q, J)
    M = (_hnf(L) * Qinv / N).change_ring(ZZ)          # LLL is scale invariant: reduce L/N instead of L
    x = _lll_shortest(M, Gram, basis)
    Nx = x.reduced_norm()
    assert RR(Nx) < _INR_CONST * RR(p * N).sqrt()
    assert Nx % p != 0, "Nx mod p = {}".format(Nx % p)

    # L = Nx * (I1 \cap x^{-1} * I2 * x)
    xc = x.conjugate()
    L = _lattice_intersection(Nx * B1, matrix(QQ, [list(xc * b * x) for b in I2.basis()]))
    M = (_hnf(L) * Qinv / Nx).change_ring(ZZ)
    Mred, Gred = _lll_reduce(M, Gram)
    beta1 = sum(c * b for c, b in zip(Mred[0], basis))
    newN = ZZ(beta1.reduced_norm() / N)
    if search_prime and _INR_CONST * RR(p * Nx).sqrt() < p:
        # Minkowski bound of this lattice is below p: last reduction step. Prefer the shortest vector whose
        # nrd/N is an odd (pseudo)prime among all vectors with nrd/N <= p (c * Gred * c^T = 2 * nrd).
        for val, c in lattice.ShortVectorsGram(Gred, 2 * p * N):
            cand = ZZ(ZZ(val) / (2 * N))      # N may be a Rational (norm of a fractional ideal)
            if cand % 2 == 0 or cand % p == 0 or not is_pseudoprime(cand):
                continue
            beta1 = sum(ci * b for ci, b in zip(vector(ZZ, c) * Mred, basis))
            assert beta1.reduced_norm() == cand * N
            newN = cand
            break
    assert RR(newN) < _INR_CONST * RR(N) * RR(p * Nx).sqrt()
    assert newN % p != 0
    beta2 = x * beta1 * xc / Nx
    if check:
        assert beta1 in I1
        assert beta2 in I2
    # EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2) with the left order cached
    J1 = A.ideal(list(_hnf(B1 * beta1.conjugate().matrix() / N)), left_order=O0, check=False)
    J2 = A.ideal(list(_hnf(B2 * beta2.conjugate().matrix() / N)), left_order=O0, check=False)
    return J1, J2, beta1, beta2, newN

def DeltaKLPT_plus(Icom, IskIchl, l, e, omega, count_iter=False):
    assert Icom.left_order() == IskIchl.left_order()
    _, qi, qj, qk = Icom.quaternion_algebra().basis()
    O = Icom.left_order()
    p = Icom.quaternion_algebra().discriminant()
    le = l**e

    # bound for the original KLPT
    B_KLPT = 96 * (log(2)/pi * omega * p * log(p))**3
    e_KLPT = ceil(log(B_KLPT, 2))

    J1, _, found1 = KLPT(Icom, 2, e_KLPT, omega)
    J2, _, found2 = KLPT(IskIchl, 2, e_KLPT, omega)
    assert found1 and found2
    assert norm(J1) == norm(J2) == 2**e_KLPT
    N = 2**e_KLPT

    iteration_count = 0
    while True:
        r = randint(0, p)
        if r == p:
            alpha = qi
        else:
            v, _ = lattice.ShortBasisDim2Euclidean(vector(ZZ, [1, r]), vector(ZZ, [0, p]))
            alpha = v[0] + v[1]*qi
        J1 = EquivalentIdeal(J1, N*alpha)
        J2 = EquivalentIdeal(J2, N*alpha.conjugate())
        N = N * ZZ(alpha.reduced_norm())

        max_iter = ceil(1/2 * log(log(RR(N/p) + 0.303, 2), 2) + 0.6) # the theoretical bound
        iter = 0
        while N > p:
            J1, J2, _, beta2, newN = IdealNormReduce(J1, J2)
            N = newN
            iteration_count += 1
            iter += 1
        assert iter <= max_iter, "Exceeded max iterations: {} > {}".format(iter, max_iter)

        def is_cyclic(nu):
            if nu / 2 in O:
                return False
            I = IskIchl.conjugate() * J1.intersection(O*nu)
            check, alpha = I.is_principal(True)
            assert check
            return (alpha / 2) not in O

        if is_pseudoprime(N):
            C, D = IdealModConstraint(O, qj, qk, SmallestGenerator(J2), SmallestGenerator(J1), N)
            if kronecker(l**e, N) == kronecker(p * (C**2 + D**2), N):
                num_vec = ceil(2*log(2) * omega * (l**2 + l - 1)/(l**2 - 1) * log(le/N**2))
                nu, found = FullStrongApproximation(O, N, C, D, le, num_vec, condition=is_cyclic)
                assert found
                if found:
                    beta2 = SmallestGenerator(J2)
                    assert beta2 * nu in J1
                    assert J1.intersection(O*nu) == J2 * nu
                    if count_iter:
                        return J1.intersection(O*nu), nu, iteration_count
                    return J1.intersection(O*nu), nu

def Qlapoti(I, e, max_tries=10**7):
    _, qi, _, _ = I.quaternion_algebra().basis()
    assert qi**2 == -1
    O = I.left_order()
    N = 2**e
    p = I.quaternion_algebra().discriminant()
    B = ceil(log(p))

    I, beta0 = SmallestEquivalentIdeal(I)
    n = norm(I)
    basis = LLLBasis(I)  # fixed for all tries; SmallGenerator only takes random combinations

    for _ in range(max_tries):
        alpha = SmallGenerator(I, B, basis)
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
        rb1, rb2 = lattice.ShortBasisDim2Euclidean(b1, b2)
        vs = lattice.EnumerateCloseVectorsDim2Euclidean(rb1, rb2, target, 10000, ceil(2 * (N - 2*r) / n))
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
            if (gamma - 1) / 2 in O or (gamma - qi) / 2 in O:
                continue
            beta1 = beta1 * beta0 / n
            beta2 = beta2 * beta0 / n
            return beta1, beta2, gamma
    raise ValueError("Qlapoti: no solution found after {} generators (nrd(I) = {}, 2^e = 2^{})".format(max_tries, n, e))
