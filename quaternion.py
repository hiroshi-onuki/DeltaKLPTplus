from xml.etree.ElementTree import PI

from sage.all import ( # type: ignore
    ZZ,
    GF,
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
    CRT,
    IntegralLattice,
    log,
    set_random_seed,
)
from sage.rings.factorint import factor_trial_division # type: ignore
from util import deterministic_sqrt_mod
import lattice

# return x, y s.t. n = x^2 + y^2, or None, None if no such x, y exist
def SumOf2Squares(n):
    if n < 0:
        return None, None
    if n == 0:
        return 0, 0
    if n == 1:
        return 1, 0

    factor = factor_trial_division(n, 100)
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
def SmallGenerator(I):
    basis = LLLBasis(I)
    a = 0
    n = 0
    N = norm(I)
    while gcd(n, N**2) != N:
        coeffs = [randint(-100, 100) for _ in range(len(basis))]
        a = sum(c * b for c, b in zip(coeffs, basis))
        n = a.reduced_norm()
    return a

def SmallestGenerator(I):
    basis = LLLBasis(I)
    a = basis[0]
    assert gcd(a.reduced_norm(), norm(I)) == norm(I)
    return a

# return I*bar(beta)/norm(I)
def EquivalentIdeal(I, beta):
    assert beta in I
    return I * (beta.conjugate() / norm(I))

# return J ~ I with nrd(J) is prime
def EquivalentPrimeIdeal(I, constraint=lambda N: True):
    O = I.left_order()
    N = norm(I)
    p = O.discriminant()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O.basis() for b2 in O.basis()])
    Q = O.basis_matrix()
    Qinv = Q.inverse()

    omega = 64
    B = ceil(2*sqrt(2)/pi * sqrt(p * omega * log(2) * log(p)/2)) * N
    L = IntegralLattice(Gram, [vector(b) * Qinv for b in I.basis()])
    v = lattice.LatticeEnumeration(L, B, condition=lambda newN: is_pseudoprime(ZZ(newN/(2*N))) and constraint(ZZ(newN/(2*N))), num_vectors=1)[0]
    alpha = sum(c * b for c, b in zip(v, O.basis()))
    return EquivalentIdeal(I, alpha), alpha, ZZ(alpha.reduced_norm() // N)

# return J ~ I with small nrd(J)
def SmallestEquivalentIdeal(I):
    O0 = I.left_order()
    basis = LLLBasis(I)
    return EquivalentIdeal(I, basis[0]), basis[0]

# return gamma in O0 s.t. nrd(gamma) = n
def FullRepresentInteger(O0, n, seed=None):
    if seed is not None:
        set_random_seed(seed)
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
def KLPT(I, l, e):
    p = I.quaternion_algebra().discriminant()
    O = I.left_order()
    _, _, qj, qk = I.quaternion_algebra().basis()
    L, alpha, N = EquivalentPrimeIdeal(I)
    beta = SmallestGenerator(L)

    omega = 64
    B_FRI = ZZ(ceil(16 * omega/pi * log(2) * p * log(p)))
    e0 = ceil(log(B_FRI / N, l))
    e1 = e - e0
    num_vectors_SA = ceil(6*log(2) * omega * log(p))

    pCD = None
    seed = 0
    while pCD is None or kronecker(l**e1, N) != kronecker(pCD, N):
        gamma = FullRepresentInteger(L.left_order(), N * l**e0, seed=seed)
        C, D = IdealModConstraint(L.left_order(), qj, qk, gamma, beta, N)
        pCD = p * (C**2 + D**2)
        seed += 1
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

def IdealNormReduce(I1, I2):
    N = norm(I1)
    assert norm(I2) == N
    O0 = I1.left_order()
    p = O0.discriminant()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
    Q = O0.basis_matrix()
    Qinv = Q.inverse()

    # L = {x in O0 | I_2 x subseteq I_1}
    L = (O0*N).intersection(I2.conjugate() * I1)
    L = IntegralLattice(Gram, [vector(b) * Qinv for b in L.basis()])
    x = L.LLL().basis()[0]
    x = sum(c * b for c, b in zip(x, O0.basis()))
    x = x / N
    Nx = x.reduced_norm()
    assert Nx < 2*sqrt(2)/pi * sqrt(p * N)
    assert Nx % p != 0, "Nx mod p = {}".format(Nx % p)

    # L = Nx * (I1 \cap x^{-1} * I2 * x)
    L = (I1 * Nx).intersection(x.conjugate() * I2 * x)
    L = IntegralLattice(Gram, [vector(b) * Qinv for b in L.basis()])
    beta1 = L.LLL().basis()[0]
    beta1 = sum(c * b for c, b in zip(beta1, O0.basis())) / Nx
    newN = ZZ(beta1.reduced_norm() / N)
    assert newN < 2*sqrt(2)/pi * N * sqrt(p * Nx)
    assert newN % p != 0
    beta2 = x * beta1 * x.conjugate() / Nx
    assert beta1 in I1
    assert beta2 in I2
    return EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2), beta1, beta2, newN

def GeneralizedDeltaKLPT_heuristic(Icom, IskIchl, l, e, norm_bound):
    assert Icom.left_order() == IskIchl.left_order()
    _, qi, qj, qk = Icom.quaternion_algebra().basis()
    O = Icom.left_order()
    p = Icom.quaternion_algebra().discriminant()
    le = l**e

    # bound for the original KLPT
    omega = 64
    B_KLPT = 96 * (log(2)/pi * omega * p * log(p))**3
    e_KLPT = ceil(log(B_KLPT, 2))

    J1, _, found1 = KLPT(Icom, 2, e_KLPT)
    J2, _, found2 = KLPT(IskIchl, 2, e_KLPT)
    assert found1 and found2
    assert norm(J1) == norm(J2) == 2**e_KLPT
    N = 2**e_KLPT

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

        while N > norm_bound:
            J1, J2, _, beta2, newN = IdealNormReduce(J1, J2)
            N = newN

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
                    return J1.intersection(O*nu), nu

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
    raise ValueError("No suitable ideals found")
