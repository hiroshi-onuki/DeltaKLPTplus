from sage.all import (
    ZZ,
    GF,
    kronecker,
    ceil,
    floor,
    gcd,
    is_prime,
    norm,
    randint,
    sqrt,
    sum_of_k_squares,
    vector,
    matrix,
    CRT,
    IntegralLattice,
)
from sage.rings.factorint import factor_trial_division

def SumOf2Squares(n):
    if n < 0:
        return None, None
    if n == 0:
        return 0, 0
    if n == 1:
        return 1, 0

    factor = factor_trial_division(n, 100)
    if factor and is_prime(factor[-1][0]):
        try:
            x, y = sum_of_k_squares(2, n)
            return x, y
        except ValueError:
            return None, None
    return None, None

def LLLBasis(I):
    B = I.quaternion_algebra()
    O = I.left_order()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace()/2 for b1 in B.basis() for b2 in B.basis()])
    L = IntegralLattice(Gram, [list(b) for b in (O(2)*I).basis()])
    return [O(b/2) for b in L.LLL().basis()]

def SmallGenerator(I):
    basis = LLLBasis(I)
    return basis[0]

def EquivalentRandomPrimeIdeal(I, constraint=lambda N: True):
    O = I.left_order()
    basis = LLLBasis(I)
    N = 0
    a = O(0)
    while not (is_prime(N) and constraint(N)):
        cs = [randint(-100, 100) for _ in range(len(basis))]
        a = sum(c * b for c, b in zip(cs, basis))
        N = ZZ(a.reduced_norm() // norm(I))
    return I * (a.conjugate() / norm(I)), a, N

# return gamma in O0 s.t. nrd(gamma) = n
def RepresentInteger(O0, n):
    p = O0.discriminant()
    assert n > p, "RepresentInteger requires n > p"

    B = floor(sqrt(n/p))
    while True:
        z = randint(-B, B)
        t = randint(-B, B)
        x, y = SumOf2Squares(n - p * (z**2 + t**2))
        if x is None or y is None:
            continue
        gamma = O0([x, y, z, t])
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

    print("IdealModConstraint: sol", M.left_kernel())
    sol = M.left_kernel()[1]

    C, D = ZZ(sol[0]), ZZ(sol[1])
    if (C == 0 and D == 0) or (C**2 + D**2) % N == 0:
        raise ValueError("IdealModConstraint: no solution")
    assert (gamma * (C*qj + D*qk)) in O0.left_ideal([alpha, N])
    return C, D

# Return a quaternion nu such that nu = Cj + Dk mod N*M and Nrd(nu) = nrd.
def StrongApproximationTwoFactors(O0, N, M, C, D, nrd, max_cnt=10):
    p = O0.discriminant()
    Nrd_mu = p * (C**2 + D**2)
    assert kronecker(Nrd_mu, N) == kronecker(nrd, N) and kronecker(Nrd_mu, M) == kronecker(nrd, M)
    lam_N = ZZ(sqrt(GF(N)(nrd)/GF(N)(Nrd_mu)))
    lam_M = ZZ(sqrt(GF(M)(nrd)/GF(M)(Nrd_mu)))
    lam = CRT([lam_N, lam_M], [N, M])
    rhs = ZZ((nrd - lam**2 * Nrd_mu) / (N*M))
    if rhs < 0:
        raise ValueError("StrongApproximationTwoFactors: no solution")
    R = ZZ.quotient_ring(N*M)

    for _ in range(max_cnt):
        c = 1
        d = ZZ((R(rhs) / R(2*p*lam) - C*c) / R(D))

        x = ZZ(-R(C)/R(D))
        L = IntegerLattice([[N*M, N*M*x], [0, (N*M)**2]])
        v = L.approximate_closest_vector([-lam*C - N*M*c, -lam*D - N*M*d])
        NMc = N*M*c + v[0]
        NMd = N*M*d + v[1]

        tmp = ZZ((nrd - p*((lam*C + NMc)**2 + (lam*D + NMd)**2)) / (N*M)**2)
        a, b = SumOf2Squares(tmp)
        if a is not None and b is not None:
            nu = O0([N*M*a, N*M*b, lam*C + NMc, lam*D + NMd])
            assert nu.reduced_norm() == nrd
            return nu

# return a left O0-ideal of norm N
# Algorithm 3 in https://eprint.iacr.org/2024/760.pdf
def RandomFixedNormIdeal(O0, N):
    N = ZZ(N)
    if N <= 0:
        raise ValueError("N must be positive")

    if N == 1:
        return O0.left_ideal([O0(1)])

    M = ceil(10*O0.discriminant() / N)
    gamma = RepresentInteger(O0, N * M)

    basis = O0.basis()
    alpha = O0(0)

    while gcd(alpha.reduced_norm(), N) != 1:
        us = [randint(0, N - 1) for _ in range(4)]
        alpha = sum(u * b for u, b in zip(us, basis))
    I = O0.left_ideal([gamma*alpha, N])
    assert norm(I) == N
    return I, gamma*alpha

def newKLPT(I, J, l, e):
    assert I.right_order() == J.left_order()
    _, qi, qj, qk = I.quaternion_algebra().basis()
    p = I.quaternion_algebra().discriminant()

    NCD = 0
    N = 3
    M = 3
    while not (kronecker(l**e, N) == kronecker(NCD, N) and kronecker(l**e, M) == kronecker(NCD, M)):
        I, alpha, N = EquivalentRandomPrimeIdeal(I)
        J = alpha * J * alpha.inverse()
        J, alpha, M = EquivalentRandomPrimeIdeal(J, constraint=lambda N: N % 4 == 1)
        x, y = SumOf2Squares(M)
        I1 = I * (x + y*qi)
        I2 = I * J
        beta1 = SmallGenerator(I1)
        beta2 = SmallGenerator(I2)
        assert norm(I1) == norm(I2) == N*M
        print("newKLPT: N = %d, M = %d, NCD = %d" % (N, M, NCD))
        C_N, D_N = IdealModConstraint(I.left_order(), qj, qk, beta2, beta1, N)
        C_M, D_M = IdealModConstraint(I.left_order(), qj, qk, beta2, beta1, M)
        C, D = CRT([C_N, C_M], [N, M]), CRT([D_N, D_M], [N, M])
        NCD = p * (C**2 + D**2)
    nu = StrongApproximationTwoFactors(I.left_order(), N, M, C, D, l**e)
    assert beta2 * nu in I1
    return I1.intersection(I1.left_order()*nu), nu