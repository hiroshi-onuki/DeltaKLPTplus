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
    proof,
)
from sage.rings.factorint import factor_trial_division
proof.all(False)

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

    M = ceil(O0.discriminant()**2 / N)
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
    N = norm(I)
    assert is_prime(N)
    assert kronecker(l**e, N) == kronecker(NCD, N)

    NCD = 0
    M = 3
    while not kronecker(l**e, M) == kronecker(NCD, M):
        J, _, M = EquivalentRandomPrimeIdeal(J, constraint=lambda N: N % 4 == 1)
        x, y = SumOf2Squares(M)
        I1 = I * (x + y*qi)
        I2 = I * J
        assert norm(I1) == norm(I2) == N*M
        beta1, beta2 = EquivalentIdealsWithSameNorm(I1, I2, N, M)
        print("newKLPT: N = %d, M = %d, NCD = %d" % (N, M, NCD))
        C_N, D_N = IdealModConstraint(I.left_order(), qj, qk, beta2, beta1, N)
        C_M, D_M = IdealModConstraint(I.left_order(), qj, qk, beta2, beta1, M)
        C, D = CRT([C_N, C_M], [N, M]), CRT([D_N, D_M], [N, M])
        NCD = p * (C**2 + D**2)
    nu = StrongApproximationTwoFactors(I.left_order(), N, M, C, D, l**e)
    assert beta2 * nu in I1
    return I1.intersection(I1.left_order()*nu), nu

# Given two O0-ideals I1, I2 with the same norm N*M, where N and M are primes,
# return beta1 in I1 and beta2 in I2 s.t. qI1(beta1) = qI2(beta2) approx p^(3/4) * N^(1/4)
def EquivalentIdealsWithSameNorm(I1, I2, N, M):
    assert I1.left_order() == I2.left_order()
    assert norm(I1) == norm(I2) == N*M
    assert is_prime(N) and is_prime(M)
    O0 = I1.left_order()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
    Q = O0.basis_matrix()
    Qinv = Q.inverse()
    ZN = ZZ.quotient_ring(ZZ(N))
    ZM = ZZ.quotient_ring(ZZ(M))

    # find vectors v1, v2, v3 corresponding to solusions of alpha2 * x * bar(alpha1) = 0 mod N*M
    alpha1 = SmallGenerator(I1)
    alpha2 = SmallGenerator(I2)
    MatN = matrix(ZN, Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv)
    MatM = matrix(ZM, Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv)
    kerN = MatN.left_kernel()
    kerM = MatM.left_kernel()
    assert kerN.dimension() == kerM.dimension() == 3
    v1N, v2N, v3N = kerN.basis()
    v1M, v2M, v3M = kerM.basis()
    v1 = vector(ZZ, [CRT([ZZ(cN), ZZ(cM)], [N, M]) for cN, cM in zip(v1N, v1M)])
    v2 = vector(ZZ, [CRT([ZZ(cN), ZZ(cM)], [N, M]) for cN, cM in zip(v2N, v2M)])
    v3 = vector(ZZ, [CRT([ZZ(cN), ZZ(cM)], [N, M]) for cN, cM in zip(v3N, v3M)])
    assert v1 * Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv % (N*M) == vector([0, 0, 0, 0])
    assert v2 * Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv % (N*M) == vector([0, 0, 0, 0])
    assert v3 * Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv % (N*M) == vector([0, 0, 0, 0])

    # (approximate) shortest solution for alpha2 * x * bar(alpha1) = 0 mod N
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
    L = IntegralLattice(Gram, [list(v) for v in [v1, v2, v3, vector([0,0,0,N*M])]])
    shortest_v = L.LLL().basis()[0]
    x = sum([c * b for c, b in zip(shortest_v, O0.basis())])

    # construct the lattice L = I1 \cap (O0 * x + Z)
    L1 = IntegralLattice(Gram, [vector(b) * Qinv for b in I1.basis()])
    Ox = IntegralLattice(Gram, [vector(b*x) * Qinv for b in O0.basis()])
    OxZ = Ox.overlattice([vector([1,0,0,0])])
    L = IntegralLattice(Gram, L1.intersection(OxZ).basis())

    bs = L.LLL().basis()
    newN = 1
    while not is_prime(newN):
        cs = [randint(-100, 100) for _ in range(len(bs))]
        v = sum(c * b for c, b in zip(cs, bs))
        beta1 = sum([c * b for c, b in zip(v, O0.basis())])
        newN = ZZ(beta1.reduced_norm() / (N*M))

    assert beta1 in I1
    beta2 = x * beta1 * x.conjugate() / x.reduced_norm()
    assert beta2 in I2

    return beta1, beta2