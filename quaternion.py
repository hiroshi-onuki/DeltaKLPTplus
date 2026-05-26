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
    log,
)
from sage.rings.factorint import factor_trial_division
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
    if factor and is_prime(factor[-1][0]):
        try:
            x, y = sum_of_k_squares(2, n)
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

# return I*bar(beta)/norm(I)
def EquivalentIdeal(I, beta):
    assert beta in I
    return I * (beta.conjugate() / norm(I))

# return J ~ I with nrd(J) is prime
def EquivalentRandomPrimeIdeal(I, constraint=lambda N: True):
    O = I.left_order()
    basis = LLLBasis(I)
    N = 0
    a = O(0)
    while not (is_prime(N) and constraint(N)):
        cs = [randint(-100, 100) for _ in range(len(basis))]
        a = sum(c * b for c, b in zip(cs, basis))
        N = ZZ(a.reduced_norm() // norm(I))
    return EquivalentIdeal(I, a), a, N

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
def StrongApproximation(O0, N, C, D, nrd, max_cnt=1000):
    p = O0.discriminant()
    Nrd_mu = p * (C**2 + D**2)
    assert kronecker(Nrd_mu, N) == kronecker(nrd, N)
    lam = ZZ(sqrt(GF(N)(nrd)/GF(N)(Nrd_mu)))
    rhs = ZZ((nrd - lam**2 * Nrd_mu) / N)
    R = ZZ.quotient_ring(N)

    c = 1
    d = ZZ((R(rhs) / R(2*p*lam) - C*c) / R(D))

    x = ZZ(-R(C)/R(D))
    b0 = vector(ZZ, [N, N*x])
    b1 = vector(ZZ, [0, N**2])
    beta1, beta0 = lattice.ShortBasisDim2Euclidean(b0, b1)
    target = vector(ZZ, [-lam*C - N*c, -lam*D - N*d])
    close = lattice.ClosestVectorDim2Euclidean(beta1, beta0, target)
    bound = ZZ(floor(nrd / p))

    for v in lattice.EnumerateCloseVectorsDim2Euclidean(beta1, beta0, target, close, max_cnt, bound):
        Nc = N*c + v[0]
        Nd = N*d + v[1]

        tmp = ZZ((nrd - p*((lam*C + Nc)**2 + (lam*D + Nd)**2)) / N**2)
        a, b = SumOf2Squares(tmp)
        if a is not None and b is not None:
            nu = O0([N*a, N*b, lam*C + Nc, lam*D + Nd])
            assert nu.reduced_norm() == nrd
            return nu
    raise ValueError("StrongApproximation: no solution found after max_cnt tries")

# return J ~ I with nrd(J) is n1*n2
def KLPT(I, n1, n2):
    p = I.quaternion_algebra().discriminant()
    O = I.left_order()
    _, _, qj, qk = I.quaternion_algebra().basis()
    assert n1 > p**(0.5)
    assert n2 > p**(2.5)
    L, _, N = EquivalentRandomPrimeIdeal(I)
    alpha = SmallGenerator(L)

    pCD = None
    while pCD is None or kronecker(n2, N) != kronecker(pCD, N):
        gamma = FullRepresentInteger(L.left_order(), n1 * N)
        C, D = IdealModConstraint(L.left_order(), qj, qk, gamma, alpha, N)
        pCD = p * (C**2 + D**2)
    nu = StrongApproximation(L.left_order(), N, C, D, n2)
    assert gamma * nu in L
    return EquivalentIdeal(L, gamma * nu)

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

def newKLPT(I, J, l, e):
    assert I.right_order() == J.left_order()
    _, qi, qj, qk = I.quaternion_algebra().basis()
    p = I.quaternion_algebra().discriminant()
    N = ZZ(norm(I))
    assert is_prime(N)

    NCD = 0
    newN = 3
    while not kronecker(l**e, newN) == kronecker(NCD, newN):
        J, _, M = EquivalentRandomPrimeIdeal(J, constraint=lambda N: N % 4 == 1)
        x, y = SumOf2Squares(M)
        I1 = I * (x + y*qi)
        I2 = I * J
        assert norm(I1) == norm(I2) == N*M
        I1, I2, newN = EquivalentIdealsWithSameNorm(I1, I2, N, M)
        assert norm(I1) == norm(I2) == newN
        beta1 = SmallGenerator(I1)
        beta2 = SmallGenerator(I2)
        C, D = IdealModConstraint(I.left_order(), qj, qk, beta2, beta1, newN)
        NCD = p * (C**2 + D**2)
    print(float(log(newN, 2)))
    nu = StrongApproximation(I.left_order(), newN, C, D, l**e)
    assert beta2 * nu in I1
    return I1.intersection(I1.left_order()*nu), nu

# Given two O0-ideals I1, I2 with the same norm N,
# return beta1 in I1 and beta2 in I2 s.t. qI1(beta1) = qI2(beta2) approx p^(3/4) * N^(1/4)
def EquivalentIdealsWithSameNorm(I1, I2, N):
    assert I1.left_order() == I2.left_order()
    assert norm(I1) == norm(I2) == N
    O0 = I1.left_order()
    p = O0.discriminant()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
    Q = O0.basis_matrix()
    Qinv = Q.inverse()
    ZN = ZZ.quotient_ring(ZZ(N))

    # construct the lattice L consisting of vectors corresponding to solutions of alpha2 * x * bar(alpha1) = 0 mod N
    alpha1 = SmallGenerator(I1)
    alpha2 = SmallGenerator(I2)
    MatN = matrix(ZN, Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv)
    w = None
    i = -1
    for col in MatN.columns():
        for j in range(4):
            if gcd(ZZ(col[j]), N) == 1:
                w = col
                i = j
                break
        if w is not None:
            break
    assert w is not None
    Lbasis = []
    for j in range(4):
        v = vector(ZZ, [0, 0, 0, 0])
        if j != i:
            v[j] = 1
            v[i] = ZZ(-w[i].inverse() * w[j])
            assert vector(ZN, v).dot_product(w) == 0
        else:
            v[i] = N
        Lbasis.append(v)
    L = IntegralLattice(Gram, Lbasis)

    # short solution for alpha2 * x * bar(alpha1) = 0 mod N with gcd(norm(x), N) = 1
    found = False
    e = 0
    while not found:
        v, _, found = lattice.LatticeEnumeration(L, ceil(log(p*N, 2)/2) + e, condition=lambda newN: gcd(newN/2, N) == 1)
        e += 1
    x = sum(c * b for c, b in zip(v, O0.basis()))
    Nx = x.reduced_norm()
    assert alpha2 * x * alpha1.conjugate() in O0 * N
    assert gcd(Nx, N) == 1

    # construct the lattice L = I1 \cap (O0 * x + Z)
    L1 = IntegralLattice(Gram, [vector(b) * Qinv for b in I1.basis()])
    Ox = IntegralLattice(Gram, [vector(b*x) * Qinv for b in O0.basis()])
    OxZ = Ox.overlattice([vector([1,0,0,0])])
    L = IntegralLattice(Gram, L1.intersection(OxZ).basis())

    # find a short vector v in L s.t. the normalized norm of the corresponding element is prime
    found = False
    e = 0
    while not found:
        v, newN, found = lattice.LatticeEnumeration(L, ceil(log(p*N**2*Nx, 2)/2) + e, condition=lambda newN: is_prime(ZZ(newN/(2*N))))
        e += 1
    beta1 = sum(c * b for c, b in zip(v, O0.basis()))
    newN = ZZ(newN / (2*N))

    assert beta1 in I1
    beta2 = x * beta1 * x.conjugate() / Nx
    assert beta2 in I2

    return EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2), newN

# Given two O0-ideals I1, I2 with the same norm N,
# return beta1 in I1 and beta2 in I2 s.t. qI1(beta1) = qI2(beta2) approx p^(1/2) * N^(1/2)
def EquivalentIdealsWithSameNormSmallN(I1, I2, N):
    assert I1.left_order() == I2.left_order()
    assert norm(I1) == norm(I2) == N
    O0 = I1.left_order()
    p = O0.discriminant()
    _, qi, _, _ = O0.quaternion_algebra().basis()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O0.basis() for b2 in O0.basis()])
    Q = O0.basis_matrix()
    Qinv = Q.inverse()
    ZN = ZZ.quotient_ring(ZZ(N))

    # construct the lattice L consisting of vectors corresponding to solutions of alpha2 * x * bar(alpha1) = 0 mod N
    alpha1 = SmallGenerator(I1)
    alpha2 = SmallGenerator(I2)
    MatN = matrix(ZN, Q * alpha1.conjugate().matrix('right') * alpha2.matrix('left') * Qinv)
    MatN = MatN[:2,:]
    b0 = None
    b1 = None
    for i in range(2):
        for j in range(4):
            if gcd(ZZ(MatN[i, j]), N) == 1:
                c = MatN[i, j].inverse() * MatN[1-i, j]
                b0 = vector(ZZ, [1, 1])
                b1 = vector(ZZ, [0, 0])
                b0[i] = -c
                b1[i] = N
                assert vector(ZN, b0) * MatN == 0
                break
        if b0 is not None:
            break
    assert b0 is not None
    L = IntegralLattice(matrix([[1, 0], [0, 1]]), [b0, b1])
    e = 0
    found = False
    while not found:
        v, Nx, found = lattice.LatticeEnumeration(L, ceil(log(N, 2)/2) + e, condition=lambda newN: gcd(newN, N) == 1)
        e += 1
    x = v[0] + v[1]*qi
    assert x.reduced_norm() == Nx
    assert alpha2 * x * alpha1.conjugate() in O0 * N

    # construct the lattice L = I1 \cap (O0 * x + Z)
    L1 = IntegralLattice(Gram, [vector(b) * Qinv for b in I1.basis()])
    Ox = IntegralLattice(Gram, [vector(b*x) * Qinv for b in O0.basis()])
    OxZ = Ox.overlattice([vector([1,0,0,0])])
    L = IntegralLattice(Gram, L1.intersection(OxZ).basis())

    # find a short vector v in L s.t. the normalized norm of the corresponding element is prime
    found = False
    e = 0
    while not found:
        v, newN, found = lattice.LatticeEnumeration(L, ceil(log(p*N**2*Nx, 2)/2) + e, condition=lambda newN: is_prime(ZZ(newN/(2*N))))
        e += 1
    beta1 = sum(c * b for c, b in zip(v, O0.basis()))
    newN = ZZ(newN / (2*N))

    assert beta1 in I1
    beta2 = x * beta1 * x.conjugate() / Nx
    assert beta2 in I2

    return EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2), newN


def deltaKLPT(I1, I2, l, e):
    assert I1.left_order() == I2.left_order()
    N = norm(I1)
    assert norm(I2) == N
    _, _, qj, qk = I1.quaternion_algebra().basis()
    O = I1.left_order()
    p = I1.quaternion_algebra().discriminant()

    C, D = 0, 0
    NCD = None
    while NCD is None or kronecker(l**e, N) != kronecker(NCD, N):
        I1, I2, N = EquivalentIdealsWithSameNorm(I1, I2, N)
        assert norm(I1) == norm(I2) == N
        beta1 = SmallGenerator(I1)
        beta2 = SmallGenerator(I2)
        C, D = IdealModConstraint(O, qj, qk, beta2, beta1, N)
        NCD = p * (C**2 + D**2)
    print(float(log(N, 2)))
    nu = StrongApproximation(O, N, C, D, l**e)
    assert beta2 * nu in I1
    return I1.intersection(O*nu), nu