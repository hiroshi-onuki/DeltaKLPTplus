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
    proof,
)
from sage.rings.factorint import factor_trial_division
proof.all(False)

def _round_div(n, d):
    n = ZZ(n)
    d = ZZ(d)
    if d < 0:
        n = -n
        d = -d
    if n >= 0:
        return (2*n + d) // (2*d)
    return -((-2*n + d) // (2*d))

def _floor_div(n, d):
    n = ZZ(n)
    d = ZZ(d)
    if d < 0:
        n = -n
        d = -d
    return n // d

def _ceil_div(n, d):
    return -_floor_div(-ZZ(n), ZZ(d))

def _floor_sqrt(n):
    n = ZZ(n)
    if n < 0:
        raise ValueError("square root of a negative integer")
    return floor(sqrt(n))

def _vec2(v):
    v = vector(ZZ, v)
    if len(v) != 2:
        raise ValueError("expected a 2-dimensional vector")
    return v

def _dot(u, v):
    return ZZ(u[0])*ZZ(v[0]) + ZZ(u[1])*ZZ(v[1])

def EuclideanNorm(v):
    v = _vec2(v)
    return _dot(v, v)

def ShortBasisEuclidean(b0, b1):
    beta0 = _vec2(b0)
    beta1 = _vec2(b1)
    if EuclideanNorm(beta0) < EuclideanNorm(beta1):
        beta0, beta1 = beta1, beta0

    gamma = beta0
    while True:
        r = _round_div(_dot(beta0, beta1), EuclideanNorm(beta1))
        gamma = beta0 - r*beta1
        if EuclideanNorm(gamma) < EuclideanNorm(beta1):
            beta0, beta1 = beta1, gamma
        else:
            break

    if EuclideanNorm(gamma) < EuclideanNorm(beta0):
        beta0 = gamma
    return beta1, beta0

def ClosestVectorEuclidean(beta1, beta0, t):
    beta1 = vector(ZZ, beta1)
    beta0 = vector(ZZ, beta0)
    t = vector(ZZ, t)
    N1 = beta1[0]**2 + beta1[1]**2
    B = beta1[0]*beta0[0] + beta1[1]*beta0[1]
    mu = N1 * beta0 - B * beta1
    Nmu = mu[0]**2 + mu[1]**2
    c = t - round(B*N1 / Nmu) * beta0
    B = beta1[0]*c[0] + beta1[1]*c[1]
    c = c - round(B / N1) * beta1
    return t - c

def EnumerateCloseVectorsEuclidean(L, t, close, m, B):
    b0, b1 = [_vec2(b) for b in L]
    t = _vec2(t)
    close = _vec2(close)
    m = ZZ(m)
    B = ZZ(B)
    if m <= 0:
        return

    d = t - close
    a = EuclideanNorm(b0)
    h = _dot(b0, b1)
    c = EuclideanNorm(b1)
    delta = a*c - h**2
    if delta <= 0:
        raise ValueError("EnumerateCloseVectorsEuclidean: degenerate lattice basis")

    det = ZZ(b0[0])*ZZ(b1[1]) - ZZ(b0[1])*ZZ(b1[0])
    y_num = ZZ(b0[0])*ZZ(d[1]) - ZZ(b0[1])*ZZ(d[0])
    y_den = det
    if y_den < 0:
        y_num = -y_num
        y_den = -y_den

    # From min_x ||d - x*b0 - y*b1||^2 = delta/a * (y-y0)^2.
    y_radius = _floor_sqrt((a*B) // delta) + 2
    y_min = _floor_div(y_num, y_den) - y_radius
    y_max = _ceil_div(y_num, y_den) + y_radius

    tries = ZZ(0)
    db0 = _dot(d, b0)
    db1 = _dot(d, b1)
    nd = EuclideanNorm(d)
    for y in range(y_min, y_max + 1):
        if tries >= m:
            break
        K = c*y**2 - 2*db1*y + nd - B
        Lx = h*y - db0
        D = Lx**2 - a*K
        if D < 0:
            continue
        x_radius = _floor_sqrt(D) + 2
        x_min = _floor_div(-Lx - x_radius, a) - 1
        x_max = _ceil_div(-Lx + x_radius, a) + 1
        for x in range(x_min, x_max + 1):
            if tries >= m:
                break
            tries += 1
            v = close + x*b0 + y*b1
            if EuclideanNorm(t - v) <= B:
                yield v

def EnumerateCloseVectors(L, t, close, m, B):
    yield from EnumerateCloseVectorsEuclidean(L, t, close, m, B)

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
    beta1, beta0 = ShortBasisEuclidean(b0, b1)
    target = vector(ZZ, [-lam*C - N*c, -lam*D - N*d])
    close = ClosestVectorEuclidean(beta1, beta0, target)
    bound = ZZ(floor(nrd / p))

    for v in EnumerateCloseVectorsEuclidean((beta1, beta0), target, close, max_cnt, bound):
        Nc = N*c + v[0]
        Nd = N*d + v[1]

        tmp = ZZ((nrd - p*((lam*C + Nc)**2 + (lam*D + Nd)**2)) / N**2)
        a, b = SumOf2Squares(tmp)
        if a is not None and b is not None:
            nu = O0([N*a, N*b, lam*C + Nc, lam*D + Nd])
            assert nu.reduced_norm() == nrd
            return nu
    raise ValueError("StrongApproximation: no solution found after max_cnt tries")

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
    nx = 0
    while gcd(nx, N*M) != 1:
        bs = L.LLL().basis()
        cs = [randint(-100, 100) for _ in range(len(bs))]
        v = sum(c * b for c, b in zip(cs, bs))
        x = sum([c * b for c, b in zip(v, O0.basis())])
        nx = x.reduced_norm()
    assert alpha2 * x * alpha1.conjugate() in O0 * (N*M)

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
    beta2 = x * beta1 * x.conjugate() / nx
    assert beta2 in I2

    return EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2), newN
