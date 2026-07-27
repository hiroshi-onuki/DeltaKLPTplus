from sage.all import ( # type: ignore
    ZZ,
    GF,
    kronecker,
    ceil,
    floor,
    gcd,
    is_pseudoprime,
    random_prime,
    norm,
    randint,
    sqrt,
    sum_of_k_squares,
    vector,
    matrix,
    CRT,
    IntegralLattice,
    log,
    RealField,
)
from sage.rings.factorint import factor_trial_division # type: ignore
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
    while not (is_pseudoprime(N) and constraint(N)):
        cs = [randint(-100, 100) for _ in range(len(basis))]
        a = sum(c * b for c, b in zip(cs, basis))
        N = ZZ(a.reduced_norm() // norm(I))
    return EquivalentIdeal(I, a), a, N

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
def StrongApproximation(O0, N, C, D, nrd, max_cnt=1000, condition=lambda nu: True):
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
    bound = ZZ(floor(nrd / p))

    for v in lattice.EnumerateCloseVectorsDim2Euclidean(beta1, beta0, target, max_cnt, bound):
        Nc = N*c + v[0]
        Nd = N*d + v[1]

        tmp = ZZ((nrd - p*((lam*C + Nc)**2 + (lam*D + Nd)**2)) / N**2)
        a, b = SumOf2Squares(tmp)
        if a is not None and b is not None:
            nu = O0([N*a, N*b, lam*C + Nc, lam*D + Nd])
            assert nu.reduced_norm() == nrd
            if condition(nu):
                return nu, True
    return None, False

# return J ~ I with nrd(J) is n1*n2
def KLPT(I, n1, n2):
    p = I.quaternion_algebra().discriminant()
    O = I.left_order()
    _, _, qj, qk = I.quaternion_algebra().basis()
    assert n1 > p**(0.5)
    assert n2 > p**(2.5)
    L, alpha, N = EquivalentRandomPrimeIdeal(I)
    beta = SmallGenerator(L)

    pCD = None
    while pCD is None or kronecker(n2, N) != kronecker(pCD, N):
        gamma = FullRepresentInteger(L.left_order(), n1 * N)
        C, D = IdealModConstraint(L.left_order(), qj, qk, gamma, beta, N)
        pCD = p * (C**2 + D**2)
    nu, found = StrongApproximation(L.left_order(), N, C, D, n2)
    if not found:
        return None, None, False
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

# the norm bound p^(3/4) * N^(1/4) * (3/4*log(p) + 1/4*log(N)) reachable by
# EquivalentIdealsWithSameNorm. Evaluated over RealField because N may exceed the
# exponent range of a double.
def EquivalentIdealsWithSameNormBound(p, N, prec=100):
    RF = RealField(prec)
    return ceil(RF(p)**(RF(3)/4) * RF(N)**(RF(1)/4) * (3*RF(p).log() + RF(N).log()) / 4)

# Given two O0-ideals I1, I2 with the same norm N,
# return beta1 in I1 and beta2 in I2 s.t. qI1(beta1) = qI2(beta2) approx p^(3/4) * N^(1/4)
def EquivalentIdealsWithSameNorm(I1, I2, N, norm_bound, num_vectors=10):
    assert I1.left_order() == I2.left_order()
    assert norm(I1) == norm(I2) == N
    O0 = I1.left_order()
    p = O0.discriminant()
    assert norm_bound >= EquivalentIdealsWithSameNormBound(p, N)
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
    # Gram is the trace form, so an enumerated value is 2*nrd; halve it before the gcd
    # test, otherwise the condition is unsatisfiable whenever N is even
    B = floor(norm_bound**2 / (p * (log(p)**2)))
    vlist = lattice.LatticeEnumeration(L, B, condition=lambda newN: gcd(newN // 2, N) == 1, num_vectors=num_vectors)
    assert vlist, "No solution found in LatticeEnumeration"
    v = vlist[randint(0, len(vlist)-1)]
    x = sum(c * b for c, b in zip(v, O0.basis()))
    Nx = x.reduced_norm()
    assert Nx < B
    assert alpha2 * x * alpha1.conjugate() in O0 * N
    assert gcd(Nx, N) == 1

    # construct the lattice L = I1 \cap (O0 * x + Z)
    L1 = IntegralLattice(Gram, [vector(b) * Qinv for b in I1.basis()])
    Ox = IntegralLattice(Gram, [vector(b*x) * Qinv for b in O0.basis()])
    OxZ = Ox.overlattice([vector(O0(1)) * Qinv])
    L = IntegralLattice(Gram, L1.intersection(OxZ).basis())

    B = N * norm_bound
    vlist = lattice.LatticeEnumeration(L, B, condition=lambda newN: is_pseudoprime(ZZ(newN/(2*N))), num_vectors=num_vectors)
    assert vlist, "No solution found in LatticeEnumeration"
    v = vlist[randint(0, len(vlist)-1)]
    beta1 = sum(c * b for c, b in zip(v, O0.basis()))
    newN = ZZ(beta1.reduced_norm() / N)

    assert beta1 in I1
    beta2 = x * beta1 * x.conjugate() / Nx
    assert beta2 in I2

    return EquivalentIdeal(I1, beta1), EquivalentIdeal(I2, beta2), beta1, beta2, newN

def deltaKLPTforSign(Icom, IskIchl, l, e, norm_bound,
                    KLPT_margin=40, EISN_vec_bound=20, SA_loop_bound=1000):
    assert Icom.left_order() == IskIchl.left_order()
    _, qi, qj, qk = Icom.quaternion_algebra().basis()
    O = Icom.left_order()
    p = Icom.quaternion_algebra().discriminant()

    # bound for the original KLPT
    B1 = ceil(p**(0.5))
    B2 = ceil(p**(2.5)*log(p))

    found = False
    while not found:
        n1 = random_prime(2**KLPT_margin*B1, lbound=B1, proof=False)
        n2 = random_prime(2**KLPT_margin*B2, lbound=B2, proof=False)
        J1, _, found = KLPT(Icom, n1, n2)
        J2, alpha2, found2 = KLPT(IskIchl, n1, n2)
        found = found and found2
    assert norm(J1) == norm(J2) == n1*n2

    # randomize the class of (J_1, J_2)
    N = n1*n2
    r = randint(0, p)
    if r < p:
        alpha = 1 + r*qi
    else:
        alpha = qi
    n = alpha.reduced_norm()
    c = 0
    if n % 2 == 0:
        c += 1
    while not is_pseudoprime(n + c**2*p):
        c += 2
    alpha += c*qj
    assert alpha.reduced_norm() == n + c**2*p == alpha.conjugate().reduced_norm()
    J1 = EquivalentIdeal(J1, N*alpha)
    J2 = EquivalentIdeal(J2, N*alpha.conjugate())
    alpha2 = alpha.conjugate() * alpha2   # keep J2 == EquivalentIdeal(IskIchl, alpha2)
    N = N * (n + c**2*p)

    C, D = 0, 0
    while True:
        NCD = None
        while NCD is None or N > norm_bound or kronecker(l**e, N) != kronecker(NCD, N):
            B = max(norm_bound, EquivalentIdealsWithSameNormBound(p, N))
            J1, J2, _, beta2, newN = EquivalentIdealsWithSameNorm(J1, J2, N, B, EISN_vec_bound)
            assert norm(J1) == norm(J2) == newN
            alpha2 = beta2*alpha2 / N
            N = newN
            beta1 = SmallGenerator(J1)
            beta2 = SmallGenerator(J2)
            C, D = IdealModConstraint(O, qj, qk, beta2, beta1, N)
            NCD = p * (C**2 + D**2)
        assert J2 == EquivalentIdeal(IskIchl, alpha2)
        
        def is_cyclic(nu):
            if nu / 2 in O:
                return False
            O1 = J1.intersection(O*nu).right_order()
            O2 = J2.right_order()
            gamma = O2.isomorphism_to(O1, conjugator=True)
            assert J1.intersection(O*nu) * gamma.inverse() * J2.conjugate() * gamma == O * gamma
            return (alpha2.conjugate() * gamma) / 2 not in O

        nu, found = StrongApproximation(O, N, C, D, l**e, SA_loop_bound, condition=is_cyclic)
        if found:
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
