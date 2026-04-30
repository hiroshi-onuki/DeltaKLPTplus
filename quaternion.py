from sage.all import (
    ZZ,
    GF,
    log,
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
)
from sage.rings.factorint import factor_trial_division
from sage.modules.free_module_integer import IntegerLattice

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

# Return a quaternion nu such that nu = Cj + Dk mod N and Nrd(nu) = l^e.
def StrongApproximation(O0, N, C, D, l, max_cnt=10):
    p = O0.discriminant()
    e = floor(3 * log(N, l) + log(p, l))
    Nrd_mu = p * (C**2 + D**2)
    if (kronecker(Nrd_mu, N) == 1 and e % 2 == 1) or (kronecker(Nrd_mu, N) == -1 and e % 2 == 0):
        e += 1

    while True:
        lam = ZZ(sqrt(GF(N)(l**e)/GF(N)(Nrd_mu)))
        rhs = ZZ((l**e - lam**2 * Nrd_mu) / N)

        for _ in range(max_cnt):
            c = GF(N).random_element()
            d = ZZ((GF(N)(rhs) / GF(N)(2*p*lam) - C*c) / GF(N)(D))
            c = ZZ(c)

            x = ZZ(-GF(N)(C)/GF(N)(D))
            L = IntegerLattice([[N, N*x], [0, N**2]])
            v = L.approximate_closest_vector([-lam*C - N*c, -lam*D - N*d])
            Nc = N*c + v[0]
            Nd = N*d + v[1]

            tmp = ZZ((l**e - p*((lam*C + N*c)**2 + (lam*D + N*d)**2)) / N**2)
            a, b = SumOf2Squares(tmp)
            if a is not None and b is not None:
                nu = O0([N*a, N*b, lam*C + N*c, lam*D + N*d])
                assert nu.reduced_norm() == l**e
                return nu
        e += 2

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
