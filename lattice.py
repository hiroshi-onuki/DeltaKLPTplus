from sage.all import ( # type: ignore
    ZZ,
    QQ,
    gcd,
    vector,
    floor,
    ceil,
    sqrt,
    matrix,
)
from sage.modules.free_module_integer import IntegerLattice # type: ignore

def EuclideanNorm(v):
    v = vector(ZZ, v)
    return v.dot_product(v)

def ShortBasisDim2Euclidean(b0, b1):
    beta0, beta1 = vector(ZZ, b0), vector(ZZ, b1)
    if EuclideanNorm(beta0) < EuclideanNorm(beta1):
        beta0, beta1 = beta1, beta0

    gamma = beta0
    while True:
        r = ZZ(floor(beta0.dot_product(beta1) / EuclideanNorm(beta1) + ZZ(1)/ZZ(2)))
        gamma = beta0 - r*beta1
        if EuclideanNorm(gamma) < EuclideanNorm(beta1):
            beta0, beta1 = beta1, gamma
        else:
            break

    if EuclideanNorm(gamma) < EuclideanNorm(beta0):
        beta0 = gamma
    return beta1, beta0

# generator: yields lattice vectors v = close + x*b0 + y*b1 with ||t - v||^2 <= B, in the same order as
# before, but lazily and with plain Python integers (the caller usually stops after a few vectors).
def EnumerateCloseVectorsDim2Euclidean(b0, b1, t, m, B):
    b0, b1 = vector(ZZ, b0), vector(ZZ, b1)
    t = vector(ZZ, t)
    L = IntegerLattice([b0, b1])
    close = L.approximate_closest_vector(t)
    m = int(m)
    B = int(B)
    if m <= 0:
        return

    d = t - close
    b00, b01 = int(b0[0]), int(b0[1])
    b10, b11 = int(b1[0]), int(b1[1])
    d0, d1 = int(d[0]), int(d[1])
    c0, c1 = int(close[0]), int(close[1])
    t0, t1 = int(t[0]), int(t[1])
    a = b00*b00 + b01*b01
    h = b00*b10 + b01*b11
    c = b10*b10 + b11*b11
    delta = a*c - h*h
    if delta <= 0:
        raise ValueError("EnumerateCloseVectorsEuclidean: degenerate lattice basis")

    det = b00*b11 - b01*b10
    y_num = b00*d1 - b01*d0
    y_den = det
    if y_den < 0:
        y_num = -y_num
        y_den = -y_den

    # From min_x ||d - x*b0 - y*b1||^2 = delta/a * (y-y0)^2.
    y_radius = int(ZZ((a*B) // delta).isqrt()) + 2
    y_min = y_num // y_den - y_radius                  # floor(y_num/y_den) - y_radius
    y_max = -((-y_num) // y_den) + y_radius            # ceil(y_num/y_den) + y_radius

    tries = 0
    db0 = d0*b00 + d1*b01
    db1 = d0*b10 + d1*b11
    nd = d0*d0 + d1*d1
    for y in range(y_min, y_max + 1):
        if tries >= m:
            return
        K = c*y*y - 2*db1*y + nd - B
        Lx = h*y - db0
        D = Lx*Lx - a*K
        if D < 0:
            continue
        x_radius = int(ZZ(D).isqrt()) + 2
        x_min = (-Lx - x_radius) // a - 1              # floor((-Lx - x_radius)/a) - 1
        x_max = -((Lx - x_radius) // a) + 1            # ceil((-Lx + x_radius)/a) + 1
        for x in range(x_min, x_max + 1):
            if tries >= m:
                return
            tries += 1
            v0 = c0 + x*b00 + y*b10
            v1 = c1 + x*b01 + y*b11
            e0 = t0 - v0
            e1 = t1 - v1
            if e0*e0 + e1*e1 <= B:
                yield vector(ZZ, [v0, v1])


# return coefficients q_i,j s.t.
# Q(sum_i x_i*basis_i) = sum_i q_i,i*(x_i + sum_{j > i} q_i,j*x_j)^2.
def MakeQuatraticForm(basis, quadratic_form):
    n = len(basis)
    C = matrix(QQ, n, n)
    q = matrix(QQ, n, n)

    for i in range(n):
        C[i, i] = quadratic_form(basis[i], basis[i])
        for j in range(i + 1, n):
            C[i, j] = quadratic_form(basis[i], basis[j])

    for i in range(n):
        q[i, i] = C[i, i] - sum(q[k, k] * q[k, i]**2 for k in range(i))
        for j in range(i + 1, n):
            q[i, j] = (C[i, j] - sum(q[k, k] * q[k, i] * q[k, j] for k in range(i))) / q[i, i]
    return q

# floor(sqrt(s) + U) for a nonnegative rational s and a rational U, computed exactly.
def _floor_sqrt_plus(s, U):
    num = s.numerator()
    den = s.denominator()
    fs = (num * den).isqrt() // den          # floor(sqrt(s))
    k = ZZ((fs + U).floor()) + 2             # guaranteed >= the true value
    while True:
        t = k - U
        if t <= 0 or t * t <= s:
            return k
        k -= 1

# enumerate vectors of the (already reduced) basis with norm < B
def LatticeEnumeration(L, B, condition, num_vectors):
    red_basis = [vector(ZZ, b) for b in L.LLL().basis()]
    n = len(red_basis)
    G = [[ZZ(red_basis[a].inner_product(red_basis[b])) for b in range(n)] for a in range(n)]
    q = MakeQuatraticForm(list(range(n)), lambda a, b: G[a][b])
    zero_vec = vector(ZZ, [0] * len(red_basis[0]))

    S = [0] * n
    U = [0] * n
    upper = [ZZ(0)] * n
    x = [ZZ(0)] * n
    S[-1] = ZZ(B)

    i = n - 1
    s = S[i] / q[i, i]
    upper[i] = _floor_sqrt_plus(s, -U[i])
    x[i] = -_floor_sqrt_plus(s, U[i]) - 1

    ret = []
    while True:
        x[i] += 1
        while i < n and x[i] > upper[i]:
            i += 1
            if i == n:
                return ret
            x[i] += 1

        if i > 0:
            S[i - 1] = S[i] - q[i, i] * (x[i] + U[i])**2
            i -= 1
            U[i] = sum(q[i, j] * x[j] for j in range(i + 1, n))
            assert q[i, i] > 0
            s = S[i] / q[i, i]
            upper[i] = _floor_sqrt_plus(s, -U[i])
            x[i] = -_floor_sqrt_plus(s, U[i]) - 1
            continue

        if any(x):
            g = gcd(x)
            coeffs = [ZZ(c // g) for c in x]
            newN = sum(coeffs[a] * G[a][b] * coeffs[b] for a in range(n) for b in range(n))
            if condition(newN):
                alpha = sum((coeffs[j] * red_basis[j] for j in range(n)), zero_vec)
                ret.append(alpha)
                if len(ret) >= num_vectors:
                    return ret
        else:
            return ret
