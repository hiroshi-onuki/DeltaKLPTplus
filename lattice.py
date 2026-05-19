from sage.all import (
    ZZ,
    QQ,
    gcd,
    vector,
    floor,
    ceil,
    sqrt,
    matrix,
    log,
)

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

def ClosestVectorDim2Euclidean(beta1, beta0, t):
    beta1, beta0, t = vector(ZZ, beta1), vector(ZZ, beta0), vector(ZZ, t)
    mu1 = EuclideanNorm(beta1)*beta0 - beta0.dot_product(beta1)*beta1
    r0 = ZZ(floor(mu1.dot_product(t)*EuclideanNorm(beta1) / EuclideanNorm(mu1) + ZZ(1)/ZZ(2)))
    residual = t - r0*beta0
    r1 = ZZ(floor(beta1.dot_product(residual) / EuclideanNorm(beta1) + ZZ(1)/ZZ(2)))
    return t - (residual - r1*beta1)

def EnumerateCloseVectorsDim2Euclidean(b0, b1, t, close, m, B):
    b0, b1 = vector(ZZ, b0), vector(ZZ, b1)
    t = vector(ZZ, t)
    close = vector(ZZ, close)
    m = ZZ(m)
    B = ZZ(B)
    if m <= 0:
        return

    d = t - close
    a = EuclideanNorm(b0)
    h = b0.dot_product(b1)
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
    y_radius = floor(sqrt((a*B) // delta)) + 2
    y_min = ZZ(floor(y_num / y_den)) - y_radius
    y_max = ZZ(ceil(y_num / y_den)) + y_radius

    tries = ZZ(0)
    db0 = d.dot_product(b0)
    db1 = d.dot_product(b1)
    nd = EuclideanNorm(d)
    ret = []
    for y in range(y_min, y_max + 1):
        if tries >= m:
            break
        K = c*y**2 - 2*db1*y + nd - B
        Lx = h*y - db0
        D = Lx**2 - a*K
        if D < 0:
            continue
        x_radius = floor(sqrt(D)) + 2
        x_min = ZZ(floor((-Lx - x_radius) / a)) - 1
        x_max = ZZ(ceil((-Lx + x_radius) / a)) + 1
        for x in range(x_min, x_max + 1):
            if tries >= m:
                break
            tries += 1
            v = close + x*b0 + y*b1
            if EuclideanNorm(t - v) <= B:
                ret.append(v)
    return ret

def lattice_inner_product(L, x, y):
    return QQ(x.inner_product(y))

def lattice_norm(L, x):
    return lattice_inner_product(L, x, x)

# return coefficients q_i,j s.t.
# Q(sum_i x_i*basis_i) = sum_i q_i,i*(x_i + sum_{j > i} q_i,j*x_j)^2.
def make_quadratic_form_coeffs(basis, quadratic_form):
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

# return alpha in L s.t. norm(alpha) < 2^a * norm(L) and condition(norm(alpha)/norm(L)).
def element_for_response(L, a, condition):
    red_basis = [vector(ZZ, b) for b in L.LLL().basis()]
    q = make_quadratic_form_coeffs(red_basis, lambda x, y: lattice_inner_product(L, x, y))

    n = len(red_basis)
    S = [0] * n
    U = [0] * n
    upper = [ZZ(0)] * n
    x = [ZZ(0)] * n
    S[-1] = ZZ(2**a)

    i = n - 1
    Z = sqrt(S[i] / q[i, i])
    upper[i] = ZZ(floor(Z - U[i]))
    x[i] = ZZ(ceil(-Z - U[i])) - 1

    while True:
        x[i] += 1
        while i < n and x[i] > upper[i]:
            i += 1
            if i == n:
                return vector(ZZ, [0] * len(red_basis[0])), ZZ(0), False
            x[i] += 1

        if i > 0:
            S[i - 1] = S[i] - q[i, i] * (x[i] + U[i])**2
            i -= 1
            U[i] = sum(q[i, j] * x[j] for j in range(i + 1, n))
            assert q[i, i] > 0
            Z = sqrt(S[i] / q[i, i])
            upper[i] = ZZ(floor(Z - U[i]))
            x[i] = ZZ(ceil(-Z - U[i])) - 1
            continue

        if any(x):
            g = gcd(x)
            coeffs = [ZZ(c // g) for c in x]
            alpha = sum((coeffs[j] * red_basis[j] for j in range(n)), vector(ZZ, [0] * len(red_basis[0])))
            newN = ZZ(lattice_norm(L, alpha))
            if condition(newN):
                return alpha, newN, True
        else:
            return vector(ZZ, [0] * len(red_basis[0])), ZZ(0), False
