from quaternion import *
from lattice import *

def ListAllEquivalentIdealsWithSameNorm(I1, I2, N):
    """
    List (J1, J2) s.t. J1 ~ I1, J2 ~ I2 and norm(J1) = norm(J2) <= N
    """
    O = I1.left_order()
    p = O.discriminant()
    Gram = matrix(ZZ, 4, 4, [(b1*b2.conjugate()).reduced_trace() for b1 in O.basis() for b2 in O.basis()])
    N1 = norm(I1)
    N2 = norm(I2)
    retList = []

    L1 = IntegralLattice(Gram, [vector(b) * O.basis_matrix().inverse() for b in I1.basis()])
    S1 = LatticeEnumeration(L1, 2*N*N1, lambda x: is_prime(ZZ(x/(2*N1))), N**2)
    S1 = [sum(c * b for c, b in zip(v, O.basis())) for v in S1]
    print("S1 size:", len(S1))
    T1 = {(alpha.reduced_norm() / N1): alpha for alpha in S1}

    L2 = IntegralLattice(Gram, [vector(b) * O.basis_matrix().inverse() for b in I2.basis()])
    S2 = LatticeEnumeration(L2, 2*N*N2, lambda x: is_prime(ZZ(x/(2*N2))), N**2)
    print("S2 size:", len(S2))
    T2 = {(alpha.reduced_norm() / N2): alpha for alpha in S2}
    for t in T2.keys():
        if T1.has_key(t):
            alpha1 = T1[t]
            alpha2 = T2[t]
            J1 = EquivalentIdeal(I1, alpha1)
            J2 = EquivalentIdeal(I2, alpha2)
            retList.append((J1, J2))
    return retList

def ListByEquivalentIdealsWithSameNorm(I1, I2, norm_bound, vec_bound, num_loops):
    # bound for the original KLPT
    p = I1.left_order().discriminant()
    B1 = ceil(p**(0.5))
    B2 = ceil(p**(2.5)*log(p))
    KLPT_margin = 2**40
    found = False
    while not found:
        n1 = random_prime(KLPT_margin*B1, lbound=B1, proof=False)
        n2 = random_prime(KLPT_margin*B2, lbound=B2, proof=False)
        J1, _, found = KLPT(I1, n1, n2)
        J2, _, found2 = KLPT(I2, n1, n2)
        found = found and found2
    assert norm(J1) == norm(J2) == n1*n2

    N = n1*n2
    retList = []
    for _ in range(num_loops):
        B = max(norm_bound, ceil(p**(3/4) * N**(1/4) * (3/4*log(p) + 1/4*log(N))))
        J1, J2, _, _, N = EquivalentIdealsWithSameNorm(J1, J2, N, B, vec_bound)
        retList.append((J1, J2))
    return retList


p = 103
B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
assert O0.is_maximal()
N = random_prime(ceil(p**4))
I1, _ = RandomFixedNormIdeal(O0, N)
I2, _ = RandomFixedNormIdeal(O0, N)

bound = 10 * p * ceil(log(p)^(4/3))
L = ListAllEquivalentIdealsWithSameNorm(I1, I2, bound)
print(len(L))
Ld = ListByEquivalentIdealsWithSameNorm(I1, I2, bound, 100, 10*len(L))
print(len(Ld))