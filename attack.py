from sage.all import (
    ZZ,
    floor,
    sqrt,
    is_prime,
    randint,
    random_prime,
    norm,
)

from quaternion import *

def simple_pullback(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
    O0 = IskIchlIrsp.left_order()
    IskIchl = IskIchlIrsp + O0 * Nsk*Nchl
    
    # find an ideal equivalent to IskIchl s.t. norm(I) is a prime <= N_bound
    reduced_basis = LLLBasis(IskIchl)
    Bs = [ZZ(floor(sqrt(N_bound * Nsk * Nchl / (b.reduced_norm())))) for b in reduced_basis]
    alpha = O0(0)
    I = IskIchl
    N = N_bound + 1
    while N > N_bound or not is_prime(N):
        cs = [randint(-B, B) for B in Bs]
        alpha = sum(c * b for c, b in zip(cs, reduced_basis))
        I = EquivalentIdeal(IskIchl, alpha)
        N = ZZ(norm(I))
    
    # compute the pullback of Irsp by I
    IIrsp = IskIchlIrsp * (alpha.conjugate() / (Nsk * Nchl))
    Irsp_pullback = IIrsp + O0 * Nrsp

    return Irsp_pullback.is_principal()

def simulate_signing(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
    O0 = IskIchlIrsp.left_order()
    _, _, qj, qk = IskIchlIrsp.quaternion_algebra().basis()
    p = IskIchlIrsp.quaternion_algebra().discriminant()
    IskIchl = IskIchlIrsp + O0 * Nsk*Nchl

    # set Icom as a short ideal equivalent to IskIchlIrsp
    reduced_basis = LLLBasis(IskIchlIrsp)
    Icom = EquivalentIdeal(IskIchlIrsp, reduced_basis[0])

    # do as the same as deltaKLPTforSign
    found = False
    while not found:
        B1 = ceil(p**(0.7))
        B2 = ceil(p**(2.8))
        found = False
        while not found:
            n1 = random_prime(B1)
            n2 = random_prime(B2)
            if n1 < p**(0.5) or n2 < p**(2.5):
                continue
            J1, _, found = KLPT(Icom, n1, n2)
            J2, alpha2, found2 = KLPT(IskIchl, n1, n2)
            found = found and found2
            B1 *= 2
            B2 *= 2
        assert norm(J1) == norm(J2) == n1*n2
        N = n1*n2
        NCD = None
        cnt = 0
        while cnt < 10 and (NCD is None or N > N_bound or kronecker(Nrsp, N) != kronecker(NCD, N)):
            J1, J2, _, beta2, newN = EquivalentIdealsWithSameNorm(J1, J2, N, 100)
            alpha2 = beta2*alpha2 / N
            N = newN
            beta1 = SmallGenerator(J1)
            beta2 = SmallGenerator(J2)
            cnt += 1
            C, D = IdealModConstraint(O0, qj, qk, beta2, beta1, N)
            NCD = p * (C**2 + D**2)
        if N > N_bound or kronecker(Nrsp, N) != kronecker(NCD, N):
            continue
        _, found = StrongApproximation(O0, N, C, D, Nrsp, 100, condition=lambda nu: not nu*alpha2.conjugate()/(2*N) in O0)
    J2Irsp = IskIchlIrsp * (alpha2.conjugate() / (Nsk * Nchl))
    Irsp_pullback = J2Irsp + O0 * Nrsp
    return Irsp_pullback.is_principal()
