from sage.all import (
    ZZ,
    floor,
    sqrt,
    is_prime,
    randint,
    valuation,
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

    # set Icom as an ideal with prime norm equivalent to IskIchlIrsp
    reduced_basis = LLLBasis(IskIchlIrsp)
    Icom, alpha, Ncom = EquivalentRandomPrimeIdeal(IskIchlIrsp, lambda N: gcd(N, norm(IskIchlIrsp)) == 1)
    IcomIrsp = O0 * alpha.conjugate() + O0 * Ncom * Nrsp
    assert IcomIrsp < Icom

    # call delta-KLPT
    e_rsp = valuation(Nrsp, 2)
    assert Nrsp == 2**e_rsp
    L, _ = deltaKLPTforSign(Icom, IskIchl, 2, e_rsp, N_bound)

    # compute the pullback of Irsp
    N = norm(L) / Nrsp
    Icom_d = L + O0 * N
    Ocom = Icom.right_order()
    Ocom_d = Icom_d.right_order()
    alpha = Ocom_d.isomorphism_to(Ocom, conjugator=True)
    I = Icom * alpha.inverse() * Icom_d.conjugate() * alpha
    assert I.is_principal()
    assert I == O0 * alpha
    assert Icom_d == EquivalentIdeal(Icom, alpha)
    IcomIrsp = IcomIrsp * (alpha.conjugate() / Ncom)
    assert IcomIrsp < Icom_d
    Irsp_pullback = IcomIrsp + O0 * Nrsp
    return Irsp_pullback.is_principal()

