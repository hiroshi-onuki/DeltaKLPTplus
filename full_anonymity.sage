from quaternion import *
from attack import *
from ring_sqisign import RingSQIsign
import util

proof.all(False)

def success_num(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound, trials=100):
    success_simple = 0
    success_signing = 0
    for _ in range(trials):
        if simple_pullback(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
            success_simple += 1
        if simulate_signing(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
            success_signing += 1
        print(f"Trial {_+1}/{trials} completed.\r", end="")
    return success_simple, success_signing

def ChallengeToIdeal(inst, sk, chl):
    """
    Return the left O0-ideal Ichl of norm 2^e_chl such that Isk ∩ Ichl corresponds to
    phi_chl ∘ phi_sk (the same computation as the beginning of SQIsign.Respond).
    """
    Isk, Msk = sk
    a, b = vector([1, chl]) * Msk.inverse()
    return Isk.intersection(inst.E0withEnd.KernelToIdeal(a, b, inst.e_chl))

def ResponseToIdeal(inst, sk, chl, rsp):
    """
    Return the left O0-ideal corresponding to phi_rsp ∘ phi_chl ∘ phi_sk, where phi_rsp ∘ phi_chl
    is the chain of five 2-power isogenies encoded by (chl, rsp) (see SQIsign.RecoverCommitment).

    The chain is walked on ideals instead of curves: at each step the current ideal J is replaced
    by an equivalent odd-norm ideal Im, whose isogeny (IdealToIsogeny) lands on the same normalized
    curve, so the kernel c*P + Q (or P + c*Q) given on the deterministic basis can be pulled back to
    E0 and pushed through J by the same beta.
    """
    E0 = inst.E0withEnd
    e = E0.e
    e0 = inst.e_chl + inst.e_rsp - 4*e
    Isk, Msk = sk
    c0without_chl, c1b, c1f, c2b, c2f, isP1b, isP1f, isP2b, isP2f = rsp
    c0 = c0without_chl * 2**inst.e_chl + chl

    # first step: kernel 2^(e-e0) * (Ppk + c0*Qpk) on Epk, written on the E0-basis via Msk
    a, b = vector([1, c0]) * Msk.inverse()
    J = Isk.intersection(E0.KernelToIdeal(a, b, e0))

    for c, isP in [(c1b, isP1b), (c1f, isP1f), (c2b, isP2b), (c2f, isP2f)]:
        Im, beta, _ = EquivalentPrimeIdeal(J, inst.sec_lambda)     # Im = J * conj(beta) / norm(J)
        Em, Pm, Qm = E0.IdealToIsogeny(Im)
        Em, (Pm, Qm) = inst._normalize_curve(Em, (Pm, Qm))
        Pmd, Qmd = inst._deterministic_torsion_basis(Em, e)
        Mm = util.BiDLP_matrix_power_two(Pm, Qm, Pmd, Qmd, e)
        v = vector([c, 1] if isP else [1, c]) * Mm.inverse()
        IK = E0.KernelToIdeal(v[0], v[1], e)                       # kernel on E0 pulled back by Im
        J = Im.intersection(IK) * (beta / norm(Im))                 # phi_K ∘ phi_J
    return J


def make_instance(e, f, lam, n_parties):
    message = b"Test message"
    p = 2**(2*lam) * f - 1
    inst = RingSQIsign(e, f, lam, n_parties)

    # generate keys and signature
    Pk, Sk = inst.Keygen()
    mPk = b''.join([util.j_invariant_to_bytes(pk) for pk in Pk])
    chl_first, Rsp = inst.Sign(Pk, Sk[0], 0, message)

    # recover commitments and challenges
    Chl = []
    Com = []
    Chl.append(chl_first)
    chl = chl_first
    for i in range(n_parties):
        com, _ = inst.RecoverCommitment(Pk[i], Chl[i], Rsp[i])
        chl = inst.Hash(message + util.j_invariant_to_bytes(com) + mPk)
        Com.append(com)
        if i < n_parties - 1:
            Chl.append(chl)

    IskIchls = [ChallengeToIdeal(inst, Sk[i], Chl[i]) for i in range(n_parties)]
    IskIchlIrsps = [ResponseToIdeal(inst, Sk[i], Chl[i], Rsp[i]) for i in range(n_parties)]
    for IskIchl, IskIchlIrsp, com in zip(IskIchls, IskIchlIrsps, Com):
        assert IskIchlIrsp.is_submodule(IskIchl)
        assert norm(IskIchlIrsp) == inst.Dmix * 2**(inst.e_chl + inst.e_rsp)
        Im, _, _ = EquivalentPrimeIdeal(IskIchlIrsp, lam)
        assert inst.E0withEnd.IdealToIsogeny(Im)[0].j_invariant() == com.j_invariant()
    return inst, Pk, Sk, IskIchl, IskIchlIrsp

def make_p(lam):
    f = 1
    p = 2**(2*lam) * f - 1
    while not is_prime(p):
        f += 1
        p = 2**(2*lam) * f - 1
    return p

lam = 5
p = make_p(lam)
print(f"Using p = {p}")

O0, IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound = make_instance(p, lam)

n_trials = 10000

# attack for actual signing
success_simple, success_signing = success_num(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound, n_trials)
print(f"Success rate of simple pullback: {success_simple}/{n_trials}")
print(f"Success rate of signing simulation: {success_signing}/{n_trials}")

# attack for random isogeny chain
Isk = IskIchlIrsp + O0 * Nsk
IchlIrsq_dummy, _ = RandomFixedNormIdeal(O0, Nchl * Nrsp)
IskIchlIrsp_dummy = Isk.intersection(IchlIrsq_dummy)
assert norm(IskIchlIrsp_dummy) == Nsk * Nchl * Nrsp
success_simple_dummy, success_signing_dummy = success_num(IskIchlIrsp_dummy, Nsk, Nchl, Nrsp, N_bound, n_trials)
print(f"Success rate of simple pullback for random isogeny chain: {success_simple_dummy}/{n_trials}")
print(f"Success rate of signing simulation for random isogeny chain: {success_signing_dummy}/{n_trials}")