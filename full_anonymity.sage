from importlib import util

from quaternion import *
from attack import *

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
        com, _ = inst.super().RecoverCommitment(Pk[i], Chl[i], Rsp[i])
        chl = inst.super().Hash(message + util.j_invariant_to_bytes(com) + mPk)
        Com.append(com)
        if i < n_parties - 1:
            Chl.append(chl)


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