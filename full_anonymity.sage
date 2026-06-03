from quaternion import *
from attack import *

proof.all(False)

def success_rate(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound, trials=100):
    success_simple = 0
    success_signing = 0
    for _ in range(trials):
        if simple_pullback(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
            success_simple += 1
        if simulate_signing(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound):
            success_signing += 1
        print(f"Trial {_+1}/{trials} completed.")
    return success_simple / trials, success_signing / trials

def make_instance(p, lam):
    assert p % 4 == 3
    B.<qi, qj, qk> = QuaternionAlgebra(QQ, -1, -p)
    O0 = B.quaternion_order([1, qi, (qi + qj)/2, (1 + qk)/2])
    assert O0.is_maximal()
    Nsk = random_prime(ceil(p**4))
    Ncom = random_prime(ceil(p**4))
    Nchl = 2**lam
    exp_rsp = ZZ(ceil(log(p, 2)) * 4 + 10)
    Nbound = 2**(ceil(log(p, 2)) + 5)

    Isk, _ = RandomFixedNormIdeal(O0, Nsk)
    Icom, _ = RandomFixedNormIdeal(O0, Ncom)
    Ichl, _ = RandomFixedNormIdeal(O0, Nchl)
    IskIchl = Isk.intersection(Ichl)

    L, nu = deltaKLPTforSign(Icom, IskIchl, 2, exp_rsp, Nbound)
    I = L + O0 * 2**exp_rsp
    assert I.is_principal()
    O = L.right_order()
    Od = IskIchl.right_order()
    beta = O.isomorphism_to(Od, conjugator=True)
    Iall = IskIchl * beta.inverse() * L.conjugate() * beta
    assert Iall.is_principal()
    IskIchlIrsp = Iall + O0 * Nsk * Nchl * 2**exp_rsp
    assert not IskIchlIrsp * B(1/2) in O0
    return IskIchlIrsp, Nsk, Nchl, 2**exp_rsp, Nbound

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
IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound = make_instance(p, lam)
success_simple, success_signing = success_rate(IskIchlIrsp, Nsk, Nchl, Nrsp, N_bound, 1000)
print(f"Success rate of simple pullback: {success_simple}")
print(f"Success rate of signing simulation: {success_signing}")