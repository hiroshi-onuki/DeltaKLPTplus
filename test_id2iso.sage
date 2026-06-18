from id2iso import IdealToIsogeny
from quaternion import RandomFixedNormIdeal, SmallGenerator
from sqisign import SQIsign

e = 248
f = 5
p = 2^e * f - 1  
lam = 128

SQIsign_instance = SQIsign(p, e, f, lam)
E0withEnd = SQIsign_instance.E0withEnd

for _ in range(10):
    N = random_prime(10*p)
    I, _ = RandomFixedNormIdeal(E0withEnd.order, N)
    print(f"Generated ideal I with norm {N}")

    EI, PI, QI = IdealToIsogeny(E0withEnd, I)
    print("Isogeny computed successfully")

    SQIsign_instance._deterministic_torsion_basis(EI, e)
    print("Deterministic torsion basis computed successfully")

    P, Q = E0withEnd.P, E0withEnd.Q
    c = randint(0, 2**e - 1)
    R = P + c*Q
    I = E0withEnd.KernelToIdeal(1, c)
    alpha = SmallGenerator(I)
    alphaP, alphaQ = E0withEnd.quaternion_action(alpha) 
    assert (alphaP + c*alphaQ).is_zero(), "KernelToIdeal failed"
    print("KernelToIdeal verified successfully")

    sk, pk = SQIsign_instance.Keygen()
    SQIsign_instance.Sign(sk, pk, b"Test message")
