from sqisign import SQIsign
from quaternion import SmallGenerator

e = 248
f = 5
p = 2^e * f - 1  
lam = 128

SQIsign_instance = SQIsign(1)
E0withEnd = SQIsign_instance.E0withEnd

for _ in range(10):
    P, Q = E0withEnd.P, E0withEnd.Q
    Pd, Qd = 2**(e-lam)*P, 2**(e-lam)*Q
    c = randint(0, 2**lam - 1)
    R = Pd + c*Qd
    I = E0withEnd.KernelToIdeal(1, c, lam)
    alpha = SmallGenerator(I)
    alphaP, alphaQ = E0withEnd.quaternion_action(alpha) 
    alphaP, alphaQ = 2**(e-lam)*alphaP, 2**(e-lam)*alphaQ
    assert (alphaP + c*alphaQ).is_zero(), "KernelToIdeal failed"
    print("KernelToIdeal verified successfully")

    sk, pk = SQIsign_instance.Keygen()
    print("Key generation successful")
    SQIsign_instance.Sign(sk, pk, b"Test message")
