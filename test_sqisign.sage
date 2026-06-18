from sqisign import SQIsign
from quaternion import SmallGenerator

e = 248
f = 5
p = 2^e * f - 1  
lam = 128

SQIsign_instance = SQIsign(p, e, f, lam)
E0withEnd = SQIsign_instance.E0withEnd

for _ in range(10):
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
